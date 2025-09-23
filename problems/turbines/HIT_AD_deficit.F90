! Template for PadeOps
!
! This problem allows for computing deficit budgets
! in the HIT + shear + AD problem (see HIT_AD_interact and
! HIT_AD_shear)
!
! Kirby Heck
! 2024 December 14

! NOTE: Initializion references HIT_shear_files
#include "HIT_shear_files/initialize.F90"
#include "HIT_shear_files/temporalHook.F90"

program HIT_deficit
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use frozen_igrid_mod, only: frozen_igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message, message_min_max, GracefulExit
    use constants, only: one, zero
    use reductions,         only: p_maxval, p_minval
    use HIT_shear_parameters
    use fof_mod, only: fof
    use budgets_time_avg_mod, only: budgets_time_avg
    use budgets_vol_avg_mod, only: budgets_vol_avg
    use budgets_time_avg_deficit_mod, only: budgets_time_avg_deficit

    implicit none

    class(igrid), allocatable, target :: hit, adsim, emptysim  ! make these polymorphic so we can freeze the turbulence
    character(len=clen) :: inputfile, HIT_InputFile, AD_InputFile, Empty_InputFile, fof_dir, filoutdir
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg, budg_tavg_empty
    type(budgets_vol_avg)  :: budg_vavg
    type(budgets_time_avg_deficit) :: budg_tavg_deficit  ! added deficit budgets, can turn off in inputfile (ON by default)
    real(rkind), dimension(:,:,:), allocatable :: utarget, vtarget, wtarget
    real(rkind) :: dt1 = one, dt2 = one, dt3 = one, dt = one
    real(rkind) :: k_bandpass_left = 10.d0, k_bandpass_right = 64.d0, TI_xloc = 0
    real(rkind) :: TI_target = -1, TI_fact = -1, KIinv_TI = 0.5d0, Kp_TI = 1.d0, integral_err = zero, time_stop_TIcont = -1

    type(fof), dimension(:), allocatable :: filt
    integer, dimension(:), allocatable :: pid
    integer :: fid, nfilters = 2, tid_FIL_FullField = 75, tid_FIL_Planes = 4, TI_xid
    integer :: aniso_x = 1
    logical :: applyFilters = .false., freeze_HIT = .false., control_TI = .false., TI_at_rotor = .true., isStratified = .false.
    logical, parameter :: synchronize_RK_substeps = .true.

    namelist /concurrent/ HIT_InputFile, AD_InputFile, Empty_InputFile, InflowSpeed, &
        k_bandpass_left, k_bandpass_right, & 
        TI_target, TI_xloc, TI_fact, freeze_HIT, advect_shear, KIinv_TI, Kp_TI, time_stop_TIcont, TI_at_rotor, &
        isStratified
    namelist /FILTER_INFO/ applyfilters, nfilters, fof_dir, tid_FIL_FullField, tid_FIL_Planes, filoutdir

    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    ! read concurrent input file
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    read(unit=ioUnit, NML=FILTER_INFO)
    close(ioUnit)

    allocate(adsim, emptysim)
    if (freeze_HIT) then
        allocate(frozen_igrid :: hit)
    else
        allocate(igrid :: hit)
    end if

    call get_is_stratified(AD_inputfile)  ! whether to load T field in initialization - glean from &PHYSICS

    ! initialize igrid objects
    simulationID = 1
    call adsim%init(AD_InputFile, .true.)  ! initialize decomposition
    call adsim%start_io(.true.)  ! don't start IO
    call adsim%printDivergence()

    call mpi_barrier(mpi_comm_world, ierr)
    call message("Initialized PRIMARY simulation")

    call emptysim%init(Empty_InputFile, .False.)  ! do not initialize decomp
    emptysim%Am_I_Primary = .false.  ! not primary
    call emptysim%start_io(.true.)
    call emptysim%printDivergence()
    call message("Initialized EMPTY PRECURSOR simulation")

    ! check to make sure the simulations have the same grid, required for deficit budgets
    if (.not. all(adsim%mesh == emptysim%mesh)) call GracefulExit("EMPTY and AD simulations mesh dimensions must match", 234)
    call mpi_barrier(mpi_comm_world, ierr)

    simulationID = 2  ! HIT box
    call hit%init(HIT_InputFile, .false.)
    hit%Am_I_Primary = .false.
    call hit%start_io(.true.)
    call hit%printDivergence()
    call message("Initialized CONCURRENT HIT simulation")
    if (freeze_HIT) call message(1, "HIT targets are FROZEN")

    ! For anisotropic PRIMARY and EMPTY domains, we will need to declare an anisotropy factor in x
    aniso_x = nint(adsim%dx / hit%dx)
    if (abs((adsim%dx / hit%dx) - real(aniso_x)) > 1e-5) then
        call GracefulExit("Anisotropy factor must be an integer >= 1.", 211)
    else if (aniso_x .ne. 1) then
        call message(0, "PRIMARY grid is anisotropic, using aniso_x factor", aniso_x)
    end if

    call make_global_zaxis(adsim)  ! allocate the global-z axis variables
    nxfringe = min(nxadsim * aniso_x, nxhitsim)  ! determine domain range from HIT to use in fringe targets

    !!!!!!!!!!!!! decide whether to turn on the TI controller !!!!!!!!!!!!!
    control_TI = .false.
    if (TI_target > 0) then
        TI_fact = one
        TI_xid = minloc(abs(adsim%mesh(:,1,1,1) - TI_xloc), 1)  ! xid corresponding to TI sampling location
        control_TI = .true.
        call message(0, "TI controller activated")
        call message(1, "TI controller parameters")
        call message(2, "Kp", Kp_TI)
        call message(2, "KIinv", KIinv_TI)
        call message(1, "tracking x-location:", adsim%mesh(TI_xid,1,1,1))
        call message(1, "target TI: ", TI_target)
        if (TI_at_rotor) call message(1, "zmid_for_TI (+/- 0.5):", zmid_for_TI)
    else if (TI_fact >= 0) then
        call message(0, "TI controller not used")
        call message(1, "Using fixed TI gain/loss: ", TI_fact)
    else
        call message(0, "No TI settings provided, superimposing HIT fluctuations")
        TI_fact = one
    end if

    !!!!!!!!!!!!! allocate target cells for the fringe !!!!!!!!!!!!!
    allocate(utarget0(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(vtarget0(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(wtarget0(adsim%gpE%xsz(1), adsim%gpE%xsz(2), adsim%gpE%xsz(3)))
    if (adsim%isStratified) allocate(Ttarget0(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    call init_fringe_targets(AD_inputfile, adsim%mesh)  ! populates utarget0, vtarget0, wtarget0

    ! allocate moving (turbulent) targets
    allocate(utarget(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(vtarget(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(wtarget(adsim%gpE%xsz(1), adsim%gpE%xsz(2), adsim%gpE%xsz(3)))

    ! initialize turbulent fluctuations as zero
    utarget = utarget0
    vtarget = vtarget0
    wtarget = wtarget0

    ! initialize bandpass filter
    call hit%spectC%init_bandpass_filter(k_bandpass_left, k_bandpass_right, hit%cbuffzC(:,:,:,1), hit%cbuffyC(:,:,:,1))

    ! now initialize turbulent fringe targets
    if (adsim%usedoublefringex) then
        call message(0, "Setting double fringe targets")
        ! first fringe is re-laminarization
        call adsim%fringe_x1%associateFringeTargets(utarget0, vtarget0, wtarget0, Ttarget0)
        call adsim%fringe_x1%associateFringeTarget_scalar(Ttarget0)
        call emptysim%fringe_x1%associateFringeTargets(utarget0, vtarget0, wtarget0, Ttarget0)
        call emptysim%fringe_x1%associateFringeTarget_scalar(Ttarget0)

        ! second fringe is turbulent
        call adsim%fringe_x2%associateFringeTargets(utarget, vtarget, wtarget, Ttarget0)
        call adsim%fringe_x2%associateFringeTarget_scalar(Ttarget0)
        call emptysim%fringe_x2%associateFringeTargets(utarget, vtarget, wtarget, Ttarget0)
        call emptysim%fringe_x2%associateFringeTarget_scalar(Ttarget0)
    else
        call message(0, "Setting fringe targets")
        ! first (only) fringe is turbulent
        call adsim%fringe_x%associateFringeTargets(utarget, vtarget, wtarget, Ttarget0)
        call adsim%fringe_x%associateFringeTarget_scalar(Ttarget0)
        call emptysim%fringe_x%associateFringeTargets(utarget, vtarget, wtarget, Ttarget0)
        call emptysim%fringe_x%associateFringeTarget_scalar(Ttarget0)
    end if

    ! phaseshift turbulent fringe targets using the laminar fringe targets
    if (control_TI) call update_TI_fact(emptysim, .true.)  ! update TI based on the EMPTY simulation
    call do_phaseshifting()

    ! initialize budgets
    call budg_tavg%init(AD_Inputfile, adsim)               !<-- Budget class initialization
    call budg_tavg_empty%init(Empty_Inputfile, emptysim)   !<-- Budget class initialization
    call budg_vavg%init(HIT_Inputfile, hit)                !<-- Budget class initialization
    call budg_tavg_deficit%init(budg_tavg_empty, inputfile, budg_tavg)

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic()
    do while (adsim%tsim < adsim%tstop)
        dt1 = adsim%get_dt(recompute=.true.)
        dt2 = emptysim%get_dt(recompute=.true.)
        if (freeze_HIT) then
            dt = min(dt1, dt2)  ! don't consider frozen_igrid dt constraints in time stepping
        else
            dt3 = hit%get_dt(recompute=.true.)
            dt = min(dt1, dt2, dt3)
        endif

        if (synchronize_RK_substeps) then
            adsim%dt = dt
            emptysim%dt = dt
            hit%dt = dt
            ! Stage 1
            call adsim%advance_SSP_RK45_Stage_1()
            call emptysim%advance_SSP_RK45_Stage_1()
            call hit%advance_SSP_RK45_Stage_1()
            ! Stage 2
            call adsim%advance_SSP_RK45_Stage_2()
            call emptysim%advance_SSP_RK45_Stage_2()
            call hit%advance_SSP_RK45_Stage_2()
            ! Stage 3
            call adsim%advance_SSP_RK45_Stage_3()
            call emptysim%advance_SSP_RK45_Stage_3()
            call hit%advance_SSP_RK45_Stage_3()
            ! Stage 4
            call adsim%advance_SSP_RK45_Stage_4()
            call emptysim%advance_SSP_RK45_Stage_4()
            call hit%advance_SSP_RK45_Stage_4()
            ! Stage 5
            call adsim%advance_SSP_RK45_Stage_5()
            call emptysim%advance_SSP_RK45_Stage_5()
            call hit%advance_SSP_RK45_Stage_5()
            ! Call wrap up
            call adsim%wrapup_timestep()
            call emptysim%wrapup_timestep()
            call hit%wrapup_timestep()

        else
            call adsim%timeAdvance(dt)
            call emptysim%timeAdvance(dt)
            call hit%timeAdvance(dt)
        end if

        call budg_tavg%doBudgets()       !<--- perform budget related operations
        call budg_vavg%doBudgets()       !<--- perform budget related operations
        call budg_tavg_empty%doBudgets()       !<--- perform budget related operations
        call budg_tavg_deficit%doBudgets()     !<--- perform budget related operations

        ! phaseshift turbulent fringe targets using the laminar fringe targets
        if (control_TI) call update_TI_fact(emptysim, .false.)
        call do_phaseshifting()

        call doTemporalStuff(adsim, 1)
        call doTemporalStuff(emptysim, 0)
        if (.not. freeze_HIT) call doTemporalStuff(hit, 2)
    end do

    ! wrapup tasks
    call budg_tavg%doBudgets(.true.)   !<--- force dump if budget calculation had started
    call budg_vavg%doBudgets(.true.)   !<--- force dump if budget calculation had started
    call budg_tavg_empty%doBudgets(.true.)   !<--- force dump if budget calculation had started
    call budg_tavg_deficit%doBudgets(.true.) !<--- force dump if budget calculation had started

    call budg_tavg%destroy()           !<-- release memory taken by the budget class
    call budg_vavg%destroy()           !<-- release memory taken by the budget class
    call budg_tavg_empty%destroy()           !<-- release memory taken by the budget class
    call budg_tavg_deficit%destroy()         !<-- release memory taken by the budget class

    if (applyfilters) then
        do fid = 1,nfilters
            call filt(fid)%destroy()
        end do
        if (allocated(pid)) deallocate(pid)
        deallocate(filt)
    end if

    call hit%finalize_io()
    call adsim%finalize_io()

    call hit%destroy()
    call adsim%destroy()

    deallocate(hit, adsim)

    ! deallocate fringe targets
    deallocate(utarget0, vtarget0, wtarget0)
    deallocate(utarget, vtarget, wtarget)
    deallocate(utarget_1d, vtarget_1d)
    deallocate(z_global)
    if (adsim%isStratified) deallocate(Ttarget0)

    call MPI_Finalize(ierr)

contains

! Do phase shifting here - program variables are still in scope
    subroutine do_phaseshifting()
        real(rkind), dimension(size(z_global,3)) :: x_shift_z, y_shift_z
        real(rkind) :: x_shift
        integer :: ad_st, hit_st, hit_en

        if (advect_shear) then
            ! need to take the full z-domain
            x_shift_z = adsim%tsim * utarget_1d(1, 1,:)
            y_shift_z = adsim%tsim * vtarget_1d(1, 1,:)

            ! if sheared, then advect the HIT flow with different freestream velocity as a function of z
            call hit%spectC%bandpassFilter_and_phaseshift_z(hit%whatC, hit%rbuffxC(:,:,:,1), x_shift_z, y_shift_z)
            call hit%interpolate_cellField_to_edgeField(hit%rbuffxC(:,:,:,1), hit%rbuffxE(:,:,:,1),0,0)
            call hit%spectC%bandpassFilter_and_phaseshift_z(hit%uhat, hit%rbuffxC(:,:,:,1), x_shift_z, y_shift_z)
            call hit%spectC%bandpassFilter_and_phaseshift_z(hit%vhat, hit%rbuffxC(:,:,:,2), x_shift_z, y_shift_z)
        else
            ! Set the true target field for AD simulation
            x_shift = adsim%tsim * InflowSpeed

            call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC, hit%rbuffxC(:,:,:,1), x_shift)
            call hit%interpolate_cellField_to_edgeField(hit%rbuffxC(:,:,:,1), hit%rbuffxE(:,:,:,1),0,0)
            call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat, hit%rbuffxC(:,:,:,1), x_shift)
            call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat, hit%rbuffxC(:,:,:,2), x_shift)

        end if
        ! Now modify rhs HIT field appropriately
        ad_st = nxADSim - nxfringe / aniso_x + 1
        hit_st = nxhitsim - nxfringe + 1
        hit_en = nxhitsim
        utarget(ad_st:nxADSim,:,:) = hit%rbuffxC(hit_st:hit_en:aniso_x,:,:,1)*TI_fact + utarget0(ad_st:nxADSim,:,:)
        vtarget(ad_st:nxADSim,:,:) = hit%rbuffxC(hit_st:hit_en:aniso_x,:,:,2)*TI_fact + vtarget0(ad_st:nxADSim,:,:)
        wtarget(ad_st:nxADSim,:,:) = hit%rbuffxE(hit_st:hit_en:aniso_x,:,:,1)*TI_fact
    end subroutine

    ! Update TI gain
    subroutine update_TI_fact(sim, first_timestep)
        use IncompressibleGrid, only : igrid
        use constants, only          : zero, one, two, three
        use reductions, only         : p_sum
        use exits, only              : message
        use HIT_shear_parameters, only: zmid_for_TI

        class(igrid), allocatable, target, intent(in) :: sim
        real(rkind), dimension(sim%gpC%xsz(2), sim%gpC%xsz(3)) :: buff1
        real(rkind) :: TI_inst, tke_avg, error
        integer :: nz_norm = 1
        logical, intent(in) :: first_timestep

        if (TI_target < 0) then
            return  ! Doesn't compute/update anything
        end if

        ! need to compute TKE, TI
        buff1 = 0.5 * ((sim%u(TI_xid,:,:)-utarget0(TI_xid,:,:))**2 + (sim%v(TI_xid,:,:)-vtarget0(TI_xid,:,:))**2 + (sim%wC(TI_xid,:,:))**2)  ! TKE
        if (TI_at_rotor) then
            where (abs(sim%mesh(TI_xid,:,:,3) - zmid_for_TI) > 0.5)
                buff1 = zero  ! mask points outside of the rotor region
            end where
            nz_norm = 1 / sim%dz  ! number of points inside masked region
        else
            nz_norm = sim%nz
        end if
        tke_avg = p_sum(buff1) / (sim%ny*nz_norm)
        TI_inst = sqrt(two / three * tke_avg) / InflowSpeed

        if (first_timestep) then
            ! try to start with a reasonable guess for the gain variable
            !!! BETTER SOLUTION FOR RESTARTS: JUST WRITE TI_FACT TO A FILE !!!
            if (TI_inst .ge. 1e-6) then  
                ! If TI_inst is not machine zero, then this is probably from restart files
                TI_fact = sqrt(three / two / hit%getMeanKE()) * TI_inst
            else
                ! If TI_inst is basically zero, then set the "guess" for TI_fact based on TI_target
                TI_fact = sqrt(three / two / hit%getMeanKE()) * TI_target
            end if
        else if (time_stop_TIcont > 0 .and. sim%tsim > time_stop_TIcont) then
            continue  !!! do not update the controller anymore !!! (but still print debug messages)
        else
            error = TI_target - TI_inst
            integral_err = integral_err + (error * sim.dt / KIinv_TI)
            TI_fact = max(zero, Kp_TI * error + integral_err)
        end if

        if (debug_TI_gain) then
            call message(1, "update_TI: TI_inst", TI_inst)
            call message(1, "update_TI: TI_fact", TI_fact)
        end if
    end subroutine

end program
