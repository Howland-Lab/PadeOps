! Template for PadeOps

#include "HIT_shear_files/initialize.F90"
#include "HIT_shear_files/temporalHook.F90"

program HIT_shear
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

    implicit none

    ! type(igrid), allocatable, target :: adsim
    class(igrid), allocatable, target :: hit, adsim  ! make these polymorphic so we can freeze the turbulence if 
    character(len=clen) :: inputfile, HIT_InputFile, AD_InputFile, fof_dir, filoutdir
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg
    type(budgets_vol_avg)  :: budg_vavg
    real(rkind), dimension(:,:,:), allocatable :: utarget, vtarget, wtarget
    real(rkind) :: dt1 = one, dt2 = one, dt = one
    real(rkind) :: k_bandpass_left = 10.d0, k_bandpass_right = 64.d0, TI_xloc = 0
    type(fof), dimension(:), allocatable :: filt
    integer, dimension(:), allocatable :: pid
    integer :: fid, nfilters = 2, tid_FIL_FullField = 75, tid_FIL_Planes = 4, TI_xid
    logical :: applyFilters = .false., control_TI = .true., freeze_HIT = .false.
    logical, parameter :: synchronize_RK_substeps = .true.

    namelist /concurrent/ HIT_InputFile, AD_InputFile, InflowSpeed, k_bandpass_left, k_bandpass_right, TI, TI_xloc, TI_fact, freeze_HIT
    namelist /FILTER_INFO/ applyfilters, nfilters, fof_dir, tid_FIL_FullField, tid_FIL_Planes, filoutdir

    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    ! read concurrent input file
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    read(unit=ioUnit, NML=FILTER_INFO)
    close(ioUnit)

    allocate(adsim)
    if (freeze_HIT) then
        allocate(frozen_igrid :: hit)
    else
        allocate(igrid :: hit)
    end if

    ! initialize igrid objects
    simulationID = 1
    call adsim%init(AD_InputFile, .true.)
    call adsim%start_io(.true.)
    call adsim%printDivergence()

    call mpi_barrier(mpi_comm_world, ierr)
    call message("Initialized PRIMARY simulation")

    simulationID = 2
    call hit%init(HIT_InputFile, .false.)
    hit%Am_I_Primary = .false.
    call hit%start_io(.true.)
    call hit%printDivergence()
    call message("Initialized CONCURRENT HIT simulation")
    if (freeze_HIT) then
        call message(1, "HIT targets are FROZEN")
    end if

    call make_global_zaxis(adsim)  ! allocate the global-z axis variables

    ! decide whether to turn on the TI controller
    if ((TI < 0) .and. (TI_fact < 0))then
        TI_fact = one
    end if

    if (TI_fact >= 0) then  ! TI_fact specified, no TI controller
        control_TI = .false.
        call message(0, "TI controller not used")
        call message(1, "Using fixed TI gain/loss: ", TI_fact)
    else  ! Use TI control
        control_TI = .true.
        TI_fact = one
        TI_xid = minloc(abs(adsim%mesh(:,1,1,1) - TI_xloc), 1)  ! xid corresponding to TI sampling location
        call message(0, "TI controller activated, tracking x-location:", adsim%mesh(TI_xid,1,1,1))
    end if

    ! allocate target cells for the fringe
    allocate(utarget0(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(vtarget0(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(wtarget0(adsim%gpE%xsz(1), adsim%gpE%xsz(2), adsim%gpE%xsz(3)))
    call init_fringe_targets(AD_inputfile, adsim%mesh)  ! populates utarget0, vtarget0, wtarget0

    ! allocate moving (turbulent) targets
    allocate(utarget(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(vtarget(adsim%gpC%xsz(1), adsim%gpC%xsz(2), adsim%gpC%xsz(3)))
    allocate(wtarget(adsim%gpE%xsz(1), adsim%gpE%xsz(2), adsim%gpE%xsz(3)))

    ! initialize turbulent fluctuations as zero
    utarget = zero
    vtarget = zero
    wtarget = zero

    ! initialize bandpass filter
    call hit%spectC%init_bandpass_filter(k_bandpass_left, k_bandpass_right, hit%cbuffzC(:,:,:,1), hit%cbuffyC(:,:,:,1))

    ! now initialize turbulent fringe targets
    if (adsim%usedoublefringex) then
        call message(0, "Setting double fringe targets")
        ! first fringe is re-laminarization
        call adsim%fringe_x1%associateFringeTargets(utarget0, vtarget0, wtarget0)

        ! second fringe is turbulent
        call adsim%fringe_x2%associateFringeTargets(utarget, vtarget, wtarget)
    else
        call message(0, "Setting fringe targets")
        ! first (only) fringe is turbulent
        call adsim%fringe_x%associateFringeTargets(utarget, vtarget, wtarget)
    end if

    ! phaseshift turbulent fringe targets using the laminar fringe targets
    call update_TI_fact(adsim, TI_xid)
    call do_phaseshifting(hit, adsim, utarget, vtarget, wtarget)

    ! initialize budgets
    call budg_tavg%init(AD_Inputfile, adsim)   !<-- Budget class initialization
    call budg_vavg%init(HIT_Inputfile, hit)    !<-- Budget class initialization

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic()
    do while (adsim%tsim < adsim%tstop)
        dt1 = adsim%get_dt(recompute=.true.)
        dt2 = hit%get_dt(recompute=.true.)
        dt = min(dt1, dt2)

        if (synchronize_RK_substeps) then
            adsim%dt = dt
            hit%dt = dt
            ! Stage 1
            call adsim%advance_SSP_RK45_Stage_1()
            call hit%advance_SSP_RK45_Stage_1()
            ! Stage 2
            call adsim%advance_SSP_RK45_Stage_2()
            call hit%advance_SSP_RK45_Stage_2()
            ! Stage 3
            call adsim%advance_SSP_RK45_Stage_3()
            call hit%advance_SSP_RK45_Stage_3()
            ! Stage 4
            call adsim%advance_SSP_RK45_Stage_4()
            call hit%advance_SSP_RK45_Stage_4()
            ! Stage 5
            call adsim%advance_SSP_RK45_Stage_5()
            call hit%advance_SSP_RK45_Stage_5()
            ! Call wrap up
            call adsim%wrapup_timestep()
            call hit%wrapup_timestep()

        else
            call adsim%timeAdvance(dt)
            call hit%timeAdvance(dt)
        end if

        call budg_tavg%doBudgets()       !<--- perform budget related operations
        call budg_vavg%doBudgets()       !<--- perform budget related operations

        ! phaseshift turbulent fringe targets using the laminar fringe targets
        call update_TI_fact(adsim, TI_xid)
        call do_phaseshifting(hit, adsim, utarget, vtarget, wtarget)

        call doTemporalStuff(adsim, 1)
        call doTemporalStuff(hit  , 2)
    end do

    ! wrapup tasks
    call budg_tavg%doBudgets(.true.)   !<--- force dump if budget calculation had started
    call budg_vavg%doBudgets(.true.)   !<--- force dump if budget calculation had started

    call budg_tavg%destroy()           !<-- release memory taken by the budget class
    call budg_vavg%destroy()           !<-- release memory taken by the budget class

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

    call MPI_Finalize(ierr)

end program
