! Template for PadeOps
#include "neutral_pbl_concurrent_files/initialize.F90"
#include "neutral_pbl_concurrent_files/temporalHook.F90"

program neutral_pbl_concurrent
    use mpi
    use decomp_2d
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message, GracefulExit
    use constants, only: one, zero
    use budgets_time_avg_mod, only: budgets_time_avg
    use budgets_time_avg_deficit_mod, only: budgets_time_avg_deficit
    use link_turbine_to_scalar, only: setup_turb_scalar_source

    implicit none

    type(igrid), allocatable, target :: primary, precursor
    character(len=clen) :: inputfile, primary_inputFile, precursor_inputFile
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg, pre_budg_tavg  ! added precursor budgets
    type(budgets_time_avg_deficit) :: budg_tavg_deficit  ! added deficit budgets, can turn off in inputfile (off by default)
    real(rkind) :: dt1, dt2, dt
    integer :: nxloc, nyloc, nzloc, tile_x, tile_y
    real(rkind), allocatable, dimension(:,:,:) :: utarget, vtarget, wtarget, Ttarget  ! fringe targets
    logical :: synchronize_RK_fringe = .true., do_deficit_budgets = .false., tile_precursor = .false., allocated_targets = .false.
    ! tiling variables: 
    type(decomp_info) :: gpC_upX, gpC_upXY, gpE_upX, gpE_upXY
    real(rkind), dimension(:,:,:), allocatable, target :: fxup_inX, fxup_inY, fxyup_inY
    real(rkind), dimension(:,:,:), allocatable, target :: fxup_inXE, fxup_inYE, fxyup_inYE

    namelist /concurrent/ primary_inputfile, precursor_inputfile, synchronize_RK_fringe, do_deficit_budgets, tile_precursor

    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    allocate(precursor)
    allocate(primary)
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)

    call primary%init(primary_inputFile, .true.)
    call primary%start_io(.false.)
    call primary%printDivergence()

    if (primary%useScalars) then
        call setup_turb_scalar_source(primary)
    end if

    call precursor%init(precursor_inputFile, .false.)
    precursor%Am_I_Primary = .false.
    call precursor%start_io(.true.)

    if (primary%usefringe) then  ! initialize fringe targets
        if (primary%nz .ne. precursor%nz) then
            ! needs to have same number of points in z
            call gracefulExit("usefringe is TRUE but vertical points do not match!", 125)
        endif

        if ((primary%nx .ne. precursor%nx) .or. (primary%ny .ne. precursor%ny)) then
            if (tile_precursor) then
                ! check if integer multiple in x, y, needs to be an integer multiple to tile periodic
                if ((mod(primary%nx, precursor%nx) .ne. 0) .or. (mod(primary%ny, precursor%ny) .ne. 0)) then
                    call gracefulExit("Tiling requires integer nx1/nx2 and ny1/ny2", 124)
                end if

                ! allocate target arrays
                call allocate_tile_targets()
                allocated_targets = .true.
                call message(0, "Using tiled fringe targets")
                call message(1, "x-fringe tiling: ", tile_x)
                call message(1, "y-fringe tiling: ", tile_y)
                call tile_targets()
                
                ! link pointers:
                call primary%fringe_x%associateFringeTargets(utarget, vtarget, wtarget, Ttarget)
                call primary%fringe_x%associateFringeTarget_scalar(Ttarget)
            else
                ! tile_precursor flag is not .TRUE.
                call gracefulExit("Number of points in x or y do not match between primary/precursor.", 123)
            end if

        else  ! initialize fringe targets normally
            call primary%fringe_x%associateFringeTargets(precursor%u, precursor%v, precursor%w, precursor%T)
            call primary%fringe_x%associateFringeTarget_scalar(precursor%T)
        end if
    end if

    call budg_tavg%init(primary_inputfile, primary)           !<-- Budget class initialization for the primary
    call pre_budg_tavg%init(precursor_inputfile, precursor)   !<-- Budget class initialization for the precursor
    if (do_deficit_budgets) then
        if (allocated_targets) then
            call gracefulExit("Deficit budgets do not currently support a tiled fringe", 126)
        end if
        call budg_tavg_deficit%init(pre_budg_tavg, primary_inputfile, budg_tavg)   !<-- Budget class initialization for the deficit
    end if

    if (primary%useWindTurbines) then  ! initialize wind farm control links
        call primary%WindTurbineArr%link_reference_domain_for_control(primary%u, primary%v, primary%rbuffyC, primary%rbuffzC, primary%gpC)
    end if

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic()
    do while (primary%tsim < primary%tstop)
        dt1 = primary%get_dt(recompute=.true.)
        dt2 = precursor%get_dt(recompute=.true.)
        dt = min(dt1, dt2)
        if (synchronize_RK_fringe) then
            primary%dt = dt
            precursor%dt = dt
            ! Stage 1
            call primary%advance_SSP_RK45_Stage_1()
            call precursor%advance_SSP_RK45_Stage_1()
            ! Stage 2
            call primary%advance_SSP_RK45_Stage_2()
            call precursor%advance_SSP_RK45_Stage_2()
            ! Stage 3
            call primary%advance_SSP_RK45_Stage_3()
            call precursor%advance_SSP_RK45_Stage_3()
            ! Stage 4
            call primary%advance_SSP_RK45_Stage_4()
            call precursor%advance_SSP_RK45_Stage_4()
            ! Stage 5
            call primary%advance_SSP_RK45_Stage_5()
            call precursor%advance_SSP_RK45_Stage_5()
            ! Call wrap up
            call primary%wrapup_timestep()
            call precursor%wrapup_timestep()
        else
            call primary%timeAdvance(dt)
            call precursor%timeAdvance(dt)
        end if

        call budg_tavg%doBudgets()
        call pre_budg_tavg%doBudgets()
        if (do_deficit_budgets) then
            call budg_tavg_deficit%doBudgets()
        end if

        if (allocated_targets) then
            call tile_targets()  ! re-tile fringe targets
        end if

        call doTemporalStuff(primary, 1)
        call doTemporalStuff(precursor, 2)

    end do

    call budg_tavg%destroy()           !<-- release memory taken by the budget class
    call pre_budg_tavg%destroy()

    call precursor%finalize_io()
    call primary%finalize_io()

    call precursor%destroy()
    call primary%destroy()

    deallocate(precursor, primary)

    if (allocated_targets) call deallocate_targets()

    call MPI_Finalize(ierr)

    ! end PadeOps simulation
contains

    subroutine allocate_tile_targets()
        ! allocates buffers for precursor simulation tiling
        nxloc = primary%gpC%xsz(1)
        nyloc = primary%gpC%xsz(2)
        nzloc = primary%gpC%xsz(3)  ! local grid, primary simulation: nx, ny, nz
        tile_x = primary%nx / precursor%nx
        tile_y = primary%ny / precursor%ny 

        ! decompositions for tiling intermediates: 
        call decomp_info_init(primary%nx, precursor%ny, primary%nz, gpC_upX)        !<-- up in X only
        call decomp_info_init(primary%nx, primary%ny, primary%nz, gpC_upXY)         !<-- up in X and Y only
        call decomp_info_init(primary%nx, precursor%ny, primary%nz+1, gpE_upX)      !<-- up in X only (edges)
        call decomp_info_init(primary%nx, primary%ny, primary%nz+1, gpE_upXY)       !<-- up in X and Y only (edges)

        ! buffer arrays for tiling intermediates: 
        allocate(utarget(nxloc, nyloc, nzloc), vtarget(nxloc, nyloc, nzloc), Ttarget(nxloc, nyloc, nzloc), wtarget(nxloc, nyloc, nzloc+1))
        allocate(fxup_inX(gpC_upX%xsz(1),gpC_upX%xsz(2),gpC_upX%xsz(3)))
        allocate(fxup_inY(gpC_upX%ysz(1),gpC_upX%ysz(2),gpC_upX%ysz(3)))
        allocate(fxyup_inY(gpC_upXY%ysz(1),gpC_upXY%ysz(2),gpC_upXY%ysz(3)))
        allocate(fxup_inXE(gpE_upX%xsz(1),gpE_upX%xsz(2),gpE_upX%xsz(3)))
        allocate(fxup_inYE(gpE_upX%ysz(1),gpE_upX%ysz(2),gpE_upX%ysz(3)))
        allocate(fxyup_inYE(gpE_upXY%ysz(1),gpE_upXY%ysz(2),gpE_upXY%ysz(3)))

    end subroutine

    subroutine deallocate_targets()
        deallocate(utarget, vtarget, Ttarget, wtarget)
        deallocate(fxup_inX, fxup_inY, fxyup_inY)
        deallocate(fxup_inXE, fxup_inYE, fxyup_inYE)
    end subroutine

    ! tile u, v, w, T fields
    subroutine tile_targets()
        call help_tile_targets(precursor%u, utarget)
        call help_tile_targets(precursor%v, vtarget)
        call help_tile_targets(precursor%w, wtarget, .true.)  ! use edge cells here
        call help_tile_targets(precursor%T, Ttarget)
    end subroutine

    ! helper function to do the tiling for one field
    subroutine help_tile_targets(arrIn, arrOut, use_edges)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        logical, optional, intent(in) :: use_edges
        real(rkind), dimension(:,:,:), pointer :: tilex_inX, tilex_inY, tilexy_inY
        type(decomp_info) :: gp_tileX, gp_tileXY
        integer :: i, nx, ny

        ! associate pointers
        if (present(use_edges)) then
            gp_tileX = gpE_upX
            gp_tileXY = gpE_upXY
            tilex_inX => fxup_inXE
            tilex_inY => fxup_inYE
            tilexy_inY => fxyup_inYE
        else
            gp_tileX = gpC_upX
            gp_tileXY = gpC_upXY
            tilex_inX => fxup_inX
            tilex_inY => fxup_inY
            tilexy_inY => fxyup_inY
        end if
        
        ! step 1: tile in x
        nx = size(arrIn, 1)
        do i=1,tile_x
            tilex_inX(1+(i-1)*nx:i*nx,:,:) = arrIn
        end do

        ! step 2: transpose x to y
        call transpose_x_to_y(tilex_inX, tilex_inY, gp_tileX)
        ny = size(tilex_inY, 2)

        ! step 3: tile in y
        do i=1,tile_y
            tilexy_inY(:, 1+(i-1)*ny:i*ny,:) = tilex_inY
        end do

        ! step 4: transpose y to x
        call transpose_y_to_x(tilexy_inY, arrOut, gp_tileXY)
    end subroutine

end program
