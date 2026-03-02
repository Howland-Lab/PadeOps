! Concurrent-precursor problem for inhomogeneous Dirichlet
! boundary conditions in the stable PBL with wind turbines.

#include "stable_pbl_concurrent_files/initialize.F90"
#include "stable_pbl_concurrent_files/temporalHook.F90"

program stable_pbl_concurrent
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use budgets_time_avg_mod, only: budgets_time_avg
    use budgets_time_avg_deficit_mod, only: budgets_time_avg_deficit
    use exits, only: message, gracefulExit

    implicit none

    type(igrid), allocatable, target :: primary, precursor
    character(len=clen) :: inputfile, primary_inputfile, precursor_inputfile
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg, pre_budg_tavg
    type(budgets_time_avg_deficit) :: budg_tavg_deficit
    real(rkind) :: dt1, dt2, dt
    logical :: synchronize_RK_fringe = .true., do_deficit_budgets = .false.

    namelist /concurrent/ primary_inputfile, precursor_inputfile, synchronize_RK_fringe, do_deficit_budgets

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(precursor)               !<-- Allocate precursor
    allocate(primary)                 !<-- Allocate primary
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)

    call compute_xdim_udim(primary_inputFile)    !<-- Reads the \PHYSICS\ namelist to compute (xdim, udim) from (Ro, Fr)
    
    ! INITIALIZE PRIMARY SIMULATION
    call primary%init(primary_inputFile, .true.)
    call primary%start_io(.false.)     ! do not dump IO fields on init (avoid overwriting turbine data)
    call primary%printDivergence()

    ! INITIALIZE PRECURSOR SIMULATION
    call precursor%init(precursor_inputFile, .false.)
    precursor%Am_I_Primary = .false.
    call precursor%start_io(.true.)

    if (primary%usefringe) then
        call primary%fringe_x%associateFringeTargets(precursor%u, precursor%v, precursor%w, precursor%T)
        call primary%fringe_x%associateFringeTarget_scalar(precursor%T)
    end if

    if (primary%useControl .AND. primary%dummy_controller)then
        if(.NOT. precursor%useControl)then
            call gracefulExit("Primary has a dummy controller, but precursor does not have a controller at all. Exiting.", 44)
        elseif(precursor%dummy_controller) then
            call gracefulExit("Both Primary and Precursor have dummy controllers. Exiting.", 44)
        else
            if(.NOT. allocated(precursor%angCont_yaw))then
                call gracefulExit("Precursor does not have an active controller, and Primary has a dummy controller. Exiting.", 44)
            end if
            ! Point to the precursor's controller
            primary%angCont_yaw_dummy => precursor%angCont_yaw
            call message(0, 'Dummy controller points to main controller.')
        end if
    end if

    call budg_tavg%init(primary_inputfile, primary)             !<-- Budget class initialization
    call pre_budg_tavg%init(precursor_inputfile, precursor)     !<-- Budget class initialization
    if (do_deficit_budgets) then                                !<-- Budget class initialization for the deficit
        call budg_tavg_deficit%init(pre_budg_tavg, primary_inputfile, budg_tavg)
    end if

    if (primary%useWindTurbines) then
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
        if (do_deficit_budgets) call budg_tavg_deficit%doBudgets()

        call doTemporalStuff(primary,  1)
        call doTemporalStuff(precursor,2)

    end do

    call budg_tavg%destroy()                !<-- release memory taken by the budget classes
    call pre_budg_tavg%destroy()
    if (do_deficit_budgets) call budg_tavg_deficit%destroy()

    call precursor%finalize_io()
    call primary%finalize_io()

    call precursor%destroy()
    call primary%destroy()

    deallocate(precursor, primary)

    call MPI_Finalize(ierr)

end program