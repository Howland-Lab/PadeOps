module budgets_multi_phase_avg_mod
    ! general imports used within function
    use kind_parameters, only: rkind
    use incompressibleGrid, only: igrid  
    use exits, only: GracefulExit
    ! import time-average type to act as parent to phase-average type
    use budgets_phase_avg_mod
    use budgets_time_avg_mod, only: time_budget_config

    implicit none 
    private
    public :: budgets_multi_phase_avg

    type :: budgets_multi_phase_avg
        logical :: do_budgets
        real(rkind), allocatable :: phases(:)
        type(budgets_phase_avg), allocatable :: phase_budgets(:)
        real(rkind) :: tol
    contains
        procedure :: init
        procedure :: doBudgets ! call doBudgets for each phase budget
    end type budgets_multi_phase_avg

contains

    subroutine init(this, inputfile, igrid_sim)
        class(budgets_multi_phase_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim
        ! values from namelist
        integer :: ioUnit, ierr
        integer:: i, nphases
        logical :: do_budgets = .false. 
        real(rkind), allocatable :: phases(:)
        real(rkind) :: tol = 0.1d0
        type(time_budget_config) :: cfg

        ! ensure using a dynamic turbine or else phase averaging doesn't make much sense
        if(.not. igrid_sim%WindTurbineArr%useDynamicTurbine) then
            call GracefulExit("Turbine isn't dynamic - right now phase averaging depends on turbine postions/speed.", 100)
        endif

        ! would be good to move all of this to a phase_budget_config once things work
        namelist /BUDGET_MULTI_PHASE_AVG/ do_budgets, phases, tol
        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_MULTI_PHASE_AVG)
        close(ioUnit)

        if (.not. allocated(phases)) then
            call GracefulExit("Phases array was not read from namelist!", 101)
        endif

        nphases = size(phases)
        allocate(this%phases(nphases))
        this%phases = phases
        this%do_budgets = do_budgets
        this%tol = tol

        ! get default time budget config values and update from the namelist
        cfg = time_budget_config()
        call cfg%update_budget_config_from_namelist(inputfile)

        allocate(this%phase_budgets(nphases))
        do i = 1, nphases
            ! each time the phase -> time budget constructor is called, the namelist is read in.
            ! this can and should be fixed... perhapes all of the values from the time budget need
            ! to be in the multi-phase cosntructor list???? perhapes we can read both in here???
            call this%phase_budgets(i)%phase_avg_init(inputfile, igrid_sim, this%phases(i), tol, cfg)
        end do
    end subroutine init

    subroutine doBudgets(this, forceDump)
        class(budgets_multi_phase_avg), intent(inout) :: this
        logical, intent(in), optional :: forceDump
        integer:: i

        do i = 1, size(this%phases)
            call this%phase_budgets(i)%doBudgets(forceDump)
        end do

    end subroutine doBudgets

end module budgets_multi_phase_avg_mod