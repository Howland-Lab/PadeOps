module budgets_multi_phase_avg_mod
    ! general imports used within function
    use kind_parameters, only: rkind
    use incompressibleGrid, only: igrid  
    ! import time-average type to act as parent to phase-average type
    use budgets_phase_avg_mod

    implicit none 
    private
    public :: budgets_multi_phase_avg

    type :: budgets_multi_phase_avg
        logical :: do_budgets
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
        logical :: do_budgets = .false. 
        real(rkind), allocatable :: phases(:)
        real(rkind) :: tol = 0.1d0

        ! would be good to move all of this to a phase_budget_config once things work
        integer :: ioUnit, ierr
        namelist /BUDGET_MULTI_PHASE_AVG/ do_budgets, phases, tol
        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_MULTI_PHASE_AVG)
        close(ioUnit)

        this%do_budgets = do_budgets
        this%phases = phases
        this%tol = tol

        ! get default time budget config values and update from the namelist
        cfg = time_budget_config()
        call update_budget_config_from_namelist(cfg, inputfile)

        allocate(this%phase_budgets(size(phases)))
        do i = 1, size(phases)
            ! each time the phase -> time budget constructor is called, the namelist is read in.
            ! this can and should be fixed... perhapes all of the values from the time budget need
            ! to be in the multi-phase cosntructor list???? perhapes we can read both in here???
            rewind(unit)
            call this%phase_budgets(i)%init(inputfile, igrid_sim, phases(i), tol, cfg)
        end do
    end subroutine init

    subroutine doBudgets(this, forceDump)
        class(budgets_multi_phase_avg), intent(inout) :: this
        logical, intent(in), optional :: forceDump

        do i = 1, size(this%phases)
            call this%phase_budgets(i)%doBudgets(forceDump)
        end do

    end subroutine doBudgets

end module budgets_multi_phase_avg_mode