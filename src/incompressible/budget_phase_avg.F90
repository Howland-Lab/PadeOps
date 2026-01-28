module budgets_phase_avg_mod
    ! general imports used within function
    use kind_parameters, only: rkind, clen
    use decomp_2d
    use incompressibleGrid, only: igrid  
    use exits, only: message, GracefulExit
    ! import time-average type to act as parent to phase-average type
    use budgets_time_avg_mod

    implicit none 
    private
    public :: budgets_phase_avg

    ! create phase-average type as child to time-average type
    type, extends(budgets_time_avg) :: budgets_phase_avg
        private
        ! phase fields needed for phase-averaged budgets
        ! note that this CANNOT be read in via a namelist -> use budget_multi_phase_avg!
        real(rkind) :: phase
        real(rkind) :: tol
    contains
        procedure :: phase_init
        procedure :: phase_doBudgets
    
    end type budgets_phase_avg

contains

    subroutine phase_init(this, inputfile, igrid_sim, phase, tol, cfg)
        class(budgets_phase_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim
        ! phase and tolerance passed in from multi-phase budget namelist
        real(rkind), intent(in)  :: phase
        real(rkind), intent(in)  :: tol
        type(time_budget_config), intent(in), optional :: cfg
        integer :: iphase
        character(len=clen) :: phase_str 

        ! ensure phase is in correct range
        if((phase < 0.0_rkind) .or. (phase > 1.0_rkind)) then
            call GracefulExit("Phase must be between 0-1.", 100)
        endif

        ! first initialize time-average values
        call this%init(inputfile, igrid_sim, cfg)
        this%phase = phase
        this%tol = tol

        ! create file name suffix for phase
        iphase = nint(this%phase * 100.0_rkind)
        write(phase_str,"(A6,I3.3)") "_phase", iphase
        this%file_suffix = phase_str
        call message(0, "Phase budget initialized: file_suffix = " // trim(this%file_suffix))
    end subroutine phase_init

    subroutine phase_doBudgets(this, sim_curr_phase, forceDump)
        class(budgets_phase_avg), intent(inout) :: this
        real(rkind), intent(in) :: sim_curr_phase
        logical, intent(in), optional :: forceDump
        logical :: startedBudget, runBudget
        ! check if we need to start (or have already started) budget calculations
        startedBudget = .FALSE.
        if ((this%tidx_budget_start>0) .and. (this%igrid_sim%step>this%tidx_budget_start)) then
            startedBudget = .TRUE.
        endif
        if ((this%time_budget_start>0) .and. (this%igrid_sim%tsim>this%time_budget_start)) then
            startedBudget = .TRUE.
        endif
        if (startedBudget) then
            ! check if we need to calculate budget due to phase
            runBudget = check_runBudget(sim_curr_phase, this%phase, this%tol)
            ! update force dump
            call this%updateForceDump(forceDump)
            ! if current phase, then call time-average doBudgets on current budget
            if (runBudget) then
                call message(0, "Updating phase budget for phase", this%phase)
                call this%doBudgets(this%forceDump) ! calls time average doBudgets
            else if (this%forceDump) then
                call this%dumpBudget()
                call message(0,"Dumped a phase budget .stt file")
                this%forceDump = .FALSE.
            end if
        end if
    end subroutine phase_doBudgets

    pure function check_runBudget(sim_curr_phase, budget_phase, tol) result(runBudget)
        real(rkind), intent(in) :: sim_curr_phase, budget_phase, tol
        logical :: runBudget
        real(rkind) :: wrapped_tol

        runBudget = .FALSE.
        ! if sim_curr_phase isn't within tol of 0 or 1, simple check
        if (abs(sim_curr_phase - budget_phase) < tol) then
            runBudget = .TRUE.
        end if
        ! if sim_curr_phase is within tol of 1
        if (budget_phase + tol > 1.0_rkind) then
            wrapped_tol = budget_phase + tol - 1.0_rkind
            if (sim_curr_phase < wrapped_tol) then
                runBudget = .TRUE.
            end if
        end if
        ! if sim_curr_phase is within tol of 0
        if (budget_phase - tol < 0.0_rkind) then
            wrapped_tol = 1.0_rkind + (budget_phase - tol)
            if (sim_curr_phase > wrapped_tol) then
                runBudget = .TRUE.
            end if
        end if
    end function check_runBudget

end module budgets_phase_avg_mod