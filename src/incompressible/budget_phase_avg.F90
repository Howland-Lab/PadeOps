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

        ! add phase and tol values
        if (abs(phase - 1.0_rkind) < tol) then ! if phase is approximately 1 -> set to 0 
            this%phase = 0.0_rkind
        else ! set phase to provided phase between [0, 1)
            this%phase = phase
        endif
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
        ! update force dump
        call this%updateForceDump(forceDump)
        ! if current phase, then call time-average doBudgets on current budget
        if (abs(sim_curr_phase - this%phase) < this%tol) then
            call message(0, "Updating phase budget for phase", this%phase)
            call this%doBudgets(forceDump) ! calls time average doBudgets
        else if (this%forceDump) then
            call this%dumpBudget()
            call message(0,"Dumped a budget .stt file")
            this%forceDump = .FALSE.
        end if
    end subroutine phase_doBudgets

    ! subroutine phase_dump_budget_field(this, field, fieldID, BudgetID)
    !     use decomp_2d_io
    !     class(budgets_phase_avg), intent(inout) :: this
    !     real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(in) :: field
    !     integer, intent(in) :: fieldID, BudgetID
    !     character(len=clen) :: fname, tempname 
    !     character(len=clen) :: budgets_dir
    !     ! write budget to file that includes the phase
    !     budgets_dir = this%get_budgets_dir()
    !     write(tempname,"(A3,I2.2,A7,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A6,I3.3,A4)") "Run", this%get_run_id(), "_budget", BudgetID, &
    !         "_term", fieldID, "_t", this%igrid_sim%step,"_n", this%get_counter(),"_phase", this%iphase,".s3D"
    !     fname = budgets_dir(:len_trim(budgets_dir))//"/"//trim(tempname)
    !     call decomp_2d_write_one(1,field,fname, this%igrid_sim%gpC)
    !     call message(0,"Dumped a phase budget file")
    ! end subroutine phase_dump_budget_field

end module budgets_phase_avg_mod