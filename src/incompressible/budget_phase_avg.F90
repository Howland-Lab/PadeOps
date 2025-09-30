module budgets_phase_avg_mod
    ! general imports used within function
    use kind_parameters, only: rkind, clen
    use decomp_2d
    use incompressibleGrid, only: igrid  
    use exits, only: GracefulExit
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
        integer :: iphase ! convert phase [0, 1) to integer [0, 100)
        real(rkind) :: tol
    contains
        procedure :: phase_avg_init
        procedure :: doBudgets
        procedure :: dump_budget_field
    
    end type budgets_phase_avg

contains

    subroutine phase_avg_init(this, inputfile, igrid_sim, phase, tol, cfg)
        class(budgets_phase_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim
        ! phase and tolerance passed in from multi-phase budget namelist
        real(rkind), intent(in)  :: phase
        real(rkind), intent(in)  :: tol
        type(time_budget_config), intent(in), optional :: cfg

        ! ensure phase is in correct range
        if((phase < 0.0_rkind) .or. (phase > 1.0_rkind)) then
            call GracefulExit("Phase must be between 0-1.", 100)
        endif
        ! first initialize time-average values
        call this%budgets_time_avg%init(inputfile, igrid_sim, cfg)
        ! add phase and tol values
        if (abs(phase - 1.0_rkind) < tol) then ! if phase is approximately 1 -> set to 0 
            this%phase = 0.0_rkind
        else ! set phase to provided phase between [0, 1)
            this%phase = phase
        endif
        this%iphase = nint(this%phase * 100.0_rkind)
        this%tol = tol
    end subroutine phase_avg_init

    subroutine doBudgets(this, forceDump)
        class(budgets_phase_avg), intent(inout) :: this
        logical, intent(in), optional :: forceDump
        real(rkind) :: delx, uturb, surge_amp, surge_freq
        real(rkind) :: sim_curr_phase
        ! get needed arguments from first wind turbine (assumes all turbines move the same)
        delx = this%igrid_sim%WindTurbineArr%dynamicArray(1)%delx
        uturb = this%igrid_sim%WindTurbineArr%dynamicArray(1)%ut
        surge_amp = this%igrid_sim%WindTurbineArr%dynamicArray(1)%surge_amplitude
        surge_freq = this%igrid_sim%WindTurbineArr%dynamicArray(1)%surge_freq
        ! TODO: right now this only uses surge!!! update to allow pitch later
        ! only do the budget if timestep is correct phase of turbine motion

        ! compute the phase of the first turbine
        sim_curr_phase = compute_phase(delx, uturb, surge_amp, surge_freq, this%tol)
        ! if current phase, then call time-average doBudgets on current budget
        if (abs(sim_curr_phase - this%phase) < this%tol) then
            call this%budgets_time_avg%doBudgets(forceDump)
        end if
    end subroutine doBudgets

    subroutine dump_budget_field(this, field, fieldID, BudgetID)
        use decomp_2d_io
        class(budgets_phase_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(in) :: field
        integer, intent(in) :: fieldID, BudgetID
        character(len=clen) :: fname, tempname 
        character(len=clen) :: budgets_dir
        ! write budget to file that includes the phase
        budgets_dir = this%get_budgets_dir()
        write(tempname,"(A3,I2.2,A7,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A6,I3.3,A4)") "Run", this%get_run_id(), "_budget", BudgetID, &
            "_term", fieldID, "_t", this%igrid_sim%step,"_n", this%get_counter(),"_phase", this%iphase,".s3D"
        fname = budgets_dir(:len_trim(budgets_dir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,field,fname, this%igrid_sim%gpC)
    end subroutine dump_budget_field

    ! TODO: right now this only uses surge!!! update to allow pitch later
    pure function compute_phase(x, U, A, f, tol) result(phase)
        !! Compute normalized phase (0–1).
        real(rkind), intent(in) :: x, U, A, f, tol
        real(rkind), parameter :: pi = acos(-1.0_rkind)
        real(rkind) :: phi, sincomp, coscomp
        real(rkind) :: phase

        ! Construct normalized sine/cosine inputs
        sincomp = (2.0_rkind * pi * f / A) * x
        coscomp = U / A

        ! Compute raw phase angle in radians
        phi = atan2(sincomp, coscomp)

        ! Normalize to [0,1)
        phase = modulo(phi / (2.0_rkind * pi), 1.0_rkind)

        if (abs(phase - 1.0_rkind) < tol) then
            phase = 0.0_rkind
        end if
    end function compute_phase

end module budgets_phase_avg_mod