module budgets_phase_avg_mod
    ! general imports used within function
    use kind_parameters, only: rkind
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
        procedure :: init
        procedure :: doBudgets
        procedure :: dump_budget_field
    
    end type budgets_phase_avg

contains

    subroutine init(this, inputfile, igrid_sim, phase, tol, cfg)
        class(budgets_phase_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim
        ! phase and tolerance passed in from multi-phase budget namelist
        real(rkind), intent(in)  :: phase
        real(rkind), intent(in)  :: tol
        type(time_budget_config), intent(in), optional :: cfg

        ! first initialize time-average values
        call this%budgets_time_avg%init(inputfile, igrid_sim, cfg)
        ! add phase values
        this%phase = phase
        this%tol = tol
        ! check if phase is between 0-1
        if((this%phase < 0.0_rkind) .or. (this%phase > 1.0_rkind)) then
            call GracefulExit("Phase must be between 0-1.", 100)
        endif
    end subroutine init

    subroutine doBudgets(this, forceDump)
        class(budgets_phase_avg), intent(inout) :: this
        logical, intent(in), optional :: forceDump
        ! only do the budget if timestep is correct phase of turbine motion
        turb1 = this%igrid_sim%WindTurbineArr%dynamicArray(1)
        ! TODO: right now this only uses surge!!! update to allow pitch later
        sim_curr_phase = compute_phase(turb1%delx, turb1%ut, turb1%surge_amplitude, turb1%surge_freq)
        if (abs(sim_curr_phase - this%phase) < tol)
            call this%budgets_time_avg%doBudgets(this, forceDump)
        end if
    end subroutine doBudgets

    subroutine dump_budget_field(this, field, fieldID, BudgetID)
        use decomp_2d_io
        class(budgets_phase_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(in) :: field
        integer, intent(in) :: fieldID, BudgetID
        character(len=clen) :: fname, tempname 
        ! Convert phase (0–1) to integer 0–100
        iphase = nint(phase * 100.0_rkind)
        write(tempname,"(A3,I2.2,A7,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A6,I3.3,A4)") "Run",this%run_id,"_budget",BudgetID,"_term",fieldID,"_t",this%igrid_sim%step,"_n",this%counter,"_phase", iphase,".s3D"
        fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)

        call decomp_2d_write_one(1,field,fname, this%igrid_sim%gpC)
    end subroutine dump_budget_field

    ! TODO: right now this only uses surge!!! update to allow pitch later
    pure function compute_phase(x, u, A, f, tol) result(phase)
        !! Compute normalized phase (0–1).
        real(rkind), intent(in) :: x, y, A, f, tol
        real(rkind), parameter :: pi = acos(-1.0_rkind)
        real(rkind) :: phi, sincomp, coscomp
        real(rkind) :: phase

        ! Construct normalized sine/cosine inputs
        sincomp = (2.0_rkind * pi * f / A) * x
        coscomp = u / A

        ! Compute raw phase angle in radians
        phi = atan2(sincomp, coscomp)

        ! Normalize to [0,1)
        phase = modulo(phi / (2.0_rkind * pi), 1.0_rkind)

        if (abs(phase - 1.0_rkind) < tol) phase = 0.0_rkind end if
    end function compute_phase

end module budgets_phase_avg_mod