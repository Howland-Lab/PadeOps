module budgets_multi_phase_avg_mod
    ! general imports used within function
    use kind_parameters, only: rkind, clen
    use incompressibleGrid, only: igrid  
    use exits, only: GracefulExit, message
    use constants, only: zero,
    ! import time-average type to act as parent to phase-average type
    use budgets_phase_avg_mod
    use budgets_time_avg_mod, only: time_budget_config

    implicit none 
    private
    public :: budgets_multi_phase_avg
    ! need to declare namelist variables at a module-level due to reading in a list
    ! this causes memory weirdness that can otherwise cause a seg-fault
    logical :: do_budgets
    integer :: nphases
    real(rkind) :: tol
    integer, parameter :: max_phases = 16 ! maximum number of phases a user could request
    real(rkind) :: phases(max_phases)
    ! need BUDGET_MULTI_PHASE_AVG to get information on multi-phase budgets
    namelist /BUDGET_MULTI_PHASE_AVG/ do_budgets, nphases, tol, phases

    type :: budgets_multi_phase_avg
        type(igrid), pointer, public :: igrid_sim => null()
        logical :: do_budgets
        real(rkind), allocatable :: phases(:)
        type(budgets_phase_avg), allocatable :: phase_budgets(:)
        real(rkind) :: tol
        integer :: nphases
        logical :: has_moving_turbine
    contains
        procedure :: init
        procedure :: doBudgets
        procedure :: destroy

    end type budgets_multi_phase_avg

contains

    subroutine init(this, inputfile, igrid_sim, amplit_inflow, freq_inflow)
        class(budgets_multi_phase_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim
        real(rkind), intent(in) optional :: amplit_inflow, freq_inflow ! if inflow is surging, these must be provided
        ! values from namelist
        integer :: ioUnit, ierr
        integer:: i, j
        type(time_budget_config) :: cfg
        character(len=clen) :: overlap_msg
        logical :: has_moving_turbine
        logical :: has_surging_inflow

        ! save pointer to igrid_sim
        this%igrid_sim => igrid_sim 

        ! read in namelist variables (declared above in module)
        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_MULTI_PHASE_AVG)
        read(unit=ioUnit, NML=AD_CoriolisINPUT)
        close(ioUnit)

        ! check if there is a dynamic turbine
        has_moving_turbine = .FALSE.
        if (igrid_sim%useWindTurbines) then
            has_moving_turbine = igrid_sim%WindTurbineArr%useDynamicTurbine
        end
        ! check if the inflow is surging
        has_surging_inflow = .FALSE.
        if (present(amplit_inflow) .and. present(freq_inflow)) then
            has_surging_inflow = (amplit_inflow /= zero) .and. (freq_inflow /= zero)
        endif
        if (.not. (has_moving_turbine .or. has_surging_inflow)) then
            call GracefulExit("Turbine isn't dynamic and inflow isn't surging. Nothing to else to perform phase budget on.", 100)
        endif
        if (has_moving_turbine .and. has_surging_inflow) then
            call GracefulExit("Turbine is dynamic and inflow is surging. Unclear what to perform phase budget on.", 100)
        endif
        ! save has_moving_turbine to budgets_multi_phase_avg object
        this%has_moving_turbine = has_moving_turbine

        ! save namelist variables to budgets_multi_phase_avg object
        this%do_budgets = do_budgets
        this%nphases = nphases
        this%phases = phases
        this%tol = tol
        if (this%do_budgets) then
            ! check to ensure users have provided phases
            if (this%nphases == 0) then
                call GracefulExit("Phase-averaged budgets turned on, but phases array is empty!", 101)
            endif
            ! check to ensure that given phases don't overlap
            do i = 1, this%nphases
                do j = 1, this%nphases
                    if ((.not. (i .eq. j)) .and. (phase_overlaps(this%phases(i), this%phases(j), this%tol))) then
                        write(overlap_msg, '(A, F10.4, A, F10.4)') &
                            "The following phases overlap with given tolerance: ", &
                            this%phases(i), " and ", this%phases(j)
                        call message(0, overlap_msg)
                    end if
                end do
            end do
            ! get default time budget config values and update from the namelist
            cfg = time_budget_config()
            call cfg%update_budget_config_from_namelist(inputfile, this%do_budgets)
            ! create one phase-average budget per requested phase (children of time-average budgets)
            allocate(this%phase_budgets(this%nphases))
            do i = 1, this%nphases
                call this%phase_budgets(i)%phase_init(inputfile, igrid_sim, this%phases(i), tol, cfg)
            end do
        end if
    end subroutine init

    ! check if two phases overlap with given tolerance
    pure function phase_overlaps(p1, p2, tol) result(overlaps)
        real(rkind), intent(in) :: p1, p2, tol
        real(rkind) :: min_phase, max_phase
        logical :: overlaps
        ! find min and max extent of p1, adjusting for wrapping around at 0 and 1
        min_phase = p1 - tol
        max_phase = p1 + tol
        if (min_phase < 0.0_rkind) min_phase = min_phase + 1.0_rkind
        if (max_phase > 1.0_rkind) max_phase = max_phase - 1.0_rkind
        ! check if phases overlap (accounting for wrap around if needed in the first branch)
        overlaps = .FALSE.
        if (min_phase > max_phase) then
             if (.not. ((p2 > max_phase) .and. (p2 < min_phase))) then
                overlaps = .TRUE.
             end if
        else if ((p2 > min_phase) .and. (p2 < max_phase)) then
            overlaps = .TRUE.
        endif
    end function phase_overlaps

    ! TODO: right now this only uses surge!!! update to allow pitch later
    pure function compute_phase(x, U, A, f) result(phase)
        !! Compute normalized phase (0–1).
        real(rkind), intent(in) :: x, U, A, f
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
    end function compute_phase ! TODO: remove compute_phase and instead use phase_turbine or phase_inflow

    subroutine doBudgets(this, forceDump, phase_inflow)
        class(budgets_multi_phase_avg), intent(inout) :: this
        logical, intent(in), optional :: forceDump
        real(rkind), intent(in), optional :: phase_inflow
        real(rkind) :: delx, uturb, surge_amp, surge_freq
        real(rkind) :: sim_curr_phase
        integer:: i
        ! call doBudgets for each phase budget
        if (this%do_budgets)  then
            ! get needed arguments from first wind turbine (assumes all turbines move the same)
            delx = this%igrid_sim%WindTurbineArr%dynamicArray(1)%delx
            uturb = this%igrid_sim%WindTurbineArr%dynamicArray(1)%ut
            surge_amp = this%igrid_sim%WindTurbineArr%dynamicArray(1)%surge_amplitude
            surge_freq = this%igrid_sim%WindTurbineArr%dynamicArray(1)%surge_freq
            ! TODO: right now this only uses surge!!! update to allow pitch later
            ! only do the budget if timestep is correct phase of turbine motion
            ! compute the phase of the first turbine
            sim_curr_phase = compute_phase(delx, uturb, surge_amp, surge_freq)
            call message(0, "Current phase ", sim_curr_phase)
            do i = 1, this%nphases
                call this%phase_budgets(i)%phase_doBudgets(sim_curr_phase, forceDump)
            end do
        end if 
    end subroutine doBudgets

    subroutine destroy(this)
        class(budgets_multi_phase_avg), intent(inout) :: this
        integer :: i
        nullify(this%igrid_sim)
        ! destroy each phase budgets
        if (this%do_budgets) then 
            do i = 1, this%nphases
                call this%phase_budgets(i)%destroy()
            end do
        end if
    end subroutine destroy

end module budgets_multi_phase_avg_mod