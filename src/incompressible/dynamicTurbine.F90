module dynamicTurbineMod
    use kind_parameters, only : rkind, clen
    use constants, only : zero, one, two, pi
    use exits, only : GracefulExit, message
    use mpi
!    use incompressibleGrid, only: igrid  ! for now, try to do this without the igrid object
    use actuatorDisk_filteredMod, only : actuatorDisk_filtered
    
    implicit none
    private
    public :: dynamicTurbine
    
    ! default module variables
    integer :: iounit

    type :: dynamicTurbine
        ! Additional degrees of freedom for turbine
        ! Modifies the forcing function by rebuilding the kernel and adjusting turbine velocity ud
        ! KSH 04/20/2024
        
        real(rkind) :: xloc, yloc, zloc  ! unperturbed turbine coordinates
        real(rkind) :: delx = zero, dely = zero, delz = zero  ! turbine position perturbation
        real(rkind) :: ut, vt, wt  ! turbine velocity direction components
        real(rkind) :: yaw, tilt, roll = zero  ! turbine angles 
        real(rkind) :: time  ! simulation time, non-dimensional
        real(rkind) :: surge_freq, surge_amplitude, pitch_amplitude
        real(rkind) :: phase_turbine = zero
        real(rkind) :: static_tilt ! mean tilt value -> tilt can vary sinusoidally around

        ! methods to implement motion: 
        logical :: use_dynamic_turbine, use_simple_periodic, use_timeseries, verbose

        ! Time series storage
        integer :: n_timeseries = 0
        real(rkind), allocatable :: ts_time(:), ts_surge(:), ts_pitch(:), ts_uturb(:)
        integer :: ts_index = 1  ! current index in time series
        character(len=clen) :: timeseries_file

        logical :: do_redraw = .false.  ! redraw turbine this timestep? 
        
        type(actuatorDisk_filtered), pointer :: turbine  ! TODO: make a generic basicTurbine class
        
    contains
        procedure :: init
        procedure :: destroy
        procedure :: time_advance 
        procedure :: sinusoid_update
        procedure :: timeseries_update
        procedure :: read_timeseries
        
    end type

contains

subroutine init(this, turbine)
    class(dynamicTurbine), intent(inout)   :: this
    class(actuatorDisk_filtered), intent(in), target :: turbine  ! TODO make generic turbine
    
    logical :: use_dynamic_turbine = .true., use_simple_periodic = .true., use_timeseries = .false.
    logical :: verbose = .false.
    character(len=clen) :: fname, timeseries_file
    real(rkind) :: surge_freq = zero, surge_amplitude = zero, pitch_amplitude = zero

    ! read namelist
    namelist /DYNAMICTURBINE/ use_simple_periodic, use_timeseries, timeseries_file, &
                              surge_freq, surge_amplitude, pitch_amplitude, &
                              verbose

    ioUnit = 55
    call turbine%get_fname(fname)  ! get inputfile name from turbine
    open(unit=ioUnit, file=trim(fname), form='FORMATTED', action="read")
    read(unit=ioUnit, NML=DYNAMICTURBINE)
    close(ioUnit)

    this%turbine => turbine
    this%verbose = verbose  ! add additional print statements
    this%time = zero  ! TODO - may need to pass a non-zero start-time in 
    call this%turbine%get_pos(this%xloc, this%yloc, this%zloc)
    call this%turbine%get_angle(this%yaw, this%tilt)  ! stored in DEGREES
    this%static_tilt = this%tilt

    ! save namelist variables
    this%use_dynamic_turbine = use_dynamic_turbine
    this%use_simple_periodic = use_simple_periodic  ! simple periodic motion given by sinusoid_update
    this%use_timeseries = use_timeseries ! motion defined by time-series
    this%timeseries_file = timeseries_file ! time-series read in if use_timeseries is true
    this%surge_freq = surge_freq            ! surge frequency, non-dimensionalized
    this%surge_amplitude = surge_amplitude  ! surge amplitude =  u_d,max/U
    this%pitch_amplitude = pitch_amplitude  ! pitch amplitude, in degrees

    ! Validate flags
    if (this%use_simple_periodic .and. this%use_timeseries) then
        call gracefulExit("Cannot specify both use_simple_periodic and use_timeseries", 424)
    endif

    ! Read time series if specified
    if (this%use_timeseries) then
        call this%read_timeseries()
    endif

    call message(1, 'Initialized dynamicTurbine module')

end subroutine

subroutine destroy(this)
    class(dynamicTurbine), intent(inout) :: this
    ! nothing to deallocate at the moment
end subroutine

! do time advancement step - consider passing time instead of dt into this function !
subroutine time_advance(this, dt)
    class(dynamicTurbine), intent(inout) :: this
    real(rkind), intent(in) :: dt
    real(rkind) :: tmp

    ! STEP 1: Update time
    this%time = this%time + dt

    ! STEP 2: first, update the position & velocity of the turbine (if not needed, skip time_advance)
    if (this%use_simple_periodic) then
        call this%sinusoid_update(dt)
    else if (this%use_timeseries) then
        call this%timeseries_update()
    else
        call gracefulExit("Unknown or missing time advance type in DYNAMICTURBINE module", 423)
    endif
    
    ! STEP 3: redraw the turbine forcing kernel
    if (this%do_redraw) then
        call this%turbine%set_pos(this%xloc + this%delx, this%yloc + this%dely, this%zloc + this%delz)
        call this%turbine%set_angle(this%yaw, this%tilt)

        ! now, redraw the turbine forcing kernel
        call this%turbine%redraw()

        ! and set do_redraw to false for the start of the next timestep
        this%do_redraw = .false.
    endif

    ! STEP 4: update turbine velocity
    call this%turbine%set_ut(this%ut, this%vt, this%wt)

    if (this%verbose) then
        tmp = this%turbine%get_udisk()
        call message(0, 'dynamicTurbine: time_advance called at t', this%time)
        call message(1, 'dynamicTurbine: position delta x', this%delx)
        call message(1, 'dynamicTurbine: velocity uturb', this%ut)
        ! call message(1, 'dynamicTurbine: velocity udisk', tmp)  ! this lags one time step
        if (this%pitch_amplitude > zero) then
            call message(1, 'dynamicTurbine: turbine tilt (deg.)', this%tilt)
        endif
        if (this%use_simple_periodic) then
            call message(1, 'dynamicTurbine: normalized turbine phase', this%phase_turbine)
        endif
    endif

end subroutine

! most basic case: sinusoidal variation
subroutine sinusoid_update(this, dt)
    class(dynamicTurbine), intent(inout) :: this
    real(rkind), intent(in) :: dt
    real(rkind) :: omega, omega_t

    if (.not. (this%surge_freq == zero)) then
        ! update the turbine phase
        this%phase_turbine = this%phase_turbine + this%surge_freq * dt
        this%phase_turbine = modulo(this%phase_turbine, one)

        ! update omega values
        omega = two * pi * this%surge_freq
        omega_t = two * pi * this%phase_turbine

        ! sinusoid needs updating every timestep as long as f != 0, A != 0
        this%do_redraw = .true.

        ! update surge velocity and position
        this%ut = this%surge_amplitude * cos(omega_t)
        this%delx = this%surge_amplitude / omega * sin(omega_t)
        ! update the pitch (tilt) as well
        this%tilt = this%pitch_amplitude * sin(omega_t) + this%static_tilt
    endif

end subroutine

subroutine read_timeseries(this)
    class(dynamicTurbine), intent(inout) :: this
    integer :: unit_ts, ios, n
    real(rkind) :: t, surge, pitch, uturb
    
    unit_ts = 56
    n = 0
    
    ! First pass: count lines
    open(unit=unit_ts, file=trim(this%timeseries_file), form='FORMATTED', action="read", iostat=ios)
    if (ios /= 0) then
        call gracefulExit("Failed to open timeseries file: " // trim(this%timeseries_file), 425)
    endif
    
    do
        read(unit_ts, *, iostat=ios) t, surge, pitch, uturb
        if (ios /= 0) exit
        n = n + 1
    enddo
    
    this%n_timeseries = n
    if (this%n_timeseries < 2) call gracefulExit("Timeseries_file needs at least 2 rows", 430)

    allocate(this%ts_time(n))
    allocate(this%ts_surge(n))
    allocate(this%ts_pitch(n))
    allocate(this%ts_uturb(n))
    
    ! Second pass: read data
    rewind(unit_ts)
    do n = 1, this%n_timeseries
        read(unit_ts, *, iostat=ios) t, surge, pitch, uturb
        if (ios /= 0) exit
        this%ts_time(n) = t
        this%ts_surge(n) = surge
        this%ts_pitch(n) = pitch
        this%ts_uturb(n) = uturb
    enddo
    close(unit_ts)

    do n = 2, this%n_timeseries
        if (this%ts_time(n) <= this%ts_time(n-1)) then
            call gracefulExit("Timeseries_file time must be strictly increasing", 431)
        end if
    end do
    
    if (this%verbose) then
        call message(1, 'Done reading timeseries!')
    endif
    ! initialize time series index
    this%ts_index = 1
end subroutine

subroutine timeseries_update(this)
    class(dynamicTurbine), intent(inout) :: this
    integer :: i
    real(rkind) :: t0, t1, a
    real(rkind), parameter :: eps = 1.0e-12_rkind

    ! Clamp outside range
    if (this%time <= this%ts_time(1)) then
        this%ts_index = 1
        this%delx = this%ts_surge(1)
        this%tilt = this%ts_pitch(1) + this%static_tilt
        this%ut   = zero
        this%do_redraw = .true.
        return
    else if (this%time >= this%ts_time(this%n_timeseries)) then
        this%ts_index = this%n_timeseries - 1
        this%delx = this%ts_surge(this%n_timeseries)
        this%tilt = this%ts_pitch(this%n_timeseries) + this%static_tilt
        this%ut   = zero
        this%do_redraw = .true.
        return
    end if

    ! Move index forward (time is monotone increasing)
    do while (this%ts_index < this%n_timeseries - 1 .and. &
              this%time > this%ts_time(this%ts_index + 1))
        this%ts_index = this%ts_index + 1
    end do

    i  = this%ts_index
    t0 = this%ts_time(i)
    t1 = this%ts_time(i+1)

    ! Exact hit on upper node: advance index and snap
    if (abs(this%time - t1) <= eps .and. i < this%n_timeseries - 1) then
        this%ts_index = i + 1
        i = this%ts_index
        this%delx = this%ts_surge(i)
        this%tilt = this%ts_pitch(i) + this%static_tilt
        this%ut   = this%ts_uturb(i)
    else
        a = (this%time - t0) / (t1 - t0)
        this%delx = this%ts_surge(i) + a * (this%ts_surge(i+1) - this%ts_surge(i))
        this%tilt = this%ts_pitch(i) + a * (this%ts_pitch(i+1) - this%ts_pitch(i)) + this%static_tilt
        this%ut   = this%ts_uturb(i) + a * (this%ts_uturb(i+1) - this%ts_uturb(i))
    end if

    this%do_redraw = .true.
end subroutine

end module
