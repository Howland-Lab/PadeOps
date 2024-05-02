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

        ! methods to implement motion: 
        logical :: use_dynamic_turbine, use_simple_periodic, verbose 

        logical :: do_redraw = .false.  ! redraw turbine this timestep? 
        
        type(actuatorDisk_filtered), pointer :: turbine  ! TODO: make a generic basicTurbine class
        
    contains
        procedure :: init
        procedure :: destroy
        procedure :: time_advance 
        procedure :: sinusoid_update
        
    end type

contains

subroutine init(this, turbine)
    class(dynamicTurbine), intent(inout)   :: this
    class(actuatorDisk_filtered), intent(in), target :: turbine  ! TODO make generic turbine
    
    logical :: use_dynamic_turbine = .true., use_simple_periodic = .true., verbose = .false.
    character(len=clen) :: fname
    real(rkind) :: surge_freq = zero, surge_amplitude = zero, pitch_amplitude = zero

    ! read namelist
    namelist /DYNAMICTURBINE/ use_simple_periodic, surge_freq, surge_amplitude, verbose, pitch_amplitude

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

    ! save namelist variables
    this%use_dynamic_turbine = use_dynamic_turbine
    this%use_simple_periodic = use_simple_periodic  ! simple periodic motion given by sinusoid_update
    this%surge_freq = surge_freq            ! surge frequency, non-dimensionalized
    this%surge_amplitude = surge_amplitude  ! surge amplitude =  u_d,max/U
    this%pitch_amplitude = pitch_amplitude  ! pitch amplitude, in degrees

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
        call this%sinusoid_update()
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
    endif

end subroutine

! most basic case: sinusoidal variation
subroutine sinusoid_update(this)
    class(dynamicTurbine), intent(inout) :: this
    real(rkind) :: omega

    omega = two * pi * this%surge_freq
    if (.not. (this%surge_freq == zero)) then
        ! sinusoid needs updating every timestep as long as f != 0, A != 0
        this%do_redraw = .true.

        ! TODO update to include other DOF
        this%ut = this%surge_amplitude * cos(omega * this%time)
        this%delx = this%surge_amplitude / omega * sin(omega * this%time)
        ! update the pitch (tilt) as well
        this%tilt = this%pitch_amplitude * sin(omega * this%time)
        
    endif

end subroutine

end module
