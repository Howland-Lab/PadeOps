module actuatorDisk_FilteredMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc
    use Gridtools, only: linspace

    implicit none

    private
    public :: actuatorDisk_filtered
     
    type :: actuatorDisk_filtered
        ! Implementation of Shapiro et al. (2019) Filtered ADM
        ! Kirby Heck 07/20/2022

        ! Actuator Disk Info
        integer :: xLoc_idx, ActutorDisk_T2ID, tInd = 1
        real(rkind) :: yaw, tilt, ut, powerBaseline, hubDirection
        real(rkind) :: xLoc, yLoc, zLoc, dx, dy, dz, dV
        real(rkind) :: diam, cT, pfactor, normfactor, OneBydelSq, Cp, thick, npts
        real(rkind) :: uface = zero, vface = zero, wface = zero  ! LES velocity, disk-averaged
        real(rkind) :: uturb, vturb, wturb  ! turbine motion vector

        real(rkind), dimension(:), allocatable :: xs, ys, zs  ! list of UNYAWED points for the ADM
        real(rkind), dimension(:), allocatable :: powerTime, uTime, vTime
        logical :: useDynamicYaw, quickDecomp

        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta, M  ! Shapiro smearing size, corr. factor M<1
        real(rkind), dimension(:), allocatable :: xline, yline, zline
        real(rkind), dimension(:,:,:), pointer :: xG, yG, zG
        
        ! Pointers to memory buffers 
        real(rkind), dimension(:,:,:), allocatable :: rbuff, blanks, speed, scalarSource

        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: color, myComm, myComm_nproc, myComm_nrank

        character(len=clen) :: fname  ! save inputfile path

    contains
        procedure :: init
        procedure :: destroy
        procedure :: get_RHS
        procedure :: get_R1
        procedure :: get_R2
        procedure :: get_R
        procedure :: get_weights
        procedure :: get_power

        ! new functions -- move to basic_turbine? 
        procedure :: get_pos
        procedure :: get_angle
        procedure :: get_fname
        procedure :: get_ut
        procedure :: get_udisk  ! velocity that the turbine sees
        procedure :: set_pos
        procedure :: set_angle
        procedure :: set_ut
        procedure :: redraw
    end type


contains

subroutine init(this, inputDir, ActuatorDisk_ID, xG, yG, zG)
    class(actuatorDisk_filtered), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDisk_ID
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit, ierr
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0
    real(rkind) :: diam=0.08d0, cT=0.65d0, yaw=0.d0, tilt=0.d0, h  !, Cp = 0.3
    real(rkind) :: thickness=1.5d0, filterWidth=0.5, time2initialize=0.d0
    logical :: useCorrection=.true., useDynamicYaw=.false., quickDecomp=.false., use_h=.false.

    ! Read input file for this turbine    
    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, cT, yaw, tilt, filterWidth, useCorrection, &
                             useDynamicYaw, thickness, quickDecomp, use_h
    write(tempname,"(A13,I4.4,A10)") "ActuatorDisk_", ActuatorDisk_ID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    this%fname = fname

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED', action="read")
    read(unit=ioUnit, NML=ACTUATOR_DISK)
    close(ioUnit)
    
    call message(0, "Initializing Actuator Disk (ADM Type=5) number", ActuatorDisk_ID)
    call tic()
    
    ! link grids and read inputs 
    this%dx=xG(2,1,1)-xG(1,1,1)
    this%dy=yG(1,2,1)-yG(1,1,1)
    this%dz=zG(1,1,2)-zG(1,1,1)
    this%dV = this%dx*this%dy*this%dz
    this%xLoc = xLoc; this%yLoc = yLoc; this%zLoc = zLoc
    this%cT = cT; this%diam = diam; this%yaw = yaw; this%tilt = tilt
    this%ut = 1.d0!; this%Cp = Cp

    this%uturb = zero; this%vturb = zero; this%wturb = zero
    
    this%nxLoc = size(xG,1); this%nyLoc = size(xG,2); this%nzLoc = size(xG,3)
    
    ! Allocate stuff
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    
    this%xG => xG; this%yG => yG; this%zG => zG
    this%xLine = xG(:,1,1); this%yLine = yG(1,:,1); this%zLine = zG(1,1,:)

    ! allocate memory buffers
    allocate(this%rbuff(this%nxLoc, this%nyLoc, this%nzLoc))
    allocate(this%blanks(this%nxLoc, this%nyLoc, this%nzLoc))
    allocate(this%speed(this%nxLoc, this%nyLoc, this%nzLoc))
    allocate(this%scalarsource(this%nxLoc, this%nyLoc, this%nzLoc))
    
    ! copied from ADM T2
    this%Am_I_Split = .TRUE. ! TODO: Fix me later by flagging where the turbine is
    if (this%Am_I_Split) then
        call MPI_COMM_SPLIT(mpi_comm_world, this%color, nrank, this%mycomm, ierr)
        call MPI_COMM_RANK(this%mycomm, this%myComm_nrank, ierr) 
        call MPI_COMM_SIZE(this%mycomm, this%myComm_nproc, ierr)
    end if 
    
    ! this ensures that only ONE turbine is keeping track of power and writing to disk
    if((this%Am_I_Split .and. this%myComm_nrank==0) .or. (.not. this%Am_I_Split)) then
!        write(*,*) "Only one allocated? YES/NO"
        allocate(this%powerTime(1000000))  
        allocate(this%uTime(1000000))  
        allocate(this%vTime(1000000))
    end if

    ! Set thickness
    this%thick = thickness*this%dx
    if (use_h) then
        ! use h to dimensionalize the filterwidth
        h = sqrt((this%xLine(2) - this%xLine(1))**2 + (this%yLine(2) - this%yLine(1))**2 + (this%zLine(2) - this%zLine(1))**2)
        this%delta = filterWidth * h
    else
        ! use the turbine diameter to dimensionalize the filterwidth
        this%delta = filterWidth * this%diam 
    endif
    ! Thickness is only used if quickDecomp = .TRUE.
    this%quickDecomp = quickDecomp
    if (quickDecomp) then
        call message(1, "ADM: using quick decomposition in x")
    else
        call message(1, "ADM: using full kernel integration")
    end if

    ! Get (unrotated) turbine location points
    call sample_on_circle(this%diam, this%yLoc, this%zLoc, this%ys, this%zs, this%dy, this%dz)
    this%npts = size(this%ys,1)
    call message(1, "NUMBER OF POINTS: ", this%npts)
    allocate(this%xs(size(this%ys)))
    this%xs = this%xLoc
     
    ! Correction factor (taylor series approx.)
    if (useCorrection) then
        this%M = one / (one + this%cT/two*this%delta/sqrt(3.d0*pi)/this%diam)
    else
        this%M = one
    end if 
    
    this%useDynamicYaw = useDynamicYaw
    if (this%useDynamicYaw) then
        call message(2, "Using Dynamic Yaw")
    else
        call message(2, "Using static turbine.")
        call this%redraw()  ! get_weights(this) 
    end if
    
    call message(2, "Smearing grid parameter, Delta", this%delta)
    call message(2, "Smearing grid parameter, Thickness", this%thick)
    call message(2, "Turbine location:")
    call message(3, "x = ", this%xLoc)
    call message(3, "y = ", this%yLoc)
    call message(3, "z = ", this%zLoc)
    call toc(mpi_comm_world, time2initialize)
    call message(2, "Time (seconds) to initialize", time2initialize)
end subroutine 

subroutine destroy(this)
    class(actuatordisk_filtered), intent(inout) :: this

    deallocate(this%rbuff, this%blanks, this%speed, this%scalarSource) 
    nullify(this%xG, this%yG, this%zG)
end subroutine 

! Convolution in x (streamwise) direction
subroutine get_R1(this, R1) 
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nxLoc), intent(out) :: R1
    real(rkind), dimension(this%nxLoc) :: xLine
    real(rkind) :: tmp

    xLine = this%xLine - this%xLoc 
    tmp = sqrt(6.d0)/this%delta
    R1 = erf(tmp*(xLine + this%thick/two)) - &
         erf(tmp*(xLine - this%thick/two))
    R1 = R1 / (two*this%thick)
end subroutine

! Convolution in the yz-plane (build numerically with Green's function)
subroutine get_R2(this, ys, zs, R2)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nyLoc, this%nzLoc), intent(out) :: R2
    real(rkind), dimension(this%nyLoc, this%nzLoc) :: yy, zz
    real(rkind), dimension(:), intent(in), allocatable :: ys, zs
    real(rkind) :: tmp
    real(rkind), dimension(this%nyLoc, this%nzLoc) :: stamp
    integer :: j
    
    R2 = zero
    stamp = zero
    yy = this%yG(1, :, :)
    zz = this%zG(1, :, :)
    tmp = -6.d0 / (this%delta**2)
    ! we need to compute the integral, iterate through points on disk face: 
    do j = 1, size(ys)
        stamp = exp(tmp*((ys(j)-yy)**2 + (zs(j)-zz)**2))
        R2 = R2 + stamp
    end do

    tmp = this%dy*this%dz*24.d0 / (pi**2 * this%diam**2 * this%delta**2)
    R2 = R2*tmp  ! weight accordingly to the sampling density in sample_on_circle()
end subroutine

! Use this generic Greens function for yawed turbines
subroutine get_R(this)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind) :: yrad, trad, xs, ys, zs, C1, xtmp, ytmp, ztmp  ! rotations, in radians
    real(rkind), dimension(this%npts) :: xi, yi, zi
    integer :: k
    
    ! First, rotate all the points with the yaw and tilt
    ! call message(1, "Building kernel for turbine yaw:", this%yaw)
    yrad = this%yaw*pi/180.d0
    trad = this%tilt*pi/180.d0
    do k = 1, this%npts
        xs = this%xs(k); ys = this%ys(k); zs = this%zs(k)
        ! apply yaw rotation, +z = positive yaw (e.g., Howland, et al. 2022)
        xtmp = (xs-this%xLoc)*cos(yrad) - (ys-this%yLoc)*sin(yrad) + this%xLoc
        ytmp = (xs-this%xLoc)*sin(yrad) + (ys-this%yLoc)*cos(yrad) + this%yLoc
        ztmp = zs
        
        ! then apply tilt rotation, +y = positive tilt (e.g., Bossuyt, et al. 2021)
        xi(k) = (xtmp-this%xLoc)*cos(trad) + (ztmp-this%zLoc)*sin(trad) + this%xLoc
        yi(k) = ytmp
        zi(k) = -(xtmp-this%xLoc)*sin(trad) + (ztmp-this%zLoc)*cos(trad) + this%zLoc
    end do
    
    ! now xi, yi, zi are the rotated coordinates, assemble w/Greens function 
    ! this may take a while... 
    ! TODO: can speed this up if only a subsection of the domain is used
    C1 = (6.d0/pi/this%delta**2)**(three/two)
    ! TODO: May need to zero scalarsource for dynamic yaw
    do k = 1, this%npts
        this%rbuff = (this%xG-xi(k))**2 + (this%yG-yi(k))**2 + (this%zG-zi(k))**2
        this%scalarsource = this%scalarsource + C1*exp(-6.d0*this%rbuff/this%delta**2) 
    end do
    
    ! scalarsource NOT necessarily normalized to integrate to 1 (yet), do this in get_weights()
end subroutine

subroutine get_weights(this)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nyLoc, this%nzLoc) :: R2
    real(rkind), dimension(this%nxLoc) :: R1
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc) :: R
        
    if ((abs(this%yaw) < 1e-3) .and. (abs(this%tilt) < 1e-3)) then
        if (this%quickDecomp) then
            !aligned with the x-direction, use the "quick" kernel creation
            call this%get_R2(this%ys, this%zs,R2)
            call this%get_R1(R1)
            
            ! Not sure why, but setting the product of R1*R2 directly leads to
            ! memory errors: (2023/07/17) 
            !this%scalarsource = spread(spread(R1,2,this%nyLoc),3,this%nzLoc) * spread(R2, 1, this%nxLoc)
            
            ! Instead, store in intermediate variable R: 
            R = spread(R2, 1, this%nxLoc) * spread(spread(R1, 2, this%nyLoc), 3,this%nzLoc)
            this%scalarsource = R
        else
            call this%get_R()  ! bypass quick kernel
        end if
    else
        ! build the kernel from 3D gaussian kernel
        call this%get_R()
    end if
     
    ! minimum threshold tolerance
    where (this%scalarsource < 1.d-10)
        this%scalarsource = 0
    end where

    ! normalize so R integrates to 1 exactly
    this%scalarsource = this%scalarsource / (p_sum(this%scalarsource)*this%dV) 
end subroutine

! sample a circle with points spaced dx, dy apart and centered at xcen, ycen
subroutine sample_on_circle(diam, xcen, ycen, xloc, yloc, dx, dy)
    use gridtools, only: linspace
    real(rkind), intent(in) :: diam, xcen, ycen, dx, dy
    real(rkind) :: R
    integer, dimension(:), allocatable :: tag
    real(rkind), dimension(:), allocatable :: xline, yline
    real(rkind), dimension(:), allocatable, intent(out) :: xloc, yloc
    real(rkind), dimension(:), allocatable :: xtmp, ytmp, rtmp
    integer :: idx, i, j, nsz, iidx, nx_per_R, ny_per_R, nx, ny, np
    
    R = diam/two
    nx_per_R = ceiling(R/dx); ny_per_R = ceiling(R/dy)
    nx = nx_per_R*2 + 1
    ny = ny_per_R*2 + 1
    np = nx*ny  ! total number of points

    allocate(xline(nx), yline(ny))
    allocate(xtmp(np),ytmp(np), rtmp(np), tag(np))
    
    ! initialize linearly-spaced arrays 
    ! this is necessary to do independently of the grid xG, yG, zG 
    ! because parallelization splits the grid up
    xline = (/(i, i=-nx_per_R, nx_per_R)/) * dx
    yline = (/(i, i=-ny_per_R, ny_per_R)/) * dy
    
    ! reshapes xline, yline: 
!    xtmp = reshape(spread(xline, 1, ny), [np])
!    ytmp = reshape(spread(yline, 2, nx), [np])  ! why doesn't reshape() work? 
    idx = 1
    do j = 1,ny
        do i = 1,nx
            xtmp(idx) = xline(i); ytmp(idx) = yline(j)
            idx = idx + 1
        end do 
    end do
    rtmp = sqrt(xtmp**2 + ytmp**2) 
    tag = 0
    where (rtmp < R) 
        tag = 1
    end where
    nsz = sum(tag) ! number of points on disk face
    
    allocate(xloc(nsz), yloc(nsz))
    iidx = 1
    do idx = 1,size(tag)
        if (tag(idx) == 1) then
            xloc(iidx) = xtmp(idx)
            yloc(iidx) = ytmp(idx)
            iidx = iidx + 1
        end if
    end do

    xloc = xloc + xcen; yloc = yloc + ycen 
    deallocate(xtmp, ytmp, rtmp, tag)  ! deallocate temporary variables
end subroutine

! Right hand side forcing term for the ADM
subroutine get_RHS(this, u, v, w, rhsxvals, rhsyvals, rhszvals)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind) :: yaw, tilt
    real(rkind) :: usp_sq, force, vface
    real(rkind), dimension(3,1) :: n=[1,0,0], tau=[0,1,0] !xn, Ft
    real(rkind), dimension(3,3) :: R, T

    ! update yaw and tilt of the turbine
    if (.not. this%useDynamicYaw .and. (this%yaw - yaw*180.d0/pi)>1.d-8) then
        call GracefulExit("Turbine prescribed yaw changed, but useDynamicYaw is OFF", 423)
    end if

    yaw = this%yaw * pi/180.d0
    tilt = this%tilt * pi/180.d0

    n = reshape([1,0,0], shape(n))  ! reset the normal vector
    tau = reshape([0, 1, 0], shape(tau))  ! also reset the tangent vector

    R = reshape([cos(yaw), -sin(yaw), 0.d0, &
                 sin(yaw), cos(yaw), 0.d0, &
                 0.d0, 0.d0, 1.d0], shape(R), order=[2, 1])
    T = reshape([cos(tilt), zero, sin(tilt), &
                 zero, one, zero, &
                 -sin(tilt), zero, cos(tilt)], shape(T), order=[2, 1])
    n = matmul(T, matmul(R, n))
    tau = matmul(T, matmul(R, tau))

    ! OLD method of computing disk velocity
    ! this%speed = u*n(1,1) + v*n(2,1) + w*n(3,1)
    ! this%ut = this%M*p_sum(this%scalarSource*this%speed)*this%dV
    ! vface = p_sum(this%scalarSource*(u*tau(1,1) + v*tau(2,1) + w*tau(3,1)))*this%dV

    ! NEW method -- requires more p_sum but results in a vector
    this%uface = p_sum(this%scalarSource * u) * this%dV
    this%vface = p_sum(this%scalarSource * v) * this%dV
    this%wface = p_sum(this%scalarSource * w) * this%dV
    this%ut = this%M * ((this%uface - this%uturb) * n(1,1) + (this%vface - this%vturb) * n(2,1) + (this%wface - this%wturb) * n(3,1))
    vface = ((this%uface - this%uturb) * tau(1,1) + (this%vface - this%vturb) * tau(2,1) + (this%wface - this%wturb) * tau(3,1))
    
    ! call message(1, 'DEBUG ActuatorDisk: this%ut', this%ut)
    ! TODO: May need to update yaw before calling get_weights()
    ! if (this%useDynamicYaw) then
    !     call this%get_weights() 
    ! end if

    ! Mean speed at the turbine, corrected with factor M
    usp_sq = (this%ut)**2
    force = -0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq

    rhsxvals = rhsxvals + force * n(1,1) * this%scalarSource
    rhsyvals = rhsyvals + force * n(2,1) * this%scalarSource
    rhszvals = rhszvals + force * n(3,1) * this%scalarSource

    if (allocated(this%powerTime)) then   ! check allocated so only one processor writes data
    !        if((this%Am_I_Split .and. this%myComm_nrank==0) .or. (.not. this%Am_I_Split)) then
        if (usp_sq /= 0.d0) then
            this%powerTime(this%tInd) = this%get_power()
            this%uTime(this%tInd) = this%ut
            this%vTime(this%tInd) = vface
            this%tInd = this%tInd + 1
        end if
    end if

end subroutine

! TODO - MOVE THESE FUNCTIONS TO A BASE CLASS

! Get power: 
function get_power(this) result(power)
    class(actuatordisk_filtered), intent(in) :: this
    real(rkind) :: power
    ! OLD: 
    !power = 0.5d0*this%Cp*(pi*(this%diam**2)/4.d0)*(this%ut)**3
    
    ! New power assumes induction theory and uses CT' = CP'
    power = 0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*this%ut**3
end function

! access turbine position
subroutine get_pos(this, x, y, z)
    class(actuatordisk_filtered), intent(in) :: this
    real(rkind), intent(out) :: x, y, z
    
    x = this%xloc; y = this%yloc; z = this%zloc
end subroutine

! accessor for yaw and tilt angles
subroutine get_angle(this, yaw, tilt)
    class(actuatorDisk_filtered), intent(in) :: this
    real(rkind), intent(out) :: yaw, tilt
    
    ! returns angles in DEGREES
    yaw = this%yaw
    tilt = this%tilt
end subroutine

! accessor for inputfile name
subroutine get_fname(this, fname)
    class(actuatordisk_filtered), intent(in) :: this
    character(len=clen), intent(out) :: fname
    
    fname = this%fname
end subroutine 

! modifier for turbine position
subroutine set_pos(this, x, y, z)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), intent(in) :: x, y, z
    
    this%xloc = x; this%yloc = y; this%zloc = z
end subroutine

! modifier for angles
subroutine set_angle(this, yaw, tilt)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), intent(in) :: yaw, tilt
    
    this%yaw = yaw; this%tilt = tilt
end subroutine

! modifier for turbine velocity
subroutine set_ut(this, ut, vt, wt)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), intent(in) :: ut, vt, wt
    
    this%uturb = ut; this%vturb = vt; this%wturb = wt
end subroutine

! accessor for turbine velocity
subroutine get_ut(this, ut, vt, wt)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), intent(out) :: ut, vt, wt
    
    ut = this%uturb; vt = this%vturb; wt = this%wturb
end subroutine

! accessor for disk velocity
function get_udisk(this) result(udisk)
    class(actuatordisk_filtered), intent(in) :: this
    real(rkind) :: udisk
    udisk = this%ut    
end function

! rebuilds forcing kernel
subroutine redraw(this)
    class(actuatordisk_filtered), intent(inout) :: this

    ! (re)sample points, this is quick
    call sample_on_circle(this%diam, this%yloc, this%zloc, this%ys, this%zs, this%dy, this%dz)
    this%npts = size(this%ys, 1)  
    this%xs = this%xloc
    
    ! (re)compute weights
    call this%get_weights()
end subroutine

end module 
