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
     
!    real(rkind), parameter :: alpha_Smooth = 1.d0 ! 0.9d0 ! Exonential smoothing constant
!    integer, parameter :: xReg = 8, yReg = 8, zReg = 8

    type :: actuatorDisk_filtered
        ! Implementation of Shapiro et al. (2019) Filtered ADM
        ! Kirby Heck 07/20/2022

        ! Actuator Disk Info
        integer :: xLoc_idx, ActutorDisk_T2ID, tInd = 1
        real(rkind) :: yaw, tilt, ut, powerBaseline, hubDirection
        real(rkind) :: xLoc, yLoc, zLoc, dx, dy, dz, dV
        real(rkind) :: diam, cT, pfactor, normfactor, OneBydelSq, Cp, thick
        real(rkind) :: uface = 0.d0, vface = 0.d0, wface = 0.d0
        real(rkind), dimension(:), allocatable :: xs, ys, zs
        real(rkind), dimension(:), allocatable :: powerTime, uTime, vTime
        logical :: useDynamicYaw

        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta, M  ! Shapiro smearing size, corr. factor M<1
        real(rkind), dimension(:), allocatable :: xline, yline, zline
        real(rkind), dimension(:,:,:), pointer :: xG, yG, zG
        real(rkind), dimension(:,:,:), allocatable :: smearing_base 
        
        ! Pointers to memory buffers 
        logical :: memory_buffers_linked = .false.
        real(rkind), dimension(:,:,:), pointer :: rbuff, blanks, speed, scalarSource
        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: color, myComm, myComm_nproc, myComm_nrank

    contains
        procedure :: init
        procedure :: destroy
        procedure :: get_RHS
        procedure :: get_RHS_old
        procedure :: get_R1
        procedure :: get_R2
        procedure :: get_weights
!        procedure :: get_RHS_withPower
        procedure :: link_memory_buffers
        procedure, private :: AD_force_point
        procedure :: get_power
        procedure :: dumpPower
!        procedure :: dumpPowerUpdate
    end type


contains

subroutine init(this, inputDir, ActuatorDisk_ID, xG, yG, zG)
    class(actuatorDisk_filtered), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDisk_ID
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0
    real(rkind) :: diam=0.08d0, cT=0.65d0, yaw=0.d0, tilt=0.d0!, Cp = 0.3
    real(rkind) :: thickness=1.5d0, filterWidth=0.5, time2initialize=0d0, tmp
    logical :: useCorrection=.TRUE., useDynamicYaw=.TRUE.

    ! Read input file for this turbine    
    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, cT, yaw, tilt, filterWidth, useCorrection, useDynamicYaw, thickness
    write(tempname,"(A13,I4.4,A10)") "ActuatorDisk_", ActuatorDisk_ID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED', action="read")
    read(unit=ioUnit, NML=ACTUATOR_DISK)
    close(ioUnit)
    
    call message(1, "Initializing Actuator Disk (ADM Type=5) number", ActuatorDisk_ID)
    call tic()
    
    ! link grids and read inputs 
    this%dx=xG(2,1,1)-xG(1,1,1)
    this%dy=yG(1,2,1)-yG(1,1,1)
    this%dz=zG(1,1,2)-zG(1,1,1)
    this%dV = this%dx*this%dy*this%dz
    this%xLoc = xLoc; this%yLoc = yLoc; this%zLoc = zLoc
    this%cT = cT; this%diam = diam; this%yaw = yaw
    this%ut = 1.d0!; this%Cp = Cp
    
    this%nxLoc = size(xG,1); this%nyLoc = size(xG,2); this%nzLoc = size(xG,3)
    
    ! Allocate stuff
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    allocate(this%rbuff(size(xG,1), size(xG,2), size(xG,3)))
    
    this%xG => xG; this%yG => yG; this%zG => zG
    this%memory_buffers_linked = .false.
    this%xLine = xG(:,1,1); this%yLine = yG(1,:,1); this%zLine = zG(1,1,:)
    
    allocate(this%powerTime(1000))  ! copied from ADM T2
    allocate(this%uTime(1000))  
    allocate(this%vTime(1000))

    ! Set thickness
    this%thick = thickness*this%dx
    this%delta = filterWidth*this%diam 
                 ! two*sqrt(this%dx**2 + this%dy**2 + this%dz**2)

    ! Get (unrotated) turbine location points
    call sample_on_circle(this%diam, this%yLoc, this%zLoc, this%ys, this%zs, this%dy, this%dz)
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
        call message(2, "Using static turbine, may not run in parallel due to memory allocation conflicts. (07/27/22)")
    end if
    
    call message(2, "Smearing grid parameter, Delta", this%delta)
    call toc(mpi_comm_world, time2initialize)
    call message(2, "Time (seconds) to initialize", time2initialize)
end subroutine 

subroutine link_memory_buffers(this, rbuff, blanks, speed, scalarSource)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(:,:,:), intent(in), target :: rbuff, blanks, speed, scalarSource
    this%rbuff => rbuff
    this%speed => speed
    this%scalarSource => scalarSource
    this%blanks => blanks
    this%memory_buffers_linked = .true. 
    
    if (.not. this%useDynamicyaw) then
        call get_weights(this, this%xs, this%ys, this%zs)!, this%scalarSource)
    end if
end subroutine 

subroutine destroy(this)
    class(actuatordisk_filtered), intent(inout) :: this

    this%memory_buffers_linked = .false.
    
    nullify(this%rbuff, this%blanks, this%speed, this%scalarSource)
    nullify(this%xG, this%yG, this%zG)
end subroutine 

! Convolution in x (streamwise) direction
subroutine get_R1(this, R1) 
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nxLoc), intent(out) :: R1
    real(rkind), dimension(this%nxLoc) :: xLine
    real(rkind) :: cbuff

    xLine = this%xLine - this%xLoc 
    cbuff = sqrt(6.d0)/this%delta
    R1 = erf(cbuff*(xLine + this%thick/two)) - &
         erf(cbuff*(xLine - this%thick/two))
    R1 = R1 / (two*this%thick)
end subroutine

! Convolution in the yz-plane (build numerically with Green's function)
subroutine get_R2(this, ys, zs, R2)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nyLoc, this%nzLoc), intent(out) :: R2
    real(rkind), dimension(this%nyLoc, this%nzLoc) :: yy, zz
    real(rkind), dimension(:), intent(in), allocatable :: ys, zs
    real(rkind) :: cbuff
    real(rkind), dimension(this%nyLoc, this%nzLoc) :: stamp
    integer :: j
    
    R2 = zero
    stamp = zero
    yy = this%yG(1, :, :)
    zz = this%zG(1, :, :)
    cbuff = -6.d0 / (this%delta**2)
    ! we need to compute the integral, iterate through points on disk face: 
    do j = 1, size(ys)
        stamp = exp(cbuff*((ys(j)-yy)**2 + (zs(j)-zz)**2))
        R2 = R2 + stamp
    end do

    cbuff = this%dy*this%dz*24.d0 / (pi**2 * this%diam**2 * this%delta**2)
    R2 = R2*cbuff  ! weight accordingly to the sampling density in sample_on_circle()
end subroutine

subroutine get_weights(this, xs, ys, zs)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(:), allocatable :: xs, ys, zs
!    real(rkind), dimension(:,:,:), intent(out), allocatable :: R
!    real(rkind), dimension(this%nxLoc,this%nyLoc,this%nzLoc), intent(out) :: R
    real(rkind), dimension(:,:,:), pointer :: R
    real(rkind), dimension(this%nyLoc, this%nzLoc) :: R2
    real(rkind), dimension(this%nxLoc) :: R1
!    real(rkind), dimension(this%nyLoc, this%nzLoc) :: yy, zz
    
    ! TODO: FIX to use xs as well
    call this%get_R2(ys,zs,R2)
    call this%get_R1(R1) 

    R => this%scalarSource
    R = spread(spread(R1,2,this%nyLoc),3,this%nzLoc) &
        * spread(R2, 1, this%nxLoc)
   
    ! minimum threshold tolerance
    where (R < 1.d-10)
        R = 0
    end where

    ! normalize so R integrates to 1 exactly
    R = R / (p_sum(R)*this%dV) !this%dx*this%dy*this%dz)
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
end subroutine

! Right hand side forcing term for the ADM
subroutine get_RHS(this, u, v, w, rhsxvals, rhsyvals, rhszvals, gamma_negative, theta)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind), dimension(:), allocatable :: xs, ys, zs
    real(rkind), intent(in) :: gamma_negative, theta
    real(rkind) :: usp_sq, force, vface!, gamma
    real(rkind), dimension(3,1) :: n=[1,0,0] !xn, Ft
    !real(rkind), dimension(3,3) :: R, T
    !real(rkind) ::  numPoints, x, y, z, scalarSource, sumVal
    !real(rkind) :: xnew,ynew,znew,cgamma,sgamma,ctheta,stheta,x2,y2,z2
    !integer :: i,j,k
    
    ! NOT IMPLEMENTED YET: YAWING (DYNAMIC OR OTHERWISE)
    
    if (this%memory_buffers_linked) then

        ! this%speed = sqrt(u**2 + v**2 + w**2)  !which one to use?     
        this%speed = u 
        if (this%useDynamicYaw) then
            call sample_on_circle(this%diam, this%yLoc, this%zLoc, ys, zs, this%dy, this%dz)
            allocate(xs(size(this%ys)))
            xs = this%xLoc
            call get_weights(this, xs, ys, zs)
        end if
        
        ! Mean speed at the turbine, corrected with factor M
        this%ut = this%M*p_sum(this%scalarSource*u)*this%dV !this%dx*this%dy*this%dz
!        ! debugging
!        write(*,*) "      ** get_RHS(): disk velocity is", this%ut
!        write(*,*) "      ** sum of disk velocities is", p_sum(this%smearing_base * u)
        usp_sq = (this%ut)**2
        force = -0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
        
        rhsxvals = rhsxvals + force * n(1,1) * this%scalarSource
        rhsyvals = rhsyvals + force * n(2,1) * this%scalarSource
        rhszvals = rhszvals + force * n(3,1) * this%scalarSource
        
        if (usp_sq /= 0.d0) then
            this%powerTime(this%tInd) = this%get_power()
            this%uTime(this%tInd) = this%ut
            this%vTime(this%tInd) = p_sum(this%scalarSource*v)*this%dV !this%dx*this%dy*this%dz
            this%tInd = this%tInd + 1
        end if
    end if

end subroutine

!subroutine get_RHS_withPower(this, u, v, w, rhsxvals, rhsyvals, rhszvals, gamma_negative, theta, wind_dir, dirType, ref_turbine)
!    class(actuatordisk_filtered), intent(inout) :: this
!    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
!    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
!    logical, intent(in) :: ref_turbine
!    real(rkind), intent(in) :: gamma_negative, theta, wind_dir
!    real(rkind) :: usp_sq, force, gamma, gamma_local
!    real(rkind), dimension(3,3) :: R, T
!    real(rkind), dimension(3,1) :: xn, Ft, n
!    real(rkind) ::  numPoints, x, y, z, scalarSource, sumVal
!    real(rkind) :: xnew, ynew, znew, cgamma, sgamma, ctheta, stheta, x2, y2, z2
!    integer :: ! i,j,k
!    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc) :: rhsxvalsg, rhsyvalsg, rhszvalsg
!    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc) :: ug, vg, wg
!    integer, intent(in) :: dirType
!
!    ug = u; vg = v; wg = w; rhsxvalsg = rhsxvals; rhsyvalsg = rhsyvals; rhszvalsg = rhszvals;
!    if (dirType==1) then
!        gamma_local = wind_dir * pi / 180.d0
!    elseif (dirType==2) then
!        gamma_local = this%hubDirection * pi / 180.d0
!    endif
!
!    ! Compute the actuator disk forcing and power
!    if (ref_turbine ) then
!        ! Use reference turbine adjacent to the leading turbine as the power
!        ! reference
!        call this%get_RHS(u, v, w, rhsxvals, rhsyvals, rhszvals, gamma_negative, theta)
!        this%powerBaseline = this%get_power() 
!    else
!        ! Use a ghost turbine with zero yaw as the reference turbine
!        ! Ghost turbine with zero yaw/tilt
!        call this%get_RHS(ug, vg, wg, rhsxvalsg, rhsyvalsg, rhszvalsg, gamma_local, theta*0.d0) 
!        this%powerBaseline = this%get_power()
!        ! Now run the model with the appropriate yaw misalignment to return the
!        ! correct values
!        call this%get_RHS(u, v, w, rhsxvals, rhsyvals, rhszvals, gamma_negative, theta) 
!    end if
!
!end subroutine

! residual from actuatordisk_yaw.f90
subroutine get_RHS_old(this, u, v, w, rhsxvals, rhsyvals, rhszvals, gamma_negative, theta)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind), intent(in) :: gamma_negative, theta
    real(rkind) :: usp_sq, force, gamma
    real(rkind), dimension(3,3) :: R, T
    real(rkind), dimension(3,1) :: xn, Ft, n
    real(rkind) ::  numPoints, x, y, z, scalarSource, sumVal
    real(rkind) :: xnew, ynew, znew, cgamma, sgamma, ctheta, stheta, x2, y2, z2
    integer :: i,j,k
     
    if (this%memory_buffers_linked) then
        ! Normal vector: gamma = yaw misalignment, theta = tilt angle
        ! Adjust the yaw misalignment angle as per the convention of
        ! Howland et al. 2019 and Shapiro et al. 2019
        gamma = -gamma_negative
        xn = reshape([1.0d0, 0.d0, 0.d0], shape(xn))
        R = reshape([cos(gamma), sin(gamma), 0.d0, -sin(gamma), cos(gamma), 0.d0, 0.d0, 0.d0, 1.d0], shape(R))
        T = reshape([1.d0, 0.d0, 0.d0, 0.d0, cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta)], shape(T)) 
        n = matmul(matmul(transpose(R), transpose(T)), xn)
        
        ! Above this is fast
        ! Translate and rotate domain
        ctheta = cos(theta); stheta = sin(theta)
        cgamma = cos(gamma); sgamma = sin(gamma)
        this%scalarSource = 0.d0

        do k = 1,size(this%xG,3)
            z = this%zG(1,1,k) - this%zLoc
            do j = 1,size(this%xG,2)
                y = this%yG(1,j,1) - this%yLoc
                !$omp simd
                do i = 1,size(this%xG,1)  
                  x = this%xG(i,1,1) - this%xLoc
                  xnew = x * cgamma - y * sgamma
                  ynew = x * sgamma + y * cgamma
                  znew = z
                  x2 = xnew
                  y2 = ynew*ctheta - znew*stheta
                  z2 = ynew*stheta + znew*ctheta
                  call this%AD_force_point(x2, y2, z2, scalarSource)
                  this%scalarSource(i,j,k) = scalarSource
                end do
            end do
        end do
        ! this part needs to be done at the end of the loops
        sumVal = p_sum(this%scalarSource) * this%dV !this%dx*this%dy*this%dz
        this%scalarSource = this%scalarSource / sumVal
        ! Get the mean velocities at the turbine face
        this%speed = u*n(1,1) + v*n(2,1) + w*n(3,1)
        !write(*,*) 'n'
        !write(*,*) n(1,1)
        !write(*,*) n(2,1)
        !write(*,*) n(3,1)
        this%blanks = 1.d0
        do k = 1,size(this%xG,3)
            do j = 1, size(this%xG,2)
                do i = 1,size(this%xG,1)
                    if (abs(this%scalarSource(i,j,k))<1D-10) then
                        this%blanks(i,j,k) = 0.d0
                    end if
                end do
            end do
        end do
        numPoints = p_sum(this%blanks)
        this%rbuff = this%blanks*this%speed
        this%ut = p_sum(this%rbuff)/numPoints    
        this%hubDirection = atan2(p_sum(this%blanks*v), p_sum(this%blanks*u)) * 180.d0 / pi
        ! Mean speed at the turbine
        usp_sq = (this%ut)**2
        force = -0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
        Ft = reshape([force, 0.d0, 0.d0], shape(Ft))
        Ft = matmul(matmul(transpose(R), transpose(T)), Ft)
        ! Can we avoid the above matmul? Just write it out
        rhsxvals = rhsxvals + Ft(1,1) * this%scalarSource 
        rhsyvals = rhsyvals + Ft(2,1) * this%scalarSource
        rhszvals = rhszvals + Ft(3,1) * this%scalarSource 

    end if 

end subroutine

! residual from actuatordisk_yaw.f90
!pure function get_power(this) result(power)
function get_power(this) result(power)
    class(actuatordisk_filtered), intent(in) :: this
    real(rkind) :: power
    ! OLD: 
    !power = 0.5d0*this%Cp*(pi*(this%diam**2)/4.d0)*(this%ut)**3
    
    ! New power assumes induction theory and uses CT' = CP'
    power = 0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*this%ut**3
end function

! residual from actuatordisk_yaw.f90
subroutine AD_force_point(this, X, Y, Z, scalarSource)
    class(actuatordisk_filtered), intent(inout) :: this
    real(rkind), intent(out) :: scalarSource
    real(rkind), intent(in) :: X,Y,Z
    real(rkind) :: delta_r = 0.8d0, delta, R!, sumVal
    real(rkind) :: tmp, tmp2
    !real(rkind) :: diamFactor = 1.75d0 ! works well for fine grids
    !real(rkind) :: diamFactor = 2.5d0, smear_x = 1.0d0 ! works well for coarse grids
    ! New parameters, 11/7/19
    real(rkind) :: diamFactor = 1.75d0, smear_x = 0.5d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Do everything in the i,j,k loop in the rhs call

    ! X,Y,Z are shifted to xc, yc, zx zero center as per AD location
    R = this%diam / 2.d0 * diamFactor ! Added to try and smear the forcing out more
    tmp = sqrt((Y/(R))**2 + (Z/(R))**2)
    tmp = (tmp-1.d0)/delta_r + 1.d0
    call Sfunc_point(tmp, tmp2)
    tmp = 1.d0 - tmp2
    delta = (this%dx)*smear_x
    scalarSource = tmp * (1.d0/(delta*sqrt(2.d0*pi))) * exp(-0.5d0*(X**2)/(delta**2))

end subroutine

! residual from actuatordisk_yaw.f90
pure subroutine Sfunc_point(x, val)
    !class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), intent(in)  :: x
    real(rkind), intent(out) :: val
    
    val = 1.d0 / (1.d0 + exp(min(1.d0/(x-1.d0+1D-18) + 1.d0/(x + 1D-18), 50.d0)))
    if (x>1.d0) then
        val = 1.d0
    end if
    if (x<0.d0) then
        val = 0.d0
    end if

end subroutine

! residual from actuatordisk_yaw.f90
subroutine dumpPower(this, outputfile, tempname)
    class(actuatordisk_filtered), intent(inout) :: this
    character(len=*),    intent(in)            :: outputfile, tempname
    integer :: fid = 1234
    character(len=clen) :: fname

    ! Write power
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname)
    !open(fid,file=trim(fname), form='unformatted',action='write',position='append')
    open(fid,file=trim(fname), form='formatted', action='write',position='append')
    write(fid, *) this%get_power()
    close(fid)

end subroutine    

! residual from actuatordisk_yaw.f90
!subroutine dumpPowerUpdate(this, outputfile, tempname, & 
!                           powerUpdate, dirUpdate, Phat, Phat_fit, yaw, yawOld, & 
!                           meanP, kw, sigma, phat_yaw, i, pBaseline, &
!                           hubDirection, Popti, stdP, alpha_m, dirStd, turbNum)
!    class(actuatordisk_filtered), intent(inout) :: this
!    character(len=*),    intent(in)            :: outputfile, tempname
!    integer :: fid = 1234
!    integer, intent(in) :: i, turbNum
!    character(len=clen) :: fname, tempname2
!    real(rkind), dimension(:), intent(in) :: powerUpdate, dirUpdate, Phat, Phat_fit, yaw, yawOld, meanP
!    real(rkind), dimension(:), intent(in) :: kw, sigma, phat_yaw, pBaseline, hubDirection
!    real(rkind), dimension(:), intent(in) :: Popti, stdP
!    real(rkind), intent(in) :: alpha_m, dirStd
!
!    ! Write power
!    write(tempname2,"(A5,I3.3,A6,I3.3,A4)") "Pvec_",i,"_turb_",turbNum,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/tdata/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) powerUpdate
!    close(fid)
!    ! Write wind direction
!    write(tempname2,"(A7,I3.3,A6,I3.3,A4)") "Dirvec_",i,"_turb_",turbNum,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/tdata/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) dirUpdate
!    close(fid)
!    ! Write Phat (this one also included wind condition distributions)
!    write(tempname2,"(A5,I3.3,A4)") "Phat_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) Phat
!    close(fid)
!    ! Write Phat_fit
!    write(tempname2,"(A5,I3.3,A4)") "Pfit_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) Phat_fit
!    close(fid)
!    ! Write yaw opti
!    write(tempname2,"(A8,I3.3,A4)") "YawOpti_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) yaw
!    close(fid)
!    ! Write yaw original in this time interval
!    write(tempname2,"(A12,I3.3,A4)") "YawInterval_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) yawOld
!    close(fid)
!    ! Write mean power in this time interval
!    write(tempname2,"(A6,I3.3,A4)") "MeanP_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) meanP / pBaseline(1)
!    close(fid)
!    ! Write std power in this time interval
!    write(tempname2,"(A5,I3.3,A4)") "stdP_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) stdP / pBaseline(1)
!    close(fid)
!    ! Write kw in this time interval
!    write(tempname2,"(A3,I3.3,A4)") "kw_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) kw
!    close(fid)
!    ! Write sigma in this time interval
!    write(tempname2,"(A6,I3.3,A4)") "sigma_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) sigma
!    close(fid)
!    ! Write phat_yaw in this time interval
!    write(tempname2,"(A8,I3.3,A4)") "phatYaw_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) phat_yaw
!    close(fid)
!    ! Write hub direction in this time interval
!    write(tempname2,"(A7,I3.3,A4)") "hubDir_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) hubDirection
!    close(fid)
!    ! Write popti in this time interval
!    write(tempname2,"(A6,I3.3,A4)") "pOpti_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) Popti
!    close(fid)
!    ! Write full power vector to file
!    !write(tempname2,"(A5,I3.3,A4)") "Pvec_",i,".txt"
!    !fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    !open(fid,file=trim(fname), form='formatted')
!    !write(fid, *) powerUpdate
!    !close(fid)
!
!    ! Data associated with the wind direction stationarity statistics
!    ! Write alpha_m in this time interval
!    write(tempname2,"(A7,I3.3,A4)") "alpham_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) alpha_m
!    close(fid)
!    ! Write dirStd in this time interval
!    write(tempname2,"(A7,I3.3,A4)") "dirstd_",i,".txt"
!    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
!    open(fid,file=trim(fname), form='formatted')
!    write(fid, *) dirStd
!    close(fid)
!
!end subroutine    

end module 
