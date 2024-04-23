module AD_Coriolis_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa, pi, one
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
    
    real(rkind), parameter :: xdim = 1000._rkind, udim = 0.45_rkind
    real(rkind), parameter :: timeDim = xdim/udim
    real(rkind), dimension(:,:,:), allocatable :: utarget, vtarget, wtarget
    real(rkind), dimension(:,:,:), allocatable :: utarget0, vtarget0, wtarget0
    ! variables set by the inputfile: 
    real(rkind) :: fringe_fact, lambdafact, freq_inflow, amplit_inflow, dt_max 
    integer :: N = 10  ! minimum number of time steps per period

    contains

subroutine init_fringe_targets(inputfile, igp)
    use exits, only: message
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max
    use IncompressibleGrid, only: igrid

    implicit none
    character(len=*),                intent(in)    :: inputfile
!    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    type(igrid), allocatable, target, intent(inout) :: igp
    real(rkind), dimension(:,:,:), pointer :: z
    real(rkind) :: Lx, Ly, Lz, uInflow, vInflow, yaw, freq, amplit
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: ioUnit
    integer :: InflowProfileType
    logical :: useGeostrophicForcing

    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw, freq, amplit

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    

    ! Initialize the velocity targets
    ! Put the same velocity profile in init fields and target
    ! Do something similar for v
    ! Compute the targets
    wtarget0 = zero 
    zMid = Lz / two
    z => igp%mesh(:,:,:,3)

    ! u,v,w targets are allocated in AD_floating.F90
    ! need to allocate u,v,w target0
    allocate(utarget0(size(utarget, 1), size(utarget, 2), size(utarget, 3)))
    allocate(vtarget0(size(vtarget, 1), size(vtarget, 2), size(vtarget, 3)))
    allocate(wtarget0(size(wtarget, 1), size(wtarget, 2), size(wtarget, 3)))
    call get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, yaw, utarget0, vtarget0) 
    
    utarget = utarget0
    vtarget = vtarget0
    wtarget = wtarget0
    
    ! get the fraction of the domain encompassed by the fringe
    call igp%fringe_x%getFringeFraction(fringe_fact)
    call igp%fringe_x%getLambdaFact(lambdafact)
    
    ! save frequency and amplitude variables: 
    freq_inflow = freq
    amplit_inflow = amplit
    if (freq_inflow == zero) then
        dt_max = 100.d0  ! don't limit dt based on 
    else
        dt_max = one / freq_inflow / N
    endif

end subroutine

subroutine update_fringe_targets(inputfile, igp)
    use exits, only: message
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use IncompressibleGrid, only: igrid

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(igrid), allocatable, target, intent(in) :: igp
    real(rkind), dimension(:,:,:), pointer :: z
    real(rkind) :: gain, lam_adjusted

    ! gain is computed from the analytical form of the system. 
    ! The gain might differ if CFL condition is used, but CFL not recommended as 
    ! the time step for numerical stability may be on the order of 1/freq.
    lam_adjusted = lambdafact * fringe_fact / igp%dt
    gain = sqrt((freq_inflow * two * pi)**two + lam_adjusted**two) / lam_adjusted
   
    ! set instantaneous u,v,w targets, modified by the transfer function 
    utarget = utarget0 * (one + amplit_inflow * sin(two*pi*freq_inflow*igp%tsim) * gain)
    vtarget = vtarget0 * (one + amplit_inflow * sin(two*pi*freq_inflow*igp%tsim) * gain)
    wtarget = wtarget0 * (one + amplit_inflow * sin(two*pi*freq_inflow*igp%tsim) * gain)

end subroutine

subroutine check_dt(igp)
    use IncompressibleGrid, only: igrid
    use exits, only: message

    implicit none
    type(igrid), allocatable, target, intent(inout) :: igp
    
    ! prevents dt from exceeding a maximum threshold defined as (1/freq)/N where
    ! N is the minimum number of samples per period 1/freq
    if (igp%dt > dt_max) then
        call message(1, "CHECK_DT: dt_max exceeded at dt=", igp%dt)
        igp%dt = dt_max
        call message(1, "CHECK_DT: Set dt to dt=", dt_max)
    end if

end subroutine 

subroutine get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, yaw, u, v)
    use kind_parameters, only: rkind
    use constants,       only: zero, one, two, pi, half
    
    implicit none
    real(rkind), dimension(:,:,:), intent(inout) :: u, v
    real(rkind), dimension(:,:,:), intent(in) :: z
    real(rkind), intent(in) :: InflowProfileAmplit, InflowProfileThick, zMid, uInflow, vInflow, yaw
    integer, intent(in) :: InflowProfileType
    integer:: i
    real(rkind) :: a_max, g_min, g_max
    real(rkind), dimension(:,:,:), allocatable :: alpha, g
    real(rkind) :: buffer=8.0d-1

    select case(InflowProfileType)
      case(0)
          u = uInflow 
          v = zero
      case(1)
          u = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          v = zero
      case(2)
          u = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          v = vInflow * tanh((z-zMid)/InflowProfileThick);
      case(3)  ! shear only (deprecated)
          u = uInflow*(one + (z-zMid)/InflowProfileThick)
          v = zero
          where (u<5.0d-1)
                u = 5.0d-1
          end where
          where (u>1.5d0)
                u = 1.5d0
          end where
      case(4)  ! veer only (deprecated)
          u = uInflow
          v = vInflow * (zMid-z)/InflowProfileThick
      case(5)  ! shear only
          u = uInflow*(one + (z-zMid)/InflowProfileThick)
          where (u<5.0d-1)
                u = 5.0d-1
          end where
          where (u>1.5d0)
                u = 1.5d0
          end where
          v = vInflow * (zMid-z)/InflowProfileThick
      case(6)  ! veer only
          alpha = (zMid-z)/InflowProfileThick*vInflow
          a_max = (pi/two)*buffer
          ! prevent any reverse flow from strong veer
          where (alpha>a_max)
                alpha = a_max
          end where
          where (alpha < -a_max)
                alpha = -a_max
          endwhere 
          u = uInflow*cos(alpha)
          v = uInflow*sin(alpha) 
      case(7)
          ! first compute the non-piecewise profiles g (vel  magnitude), alpha
          ! (vel dir)
          g = uInflow*((z-zMid)/InflowProfileThick + one)
          alpha = (zMid-z)/InflowProfileThick*vInflow

          ! limit by buffer (default 0.8)
          a_max = (pi/two)*buffer
          g_min = max(one-buffer, -abs(a_max/(vInflow+1d-18)) + one)
          g_max = min(one+buffer, abs(a_max/(vInflow+1d-18)) + one)

          ! ensure there is no reverse flow in the domain: enforce the bounds on
          ! g=|u| and alpha set by buffer
          do i=1, size(z,3)
              if (g(1,1,i) < g_min) then
                  g(:,:,i) = g_min
                  alpha(:,:,i) = vInflow*(one-g_min)
               else if (g(1,1,i) > g_max) then
                   g(:,:,i) = g_max
                   alpha(:,:,i) = vInflow*(one-g_max)
              end if    
          end do
          ! set the velocity components u and v: 
          u = g*cos(alpha)
          v = g*sin(alpha)
      case(8)
          ! Uniform yawed inflow
          u = uInflow*cos(yaw*pi/180.d0)
          v = -uInflow*sin(yaw*pi/180.d0)
    end select
end subroutine

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use AD_Coriolis_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one,two
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = one, Ly = one, Lz = one, yaw, amplit, freq
    real(rkind) :: uInflow, vInflow  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick
    integer :: InflowProfileType
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw, amplit, freq

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    

    !Lx = two*pi; Ly = two*pi; Lz = one

    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nxg,rkind)
        dy = Ly/real(nyg,rkind)
        dz = Lz/real(nzg,rkind)

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
                end do
            end do
        end do

        ! Shift everything to the origin 
        x = x - dx
        y = y - dy
        z = z - dz 

    end associate

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use AD_Coriolis_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind), dimension(:,:,:), allocatable :: randArr, ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE
    real(rkind)  :: Lx = one, Ly = one, Lz = one, G_alpha, yaw, amplit, freq
    real(rkind) :: uInflow, vInflow  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw, amplit, freq

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)

    wC = zero
    zMid = Lz / 2.d0
    
    ! initialize inflow profile
    call get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, yaw, u, v)
    
    !allocate(randArr(size(u,1),size(u,2),size(u,3)))
    !call gaussian_random(randArr,-one,one,seedu + 10*nrank)
    !!do k = 1,size(randArr,3)
    !!     u(:,:,k) = u(:,:,k) + randscale*randArr(:,:,k)
    !!end do
    !deallocate(randArr)

    call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
    call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
    call message_min_max(1,"Bounds for w:", p_minval(minval(w)), p_maxval(maxval(w)))

    ! Interpolate wC to w
    allocate(ybuffC(decompC%ysz(1),decompC%ysz(2), decompC%ysz(3)))
    allocate(ybuffE(decompE%ysz(1),decompE%ysz(2), decompE%ysz(3)))

    allocate(zbuffC(decompC%zsz(1),decompC%zsz(2), decompC%zsz(3)))
    allocate(zbuffE(decompE%zsz(1),decompE%zsz(2), decompE%zsz(3)))
   
    nz = decompC%zsz(3)
    nzE = nz + 1

    call transpose_x_to_y(wC,ybuffC,decompC)
    call transpose_y_to_z(ybuffC,zbuffC,decompC)
    zbuffE = zero
    zbuffE(:,:,2:nzE-1) = half*(zbuffC(:,:,1:nz-1) + zbuffC(:,:,2:nz))
    call transpose_z_to_y(zbuffE,ybuffE,decompE)
    call transpose_y_to_x(ybuffE,w,decompE) 
    
    

    deallocate(ybuffC,ybuffE,zbuffC, zbuffE) 
  
      
    nullify(u,v,w,x,y,z)
   

    call message(0,"Velocity Field Initialized")

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [308, 512]
    yplanes = [512]
    zplanes = [256]

end subroutine

subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

end subroutine

subroutine setInhomogeneousNeumannBC_Temp(inputfile, wTh_surf)
    use kind_parameters,    only: rkind
    use constants, only: one, zero 
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: wTh_surf
    integer :: ioUnit 
    real(rkind) :: ThetaRef, Lx, Ly, Lz, yaw, amplit, freq
    logical :: initPurturbations = .false. 
    real(rkind) :: uInflow, vInflow  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw, amplit, freq
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, G_alpha, yaw, amplit, freq
    integer :: iounit
    real(rkind) :: uInflow, vInflow  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw, amplit, freq
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, G_alpha, yaw, amplit, freq
    integer :: iounit
    real(rkind) :: uInflow, vInflow  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw, amplit, freq

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    
     
    Tref = 0.d0
    
    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
    character(len=*),                intent(in)    :: inputfile
    integer, parameter :: nprobes = 2
    
    ! IMPORTANT : Convention is to allocate probe_locs(3,nprobes)
    ! Example: If you have at least 3 probes:
    ! probe_locs(1,3) : x -location of the third probe
    ! probe_locs(2,3) : y -location of the third probe
    ! probe_locs(3,3) : z -location of the third probe


    ! Add probes here if needed
    ! Example code: The following allocates 2 probes at (0.1,0.1,0.1) and
    ! (0.2,0.2,0.2)  
    print*, inputfile
    allocate(probe_locs(3,nprobes))
    probe_locs(1,1) = 0.1d0; probe_locs(2,1) = 0.1d0; probe_locs(3,1) = 0.1d0;
    probe_locs(1,2) = 0.2d0; probe_locs(2,2) = 0.2d0; probe_locs(3,2) = 0.2d0;


end subroutine

subroutine initScalar(decompC, inpDirectory, mesh, scalar_id, scalarField)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    type(decomp_info),                                          intent(in)    :: decompC
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarField

    scalarField = 0.d0
end subroutine 

subroutine setScalar_source(decompC, inpDirectory, mesh, scalar_id, scalarSource)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    type(decomp_info),                                          intent(in)    :: decompC
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarSource

    scalarSource = 0.d0
end subroutine 


