module convective_igrid_parameters

      ! TAKE CARE OF TIME NON-DIMENSIONALIZATION IN THIS MODULE

    use exits, only: message, GraceFulExit
    use kind_parameters,  only: rkind
    use constants, only: zero, kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
    
    real(rkind), parameter :: xDim = 1000._rkind, uDim = 12._rkind !sqrt(3.0_rkind**2+9._rkind**2) ! xDim is set to 100 so that any turbine can be used
    real(rkind), parameter :: timeDim = xDim/uDim

end module     

module stratified_pblBCsmod
    use kind_parameters, only: rkind, clen
    use basic_io, only: read_2d_ascii 
    use interpolation, only: spline, ispline, binarysearch 
    real(rkind), dimension(:), allocatable :: t_geo, G_geo, t_flux, wT_flux, a_geo, b_geo, c_geo, a_flux, b_flux, c_flux
    real(rkind), dimension(:), allocatable :: t_galpha, G_alpha, a_galpha, b_galpha, c_galpha
    real(rkind) :: G_tolerance = 0.1d0
contains


    subroutine setup_stratified_pblBCs(inputfile)
        character(len=*),                intent(in)    :: inputfile
        real(rkind), dimension(:,:), allocatable :: data2read
        character(len=clen) :: fname_G, fname_wtheta, fname_galpha
        integer :: ioUnit
        real(rkind), parameter :: xDim = 1000._rkind, uDim = 12._rkind !sqrt(3.0_rkind**2+9._rkind**2)
        real(rkind) :: G_BC = 8.0, wtheta_BC = 0.0, g_alpha_BC = 0.0
        logical :: read_BCs = .true.
        namelist /stratified_pbl_BCS/ fname_G, fname_wtheta, fname_galpha, G_tolerance, read_BCs, G_BC, wtheta_BC, g_alpha_BC

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=STRATIFIED_PBL_BCS)
        close(ioUnit)    

        if (read_BCs) then
            ! Geostrophic wind speed
            call read_2d_ascii(data2read, fname_G)
            allocate(t_geo(size(data2read,1)))
            allocate(G_geo(size(data2read,1)))
            allocate(a_geo(size(data2read,1)))
            allocate(b_geo(size(data2read,1)))
            allocate(c_geo(size(data2read,1)))
            t_geo = data2read(:,1) * 60.d0*60.d0/xDim
            G_geo = data2read(:,2)
            call spline(t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
            deallocate(data2read)

            ! Surface sensible heat flux <w'theta'>
            call read_2d_ascii(data2read, fname_wtheta)
            allocate(t_flux(size(data2read,1)))
            allocate(wT_flux(size(data2read,1)))
            allocate(a_flux(size(data2read,1)))
            allocate(b_flux(size(data2read,1)))
            allocate(c_flux(size(data2read,1)))
            t_flux = data2read(:,1) * 60.d0*60.d0/xDim
            wT_flux = data2read(:,2)
            call spline(t_flux, wT_flux, a_flux, b_flux, c_flux, size(t_flux))
            deallocate(data2read)

            ! Geostrophic wind direction in computational domain g_alpha 
            ! NOTE: This is different from the frameAngle variable which describes
            ! the geostrophic wind direction in cardinal directions
            call read_2d_ascii(data2read, fname_galpha)
            allocate(t_galpha(size(data2read,1)))
            allocate(G_alpha(size(data2read,1)))
            allocate(a_galpha(size(data2read,1)))
            allocate(b_galpha(size(data2read,1)))
            allocate(c_galpha(size(data2read,1)))
            t_galpha = data2read(:,1) * 60.d0*60.d0/xDim
            G_alpha = data2read(:,2)
            call spline(t_galpha, G_alpha, a_galpha, b_galpha, c_galpha, size(t_galpha))
            deallocate(data2read)
        else
            ! Geostrophic wind speed
            allocate(t_geo(2))
            allocate(G_geo(2))
            allocate(a_geo(2))
            allocate(b_geo(2))
            allocate(c_geo(2))
            t_geo(1) = 0.d0 ! TODO this is a workaroung when the BCs are not read from a file 
            t_geo(2) = 100000000.d0           
            G_geo = G_BC
            call spline(t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))

            
            ! Surface sensible heat flux <w'theta'>
            allocate(t_flux(2))
            allocate(wT_flux(2))
            allocate(a_flux(2))
            allocate(b_flux(2))
            allocate(c_flux(2))
            t_flux(1) = 0.d0 ! TODO this is a workaroung when the BCs are not read from a file 
            t_flux(2) = 100000000.d0           
            wT_flux = wtheta_BC
            call spline(t_flux, wT_flux, a_flux, b_flux, c_flux, size(t_flux))


            ! Geostrophic wind direction in computational domain g_alpha 
            ! NOTE: This is different from the frameAngle variable which describes
            ! the geostrophic wind direction in cardinal directions
            allocate(t_galpha(2))
            allocate(G_alpha(2))
            allocate(a_galpha(2))
            allocate(b_galpha(2))
            allocate(c_galpha(2))
            t_galpha(1) = 0.d0 ! TODO this is a workaroung when the BCs are not read from a file 
            t_galpha(2) = 10000000.d0
            G_alpha = g_alpha_BC
            call spline(t_galpha, G_alpha, a_galpha, b_galpha, c_galpha, size(t_galpha))


        end if


    end subroutine

    subroutine stratified_pblBCs_CorrectnessCheck(time, G)
        use exits, only: GracefulExit
        real(rkind), intent(in) :: time, G
        real(rkind) :: Gtrue

        ! Test function
        !call linear_interp(size(t_geo,1),t_geo,G_geo,10.8d0,Gtrue)
        !Gtrue = ispline(10.8d0, t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
        !Gtrue = ispline(10.8d0, t_flux, wT_flux, a_flux, b_flux, c_flux, size(t_flux))
        !call linear_interp(size(t_flux,1),t_flux,wT_flux,10.8d0,Gtrue)
 
        ! Interpolate G: 
        !Gtrue = ispline(time, t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
        call linear_interp(size(t_geo,1),t_geo,G_geo,time,Gtrue)
        
        if (abs(Gtrue - G) .ge. G_tolerance) then 
            print*, "time:", time
            print*, "G entered:", G
            print*, "G from BC:", Gtrue
            call gracefulExit("Incorrect Geostrophic inputs",45)
        end if 

    end subroutine     
    
    subroutine get_stratified_pblBCs(time, G, wTheta, Gangle, useControl)
        real(rkind), intent(in) :: time
        logical, intent(in) :: useControl
        real(rkind), intent(out) :: G, wTheta, Gangle 

        ! Convert units for time? 
        ! This assumes that the input files' time vector is appropriately
        ! nondimensionalized       
         
        ! Interpolate G: 
        ! Only works for the stratified_pbl case based on 1 hour of unstable
        ! initialization
        ! I directly modified the RESTART file to be zero time!!!
        !G = ispline(time, t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
        ! Interpolate wTheta: 
        !wTheta = ispline(time, t_flux, wT_flux, a_flux, b_flux, c_flux, size(t_flux))
 
        ! Linear inteprolation
        call linear_interp(size(t_geo,1),t_geo,G_geo,time,G)
        call linear_interp(size(t_flux,1),t_flux,wT_flux,time,wTheta)

        if (.not. usecontrol) then
            call linear_interp(size(t_galpha,1),t_galpha,G_alpha,time,Gangle)
        end if  

    end subroutine 

    
    subroutine linear_interp(xlen,x,y,xv,yv)
        use kind_parameters, only: rkind
        use decomp_2d,        only: decomp_info
        implicit none
        integer, intent(in)                   :: xlen
        real(rkind), dimension(xlen), intent(in) :: x, y
        real(rkind), intent(in)               :: xv
        real(rkind), intent(out)              :: yv
        integer                               :: ind
     
        ! binary search
        call binarysearch(xlen, x, xv, ind, 1.d-6)
    
        ! Linear interpolation
        yv = y(ind) + (xv-x(ind)) * (y(ind+1)-y(ind)) / (x(ind+1)-x(ind))
    
    end subroutine

end module 

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use convective_igrid_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half, three
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random, uniform_random
    use decomp_2d          
    use reductions,         only: p_maxval
    use constants,          only: pi
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, T, x, y, z
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE, ztmp
    integer :: nz, nzE, k, init_T_profile
    real(rkind) :: sig, ebar, zH
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, dTsurf_dt = -0.05d0, z0init = 1.d-4, Tsurf0 = 290.d0, frameAngle = 0.d0
    real(rkind), dimension(:,:,:), allocatable :: randArr1, randArr2, randArr3
    
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0, init_T_profile, frameAngle

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    
    
    !!!!!!!!!!!!!!!!!!!!! DON'T CHANGE THE POINTERS / ALLOCATIONS !!!!!!!!!!!!!!!!!!!!!!
    u  => fieldsC(:,:,:,1); v  => fieldsC(:,:,:,2); wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1); T  => fieldsC(:,:,:,7) 
    z => mesh(:,:,:,3); y => mesh(:,:,:,2); x => mesh(:,:,:,1)
    !allocate(randArr(size(T,1),size(T,2),size(T,3)))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Add random numbers
    allocate(randArr1(size(T,1),size(T,2),size(T,3)))
    allocate(randArr2(size(T,1),size(T,2),size(T,3)))
    allocate(randArr3(size(T,1),size(T,2),size(T,3)))
    call uniform_random(randArr1,-half,half,seedu + 10*nrank)
    call uniform_random(randArr2,-half,half,seedv + 10*nrank)
    call uniform_random(randArr3,-half,half,seedw + 10*nrank)

    u = uDim
    v = zero/uDim
    wC = zero/uDim

    ! Added to account for frame angle    
    u = u * cos(frameAngle * pi / 180.d0)
    v = v * sin(frameAngle * pi / 180.d0) 

    allocate(ztmp(decompC%xsz(1),decompC%xsz(2),decompC%xsz(3)))
    ztmp = z*xDim

    select case(init_T_profile)
    
    ! Linear profile (best for neutral case)
    case(0)
        T = 0.003d0*ztmp + Tsurf0

    ! Neutral profile for the unstable case
    case(1)
        do k = 1, decompC%xsz(3)
            if (ztmp(1,1,k)<700.d0) then
                T(:,:,k) = Tsurf0
            elseif (ztmp(1,1,k)>=700.d0 .and. ztmp(1,1,k)<800.d0) then
                T(:,:,k) = Tsurf0 + (ztmp(1,1,k)-700.d0) * 3.d0 / 100.d0
            elseif (ztmp(1,1,k)>=800.d0) then
                T(:,:,k) = Tsurf0 + (800.d0 - 700.d0) * 3.d0 / 100.d0 + (ztmp(1,1,k)-800.d0) * 0.003d0
            end if
        end do
    
    ! Gabls modified initial profile for the stable case
    case(2)
        T = 0.003d0*(ztmp - 100.d0) + Tsurf0
        where(ztmp < 100.d0)
            T = Tsurf0
        end where
        ! T = T + 0.0001d0*ztmp

        if(ztmp(1,1,k) .le. zH) then 
            ebar = 0.5d0*(1.0d0-ztmp(1,1,k)/zH)
         else
            ebar = 0.0d0
         endif
         sig = sqrt(two*ebar/three) * sqrt(12.0d0) ! 1/sqrt(12) is stdev of uniform random distribution of width 1.
         sig = sig/uDim
         T(:,:,k) =  T(:,:,k) + sig*randArr1(:,:,k)
    
    ! Unstable initialization from Moeng and Sullivan 1994
    case(3)
        u = 10.d0
        T = 300.d0
        do k = 1, decompC%xsz(3)
            if (ztmp(1,1,k)>937.d0 .and. ztmp(1,1,k)<=1062) then
                T(:,:,k) = 300.d0 + (8.d0 / 125.d0) * (ztmp(1,1,k) - 937)
            elseif (ztmp(1,1,k)>=1062.d0) then
                T(:,:,k) = 308.d0 + 0.003d0 * ztmp(1,1,k)
            end if
        end do 
    
    ! Neutral initialization from neutral_pbl
    case(4)
        T = 0.003d0*(ztmp - 700.d0) + 300.d0
        where(ztmp < 700.d0)
            T = 300.d0
        end where

    ! Neutral initialization from neutral_pbl scaled for smaller z0
    case(5)
        T = 0.003d0*(ztmp - 400.d0) + Tsurf0
        where(ztmp < 400.d0)
            T = Tsurf0
        end where

    ! This is the original from igridSGS-MFH
    ! From Kumar et al. (2006)
    case(6)
    do k = 1, decompC%xsz(3)
       if (ztmp(1,1,k)<800.d0) then
           T(:,:,k) = Tsurf0
       elseif (ztmp(1,1,k)>=800.d0 .and. ztmp(1,1,k)<1000.d0) then
           T(:,:,k) = Tsurf0 + (ztmp(1,1,k)-800.d0) * 6.d0 / 200.d0
       elseif (ztmp(1,1,k)>=1000.d0) then
           T(:,:,k) = Tsurf0 + (1000.d0 - 800.d0) * 6.d0 / 200.d0 + (ztmp(1,1,k)-1000.d0) * 1.d0 / 100.d0
       end if
    end do

    ! Neutral profile for the unstable case from Pino and VILÀ-GUERAU DE ARELLANO
    case(7)
    do k = 1, decompC%xsz(3)
        if (ztmp(1,1,k)<625.d0) then
            T(:,:,k) = Tsurf0
        elseif (ztmp(1,1,k)>=625.d0 .and. ztmp(1,1,k)<825.d0) then
            T(:,:,k) = Tsurf0 + (ztmp(1,1,k)-625.d0) * 6.d0 / 200.d0
        elseif (ztmp(1,1,k)>=825.d0) then
            T(:,:,k) = Tsurf0 + 6.d0 + (ztmp(1,1,k)-825.d0) * 0.003d0
        end if
    end do
    ! Neutral profile for the unstable case from Pino and VILÀ-GUERAU DE ARELLANO with higher BL
    case(8)
    do k = 1, decompC%xsz(3)
        if (ztmp(1,1,k)<800.d0) then
            T(:,:,k) = Tsurf0
        elseif (ztmp(1,1,k)>=800.d0 .and. ztmp(1,1,k)<1000.d0) then
            T(:,:,k) = Tsurf0 + (ztmp(1,1,k)-800.d0) * 6.d0 / 200.d0
        elseif (ztmp(1,1,k)>=1000.d0) then
            T(:,:,k) = Tsurf0 + 6.d0 + (ztmp(1,1,k)-1000.d0) * 0.003d0
        end if
    end do

    case default
        call gracefulExit("Invalid choice for initial potential temperature profile",423)
    end select


    !do k = 1, decompC%xsz(3)
    !  if(ztmp(1,1,k) < 800.d0) then
    !    T(:,:,k) = 301.D0
    !  elseif(ztmp(1,1,k) < 850.d0) then
    !    T(:,:,k) = 286.0d0
    !  elseif(ztmp(1,1,k) < 900.d0) then
    !    T(:,:,k) = 286.0d0 + (ztmp(1,1,k)-850.d0)/50.0d0*2.0d0
    !  elseif(ztmp(1,1,k) < 1000.d0) then
    !    T(:,:,k) = 288.0d0 + (ztmp(1,1,k)-900.d0)/100.0d0*4.0d0
    !  else
    !    T(:,:,k) = 292.0d0 + (ztmp(1,1,k)-1000.d0)/1000.d0*8.0d0
    !  endif
    !enddo

    select case(init_T_profile)
    
        ! Linear profile (best for neutral case)
        case(0)
            zH = 700.d0
    
        ! Neutral profile for the unstable case
        case(1)
            zH = 700.d0
        
        ! Gabls modified initial profile for the stable case
        case(2)
            zH = 500.d0
        
        ! Unstable initialization from Moeng and Sullivan 1994
        case(3)
            zH = 937.d0
        
        ! Neutral initialization from neutral_pbl
        case(4)
            zH = 700.d0
    
        ! Neutral initialization from neutral_pbl scaled for smaller z0
        case(5)
            zH = 400.d0

        ! This is the original from igridSGS-MFH
        ! From Kumar et al. (2006)
        case(6)
            zH = 800.d0

        ! From Pino
        case(7)
            zH = 825.d0

        case(8)
            zH = 1000.d0

        case default
            call gracefulExit("Invalid choice for initial potential temperature profile",423)
        end select
    

    do k = 1,size(u,3)
        if(ztmp(1,1,k) .le. zH) then 
           ebar = 0.5d0*(1.0d0-ztmp(1,1,k)/zH)
        else
           ebar = 0.0d0
        endif
        sig = sqrt(two*ebar/three) * sqrt(12.0d0) ! 1/sqrt(12) is stdev of uniform random distribution of width 1.
        sig = sig/uDim
        !write(*,'(5(e19.12,1x))') ztmp(1,1,k), ebar, sig, maxval(randArr1), minval(randArr1)
         u(:,:,k) =  u(:,:,k) + sig*randArr1(:,:,k)
         v(:,:,k) =  v(:,:,k) + sig*randArr2(:,:,k)
        wC(:,:,k) = wC(:,:,k) + sig*randArr3(:,:,k)
    end do
    deallocate(randArr3)
    deallocate(randArr2)
    deallocate(randArr1)
    deallocate(ztmp)

    !!!!!!!!!!!!!!!!!!!!! DON'T CHANGE ANYTHING UNDER THIS !!!!!!!!!!!!!!!!!!!!!!
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
    ! Deallocate local memory 
    deallocate(ybuffC,ybuffE,zbuffC, zbuffE) 
    nullify(u,v,w,x,y,z)
    call message(0,"Velocity Field Initialized")
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine

subroutine setInhomogeneousNeumannBC_Temp(inputfile, wTh_surf)
    use kind_parameters,    only: rkind
    use convective_igrid_parameters
    use constants, only: one, zero 
    use stratified_pblBCsmod
    implicit none
    real(rkind), intent(inout) :: wTh_surf
    character(len=*),                intent(in)    :: inputfile
    real(rkind) :: tmp, tmp2, time

    time = wTh_surf ! temporary workaround

    call get_stratified_pblBCs(time, tmp, wTh_surf, tmp2, .false.)  ! << will be an issue for restarts

end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use convective_igrid_parameters
    use constants, only: one, zero 
    implicit none
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    character(len=*),                intent(in)    :: inputfile
    integer :: ioUnit,  init_T_profile
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, z0init = 1.d-4, Tsurf0 = 290.0d0, frameAngle = 0.d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0, init_T_profile, frameAngle
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    

    ! Do nothing here
    Tsurf = Tsurf0
    dTsurf_dt = zero
end subroutine

subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [64]
    yplanes = [64]
    zplanes = [256]

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
    allocate(probe_locs(3,nprobes))
    probe_locs(1,1) = 0.1d0; probe_locs(2,1) = 0.1d0; probe_locs(3,1) = 0.1d0;
    probe_locs(1,2) = 0.2d0; probe_locs(2,2) = 0.2d0; probe_locs(3,2) = 0.2d0;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! THE SUBROUTINES UNDER THIS DON'T TYPICALLY NEED TO BE CHANGED !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use convective_igrid_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one,two
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn, init_T_profile
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, dTsurf_dt = -0.05d0, z0init = 1.d-4, Tsurf0 = 290.d0, frameAngle = 0.d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0, init_T_profile, frameAngle

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
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

subroutine set_Reference_Temperature(inputfile, Thetaref)
    use kind_parameters,    only: rkind
    use constants, only: one, zero 
    implicit none
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Thetaref
    integer :: ioUnit, init_T_profile
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, dTsurf_dt = -0.05d0, z0init = 1.d-4, Tsurf0 = 290.d0, frameAngle = 0.d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0, init_T_profile, frameAngle
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    

    Thetaref = Tref
    ! This will set the value of Tref.     

end subroutine


subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

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






