module HIT_shear_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa, zero
    implicit none

    ! I realize it is probably bad practice to store information here, but it is the easiest way I've found
    ! -KSH 10/02/2024

    integer :: simulationID = 0
    integer :: nxADSim, nyADSim, nzADSim, nxHITSim, nyHITSim, nzHITSim
    real(rkind), dimension(:,:,:), allocatable :: utarget0, vtarget0, wtarget0   ! u, v, w laminar fringe targets
    real(rkind) :: InflowSpeed = 1.d0, TI = -1, TI_fact = -1, Kp_TI = 0.1d0
    real(rkind), dimension(:,:,:), allocatable :: z_global, utarget_1d, vtarget_1d, wtarget_1d  ! global z-axis of shape (1,1,nz)
    logical :: inflow_varies_in_z = .false., debug_TI_gain = .true.
contains

! build the velocity profiles
    subroutine get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, yaw, u, v)
        use kind_parameters, only: rkind
        use constants,       only: zero, one, two, pi, half
        use exits,           only: gracefulExit

        implicit none
        real(rkind), dimension(:,:,:), intent(inout) :: u, v
        real(rkind), dimension(:,:,:), intent(in) :: z
        real(rkind), intent(in) :: InflowProfileAmplit, InflowProfileThick, zMid, uInflow, vInflow, yaw
        integer, intent(in) :: InflowProfileType
        integer:: i
        real(rkind) :: a_max, g_min, g_max
        real(rkind), dimension(:,:,:), allocatable :: alpha, g
        real(rkind) :: buffer=8.0d-1  ! buffer value = 1 - umin

        select case(InflowProfileType)
          case(0)
            u = uInflow
            v = zero
          case(1)  ! tanh shear and veer
            u = uInflow
            v = uInflow * buffer * tanh(vinflow * InflowProfileAmplit * (z-zMid) / buffer)
          case(2)
            u = uInflow*(one  + buffer * tanh(InflowProfileAmplit * (z-zMid) / buffer))
            v = uInflow * buffer * tanh(vinflow * InflowProfileAmplit * (z-zMid) / buffer)
          case(3)  ! shear only (deprecated)
            call GracefulExit("inflow 3 deprecated", 999)
          case(4)  ! veer only (deprecated)
            call GracefulExit("inflow 4 deprecated", 999)
          case(5)  ! shear only
            call GracefulExit("inflow 5 deprecated", 999)
          case(6)  ! veer only - FIX THIS
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
          case(7)  ! shear and veer - FIX THIS
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

! makes global fringe targets necessary for do_phaseshifting (spectral z-decomp workaround)
    subroutine make_global_zaxis(adsim)
        use IncompressibleGrid, only : igrid
        use constants, only : zero, half

        type(igrid), allocatable, target :: adsim
        integer :: i

        allocate(z_global(1,1,adsim%nz), utarget_1d(1,1,adsim%nz), vtarget_1d(1,1,adsim%nz), wtarget_1d(1,1,adsim%nz))  ! all use the "global" z-axis
        z_global = reshape((/(real(i) - half, i=1, adsim%nz)/) * adsim%dz, (/1,1,adsim%nz/))  ! (1,1,nz) array
        utarget_1d = zero
        vtarget_1d = zero
        wtarget_1d = zero
    end subroutine

! fringe function
    pure subroutine Sfunc(x, val)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind), dimension(:,:,:), intent(out) :: val

        val = 0.d0
        where (x>0.d0)
            val = 1.d0/(1.d0 + exp(min(1.d0/(x - 1.d0 + 1.d-18) + 1.d0/(x + 1.d-18),50.d0)))
        end where

        where (x>1.d0)
            val = 1.d0
        end where

    end subroutine

! initialize fringe targets with (laminar) flow
    subroutine init_fringe_targets(inputfile, mesh)
        use exits, only: message
        use kind_parameters,    only: rkind
        use constants,          only: zero, one, two, pi, half
        use random,             only: gaussian_random

        implicit none
        character(len=*),                intent(in)    :: inputfile
        real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
        real(rkind), dimension(:,:,:), pointer :: z
        integer :: ioUnit
        real(rkind) :: Lx, Ly, Lz, uInflow = one, vInflow = zero, yaw = zero
        real(rkind) :: InflowProfileAmplit = one, InflowProfileThick = zero, zmid=-1
        integer :: InflowProfileType = 1

        namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, zmid, &
            InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=AD_CoriolisINPUT)
        close(ioUnit)

        ! Initialize the velocity targets for the fringe (without HIT)
        wtarget0 = zero
        if (zmid < 0) then
            zMid = Lz / two
        end if
        z => mesh(:,:,:,3)
        call get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, yaw, utarget0, vtarget0)
        call get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z_global, zMid, InflowProfileType, yaw, utarget_1d, vtarget_1d)

        InflowSpeed = uInflow  ! set the "linear advection" velocity
        if (InflowProfileType == 0) then
            inflow_varies_in_z = .false.
        else
            inflow_varies_in_z = .true.
        endif

        ! The velocity profile in z needs to go to slip wall at the top
        ! Both u and v need slip conditions
    end subroutine

! Do phase shifting here
    subroutine do_phaseshifting(hitsim, adsim, u, v, w)
        use kind_parameters,  only: rkind
        use IncompressibleGrid, only: igrid

        type(igrid), allocatable, target :: hitsim, adsim
        real(rkind), dimension(:,:,:), allocatable, intent(inout) :: u, v, w
        real(rkind), dimension(size(z_global,3)) :: x_shift_z, y_shift_z
        real(rkind) :: x_shift

        if (inflow_varies_in_z) then
            ! need to take the full z-domain
            x_shift_z = adsim%tsim * utarget_1d(1, 1,:)
            y_shift_z = adsim%tsim * vtarget_1d(1, 1,:)

            call hitsim%spectC%bandpassFilter_and_phaseshift_z(hitsim%whatC, u(nxADSim-nxHITSim+1:nxADSim,:,:), x_shift_z, y_shift_z)
            call hitsim%interpolate_cellField_to_edgeField(u(nxADSim-nxHITSim+1:nxADSim,:,:), w(nxADSim-nxHITSim+1:nxADSim,:,:),0,0)
            call hitsim%spectC%bandpassFilter_and_phaseshift_z(hitsim%uhat, u(nxADSim-nxHITSim+1:nxADSim,:,:), x_shift_z, y_shift_z)
            call hitsim%spectC%bandpassFilter_and_phaseshift_z(hitsim%vhat, v(nxADSim-nxHITSim+1:nxADSim,:,:), x_shift_z, y_shift_z)
        else
            ! Set the true target field for AD simulation
            x_shift = adsim%tsim * InflowSpeed

            ! old code:
            call hitsim%spectC%bandpassFilter_and_phaseshift(hitsim%whatC, u(nxADSim-nxHITSim+1:nxADSim,:,:), x_shift)
            call hitsim%interpolate_cellField_to_edgeField(u(nxADSim-nxHITSim+1:nxADSim,:,:), w(nxADSim-nxHITSim+1:nxADSim,:,:),0,0)
            call hitsim%spectC%bandpassFilter_and_phaseshift(hitsim%uhat, u(nxADSim-nxHITSim+1:nxADSim,:,:), x_shift)
            call hitsim%spectC%bandpassFilter_and_phaseshift(hitsim%vhat, v(nxADSim-nxHITSim+1:nxADSim,:,:), x_shift)

        end if
        ! Now scale rhw HIT field appropriately
        u(nxADSim-nxHITSim+1:nxADSim,:,:) = u(nxADSim-nxHITSim+1:nxADSim,:,:)*TI_fact + utarget0(nxADSim-nxHITSim+1:nxADSim,:,:)
        v(nxADSim-nxHITSim+1:nxADSim,:,:) = v(nxADSim-nxHITSim+1:nxADSim,:,:)*TI_fact + vtarget0(nxADSim-nxHITSim+1:nxADSim,:,:)
        w(nxADSim-nxHITSim+1:nxADSim,:,:) = w(nxADSim-nxHITSim+1:nxADSim,:,:)*TI_fact
    end subroutine

! Update TI gain
    subroutine update_TI_fact(adsim, xid)
        use IncompressibleGrid, only : igrid
        use constants, only          : zero, one, two, three
        use reductions, only         : p_sum
        use exits, only              : message

        type(igrid), allocatable, target :: adsim
        integer, intent(in) :: xid
        real(rkind), dimension(adsim%gpC%xsz(2), adsim%gpC%xsz(3)) :: buff1, buff2
        real(rkind) :: TI_inst

        if (TI < 0) then
            TI_fact = one
            return  ! Doesn't compute anything
        end if

        ! need to compute TKE, TI
        buff1 = 0.5 * ((adsim%u(xid,:,:)-utarget0(xid,:,:))**2 + (adsim%v(xid,:,:)-vtarget0(xid,:,:))**2 + (adsim%wC(xid,:,:))**2)  ! TKE
        buff2 = sqrt(utarget0(xid,:,:)**2 + vtarget0(xid,:,:)**2)  ! U_inf velocity
        buff1 = sqrt(two / three * buff1) / buff2
        TI_inst = p_sum(buff1) / (adsim%ny*adsim%nz)  ! mean TI at the given xid
        TI_fact = max(zero, TI_fact + ((TI - TI_inst) * Kp_TI))

        if (debug_TI_gain) then
            call message(1, "update_TI: TI_inst", TI_inst)
            call message(1, "update_TI: TI_fact", TI_fact)
        end if
    end subroutine

end module  ! end module functions


subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use HIT_shear_parameters
    use kind_parameters,  only: rkind
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    character(len=*),                intent(in)    :: inputfile
    integer :: i,j,k, ioUnit
    integer :: nxg, nyg, nzg
    integer :: ix1, ixn, iy1, iyn, iz1, izn

    real(rkind) :: Lx, Ly, Lz, uInflow, vInflow, yaw
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zmid=-1
    integer :: InflowProfileType

    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, zmid, &
        InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw

    select case (simulationID)
      case (1)
        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=AD_CoriolisINPUT)
        close(ioUnit)
      case (2)
        Lx = two*pi
        Ly = two*pi
        Lz = two*pi  ! this is fixed as a (2*pi)^3 box  TODO
    end select

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

    if (simulationID == 1) then  ! primary simulation (AD)
        nxADSim = nxg; nyADSim = nyg; nzADSim = nzg
    elseif (simulationID == 2) then  ! HIT simulation
        nxHITSim = nxg; nyHITSim = nyg; nzHITSim = nzg
    end if

    call message(0, "meshgen_wallM: initialized grid")
    call message(1, "nx", nxg)
    call message(1, "ny", nyg)
    call message(1, "nz", nzg)

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use HIT_shear_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max, gracefulExit, message
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    integer :: ioUnit

    real(rkind) :: Lx, Ly, Lz, uInflow = one, vInflow = zero, yaw = zero
    real(rkind) :: InflowProfileAmplit = zero, InflowProfileThick = zero, zmid=-1
    integer :: InflowProfileType = 0

    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, zmid, &
        InflowProfileAmplit, InflowProfileThick, InflowProfileType, yaw

    if (simulationID == 1) then ! for adsim only

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=AD_CoriolisINPUT)
        close(ioUnit)

        ! link pointers
        u  => fieldsC(:,:,:,1)
        v  => fieldsC(:,:,:,2)
        wC => fieldsC(:,:,:,3)
        w  => fieldsE(:,:,:,1)

        x => mesh(:,:,:,1)
        y => mesh(:,:,:,2)
        z => mesh(:,:,:,3)

        if (zmid < 0) then
            zmid = Lz * half
        end if

        ! initialize velocity fields
        call get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, yaw, u, v)
        wC= zero
        w = zero

        call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
        call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
        call message_min_max(1,"Bounds for w:", p_minval(minval(w)), p_maxval(maxval(w)))


        nullify(u,v,w,x,y,z)
        call message(0,"Velocity Field for Simulation 1 Initialized")
    else
        call gracefulExit("Only the Actuator disk simulation can be initialized like this. &
        & Check input file for HIT to ensure that it is restarted.",13)
    end if

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    use HIT_shear_parameters
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    if (simulationID == 1) then
        !allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))
        allocate(yplanes(nyplanes))
        allocate(xplanes(nxplanes))
        !allocate(zplanes(nzplanes))
        !xplanes = [300,400,500,600,700]
        yplanes = [nyADSim/2]
        xPlanes = [5*nxADSim/8] !800
        !zplanes = [128]
    end if
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
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: wTh_surf
    real(rkind) :: ThetaRef, Lx, Ly, Lz, uInflow = 1.d0
    integer :: iounit
    ! namelist /HIT_AD_interactINPUT/ Lx, Ly, Lz, uInflow

    wTh_surf = zero

    ! ioUnit = 11
    ! open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    ! read(unit=ioUnit, NML=HIT_AD_interactINPUT)
    ! close(ioUnit)

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, uInflow = 1.d0
    ! integer :: iounit
    ! namelist /HIT_AD_interactINPUT/ Lx, Ly, Lz, uInflow

    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one


    ! ioUnit = 11
    ! open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    ! read(unit=ioUnit, NML=HIT_AD_interactINPUT)
    ! close(ioUnit)

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, uInflow = 1.d0
    ! integer :: iounit

    ! namelist /HIT_AD_interactINPUT/ Lx, Ly, Lz, uInflow

    ! ioUnit = 11
    ! open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    ! read(unit=ioUnit, NML=HIT_AD_interactINPUT)
    ! close(ioUnit)

    Tref = 0.d0

    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
    character(len=*),                intent(in)    :: inputfile
    integer, parameter :: nprobes = 1

    ! IMPORTANT : Convention is to allocate probe_locs(3,nprobes)
    ! Example: If you have at least 3 probes:
    ! probe_locs(1,3) : x -location of the third probe
    ! probe_locs(2,3) : y -location of the third probe
    ! probe_locs(3,3) : z -location of the third probe


    ! Add probes here if needed
    ! Example code: The following allocates 2 probes at (0.1,0.1,0.1) and
    ! (0.2,0.2,0.2)
    allocate(probe_locs(3,nprobes))

    ! Probe 1
    probe_locs(1,1) = 0.1d0;
    probe_locs(2,1) = 3.141592653589d0;
    probe_locs(3,1) = 3.141592653589d0;

    ! add more probes below...
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


subroutine setScalar_source(decompC, inputfile, mesh, scalar_id, scalarSource)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    use HIT_shear_parameters, only: Sfunc
    use constants, only: pi
    use exits, only: message
    use reductions, only: p_sum

    type(decomp_info),                                          intent(in)    :: decompC
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarSource
    real(rkind), dimension(:,:,:), allocatable :: r, lambda, tmp
    real(rkind), dimension(:,:,:), pointer :: x, y, z
    real(rkind) :: xc = pi, yc = pi, zc = pi, rin = 0.75d0, rout = 1.25d0, delta_r = 0.22d0
    real(rkind) :: smear_x = 2.5d0, delta
    real(rkind) :: sumVal

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)


    allocate(r(size(x,1),size(x,2),size(x,3)))
    allocate(lambda(size(x,1),size(x,2),size(x,3)))
    allocate(tmp(size(x,1),size(x,2),size(x,3)))

    r = sqrt((y - yc)**2 + (z - zc)**2)

    select case (scalar_id)
      case (1)
        tmp = (r - rout)/delta_r + 1
        call Sfunc(tmp, lambda)
        lambda = -lambda
      case (2)
        tmp = (r - rin)/delta_r
        call Sfunc(tmp, lambda)
        lambda = 1.d0 - lambda
    end select

    r = x - xc
    delta = (x(2,1,1) - x(1,1,1))*smear_x
    tmp = (1.d0/(delta*sqrt(2.d0*pi)))*exp(-0.5d0*(r**2)/(delta**2))
    scalarSource = tmp*lambda
    sumVal = p_sum(sum(scalarSource))*((x(2,1,1) - x(1,1,1))**3)

    call message(2,"Scalar source initialized with domain integrated value", sumVal)
    deallocate(r, lambda, tmp)

end subroutine
