module streamtube_mod
   use mpi
   use exits, only: message, gracefulExit
   use constants, only: one, two, zero, half
   use kind_parameters, only: rkind, clen
   use PadeDerOps, only: Pade6stagg
   use spectralMod, only: spectral
   use decomp_2d
   use decomp_2d_io

   implicit none

   integer :: nx = 0, ny = 0, nz = 0
   integer :: prow = 0, pcol = 0
   integer :: NumericalSchemeVert = 1

   real(rkind) :: Lx = one, Ly = one, Lz = one
   real(rkind) :: dx = one, dy = one, dz = one

   logical :: PeriodicInZ = .false.
   integer :: botWall = 3, topWall = 2
   integer :: uBC_bottom, uBC_top
   integer :: vBC_bottom, vBC_top
   integer :: wBC_bottom, wBC_top

   type(spectral), target :: spectC, spectE
   type(decomp_info) :: gpC, gpE
   type(decomp_info), pointer :: sp_gpC, sp_gpE
   type(Pade6stagg) :: Pade6opZ

   ! Local x-pencil mesh: mesh(:,:,:,1:3) = x,y,z.
   real(rkind), dimension(:,:,:,:), allocatable, target :: mesh

   ! Buffers for derivatives/interpolation.
   real(rkind), dimension(:,:,:), allocatable :: rbuffxC
   real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
   complex(rkind), dimension(:,:,:), allocatable :: cbuffyC
   complex(rkind), dimension(:,:,:,:), allocatable, target :: cbuffzC

   ! Velocity-gradient pointers.
   real(rkind), dimension(:,:,:), pointer :: dudx => null(), dudy => null(), dudz => null()
   real(rkind), dimension(:,:,:), pointer :: dvdx => null(), dvdy => null(), dvdz => null()
   real(rkind), dimension(:,:,:), pointer :: dwdx => null(), dwdy => null(), dwdz => null()

   ! Numerical controls.
   real(rkind) :: umin_streamtube = 1.0e-12_rkind
   real(rkind) :: invalid_value = huge(one)
   integer :: streamtube_cutcell_nsamp_y = 11
   integer :: streamtube_cutcell_nsamp_z = 11

contains

!===============================================================================
! Initialization / finalization
!===============================================================================

   subroutine initialize_streamtube_module(nx_in, ny_in, nz_in, Lx_in, Ly_in, Lz_in, prow_in, pcol_in)
      implicit none
      integer, intent(in) :: nx_in, ny_in, nz_in
      integer, intent(in) :: prow_in, pcol_in
      real(rkind), intent(in) :: Lx_in, Ly_in, Lz_in

      logical :: periodicbcs(3)

      nx = nx_in; ny = ny_in; nz = nz_in
      Lx = Lx_in; Ly = Ly_in; Lz = Lz_in
      prow = prow_in; pcol = pcol_in

      dx = Lx / real(nx, rkind)
      dy = Ly / real(ny, rkind)
      dz = Lz / real(nz, rkind)

      periodicbcs(1) = .true.
      periodicbcs(2) = .true.
      periodicbcs(3) = PeriodicInZ

      call decomp_2d_init(nx, ny, nz, prow, pcol, periodicbcs)
      call get_decomp_info(gpC)
      call decomp_info_init(nx, ny, nz + 1, gpE)

      call message(0, 'Initializing streamtube spectral operators ...')
      call spectC%init('x', nx, ny, nz, dx, dy, dz, 'FOUR', '2/3rd', &
                       dimTransform=2, fixOddball=.false., init_periodicInZ=PeriodicInZ)
      call spectE%init('x', nx, ny, nz + 1, dx, dy, dz, 'FOUR', '2/3rd', &
                       dimTransform=2, fixOddball=.false., init_periodicInZ=PeriodicInZ)

      sp_gpC => spectC%spectdecomp
      sp_gpE => spectE%spectdecomp

      call allocate_streamtube_memory()
      call create_local_mesh()

      call Pade6opZ%init(gpC, sp_gpC, gpE, sp_gpE, dz, NumericalSchemeVert, PeriodicInZ, spectC)
      call get_boundary_conditions_stencil()
      call associate_gradient_pointers()

      call message(0, 'Streamtube module initialized.')
   end subroutine initialize_streamtube_module


   subroutine allocate_streamtube_memory()
      implicit none

      call message(0, 'Allocating streamtube work arrays ...')

      allocate(mesh(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3), 3))
      allocate(duidxj(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3), 9))
      allocate(cbuffyC(sp_gpC%ysz(1), sp_gpC%ysz(2), sp_gpC%ysz(3)))
      allocate(cbuffzC(sp_gpC%zsz(1), sp_gpC%zsz(2), sp_gpC%zsz(3), 2))

      mesh = zero
      duidxj = zero
      cbuffyC = cmplx(zero, zero, kind=rkind)
      cbuffzC = cmplx(zero, zero, kind=rkind)
   end subroutine allocate_streamtube_memory


   subroutine create_local_mesh()
      implicit none
      integer :: i, j, k
      integer :: ix1, iy1, iz1

      ix1 = gpC%xst(1)
      iy1 = gpC%xst(2)
      iz1 = gpC%xst(3)

      do k = 1, size(mesh,3)
         do j = 1, size(mesh,2)
            do i = 1, size(mesh,1)
               mesh(i,j,k,1) = real(ix1 + i - 2, rkind) * dx
               mesh(i,j,k,2) = real(iy1 + j - 2, rkind) * dy
               mesh(i,j,k,3) = real(iz1 + k - 2, rkind) * dz
            end do
         end do
      end do

      ! Cell-center vertical coordinate
      mesh(:,:,:,3) = mesh(:,:,:,3) + half * dz
   end subroutine create_local_mesh


   subroutine associate_gradient_pointers()
      implicit none

      dudx => duidxj(:,:,:,1)
      dudy => duidxj(:,:,:,2)
      dudz => duidxj(:,:,:,3)

      dvdx => duidxj(:,:,:,4)
      dvdy => duidxj(:,:,:,5)
      dvdz => duidxj(:,:,:,6)

      dwdx => duidxj(:,:,:,7)
      dwdy => duidxj(:,:,:,8)
      dwdz => duidxj(:,:,:,9)
   end subroutine associate_gradient_pointers


   subroutine release_streamtube_module()
      implicit none

      if (allocated(mesh)) deallocate(mesh)
      if (allocated(duidxj)) deallocate(duidxj)
      if (allocated(cbuffyC)) deallocate(cbuffyC)
      if (allocated(cbuffzC)) deallocate(cbuffzC)
      if (allocated(rbuffxC)) deallocate(rbuffxC)

      nullify(dudx, dudy, dudz)
      nullify(dvdx, dvdy, dvdz)
      nullify(dwdx, dwdy, dwdz)

      call spectC%destroy()
      call spectE%destroy()
      call Pade6opZ%destroy()
      call decomp_info_finalize(gpC)
      call decomp_info_finalize(gpE)
      call decomp_2d_finalize()
   end subroutine release_streamtube_module

!===============================================================================
! Main user-facing driver
!===============================================================================

   subroutine march_streamtube_xpencil(u, v, w, x_march, y0, z0, tube_y, tube_z)
      implicit none

      real(rkind), intent(in) :: u(:,:,:), v(:,:,:), w(:,:,:)
      real(rkind), intent(in) :: x_march(:)
      real(rkind), intent(in) :: y0(:), z0(:)
      real(rkind), intent(out) :: tube_y(:,:), tube_z(:,:)

      integer :: nxm, np, n
      real(rkind) :: dxm
      real(rkind), allocatable :: y_old(:), z_old(:)
      real(rkind), allocatable :: y_new(:), z_new(:)

      nxm = size(x_march)
      np = size(y0)

      if (size(z0) /= np) call gracefulExit('march_streamtube_xpencil: z0 size mismatch.', 201)
      if (size(tube_y,1) /= nxm .or. size(tube_y,2) /= np) call gracefulExit('march_streamtube_xpencil: tube_y size mismatch.', 202)
      if (size(tube_z,1) /= nxm .or. size(tube_z,2) /= np) call gracefulExit('march_streamtube_xpencil: tube_z size mismatch.', 203)

      allocate(y_old(np), z_old(np), y_new(np), z_new(np))

      call message(0, 'Computing velocity gradients for streamtube interpolation ...')
      call compute_velocity_gradients(u, v, w)

      y_old = y0
      z_old = z0
      tube_y(1,:) = y_old
      tube_z(1,:) = z_old

      do n = 1, nxm - 1
         call message(1, 'At x station: ', x_march(n))
         dxm = x_march(n+1) - x_march(n)

         call rk4_step_all_points(x_march(n), dxm, y_old, z_old, y_new, z_new, u, v, w)

         ! y_new/z_new are identical on all ranks because every RK-stage RHS is
         ! globally reduced with MPI_Allreduce. Thus a point can leave one rank's
         ! y-z pencil and be picked up by another rank at the next station.
         y_old = y_new
         z_old = z_new

         tube_y(n+1,:) = y_old
         tube_z(n+1,:) = z_old
      end do

      deallocate(y_old, z_old, y_new, z_new)
   end subroutine march_streamtube_xpencil

!===============================================================================
! RK4 integration in x
!===============================================================================

   subroutine rk4_step_all_points(x0, dxm, y0, z0, y1, z1, u, v, w)
      implicit none

      real(rkind), intent(in) :: x0, dxm
      real(rkind), intent(in) :: y0(:), z0(:)
      real(rkind), intent(out) :: y1(:), z1(:)
      real(rkind), intent(in) :: u(:,:,:), v(:,:,:), w(:,:,:)

      integer :: np, p
      real(rkind), allocatable :: k1y(:), k1z(:), k2y(:), k2z(:), k3y(:), k3z(:), k4y(:), k4z(:)
      real(rkind), allocatable :: yt(:), zt(:)
      logical, allocatable :: ok1(:), ok2(:), ok3(:), ok4(:)

      np = size(y0)

      allocate(k1y(np), k1z(np), k2y(np), k2z(np), k3y(np), k3z(np), k4y(np), k4z(np))
      allocate(yt(np), zt(np))
      allocate(ok1(np), ok2(np), ok3(np), ok4(np))

      call rhs_all_points(x0, y0, z0, k1y, k1z, ok1, u, v, w)

      yt = y0 + half * dxm * k1y
      zt = z0 + half * dxm * k1z
      call rhs_all_points(x0 + half*dxm, yt, zt, k2y, k2z, ok2, u, v, w)

      yt = y0 + half * dxm * k2y
      zt = z0 + half * dxm * k2z
      call rhs_all_points(x0 + half*dxm, yt, zt, k3y, k3z, ok3, u, v, w)

      yt = y0 + dxm * k3y
      zt = z0 + dxm * k3z
      call rhs_all_points(x0 + dxm, yt, zt, k4y, k4z, ok4, u, v, w)

      do p = 1, np
         if (ok1(p) .and. ok2(p) .and. ok3(p) .and. ok4(p)) then
            y1(p) = y0(p) + dxm * (k1y(p) + two*k2y(p) + two*k3y(p) + k4y(p)) / 6.0_rkind
            z1(p) = z0(p) + dxm * (k1z(p) + two*k2z(p) + two*k3z(p) + k4z(p)) / 6.0_rkind
         else
            y1(p) = invalid_value
            z1(p) = invalid_value
         end if
      end do

      deallocate(k1y, k1z, k2y, k2z, k3y, k3z, k4y, k4z)
      deallocate(yt, zt)
      deallocate(ok1, ok2, ok3, ok4)
   end subroutine rk4_step_all_points


   subroutine rhs_all_points(xp, yp, zp, rhs_y, rhs_z, ok, u, v, w)
      implicit none

      real(rkind), intent(in) :: xp
      real(rkind), intent(in) :: yp(:), zp(:)
      real(rkind), intent(out) :: rhs_y(:), rhs_z(:)
      logical, intent(out) :: ok(:)
      real(rkind), intent(in) :: u(:,:,:), v(:,:,:), w(:,:,:)

      integer :: np, p, ierr
      real(rkind), allocatable :: rhs_y_local(:), rhs_z_local(:)
      real(rkind), allocatable :: rhs_y_global(:), rhs_z_global(:)
      real(rkind), allocatable :: found_local(:), found_global(:)
      real(rkind) :: up, vp, wp
      logical :: found

      np = size(yp)

      allocate(rhs_y_local(np), rhs_z_local(np))
      allocate(rhs_y_global(np), rhs_z_global(np))
      allocate(found_local(np), found_global(np))

      rhs_y_local = zero
      rhs_z_local = zero
      found_local = zero

      do p = 1, np
         if (.not. is_valid_point(yp(p), zp(p))) cycle

         call interp_velocity_taylor(xp, yp(p), zp(p), u, v, w, up, vp, wp, found)

         if (found) then
            if (abs(up) > umin_streamtube) then
               rhs_y_local(p) = vp / up
               rhs_z_local(p) = wp / up
               found_local(p) = one
            end if
         end if
      end do

      call mpi_allreduce(rhs_y_local, rhs_y_global, np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call mpi_allreduce(rhs_z_local, rhs_z_global, np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call mpi_allreduce(found_local, found_global, np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      do p = 1, np
         if (found_global(p) > half) then
            rhs_y(p) = rhs_y_global(p) / found_global(p)
            rhs_z(p) = rhs_z_global(p) / found_global(p)
            ok(p) = .true.
         else
            rhs_y(p) = zero
            rhs_z(p) = zero
            ok(p) = .false.
         end if
      end do

      deallocate(rhs_y_local, rhs_z_local)
      deallocate(rhs_y_global, rhs_z_global)
      deallocate(found_local, found_global)
   end subroutine rhs_all_points

!===============================================================================
! Velocity interpolation
!===============================================================================

   subroutine interp_velocity_taylor(xp, yp, zp, u, v, w, up, vp, wp, found)
      implicit none

      real(rkind), intent(in) :: xp, yp, zp
      real(rkind), intent(in) :: u(:,:,:), v(:,:,:), w(:,:,:)
      real(rkind), intent(out) :: up, vp, wp
      logical, intent(out) :: found

      integer :: i, j, k
      real(rkind) :: dxp, dyp, dzp

      up = zero
      vp = zero
      wp = zero
      found = .false.

      ! Locate the owning cell once, then reuse the same local index and Taylor
      ! displacement for u, v, and w.  The passed fields are all cell-centered.
      call find_nearest_owned_cell(xp, yp, zp, i, j, k, dxp, dyp, dzp, found)
      if (.not. found) return

      call interp_scalar_taylor_1st(i, j, k, dxp, dyp, dzp, u, dudx, dudy, dudz, up)
      call interp_scalar_taylor_1st(i, j, k, dxp, dyp, dzp, v, dvdx, dvdy, dvdz, vp)
      call interp_scalar_taylor_1st(i, j, k, dxp, dyp, dzp, w, dwdx, dwdy, dwdz, wp)

      found = .true.
   end subroutine interp_velocity_taylor


   subroutine interp_scalar_taylor_1st(i, j, k, dxp, dyp, dzp, f, dfdx, dfdy, dfdz, fp)
      implicit none

      integer, intent(in) :: i, j, k
      real(rkind), intent(in) :: dxp, dyp, dzp
      real(rkind), intent(in) :: f(:,:,:)
      real(rkind), intent(in) :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
      real(rkind), intent(out) :: fp

      fp = f(i,j,k) + dxp*dfdx(i,j,k) + dyp*dfdy(i,j,k) + dzp*dfdz(i,j,k)
   end subroutine interp_scalar_taylor_1st


   subroutine find_nearest_owned_cell(xp, yp, zp, i, j, k, dxp, dyp, dzp, found)
      implicit none

      real(rkind), intent(in) :: xp, yp, zp
      integer, intent(out) :: i, j, k
      real(rkind), intent(out) :: dxp, dyp, dzp
      logical, intent(out) :: found

      integer :: ig, jg, kg
      real(rkind) :: xc, yc, zc

      found = .false.
      i = -1; j = -1; k = -1
      dxp = zero; dyp = zero; dzp = zero

      if (zp < zero .or. zp > Lz) then
         found = .false.
         return
      end if

      ! x and y are periodic; z is wall-bounded and remains clamped.
      call nearest_global_index_uniform_periodic(xp, dx, nx, ig)
      call nearest_global_index_uniform_periodic(yp, dy, ny, jg)
      call nearest_global_index_uniform_zcell(zp, dz, nz, kg)

      if (ig < gpC%xst(1) .or. ig > gpC%xen(1)) return
      if (jg < gpC%xst(2) .or. jg > gpC%xen(2)) return
      if (kg < gpC%xst(3) .or. kg > gpC%xen(3)) return

      i = ig - gpC%xst(1) + 1
      j = jg - gpC%xst(2) + 1
      k = kg - gpC%xst(3) + 1

      if (i < 1 .or. i > size(mesh,1)) return
      if (j < 1 .or. j > size(mesh,2)) return
      if (k < 1 .or. k > size(mesh,3)) return

      xc = mesh(i,j,k,1)
      yc = mesh(i,j,k,2)
      zc = mesh(i,j,k,3)

      ! Use the nearest periodic image for the Taylor displacement.  Wrapping the
      ! index alone is not sufficient: xp-xc or yp-yc can otherwise be O(L).
      dxp = periodic_delta(xp - xc, Lx)
      dyp = periodic_delta(yp - yc, Ly)
      dzp = zp - zc

      found = .true.
   end subroutine find_nearest_owned_cell


   subroutine nearest_global_index_uniform_periodic(xp, dd, nn, ig)
      implicit none
      real(rkind), intent(in) :: xp, dd
      integer, intent(in) :: nn
      integer, intent(out) :: ig

      ig = modulo(nint(xp / dd), nn) + 1
   end subroutine nearest_global_index_uniform_periodic


   subroutine nearest_global_index_uniform_zcell(zp, dd, nn, kg)
      implicit none
      real(rkind), intent(in) :: zp, dd
      integer, intent(in) :: nn
      integer, intent(out) :: kg

      ! z cell centers are at (k-1)*dz + dz/2.
      kg = nint((zp - half*dd) / dd) + 1
      kg = max(1, min(nn, kg))
   end subroutine nearest_global_index_uniform_zcell


   pure real(rkind) function periodic_delta(delta, period)
      implicit none
      real(rkind), intent(in) :: delta, period

      if (period > zero) then
         periodic_delta = delta - period * anint(delta / period)
      else
         periodic_delta = delta
      end if
   end function periodic_delta


   logical function is_valid_point(y, z)
      implicit none
      real(rkind), intent(in) :: y, z

      is_valid_point = .true.
      if (abs(y) >= half*invalid_value) is_valid_point = .false.
      if (abs(z) >= half*invalid_value) is_valid_point = .false.
   end function is_valid_point

!===============================================================================
! Velocity gradients used by Taylor interpolation
!===============================================================================

   subroutine compute_velocity_gradients(u, v, w)
      implicit none
      real(rkind), intent(in) :: u(:,:,:), v(:,:,:), w(:,:,:)

      duidxj = zero

      call ddx_R2R(u, dudx)
      call ddy_R2R(u, dudy)
      call ddz_R2R(u, dudz, uBC_bottom, uBC_top)

      call ddx_R2R(v, dvdx)
      call ddy_R2R(v, dvdy)
      call ddz_R2R(v, dvdz, vBC_bottom, vBC_top)

      call ddx_R2R(w, dwdx)
      call ddy_R2R(w, dwdy)
      call ddz_R2R(w, dwdz, wBC_bottom, wBC_top)
   end subroutine compute_velocity_gradients


   subroutine ddx_R2R(f, dfdx)
      implicit none
      real(rkind), dimension(:,:,:), intent(in) :: f
      real(rkind), dimension(:,:,:), intent(out) :: dfdx

      call spectC%fft(f, cbuffyC)
      call spectC%mtimes_ik1_ip(cbuffyC)
      call spectC%dealias(cbuffyC)
      call spectC%ifft(cbuffyC, dfdx)
   end subroutine ddx_R2R


   subroutine ddy_R2R(f, dfdy)
      implicit none
      real(rkind), dimension(:,:,:), intent(in) :: f
      real(rkind), dimension(:,:,:), intent(out) :: dfdy

      call spectC%fft(f, cbuffyC)
      call spectC%mtimes_ik2_ip(cbuffyC)
      call spectC%dealias(cbuffyC)
      call spectC%ifft(cbuffyC, dfdy)
   end subroutine ddy_R2R


   subroutine ddz_R2R(f, dfdz, n1, n2)
      implicit none
      real(rkind), dimension(:,:,:), intent(in) :: f
      real(rkind), dimension(:,:,:), intent(out) :: dfdz
      integer, intent(in) :: n1, n2

      call spectC%fft(f, cbuffyC)
      call transpose_y_to_z(cbuffyC, cbuffzC(:,:,:,1), sp_gpC)
      call Pade6opZ%ddz_C2C(cbuffzC(:,:,:,1), cbuffzC(:,:,:,2), n1, n2)
      call transpose_z_to_y(cbuffzC(:,:,:,2), cbuffyC, sp_gpC)
      call spectC%dealias(cbuffyC)
      call spectC%ifft(cbuffyC, dfdz)
   end subroutine ddz_R2R


   subroutine get_boundary_conditions_stencil()
      implicit none

      wBC_bottom = -1
      wBC_top = -1

      select case (botWall)
      case (1)
         call message(1, 'Bottom wall: no-slip wall')
         uBC_bottom = 0
         vBC_bottom = 0
         wBC_bottom = 1
      case (2)
         call message(1, 'Bottom wall: slip wall')
         uBC_bottom = 1
         vBC_bottom = 1
      case (3)
         call message(1, 'Bottom wall: wall model')
         uBC_bottom = 0
         vBC_bottom = 0
      case default
         call gracefulExit('Invalid choice for bottom wall BCs.', 301)
      end select

      select case (topWall)
      case (1)
         call message(1, 'Top wall: no-slip wall')
         uBC_top = 0
         vBC_top = 0
         wBC_top = 1
      case (2)
         call message(1, 'Top wall: slip wall')
         uBC_top = 1
         vBC_top = 1
      case (3)
         call message(1, 'Top wall: wall model')
         uBC_top = 0
         vBC_top = 0
      case default
         call gracefulExit('Invalid choice for top wall BCs.', 302)
      end select
   end subroutine get_boundary_conditions_stencil

!===============================================================================
! Output helpers
!===============================================================================

   subroutine export_streamtube_csv(filename, x_march, tube_y, tube_z)
      implicit none
      character(len=*), intent(in) :: filename
      real(rkind), intent(in) :: x_march(:)
      real(rkind), intent(in) :: tube_y(:,:), tube_z(:,:)

      integer :: unit, n, p, nxm, np

      nxm = size(x_march)
      np = size(tube_y, 2)

      open(newunit=unit, file=trim(filename), status='replace', action='write', form='formatted')

      write(unit, '(A)', advance='no') 'x'
      do p = 1, np
         write(unit, '(A,I0)', advance='no') ',y', p
         write(unit, '(A,I0)', advance='no') ',z', p
      end do
      write(unit, *)

      do n = 1, nxm
         write(unit, '(ES16.8)', advance='no') x_march(n)
         do p = 1, np
            write(unit, '(A,ES16.8,A,ES16.8)', advance='no') ',', tube_y(n,p), ',', tube_z(n,p)
         end do
         write(unit, *)
      end do

      close(unit)
   end subroutine export_streamtube_csv


!===============================================================================
! Optional geometry diagnostics
!===============================================================================

subroutine compute_contour_extents(tube_y, tube_z, valid, width_y, height_z)
implicit none

real(rkind), intent(in)  :: tube_y(:,:), tube_z(:,:)
logical, intent(in)      :: valid(:,:)
real(rkind), intent(out) :: width_y(:), height_z(:)

integer :: n

if (size(tube_z,1) /= size(tube_y,1)) call gracefulExit('compute_contour_extents: tube_z station mismatch.', 551)
if (size(tube_z,2) /= size(tube_y,2)) call gracefulExit('compute_contour_extents: tube_z point mismatch.', 552)
if (size(valid,1) /= size(tube_y,1)) call gracefulExit('compute_contour_extents: valid station mismatch.', 553)
if (size(valid,2) /= size(tube_y,2)) call gracefulExit('compute_contour_extents: valid point mismatch.', 554)

do n = 1, size(tube_y,1)

   if (all(valid(n,:))) then
      width_y(n) = maxval(tube_y(n,:), mask=valid(n,:)) - &
                     minval(tube_y(n,:), mask=valid(n,:))

      height_z(n) = maxval(tube_z(n,:), mask=valid(n,:)) - &
                     minval(tube_z(n,:), mask=valid(n,:))
   else
      width_y(n) = zero
      height_z(n) = zero
   end if

end do

end subroutine compute_contour_extents


subroutine compute_contour_area(tube_y, tube_z, valid, area)
implicit none

real(rkind), intent(in) :: tube_y(:,:), tube_z(:,:)
logical, intent(in) :: valid(:,:)
real(rkind), intent(out) :: area(:)

integer :: n, p, pp, np, nv
real(rkind) :: accum
real(rkind), allocatable :: yv(:), zv(:)

np = size(tube_y,2)

if (size(tube_z,1) /= size(tube_y,1)) call gracefulExit('compute_contour_area: tube_z station mismatch.', 541)
if (size(tube_z,2) /= np) call gracefulExit('compute_contour_area: tube_z point mismatch.', 542)
if (size(valid,1) /= size(tube_y,1)) call gracefulExit('compute_contour_area: valid station mismatch.', 543)
if (size(valid,2) /= np) call gracefulExit('compute_contour_area: valid point mismatch.', 544)
if (size(area) /= size(tube_y,1)) call gracefulExit('compute_contour_area: area size mismatch.', 545)

allocate(yv(np), zv(np))

do n = 1, size(tube_y,1)

   if (.not. all(valid(n,:))) then
      area(n) = zero
      cycle
   end if

   call compact_valid_contour(tube_y(n,:), tube_z(n,:), valid(n,:), yv, zv, nv)

   if (nv < 3) then
      area(n) = zero
      cycle
   end if

   accum = zero

   do p = 1, nv
      pp = p + 1
      if (pp > nv) pp = 1
      accum = accum + yv(p)*zv(pp) - zv(p)*yv(pp)
   end do

   area(n) = half * abs(accum)

end do

deallocate(yv, zv)

end subroutine compute_contour_area

   subroutine build_x_march(xa, xb, xarr)
      implicit none
      real(rkind), intent(in) :: xa, xb
      real(rkind), intent(out) :: xarr(:)
      integer :: n, nloc

      nloc = size(xarr)

      if (nloc < 2) then
         call gracefulExit('build_x_march: need at least 2 stations.', 570)
      end if

      do n = 1, nloc
         xarr(n) = xa + (xb - xa) * real(n - 1, rkind) / real(nloc - 1, rkind)
      end do
   end subroutine build_x_march


   subroutine count_initial_contour_points(filename, npts)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: npts

      integer :: unit, ios
      character(len=4096) :: line
      real(rkind) :: ytmp, ztmp
      logical :: ok

      npts = 0
      open(newunit=unit, file=trim(filename), status='old', action='read', form='formatted')

      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         call parse_yz_line(line, ytmp, ztmp, ok)
         if (ok) npts = npts + 1
      end do

      close(unit)
   end subroutine count_initial_contour_points


   subroutine read_initial_contour_csv(filename, yy, zz)
      implicit none
      character(len=*), intent(in) :: filename
      real(rkind), intent(out) :: yy(:), zz(:)

      integer :: unit, ios, p, npts
      character(len=4096) :: line
      real(rkind) :: ytmp, ztmp
      logical :: ok

      npts = size(yy)
      if (size(zz) /= npts) call gracefulExit('Internal contour array size mismatch.', 201)

      p = 0
      open(newunit=unit, file=trim(filename), status='old', action='read', form='formatted')

      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit

         call parse_yz_line(line, ytmp, ztmp, ok)
         if (.not. ok) cycle

         p = p + 1
         if (p > npts) then
            call gracefulExit('Initial contour file has more valid points than num_contour_points.', 202)
         end if

         yy(p) = ytmp
         zz(p) = ztmp
      end do

      close(unit)

      if (p /= npts) then
         call gracefulExit('Initial contour file has fewer valid points than expected.', 203)
      end if
   end subroutine read_initial_contour_csv


   subroutine parse_yz_line(line_in, yy, zz, ok)
      implicit none
      character(len=*), intent(in) :: line_in
      real(rkind), intent(out) :: yy, zz
      logical, intent(out) :: ok

      character(len=len(line_in)) :: line
      integer :: i, ios

      line = adjustl(line_in)

      ok = .false.
      yy = zero
      zz = zero

      if (len_trim(line) == 0) return
      if (line(1:1) == '#') return

      do i = 1, len(line)
         if (line(i:i) == ',') line(i:i) = ' '
      end do

      read(line, *, iostat=ios) yy, zz
      ok = (ios == 0)
   end subroutine parse_yz_line


   subroutine export_geometry_csv(filename, xarr, area_in, width_in, height_in, flow_rate_in)
      implicit none
      character(len=*), intent(in) :: filename
      real(rkind), intent(in) :: xarr(:)
      real(rkind), intent(in) :: area_in(:), width_in(:), height_in(:)
      real(rkind), intent(in) :: flow_rate_in(:)

      integer :: unit, n

      call message(0, 'Writing streamtube geometry to '//trim(filename))

      open(newunit=unit, file=trim(filename), status='replace', action='write', form='formatted')
      write(unit, '(A)') 'x,area,width_y,height_z,flow_rate'

      do n = 1, size(xarr)
         write(unit, '(ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') &
            xarr(n), ',', area_in(n), ',', width_in(n), ',', height_in(n), ',', flow_rate_in(n)
      end do

      close(unit)
   end subroutine export_geometry_csv

subroutine export_streamtube_stl(filename, x_march, tube_y, tube_z, valid)
   implicit none

   character(len=*), intent(in) :: filename
   real(rkind), intent(in) :: x_march(:)
   real(rkind), intent(in) :: tube_y(:,:), tube_z(:,:)
   logical, intent(in) :: valid(:,:)

   integer :: unit
   integer :: n, p, pp
   integer :: nxm, np
   real(rkind) :: v1(3), v2(3), v3(3), v4(3)
   logical :: ok

   nxm = size(x_march)
   np  = size(tube_y, 2)

   if (size(tube_y,1) /= nxm) call gracefulExit('export_streamtube_stl: tube_y station mismatch.', 401)
   if (size(tube_z,1) /= nxm) call gracefulExit('export_streamtube_stl: tube_z station mismatch.', 402)
   if (size(tube_z,2) /= np ) call gracefulExit('export_streamtube_stl: tube_z point mismatch.', 403)

   call message(0, 'Writing streamtube STL to '//trim(filename))

   open(newunit=unit, file=trim(filename), status='replace', action='write', form='formatted')

   write(unit,'(A)') 'solid streamtube'

   do n = 1, nxm - 1
      do p = 1, np

         pp = p + 1
         if (pp > np) pp = 1

         ok = all(valid(n,:)) .and. all(valid(n+1,:))

         if (.not. ok) cycle

         v1 = (/ x_march(n),   tube_y(n,  p ), tube_z(n,  p ) /)
         v2 = (/ x_march(n+1), tube_y(n+1,p ), tube_z(n+1,p ) /)
         v3 = (/ x_march(n+1), tube_y(n+1,pp), tube_z(n+1,pp) /)
         v4 = (/ x_march(n),   tube_y(n,  pp), tube_z(n,  pp) /)

         ! Quad split:
         ! Triangle 1: v1 -> v2 -> v3
         ! Triangle 2: v1 -> v3 -> v4
         call write_stl_triangle(unit, v1, v2, v3)
         call write_stl_triangle(unit, v1, v3, v4)

      end do
   end do

   write(unit,'(A)') 'endsolid streamtube'

   close(unit)

end subroutine export_streamtube_stl

subroutine write_stl_triangle(unit, a, b, c)
   implicit none

   integer, intent(in) :: unit
   real(rkind), intent(in) :: a(3), b(3), c(3)

   real(rkind) :: nvec(3)

   call triangle_normal(a, b, c, nvec)

   write(unit,'(A,3(1X,ES16.8))') '  facet normal', nvec(1), nvec(2), nvec(3)
   write(unit,'(A)')              '    outer loop'
   write(unit,'(A,3(1X,ES16.8))') '      vertex', a(1), a(2), a(3)
   write(unit,'(A,3(1X,ES16.8))') '      vertex', b(1), b(2), b(3)
   write(unit,'(A,3(1X,ES16.8))') '      vertex', c(1), c(2), c(3)
   write(unit,'(A)')              '    endloop'
   write(unit,'(A)')              '  endfacet'

end subroutine write_stl_triangle

subroutine triangle_normal(a, b, c, nvec)
   implicit none

   real(rkind), intent(in)  :: a(3), b(3), c(3)
   real(rkind), intent(out) :: nvec(3)

   real(rkind) :: ab(3), ac(3)
   real(rkind) :: mag

   ab = b - a
   ac = c - a

   nvec(1) = ab(2)*ac(3) - ab(3)*ac(2)
   nvec(2) = ab(3)*ac(1) - ab(1)*ac(3)
   nvec(3) = ab(1)*ac(2) - ab(2)*ac(1)

   mag = sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)

   if (mag > tiny(one)) then
      nvec = nvec / mag
   else
      nvec = zero
   end if

end subroutine triangle_normal

subroutine compute_streamtube_flow_rate(u, x_march, tube_y, tube_z, valid, flow_rate)
   implicit none

   real(rkind), intent(in)  :: u(:,:,:)
   real(rkind), intent(in)  :: x_march(:)
   real(rkind), intent(in)  :: tube_y(:,:), tube_z(:,:)
   logical, intent(in)      :: valid(:,:)
   real(rkind), intent(out) :: flow_rate(:)
   real(rkind), allocatable :: q_local_arr(:), q_global_arr(:)

   integer :: n, ierr
   integer :: nxm, np
   integer :: ig0, ig1
   integer :: il0, il1
   integer :: nv
   real(rkind) :: wx0, wx1, q0_local, q1_local
   logical :: have_i0, have_i1
   real(rkind), allocatable :: yv(:), zv(:)

   nxm = size(x_march)
   np  = size(tube_y, 2)

   if (size(tube_y,1) /= nxm) call gracefulExit('compute_streamtube_flow_rate: tube_y station mismatch.', 501)
   if (size(tube_z,1) /= nxm) call gracefulExit('compute_streamtube_flow_rate: tube_z station mismatch.', 502)
   if (size(tube_z,2) /= np ) call gracefulExit('compute_streamtube_flow_rate: tube_z point mismatch.', 503)
   if (size(valid,1) /= nxm .or. size(valid,2) /= np) call gracefulExit('compute_streamtube_flow_rate: valid size mismatch.', 504)
   if (size(flow_rate) /= nxm) call gracefulExit('compute_streamtube_flow_rate: flow_rate size mismatch.', 505)

   allocate(yv(np), zv(np))
   allocate(q_local_arr(nxm), q_global_arr(nxm))
   q_local_arr = zero
   q_global_arr = zero

   flow_rate = zero

   do n = 1, nxm

      ! For a closed streamtube contour, do not compact a partially invalid
      ! polygon. That would close the contour with artificial chords.
      if (.not. all(valid(n,:))) then
         q_local_arr(n) = zero
         cycle
      end if

      call compact_valid_contour(tube_y(n,:), tube_z(n,:), valid(n,:), yv, zv, nv)

      if (nv < 3) then
         flow_rate(n) = zero
         cycle
      end if

      call bracket_x_indices(x_march(n), ig0, ig1, wx0, wx1)

      call global_x_to_local(ig0, il0, have_i0)
      call global_x_to_local(ig1, il1, have_i1)

      q0_local = zero
      q1_local = zero

      if (have_i0) then
         call integrate_u_on_xplane_cutcells(u, il0, yv, zv, nv, q0_local)
      end if

      if (have_i1) then
         call integrate_u_on_xplane_cutcells(u, il1, yv, zv, nv, q1_local)
      end if

      q_local_arr(n) = wx0*q0_local + wx1*q1_local

   end do

   call mpi_allreduce(q_local_arr, q_global_arr, nxm, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, MPI_COMM_WORLD, ierr)

   flow_rate = q_global_arr

   deallocate(q_local_arr, q_global_arr)
   deallocate(yv, zv)

end subroutine compute_streamtube_flow_rate

logical function valid_mask_has_gap(valid)
   implicit none

   logical, intent(in) :: valid(:)

   integer :: p, np
   integer :: ntrans

   np = size(valid)
   ntrans = 0

   do p = 1, np
      if (valid(p) .neqv. valid(merge(p+1,1,p<np))) then
         ntrans = ntrans + 1
      end if
   end do

   ! For a circular list:
   ! ntrans = 0 means all valid or all invalid.
   ! ntrans = 2 means one contiguous valid block and one invalid block.
   ! ntrans > 2 means multiple gaps.
   valid_mask_has_gap = ntrans > 2

end function valid_mask_has_gap

subroutine bracket_x_indices(xp_in, ig0, ig1, wx0, wx1)
   implicit none

   real(rkind), intent(in)  :: xp_in
   integer, intent(out)    :: ig0, ig1
   real(rkind), intent(out) :: wx0, wx1

   real(rkind) :: xp
   real(rkind) :: s
   integer :: i0

   ! x is periodic in this module.
   xp = modulo(xp_in, Lx)

   s  = xp / dx
   i0 = floor(s)

   wx1 = s - real(i0, rkind)
   wx0 = one - wx1

   ig0 = modulo(i0, nx) + 1
   ig1 = modulo(i0 + 1, nx) + 1

   ! Exact or near-exact grid plane.
   if (abs(wx1) < 10.0_rkind*epsilon(one)) then
      wx0 = one
      wx1 = zero
      ig1 = ig0
   end if

end subroutine bracket_x_indices

subroutine global_x_to_local(ig, iloc, have_i)
   implicit none

   integer, intent(in)  :: ig
   integer, intent(out) :: iloc
   logical, intent(out) :: have_i

   have_i = .false.
   iloc = -1

   if (ig < gpC%xst(1) .or. ig > gpC%xen(1)) return

   iloc = ig - gpC%xst(1) + 1

   if (iloc < 1 .or. iloc > size(mesh,1)) then
      have_i = .false.
   else
      have_i = .true.
   end if

end subroutine global_x_to_local

subroutine integrate_u_on_xplane_cutcells(u, iloc, yy, zz, nv, q_local)
   implicit none

   real(rkind), intent(in)  :: u(:,:,:)
   integer, intent(in)      :: iloc
   real(rkind), intent(in)  :: yy(:), zz(:)
   integer, intent(in)      :: nv
   real(rkind), intent(out) :: q_local

   integer :: j, k
   real(rkind) :: yc, zc
   real(rkind) :: ylo, yhi, zlo, zhi
   real(rkind) :: frac
   real(rkind) :: cell_area
   real(rkind) :: ymin_poly, ymax_poly
   real(rkind) :: zmin_poly, zmax_poly

   q_local = zero
   cell_area = dy * dz

   call valid_contour_bounds(yy, zz, nv, ymin_poly, ymax_poly, zmin_poly, zmax_poly)

   if (ymax_poly <= ymin_poly) return
   if (zmax_poly <= zmin_poly) return

   do k = 1, size(u,3)
      do j = 1, size(u,2)

         yc = mesh(iloc,j,k,2)
         zc = mesh(iloc,j,k,3)

         ylo = yc - half*dy
         yhi = yc + half*dy
         zlo = zc - half*dz
         zhi = zc + half*dz

         if (yhi < ymin_poly) cycle
         if (ylo > ymax_poly) cycle
         if (zhi < zmin_poly) cycle
         if (zlo > zmax_poly) cycle

         call cell_area_fraction_by_sampling(yy, zz, nv, ylo, yhi, zlo, zhi, frac)

         if (frac > zero) then
            q_local = q_local + u(iloc,j,k) * frac * cell_area
         end if

      end do
   end do

end subroutine integrate_u_on_xplane_cutcells

subroutine valid_contour_bounds(yy, zz, nv, ymin_poly, ymax_poly, zmin_poly, zmax_poly)
   implicit none

   real(rkind), intent(in)  :: yy(:), zz(:)
   integer, intent(in)      :: nv
   real(rkind), intent(out) :: ymin_poly, ymax_poly
   real(rkind), intent(out) :: zmin_poly, zmax_poly

   integer :: p

   if (nv < 1) then
      ymin_poly = zero
      ymax_poly = zero
      zmin_poly = zero
      zmax_poly = zero
      return
   end if

   ymin_poly = yy(1)
   ymax_poly = yy(1)
   zmin_poly = zz(1)
   zmax_poly = zz(1)

   do p = 2, nv
      ymin_poly = min(ymin_poly, yy(p))
      ymax_poly = max(ymax_poly, yy(p))
      zmin_poly = min(zmin_poly, zz(p))
      zmax_poly = max(zmax_poly, zz(p))
   end do

end subroutine valid_contour_bounds

subroutine cell_area_fraction_by_sampling(yy, zz, nv, ylo, yhi, zlo, zhi, frac)
   implicit none

   real(rkind), intent(in)  :: yy(:), zz(:)
   integer, intent(in)      :: nv
   real(rkind), intent(in)  :: ylo, yhi, zlo, zhi
   real(rkind), intent(out) :: frac

   logical :: c1, c2, c3, c4
   integer :: ninside

   c1 = point_in_polygon(yy, zz, nv, ylo, zlo)
   c2 = point_in_polygon(yy, zz, nv, yhi, zlo)
   c3 = point_in_polygon(yy, zz, nv, yhi, zhi)
   c4 = point_in_polygon(yy, zz, nv, ylo, zhi)

   ninside = 0
   if (c1) ninside = ninside + 1
   if (c2) ninside = ninside + 1
   if (c3) ninside = ninside + 1
   if (c4) ninside = ninside + 1

   if (ninside == 4) then
      frac = one
   else
      call partial_cell_area_fraction_sampling(yy, zz, nv, ylo, yhi, zlo, zhi, frac)
   end if

end subroutine cell_area_fraction_by_sampling

subroutine partial_cell_area_fraction_sampling(yy, zz, nv, ylo, yhi, zlo, zhi, frac)
   implicit none

   real(rkind), intent(in)  :: yy(:), zz(:)
   integer, intent(in)      :: nv
   real(rkind), intent(in)  :: ylo, yhi, zlo, zhi
   real(rkind), intent(out) :: frac

   integer :: iy, iz
   integer :: ny_samp, nz_samp
   integer :: ninside, ntotal
   real(rkind) :: yp, zp
   real(rkind) :: sy, sz

   ny_samp = max(3, streamtube_cutcell_nsamp_y)
   nz_samp = max(3, streamtube_cutcell_nsamp_z)

   ninside = 0
   ntotal  = ny_samp * nz_samp

   do iz = 1, nz_samp

      sz = real(iz - 1, rkind) / real(nz_samp - 1, rkind)
      zp = zlo + sz * (zhi - zlo)

      do iy = 1, ny_samp

         sy = real(iy - 1, rkind) / real(ny_samp - 1, rkind)
         yp = ylo + sy * (yhi - ylo)

         if (point_in_polygon(yy, zz, nv, yp, zp)) then
            ninside = ninside + 1
         end if

      end do
   end do

   frac = real(ninside, rkind) / real(ntotal, rkind)

end subroutine partial_cell_area_fraction_sampling

subroutine build_streamtube_valid_mask(tube_y, tube_z, valid)
   implicit none

   real(rkind), intent(in)  :: tube_y(:,:), tube_z(:,:)
   logical, intent(out)     :: valid(:,:)

   integer :: n, p

   if (size(valid,1) /= size(tube_y,1)) call gracefulExit('build_streamtube_valid_mask: station mismatch.', 531)
   if (size(valid,2) /= size(tube_y,2)) call gracefulExit('build_streamtube_valid_mask: point mismatch.', 532)
   if (size(tube_z,1) /= size(tube_y,1)) call gracefulExit('build_streamtube_valid_mask: tube_z station mismatch.', 533)
   if (size(tube_z,2) /= size(tube_y,2)) call gracefulExit('build_streamtube_valid_mask: tube_z point mismatch.', 534)

   do n = 1, size(tube_y,1)
      do p = 1, size(tube_y,2)
         valid(n,p) = is_valid_point(tube_y(n,p), tube_z(n,p))
      end do
   end do

end subroutine build_streamtube_valid_mask

logical function point_in_polygon(yy, zz, nv, yp, zp)
!===============================================================================
! Point-in-polygon classification
!
! Implements the dual-perspective point-in-polygon method of Ali and Guaily.
! The method classifies a point by identifying the closest polygon vertex and
! combining the inside/outside perspectives of the two parent edges attached to
! that vertex. Points lying exactly on polygon edges are treated as inside for
! the present streamtube area-fraction calculation.
!
! Reference:
!   Ali, K. M. and Guaily, A. (2020). "Dual perspective method for solving the
!   point in a polygon problem." arXiv:2012.05001.
!
! Notes:
!   - The polygon vertices must be ordered consistently around the contour.
!   - The polygon must be non-self-intersecting.
!   - This implementation assumes compacted valid vertices are passed in; invalid
!     streamtube points should be removed before calling point_in_polygon.
!===============================================================================
   implicit none

   real(rkind), intent(in) :: yy(:), zz(:)
   integer, intent(in) :: nv
   real(rkind), intent(in) :: yp, zp

   integer :: ic, iprev, inext
   integer :: s1, s2
   logical :: L1, L2
   logical :: dismiss1, dismiss2
   real(rkind) :: c
   real(rkind) :: A1y, A1z, A2y, A2z
   real(rkind) :: tol
   real(rkind) :: area_signed
   logical :: ccw

   point_in_polygon = .false.

   if (nv < 3) return

   tol = point_polygon_tolerance(yy, zz, nv)

   call polygon_signed_area(yy, zz, nv, area_signed)

   if (abs(area_signed) <= tol) then
      point_in_polygon = .false.
      return
   end if

   ccw = area_signed > zero

   call closest_vertex(yy, zz, nv, yp, zp, ic)

   if (ic < 1) then
      point_in_polygon = .false.
      return
   end if

   if ((yp - yy(ic))**2 + (zp - zz(ic))**2 <= tol**2) then
      point_in_polygon = .true.
      return
   end if

   iprev = ic - 1
   if (iprev < 1) iprev = nv

   inext = ic + 1
   if (inext > nv) inext = 1

   A1y = yy(inext) - yy(ic)
   A1z = zz(inext) - zz(ic)

   A2y = yy(iprev) - yy(ic)
   A2z = zz(iprev) - zz(ic)

   c = A1y*A2z - A1z*A2y

   call edge_perspective_sign(yy(ic), zz(ic), yy(inext), zz(inext), &
                              yp, zp, ccw, tol, s1, L1, dismiss1)

   call edge_perspective_sign(yy(iprev), zz(iprev), yy(ic), zz(ic), &
                              yp, zp, ccw, tol, s2, L2, dismiss2)

   if ((s1 == 0 .and. L1) .or. (s2 == 0 .and. L2)) then
      point_in_polygon = .true.
      return
   end if

   if (dismiss1 .and. .not. dismiss2) then
      point_in_polygon = s2 <= 0
      return
   end if

   if (dismiss2 .and. .not. dismiss1) then
      point_in_polygon = s1 <= 0
      return
   end if

   if (dismiss1 .and. dismiss2) then
      point_in_polygon = .false.
      return
   end if

   if (c > tol) then

      if (s1 == 1 .and. s2 == 1) then
         point_in_polygon = .false.
      else
         point_in_polygon = .true.
      end if

   else if (abs(c) <= tol) then

      if (s1 == 1) then
         point_in_polygon = .false.
      else
         point_in_polygon = .true.
      end if

   else

      if (s1 == -1 .and. s2 == -1) then
         point_in_polygon = .true.
      else
         point_in_polygon = .false.
      end if

   end if

end function point_in_polygon

subroutine closest_vertex(yy, zz, nv, yp, zp, ic)
   implicit none

   real(rkind), intent(in) :: yy(:), zz(:)
   integer, intent(in) :: nv
   real(rkind), intent(in) :: yp, zp
   integer, intent(out) :: ic

   integer :: p
   real(rkind) :: d2, d2min

   ic = 1
   d2min = (yp - yy(1))**2 + (zp - zz(1))**2

   do p = 2, nv
      d2 = (yp - yy(p))**2 + (zp - zz(p))**2

      if (d2 < d2min) then
         d2min = d2
         ic = p
      end if
   end do

end subroutine closest_vertex

subroutine edge_perspective_sign(ya, za, yb, zb, yp, zp, ccw, tol, s, L, dismiss)
   implicit none

   real(rkind), intent(in) :: ya, za
   real(rkind), intent(in) :: yb, zb
   real(rkind), intent(in) :: yp, zp
   logical, intent(in) :: ccw
   real(rkind), intent(in) :: tol

   integer, intent(out) :: s
   logical, intent(out) :: L
   logical, intent(out) :: dismiss

   real(rkind) :: ey, ez
   real(rkind) :: ny, nz
   real(rkind) :: my, mz
   real(rkind) :: ry, rz
   real(rkind) :: d

   ey = yb - ya
   ez = zb - za

   ! Outward normal.
   !
   ! For a counter-clockwise polygon, the interior is on the left of each
   ! directed edge, so the outward normal is the right normal.
   !
   ! Edge vector:       e = (ey, ez)
   ! Right normal:      n = ( ez, -ey)
   ! Left normal:       n = (-ez,  ey)
   if (ccw) then
      ny =  ez
      nz = -ey
   else
      ny = -ez
      nz =  ey
   end if

   my = half * (ya + yb)
   mz = half * (za + zb)

   ry = yp - my
   rz = zp - mz

   d = ny*ry + nz*rz

   L = point_on_segment(ya, za, yb, zb, yp, zp, tol)

   if (d > tol) then
      s = 1
      dismiss = .false.
   else if (d < -tol) then
      s = -1
      dismiss = .false.
   else
      s = 0

      if (L) then
         dismiss = .false.
      else
         ! Point is on the extension of the edge, but not on the edge itself.
         ! This parent edge perspective is dismissed.
         dismiss = .true.
      end if
   end if

end subroutine edge_perspective_sign

logical function point_on_segment(ya, za, yb, zb, yp, zp, tol)
   implicit none

   real(rkind), intent(in) :: ya, za
   real(rkind), intent(in) :: yb, zb
   real(rkind), intent(in) :: yp, zp
   real(rkind), intent(in) :: tol

   real(rkind) :: ey, ez
   real(rkind) :: py, pz
   real(rkind) :: crossp
   real(rkind) :: dotp
   real(rkind) :: len2, elen

   ey = yb - ya
   ez = zb - za

   py = yp - ya
   pz = zp - za

   crossp = ey*pz - ez*py
   len2 = ey*ey + ez*ez
   elen = sqrt(len2)

   if (len2 <= tol**2) then
      point_on_segment = (py*py + pz*pz <= tol**2)
      return
   end if

   if (abs(crossp) > tol * elen) then
      point_on_segment = .false.
      return
   end if

   dotp = py*ey + pz*ez

   point_on_segment = dotp >= -tol*elen .and. dotp <= len2 + tol*elen

end function point_on_segment

subroutine polygon_signed_area(yy, zz, nv, area_signed)
   implicit none

   real(rkind), intent(in) :: yy(:), zz(:)
   integer, intent(in) :: nv
   real(rkind), intent(out) :: area_signed

   integer :: p, pp
   real(rkind) :: accum
   real(rkind) :: yc, zc
   real(rkind) :: yp, zp, ypp, zpp

   if (nv < 3) then
      area_signed = zero
      return
   end if

   yc = sum(yy(1:nv)) / real(nv, rkind)
   zc = sum(zz(1:nv)) / real(nv, rkind)

   accum = zero

   do p = 1, nv
      pp = p + 1
      if (pp > nv) pp = 1

      yp  = yy(p)  - yc
      zp  = zz(p)  - zc
      ypp = yy(pp) - yc
      zpp = zz(pp) - zc

      accum = accum + yp*zpp - zp*ypp
   end do

   area_signed = half * accum

end subroutine polygon_signed_area

real(rkind) function point_polygon_tolerance(yy, zz, nv)
   implicit none

   real(rkind), intent(in) :: yy(:), zz(:)
   integer, intent(in) :: nv

   real(rkind) :: ymin, ymax
   real(rkind) :: zmin, zmax
   real(rkind) :: scale

   if (nv < 1) then
      point_polygon_tolerance = 100.0_rkind * epsilon(one)
      return
   end if

   ymin = minval(yy(1:nv))
   ymax = maxval(yy(1:nv))
   zmin = minval(zz(1:nv))
   zmax = maxval(zz(1:nv))

   scale = max(ymax - ymin, zmax - zmin)
   scale = max(scale, one)

   point_polygon_tolerance = 100.0_rkind * epsilon(one) * scale

end function point_polygon_tolerance

subroutine compact_valid_contour(yy, zz, valid, yv, zv, nv)
   implicit none

   real(rkind), intent(in) :: yy(:), zz(:)
   logical, intent(in) :: valid(:)
   real(rkind), intent(out) :: yv(:), zv(:)
   integer, intent(out) :: nv

   integer :: p, np

   np = size(yy)

   if (size(zz) /= np) call gracefulExit('compact_valid_contour: zz size mismatch.', 520)
   if (size(valid) /= np) call gracefulExit('compact_valid_contour: valid size mismatch.', 521)
   if (size(yv) < np) call gracefulExit('compact_valid_contour: yv too small.', 522)
   if (size(zv) < np) call gracefulExit('compact_valid_contour: zv too small.', 523)

   nv = 0

   do p = 1, np
      if (valid(p)) then
         nv = nv + 1
         yv(nv) = yy(p)
         zv(nv) = zz(p)
      end if
   end do

end subroutine compact_valid_contour

subroutine check_streamtube_valid_topology(valid)
   implicit none

   logical, intent(in) :: valid(:,:)

   integer :: n

   do n = 1, size(valid,1)
      if (valid_mask_has_gap(valid(n,:))) then
         call gracefulExit( &
            'Streamtube contour has multiple invalid gaps; polygon topology is ambiguous.', &
            560)
      end if
   end do

end subroutine check_streamtube_valid_topology

end module streamtube_mod


program constructStreamtube
   use streamtube_mod

   implicit none

   integer :: ioUnit, ierr
   character(len=clen) :: inputfile, ers, filename

   character(len=clen) :: outputdir = '.'
   character(len=clen) :: inputdir = '.'

   ! User-facing namelist controls.
   character(len=clen) :: initial_contour_file = ''
   character(len=clen) :: output_streamtube_file = 'streamtube.csv'
   character(len=clen) :: output_geometry_file = 'streamtube_geometry.csv'
   character(len=clen) :: output_stl_file = 'streamtube.stl'
   character(len=2) :: crid, cbrid
   character(len=6) :: ctid, ccount

   integer :: RunID, BaseRunID, TID, counter
   integer :: num_x_stations = 101
   integer :: num_contour_points = 0
   real(rkind) :: x_start = zero
   real(rkind) :: x_end = one
   logical :: write_geometry = .true., write_stl = .true., exists

   ! Internal names retained for consistency with the streamtube module.
   integer :: nxm, np

   real(rkind), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   real(rkind), allocatable :: x_march(:)
   real(rkind), allocatable :: y0(:), z0(:)
   real(rkind), allocatable :: tube_y(:,:), tube_z(:,:)
   real(rkind), allocatable :: area(:), width_y(:), height_z(:), flow_rate(:)
   logical, allocatable :: valid(:,:)

   namelist /INPUT/ nx, ny, nz, Lx, Ly, Lz, prow, pcol, NumericalSchemeVert, &
                    PeriodicInZ, botWall, topWall, &
                    RunID, BaseRunID, TID, counter,&
                    initial_contour_file, inputdir, outputdir, &
                    x_start, x_end, num_x_stations, num_contour_points, &
                    output_streamtube_file, write_geometry, output_geometry_file, &
                    output_stl_file, write_stl

   call MPI_Init(ierr)
   call get_command_argument(1, inputfile)

   if (len_trim(inputfile) == 0) then
      call gracefulExit('Usage: constructStreamtube input.dat', 100)
   end if

   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', status='old', action='read')
   read(unit=ioUnit, NML=INPUT, IOSTAT=ierr)
   if (ierr /= 0) then
      write(ers,'(I0)') ierr
      call gracefulExit('Reading failed for INPUT with error '//trim(ers), 101)
   end if
   close(ioUnit)

   if (len_trim(initial_contour_file) == 0) call gracefulExit('INPUT requires initial_contour_file.', 105)
   if (num_x_stations < 2) call gracefulExit('num_x_stations must be at least 2.', 106)

   nxm = num_x_stations

   ! Count the initial contour points unless the user supplies the expected count.
   if (num_contour_points > 0) then
      np = num_contour_points
   else
      call count_initial_contour_points(trim(initial_contour_file), np)
   end if

   if (np < 3) call gracefulExit('The initial contour must contain at least 3 valid points.', 107)

   call initialize_streamtube_module(nx, ny, nz, Lx, Ly, Lz, prow, pcol)

   allocate(u(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)), source=zero)
   allocate(v(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)), source=zero)
   allocate(w(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)), source=zero)
   allocate(rbuffxC(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)), source=zero)

   allocate(x_march(nxm))
   allocate(y0(np), z0(np))
   allocate(tube_y(nxm,np), tube_z(nxm,np))

   call build_x_march(x_start, x_end, x_march)
   call read_initial_contour_csv(trim(initial_contour_file), y0, z0)

   write(crid, '(I2.2)') RunID
   write(cbrid, '(I2.2)') BaseRunID
   write(ctid, '(I6.6)') TID
   write(ccount, '(I6.6)') counter

   ! Read the velocity fields
   
   ! Delta u
   filename = trim(inputdir)//'/Run'//trim(crid)//'_comp_deficit_budget0_term01_t'//trim(ctid)//'_n'//trim(ccount)//'.s3D'
   inquire(file=trim(filename), exist=exists)
   if (.not. exists) call gracefulExit('Missing: '//trim(filename), 707)
   call message(0, 'Reading '//trim(filename))
   call decomp_2d_read_one(1, rbuffxC, trim(filename), gpC)
   u = rbuffxC

   ! Base u
   filename = trim(inputdir)//'/Run'//trim(cbrid)//'_budget0_term01_t'//trim(ctid)//'_n'//trim(ccount)//'.s3D'
   inquire(file=trim(filename), exist=exists)
   if (.not. exists) call gracefulExit('Missing: '//trim(filename), 708)
   call message(0, 'Reading '//trim(filename))
   call decomp_2d_read_one(1, rbuffxC, trim(filename), gpC)
   u = u + rbuffxC

   ! Delta v
   filename = trim(inputdir)//'/Run'//trim(crid)//'_comp_deficit_budget0_term02_t'//trim(ctid)//'_n'//trim(ccount)//'.s3D'
   inquire(file=trim(filename), exist=exists)
   if (.not. exists) call gracefulExit('Missing: '//trim(filename), 709)
   call message(0, 'Reading '//trim(filename))
   call decomp_2d_read_one(1, rbuffxC, trim(filename), gpC)
   v = rbuffxC

   ! Base v
   filename = trim(inputdir)//'/Run'//trim(cbrid)//'_budget0_term02_t'//trim(ctid)//'_n'//trim(ccount)//'.s3D'
   inquire(file=trim(filename), exist=exists)
   if (.not. exists) call gracefulExit('Missing: '//trim(filename), 710)
   call message(0, 'Reading '//trim(filename))
   call decomp_2d_read_one(1, rbuffxC, trim(filename), gpC)
   v = v + rbuffxC

   ! Delta w
   filename = trim(inputdir)//'/Run'//trim(crid)//'_comp_deficit_budget0_term03_t'//trim(ctid)//'_n'//trim(ccount)//'.s3D'
   inquire(file=trim(filename), exist=exists)
   if (.not. exists) call gracefulExit('Missing: '//trim(filename), 711)
   call message(0, 'Reading '//trim(filename))
   call decomp_2d_read_one(1, rbuffxC, trim(filename), gpC)
   w = rbuffxC

   ! Base w
   filename = trim(inputdir)//'/Run'//trim(cbrid)//'_budget0_term03_t'//trim(ctid)//'_n'//trim(ccount)//'.s3D'
   inquire(file=trim(filename), exist=exists)
   if (.not. exists) call gracefulExit('Missing: '//trim(filename), 712)
   call message(0, 'Reading '//trim(filename))
   call decomp_2d_read_one(1, rbuffxC, trim(filename), gpC)
   w = w + rbuffxC
   
   call march_streamtube_xpencil(u, v, w, x_march, y0, z0, tube_y, tube_z)

   allocate(valid(nxm,np))
   call build_streamtube_valid_mask(tube_y, tube_z, valid)

   call check_streamtube_valid_topology(valid)

   if (write_geometry) then
      call message(0, "Post Processing the streamtube ...")
      allocate(area(nxm), width_y(nxm), height_z(nxm), flow_rate(nxm))

      call compute_contour_area(tube_y, tube_z, valid, area)
      call compute_contour_extents(tube_y, tube_z, valid, width_y, height_z)
      call compute_streamtube_flow_rate(u, x_march, tube_y, tube_z, valid, flow_rate)
   end if

   if(nrank == 0) then
      call export_streamtube_csv(trim(outputdir)//'/'//trim(output_streamtube_file), x_march, tube_y, tube_z)

      if (write_geometry) then
         call export_geometry_csv(trim(outputdir)//'/'//trim(output_geometry_file), &
                                  x_march, area, width_y, height_z, flow_rate)
      end if

      if(write_stl)then
         call export_streamtube_stl(trim(outputdir)//'/'//trim(output_stl_file), &
                                    x_march, tube_y, tube_z, valid)
      end if
   end if

   if (write_geometry) then
      deallocate(area, width_y, height_z, flow_rate)
   end if

   deallocate(tube_y, tube_z)
   deallocate(y0, z0)
   deallocate(x_march)
   deallocate(u, v, w)
   deallocate(valid)

   call release_streamtube_module()
   call MPI_Finalize(ierr)

end program constructStreamtube
