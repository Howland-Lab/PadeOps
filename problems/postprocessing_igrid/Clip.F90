module clip_mod
   use mpi
   use decomp_2d
   use decomp_2d_io
   use exits,          only: message, gracefulExit
   use constants,      only: zero, half
   use kind_parameters, only: rkind, clen

   implicit none

   private

   public :: initialize_clipper
   public :: clip_and_write_file
   public :: finalize_clipper

   ! Full-domain descriptors
   type(decomp_info), target :: gpC, gpE

   ! Clipped-domain descriptor
   type(decomp_info), target :: gpClip

   ! MPI
   integer :: myrank = 0, nprocs = 1

   ! Global full-domain sizes
   integer :: nx = 0, ny = 0, nz = 0
   integer :: prow = 0, pcol = 0

   ! Physical dimensions
   real(rkind) :: Lx = zero, Ly = zero, Lz = zero
   real(rkind) :: dx = zero, dy = zero, dz = zero

   ! User bounds
   real(rkind) :: x1 = zero, x2 = zero
   real(rkind) :: y1 = zero, y2 = zero
   real(rkind) :: z1 = zero, z2 = zero

   ! Clipped global old-index bounds
   integer :: ix1 = 0, ix2 = 0
   integer :: iy1 = 0, iy2 = 0
   integer :: iz1 = 0, iz2 = 0

   ! Clipped global sizes
   integer :: nx_clip = 0
   integer :: ny_clip = 0
   integer :: nz_clip = 0

   ! Field grid type:
   !   'C' -> cell-centered, old global size nx,ny,nz
   !   'E' -> vertically staggered, old global size nx,ny,nz+1
   character(len=1) :: field_grid = 'C'

   ! Owner lookup for gpClip
   integer, allocatable :: clip_xst_all(:,:), clip_xen_all(:,:)

   ! Cached ownership of the old descriptor (set once in initialize_clipper).
   integer, allocatable :: old_xst_all(:,:), old_xen_all(:,:)

contains

!=======================================================================
! Initialize full-domain descriptors and clipped descriptor.
!=======================================================================
subroutine initialize_clipper(nx_in, ny_in, nz_in, Lx_in, Ly_in, Lz_in, &
                              prow_in, pcol_in, periodicbcs,             &
                              x1_in, x2_in, y1_in, y2_in, z1_in, z2_in, &
                              field_grid_in)
   implicit none

   integer, intent(in) :: nx_in, ny_in, nz_in
   integer, intent(in) :: prow_in, pcol_in
   real(rkind), intent(in) :: Lx_in, Ly_in, Lz_in
   real(rkind), intent(in) :: x1_in, x2_in
   real(rkind), intent(in) :: y1_in, y2_in
   real(rkind), intent(in) :: z1_in, z2_in
   logical, intent(in) :: periodicbcs(3)
   character(len=*), intent(in) :: field_grid_in

   integer :: ierr

   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

   nx = nx_in
   ny = ny_in
   nz = nz_in

   Lx = Lx_in
   Ly = Ly_in
   Lz = Lz_in

   prow = prow_in
   pcol = pcol_in

   x1 = x1_in
   x2 = x2_in
   y1 = y1_in
   y2 = y2_in
   z1 = z1_in
   z2 = z2_in

   field_grid = adjustl(field_grid_in)
   
   dx = Lx / real(nx, rkind)
   dy = Ly / real(ny, rkind)
   dz = Lz / real(nz, rkind)

   if (field_grid /= 'C' .and. field_grid /= 'c' .and. &
       field_grid /= 'E' .and. field_grid /= 'e') then
      call gracefulExit("field_grid must be either 'C' or 'E'.", 501)
   end if

   ! Full-domain decomp.
   call decomp_2d_init(nx, ny, nz, prow, pcol, periodicbcs)
   call get_decomp_info(gpC)

   ! Vertically staggered / edge descriptor.
   call decomp_info_init(nx, ny, nz + 1, gpE)

   call compute_clip_indices()

   ! New clipped descriptor.
   call decomp_info_init(nx_clip, ny_clip, nz_clip, gpClip)
   
   call build_clip_owner_lookup()
   ! Cache old-descriptor ownership once.
   ! field_grid is already validated and fixed at this point.
   if (field_grid == 'C' .or. field_grid == 'c') then
      call gather_descriptor_ownership(gpC, old_xst_all, old_xen_all)
   else
      call gather_descriptor_ownership(gpE, old_xst_all, old_xen_all)
   end if

   if (myrank == 0) then
      call message(0, 'Clipper initialized.')
      call message(0, '  field_grid = '//field_grid)
      call message(0, '  clipped old-index bounds:')
      call message(0, '    ix = '//trim(to_string(ix1))//' : '//trim(to_string(ix2)))
      call message(0, '    iy = '//trim(to_string(iy1))//' : '//trim(to_string(iy2)))
      call message(0, '    iz = '//trim(to_string(iz1))//' : '//trim(to_string(iz2)))
      call message(0, '  clipped size:')
      call message(0, '    nx_clip = '//trim(to_string(nx_clip)))
      call message(0, '    ny_clip = '//trim(to_string(ny_clip)))
      call message(0, '    nz_clip = '//trim(to_string(nz_clip)))
   end if

end subroutine initialize_clipper

!=======================================================================
! Compute index bounds for the selected grid.
!
! C grid:
!   x_i = (i - 1) dx
!   y_j = (j - 1) dy
!   z_k = (k - 1/2) dz
!
! E grid:
!   x_i = (i - 1) dx
!   y_j = (j - 1) dy
!   z_k = (k - 1) dz, k = 1,...,nz+1
!=======================================================================
subroutine compute_clip_indices()
   implicit none

   real(rkind) :: xmin, xmax
   real(rkind) :: ymin, ymax
   real(rkind) :: zmin, zmax
   real(rkind), parameter :: clip_eps = 100.0_rkind * epsilon(1.0_rkind)

   xmin = min(x1, x2)
   xmax = max(x1, x2)

   ymin = min(y1, y2)
   ymax = max(y1, y2)

   zmin = min(z1, z2)
   zmax = max(z1, z2)

   ix1 = max(1,  int(ceiling(xmin/dx + 1.0_rkind - clip_eps * (abs(xmin/dx) + 1.0_rkind))))
   ix2 = min(nx, int(floor  (xmax/dx + 1.0_rkind + clip_eps * (abs(xmax/dx) + 1.0_rkind))))
   iy1 = max(1,  int(ceiling(ymin/dy + 1.0_rkind - clip_eps * (abs(ymin/dy) + 1.0_rkind))))
   iy2 = min(ny, int(floor  (ymax/dy + 1.0_rkind + clip_eps * (abs(ymax/dy) + 1.0_rkind))))

   select case(field_grid)
   case('C','c')
      ! z_k = (k - 1/2) dz
      iz1 = max(1,  int(ceiling(zmin/dz + half - clip_eps* (abs(zmin/dz) + 1.0_rkind))))
      iz2 = min(nz, int(floor  (zmax/dz + half + clip_eps* (abs(zmax/dz) + 1.0_rkind))))

   case('E','e')
      ! z_k = (k - 1) dz, k = 1,...,nz+1
      iz1 = max(1,      int(ceiling(zmin/dz + 1.0_rkind - clip_eps* (abs(zmin/dz) + 1.0_rkind))))
      iz2 = min(nz + 1, int(floor  (zmax/dz + 1.0_rkind + clip_eps* (abs(zmax/dz) + 1.0_rkind))))
   end select

   if (ix2 < ix1 .or. iy2 < iy1 .or. iz2 < iz1) then
      call gracefulExit('Clip bounds do not intersect the selected grid.', 502)
   end if

   nx_clip = ix2 - ix1 + 1
   ny_clip = iy2 - iy1 + 1
   nz_clip = iz2 - iz1 + 1

end subroutine compute_clip_indices

subroutine build_clip_owner_lookup()
   implicit none

   integer :: ierr
   integer :: local_bounds(6)
   integer, allocatable :: all_bounds(:)

   integer :: r
   integer :: p
   
   if (allocated(clip_xst_all)) deallocate(clip_xst_all)
   if (allocated(clip_xen_all)) deallocate(clip_xen_all)
   
   allocate(clip_xst_all(3,nprocs))
   allocate(clip_xen_all(3,nprocs))
   allocate(all_bounds(6*nprocs))

   local_bounds(1:3) = gpClip%xst(1:3)
   local_bounds(4:6) = gpClip%xen(1:3)

   call MPI_Allgather(local_bounds, 6, MPI_INTEGER, &
                      all_bounds,    6, MPI_INTEGER, &
                      MPI_COMM_WORLD, ierr)

   if (ierr /= MPI_SUCCESS) then
      call gracefulExit('build_clip_owner_lookup: MPI_Allgather failed.', 602)
   end if

   do r = 1, nprocs
      p = 6*(r - 1)

      clip_xst_all(1,r) = all_bounds(p + 1)
      clip_xst_all(2,r) = all_bounds(p + 2)
      clip_xst_all(3,r) = all_bounds(p + 3)

      clip_xen_all(1,r) = all_bounds(p + 4)
      clip_xen_all(2,r) = all_bounds(p + 5)
      clip_xen_all(3,r) = all_bounds(p + 6)
   end do

   deallocate(all_bounds)

end subroutine build_clip_owner_lookup

subroutine gather_descriptor_ownership(gp, xst_all, xen_all)
   implicit none

   type(decomp_info), intent(in) :: gp

   integer, allocatable, intent(out) :: xst_all(:,:)
   integer, allocatable, intent(out) :: xen_all(:,:)

   integer :: ierr
   integer :: local_bounds(6)
   integer, allocatable :: all_bounds(:)

   integer :: r
   integer :: p

   if (allocated(xst_all)) deallocate(xst_all)
   if (allocated(xen_all)) deallocate(xen_all)

   allocate(xst_all(3,nprocs))
   allocate(xen_all(3,nprocs))
   allocate(all_bounds(6*nprocs))

   local_bounds(1:3) = gp%xst(1:3)
   local_bounds(4:6) = gp%xen(1:3)

   call MPI_Allgather(local_bounds, 6, MPI_INTEGER, &
                      all_bounds,    6, MPI_INTEGER, &
                      MPI_COMM_WORLD, ierr)

   if (ierr /= MPI_SUCCESS) then
      call gracefulExit('gather_descriptor_ownership: MPI_Allgather failed.', 605)
   end if

   do r = 1, nprocs
      p = 6*(r - 1)

      xst_all(1,r) = all_bounds(p + 1)
      xst_all(2,r) = all_bounds(p + 2)
      xst_all(3,r) = all_bounds(p + 3)

      xen_all(1,r) = all_bounds(p + 4)
      xen_all(2,r) = all_bounds(p + 5)
      xen_all(3,r) = all_bounds(p + 6)
   end do

   deallocate(all_bounds)

end subroutine gather_descriptor_ownership

!=======================================================================
! Clip one scalar 3D file.
!=======================================================================
subroutine clip_and_write_file(inputdir, infile, basefile, outputdir, outfile)
    implicit none

    character(len=*), intent(in) :: inputdir
    character(len=*), intent(in) :: infile
    character(len=*), intent(in) :: basefile
    character(len=*), intent(in) :: outputdir
    character(len=*), intent(in) :: outfile

    real(rkind), allocatable :: field_old(:,:,:), buffer(:,:,:)
    real(rkind), allocatable :: field_clip(:,:,:)
    character(len=clen) :: filename, outfile_local

    select case(field_grid)
    case('C','c')
        allocate(field_old(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
        allocate(buffer(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))

    case('E','e')
        allocate(field_old(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
        allocate(buffer(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
    end select

    allocate(field_clip(gpClip%xsz(1), gpClip%xsz(2), gpClip%xsz(3)))

    filename = trim(inputdir)//'/'//trim(infile)
    if (myrank == 0) then
        call message(0, 'Reading '//trim(filename))
    end if
    select case(field_grid)
    case('C','c')
        call decomp_2d_read_one(1, field_old, trim(filename), gpC)

    case('E','e')
        call decomp_2d_read_one(1, field_old, trim(filename), gpE)
    end select

    if(trim(basefile) /= '') then
        filename = trim(inputdir)//'/'//trim(basefile)
        if (myrank == 0) then
            call message(0, 'Reading '//trim(filename))
        end if
        select case(field_grid)
        case('C','c')
            call decomp_2d_read_one(1, buffer, trim(filename), gpC)
        case('E','e')
            call decomp_2d_read_one(1, buffer, trim(filename), gpE)
        end select

        field_old = field_old + buffer
    end if

    field_clip = zero
    call redistribute_clip_xpencil(field_old, field_clip)

    if(trim(outfile) == '') then
        outfile_local = 'clipped_'//trim(infile)
    else
        outfile_local = trim(outfile)
    end if
    filename = trim(outputdir)//'/'//trim(outfile_local)
    if (myrank == 0) then
        call message(0, 'Writing clipped field to '//trim(filename))  ! outfile may be ''
    end if
    call decomp_2d_write_one(1, field_clip, trim(filename), gpClip)

    deallocate(field_old)
    deallocate(field_clip)
    deallocate(buffer)

end subroutine clip_and_write_file

subroutine redistribute_clip_xpencil(field_old, field_clip)
   implicit none

   real(rkind), intent(in)  :: field_old(:,:,:)
   real(rkind), intent(out) :: field_clip(:,:,:)

   type(decomp_info), pointer :: gpOld

   integer :: ierr

   ! MPI source/destination ranks.
   integer :: src_rank, dst_rank

   ! Loop indices in clipped global coordinates.
   integer :: ig_clip, jg_clip, kg_clip

   ! Local indices in old and clipped arrays.
   integer :: i_old, j_old, k_old
   integer :: i_clip, j_clip, k_clip

   ! Intersection box in clipped global coordinates.
   integer :: ixlo, ixhi
   integer :: iylo, iyhi
   integer :: izlo, izhi
   logical :: has_box

   ! Communication metadata.
   integer, allocatable :: send_counts(:), recv_counts(:)
   integer, allocatable :: send_displs(:), recv_displs(:)

   ! Send/receive buffers.
   real(rkind), allocatable :: send_val(:), recv_val(:)

   ! Buffer positions.
   integer :: pos
   integer :: total_send, total_recv

   ! Optional consistency counters.
   integer :: expected_count
   integer :: packed_count
   integer :: unpacked_count

   !------------------------------------------------------------
   ! Select old descriptor.
   !
   ! C: cell-centered descriptor gpC, global size nx,ny,nz.
   ! E: vertically staggered descriptor gpE, global size nx,ny,nz+1.
   !------------------------------------------------------------
   if (field_grid == 'C' .or. field_grid == 'c') then
      gpOld => gpC
   else
      gpOld => gpE
   end if

   !------------------------------------------------------------
   ! Allocate communication arrays.
   !------------------------------------------------------------
   allocate(send_counts(nprocs), recv_counts(nprocs))
   allocate(send_displs(nprocs), recv_displs(nprocs))

   send_counts = 0
   recv_counts = 0
   send_displs = 0
   recv_displs = 0

   !------------------------------------------------------------
   ! Count how many values this rank sends to each destination.
   !
   ! For destination rank dst_rank, the send block is:
   !
   !   old ownership of myrank, mapped to clipped coordinates
   !             intersected with
   !   gpClip ownership of dst_rank
   !
   ! All bounds returned by get_source_dest_intersection are in
   ! clipped global coordinates.
   !------------------------------------------------------------
   do dst_rank = 0, nprocs - 1

      call get_source_dest_intersection(myrank, dst_rank, old_xst_all, old_xen_all, &
                                        ixlo, ixhi, iylo, iyhi, izlo, izhi, has_box)

      if (has_box) then
         send_counts(dst_rank + 1) = get_intersection_count(ixlo, ixhi, &
                                                            iylo, iyhi, &
                                                            izlo, izhi)
      else
         send_counts(dst_rank + 1) = 0
      end if

   end do

   !------------------------------------------------------------
   ! Exchange counts.
   !------------------------------------------------------------
   call MPI_Alltoall(send_counts, 1, MPI_INTEGER, &
                     recv_counts, 1, MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

   if (ierr /= MPI_SUCCESS) then
      call gracefulExit('redistribute_clip_xpencil: MPI_Alltoall counts failed.', 701)
   end if

   call prefix_displs(send_counts, send_displs)
   call prefix_displs(recv_counts, recv_displs)

   total_send = sum(send_counts)
   total_recv = sum(recv_counts)

   allocate(send_val(max(total_send,1)))
   allocate(recv_val(max(total_recv,1)))

   send_val = zero
   recv_val = zero

   !------------------------------------------------------------
   ! Pack values only.
   !
   ! Critical ordering:
   !
   !   do kg_clip = izlo, izhi
   !      do jg_clip = iylo, iyhi
   !         do ig_clip = ixlo, ixhi
   !
   ! The receiver reconstructs the same intersection box and uses
   ! the same loop order during unpacking.
   !------------------------------------------------------------
   do dst_rank = 0, nprocs - 1

      call get_source_dest_intersection(myrank, dst_rank, old_xst_all, old_xen_all, &
                                        ixlo, ixhi, iylo, iyhi, izlo, izhi, has_box)

      if (.not. has_box) cycle

      expected_count = get_intersection_count(ixlo, ixhi, iylo, iyhi, izlo, izhi)
      packed_count   = 0

      pos = send_displs(dst_rank + 1)

      do kg_clip = izlo, izhi
         k_old  = kg_clip + iz1 - gpOld%xst(3)

         do jg_clip = iylo, iyhi
            j_old  = jg_clip + iy1 - gpOld%xst(2)

            do ig_clip = ixlo, ixhi
               i_old  = ig_clip + ix1 - gpOld%xst(1)

               if (i_old < 1 .or. i_old > gpOld%xsz(1) .or. &
                   j_old < 1 .or. j_old > gpOld%xsz(2) .or. &
                   k_old < 1 .or. k_old > gpOld%xsz(3)) then
                  call gracefulExit('redistribute_clip_xpencil: old local index out of range during packing.', 702)
               end if

               pos = pos + 1
               packed_count = packed_count + 1

               send_val(pos) = field_old(i_old,j_old,k_old)

            end do
         end do
      end do

      if (packed_count /= expected_count) then
         call gracefulExit('redistribute_clip_xpencil: packed count mismatch.', 703)
      end if

      if (pos /= send_displs(dst_rank + 1) + send_counts(dst_rank + 1)) then
         call gracefulExit('redistribute_clip_xpencil: send displacement mismatch after packing.', 704)
      end if

   end do

   !------------------------------------------------------------
   ! Exchange values only.
   !------------------------------------------------------------
   call MPI_Alltoallv(send_val, send_counts, send_displs, MPI_DOUBLE_PRECISION, &
                      recv_val, recv_counts, recv_displs, MPI_DOUBLE_PRECISION, &
                      MPI_COMM_WORLD, ierr)

   if (ierr /= MPI_SUCCESS) then
      call gracefulExit('redistribute_clip_xpencil: MPI_Alltoallv values failed.', 705)
   end if

   !------------------------------------------------------------
   ! Unpack into local gpClip array.
   !
   ! For each source rank src_rank, reconstruct:
   !
   !   old ownership of src_rank, mapped to clipped coordinates
   !             intersected with
   !   gpClip ownership of myrank
   !
   ! Then unpack using the same kg-jg-ig ordering used by sender.
   !------------------------------------------------------------
   do src_rank = 0, nprocs - 1

      call get_source_dest_intersection(src_rank, myrank, old_xst_all, old_xen_all, &
                                        ixlo, ixhi, iylo, iyhi, izlo, izhi, has_box)

      if (.not. has_box) cycle

      expected_count = get_intersection_count(ixlo, ixhi, iylo, iyhi, izlo, izhi)
      unpacked_count = 0

      pos = recv_displs(src_rank + 1)

      do kg_clip = izlo, izhi

         k_clip = kg_clip - gpClip%xst(3) + 1

         do jg_clip = iylo, iyhi

            j_clip = jg_clip - gpClip%xst(2) + 1

            do ig_clip = ixlo, ixhi

               i_clip = ig_clip - gpClip%xst(1) + 1

               if (i_clip < 1 .or. i_clip > gpClip%xsz(1) .or. &
                   j_clip < 1 .or. j_clip > gpClip%xsz(2) .or. &
                   k_clip < 1 .or. k_clip > gpClip%xsz(3)) then
                  call gracefulExit('redistribute_clip_xpencil: clipped local index out of range during unpacking.', 706)
               end if

               pos = pos + 1
               unpacked_count = unpacked_count + 1

               field_clip(i_clip,j_clip,k_clip) = recv_val(pos)

            end do
         end do
      end do

      if (unpacked_count /= expected_count) then
         call gracefulExit('redistribute_clip_xpencil: unpacked count mismatch.', 707)
      end if

      if (pos /= recv_displs(src_rank + 1) + recv_counts(src_rank + 1)) then
         call gracefulExit('redistribute_clip_xpencil: receive displacement mismatch after unpacking.', 708)
      end if

   end do

   !------------------------------------------------------------
   ! Cleanup.
   !------------------------------------------------------------
   deallocate(send_counts, recv_counts)
   deallocate(send_displs, recv_displs)
   deallocate(send_val, recv_val)
   nullify(gpOld)

end subroutine redistribute_clip_xpencil

pure integer function get_intersection_count(ixlo, ixhi, iylo, iyhi, izlo, izhi) result(npts)
   implicit none

   integer, intent(in) :: ixlo, ixhi
   integer, intent(in) :: iylo, iyhi
   integer, intent(in) :: izlo, izhi

   if (ixhi < ixlo .or. iyhi < iylo .or. izhi < izlo) then
      npts = 0
   else
      npts = (ixhi - ixlo + 1) * &
             (iyhi - iylo + 1) * &
             (izhi - izlo + 1)
   end if

end function get_intersection_count

subroutine get_source_dest_intersection(src_rank, dst_rank, old_xst_all, old_xen_all, &
                                        ixlo, ixhi, iylo, iyhi, izlo, izhi, has_box)
   implicit none

   integer, intent(in) :: src_rank
   integer, intent(in) :: dst_rank

   integer, intent(in) :: old_xst_all(:,:)
   integer, intent(in) :: old_xen_all(:,:)

   integer, intent(out) :: ixlo, ixhi
   integer, intent(out) :: iylo, iyhi
   integer, intent(out) :: izlo, izhi

   logical, intent(out) :: has_box

   integer :: src_i1, src_i2
   integer :: src_j1, src_j2
   integer :: src_k1, src_k2

   integer :: dst_i1, dst_i2
   integer :: dst_j1, dst_j2
   integer :: dst_k1, dst_k2

   integer :: src
   integer :: dst

   ! src_rank and dst_rank are MPI ranks: 0,...,nprocs-1.
   ! Array storage uses Fortran indexing: 1,...,nprocs.
   src = src_rank + 1
   dst = dst_rank + 1

   if (src < 1 .or. src > nprocs) then
      has_box = .false.
      ixlo = 1; ixhi = 0
      iylo = 1; iyhi = 0
      izlo = 1; izhi = 0
      return
   end if

   if (dst < 1 .or. dst > nprocs) then
      has_box = .false.
      ixlo = 1; ixhi = 0
      iylo = 1; iyhi = 0
      izlo = 1; izhi = 0
      return
   end if

   !------------------------------------------------------------
   ! Source old ownership mapped into clipped global coordinates.
   !
   ! Old global index:
   !   ig_old
   !
   ! Clipped global index:
   !   ig_clip = ig_old - ix1 + 1
   !
   ! Therefore:
   !   ig_old = ix1 maps to ig_clip = 1
   !------------------------------------------------------------
   src_i1 = max(1,       old_xst_all(1,src) - ix1 + 1)
   src_i2 = min(nx_clip, old_xen_all(1,src) - ix1 + 1)

   src_j1 = max(1,       old_xst_all(2,src) - iy1 + 1)
   src_j2 = min(ny_clip, old_xen_all(2,src) - iy1 + 1)

   src_k1 = max(1,       old_xst_all(3,src) - iz1 + 1)
   src_k2 = min(nz_clip, old_xen_all(3,src) - iz1 + 1)

   !------------------------------------------------------------
   ! Destination gpClip ownership is already in clipped global
   ! coordinates.
   !------------------------------------------------------------
   dst_i1 = clip_xst_all(1,dst)
   dst_i2 = clip_xen_all(1,dst)

   dst_j1 = clip_xst_all(2,dst)
   dst_j2 = clip_xen_all(2,dst)

   dst_k1 = clip_xst_all(3,dst)
   dst_k2 = clip_xen_all(3,dst)

   !------------------------------------------------------------
   ! Source-destination intersection in clipped coordinates.
   !------------------------------------------------------------
   ixlo = max(src_i1, dst_i1)
   ixhi = min(src_i2, dst_i2)

   iylo = max(src_j1, dst_j1)
   iyhi = min(src_j2, dst_j2)

   izlo = max(src_k1, dst_k1)
   izhi = min(src_k2, dst_k2)

   has_box = (ixhi >= ixlo .and. iyhi >= iylo .and. izhi >= izlo)

end subroutine get_source_dest_intersection

subroutine prefix_displs(counts, displs)
   implicit none

   integer, intent(in)  :: counts(:)
   integer, intent(out) :: displs(:)

   integer :: r

   if (size(displs) /= size(counts)) then
      call gracefulExit('prefix_displs: counts and displs have inconsistent sizes.', 601)
   end if

   if (size(counts) < 1) return

   displs(1) = 0

   do r = 2, size(counts)
      displs(r) = displs(r-1) + counts(r-1)
   end do

end subroutine prefix_displs

!=======================================================================
! Cleanup.
!=======================================================================
subroutine finalize_clipper()
    implicit none

    if (allocated(clip_xst_all)) deallocate(clip_xst_all)
    if (allocated(clip_xen_all)) deallocate(clip_xen_all)
    if (allocated(old_xst_all))  deallocate(old_xst_all)
    if (allocated(old_xen_all))  deallocate(old_xen_all)

    call decomp_info_finalize(gpClip)
    call decomp_info_finalize(gpE)
    call decomp_info_finalize(gpC)
    call decomp_2d_finalize()

end subroutine finalize_clipper

!=======================================================================
! Small integer-to-string helper.
!=======================================================================
pure function to_string(i) result(str)
   implicit none

   integer, intent(in) :: i
   character(len=32) :: str

   write(str, '(I0)') i

end function to_string

end module clip_mod

program clip
   use kind_parameters, only: rkind, clen
   use exits,           only: gracefulExit
   use mpi
   use clip_mod
   
   implicit none

   integer :: ierr
   integer :: ioUnit
   character(len=clen) :: inputfile=''
   character(len=clen) :: infile='', outfile='', basefile=''
   character(len=clen) :: ers
   character(len=clen) :: inputdir='.', outputdir='.'

   integer :: nx=0, ny=0, nz=0
   integer :: prow=0, pcol=0
   real(rkind) :: Lx=0.0_rkind, Ly=0.0_rkind, Lz=0.0_rkind
   real(rkind) :: x1=0.0_rkind, x2=0.0_rkind, y1=0.0_rkind, y2=0.0_rkind, z1=0.0_rkind, z2=0.0_rkind
   logical :: periodicbcs(3)
   logical :: periodic_x = .true., periodic_y=.true., periodic_z = .false.
   character(len=1) :: field_grid='C'

   namelist /INPUT/ infile, outfile, basefile, nx, ny, nz, Lx, Ly, Lz,&
                    prow, pcol, field_grid, &
                    x1, x2, y1, y2, z1, z2, periodic_x, periodic_y, periodic_z,&
                    inputdir, outputdir

   call MPI_Init(ierr)

   call GET_COMMAND_ARGUMENT(1, inputfile)

   if (len_trim(inputfile) == 0) then
      call gracefulExit('clip3d: cannot read input file', 100)
   end if

   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='formatted', status='old', action='read')

   read(unit=ioUnit, nml=INPUT, iostat=ierr)
   if (ierr /= 0) then
      write(ers,'(I0)') ierr
      call gracefulExit('Reading failed for INPUT with error '//trim(ers), 101)
   end if

   close(ioUnit)

   periodicbcs(1) = periodic_x
   periodicbcs(2) = periodic_y
   periodicbcs(3) = periodic_z

   call initialize_clipper(nx, ny, nz, Lx, Ly, Lz, &
                           prow, pcol, periodicbcs, &
                           x1, x2, y1, y2, z1, z2, &
                           field_grid)

   call clip_and_write_file(trim(inputdir), trim(infile), trim(basefile), trim(outputdir), trim(outfile))

   call finalize_clipper()

   call MPI_Finalize(ierr)

end program clip