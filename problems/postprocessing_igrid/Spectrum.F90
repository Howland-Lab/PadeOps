module spectrum_mod
   use mpi
   use decomp_2d
   use decomp_2d_io
   use kind_parameters, only: rkind, clen, mpirkind
   use constants, only: zero, one, two, half, pi
   use exits, only: message, gracefulExit
   use reductions, only: p_sum
   use spectralMod, only: spectral

   implicit none

   character(len=clen) :: inputdir, outputdir, fields
   integer :: nx, ny, nz, prow=0, pcol=0
   real(rkind) :: Lx=one, Ly=one, Lz=one, dx=one, dy=one, dz=one
   real(rkind) :: npts=one, nhorz=one
   logical :: remove_spatial_mean=.false.
   logical :: remove_horizontal_mean=.false.
   logical :: include_one_half=.false.
   logical :: write_density=.false.
   logical :: write_x_spectrum=.true.
   logical :: write_y_spectrum=.true.
   logical :: write_vertical_summary=.true.
   logical :: write_height_spectra=.false.
   logical :: write_yz_plane_spectra=.false.

   type(decomp_info) :: gpC
   type(decomp_info) :: gpYSpec
   type(spectral) :: spectC
   real(rkind), dimension(:,:,:), allocatable :: field, rbuffxC
   real(rkind), dimension(:,:,:), allocatable :: field_y
   complex(rkind), dimension(:,:,:), allocatable :: fhat
   complex(rkind), dimension(:,:,:), allocatable :: field_yhat
   real(rkind), dimension(:), allocatable :: kbin, spectrum_local, spectrum_global
   integer, dimension(:), allocatable :: counts_local, counts_global
   real(rkind), dimension(:), allocatable :: kxbin, kybin, spectrum_x_local, spectrum_x_global
   real(rkind), dimension(:), allocatable :: spectrum_y_local, spectrum_y_global
   integer, dimension(:), allocatable :: counts_x_local, counts_x_global
   integer, dimension(:), allocatable :: counts_y_local, counts_y_global
   real(rkind), dimension(:,:), allocatable :: spectrum_height_local, spectrum_height_global
   integer, dimension(:,:), allocatable :: counts_height_local, counts_height_global
   real(rkind), dimension(:,:), allocatable :: spectrum_yzplane_local, spectrum_yzplane_global
   integer :: nbins, ierr
   integer :: nbins_x, nbins_y
   real(rkind) :: dk, dkx, dky, normfact, normfact_yzplane
   integer(kind=8) :: plan_r2c_y_plane = 0

   include "fftw3.f"

   type field_component
      character(len=clen) :: filename = ''
      real(rkind) :: scale = one
   end type field_component

   type spectrum_field
      character(len=clen) :: name = ''
      integer :: ncomponents = 0
      type(field_component), allocatable :: components(:)
   end type spectrum_field

contains

   subroutine read_field_specs(filename, specs)
      implicit none

      character(len=*), intent(in) :: filename
      type(spectrum_field), allocatable, intent(out) :: specs(:)
      integer :: unit, ios, nfields, ifield, icomp
      character(len=2048) :: line, clean_line, field_name
      logical :: in_block, waiting_for_open
      integer, allocatable :: component_counts(:)

      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) call gracefulExit("Could not open field-spec file.", 112)

      ! First pass: count fields. Braces may be on the field-name line or
      ! on the following line.
      nfields = 0
      in_block = .false.
      waiting_for_open = .false.
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit

         call strip_comments(line, clean_line)
         if (len_trim(clean_line) == 0) cycle

         if (waiting_for_open) then
            if (trim(clean_line) /= '{') call gracefulExit("Field-spec parser expected '{'.", 113)
            waiting_for_open = .false.
            in_block = .true.
            cycle
         end if

         if (in_block) then
            if (trim(clean_line) == '}') then
               in_block = .false.
            else
               component_counts(nfields) = component_counts(nfields) + 1
            end if
            cycle
         end if

         nfields = nfields + 1
         if (.not. allocated(component_counts)) allocate(component_counts(1024))
         if (nfields > size(component_counts)) call grow_integer_array(component_counts)
         component_counts(nfields) = 0
         if (index(clean_line, '{') > 0) then
            in_block = .true.
         else
            waiting_for_open = .true.
         end if
      end do

      close(unit)

      if (ios > 0) call gracefulExit("Error while reading field-spec file.", 556)
      if (nfields == 0) call gracefulExit("Field-spec file contains no fields.", 114)
      if (in_block .or. waiting_for_open) call gracefulExit("Field-spec file ended before a field block was closed.", 115)

      allocate(specs(nfields))
      do ifield = 1, nfields
         if (component_counts(ifield) == 0) call gracefulExit("Field block contains no components.", 119)
         allocate(specs(ifield)%components(component_counts(ifield)))
      end do

      ! Second pass: read field names and component lines.
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) call gracefulExit("Could not reopen field-spec file.", 116)

      ifield = 0
      icomp = 0
      in_block = .false.
      waiting_for_open = .false.
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit

         call strip_comments(line, clean_line)
         if (len_trim(clean_line) == 0) cycle

         if (waiting_for_open) then
            if (trim(clean_line) /= '{') call gracefulExit("Field-spec parser expected '{'.", 117)
            waiting_for_open = .false.
            in_block = .true.
            icomp = 0
            cycle
         end if

         if (.not. in_block) then
            ifield = ifield + 1
            call parse_field_name(clean_line, field_name, in_block)
            specs(ifield)%name = trim(field_name)
            specs(ifield)%ncomponents = component_counts(ifield)
            icomp = 0
            if (.not. in_block) waiting_for_open = .true.
            cycle
         end if

         if (trim(clean_line) == '}') then
            in_block = .false.
         else
            icomp = icomp + 1
            call parse_component(clean_line, &
               specs(ifield)%components(icomp)%filename, &
               specs(ifield)%components(icomp)%scale)
         end if
      end do
      close(unit)

      if (ios > 0) call gracefulExit("Error while reading field-spec file.", 557)
   end subroutine read_field_specs

   subroutine grow_integer_array(values)
      implicit none
      integer, allocatable, intent(inout) :: values(:)
      integer, allocatable :: old_values(:)
      integer :: old_size

      old_size = size(values)
      allocate(old_values(old_size))
      old_values = values
      deallocate(values)
      allocate(values(2*old_size))
      values(1:old_size) = old_values
      values(old_size+1:) = 0
      deallocate(old_values)
   end subroutine grow_integer_array

   subroutine parse_field_name(line, field_name, block_is_open)
      implicit none
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: field_name
      logical, intent(out) :: block_is_open
      integer :: brace_pos

      brace_pos = index(line, '{')
      block_is_open = brace_pos > 0
      if (block_is_open) then
         field_name = adjustl(line(:brace_pos-1))
      else
         field_name = adjustl(line)
      end if
      if (len_trim(field_name) == 0) call gracefulExit("Field-spec parser found an empty field name.", 120)
   end subroutine parse_field_name

   subroutine parse_component(line, filename, scale)
      implicit none
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: filename
      real(rkind), intent(out) :: scale
      character(len=2048) :: normalized
      character(len=clen) :: token1, token2, token3, token4
      integer :: ios

      call normalize_component_line(line, normalized)
      token1 = ''
      token2 = ''
      token3 = ''
      token4 = ''
      read(normalized, *, iostat=ios) token1, token2
      if (ios /= 0) call gracefulExit("Malformed field component line: "//trim(line), 121)

      if (trim(token1) == 'file') then
         read(normalized, *, iostat=ios) token1, token2, token3, token4
         if (ios /= 0) call gracefulExit("Malformed keyed field component line: "//trim(line), 121)
         filename = trim(token2)
         if (trim(token3) /= 'scale') call gracefulExit("Component line missing scale key: "//trim(line), 122)
         read(token4, *, iostat=ios) scale
      else
         filename = trim(token1)
         read(token2, *, iostat=ios) scale
      end if

      if (ios /= 0) call gracefulExit("Could not parse component scale: "//trim(line), 123)
      if (len_trim(filename) == 0) call gracefulExit("Component line has an empty filename.", 124)
   end subroutine parse_component

   subroutine normalize_component_line(line, normalized)
      implicit none
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: normalized
      integer :: i

      normalized = line
      do i = 1, len_trim(normalized)
         select case (normalized(i:i))
         case ('{', '}', ',', ':')
            normalized(i:i) = ' '
         end select
      end do
      normalized = adjustl(normalized)
   end subroutine normalize_component_line

   subroutine strip_comments(line, clean_line)
      implicit none
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: clean_line
      integer :: bang_pos, hash_pos, comment_pos

      clean_line = line
      bang_pos = index(clean_line, '!')
      hash_pos = index(clean_line, '#')
      comment_pos = 0
      if (bang_pos > 0) comment_pos = bang_pos
      if ((hash_pos > 0) .and. ((comment_pos == 0) .or. (hash_pos < comment_pos))) comment_pos = hash_pos
      if (comment_pos > 0) clean_line(comment_pos:) = ' '
      clean_line = adjustl(clean_line)
   end subroutine strip_comments

   subroutine read_field(spec)
      implicit none
      type(spectrum_field), intent(in) :: spec
      character(len=clen) :: infile
      integer :: i
      
      field = zero
      rbuffxC = zero
      do i = 1, spec%ncomponents
         infile = resolve_input_path(spec%components(i)%filename)
         call message(1, 'Reading '//trim(infile))
         call decomp_2d_read_one(1, rbuffxC, trim(infile), gpC)
         field = field + spec%components(i)%scale*rbuffxC
      end do
      
   end subroutine read_field

   function resolve_input_path(filename) result(path)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=clen) :: path

      if (filename(1:1) == '/') then
         path = trim(filename)
      else
         path = trim(inputdir)//'/'//trim(filename)
      end if
   end function resolve_input_path

   subroutine remove_requested_means()
      real(rkind) :: mean_value

      if (remove_spatial_mean) then
         mean_value = p_sum(sum(field))/npts
         field = field - mean_value
         call message(0, 'Removed volume mean:', mean_value)
      end if

      if (remove_horizontal_mean) then
         call remove_horizontal_profile()
      end if
   end subroutine remove_requested_means

   subroutine remove_horizontal_profile()
      real(rkind), dimension(:), allocatable :: profile_local, profile_global
      integer :: i, j, k, kg

      allocate(profile_local(nz), profile_global(nz))
      profile_local = zero
      profile_global = zero

      do k = 1, gpC%xsz(3)
         kg = gpC%xst(3) + k - 1
         do j = 1, gpC%xsz(2)
            do i = 1, gpC%xsz(1)
               profile_local(kg) = profile_local(kg) + field(i,j,k)
            end do
         end do
      end do

      call MPI_Allreduce(profile_local, profile_global, nz, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
      profile_global = profile_global/nhorz

      do k = 1, gpC%xsz(3)
         kg = gpC%xst(3) + k - 1
         field(:,:,k) = field(:,:,k) - profile_global(kg)
      end do

      call message(0, 'Removed horizontally averaged mean profile')
      deallocate(profile_local, profile_global)
   end subroutine remove_horizontal_profile

   subroutine init_bins()
      integer :: b
      real(rkind) :: kmax

      dk = min(two*pi/Lx, two*pi/Ly)
      dkx = two*pi/Lx
      dky = two*pi/Ly
      kmax = sqrt((pi/dx)**2 + (pi/dy)**2)
      nbins = int(kmax/dk) + 2
      nbins_x = nx/2 + 1
      nbins_y = ny/2 + 1
      normfact = one/(nhorz*nhorz*real(nz,rkind))
      normfact_yzplane = one/(real(ny,rkind)*real(ny,rkind)*real(nz,rkind))

      allocate(kbin(nbins), spectrum_local(nbins), spectrum_global(nbins))
      allocate(counts_local(nbins), counts_global(nbins))
      allocate(kxbin(nbins_x), kybin(nbins_y))
      allocate(spectrum_x_local(nbins_x), spectrum_x_global(nbins_x))
      allocate(spectrum_y_local(nbins_y), spectrum_y_global(nbins_y))
      allocate(counts_x_local(nbins_x), counts_x_global(nbins_x))
      allocate(counts_y_local(nbins_y), counts_y_global(nbins_y))
      allocate(spectrum_height_local(nbins,nz), spectrum_height_global(nbins,nz))
      allocate(counts_height_local(nbins,nz), counts_height_global(nbins,nz))
      allocate(spectrum_yzplane_local(nbins_y,nx), spectrum_yzplane_global(nbins_y,nx))

      do b = 1, nbins
         kbin(b) = (real(b,rkind) - half)*dk
      end do

      do b = 1, nbins_x
         kxbin(b) = real(b-1,rkind)*dkx
      end do

      do b = 1, nbins_y
         kybin(b) = real(b-1,rkind)*dky
      end do

      spectrum_local = zero
      spectrum_global = zero
      counts_local = 0
      counts_global = 0
      spectrum_x_local = zero
      spectrum_x_global = zero
      spectrum_y_local = zero
      spectrum_y_global = zero
      counts_x_local = 0
      counts_x_global = 0
      counts_y_local = 0
      counts_y_global = 0
      spectrum_height_local = zero
      spectrum_height_global = zero
      counts_height_local = 0
      counts_height_global = 0
      spectrum_yzplane_local = zero
      spectrum_yzplane_global = zero
   end subroutine init_bins

   subroutine init_yz_plane_spectra()
      integer :: ierr_local
      integer :: n_sizeact, n_sizeinput, n_sizeoutput, n_howmany, n_jump, n_chunk
      real(rkind), dimension(:,:), allocatable :: real_arr_2d
      complex(rkind), dimension(:,:), allocatable :: cmplx_arr_2d

      if (.not. write_yz_plane_spectra) return

      call decomp_info_init(nx, ny/2 + 1, nz, gpYSpec)
      allocate(field_y(gpC%ysz(1), gpC%ysz(2), gpC%ysz(3)), stat=ierr_local)
      if (ierr_local /= 0) call gracefulExit("Could not allocate y-pencil field buffer.", 130)
      allocate(field_yhat(gpYSpec%ysz(1), gpYSpec%ysz(2), gpYSpec%ysz(3)), stat=ierr_local)
      if (ierr_local /= 0) call gracefulExit("Could not allocate y-plane spectrum buffer.", 131)

      allocate(real_arr_2d(gpC%ysz(1), gpC%ysz(2)), stat=ierr_local)
      if (ierr_local /= 0) call gracefulExit("Could not allocate y FFT planning real buffer.", 132)
      allocate(cmplx_arr_2d(gpYSpec%ysz(1), gpYSpec%ysz(2)), stat=ierr_local)
      if (ierr_local /= 0) call gracefulExit("Could not allocate y FFT planning complex buffer.", 133)

      n_sizeact = gpC%ysz(2)
      n_sizeinput = gpC%ysz(2)
      n_sizeoutput = gpYSpec%ysz(2)
      n_howmany = gpC%ysz(1)
      n_jump = gpC%ysz(1)
      n_chunk = 1
      call dfftw_plan_many_dft_r2c(plan_r2c_y_plane, 1, n_sizeact, &
         n_howmany, real_arr_2d, n_sizeinput, n_jump, n_chunk, &
         cmplx_arr_2d, n_sizeoutput, n_jump, n_chunk, FFTW_MEASURE)

      deallocate(real_arr_2d, cmplx_arr_2d)
      field_y = zero
      field_yhat = cmplx(zero, zero, kind=rkind)
   end subroutine init_yz_plane_spectra

   subroutine compute_spectrum()
      integer :: i, j, k, ig, jg, kg, ibin, ixbin, iybin, multiplicity
      real(rkind) :: kmag, amp2, factor

      factor = one
      if (include_one_half) factor = half

      spectrum_local = zero
      counts_local = 0
      spectrum_x_local = zero
      spectrum_y_local = zero
      counts_x_local = 0
      counts_y_local = 0
      spectrum_height_local = zero
      counts_height_local = 0

      do k = 1, size(fhat,3)
         kg = spectC%spectdecomp%yst(3) + k - 1
         do j = 1, size(fhat,2)
            jg = spectC%spectdecomp%yst(2) + j - 1
            iybin = y_abs_bin(jg)
            do i = 1, size(fhat,1)
               kmag = sqrt(spectC%kabs_sq(i,j,k))
               ibin = int(kmag/dk) + 1

               ig = spectC%spectdecomp%yst(1) + i - 1
               ixbin = ig
               multiplicity = hermitian_multiplicity(ig)
               amp2 = factor*real(multiplicity,rkind)* &
                  real(fhat(i,j,k)*conjg(fhat(i,j,k)), rkind)*normfact

               if ((ibin >= 1) .and. (ibin <= nbins)) then
                  spectrum_local(ibin) = spectrum_local(ibin) + amp2
                  counts_local(ibin) = counts_local(ibin) + multiplicity
                  spectrum_height_local(ibin,kg) = spectrum_height_local(ibin,kg) + amp2*real(nz,rkind)
                  counts_height_local(ibin,kg) = counts_height_local(ibin,kg) + multiplicity
               end if

               if ((ixbin >= 1) .and. (ixbin <= nbins_x)) then
                  spectrum_x_local(ixbin) = spectrum_x_local(ixbin) + amp2
                  counts_x_local(ixbin) = counts_x_local(ixbin) + multiplicity
               end if

               if ((iybin >= 1) .and. (iybin <= nbins_y)) then
                  spectrum_y_local(iybin) = spectrum_y_local(iybin) + amp2
                  counts_y_local(iybin) = counts_y_local(iybin) + multiplicity
               end if
            end do
         end do
      end do

      call MPI_Reduce(spectrum_local, spectrum_global, nbins, mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(counts_local, counts_global, nbins, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(spectrum_x_local, spectrum_x_global, nbins_x, mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(spectrum_y_local, spectrum_y_global, nbins_y, mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(counts_x_local, counts_x_global, nbins_x, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(counts_y_local, counts_y_global, nbins_y, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(spectrum_height_local, spectrum_height_global, nbins*nz, &
         mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(counts_height_local, counts_height_global, nbins*nz, &
         MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if (write_density .and. nrank == 0) then
         spectrum_global = spectrum_global/dk
         spectrum_x_global = spectrum_x_global/dkx
         spectrum_y_global = spectrum_y_global/dky
         spectrum_height_global = spectrum_height_global/dk
      end if
   end subroutine compute_spectrum

   integer function y_abs_bin(jg)
      integer, intent(in) :: jg

      if (jg <= ny/2 + 1) then
         y_abs_bin = jg
      else
         y_abs_bin = ny - jg + 2
      end if
   end function y_abs_bin

   integer function hermitian_multiplicity(ig)
      integer, intent(in) :: ig

      if ((ig == 1) .or. (ig == nx/2 + 1)) then
         hermitian_multiplicity = 1
      else
         hermitian_multiplicity = 2
      end if
   end function hermitian_multiplicity

   integer function y_hermitian_multiplicity(iy)
      integer, intent(in) :: iy

      if ((iy == 1) .or. (iy == ny/2 + 1)) then
         y_hermitian_multiplicity = 1
      else
         y_hermitian_multiplicity = 2
      end if
   end function y_hermitian_multiplicity

   subroutine compute_yz_plane_spectrum()
      integer :: i, j, k, ig, multiplicity
      real(rkind) :: amp2, factor

      if (.not. write_yz_plane_spectra) return

      factor = one
      if (include_one_half) factor = half

      spectrum_yzplane_local = zero
      field_y = zero
      field_yhat = cmplx(zero, zero, kind=rkind)

      call transpose_x_to_y(field, field_y, gpC)

      do k = 1, gpC%ysz(3)
         call dfftw_execute_dft_r2c(plan_r2c_y_plane, field_y(:,:,k), field_yhat(:,:,k))
      end do

      do k = 1, size(field_yhat,3)
         do j = 1, size(field_yhat,2)
            multiplicity = y_hermitian_multiplicity(j)
            do i = 1, size(field_yhat,1)
               ig = gpYSpec%yst(1) + i - 1
               amp2 = factor*real(multiplicity,rkind)* &
                  real(field_yhat(i,j,k)*conjg(field_yhat(i,j,k)), rkind)*normfact_yzplane
               if ((ig >= 1) .and. (ig <= nx)) then
                  spectrum_yzplane_local(j,ig) = spectrum_yzplane_local(j,ig) + amp2
               end if
            end do
         end do
      end do

      call MPI_Reduce(spectrum_yzplane_local, spectrum_yzplane_global, nbins_y*nx, &
         mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if (write_density .and. nrank == 0) spectrum_yzplane_global = spectrum_yzplane_global/dky
   end subroutine compute_yz_plane_spectrum

   subroutine parseval_check()
      real(rkind) :: physical_energy, spectral_energy, spectral_x_energy, spectral_y_energy, factor

      factor = one
      if (include_one_half) factor = half

      physical_energy = factor*p_sum(sum(field*field))/npts
      if (nrank == 0) then
         if (write_density) then
            spectral_energy = sum(spectrum_global)*dk
         else
            spectral_energy = sum(spectrum_global)
         end if
         call message(0, 'Physical-space variance/energy:', physical_energy)
         call message(0, 'Horizontal spectrum-integrated energy:', spectral_energy)
         call message(0, 'Horizontal spectrum Parseval absolute error:', abs(spectral_energy - physical_energy))
         if (write_density) then
            spectral_x_energy = sum(spectrum_x_global)*dkx
            spectral_y_energy = sum(spectrum_y_global)*dky
         else
            spectral_x_energy = sum(spectrum_x_global)
            spectral_y_energy = sum(spectrum_y_global)
         end if
         call message(0, 'Streamwise spectrum-integrated energy:', spectral_x_energy)
         call message(0, 'Spanwise spectrum-integrated energy:', spectral_y_energy)
      end if
   end subroutine parseval_check

   subroutine export_csv(field_name)
      character(len=*), intent(in) :: field_name
      character(len=clen) :: outfile
      integer :: unit, b

      if (nrank /= 0) return

      outfile = trim(outputdir)//'/spectrum_'//trim(sanitize_field_name(field_name))//'.csv'
      call message(0, 'Writing vertically averaged horizontal spectrum to '//trim(outfile))

      open(newunit=unit, file=trim(outfile), status='replace', action='write', form='formatted')
      write(unit, '(A)') 'k,E'
      do b = 1, nbins
         if (counts_global(b) > 0) write(unit, '(ES24.16,",",ES24.16)') kbin(b), spectrum_global(b)
      end do
      close(unit)
   end subroutine export_csv

   subroutine export_directional_csv(field_name)
      character(len=*), intent(in) :: field_name
      character(len=clen) :: outfile, clean_name
      integer :: unit, b

      if (nrank /= 0) return

      clean_name = sanitize_field_name(field_name)

      if (write_x_spectrum) then
         outfile = trim(outputdir)//'/spectrum_x_'//trim(clean_name)//'.csv'
         call message(0, 'Writing streamwise spectrum to '//trim(outfile))
         open(newunit=unit, file=trim(outfile), status='replace', action='write', form='formatted')
         write(unit, '(A)') 'kx,E'
         do b = 1, nbins_x
            if (counts_x_global(b) > 0) write(unit, '(ES24.16,",",ES24.16)') kxbin(b), spectrum_x_global(b)
         end do
         close(unit)
      end if

      if (write_y_spectrum) then
         outfile = trim(outputdir)//'/spectrum_y_'//trim(clean_name)//'.csv'
         call message(0, 'Writing spanwise spectrum to '//trim(outfile))
         open(newunit=unit, file=trim(outfile), status='replace', action='write', form='formatted')
         write(unit, '(A)') 'ky,E'
         do b = 1, nbins_y
            if (counts_y_global(b) > 0) write(unit, '(ES24.16,",",ES24.16)') kybin(b), spectrum_y_global(b)
         end do
         close(unit)
      end if
   end subroutine export_directional_csv

   subroutine export_vertical_summary_csv(field_name)
      character(len=*), intent(in) :: field_name
      character(len=clen) :: outfile, clean_name
      integer :: unit, b, kg
      real(rkind) :: energy_sum, zc, zspread, zg

      if (nrank /= 0) return
      if (.not. write_vertical_summary) return

      clean_name = sanitize_field_name(field_name)
      outfile = trim(outputdir)//'/spectrum_zsummary_'//trim(clean_name)//'.csv'
      call message(0, 'Writing vertical spectrum summary to '//trim(outfile))

      open(newunit=unit, file=trim(outfile), status='replace', action='write', form='formatted')
      write(unit, '(A)') 'k,E,z_centroid,z_spread'
      do b = 1, nbins
         if (counts_global(b) <= 0) cycle
         energy_sum = sum(spectrum_height_global(b,:))
         if (energy_sum <= zero) cycle

         zc = zero
         do kg = 1, nz
            zg = (real(kg,rkind) - half)*dz
            zc = zc + zg*spectrum_height_global(b,kg)
         end do
         zc = zc/energy_sum

         zspread = zero
         do kg = 1, nz
            zg = (real(kg,rkind) - half)*dz
            zspread = zspread + (zg - zc)**2*spectrum_height_global(b,kg)
         end do
         zspread = sqrt(zspread/energy_sum)

         write(unit, '(ES24.16,",",ES24.16,",",ES24.16,",",ES24.16)') &
            kbin(b), spectrum_global(b), zc, zspread
      end do
      close(unit)
   end subroutine export_vertical_summary_csv

   subroutine export_height_spectra_csv(field_name)
      character(len=*), intent(in) :: field_name
      character(len=clen) :: outfile, clean_name
      integer :: unit, b, kg
      real(rkind) :: zg

      if (nrank /= 0) return
      if (.not. write_height_spectra) return

      clean_name = sanitize_field_name(field_name)
      outfile = trim(outputdir)//'/spectrum_height_'//trim(clean_name)//'.csv'
      call message(0, 'Writing height-resolved horizontal spectra to '//trim(outfile))

      open(newunit=unit, file=trim(outfile), status='replace', action='write', form='formatted')
      write(unit, '(A)') 'z,k,E'
      do kg = 1, nz
         zg = (real(kg,rkind) - half)*dz
         do b = 1, nbins
            if (counts_height_global(b,kg) > 0) write(unit, '(ES24.16,",",ES24.16,",",ES24.16)') &
               zg, kbin(b), spectrum_height_global(b,kg)
         end do
      end do
      close(unit)
   end subroutine export_height_spectra_csv

   subroutine export_yz_plane_spectra_csv(field_name)
      character(len=*), intent(in) :: field_name
      character(len=clen) :: outfile, clean_name
      integer :: unit, b, ig
      real(rkind) :: xg

      if (nrank /= 0) return
      if (.not. write_yz_plane_spectra) return

      clean_name = sanitize_field_name(field_name)
      outfile = trim(outputdir)//'/spectrum_yzplane_'//trim(clean_name)//'.csv'
      call message(0, 'Writing y spectra for each y-z plane to '//trim(outfile))

      open(newunit=unit, file=trim(outfile), status='replace', action='write', form='formatted')
      write(unit, '(A)') 'x,ky,E'
      do ig = 1, nx
         xg = (real(ig,rkind) - half)*dx
         do b = 1, nbins_y
            write(unit, '(ES24.16,",",ES24.16,",",ES24.16)') xg, kybin(b), spectrum_yzplane_global(b,ig)
         end do
      end do
      close(unit)
   end subroutine export_yz_plane_spectra_csv

   function sanitize_field_name(field_name) result(clean_name)
      implicit none
      character(len=*), intent(in) :: field_name
      character(len=clen) :: clean_name
      integer :: i

      clean_name = adjustl(field_name)
      do i = 1, len_trim(clean_name)
         select case (clean_name(i:i))
         case (' ', '/', '\', ':', ',', ';', '{', '}', '(', ')', '[', ']')
            clean_name(i:i) = '_'
         end select
      end do
   end function sanitize_field_name

end module spectrum_mod

program spectrum
   use mpi
   use decomp_2d
   use decomp_2d_io
   use spectrum_mod
   use exits, only: message, gracefulExit

   implicit none

   character(len=clen) :: inputfile, infile
   integer :: ioUnit, ierr_local, ispec
   logical :: exists
   type(spectrum_field), allocatable :: specs(:)
      
   namelist /INPUT/ inputdir, outputdir, nx, ny, nz, Lx, Ly, Lz, prow, pcol, fields, &
                    remove_spatial_mean, remove_horizontal_mean, include_one_half, write_density, &
                    write_x_spectrum, write_y_spectrum, write_vertical_summary, write_height_spectra, &
                    write_yz_plane_spectra

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   call GET_COMMAND_ARGUMENT(1, inputfile)
   if (len_trim(inputfile) == 0) call gracefulExit('Usage: Spectrum.x input.dat', 100)

   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', status='old', iostat=ierr_local)
   if (ierr_local /= 0) call gracefulExit('Could not open input namelist file.', 101)
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)

   if (mod(nx,2) /= 0) call gracefulExit('nx must be even for PadeOps real-to-complex FFT storage.', 102)
   if (mod(ny,2) /= 0) call gracefulExit('ny must be even.', 103)
   if (mod(nz,2) /= 0) call gracefulExit('nz must be even.', 104)

   dx = Lx/real(nx,rkind)
   dy = Ly/real(ny,rkind)
   dz = Lz/real(nz,rkind)
   npts = real(nx,rkind)*real(ny,rkind)*real(nz,rkind)
   nhorz = real(nx,rkind)*real(ny,rkind)

   call decomp_2d_init(nx, ny, nz, prow, pcol)
   call decomp_info_init(nx, ny, nz, gpC)

   call spectC%init('x', nx, ny, nz, dx, dy, dz, 'FOUR', '2/3rd', &
      dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
   
   allocate(field(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   allocate(rbuffxC(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   call spectC%alloc_r2c_out(fhat)

   if (len_trim(fields) == 0) call gracefulExit('No field-spec file was provided in fields.', 105)
   if (fields(1:1) == '/') then
      infile = trim(fields)
   else
      infile = trim(inputdir)//'/'//trim(fields)
   end if
   inquire(file=trim(infile), exist=exists)
   if (.not. exists) call gracefulExit('Input field-spec file not found: '//trim(infile), 106)
   call init_bins()
   call init_yz_plane_spectra()
   call read_field_specs(trim(infile), specs)

   do ispec = 1, size(specs)
      call message(0, 'Computing spectrum for '//trim(specs(ispec)%name))
      call read_field(specs(ispec))
      call remove_requested_means()
      call compute_yz_plane_spectrum()
      call spectC%fft(field, fhat)
      call compute_spectrum()
      call parseval_check()
      call export_csv(specs(ispec)%name)
      call export_directional_csv(specs(ispec)%name)
      call export_vertical_summary_csv(specs(ispec)%name)
      call export_height_spectra_csv(specs(ispec)%name)
      call export_yz_plane_spectra_csv(specs(ispec)%name)
   end do

   deallocate(field, rbuffxC, fhat, kbin, spectrum_local, spectrum_global, counts_local, counts_global)
   deallocate(kxbin, kybin, spectrum_x_local, spectrum_x_global, spectrum_y_local, spectrum_y_global)
   deallocate(counts_x_local, counts_x_global, counts_y_local, counts_y_global)
   deallocate(spectrum_height_local, spectrum_height_global, counts_height_local, counts_height_global)
   deallocate(spectrum_yzplane_local, spectrum_yzplane_global)
   if (write_yz_plane_spectra) then
      call dfftw_destroy_plan(plan_r2c_y_plane)
      deallocate(field_y, field_yhat)
      call decomp_info_finalize(gpYSpec)
   end if
   call spectC%destroy()
   call decomp_info_finalize(gpC)
   call decomp_2d_finalize()
   call MPI_Finalize(ierr)

end program spectrum
