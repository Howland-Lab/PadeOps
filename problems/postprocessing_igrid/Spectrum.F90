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

   type(decomp_info) :: gpC
   type(spectral) :: spectC
   real(rkind), dimension(:,:,:), allocatable :: field, rbuffxC
   complex(rkind), dimension(:,:,:), allocatable :: fhat
   real(rkind), dimension(:), allocatable :: kbin, spectrum_local, spectrum_global
   integer, dimension(:), allocatable :: counts_local, counts_global
   integer :: nbins, ierr
   real(rkind) :: dk, normfact

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
      kmax = sqrt((pi/dx)**2 + (pi/dy)**2)
      nbins = int(kmax/dk) + 2
      normfact = one/(nhorz*nhorz*real(nz,rkind))

      allocate(kbin(nbins), spectrum_local(nbins), spectrum_global(nbins))
      allocate(counts_local(nbins), counts_global(nbins))

      do b = 1, nbins
         kbin(b) = (real(b,rkind) - half)*dk
      end do

      spectrum_local = zero
      spectrum_global = zero
      counts_local = 0
      counts_global = 0
   end subroutine init_bins

   subroutine compute_spectrum()
      integer :: i, j, k, ig, ibin, multiplicity
      real(rkind) :: kmag, amp2, factor

      factor = one
      if (include_one_half) factor = half

      spectrum_local = zero
      counts_local = 0

      do k = 1, size(fhat,3)
         do j = 1, size(fhat,2)
            do i = 1, size(fhat,1)
               kmag = sqrt(spectC%kabs_sq(i,j,k))
               ibin = int(kmag/dk) + 1

               if ((ibin >= 1) .and. (ibin <= nbins)) then
                  ig = spectC%spectdecomp%yst(1) + i - 1
                  multiplicity = hermitian_multiplicity(ig)
                  amp2 = real(fhat(i,j,k)*conjg(fhat(i,j,k)), rkind)
                  spectrum_local(ibin) = spectrum_local(ibin) + &
                     factor*real(multiplicity,rkind)*amp2*normfact
                  counts_local(ibin) = counts_local(ibin) + multiplicity
               end if
            end do
         end do
      end do

      call MPI_Reduce(spectrum_local, spectrum_global, nbins, mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(counts_local, counts_global, nbins, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if (write_density .and. nrank == 0) spectrum_global = spectrum_global/dk
   end subroutine compute_spectrum

   integer function hermitian_multiplicity(ig)
      integer, intent(in) :: ig

      if ((ig == 1) .or. (ig == nx/2 + 1)) then
         hermitian_multiplicity = 1
      else
         hermitian_multiplicity = 2
      end if
   end function hermitian_multiplicity

   subroutine parseval_check()
      real(rkind) :: physical_energy, spectral_energy, factor

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
         call message(0, 'Spectrum-integrated energy:', spectral_energy)
         call message(0, 'Parseval absolute error:', abs(spectral_energy - physical_energy))
      end if
   end subroutine parseval_check

   subroutine export_csv(field_name)
      character(len=*), intent(in) :: field_name
      character(len=clen) :: outfile
      integer :: unit, b

      if (nrank /= 0) return

      outfile = trim(outputdir)//'/spectrum_'//trim(sanitize_field_name(field_name))//'.csv'
      call message(0, 'Writing spectrum to '//trim(outfile))

      open(newunit=unit, file=trim(outfile), status='replace', action='write', form='formatted')
      write(unit, '(A)') 'k,E'
      do b = 1, nbins
         if (counts_global(b) > 0) write(unit, '(ES24.16,",",ES24.16)') kbin(b), spectrum_global(b)
      end do
      close(unit)
   end subroutine export_csv

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
                    remove_spatial_mean, remove_horizontal_mean, include_one_half, write_density

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

   call spectC%init('x', nx, ny, nz, dx, dy, dz, 'FOUR', '2/3rd',dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
   
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
   call read_field_specs(trim(infile), specs)

   do ispec = 1, size(specs)
      call message(0, 'Computing spectrum for '//trim(specs(ispec)%name))
      call read_field(specs(ispec))
      call remove_requested_means()
      call spectC%fft(field, fhat)
      call compute_spectrum()
      call parseval_check()
      call export_csv(specs(ispec)%name)
   end do

   deallocate(field, rbuffxC, fhat, kbin, spectrum_local, spectrum_global, counts_local, counts_global)
   call spectC%destroy()
   call decomp_info_finalize(gpC)
   call decomp_2d_finalize()
   call MPI_Finalize(ierr)

end program spectrum
