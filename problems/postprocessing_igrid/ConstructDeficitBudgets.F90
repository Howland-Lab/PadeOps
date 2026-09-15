module constructDeficitBudgets_mod
   use mpi
   use exits, only: message, gracefulExit
   use constants, only: one, two, zero, half
   use kind_parameters,  only: rkind, clen
   use timer, only: tic, toc
   use PadeDerOps, only: Pade6stagg
   use spectralMod, only: spectral
   use decomp_2d 
   use decomp_2d_io

   implicit none

   character(len=clen) :: inputdir, outputdir
   real(rkind) :: Lx = one, Ly = one, Lz = one
   integer :: botWall=3, topWall=2, botBC_temp=0
   logical :: PeriodicInZ=.false.
   integer :: prow=0, pcol=0, nx, ny, nz, RID, BRID, NumericalSchemeVert=1
   integer :: startIDX=-1, endIDX=999999
   logical :: writeDependentVariables = .false.
   logical :: do_box_averaging=.true.
   logical :: do_x_budget=.true.
   logical :: do_y_budget=.false.
   logical :: do_z_budget=.false.
   logical :: do_TKE_budget=.false.
   logical :: do_MKE_budget=.false.
   logical :: do_TMP_budget=.false.
   real(rkind) :: x1=zero, x2=zero, y1=zero, y2=zero, z1=zero, z2=zero
   ! ------------------------------------------------------------------ !
   
   type(spectral), target  :: spectE, spectC
   type(decomp_info) :: gpC, gpE
   type(decomp_info), pointer :: sp_gpC, sp_gpE
   type(Pade6stagg) :: Pade6opZ
   real(rkind) :: dx, dy, dz
   real(rkind), dimension(:,:,:,:), allocatable, target :: mesh
   integer :: uBC_bottom, uBC_top, vBC_bottom, vBC_top, wBC_bottom, wBC_top
   integer :: nx_box, ix1g, ix2g
   real(rkind), dimension(:), allocatable :: xstations
   character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)
   integer :: num_profiles
   real(rkind), dimension(:,:), allocatable :: xprofiles
   real(rkind), dimension(:,:), allocatable :: yprofiles
   real(rkind), dimension(:,:), allocatable :: zprofiles
   real(rkind), dimension(:,:), allocatable :: TKEprofiles
   real(rkind), dimension(:,:), allocatable :: MKEprofiles
   real(rkind), dimension(:,:), allocatable :: TMPprofiles
   real(rkind), dimension(:,:,:,:), allocatable, target :: rbuffxC
   real(rkind), dimension(:,:,:), pointer :: du, dv, dw
   real(rkind), dimension(:,:,:), pointer :: ubase, vbase, wbase
   real(rkind), dimension(:,:,:), pointer :: buffer, bf
   complex(rkind), dimension(:,:,:), allocatable :: cbuffyC
   complex(rkind), dimension(:,:,:,:), allocatable :: cbuffzC
   
   contains

   function csv_file_name(key, stamp, name) result(filename)
      implicit none
      character(len=*), intent(in) :: key, stamp, name
      character(len=clen) :: filename
      character(len=3) :: crid

      write(crid, '(I2.2)') RID
      filename = trim(outputdir)//'/Run'//trim(crid)//'_t'//trim(key)//'_n'//trim(stamp)//'_'//trim(name)//'_Budgets_XProfile_'//'.csv'
   end function

   subroutine export_csv(filename, profiles)
      implicit none
      character(len=*), intent(in) :: filename
      real(rkind), dimension(:,:), intent(in) :: profiles
      character(len=3) :: tid
      integer :: i, j
      integer :: nx, ny
      integer :: unit

      nx = size(profiles, 1)
      ny = size(profiles, 2)
      unit =1045

      call message(1, 'Exporting profiles to '//trim(filename))

      ! Open file
      open(newunit=unit, file=filename, status='replace', action='write', form='formatted')

      write(unit, '(A1,",")', advance='no') 'x'
      do j = 1, ny
         write(tid, '(I3.3)') j
         if (j < ny) then
            write(unit, '(A,",")', advance='no') 'T'//trim(tid)
         else
            write(unit, '(A)') 'T'//trim(tid)
         end if
      end do

      ! Write data row by row
      do i = 1, nx
         write(unit, '(ES16.8,",")', advance='no') xstations(i)
         do j = 1, ny
               if (j < ny) then
                  write(unit, '(ES16.8,",")', advance='no') profiles(i,j)
               else
                  write(unit, '(ES16.8)') profiles(i,j)
               end if
         end do
      end do

      close(unit)
   end subroutine export_csv

   subroutine dump_budget_field(field, fieldID, BudgetID, key, stamp)
        real(rkind), dimension(:,:,:), intent(in) :: field
        character(len=*), intent(in) :: key, stamp, fieldID, BudgetID
        character(len=clen) :: fname, tempname 
        character(len=2) :: crid

        write(crid, '(I2.2)') RID
        write(tempname,"(A)") "Run"//crid//"_comp_deficit_budget"//BudgetID//"_term"//fieldID//"_t"//trim(key)//"_n"//trim(stamp)//".s3D"
        fname = trim(outputdir)//"/"//trim(tempname)

        call message(2, 'Writing a budget field to '//trim(fname))
        call decomp_2d_write_one(1,field, trim(fname), gpC)
   end subroutine

   function depedent_variable(budgetType, idx)
      implicit none
      integer, intent(in) :: idx, budgetType
      logical :: depedent_variable

      depedent_variable = .false.
      if((budgetType == 1) .or. (budgetType == 2) .or. (budgetType == 3))then
         ! X, Y, or Z momentum equation
         if((idx < 10) .or. (idx > 15)) depedent_variable = .true.
      elseif(budgetType == 4)then
         ! TKE equation
         if(idx <= 12) depedent_variable = .true.
      elseif(budgetType == 5)then
         ! MKE equation
         depedent_variable = .true.
      elseif(budgetType == 6)then
         ! TMP
         depedent_variable = .true.
      end if
   end function depedent_variable

   subroutine intersectBoxAndMesh()
      implicit none

      integer :: iL
      integer :: ig
      real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax, xplane
      character(len=4) :: ix1gc, ix2gc 
      integer, parameter :: HUGE_I = huge(1)
      integer :: ibox

      ! We use x-pencils. All ranks see the whole x range. All calculations here are local.

      !----------------------------
      ! Bounds (make robust to x1>x2 etc.)
      !----------------------------
      xmin = min(x1, x2);  xmax = max(x1, x2)
      ymin = min(y1, y2);  ymax = max(y1, y2)
      zmin = min(z1, z2);  zmax = max(z1, z2)

      ix1g =  HUGE_I
      ix2g = -HUGE_I

      do iL = 1, size(mesh,1)
         ig = gpC%xst(1) + (iL - 1)  ! local-to-global x index

         ! x is constant on an x-plane for structured meshes; sample one point on that plane
         xplane = mesh(iL, 1, 1, 1)

         if (xplane >= xmin .and. xplane <= xmax) then
            ix1g = min(ix1g, ig)
            ix2g = max(ix2g, ig)
         end if
      end do

      ! Handle: box does not intersect any x-plane anywhere
      if (ix2g < ix1g .or. ix1g == HUGE_I .or. ix2g == -HUGE_I) then
         call gracefulExit('Invalid box bounds.', 124)
      end if

      nx_box = ix2g - ix1g + 1
      write(ix1gc, '(I4.4)')ix1g
      write(ix2gc, '(I4.4)')ix2g
      call message(0,'Box intersects X dimension between indices '//trim(ix1gc)//' and '//trim(ix2gc))

      allocate(xstations(nx_box))
      do ibox = 1, nx_box
         iL = ix1g + ibox - 1
         xstations(ibox) = mesh(iL,1,1,1)
      end do
   end subroutine

   subroutine integrate_box_yz(f, prof)
      implicit none
      real(rkind), intent(in)        :: f(:,:,:)          ! local field: (xsz1,xsz2,xsz3)
      real(rkind), dimension(:), intent(out) :: prof

      ! Locals
      integer :: ierr
      integer :: iL
      integer :: ig
      real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax, xplane
      real(rkind), allocatable :: prof_local(:)
      logical, allocatable :: mask_yz(:,:)
      
      prof = zero
      allocate(prof_local(nx_box))
      prof_local = zero

      ! mask over local y-z plane
      allocate(mask_yz(size(f,2), size(f,3)))

      !----------------------------
      ! Bounds (make robust to x1>x2 etc.)
      !----------------------------
      xmin = min(x1, x2);  xmax = max(x1, x2)
      ymin = min(y1, y2);  ymax = max(y1, y2)
      zmin = min(z1, z2);  zmax = max(z1, z2)

      !----------------------------
      ! Local contribution: for each local x-plane that lies in [xmin,xmax],
      ! sum f over (y,z) points whose (y,z) are within box bounds.
      ! Accumulate into prof_local at the position corresponding to global x-index.
      !----------------------------
      do iL = 1, size(f,1)
         ig = gpC%xst(1) + (iL - 1)
         xplane = mesh(iL, 1, 1, 1)

         if (xplane >= xmin .and. xplane <= xmax) then
            ! mask for this x-plane in y-z
            mask_yz = (mesh(iL, :, :, 2) >= ymin .and. mesh(iL, :, :, 2) <= ymax) .and. &
                        (mesh(iL, :, :, 3) >= zmin .and. mesh(iL, :, :, 3) <= zmax)

            prof_local(ig - ix1g + 1) = prof_local(ig - ix1g + 1) + sum(f(iL, :, :), mask=mask_yz)
         end if
      end do
      prof_local = prof_local * dy*dz ! Area element

      !----------------------------
      ! Global reduction: sum contributions from all ranks
      !----------------------------
      call mpi_allreduce(prof_local, prof, nx_box, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      deallocate(mask_yz)
      deallocate(prof_local)
   end subroutine integrate_box_yz

   logical function TimeWithinRange(tidx, istart, iend)
      implicit none
      character(*), intent(in) :: tidx
      integer, intent(in) :: istart, iend
      integer :: itime
      integer :: ios

      read(tidx, '(I6)', iostat=ios) itime
      if (ios /= 0) then
         TimeWithinRange = .false.
         return
      end if
      TimeWithinRange = (itime >= istart .and. itime <= iend)
   end function TimeWithinRange

   subroutine get_boundary_conditions_stencil()
      implicit none

      wBC_bottom = -1
      wBC_top = -1  

      !! Bottom wall 
      call message(0,"Bottom Wall Boundary Condition is:")
      select case (botWall)
      case(1)
         call message(1,"No-Slip Wall")
         ! NOTE: no-slip wall requires both w = 0 and dwdz = 0. Therefore, w
         ! is an even extension, which also satisfies w = 0.
         uBC_bottom = 0
         vBC_bottom = 0
         wBC_bottom = 1
      case(2) 
         call message(1,"Slip Wall")
         uBC_bottom = 1
         vBC_bottom = 1 
      case(3) 
         call message(1,"Wall Model")
         uBC_bottom = 0
         vBC_bottom = 0
      case default
         call gracefulExit("Invalid choice for BOTTOM WALL BCs",423)
      end select

      !! Top wall 
      call message(0,"Top Wall Boundary Condition is:")
      select case (TopWall)
      case(1)
         call message(1,"No-Slip Wall")
         ! NOTE: no-slip wall requires both w = 0 and dwdz = 0. Therefore, w
         ! is an even extension, which also satisfies w = 0.
         uBC_top = 0
         vBC_top = 0
         wBC_top = 1
      case(2) 
         call message(1,"Slip Wall")
         uBC_top = 1 
         vBC_top = 1 
      case(3) 
         call message(1,"Wall Model")
         uBC_top = 0
         vBC_top = 0
      case default
         call gracefulExit("Invalid choice for TOP WALL BCs",13)
      end select

   end subroutine

   subroutine read_file(filename, field)
      implicit none
      character(*), intent(in) :: filename
      real(rkind), dimension(:,:,:), intent(out) :: field
      logical :: exists
      character(len=clen) :: file

      file = trim(inputdir)//'/'//trim(filename)
      inquire(file=trim(file), exist=exists)
      if(exists)then
         call message(1, 'Reading '//trim(file))
         call decomp_2d_read_one(1, field, trim(file), gpC)
      else
         call gracefulExit('Not found: '//trim(file), 915)
      end if
   end subroutine

   function getPattern(rid, budgetid, termid, key, stamp, isBase)
      implicit none
      character(len=clen) :: getPattern
      integer, intent(in) :: rid, budgetid, termid
      character(len=1) :: cbdgtid
      character(len=2) :: ctermid, crid
      character(len=*), optional :: key, stamp
      logical, optional :: isBase

      write(crid, '(I2.2)') rid
      write(cbdgtid, '(I1)') budgetid
      write(ctermid, '(I2.2)') termid

      getPattern = 'Run'//trim(crid)//'_comp_deficit_budget'//cbdgtid//'_term'//ctermid
      if(present(isBase))then
         if(isBase)then
            getPattern = 'Run'//trim(crid)//'_budget'//cbdgtid//'_term'//ctermid
         end if
      end if      

      if (present(key))then
         getPattern = trim(getPattern)//'_t'//trim(key)
      else
         getPattern = trim(getPattern)//'_t*'
      end if

      if (present(stamp))then
         getPattern = trim(getPattern)//'_n'//trim(stamp)
      else
         getPattern = trim(getPattern)//'_n~'
      end if

      getPattern = trim(getPattern)//'.s3D'
   end function

   ! Small helpers using ISO_C_BINDING to getpid()
   function getpid() result(pid)
      use iso_c_binding, only: c_int
      implicit none
      integer :: pid
      interface
         function c_getpid() bind(C, name="getpid") result(c_pid)
         import :: c_int
         integer(c_int) :: c_pid
         end function c_getpid
      end interface
      pid = int(c_getpid(), kind(pid))
   end function getpid

   pure function to_string(i) result(str)
      integer, intent(in) :: i
      character(len=32) :: str
      write(str, '(I0)') i
   end function to_string

   ! String utility functions
   logical pure function starts_with(s, pre) result(ok)
      character(*), intent(in) :: s, pre
      integer :: lp
      lp = len_trim(pre)
      if (lp == 0) then
         ok = .true.
      else
         ok = (len_trim(s) >= lp) .and. (s(1:lp) == pre(1:lp))
      end if
   end function starts_with

   logical pure function ends_with(s, suf) result(ok)
      character(*), intent(in) :: s, suf
      integer :: ls, ts
      ls = len_trim(suf); ts = len_trim(s)
      if (ls == 0) then
         ok = .true.
      else
         ok = (ts >= ls) .and. (s(ts-ls+1:ts) == suf(1:ls))
      end if
   end function ends_with

   ! Check if VAL is in LIST
   logical pure function in_list(list, n, val) result(found)
      character(len=*), intent(in) :: list(:)
      integer,          intent(in) :: n
      character(len=*), intent(in) :: val
      integer :: i
      found = .false.
      do i = 1, n
         if (list(i) == val) then
         found = .true.; return
         end if
      end do
   end function in_list

   subroutine get_keys_stamps()
      implicit none
      character(len=:), allocatable :: keys(:), stamps(:)
      character(len=clen) :: pattern
      integer :: k

      pattern = getPattern(rid, 0, 1)
      call message(0, 'Extracting time stamps with a pattern: '//trim(pattern))

      call list_matching_keys_budget(trim(inputdir), trim(pattern), keys, stamps)
      call sort_keys_and_stamps_numeric(keys, stamps, sorted_keys, sorted_stamps)

      call message(0, 'Found time stamps are: ')
      do k=1, size(sorted_keys)
         if(TimeWithinRange(trim(sorted_keys(k)), startIDX, endIDX))then
            call message(1, 'Time: '//trim(sorted_keys(k))//', # Frames: '//trim(sorted_stamps(k))//'  (within range)')
         else
            call message(1, 'Time: '//trim(sorted_keys(k))//', # Frames: '//trim(sorted_stamps(k))//'  (out of range)')
      end if
      end do
      call message(0, '  ')
  end subroutine

   subroutine sort_keys_and_stamps_numeric(keys, stamps, sorted_keys, sorted_stamps)
      !! Sort KEYS (time stamps) by their integer value (ascending),
      !! and apply the same ordering to STAMPS.
      !!
      !! Input:
      !!   keys(:)   - character time stamps, e.g. "000900", "001050"
      !!   stamps(:) - corresponding "~" stamps, e.g. "123456", "654321"
      !!
      !! Output (allocatable):
      !!   sorted_keys(:), sorted_stamps(:) - reordered copies
      !!
      character(len=*), intent(in)  :: keys(:)
      character(len=*), intent(in)  :: stamps(:)
      character(len=:), allocatable, intent(out) :: sorted_keys(:)
      character(len=:), allocatable, intent(out) :: sorted_stamps(:)

      integer :: n, i, j, ios, val
      integer, allocatable :: vals(:), idx(:)
      integer :: maxlen_k, maxlen_s
      character(len=:), allocatable :: s

      ! Basic checks
      n = size(keys)
      if (n == 0 .or. size(stamps) /= n) then
         allocate(character(len=1) :: sorted_keys(0))
         allocate(character(len=1) :: sorted_stamps(0))
         return
      end if

      allocate(vals(n), idx(n))

      ! Parse integers from KEYS; non-numeric => sent to the end
      do i = 1, n
         s = trim(keys(i))
         read(s, *, iostat=ios) val
         if (ios == 0) then
         vals(i) = val
         else
         vals(i) = huge(1)    ! put non-numeric keys after numeric ones
         end if
         idx(i) = i
      end do

      ! Simple O(n^2) indirect sort of idx by vals
      do i = 1, n-1
         do j = i+1, n
         if (vals(idx(j)) < vals(idx(i))) then
            call swap(idx(i), idx(j))   ! your existing swap(int,int)
         end if
         end do
      end do

      ! Decide output lengths
      maxlen_k = 0
      maxlen_s = 0
      do i = 1, n
         maxlen_k = max(maxlen_k, len_trim(keys(i)))
         maxlen_s = max(maxlen_s, len_trim(stamps(i)))
      end do
      if (maxlen_k <= 0) maxlen_k = 1
      if (maxlen_s <= 0) maxlen_s = 1

      ! Allocate outputs with trimmed lengths
      allocate(character(len=maxlen_k) :: sorted_keys(n))
      allocate(character(len=maxlen_s) :: sorted_stamps(n))

      ! Fill outputs according to permutation idx
      do i = 1, n
         sorted_keys(i)   = adjustl(keys(idx(i))(1:maxlen_k))
         sorted_stamps(i) = adjustl(stamps(idx(i))(1:maxlen_s))
      end do

      deallocate(vals, idx)

   end subroutine sort_keys_and_stamps_numeric

   pure subroutine swap(a, b)
      integer, intent(inout) :: a, b
      integer :: t
      t = a; a = b; b = t
   end subroutine swap

   ! Split pattern with one '*' into prefix and suffix
  subroutine split_one_star(pattern, prefix, suffix, ok)
      character(*), intent(in)  :: pattern
      character(len=:), allocatable, intent(out) :: prefix, suffix
      logical, intent(out) :: ok
      integer :: p, q, n
      n = len_trim(pattern)
      p = index(pattern(:n), '*')
      if (p == 0) then
         ok = .false.; prefix = ''; suffix = ''; return
      end if
      q = index(pattern(p+1:n), '*')
      if (q /= 0) then
         ok = .false.; prefix = ''; suffix = ''; return
      end if
      prefix = pattern(:p-1)
      suffix = pattern(p+1:n)
      ok = .true.
   end subroutine split_one_star

   subroutine list_matching_keys_budget(dir, pattern, keys, stamps)
      ! To handle files like:
      !   Run06_budget0_term13_t*_n~.s3D
      !
      ! where:
      !   *  -> time stamp (returned in KEYS)
      !   ~  -> 6-digit stamp (returned in STAMPS)
      !
      ! Example filenames:
      !   Run06_budget0_term13_t000900_n123456.s3D
      !   Run06_budget0_term13_t001050_n654321.s3D
      !
      ! Result:
      !   keys   = ["000900","001050",...]
      !   stamps = ["123456","654321",...]
      !
      character(*), intent(in) :: dir
      character(*), intent(in) :: pattern
      character(len=:), allocatable, intent(out) :: keys(:)
      character(len=:), allocatable, intent(out) :: stamps(:)

      character(len=:), allocatable :: pre, suf
      character(len=:), allocatable :: d_esc, p_glob, tmpfile, cmd
      character(len=4096) :: line
      integer :: istat, u, nlines, maxlen_k, maxlen_s, klen
      integer :: ts, lp, pos_n, extpos
      logical :: ok, ex

      ! Default empty result
      allocate(keys(0),   mold='     ')
      allocate(stamps(0), mold='     ')

      ! Split pattern around the single '*' to get prefix PRE (up to 't')
      call split_one_star(pattern, pre, suf, ok)
      if (.not. ok) then
         ! either no '*' or more than one '*'
         return
      end if

      ! Escape directory name
      d_esc  = escape_single_quotes(trim(dir))

      ! Build a glob pattern for 'find':
      !   original:  Run06_budget0_term13_t*_n~.s3D
      !   glob:      Run06_budget0_term13_t*_n*.s3D
      !
      ! i.e. replace '~' with '*' so we ignore the 6-digit stamp in the shell.
      block
         integer :: i, L
         character(len=:), allocatable :: tmp
         L = len_trim(pattern)
         allocate(character(len=L) :: tmp)
         tmp = pattern
         do i = 1, L
         if (tmp(i:i) == '~') tmp(i:i) = '*'
         end do
         p_glob = escape_single_quotes(trim(tmp))
      end block

      tmpfile = '/tmp/fortran_glob_'//to_string(getpid())//'_keys.txt'

      cmd = "find '"//d_esc//"' -maxdepth 1 -type f -name '"//p_glob// &
            "' -printf '%f\n' > '"//tmpfile//"' 2>/dev/null"
      call execute_command_line(cmd, exitstat=istat)
      if (istat /= 0) return

      inquire(file=tmpfile, exist=ex); if (.not. ex) return

      ! Count matches first
      nlines = 0
      open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
      if (istat /= 0) return
      do
         read(u,'(A)', iostat=istat) line
         if (istat /= 0) exit
         nlines = nlines + 1
      end do
      close(u)

      if (nlines == 0) then
         call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
         return
      end if

      ! Temp store (over-allocated), we'll dedupe then shrink
      if (allocated(keys))   deallocate(keys)
      if (allocated(stamps)) deallocate(stamps)
      allocate(character(len=clen) :: keys(nlines))
      allocate(character(len=clen) :: stamps(nlines))
      klen      = 0
      maxlen_k  = 0
      maxlen_s  = 0

      open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
      if (istat /= 0) then
         deallocate(keys);   allocate(keys(0),   mold='     ')
         deallocate(stamps); allocate(stamps(0), mold='     ')
         call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
         return
      end if

      lp = len_trim(pre)

      do
         read(u,'(A)', iostat=istat) line
         if (istat /= 0) exit
         ts = len_trim(line)
         if (ts <= 0) cycle

         ! Must start with PRE (e.g. "Run06_budget0_term13_t")
         if (.not. starts_with(line(:ts), pre)) cycle

         ! Find the "_n" that comes after the timestamp
         pos_n = index(line(:ts), '_n')
         if (pos_n <= 0) cycle   ! no "_n" -> not our file

         ! Check extension ".s3D"
         if (ts < 4) cycle
         extpos = ts - 3          ! position of '.' in ".s3D"
         if (line(extpos:ts) /= '.s3D') cycle

         ! Extract timestamp between PRE and "_n"
         if (pos_n <= lp+1) cycle   ! nothing between prefix and "_n"
         ! time stamp (*)
         block
         character(len=clen) :: tstamp, sstamp
         integer :: lt, ls

         tstamp = line(lp+1 : pos_n-1)

         ! Extract the 6-digit stamp (~) between "n" and ".s3D"
         ! line: "..._n123456.s3D"
         ! pos_n: index of "_"
         ! 'n' is pos_n+1, stamp starts at pos_n+2, ends at extpos-1
         if (extpos <= pos_n+2) cycle
         sstamp = line(pos_n+2 : extpos-1)

         ! Deduplicate based on time stamp; if same time stamp appears twice
         ! we'll ignore duplicates (assuming 1-to-1 as you said).
         if (.not. in_list(keys, klen, trim(tstamp))) then
            klen = klen + 1
            keys(klen)   = trim(tstamp)
            stamps(klen) = trim(sstamp)
            lt = len_trim(tstamp)
            ls = len_trim(sstamp)
            maxlen_k = max(maxlen_k, lt)
            maxlen_s = max(maxlen_s, ls)
         end if
         end block
      end do

      close(u)
      call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)

      ! Resize KEYS and STAMPS to exactly klen and appropriate lengths
      if (klen == 0) then
         deallocate(keys);   allocate(keys(0),   mold='     ')
         deallocate(stamps); allocate(stamps(0), mold='     ')
      else
         block
         character(len=:), allocatable :: tmpk(:), tmps(:)
         integer :: j

         allocate(character(len=maxlen_k) :: tmpk(klen))
         allocate(character(len=maxlen_s) :: tmps(klen))

         do j = 1, klen
            tmpk(j) = adjustl(keys(j)(:maxlen_k))
            tmps(j) = adjustl(stamps(j)(:maxlen_s))
         end do

         call move_alloc(tmpk, keys)
         call move_alloc(tmps, stamps)
         end block
      end if
   end subroutine list_matching_keys_budget

   pure function escape_single_quotes(s) result(t)
      character(*), intent(in) :: s
      character(len=:), allocatable :: t
      integer :: i, n, extra, pos
      n = len_trim(s)
      extra = 0
      do i = 1, n
         if (s(i:i) == "'") extra = extra + 3  ! "'" -> '\'' (3 extra chars)
      end do
      t = repeat(' ', n + extra)
      pos = 1
      do i = 1, n
         if (s(i:i) == "'") then
         t(pos:pos) = "'"; pos = pos + 1
         t(pos:pos) = "\"; pos = pos + 1
         t(pos:pos) = "'"; pos = pos + 1
         t(pos:pos) = "'"; pos = pos + 1
         else
         t(pos:pos) = s(i:i); pos = pos + 1
         end if
      end do
      if (pos <= len(t)) t = t(:pos-1)
   end function escape_single_quotes

   subroutine ddx_R2R(f, dfdx)
        real(rkind), dimension(:,:,:), intent(in) :: f
        real(rkind), dimension(:,:,:), intent(out) :: dfdx

        call spectC%fft(f, cbuffyC)
        call spectC%mtimes_ik1_ip(cbuffyC)
        call spectC%dealias(cbuffyC)
        call spectC%ifft(cbuffyC, dfdx)
    end subroutine 

    subroutine ddy_R2R(f, dfdy)
        real(rkind), dimension(:,:,:), intent(in) :: f
        real(rkind), dimension(:,:,:), intent(out) :: dfdy

        call spectC%fft(f, cbuffyC)
        call spectC%mtimes_ik2_ip(cbuffyC)
        call spectC%dealias(cbuffyC)
        call spectC%ifft(cbuffyC, dfdy)
    end subroutine 

    subroutine ddz_R2R(f, dfdz, n1, n2)
        real(rkind), dimension(:,:,:), intent(in) :: f
        real(rkind), dimension(:,:,:), intent(out) :: dfdz
        integer, intent(in) :: n1, n2

        call spectC%fft(f, cbuffyC)
        call transpose_y_to_z(cbuffyC, cbuffzC(:,:,:,1), sp_gpC)
        call Pade6opZ%ddz_C2C(cbuffzC(:,:,:,1), cbuffzC(:,:,:,2), n1, n2)
        call transpose_z_to_y(cbuffzC(:,:,:,2), cbuffyC, sp_gpC)
        call spectC%dealias(cbuffyC)
        call spectC%ifft(cbuffyC, dfdz)
    end subroutine

   subroutine initializeEverything()
      implicit none
      integer :: ix1, iy1, iz1
      integer :: ixn, iyn, izn
      integer :: i,j,k

      ! Allocate memory
      call message(0,'Allocating memory ...')
      allocate(mesh(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 3))

      if(allocated(rbuffxC))deallocate(rbuffxC)
      allocate(rbuffxC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 8))
      du => rbuffxC(:,:,:,1)
      dv => rbuffxC(:,:,:,2)
      dw => rbuffxC(:,:,:,3)
      ubase => rbuffxC(:,:,:,4)
      vbase => rbuffxC(:,:,:,5)
      wbase => rbuffxC(:,:,:,6)
      bf => rbuffxC(:,:,:,7)
      buffer => rbuffxC(:,:,:,8)
      
      allocate(cbuffyC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3)))
      allocate(cbuffzC(sp_gpC%zsz(1),sp_gpC%zsz(2),sp_gpC%zsz(3),2)) 

      ! Create Mesh
      ix1 = gpC%xst(1); iy1 = gpC%xst(2); iz1 = gpC%xst(3)
      ixn = gpC%xen(1); iyn = gpC%xen(2); izn = gpC%xen(3)
      do k=1,size(mesh,3)
          do j=1,size(mesh,2)
              do i=1,size(mesh,1)
                  mesh(i,j,k,1) = real( ix1 + i - 1, rkind ) * dx
                  mesh(i,j,k,2) = real( iy1 + j - 1, rkind ) * dy
                  mesh(i,j,k,3) = real( iz1 + k - 1, rkind ) * dz + dz/two
              end do
          end do
      end do
      mesh(:,:,:,1) = mesh(:,:,:,1) - dx; mesh(:,:,:,2) = mesh(:,:,:,2) - dy; mesh(:,:,:,3) = mesh(:,:,:,3) - dz 
      call message(0,'Created mesh.')

      ! Initialize Padeder
      call Pade6opz%init(gpC, sp_gpC, gpE, sp_gpE, dz, NumericalSchemeVert, PeriodicInZ, spectC)
      call message(0,'Pade operations initialized')

      ! BCs for ddz
      call get_boundary_conditions_stencil()
      call message(0,'Identified boundary condition stenciles')

      ! Intersect the box with the mesh
      if(do_box_averaging)then
         call intersectBoxAndMesh()
         call message(0,'Control volume box intersected with the mesh')
      end if
      
      if(do_box_averaging) then
         if(do_x_budget)allocate(xprofiles(nx_box, 24))
         if(do_y_budget)allocate(yprofiles(nx_box, 24))
         if(do_z_budget)allocate(zprofiles(nx_box, 24))
         if(do_TKE_budget)allocate(TKEprofiles(nx_box, 31))
         if(do_MKE_budget)allocate(MKEprofiles(nx_box, 26))
         if(do_TMP_budget)allocate(TMPprofiles(nx_box, 3))
      end if
   end subroutine

   subroutine release_memory()
    implicit none

      if(allocated(mesh)) deallocate(mesh)
      if(allocated(xprofiles)) deallocate(xprofiles)
      if(allocated(yprofiles)) deallocate(yprofiles)
      if(allocated(zprofiles)) deallocate(zprofiles)
      if(allocated(TKEprofiles)) deallocate(TKEprofiles)
      if(allocated(MKEprofiles)) deallocate(MKEprofiles)
      if(allocated(TMPprofiles)) deallocate(TMPprofiles)
      if(allocated(xstations)) deallocate(xstations)
      if(allocated(rbuffxC)) deallocate(rbuffxC)
      nullify(du, dv, dw, ubase, vbase, wbase, bf, buffer)
      
      call spectC%destroy()
      call spectE%destroy()
      call Pade6opZ%destroy()
      call decomp_info_finalize(gpC)
      call decomp_info_finalize(gpE)
      call decomp_2d_finalize()
  end subroutine

  subroutine read_velocity(crid, cbrid, key, stamp)
   implicit none
   character(len=*), intent(in) :: crid, cbrid, key, stamp
   
   call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', du) !du
   call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', dv) !dv
   call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', dw) !dw
   call read_file('Run'//trim(cbrid)//'_budget0_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', ubase) !ubase
   call read_file('Run'//trim(cbrid)//'_budget0_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', vbase) !vbase
   call read_file('Run'//trim(cbrid)//'_budget0_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', wbase) !wbase
  end subroutine

  subroutine zbudget(crid, key, stamp)
   implicit none
   character(len=*), intent(in) :: crid, key, stamp
   character(len=2) :: idx_str
   integer :: idx

   do idx = 1, 24
      select case(idx)
      case(1)
         ! Advection: delta u_1 * partial_1 (delta u_3)
         call ddx_R2R(dw, bf); buffer = du * bf
      case(2)
         ! Advection: delta u_2 * partial_2 (delta u_3)
         call ddy_R2R(dw, bf); buffer = dv * bf
      case(3)
         ! Advection: delta u_3 * partial_3 (delta u_3)
         call ddz_R2R(dw, bf, wBC_bottom, wBC_top); buffer = dw * bf
      case(4)
         ! Advection: delta u_1 * partial_1 (base u_3)
         call ddx_R2R(wbase, bf); buffer = du * bf
      case(5)
         ! Advection: delta u_2 * partial_2 (base u_3)
         call ddy_R2R(wbase, bf); buffer = dv * bf
      case(6)
         ! Advection: delta u_3 * partial_3 (base u_3)
         call ddz_R2R(wbase, bf, wBC_bottom, wBC_top); buffer = dw * bf
      case(7)
         ! Advection: base u_1 * partial_1 (delta u_3)
         call ddx_R2R(dw, bf); buffer = ubase * bf
      case(8)
         ! Advection: base u_2 * partial_2 (delta u_3)
         call ddy_R2R(dw, bf); buffer = vbase * bf
      case(9)
         ! Advection: base u_3 * partial_3 (delta u_3)
         call ddz_R2R(dw, bf, wBC_bottom, wBC_top); buffer = wbase * bf
      case(10)
         ! pressure gradient: partial_2 (delta p)
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term20_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
       case(11)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_3' delta u_j')
         ! partial_j mean(delta u_3' delta u_j') = mean(delta u_j' partial_j delta u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(12)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_3' base u_j')
         ! partial_j mean(delta u_3' base u_j') = mean(base u_j' partial_j delta u_3')
        call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(13)
         ! Divergence of Reynolds stresses: partial_j mean(base u_3' delta u_j')
         ! partial_j mean(base u_3' delta u_j') = mean(delta u_j' partial_j base u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(14)
         ! w_sgs
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(15)
         ! wbuoyancy
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term17_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         buffer = - buffer ! Put it on the LHS
      case(16)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_3' delta u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(17)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_3' delta u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(18)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_3' delta u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, 1, 1)
      case(19)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_3' base u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(20)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_3' base u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(21)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_3' base u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, 1, 1)
      case(22)
         ! Divergence of Reynolds stresses: partial_1 mean(base u_3' delta u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(23)
         ! Divergence of Reynolds stresses: partial_2 mean(base u_3' delta u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(24)
         ! Divergence of Reynolds stresses: partial_3 mean(base u_3' delta u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, 1, 1)
      end select

      ! Average this budget term across the box
      if(do_box_averaging) call integrate_box_yz(buffer, zprofiles(:,idx))

      ! Write to file calculated dependent variables if requested
      if(writeDependentVariables .and. depedent_variable(3, idx))then
         write(idx_str, '(I2.2)') idx
         call dump_budget_field(buffer, idx_str, '7', trim(key), trim(stamp))
      end if
   end do

  end subroutine

  subroutine ybudget(crid, key, stamp)
   implicit none
   character(len=*), intent(in) :: crid, key, stamp
   character(len=2) :: idx_str
   integer :: idx

   do idx = 1, 24
      select case(idx)
      case(1)
         ! Advection: delta u_1 * partial_1 (delta u_2)
         call ddx_R2R(dv, bf); buffer = du * bf
      case(2)
         ! Advection: delta u_2 * partial_2 (delta u_2)
         call ddy_R2R(dv, bf); buffer = dv * bf
      case(3)
         ! Advection: delta u_3 * partial_3 (delta u_2)
         call ddz_R2R(dv, bf, vBC_bottom, vBC_top); buffer = dw * bf
      case(4)
         ! Advection: delta u_1 * partial_1 (base u_2)
         call ddx_R2R(vbase, bf); buffer = du * bf
      case(5)
         ! Advection: delta u_2 * partial_2 (base u_2)
         call ddy_R2R(vbase, bf); buffer = dv * bf
      case(6)
         ! Advection: delta u_3 * partial_3 (base u_2)
         call ddz_R2R(vbase, bf, vBC_bottom, vBC_top); buffer = dw * bf
      case(7)
         ! Advection: base u_1 * partial_1 (delta u_2)
         call ddx_R2R(dv, bf); buffer = ubase * bf
      case(8)
         ! Advection: base u_2 * partial_2 (delta u_2)
         call ddy_R2R(dv, bf); buffer = vbase * bf
      case(9)
         ! Advection: base u_3 * partial_3 (delta u_2)
         call ddz_R2R(dv, bf, vBC_bottom, vBC_top); buffer = wbase * bf
      case(10)
         ! pressure gradient: partial_2 (delta p)
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term19_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
       case(11)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_2' delta u_j')
         ! partial_j mean(delta u_2' delta u_j') = mean(delta u_j' partial_j delta u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(12)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_2' base u_j')
         ! partial_j mean(delta u_2' base u_j') = mean(base u_j' partial_j delta u_2')
        call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(13)
         ! Divergence of Reynolds stresses: partial_j mean(base u_2' delta u_j')
         ! partial_j mean(base u_2' delta u_j') = mean(delta u_j' partial_j base u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(14)
         ! v_sgs
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(15)
         ! v_cor
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term16_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         buffer = - buffer ! Put it on the LHS
      case(16)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_2' delta u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(17)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_2' delta u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(18)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_2' delta u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, -1, -1)
      case(19)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_2' base u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(20)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_2' base u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(21)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_2' base u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, -1, -1)
      case(22)
         ! Divergence of Reynolds stresses: partial_1 mean(base u_2' delta u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(23)
         ! Divergence of Reynolds stresses: partial_2 mean(base u_2' delta u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(24)
         ! Divergence of Reynolds stresses: partial_3 mean(base u_2' delta u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, -1, -1)
      end select

      ! Average this budget term across the box
      if(do_box_averaging) call integrate_box_yz(buffer, yprofiles(:,idx))

      ! Write to file calculated dependent variables if requested
      if(writeDependentVariables .and. depedent_variable(2, idx))then
         write(idx_str, '(I2.2)') idx
         call dump_budget_field(buffer, idx_str, '6', trim(key), trim(stamp))
      end if
   end do

  end subroutine

  subroutine xbudget(crid, key, stamp)
   implicit none
   character(len=*), intent(in) :: crid, key, stamp
   character(len=2) :: idx_str
   integer :: idx

   do idx = 1, 24
      select case(idx)
      case(1)
         ! Advection: delta u_1 * partial_1 (delta u_1)
         call ddx_R2R(du, bf); buffer = du * bf
      case(2)
         ! Advection: delta u_2 * partial_2 (delta u_1)
         call ddy_R2R(du, bf); buffer = dv * bf
      case(3)
         ! Advection: delta u_3 * partial_3 (delta u_1)
         call ddz_R2R(du, bf, uBC_bottom, uBC_top); buffer = dw * bf
      case(4)
         ! Advection: delta u_1 * partial_1 (base u_1)
         call ddx_R2R(ubase, bf); buffer = du * bf
      case(5)
         ! Advection: delta u_2 * partial_2 (base u_1)
         call ddy_R2R(ubase, bf); buffer = dv * bf
      case(6)
         ! Advection: delta u_3 * partial_3 (base u_1)
         call ddz_R2R(ubase, bf, uBC_bottom, uBC_top); buffer = dw * bf
      case(7)
         ! Advection: base u_1 * partial_1 (delta u_1)
         call ddx_R2R(du, bf); buffer = ubase * bf
      case(8)
         ! Advection: base u_2 * partial_2 (delta u_1)
         call ddy_R2R(du, bf); buffer = vbase * bf
      case(9)
         ! Advection: base u_3 * partial_3 (delta u_1)
         call ddz_R2R(du, bf, uBC_bottom, uBC_top); buffer = wbase * bf
      case(10)
         ! pressure gradient: partial_1 (delta p)
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term18_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
       case(11)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_1' delta u_j')
         ! partial_j mean(delta u_1' delta u_j') = mean(delta u_j' partial_j delta u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(12)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_1' base u_j')
         ! partial_j mean(delta u_1' base u_j') = mean(base u_j' partial_j delta u_1')
        call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(13)
         ! Divergence of Reynolds stresses: partial_j mean(base u_1' delta u_j')
         ! partial_j mean(base u_1' delta u_j') = mean(delta u_j' partial_j base u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(14)
         ! u_sgs
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      case(15)
         ! u_cor
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         buffer = - buffer ! Put it on the LHS
      case(16)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_1' delta u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(17)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_1' delta u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(18)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_1' delta u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, -1, -1)
      case(19)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_1' base u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(20)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_1' base u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(21)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_1' base u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, -1, -1)
      case(22)
         ! Divergence of Reynolds stresses: partial_1 mean(base u_1' delta u_1')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(bf, buffer)
      case(23)
         ! Divergence of Reynolds stresses: partial_2 mean(base u_1' delta u_2')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(bf, buffer)
      case(24)
         ! Divergence of Reynolds stresses: partial_3 mean(base u_1' delta u_3')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(bf, buffer, -1, -1)
      end select

      ! Average this budget term across the box
      if(do_box_averaging) call integrate_box_yz(buffer, xprofiles(:,idx))

      ! Write to file calculated dependent variables if requested
      if(writeDependentVariables .and. depedent_variable(1, idx))then
         write(idx_str, '(I2.2)') idx
         call dump_budget_field(buffer, idx_str, '5', trim(key), trim(stamp))
      end if
   end do

  end subroutine

  subroutine TKEbudget(crid, cbrid, key, stamp)
   implicit none
   character(len=*), intent(in) :: crid, cbrid, key, stamp
   character(len=2) :: idx_str
   integer :: idx
   real(rkind), dimension(:,:,:), allocatable :: bf1

   allocate(bf1(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   do idx = 1, 31
      select case(idx)
      case(1)
         ! Advection: delta u_j * partial_j (delta u_i' delta u_i')/2 
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         bf = half*buffer; buffer=zero
         call ddx_R2R(bf, bf1); buffer = buffer + bf1*du
         call ddy_R2R(bf, bf1); buffer = buffer + bf1*dv
         call ddz_R2R(bf, bf1, 1, 1); buffer = buffer + bf1*dw ! even

      case(2)
         ! Advection: delta u_j * partial_j (delta u_i' base u_i') 
         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         buffer=zero
         call ddx_R2R(bf1, bf); buffer = buffer + bf*du
         call ddy_R2R(bf1, bf); buffer = buffer + bf*dv
         call ddz_R2R(bf1, bf, 1, 1); buffer = buffer + bf*dw ! even

      case(3)
         ! Advection: delta u_j * partial_j (base u_i' base u_i')/2 
         buffer=zero
         call read_file('Run'//trim(cbrid)//'_budget0_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         call read_file('Run'//trim(cbrid)//'_budget0_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         call read_file('Run'//trim(cbrid)//'_budget0_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         bf = half*buffer; buffer=zero
         call ddx_R2R(bf, bf1); buffer = buffer + bf1*du
         call ddy_R2R(bf, bf1); buffer = buffer + bf1*dv
         call ddz_R2R(bf, bf1, 1, 1); buffer = buffer + bf1*dw ! even

      case(4)
         ! Advection: base u_j * partial_j (delta u_i' delta u_i')/2
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); buffer=buffer+bf
         bf = half*buffer; buffer=zero
         call ddx_R2R(bf, bf1); buffer = buffer + bf1*ubase
         call ddy_R2R(bf, bf1); buffer = buffer + bf1*vbase
         call ddz_R2R(bf, bf1, 1, 1); buffer = buffer + bf1*wbase ! even

      case(5)
         ! Advection: base u_j * partial_j (delta u_i' base u_i') 
         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         buffer=zero
         call ddx_R2R(bf1, bf); buffer = buffer + bf*ubase
         call ddy_R2R(bf1, bf); buffer = buffer + bf*vbase
         call ddz_R2R(bf1, bf, 1, 1); buffer = buffer + bf*wbase ! even

      case(6)
         ! Production: mean(delta u_i' delta u_j') partial_j mean(delta u_i)
         buffer=zero

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(du, bf1, uBC_bottom, uBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dv, bf1, vBC_bottom, vBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dw, bf1, wBC_bottom, wBC_top); buffer = buffer + bf * bf1
         
      case(7)
         ! Production: mean(delta u_i' base u_j') partial_j mean(delta u_i)
         buffer=zero

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(du, bf1, uBC_bottom, uBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dv, bf1, vBC_bottom, vBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dw, bf1, wBC_bottom, wBC_top); buffer = buffer + bf * bf1

      case(8)
         ! Production: mean(base u_i' delta u_j') partial_j mean(delta u_i)
         buffer=zero

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(du, bf1, uBC_bottom, uBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dv, bf1, vBC_bottom, vBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dw, bf1, wBC_bottom, wBC_top); buffer = buffer + bf * bf1

      case(9)
         ! Production: mean(base u_i' base u_j') partial_j mean(delta u_i)
         buffer=zero

         call read_file('Run'//trim(cbrid)//'_budget0_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(cbrid)//'_budget0_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(du, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(cbrid)//'_budget0_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(du, bf1, uBC_bottom, uBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(cbrid)//'_budget0_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(cbrid)//'_budget0_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dv, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(cbrid)//'_budget0_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dv, bf1, vBC_bottom, vBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(cbrid)//'_budget0_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(cbrid)//'_budget0_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(dw, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(cbrid)//'_budget0_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(dw, bf1, wBC_bottom, wBC_top); buffer = buffer + bf * bf1

      case(10)
         ! Production: mean(delta u_i' delta u_j') partial_j mean(base u_i)
         buffer=zero

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(ubase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(ubase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(ubase, bf1, uBC_bottom, uBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(vbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(vbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(vbase, bf1, vBC_bottom, vBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(wbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(wbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(wbase, bf1, wBC_bottom, wBC_top); buffer = buffer + bf * bf1

      case(11)
         ! Production: mean(delta u_i' base u_j') partial_j mean(base u_i)
         buffer=zero

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(ubase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(ubase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(ubase, bf1, uBC_bottom, uBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(vbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(vbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(vbase, bf1, vBC_bottom, vBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(wbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(wbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(wbase, bf1, wBC_bottom, wBC_top); buffer = buffer + bf * bf1

      case(12)
         ! Production: mean(base u_i' delta u_j') partial_j mean(base u_i)
         buffer=zero

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(ubase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(ubase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(ubase, bf1, uBC_bottom, uBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(vbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(vbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(vbase, bf1, vBC_bottom, vBC_top); buffer = buffer + bf * bf1

         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddx_R2R(wbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddy_R2R(wbase, bf1); buffer = buffer + bf * bf1
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         call ddz_R2R(wbase, bf1, wBC_bottom, wBC_top); buffer = buffer + bf * bf1

      case(13)
         ! Buoyancy: mean(delta w' delta wb')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = - bf

      case(14)
         ! Buoyancy: mean(delta w' base wb')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = - bf

      case(15)
         ! Buoyancy covariance: mean(base w' delta wb')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = - bf

      case(16)
         ! Pressure covariance: mean(delta u_j' partial_j delta p')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         
      case(17)
         ! Pressure covariance: mean(base u_j' partial_j delta p')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         
      case(18)
         ! Pressure covariance: mean(delta u_j' partial_j base p')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         
      case(19)
         ! Transport: mean(delta u_i' delta u_j' partial_j delta u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term19_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         
      case(20)
         ! Transport: mean(delta u_i' base u_j' partial_j delta u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term18_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)

      case(21)
         ! Transport: mean(delta u_i' delta u_j' partial_j base u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term17_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
      
      case(22)
         ! Transport: mean(base u_i' delta u_j' partial_j delta u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term16_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         
      case(23)
         ! Transport: mean(delta u_i' base u_j' partial_j base u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)

      case(24)
         ! Transport: mean(base u_i' base u_j' partial_j delta u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)

      case(25)
         ! Transport: mean(base u_i' delta u_j' partial_j base u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         
      case(26)
         ! SGS transport: partial_j mean(base u_i' delta tau_ij')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
        
      case(27)
         ! SGS transport: partial_j mean(delta u_i' base tau_ij')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)

      case(28)
         ! SGS transport: partial_j mean(delta u_i' delta tau_ij')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', buffer)
         
      case(29)
         ! SGS Dissipation: mean(delta tau_ij' partial_j base u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = -bf

      case(30)
         ! SGS Dissipation: mean(base tau_ij' partial_j delta u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = -bf
         
      case(31)
         ! SGS Dissipation: mean(delta tau_ij' partial_j delta u_i')
         call read_file('Run'//trim(crid)//'_comp_deficit_budget3_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = -bf
      end select

      ! Average this budget term across the box
      if(do_box_averaging) call integrate_box_yz(buffer, TKEprofiles(:,idx))

      ! Write to file calculated dependent variables if requested
      if(writeDependentVariables .and. depedent_variable(4, idx))then
         write(idx_str, '(I2.2)') idx
         call dump_budget_field(buffer, idx_str, '4', trim(key), trim(stamp))
      end if
   end do
   deallocate(bf1)
  end subroutine

  subroutine MKEbudget(crid, cbrid, key, stamp)
   implicit none
   character(len=*), intent(in) :: crid, cbrid, key, stamp
   character(len=2) :: idx_str
   integer :: idx
   
   do idx = 1, 26
      select case(idx)
      case(1)
         ! Advection: delta u_i base u_j partial_j base u_i
         buffer=zero
         call ddx_R2R(ubase, bf); buffer = buffer + du * ubase * bf
         call ddy_R2R(ubase, bf); buffer = buffer + du * vbase * bf
         call ddz_R2R(ubase, bf, uBC_bottom, uBC_top); buffer = buffer + du * wbase * bf
         call ddx_R2R(vbase, bf); buffer = buffer + dv * ubase * bf
         call ddy_R2R(vbase, bf); buffer = buffer + dv * vbase * bf
         call ddz_R2R(vbase, bf, vBC_bottom, vBC_top); buffer = buffer + dv * wbase * bf
         call ddx_R2R(wbase, bf); buffer = buffer + dw * ubase * bf
         call ddy_R2R(wbase, bf); buffer = buffer + dw * vbase * bf
         call ddz_R2R(wbase, bf, wBC_bottom, wBC_top); buffer = buffer + dw * wbase * bf

      case(2)
         ! Advection: base u_i base u_j partial_j delta u_i
         buffer=zero
         call ddx_R2R(du, bf); buffer = buffer + ubase * ubase * bf
         call ddy_R2R(du, bf); buffer = buffer + ubase * vbase * bf
         call ddz_R2R(du, bf, uBC_bottom, uBC_top); buffer = buffer + ubase * wbase * bf
         call ddx_R2R(dv, bf); buffer = buffer + vbase * ubase * bf
         call ddy_R2R(dv, bf); buffer = buffer + vbase * vbase * bf
         call ddz_R2R(dv, bf, vBC_bottom, vBC_top); buffer = buffer + vbase * wbase * bf
         call ddx_R2R(dw, bf); buffer = buffer + wbase * ubase * bf
         call ddy_R2R(dw, bf); buffer = buffer + wbase * vbase * bf
         call ddz_R2R(dw, bf, wBC_bottom, wBC_top); buffer = buffer + wbase * wbase * bf

      case(3)
         ! Advection: delta u_i base u_j partial_j delta u_i
         buffer=zero
         call ddx_R2R(du, bf); buffer = buffer + du * ubase * bf
         call ddy_R2R(du, bf); buffer = buffer + du * vbase * bf
         call ddz_R2R(du, bf, uBC_bottom, uBC_top); buffer = buffer + du * wbase * bf
         call ddx_R2R(dv, bf); buffer = buffer + dv * ubase * bf
         call ddy_R2R(dv, bf); buffer = buffer + dv * vbase * bf
         call ddz_R2R(dv, bf, vBC_bottom, vBC_top); buffer = buffer + dv * wbase * bf
         call ddx_R2R(dw, bf); buffer = buffer + dw * ubase * bf
         call ddy_R2R(dw, bf); buffer = buffer + dw * vbase * bf
         call ddz_R2R(dw, bf, wBC_bottom, wBC_top); buffer = buffer + dw * wbase * bf

      case(4)
         ! Advection: base u_i delta u_j partial_j base u_i
         buffer=zero
         call ddx_R2R(ubase, bf); buffer = buffer + ubase * du * bf
         call ddy_R2R(ubase, bf); buffer = buffer + ubase * dv * bf
         call ddz_R2R(ubase, bf, uBC_bottom, uBC_top); buffer = buffer + ubase * dw * bf
         call ddx_R2R(vbase, bf); buffer = buffer + vbase * du * bf
         call ddy_R2R(vbase, bf); buffer = buffer + vbase * dv * bf
         call ddz_R2R(vbase, bf, vBC_bottom, vBC_top); buffer = buffer + vbase * dw * bf
         call ddx_R2R(wbase, bf); buffer = buffer + wbase * du * bf
         call ddy_R2R(wbase, bf); buffer = buffer + wbase * dv * bf
         call ddz_R2R(wbase, bf, wBC_bottom, wBC_top); buffer = buffer + wbase * dw * bf

      case(5)
         ! Advection: delta u_i delta u_j partial_j base u_i
         buffer=zero
         call ddx_R2R(ubase, bf); buffer = buffer + du * du * bf
         call ddy_R2R(ubase, bf); buffer = buffer + du * dv * bf
         call ddz_R2R(ubase, bf, uBC_bottom, uBC_top); buffer = buffer + du * dw * bf
         call ddx_R2R(vbase, bf); buffer = buffer + dv * du * bf
         call ddy_R2R(vbase, bf); buffer = buffer + dv * dv * bf
         call ddz_R2R(vbase, bf, vBC_bottom, vBC_top); buffer = buffer + dv * dw * bf
         call ddx_R2R(wbase, bf); buffer = buffer + dw * du * bf
         call ddy_R2R(wbase, bf); buffer = buffer + dw * dv * bf
         call ddz_R2R(wbase, bf, wBC_bottom, wBC_top); buffer = buffer + dw * dw * bf

      case(6)
         ! Advection: base u_i delta u_j partial_j delta u_i
         buffer=zero
         call ddx_R2R(du, bf); buffer = buffer + ubase * du * bf
         call ddy_R2R(du, bf); buffer = buffer + ubase * dv * bf
         call ddz_R2R(du, bf, uBC_bottom, uBC_top); buffer = buffer + ubase * dw * bf
         call ddx_R2R(dv, bf); buffer = buffer + vbase * du * bf
         call ddy_R2R(dv, bf); buffer = buffer + vbase * dv * bf
         call ddz_R2R(dv, bf, vBC_bottom, vBC_top); buffer = buffer + vbase * dw * bf
         call ddx_R2R(dw, bf); buffer = buffer + wbase * du * bf
         call ddy_R2R(dw, bf); buffer = buffer + wbase * dv * bf
         call ddz_R2R(dw, bf, wBC_bottom, wBC_top); buffer = buffer + wbase * dw * bf

      case(7)
         ! Advection: delta u_i delta u_j partial_j delta u_i
         buffer=zero
         call ddx_R2R(du, bf); buffer = buffer + du * du * bf
         call ddy_R2R(du, bf); buffer = buffer + du * dv * bf
         call ddz_R2R(du, bf, uBC_bottom, uBC_top); buffer = buffer + du * dw * bf
         call ddx_R2R(dv, bf); buffer = buffer + dv * du * bf
         call ddy_R2R(dv, bf); buffer = buffer + dv * dv * bf
         call ddz_R2R(dv, bf, vBC_bottom, vBC_top); buffer = buffer + dv * dw * bf
         call ddx_R2R(dw, bf); buffer = buffer + dw * du * bf
         call ddy_R2R(dw, bf); buffer = buffer + dw * dv * bf
         call ddz_R2R(dw, bf, wBC_bottom, wBC_top); buffer = buffer + dw * dw * bf

      case(8)
         ! Buoyancy: delta wb * delta w
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term17_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = - bf * dw

      case(9)
         ! Buoyancy: delta wb * base w
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term17_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = - bf * wbase

      case(10)
         ! Buoyancy: base wb * delta w
         ! Make sure that squeeze was .true. in the main simulation
         ! This is overloading an existing budget term in budget0.
         call read_file('Run'//trim(cbrid)//'_budget0_term25_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = - bf * dw

      case(11)
         ! Pressure gradient: delta u_i * d_i delta p
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term18_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term19_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term20_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf

      case(12)
         ! Pressure gradient: base u_i * d_i delta p
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term18_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + ubase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term19_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + vbase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term20_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + wbase * bf

      case(13)
         ! Pressure gradient: delta u_i * d_i base p
         ! Make sure that squeeze was .true. in the main simulation
         buffer=zero
         call read_file('Run'//trim(cbrid)//'_budget0_term17_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(cbrid)//'_budget0_term18_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(cbrid)//'_budget0_term19_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf

      case(14)
         ! SGS stresses: delta u_i * d_j delta tau_ij
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf

      case(15)
         ! SGS stresses: base u_i * d_j delta tau_ij
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + ubase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term13_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + vbase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term14_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + wbase * bf

      case(16)
         ! SGS stresses: delta u_i * d_j base tau_ij
         ! Make sure that squeeze was .true. in the main simulation
         buffer=zero
         call read_file('Run'//trim(cbrid)//'_budget0_term20_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(cbrid)//'_budget0_term21_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(cbrid)//'_budget0_term22_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf

      case(17)
         ! Production: delta u_i * mean(delta u_j' d_j delta u_i')
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf
         
      case(18)
         ! Production: delta u_i * mean(delta u_j' d_j base u_i')
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf

      case(19)
         ! Production: delta u_i * mean(base u_j' d_j delta u_i')
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf

      case(20)
         ! Production: delta u_i * mean(base u_j' d_j base u_i')
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term10_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + du * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term11_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dv * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + dw * bf

      case(21)
         ! Production: base u_i * mean(delta u_j' d_j delta u_i')
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + ubase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term02_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + vbase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term03_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + wbase * bf

      case(22)
         ! Production: base u_i * mean(delta u_j' d_j base u_i')
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + ubase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term05_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + vbase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + wbase * bf
         
      case(23)
         ! Production: base u_i * mean(base u_j' d_j delta u_i')
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + ubase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term08_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + vbase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget2_term09_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer=buffer + wbase * bf
         
      case(24)
         ! Coriolis: delta u_i * delta ucor_i
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = buffer - du * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term16_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = buffer - dv * bf
         
      case(25)
         ! Coriolis: delta u_i * base ucor_i
         ! Make sure that squeeze was .true. in the main simulation
         buffer=zero
         call read_file('Run'//trim(cbrid)//'_budget0_term23_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = buffer - du * bf
         call read_file('Run'//trim(cbrid)//'_budget0_term24_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = buffer - dv * bf

      case(26)
         ! Coriolis: base u_i * delta ucor_i
         buffer=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = buffer - ubase * bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget0_term16_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf)
         buffer = buffer - vbase * bf
      end select

      ! Average this budget term across the box
      if(do_box_averaging) call integrate_box_yz(buffer, MKEprofiles(:,idx))

      ! Write to file calculated dependent variables if requested
      if(writeDependentVariables .and. depedent_variable(5, idx))then
         write(idx_str, '(I2.2)') idx
         call dump_budget_field(buffer, idx_str, '8', trim(key), trim(stamp))
      end if
   end do

  end subroutine

  subroutine TMPbudget(crid, key, stamp)
   implicit none
   character(len=*), intent(in) :: crid, key, stamp
   character(len=2) :: idx_str
   integer :: idx
   real(rkind), dimension(:,:,:), allocatable :: bf1

   allocate(bf1(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   do idx = 1, 3
      select case(idx)
      case(1)
         ! Advection: base u_j * partial_j (delta u_i' delta u_i')/2
         buffer=zero
         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         bf1 = half*bf1
         call ddx_R2R(bf1, bf); buffer = buffer + bf * ubase

         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call ddx_R2R(bf1, bf); buffer = buffer + bf * ubase
         
      case(2)
         ! Advection: base u_j * partial_j (delta u_i' delta u_i')/2
         buffer=zero
         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         bf1 = half*bf1
         call ddy_R2R(bf1, bf); buffer = buffer + bf * vbase

         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call ddy_R2R(bf1, bf); buffer = buffer + bf * vbase
         
      case(3)
         ! Advection: base u_j * partial_j (delta u_i' delta u_i')/2
         buffer=zero
         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term01_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term04_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term06_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         bf1 = half*bf1
         call ddz_R2R(bf1, bf, 1, 1); buffer = buffer + bf * wbase

         bf1=zero
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term07_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term12_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call read_file('Run'//trim(crid)//'_comp_deficit_budget1_term15_t'//trim(key)//'_n'//trim(stamp)//'.s3D', bf); bf1=bf1+bf
         call ddz_R2R(bf1, bf, 1, 1); buffer = buffer + bf * wbase
      end select

      ! Average this budget term across the box
      if(do_box_averaging) call integrate_box_yz(buffer, TMPprofiles(:,idx))

      ! Write to file calculated dependent variables if requested
      if(writeDependentVariables .and. depedent_variable(6, idx))then
         write(idx_str, '(I2.2)') idx
         call dump_budget_field(buffer, idx_str, '9', trim(key), trim(stamp))
      end if
   end do
   deallocate(bf1)
  end subroutine

end module constructDeficitBudgets_mod

program constructDeficitBudgets
   use constructDeficitBudgets_mod

   implicit none
   integer :: ioUnit, ierr, k
   logical :: periodicbcs(3)
   character(len=clen) :: inputfile, ers
   character(len=2) :: crid, cbrid

   namelist /INPUT/ inputdir, outputdir, nx, ny, nz, Lx, Ly, Lz, prow, pcol, RID, &
                    BRID, writeDependentVariables, startIDX, endIDX, &
                    do_box_averaging, NumericalSchemeVert, &
                    PeriodicInZ, botWall, topWall, botBC_temp, &
                    do_x_budget, do_y_budget, do_z_budget, do_TKE_budget, do_MKE_budget, &
                    do_TMP_budget
   namelist /BOX/ x1, x2, y1, y2, z1, z2

   ! Do MPI stuff
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)

   ! Do file IO - input file
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', status='old', action='read')
   read(unit=ioUnit, NML=INPUT, IOSTAT=ierr)
   if (ierr/=0)then
      write(ers,'(I0)')ierr
      call gracefulExit("Reading failed for INPUT with error "//trim(ers), 101)
   end if
   read(unit=ioUnit, NML=BOX, IOSTAT=ierr)
   if (ierr/=0)then
      write(ers,'(I0)')ierr
      call gracefulExit("Reading failed for BOX with error "//trim(ers), 104)
   end if
   close(ioUnit)    

   periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = .false.
   call decomp_2d_init(nx, ny, nz, prow, pcol, periodicbcs)
   call get_decomp_info(gpC)
   call decomp_info_init(nx, ny, nz + 1, gpE)

   ! Initialize spectral
   dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
   call spectC%init("x",nx,ny,nz, dx, dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
   call spectE%init("x",nx,ny,nz + 1,dx,dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
   sp_gpC => spectC%spectdecomp
   sp_gpE => spectE%spectdecomp

   call initializeEverything()

   ! Get file list and sort by time  
   call get_keys_stamps()
   write(crid, '(I2.2)') RID
   write(cbrid, '(I2.2)') BRID

   ! Loop through time frames
   do k = 1, size(sorted_keys)
      call tic()

      if(.not. TimeWithinRange(trim(sorted_keys(k)), startIDX, endIDX)) cycle
      call message(0, 'Time Index: '//trim(sorted_keys(k))//', # Frames: '//trim(sorted_stamps(k)))
      call read_velocity(crid, cbrid, trim(sorted_keys(k)), trim(sorted_stamps(k)))

      if(do_x_budget)then 
         call xbudget(crid, trim(sorted_keys(k)), trim(sorted_stamps(k)))
         if((nrank == 0) .and. do_box_averaging)then
            call export_csv(csv_file_name(trim(sorted_keys(k)), trim(sorted_stamps(k)), "X"), xprofiles)
         end if
      end if

      if(do_y_budget)then 
         call ybudget(crid, trim(sorted_keys(k)), trim(sorted_stamps(k)))
         if((nrank == 0) .and. do_box_averaging)then
            call export_csv(csv_file_name(trim(sorted_keys(k)), trim(sorted_stamps(k)), "Y"), yprofiles)
         end if
      end if

      if(do_z_budget)then 
         call zbudget(crid, trim(sorted_keys(k)), trim(sorted_stamps(k)))
         if((nrank == 0) .and. do_box_averaging)then
            call export_csv(csv_file_name(trim(sorted_keys(k)), trim(sorted_stamps(k)), "Z"), zprofiles)
         end if
      end if

      if(do_TKE_budget)then 
         call TKEbudget(crid, cbrid, trim(sorted_keys(k)), trim(sorted_stamps(k)))
         if((nrank == 0) .and. do_box_averaging)then
            call export_csv(csv_file_name(trim(sorted_keys(k)), trim(sorted_stamps(k)), "TKE"), TKEprofiles)
         end if
      end if

      if(do_MKE_budget)then 
         call MKEbudget(crid, cbrid, trim(sorted_keys(k)), trim(sorted_stamps(k)))
         if((nrank == 0) .and. do_box_averaging)then
            call export_csv(csv_file_name(trim(sorted_keys(k)), trim(sorted_stamps(k)), "MKE"), MKEprofiles)
         end if
      end if

      if(do_TMP_budget)then 
         call TMPbudget(crid, trim(sorted_keys(k)), trim(sorted_stamps(k)))
         if((nrank == 0) .and. do_box_averaging)then
            call export_csv(csv_file_name(trim(sorted_keys(k)), trim(sorted_stamps(k)), "TMP"), TMPprofiles)
         end if
      end if

      call message(0, ' ')
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call toc()
   end do  

   call release_memory()  
   call MPI_FINALIZE(ierr) 

end program constructDeficitBudgets
