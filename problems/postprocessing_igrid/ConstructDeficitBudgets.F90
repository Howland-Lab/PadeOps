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

   external :: mpi_allreduce

   character(len=clen) :: inputdir, outputdir, tag='notag' 
   real(rkind) :: Lx = one, Ly = one, Lz = one
   integer :: botWall=3, topWall=2, botBC_temp=0
   logical :: PeriodicInZ=.false.
   type(spectral), target  :: spectE, spectC
   type(decomp_info) :: gpC, gpE
   type(decomp_info), pointer :: sp_gpC, sp_gpE
   type(Pade6stagg) :: Pade6opZ
   real(rkind) :: dx, dy, dz
   real(rkind), dimension(:,:), allocatable, target :: profiles
   real(rkind), dimension(:,:,:,:), allocatable, target :: mesh, Budget0, Budget1, Budget2, Budget3, duidxj, baseBudget0, duidxj_base
   real(rkind), dimension(:,:,:,:), allocatable, target :: rbuffxC
   complex(rkind), dimension(:,:,:), allocatable :: cbuffyC
   complex(rkind), dimension(:,:,:,:), allocatable, target :: cbuffzC
   integer :: prow=0, pcol=0, nx, ny, nz, RID, BRID, NumericalSchemeVert=1
   integer :: startIDX=-1, endIDX=999999
   integer :: uBC_bottom, uBC_top, vBC_bottom, vBC_top, wBC_bottom, wBC_top
   integer :: nx_box, ix1g, ix2g
   real(rkind) :: x1=zero, x2=zero, y1=zero, y2=zero, z1=zero, z2=zero
   integer :: num_profiles
   real(rkind), dimension(:), allocatable :: xstations
   logical :: writeDependentVariables = .false.
   integer :: budgettype=1 ! 1: x-Momentum, 2: y-Momentum, 3: z-Momentum, 4: TKE
   real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz
   real(rkind), dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
   real(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
   real(rkind), dimension(:,:,:), pointer :: dudx_base, dudy_base, dudz_base
   real(rkind), dimension(:,:,:), pointer :: dvdx_base, dvdy_base, dvdz_base
   real(rkind), dimension(:,:,:), pointer :: dwdx_base, dwdy_base, dwdz_base
   character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)
   logical :: do_box_averaging=.true.

   contains

   subroutine export_csv(key, stamp)
      implicit none
      character(len=*), intent(in) :: key, stamp
      character(clen) :: filename
      character(len=3) :: crid, tid
      integer :: i, j
      integer :: nx, ny
      integer :: unit
      character(len=3) :: name

      nx = size(profiles, 1)
      ny = size(profiles, 2)
      unit =1045

      select case(budgettype)
      case(1)
         name = 'X'
      case(2)
         name = 'Y'
      case(3)
         name = 'Z'
      case(4)
         name = 'TKE'
      end select

      write(crid, '(I2.2)') RID
      filename = trim(outputdir)//'/Run'//trim(crid)//'_t'//trim(key)//'_n'//trim(stamp)//'_'//trim(name)//'_Budgets_XProfile_'//trim(tag)//'.csv'

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

   subroutine compute_budgets(key, stamp)
      implicit none
      character(len=*), intent(in) :: key, stamp
      integer :: idx
      real(rkind), dimension(:,:,:), pointer :: buffer
      character(len=2) :: idx_str
      character(1) :: additional

      buffer => rbuffxC(:,:,:,3)

      do idx=1,num_profiles

         call message(1, 'Computing budget profile with index ', idx)

         select case(budgettype)
         case(1)
            call compute_X_budget_component(idx, buffer)
            additional = '5'
         case(2)
            call compute_Y_budget_component(idx, buffer)
            additional = '6'
         case(3)
            call compute_Z_budget_component(idx, buffer)
            additional = '7'
         case(4)
            call compute_TKE_budget_component(idx, buffer)
            additional = '4'
         end select

         ! Average this budget term across the box
         if(do_box_averaging) call integrate_box_yz(buffer, profiles(:,idx))

         ! Write to file calculated dependent variables if requested
         if(writeDependentVariables .and. depedent_variable(idx))then
            write(idx_str, '(I2.2)') idx
            call dump_budget_field(buffer, idx_str, additional, trim(key), trim(stamp))
         end if
      end do

      nullify(buffer)
   end subroutine

   function depedent_variable(idx)
      implicit none
      integer, intent(in) :: idx
      logical :: depedent_variable

      depedent_variable = .false.
      if((budgettype == 1) .or. (budgettype == 2) .or. (budgettype == 3))then
         ! X, Y, or Z momentum equation
         if((idx < 10) .or. (idx > 15)) depedent_variable = .true.
      elseif(budgettype == 4)then
         ! TKE equation
         if(idx <= 12) depedent_variable = .true.
      end if
   end function depedent_variable

   subroutine compute_X_budget_component(idx, buffer)
      implicit none
      integer, intent(in) :: idx
      real(rkind), dimension(:,:,:), intent(out) :: buffer
      real(rkind), dimension(:,:,:), pointer :: BF1, BF2

      BF1 => rbuffxC(:,:,:,1)
      BF2 => rbuffxC(:,:,:,2)

      buffer = zero
      select case(idx)
      case(1)
         ! Advection: delta u_1 * partial_1 (delta u_1)
         buffer = budget0(:,:,:,1) * dudx
      case(2)
         ! Advection: delta u_2 * partial_2 (delta u_1)
         buffer = budget0(:,:,:,2) * dudy
      case(3)
         ! Advection: delta u_3 * partial_3 (delta u_1)
         buffer = budget0(:,:,:,3) * dudz
      case(4)
         ! Advection: delta u_1 * partial_1 (base u_1)
         buffer = budget0(:,:,:,1) * dudx_base
      case(5)
         ! Advection: delta u_2 * partial_2 (base u_1)
         buffer = budget0(:,:,:,2) * dudy_base
      case(6)
         ! Advection: delta u_3 * partial_3 (base u_1)
         buffer = budget0(:,:,:,3) * dudz_base
      case(7)
         ! Advection: base u_1 * partial_1 (delta u_1)
         buffer = baseBudget0(:,:,:,1) * dudx
      case(8)
         ! Advection: base u_2 * partial_2 (delta u_1)
         buffer = baseBudget0(:,:,:,2) * dudy
      case(9)
         ! Advection: base u_3 * partial_3 (delta u_1)
         buffer = baseBudget0(:,:,:,3) * dudz
      case(10)
         ! pressure gradient: partial_1 (delta p)
         buffer = budget0(:,:,:,18)
      case(11)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_1' delta u_j')
         ! partial_j mean(delta u_1' delta u_j') = mean(delta u_j' partial_j delta u_1')
         buffer = budget2(:,:,:,1)
      case(12)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_1' base u_j')
         ! partial_j mean(delta u_1' base u_j') = mean(base u_j' partial_j delta u_1')
         buffer = budget2(:,:,:,7)
      case(13)
         ! Divergence of Reynolds stresses: partial_j mean(base u_1' delta u_j')
         ! partial_j mean(base u_1' delta u_j') = mean(delta u_j' partial_j base u_1')
         buffer = budget2(:,:,:,4)
      case(14)
         ! u_sgs
         buffer = budget0(:,:,:,12)
      case(15)
         ! u_cor
         buffer = budget0(:,:,:,15)
      case(16)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_1' delta u_1')
         call ddx_R2R(budget1(:,:,:,1), buffer)
      case(17)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_1' delta u_2')
         call ddy_R2R(budget1(:,:,:,2), buffer)
      case(18)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_1' delta u_3')
         call ddz_R2R(budget1(:,:,:,3), buffer, -1, -1) ! budget1(:,:,:,3) is odd
      case(19)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_1' base u_1')
         call ddx_R2R(budget1(:,:,:,7), buffer)
      case(20)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_1' base u_2')
         call ddy_R2R(budget1(:,:,:,8), buffer)
      case(21)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_1' base u_3')
         call ddz_R2R(budget1(:,:,:,10), buffer, -1, -1)
      case(22)
         ! Divergence of Reynolds stresses: partial_1 mean(base u_1' delta u_1')
         call ddx_R2R(budget1(:,:,:,7), buffer)
      case(23)
         ! Divergence of Reynolds stresses: partial_2 mean(base u_1' delta u_2')
         call ddy_R2R(budget1(:,:,:,9), buffer)
      case(24)
         ! Divergence of Reynolds stresses: partial_3 mean(base u_1' delta u_3')
         call ddz_R2R(budget1(:,:,:,11), buffer, -1, -1)
      end select
   end subroutine

   subroutine compute_Y_budget_component(idx, buffer)
      implicit none
      integer, intent(in) :: idx
      real(rkind), dimension(:,:,:), intent(out) :: buffer
      real(rkind), dimension(:,:,:), pointer :: BF1, BF2

      BF1 => rbuffxC(:,:,:,1)
      BF2 => rbuffxC(:,:,:,2)

      buffer = zero
      select case(idx)
      case(1)
         ! Advection: delta u_1 * partial_1 (delta u_2)
         buffer = budget0(:,:,:,1) * dvdx
      case(2)
         ! Advection: delta u_2 * partial_2 (delta u_2)
         buffer = budget0(:,:,:,2) * dvdy
      case(3)
         ! Advection: delta u_3 * partial_3 (delta u_2)
         buffer = budget0(:,:,:,3) * dvdz
      case(4)
         ! Advection: delta u_1 * partial_1 (base u_2)
         buffer = budget0(:,:,:,1) * dvdx_base
      case(5)
         ! Advection: delta u_2 * partial_2 (base u_2)
         buffer = budget0(:,:,:,2) * dvdy_base
      case(6)
         ! Advection: delta u_3 * partial_3 (base u_2)
         buffer = budget0(:,:,:,3) * dvdz_base
      case(7)
         ! Advection: base u_1 * partial_1 (delta u_2)
         buffer = baseBudget0(:,:,:,1) * dvdx
      case(8)
         ! Advection: base u_2 * partial_2 (delta u_2)
         buffer = baseBudget0(:,:,:,2) * dvdy
      case(9)
         ! Advection: base u_3 * partial_3 (delta u_2)
         buffer = baseBudget0(:,:,:,3) * dvdz
      case(10)
         ! pressure gradient: partial_2 (delta p)
         buffer = budget0(:,:,:,19)
      case(11)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_2' delta u_j')
         ! partial_j mean(delta u_2' delta u_j') = mean(delta u_j' partial_j delta u_2')
         buffer = budget2(:,:,:,2)
      case(12)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_2' base u_j')
         ! partial_j mean(delta u_2' base u_j') = mean(base u_j' partial_j delta u_2')
         buffer = budget2(:,:,:,8)
      case(13)
         ! Divergence of Reynolds stresses: partial_j mean(base u_2' delta u_j')
         ! partial_j mean(base u_2' delta u_j') = mean(delta u_j' partial_j base u_2')
         buffer = budget2(:,:,:,5)
      case(14)
         ! v_sgs
         buffer = budget0(:,:,:,13)
      case(15)
         ! v_cor
         buffer = budget0(:,:,:,16)
      case(16)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_2' delta u_1')
         call ddx_R2R(budget1(:,:,:,2), buffer)
      case(17)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_2' delta u_2')
         call ddy_R2R(budget1(:,:,:,4), buffer)
      case(18)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_2' delta u_3')
         call ddz_R2R(budget1(:,:,:,5), buffer, -1, -1) ! budget1(:,:,:,5) is odd
      case(19)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_2' base u_1')
         call ddx_R2R(budget1(:,:,:,9), buffer)
      case(20)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_2' base u_2')
         call ddy_R2R(budget1(:,:,:,12), buffer)
      case(21)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_2' base u_3')
         call ddz_R2R(budget1(:,:,:,13), buffer, -1, -1)
      case(22)
         ! Divergence of Reynolds stresses: partial_1 mean(base u_2' delta u_1')
         call ddx_R2R(budget1(:,:,:,8), buffer)
      case(23)
         ! Divergence of Reynolds stresses: partial_2 mean(base u_2' delta u_2')
         call ddy_R2R(budget1(:,:,:,12), buffer)
      case(24)
         ! Divergence of Reynolds stresses: partial_3 mean(base u_2' delta u_3')
         call ddz_R2R(budget1(:,:,:,14), buffer, -1, -1)
      end select
   end subroutine

   subroutine compute_Z_budget_component(idx, buffer)
      implicit none
      integer, intent(in) :: idx
      real(rkind), dimension(:,:,:), intent(out) :: buffer
      real(rkind), dimension(:,:,:), pointer :: BF1, BF2

      BF1 => rbuffxC(:,:,:,1)
      BF2 => rbuffxC(:,:,:,2)

      buffer = zero
      select case(idx)
      case(1)
         ! Advection: delta u_1 * partial_1 (delta u_3)
         buffer = budget0(:,:,:,1) * dwdx
      case(2)
         ! Advection: delta u_2 * partial_2 (delta u_3)
         buffer = budget0(:,:,:,2) * dwdy
      case(3)
         ! Advection: delta u_3 * partial_3 (delta u_3)
         buffer = budget0(:,:,:,3) * dwdz
      case(4)
         ! Advection: delta u_1 * partial_1 (base u_3)
         buffer = budget0(:,:,:,1) * dwdx_base
      case(5)
         ! Advection: delta u_2 * partial_2 (base u_3)
         buffer = budget0(:,:,:,2) * dwdy_base
      case(6)
         ! Advection: delta u_3 * partial_3 (base u_3)
         buffer = budget0(:,:,:,3) * dwdz_base
      case(7)
         ! Advection: base u_1 * partial_1 (delta u_3)
         buffer = baseBudget0(:,:,:,1) * dwdx
      case(8)
         ! Advection: base u_2 * partial_2 (delta u_3)
         buffer = baseBudget0(:,:,:,2) * dwdy
      case(9)
         ! Advection: base u_3 * partial_3 (delta u_3)
         buffer = baseBudget0(:,:,:,3) * dwdz
      case(10)
         ! pressure gradient: partial_3 (delta p)
         buffer = budget0(:,:,:,20)
      case(11)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_3' delta u_j')
         ! partial_j mean(delta u_3' delta u_j') = mean(delta u_j' partial_j delta u_3')
         buffer = budget2(:,:,:,3)
      case(12)
         ! Divergence of Reynolds stresses: partial_j mean(delta u_3' base u_j')
         ! partial_j mean(delta u_3' base u_j') = mean(base u_j' partial_j delta u_3')
         buffer = budget2(:,:,:,9)
      case(13)
         ! Divergence of Reynolds stresses: partial_j mean(base u_2' delta u_j')
         ! partial_j mean(base u_2' delta u_j') = mean(delta u_j' partial_j base u_2')
         buffer = budget2(:,:,:,6)
      case(14)
         ! w_sgs
         buffer = budget0(:,:,:,14)
      case(15)
         ! wb
         buffer = budget0(:,:,:,17)
      case(16)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_3' delta u_1')
         call ddx_R2R(budget1(:,:,:,3), buffer)
      case(17)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_3' delta u_2')
         call ddy_R2R(budget1(:,:,:,5), buffer)
      case(18)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_3' delta u_3')
         call ddz_R2R(budget1(:,:,:,6), buffer, -1, -1) ! budget1(:,:,:,6) is odd
      case(19)
         ! Divergence of Reynolds stresses: partial_1 mean(delta u_3' base u_1')
         call ddx_R2R(budget1(:,:,:,11), buffer)
      case(20)
         ! Divergence of Reynolds stresses: partial_2 mean(delta u_3' base u_2')
         call ddy_R2R(budget1(:,:,:,14), buffer)
      case(21)
         ! Divergence of Reynolds stresses: partial_3 mean(delta u_3' base u_3')
         call ddz_R2R(budget1(:,:,:,15), buffer, -1, -1)
      case(22)
         ! Divergence of Reynolds stresses: partial_1 mean(base u_3' delta u_1')
         call ddx_R2R(budget1(:,:,:,10), buffer)
      case(23)
         ! Divergence of Reynolds stresses: partial_2 mean(base u_3' delta u_2')
         call ddy_R2R(budget1(:,:,:,13), buffer)
      case(24)
         ! Divergence of Reynolds stresses: partial_3 mean(base u_3' delta u_3')
         call ddz_R2R(budget1(:,:,:,15), buffer, -1, -1)
      end select
   end subroutine

   subroutine compute_TKE_budget_component(idx, buffer)
      implicit none
      integer, intent(in) :: idx
      real(rkind), dimension(:,:,:), intent(out) :: buffer
      real(rkind), dimension(:,:,:), pointer :: BF1, BF2

      BF1 => rbuffxC(:,:,:,1)
      BF2 => rbuffxC(:,:,:,2)
      
      buffer = zero
      select case(idx)
      case(1)
         ! Advection: delta u_j * partial_j (delta u_i' delta u_i')/2 
         BF1 = half*(budget1(:,:,:,1) + budget1(:,:,:,4) + budget1(:,:,:,6))
         call ddx_R2R(BF1, BF2); buffer = buffer + BF2*budget0(:,:,:,1)
         call ddy_R2R(BF1, BF2); buffer = buffer + BF2*budget0(:,:,:,2)
         call ddz_R2R(BF1, BF2, 1, 1); buffer = buffer + BF2*budget0(:,:,:,3) ! BF1 is even
      
      case(2)
         ! Advection: delta u_j * partial_j (delta u_i' base u_i') 
         BF1 = (budget1(:,:,:,7) + budget1(:,:,:,12) + budget1(:,:,:,15))
         call ddx_R2R(BF1, BF2); buffer = buffer + BF2*budget0(:,:,:,1)
         call ddy_R2R(BF1, BF2); buffer = buffer + BF2*budget0(:,:,:,2)
         call ddz_R2R(BF1, BF2, 1, 1); buffer = buffer + BF2*budget0(:,:,:,3) ! BF1 is even
      
      case(3)
         ! Advection: delta u_j * partial_j (base u_i' base u_i')/2 
         BF1 = half*(baseBudget0(:,:,:,4) + baseBudget0(:,:,:,7) + baseBudget0(:,:,:,9))
         call ddx_R2R(BF1, BF2); buffer = buffer + BF2*budget0(:,:,:,1)
         call ddy_R2R(BF1, BF2); buffer = buffer + BF2*budget0(:,:,:,2)
         call ddz_R2R(BF1, BF2, 1, 1); buffer = buffer + BF2*budget0(:,:,:,3) ! BF1 is even
      
      case(4)
         ! Advection: base u_j * partial_j (delta u_i' delta u_i')/2
         BF1 = half*(budget1(:,:,:,1) + budget1(:,:,:,4) + budget1(:,:,:,6))
         call ddx_R2R(BF1, BF2); buffer = buffer + BF2*baseBudget0(:,:,:,1)
         call ddy_R2R(BF1, BF2); buffer = buffer + BF2*baseBudget0(:,:,:,2)
         call ddz_R2R(BF1, BF2, 1, 1); buffer = buffer + BF2*baseBudget0(:,:,:,3) ! BF1 is even
      
      case(5)
         ! Advection: base u_j * partial_j (delta u_i' base u_i') 
         BF1 = (budget1(:,:,:,7) + budget1(:,:,:,12) + budget1(:,:,:,15))
         call ddx_R2R(BF1, BF2); buffer = buffer + BF2*baseBudget0(:,:,:,1)
         call ddy_R2R(BF1, BF2); buffer = buffer + BF2*baseBudget0(:,:,:,2)
         call ddz_R2R(BF1, BF2, 1, 1); buffer = buffer + BF2*baseBudget0(:,:,:,3) ! BF1 is even
      
      case(6)
         ! Production: mean(delta u_i' delta u_j') partial_j mean(delta u_i)
         buffer = dudx * budget1(:,:,:,1) + dudy * budget1(:,:,:,2) + dudz * budget1(:,:,:,3) + &
                  dvdx * budget1(:,:,:,2) + dvdy * budget1(:,:,:,4) + dvdz * budget1(:,:,:,5) + &
                  dwdx * budget1(:,:,:,3) + dwdy * budget1(:,:,:,5) + dwdz * budget1(:,:,:,6)
      case(7)
         ! Production: mean(delta u_i' base u_j') partial_j mean(delta u_i)
         buffer = dudx * budget1(:,:,:,7) + dudy * budget1(:,:,:,8)  + dudz * budget1(:,:,:,10) + &
                  dvdx * budget1(:,:,:,9) + dvdy * budget1(:,:,:,12) + dvdz * budget1(:,:,:,13) + &
                  dwdx * budget1(:,:,:,11)+ dwdy * budget1(:,:,:,14) + dwdz * budget1(:,:,:,15)
      case(8)
         ! Production: mean(base u_i' delta u_j') partial_j mean(delta u_i)
         buffer = dudx * budget1(:,:,:,7) + dudy * budget1(:,:,:,9)  + dudz * budget1(:,:,:,11) + &
                  dvdx * budget1(:,:,:,8) + dvdy * budget1(:,:,:,12) + dvdz * budget1(:,:,:,14) + &
                  dwdx * budget1(:,:,:,10)+ dwdy * budget1(:,:,:,13) + dwdz * budget1(:,:,:,15)
      case(9)
         ! Production: mean(base u_i' base u_j') partial_j mean(delta u_i)
         buffer = dudx * baseBudget0(:,:,:,4) + dudy * baseBudget0(:,:,:,5) + dudz * baseBudget0(:,:,:,6) + &
                  dvdx * baseBudget0(:,:,:,5) + dvdy * baseBudget0(:,:,:,7) + dvdz * baseBudget0(:,:,:,8) + &
                  dwdx * baseBudget0(:,:,:,6) + dwdy * baseBudget0(:,:,:,8) + dwdz * baseBudget0(:,:,:,9)
      case(10)
         ! Production: mean(delta u_i' delta u_j') partial_j mean(base u_i)
         buffer = dudx_base * budget1(:,:,:,1) + dudy_base * budget1(:,:,:,2) + dudz_base * budget1(:,:,:,3) + &
                  dvdx_base * budget1(:,:,:,2) + dvdy_base * budget1(:,:,:,4) + dvdz_base * budget1(:,:,:,5) + &
                  dwdx_base * budget1(:,:,:,3) + dwdy_base * budget1(:,:,:,5) + dwdz_base * budget1(:,:,:,6)
      case(11)
         ! Production: mean(delta u_i' base u_j') partial_j mean(base u_i)
         buffer = dudx_base * budget1(:,:,:,7) + dudy_base * budget1(:,:,:,8)  + dudz_base * budget1(:,:,:,10) + &
                  dvdx_base * budget1(:,:,:,9) + dvdy_base * budget1(:,:,:,12) + dvdz_base * budget1(:,:,:,13) + &
                  dwdx_base * budget1(:,:,:,11)+ dwdy_base * budget1(:,:,:,14) + dwdz_base * budget1(:,:,:,15)
      case(12)
         ! Production: mean(base u_i' delta u_j') partial_j mean(base u_i)
         buffer = dudx_base * budget1(:,:,:,7) + dudy_base * budget1(:,:,:,9)  + dudz_base * budget1(:,:,:,11) + &
                  dvdx_base * budget1(:,:,:,8) + dvdy_base * budget1(:,:,:,12) + dvdz_base * budget1(:,:,:,14) + &
                  dwdx_base * budget1(:,:,:,10)+ dwdy_base * budget1(:,:,:,13) + dwdz_base * budget1(:,:,:,15)

      case(13)
         ! Buoyancy: mean(delta w' delta wb')
         buffer = - budget3(:,:,:,10)
      
      case(14)
         ! Buoyancy: mean(delta w' base wb')
         buffer = - budget3(:,:,:,11)

      case(15)
         ! Buoyancy covariance: mean(base w' delta wb')
         buffer = - budget3(:,:,:,12)

      case(16)
         ! Pressure covariance: mean(delta u_j' partial_j delta p')
         buffer = budget3(:,:,:,1)

      case(17)
         ! Pressure covariance: mean(base u_j' partial_j delta p')
         buffer = budget3(:,:,:,2)

      case(18)
         ! Pressure covariance: mean(delta u_j' partial_j base p')
         buffer = budget3(:,:,:,3)

      case(19)
         ! Transport: mean(delta u_i' delta u_j' partial_j delta u_i')
         buffer = budget3(:,:,:,19)

      case(20)
         ! Transport: mean(delta u_i' base u_j' partial_j delta u_i')
         buffer = budget3(:,:,:,18)

      case(21)
         ! Transport: mean(delta u_i' delta u_j' partial_j base u_i')
         buffer = budget3(:,:,:,17)

      case(22)
         ! Transport: mean(base u_i' delta u_j' partial_j delta u_i')
         buffer = budget3(:,:,:,16)

      case(23)
         ! Transport: mean(delta u_i' base u_j' partial_j base u_i')
         buffer = budget3(:,:,:,15)

      case(24)
         ! Transport: mean(base u_i' base u_j' partial_j delta u_i')
         buffer = budget3(:,:,:,14)

      case(25)
         ! Transport: mean(base u_i' delta u_j' partial_j base u_i')
         buffer = budget3(:,:,:,13)

      case(26)
         ! SGS transport: partial_j mean(base u_i' delta tau_ij')
         buffer = budget3(:,:,:,4)

      case(27)
         ! SGS transport: partial_j mean(delta u_i' base tau_ij')
         buffer = budget3(:,:,:,5)

      case(28)
         ! SGS transport: partial_j mean(delta u_i' delta tau_ij')
         buffer = budget3(:,:,:,6)

      case(29)
         ! SGS Dissipation: mean(delta tau_ij' partial_j base u_i')
         buffer = -budget3(:,:,:,7)

      case(30)
         ! SGS Dissipation: mean(base tau_ij' partial_j delta u_i')
         buffer = -budget3(:,:,:,8)

      case(31)
         ! SGS Dissipation: mean(delta tau_ij' partial_j delta u_i')
         buffer = -budget3(:,:,:,9)         
      end select 
  
      nullify(BF1, BF2)
   end subroutine

   subroutine resetEverything()
      implicit none

      if(allocated(budget0)) budget0 = zero
      if(allocated(budget1)) budget1 = zero
      if(allocated(budget2)) budget2 = zero
      if(allocated(budget3)) budget3 = zero
      if(allocated(baseBudget0)) baseBudget0 = zero
      if(allocated(duidxj)) duidxj = zero
      if(allocated(duidxj_base)) duidxj_base = zero
      if(allocated(profiles)) profiles = zero
   end subroutine

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
      integer, parameter :: HUGE_I = huge(1)

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

   subroutine compute_duidxj()
      implicit none
      call message(1, 'Computing velocity gradients ...')

      call ddx_R2R(budget0(:,:,:,1), dudx)
      call ddy_R2R(budget0(:,:,:,1), dudy)
      call ddz_R2R(budget0(:,:,:,1), dudz, uBC_bottom, uBC_top)
      call ddx_R2R(budget0(:,:,:,2), dvdx)
      call ddy_R2R(budget0(:,:,:,2), dvdy)
      call ddz_R2R(budget0(:,:,:,2), dvdz, vBC_bottom, vBC_top)
      call ddx_R2R(budget0(:,:,:,3), dwdx)
      call ddy_R2R(budget0(:,:,:,3), dwdy)
      call ddz_R2R(budget0(:,:,:,3), dwdz, wBC_bottom, wBC_top)

      call ddx_R2R(baseBudget0(:,:,:,1), dudx_base)
      call ddy_R2R(baseBudget0(:,:,:,1), dudy_base)
      call ddz_R2R(baseBudget0(:,:,:,1), dudz_base, uBC_bottom, uBC_top)
      call ddx_R2R(baseBudget0(:,:,:,2), dvdx_base)
      call ddy_R2R(baseBudget0(:,:,:,2), dvdy_base)
      call ddz_R2R(baseBudget0(:,:,:,2), dvdz_base, vBC_bottom, vBC_top)
      call ddx_R2R(baseBudget0(:,:,:,3), dwdx_base)
      call ddy_R2R(baseBudget0(:,:,:,3), dwdy_base)
      call ddz_R2R(baseBudget0(:,:,:,3), dwdz_base, wBC_bottom, wBC_top)
   end subroutine

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

   subroutine readBudgets(key, stamp)
   implicit none
   character(*), intent(in) :: key, stamp
   integer :: idx, budgetid
   character(len=clen) :: pattern, filename
   logical :: exists
   real(rkind), dimension(:,:,:,:), pointer :: budget
   
   do budgetid=0,3
      select case(budgetid)
      case(0)
         budget => budget0
      case(1)
         budget => budget1
      case(2)
         budget => budget2
      case(3)
         budget => budget3
      end select

      if((budgetid == 3) .and. (budgettype /= 4)) cycle ! budget3 is only relevant for TKE budgets

      do idx = 1, size(budget, 4)
         pattern  = getPattern(RID, budgetid, idx, key=key, stamp=stamp)
         filename = trim(inputdir)//'/'//trim(pattern)
         inquire(file=trim(filename), exist=exists)
         if(exists)then
            call message(1, 'Reading '//trim(filename))
            call decomp_2d_read_one(1, budget(:,:,:,idx), trim(filename), gpC)
         else
            call message(1, 'Not found: '//trim(filename)//' ... skipping')
            cycle
         end if
      end do
   end do

   call message(1, 'Reading base flow budget 0')
   do idx = 1,9
      pattern  = getPattern(BRID, 0, idx, key=key, stamp=stamp, isBase=.True.)
      filename = trim(inputdir)//'/'//trim(pattern)
      inquire(file=trim(filename), exist=exists)
      if(exists)then
         call message(1, 'Reading '//trim(filename))
         call decomp_2d_read_one(1, baseBudget0(:,:,:,idx), trim(filename), gpC)
      else
         call message(1, 'Not found: '//trim(filename)//' ... skipping')
         cycle
      end if
   end do

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
      allocate(duidxj(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 9))
      allocate(duidxj_base(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 9))
      allocate(Budget0(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 20))
      allocate(Budget1(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 15))
      allocate(Budget2(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 15))
      if(budgettype == 4) allocate( Budget3(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 19))
      allocate(baseBudget0(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 9))

      ! Allocate Buffers
      allocate(rbuffxC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3), 3))
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
      call message(0,'All memory allocated.')

      ! Initialize Padeder
      call Pade6opz%init(gpC, sp_gpC, gpE, sp_gpE, dz, NumericalSchemeVert,PeriodicInZ,spectC)
      call message(0,'Pade operations initialized')

      ! BCs for ddz
      call get_boundary_conditions_stencil()
      call message(0,'Identified boundary condition stenciles')

      ! Intersect the box with the mesh
      if(do_box_averaging)then
         call intersectBoxAndMesh()
         call message(0,'Control volume box intersected with the mesh')
      end if

      ! Allocate holder of x-profiles
      select case (budgettype)
      case(1)
         num_profiles = 24
      case(2)
         num_profiles = 24
      case(3)
         num_profiles = 24
      case(4)
         num_profiles = 31
      end select
      allocate(profiles(nx_box, num_profiles))

      ! Associate pointer
      dudx => duidxj(:,:,:,1)
      dudy => duidxj(:,:,:,2)
      dudz => duidxj(:,:,:,3)
      dvdx => duidxj(:,:,:,4)
      dvdy => duidxj(:,:,:,5)
      dvdz => duidxj(:,:,:,6)
      dwdx => duidxj(:,:,:,7)
      dwdy => duidxj(:,:,:,8)
      dwdz => duidxj(:,:,:,9)  

      dudx_base => duidxj_base(:,:,:,1)
      dudy_base => duidxj_base(:,:,:,2)
      dudz_base => duidxj_base(:,:,:,3)
      dvdx_base => duidxj_base(:,:,:,4)
      dvdy_base => duidxj_base(:,:,:,5)
      dvdz_base => duidxj_base(:,:,:,6)
      dwdx_base => duidxj_base(:,:,:,7)
      dwdy_base => duidxj_base(:,:,:,8)
      dwdz_base => duidxj_base(:,:,:,9)  

      call resetEverything()
   end subroutine

   subroutine release_memory()
    implicit none

      deallocate(mesh, duidxj, duidxj_base, Budget0, Budget1, Budget2, baseBudget0)
      if(allocated(Budget3)) deallocate(Budget3)
      if(allocated(rbuffxC)) deallocate(rbuffxC)
      if(allocated(cbuffyC)) deallocate(cbuffyC)
      if(allocated(cbuffzC)) deallocate(cbuffzC)
      if(allocated(profiles)) deallocate(profiles)
      if(allocated(xstations)) deallocate(xstations)
      
      nullify(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
      nullify(dudx_base, dudy_base, dudz_base, dvdx_base, dvdy_base, dvdz_base, dwdx_base, dwdy_base, dwdz_base)

      call spectC%destroy()
      call spectE%destroy()
      call Pade6opZ%destroy()
      call decomp_info_finalize(gpC)
      call decomp_info_finalize(gpE)
      call decomp_2d_finalize()
  end subroutine

end module constructDeficitBudgets_mod

program constructDeficitBudgets
   use constructDeficitBudgets_mod
   
   implicit none
   integer :: ioUnit, ierr, k
   logical :: periodicbcs(3)
   character(len=clen) :: inputfile, ers
      
   namelist /INPUT/ inputdir, outputdir, nx, ny, nz, Lx, Ly, Lz, prow, pcol, RID, &
                    BRID, budgettype, writeDependentVariables, startIDX, endIDX, tag, &
                    do_box_averaging
   namelist /NUMERICS/ NumericalSchemeVert
   namelist /BCs/ PeriodicInZ, botWall, topWall, botBC_temp
   namelist /BOX/ x1, x2, y1, y2, z1, z2

   ! Do MPI stuff
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)

   ! Do file IO - input file
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', status='old', action='read')
   read(unit=ioUnit, NML=INPUT, IOSTAT=ierr)
   if (ierr/=0)then
      write(ers,'(I)')ierr
      call gracefulExit("Reading failed for INPUT with error "//trim(ers), 101)
   end if
   read(unit=ioUnit, NML=NUMERICS, IOSTAT=ierr)
   if (ierr/=0)then
      write(ers,'(I)')ierr
      call gracefulExit("Reading failed for NUMERICS with error "//trim(ers), 102)
   end if
   read(unit=ioUnit, NML=BCs, IOSTAT=ierr)
   if (ierr/=0)then
      write(ers,'(I)')ierr
      call gracefulExit("Reading failed for BCs with error "//trim(ers), 103)
   end if
   read(unit=ioUnit, NML=BOX, IOSTAT=ierr)
   if (ierr/=0)then
      write(ers,'(I)')ierr
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
   
   ! Loop through time frames
   do k = 1, size(sorted_keys)
      call tic()

      if(.not. TimeWithinRange(trim(sorted_keys(k)), startIDX, endIDX)) cycle

      call message(0, 'Time Index: '//trim(sorted_keys(k))//', # Frames: '//trim(sorted_stamps(k)))

      ! Read Budgets
      call readBudgets(trim(sorted_keys(k)), trim(sorted_stamps(k)))

      ! Compute gradients
      call compute_duidxj()

      ! Compute Budgets
      call compute_budgets(trim(sorted_keys(k)), trim(sorted_stamps(k)))  

      ! Export profiles
      if((nrank == 0) .and. do_box_averaging)then
         call export_csv(trim(sorted_keys(k)), trim(sorted_stamps(k)))
      end if

      call resetEverything()
      call message(0, ' ')
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call toc()
   end do  

   call release_memory()  
   call MPI_FINALIZE(ierr) 

end program constructDeficitBudgets