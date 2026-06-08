module fringeMethod
   use kind_parameters, only: rkind, clen
   use decomp_2d
   use spectralMod, only: spectral  
   use exits, only: message
   implicit none
   private
   public :: fringe

   type :: fringe 
      private
      logical, public                                       :: TargetsAssociated = .false. 
      real(rkind), public, dimension(:,:,:), pointer        :: u_target, v_target, w_target, T_target, F_target
      real(rkind), public, dimension(:,:,:), allocatable    :: u_for_shifts, v_for_shifts, w_for_shifts, T_for_shifts, F_for_shifts
      complex(rkind), dimension(:,:,:), pointer        :: uhat, vhat, what, That
      real(rkind), dimension(:,:,:), allocatable    :: Fringe_kernel_cells, Fringe_kernel_edges
      type(spectral),    pointer                    :: spectC, spectE
      type(decomp_info), pointer                    :: gpC, gpE, sp_gpC, sp_gpE
      real(rkind),    dimension(:,:,:,:), pointer   :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), pointer   :: cbuffyC, cbuffyE
      real(rkind)                                   :: LambdaFact, LambdaFactPotTemp, xshift=0, yshift=0
      integer :: myFringeID = 1
      logical :: useTwoFringex = .false., T_linked_for_shifts = .false.
      logical, public :: useFringeAsSponge_Scalar = .true., do_shifts = .false.
      logical :: firstCallComplete = .false.
      logical :: firstCallCompleteScalar = .false.
      contains
         procedure :: init
         procedure :: destroy
         procedure :: addFringeRHS
         procedure :: associateFringeTargets
         procedure :: allocateTargetArray_Cells
         procedure :: allocateTargetArray_Edges
         procedure :: addFringeRHS_scalar
         procedure :: associateFringeTarget_scalar
         procedure :: getFringeFraction
         procedure :: getLambdaFact
         procedure :: link_igrid_pointers
         procedure :: update_fringe_shifts
         procedure, private :: phaseshift_cell
         procedure, private :: phaseshift_edge
   end type
    
contains
   
   subroutine allocateTargetArray_Cells(this, array)
      class(fringe), intent(in)  :: this
      real(rkind), dimension(:,:,:), allocatable, intent(out) :: array

      allocate(array(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))

   end subroutine
   
   subroutine allocateTargetArray_Edges(this, array)
      class(fringe), intent(in)  :: this
      real(rkind), dimension(:,:,:), allocatable, intent(out) :: array

      allocate(array(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))

   end subroutine

   subroutine addFringeRHS(this, dt, urhs, vrhs, wrhs, uC, vC, wE, urhsF,vrhsF,wrhsF,addF)
      class(fringe),                                                                        intent(inout)  :: this
      real(rkind),                                                                         intent(in)     :: dt
      real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),          intent(in)     :: uC, vC 
      real(rkind),    dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)),          intent(in)     :: wE
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout)  :: urhs, vrhs
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout)  :: wrhs
      logical, intent(in), optional :: addF
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout), optional  :: urhsF, vrhsF
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout), optional  :: wrhsF
      
      logical :: AllOptionalsPresent

      if (present(urhsF) .and. present(vrhsF) .and. present(wrhsF) .and. present(addF)) then
         AllOptionalsPresent = .true. 
      else
         AllOptionalsPresent = .false. 
      end if


      if (this%targetsAssociated) then
         if (this%do_shifts) then
            call this%update_fringe_shifts()
         end if

         ! u velocity source term 
         this%rbuffxC(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_cells)*(this%u_target - uC)
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         urhs = urhs + this%cbuffyC(:,:,:,1)
         if (allOptionalsPresent) then
            if (addF) then
               urhsF = urhsF +  this%cbuffyC(:,:,:,1)
            else
               urhsF = this%cbuffyC(:,:,:,1)
            end if
         end if

         ! v velocity source term 
         this%rbuffxC(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_cells)*(this%v_target - vC)
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         vrhs = vrhs + this%cbuffyC(:,:,:,1)
         if (allOptionalsPresent) then
            if (addF) then
               vrhsF = vrhsF +  this%cbuffyC(:,:,:,1)
            else
               vrhsF = this%cbuffyC(:,:,:,1)
            end if
         end if
         
         ! w velocity source term 
         this%rbuffxE(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_edges)*(this%w_target - wE)
         call this%spectE%fft(this%rbuffxE(:,:,:,1), this%cbuffyE(:,:,:,1))      
         wrhs = wrhs + this%cbuffyE(:,:,:,1)
         if (allOptionalsPresent) then
            if (addF) then
               wrhsF = wrhsF +  this%cbuffyE(:,:,:,1)
            else
               wrhsF = this%cbuffyE(:,:,:,1)
            end if
         end if
      end if 

   end subroutine

   subroutine addFringeRHS_scalar(this, dt, Frhs, F)
      class(fringe),                                                                      intent(inout)        :: this
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)),intent(inout)        :: Frhs  
      real(rkind),                                                                        intent(in)           :: dt
      real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),         intent(in)           :: F 

      ! call message(0, 'DEBUG addFringeRHS_scalar: entering if statement')
      ! if (this%useFringeAsSponge_Scalar) then
      !    call message(1, 'DEBUG use fringe as sponge scalar')
      ! else
      !     call message(1, 'DEBUG no fringe as sponge scalar')
      ! end if
      ! if (associated(this%F_target)) then
      !     call message(1, 'DEBUG associated F target')
      ! else
      !     call message(1, 'DEBUG no associated F target')
      ! end if 
      if (this%firstCallCompleteScalar) then      
          if (this%useFringeAsSponge_Scalar) then
                this%rbuffxC(:,:,:,1) = -(this%LambdafactPotTemp/dt)*(this%Fringe_kernel_cells)*(F)
                call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
                Frhs = Frhs + this%cbuffyC(:,:,:,1)
                ! call message(1, 'DEBUG firstcallcompletescalar TRUE')
          else
             if (associated(this%F_target)) then
                this%rbuffxC(:,:,:,1) = (this%LambdafactPotTemp/dt)*(this%Fringe_kernel_cells)*(this%F_target - F)
                call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
                Frhs = Frhs + this%cbuffyC(:,:,:,1)
                ! call message(1, 'DEBUG ftarget TRUE')
             end if 
          end if 
      else
          this%firstCallCompleteScalar = .true.
      end if
   end subroutine

   subroutine destroy(this)
      class(fringe), intent(inout) :: this

      deallocate(this%Fringe_kernel_cells)
      deallocate(this%Fringe_kernel_edges)
      this%TargetsAssociated = .false.

      nullify(this%uhat, this%vhat, this%what, this%u_target, this%v_target, this%w_target)

      if (this%do_shifts) then
         deallocate(this%u_for_shifts, this%v_for_shifts, this%w_for_shifts)
         if (this%T_linked_for_shifts) then
            deallocate(this%T_for_shifts)
            nullify(this%That, this%T_target)
         end if
      end if
   end subroutine

   subroutine associateFringeTarget_scalar(this, Ftarget)
      class(fringe), intent(inout) :: this
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), target :: Ftarget

      this%F_target => Ftarget
      this%useFringeAsSponge_Scalar = .false. 

   end subroutine 

   subroutine associateFringeTargets(this, utarget, vtarget, wtarget, Ttarget)
      class(fringe), intent(inout) :: this
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), target           :: utarget, vtarget
      real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in), target           :: wtarget 
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), optional, target :: Ttarget 

      this%u_target => utarget
      this%v_target => vtarget
      this%w_target => wtarget

      if (present(Ttarget)) then
         this%T_target => Ttarget
         this%useFringeAsSponge_Scalar = .false. 
      end if
      this%TargetsAssociated = .true.

      call message(0, "Fringe targets successfully associated.")
   end subroutine

   subroutine init(this, inputfile, dx, x, dy, y, spectC, spectE, gpC, gpE, rbuffxC, rbuffxE, cbuffyC, cbuffyE, fringeID)
      use reductions, only: p_minval, p_maxval
      use exits, only: message_min_max
      use decomp_2d_io
      use mpi
      use constants, only : zero
      class(fringe), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile 
      type(decomp_info), intent(in), target :: gpC, gpE
      real(rkind), dimension(gpC%xsz(1)), intent(in) :: x
      real(rkind), dimension(gpC%xsz(2)), intent(in) :: y
      real(rkind), intent(in) :: dx, dy
      type(spectral), intent(in), target :: spectC, spectE
      real(rkind),    dimension(:,:,:,:), target, intent(in) :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), target, intent(in) :: cbuffyC, cbuffyE
      integer, intent(in), optional :: fringeID

      real(rkind) :: Lx, Ly, LambdaFact = 2.45d0, LambdaFact2 = 2.45d0, LambdaFactPotTemp = 2.45d0
      real(rkind) :: Fringe_yst = 1.d0, Fringe_yen = 1.d0
      real(rkind) :: Fringe_xst = 0.75d0, Fringe_xen = 1.d0
      real(rkind) :: Fringe_delta_st_x = 1.d0, Fringe_delta_st_y = 1.d0, Fringe_delta_en_x = 1.d0, Fringe_delta_en_y = 1.d0
      
      real(rkind) :: Fringe1_delta_st_x = 1.d0, Fringe1_delta_en_x = 1.d0
      real(rkind) :: Fringe2_delta_st_x = 1.d0, Fringe2_delta_en_x = 1.d0
      real(rkind) :: Fringe1_xst = 0.75d0, Fringe1_xen = 1.d0
      real(rkind) :: Fringe2_xst = 0.75d0, Fringe2_xen = 1.d0
      real(rkind) :: xshift = zero, yshift = zero

      real(rkind) :: small, big

      integer :: ioUnit = 10, i, j, k, nx, ierr
      real(rkind), dimension(:), allocatable :: x1, x2, Fringe_func, S1, S2, y1, y2
      logical :: Apply_x_fringe = .true., Apply_y_fringe = .false., do_shifts = .false.
      namelist /FRINGE/ Apply_x_fringe, Apply_y_fringe, Fringe_xst, Fringe_xen, Fringe_delta_st_x, Fringe_delta_en_x, &
                        Fringe_delta_st_y, Fringe_delta_en_y, LambdaFact, Fringe_yen, Fringe_yst, LambdaFactPotTemp, LambdaFact2, &
                        Fringe1_delta_st_x, Fringe1_delta_en_x, Fringe1_xst, Fringe1_xen, &  ! consider switching to using Fringe_xst, etc. instead of Fringe1_xst
                        Fringe2_delta_st_x, Fringe2_delta_en_x, Fringe2_xst, Fringe2_xen, &
                        do_shifts, xshift, yshift
    
      if (present(fringeID)) then
         this%myFringeID = fringeID
         this%useTwoFringex = .true. 
      end if
      nx = gpC%xsz(1)
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=FRINGE)
      close(ioUnit)

      Lx = maxval(x) + dx
      Ly = p_maxval(maxval(y)) + dy
      this%gpC => gpC
      this%gpE => gpE
      this%spectC => spectC
      this%spectE => spectE
      this%sp_gpC => spectC%spectdecomp
      this%sp_gpE => spectE%spectdecomp
      this%rbuffxC => rbuffxC 
      this%rbuffxE => rbuffxE
      this%cbuffyC => cbuffyC 
      this%cbuffyE => cbuffyE
      this%do_shifts = do_shifts  ! apply shifted boundary conditions - Munters, Meneveau, Meyers (2016)
      this%xshift = xshift        ! amount to shift x-boundary conditions (non-dim to L, not Lx)
      this%yshift = yshift        ! amount to shift y-boundary conditions (non-dim to L, not Ly)

      if (this%do_shifts) then
         call message(1, "Applying shifted boundary conditions to fringe targets")
         call message(2, "x-shift", this%xshift)
         call message(2, "y-shift", this%yshift)
      end if

      this%useFringeAsSponge_Scalar = .true. 
      
      allocate(this%Fringe_kernel_cells(nx, gpC%xsz(2), gpC%xsz(3)))
      allocate(this%Fringe_kernel_edges(nx, gpE%xsz(2), gpE%xsz(3)))
      
      this%Fringe_kernel_cells = 0.d0
      this%Fringe_kernel_edges = 0.d0
      ! Maybe need a flag to ensure that it is define and the problem is
      ! stratified
      this%LambdaFactPotTemp = LambdaFactPotTemp
      
      ! call message(0, 'DEBUG: lambdafactpottemp', this%LambdaFactPotTemp)  
      if (this%usetwoFringex) then
         select case (this%myFringeID)
         case(1)
            this%LambdaFact   = LambdaFact
            Fringe_xst        = Fringe1_xst*Lx
            Fringe_xen        = Fringe1_xen*Lx
            Fringe_delta_st_x = Fringe1_delta_st_x*Lx
            Fringe_delta_en_x = Fringe1_delta_en_x*Lx
         case(2)
            this%LambdaFact   = LambdaFact2
            Fringe_xst        = Fringe2_xst*Lx
            Fringe_xen        = Fringe2_xen*Lx
            Fringe_delta_st_x = Fringe2_delta_st_x*Lx
            Fringe_delta_en_x = Fringe2_delta_en_x*Lx
         end select
      else
            this%LambdaFact    = LambdaFact
            Fringe_xst        = Fringe_xst*Lx
            Fringe_xen        = Fringe_xen*Lx
            Fringe_delta_st_x = Fringe_delta_st_x*Lx
            Fringe_delta_en_x = Fringe_delta_en_x*Lx
      end if 
      
      if (Apply_x_fringe) then
         ! x - direction fringe
         allocate(x1         (nx))
         allocate(x2         (nx))
         allocate(S1         (nx))
         allocate(S2         (nx))
         allocate(Fringe_func(nx))
     
         x1 = ((x -  Fringe_xst)/Fringe_delta_st_x)
         x2 = ((x -  Fringe_xen)/Fringe_delta_en_x) + 1.d0
         call S_fringe(x1, S1)
         call S_fringe(x2, S2)
         Fringe_func = S1 - S2

         do k = 1,this%gpC%xsz(3)
            do j = 1,this%gpC%xsz(2)
                this%Fringe_kernel_cells(:,j,k) = Fringe_func    
            end do 
         end do

         do k = 1,this%gpE%xsz(3)
            do j = 1,this%gpE%xsz(2)
                this%Fringe_kernel_edges(:,j,k) = Fringe_func    
            end do 
         end do
         deallocate(x1, x2, S1, S2, Fringe_func)
      end if 
       
      small = p_minval(minval(this%Fringe_kernel_cells))
      big = p_maxval(maxval(this%Fringe_kernel_cells))
      call message_min_max(1,"Bounds for Fringe_funcC:", small, big)

      small = p_minval(minval(this%Fringe_kernel_edges))
      big = p_maxval(maxval(this%Fringe_kernel_edges))
      call message_min_max(1,"Bounds for Fringe_funcE:", small, big)

      if (Apply_y_fringe) then
         Fringe_yst        = Fringe_yst*Ly
         Fringe_yen        = Fringe_yen*Ly
         Fringe_delta_st_y = Fringe_delta_st_y*Ly
         Fringe_delta_en_y = Fringe_delta_en_y*Ly
         ! y direction fringe 1
         allocate(y1         (this%gpC%xsz(2)))
         allocate(y2         (this%gpC%xsz(2)))
         allocate(S1         (this%gpC%xsz(2)))
         allocate(S2         (this%gpC%xsz(2)))
         allocate(Fringe_func(this%gpC%xsz(2)))
     
         y1 = ((y -  Fringe_yst)/Fringe_delta_st_y)
         y2 = ((y -  Fringe_yen)/Fringe_delta_en_y) + 1.d0
         call S_fringe(y1, S1)
         call S_fringe(y2, S2)
         Fringe_func = S1 - S2

         do k = 1,this%gpC%xsz(3)
            do j = 1,this%gpC%xsz(2)
               do i = 1,nx
                  this%Fringe_kernel_cells(i,j,k) = this%Fringe_kernel_cells(i,j,k) + Fringe_func(j)
               end do 
            end do 
         end do

         do k = 1,this%gpE%xsz(3)
            do j = 1,this%gpE%xsz(2)
               do i = 1,nx
                  this%Fringe_kernel_edges(i,j,k) = this%Fringe_kernel_edges(i,j,k) + Fringe_func(j)
               end do 
            end do 
         end do
         deallocate(y1, y2, S1, S2, Fringe_func)
      
         small = p_minval(minval(this%Fringe_kernel_cells))
         big = p_maxval(maxval(this%Fringe_kernel_cells))
         call message_min_max(1,"Bounds for Fringe_funcC:", small, big)
      
         small = p_minval(minval(this%Fringe_kernel_edges))
         big = p_maxval(maxval(this%Fringe_kernel_edges))
         call message_min_max(1,"Bounds for Fringe_funcE:", small, big)

      end if 
      
      where (this%Fringe_kernel_cells > 1) this%Fringe_kernel_cells = 1.d0 
      where (this%Fringe_kernel_edges > 1) this%Fringe_kernel_edges = 1.d0 

      this%firstCallComplete = .false.
      this%firstCallCompleteScalar = .false.

      call message(0, "Fringe initialized successfully.")

   end subroutine


   subroutine link_igrid_pointers(this, uhat, vhat, what, That)
      class(fringe), intent(inout) :: this
      ! These are Fourier-space y-pencil arrays.  Their x extent is nx/2+1,
      ! so the physical gpC/gpE descriptors are not shape-compatible.
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in), target           :: uhat, vhat
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(in), target           :: what
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in), optional, target :: That

      this%uhat => uhat
      this%vhat => vhat
      this%what => what
      if (present(That)) then
         this%That => That
         this%T_linked_for_shifts = .true.
      end if

   end subroutine


   pure subroutine S_fringe(x, output)
      real(rkind), dimension(:), intent(in)    :: x
      real(rkind), dimension(:), intent(out)   :: output
      integer :: i
      real(rkind) :: exparg

      do i = 1,size(x)
        if (x(i) .le. 0.d0) then
           output(i) = 0.d0
        else if (x(i) .ge. 1.d0) then
           output(i) = 1.d0
        else
           exparg = 1.d0/(x(i) - 1.d0 + 1.0D-32) + 1.d0/(x(i) + 1.0D-32)
           exparg = min(exparg,708.0d0) ! overflows if exparg > 709. need a better fix for this
           output(i) = 1.d0/(1.d0 + exp(exparg))
        end if
      end do

   end subroutine

   subroutine getFringeFraction(this, output)
      use reductions, only: p_mean, p_sum
      class(fringe), intent(inout) :: this
      real(rkind), intent(out) :: output
      real(rkind) :: pcount

      ! computes the mean of the fringe kernel and returns
      ! note: p_mean doesn't always work if partitions are different sizes
      pcount = p_sum(this%gpC%xsz(1) * this%gpC%xsz(2) * this%gpC%xsz(3))
      output = p_sum(this%Fringe_kernel_cells) / pcount
   end subroutine

   subroutine getLambdaFact(this, output)
      class(fringe), intent(inout) :: this
      real(rkind), intent(out) :: output

      output = this%LambdaFact
   end subroutine


   subroutine update_fringe_shifts(this)
      class(fringe), intent(inout) :: this
      ! Perform lateral shifting here
      call this%phaseshift_cell(this%uhat, this%u_for_shifts, this%xshift, this%yshift)
      call this%phaseshift_cell(this%vhat, this%v_for_shifts, this%xshift, this%yshift)
      call this%phaseshift_edge(this%what, this%w_for_shifts, this%xshift, this%yshift)
      ! Temperature storage is absent in unstratified runs.  The association
      ! flag is a more reliable guard than passing another flow-state logical.
      if (this%T_linked_for_shifts) then
         call this%phaseshift_cell(this%That, this%T_for_shifts, this%xshift, this%yshift)
      end if
   end subroutine


   subroutine phaseshift_cell(this, uhat, uFilt, xshift, yshift)
      use constants, only: imi
      class(fringe), intent(inout) :: this
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(in) :: uhat
      complex(rkind), dimension(size(uhat,1),size(uhat,2),size(uhat,3)) :: tmp
      real(rkind),    dimension(this%spectC%physdecomp%xsz(1),this%spectC%physdecomp%xsz(2), this%spectC%physdecomp%xsz(3)), intent(out) :: uFilt
      real(rkind), intent(in) :: xshift, yshift
      integer :: i,j,k

      tmp = uhat
      do k = 1,size(tmp,3)             ! loop through all z-levels
         do j = 1,size(tmp,2)          ! loop through all y-wavenumbers
               do i = 1,size(tmp,1)    ! loop through all x-wavenumbers
                  tmp(i,j,k) = uhat(i,j,k)*exp(-imi*(this%spectC%k1(i,1,1)*xshift+this%spectC%k2(1,j,1)*yshift))
               end do
         end do
      end do
      call this%spectC%ifft(tmp, uFilt)  ! inverse FFT back to real space

   end subroutine

   subroutine phaseshift_edge(this, what, wFilt, xshift, yshift)
      use constants, only: imi
      class(fringe), intent(inout) :: this
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(in) :: what
      complex(rkind), dimension(size(what,1),size(what,2),size(what,3)) :: tmp
      real(rkind), dimension(this%spectE%physdecomp%xsz(1),this%spectE%physdecomp%xsz(2), &
                             this%spectE%physdecomp%xsz(3)), intent(out) :: wFilt
      real(rkind), intent(in) :: xshift, yshift
      integer :: i, j, k

      tmp = what
      do k = 1,size(tmp,3)
         do j = 1,size(tmp,2)
            do i = 1,size(tmp,1)
               tmp(i,j,k) = what(i,j,k)*exp(-imi*(this%spectE%k1(i,1,1)*xshift + &
                                                  this%spectE%k2(1,j,1)*yshift))
            end do
         end do
      end do
      call this%spectE%ifft(tmp, wFilt)
   end subroutine

end module 
