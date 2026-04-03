module fringeADMethod
    use kind_parameters, only: rkind, clen
    use decomp_2d
    use exits, only: message
    use constants, only: zero, one, half, two
    implicit none
    private
    public :: fringeAD

    type :: fringeAD
    real(rkind), dimension(:,:,:), allocatable :: Fringe_kernel
    contains
        procedure :: init
        procedure :: destroy
        procedure :: S_fringe
    end type fringeAD

    contains

    subroutine destroy(this)
      class(fringeAD), intent(inout) :: this
      if(allocated(this%Fringe_kernel))deallocate(this%Fringe_kernel)
    end subroutine

    subroutine init(this, inputfile, nx, ny, nz, x, z, Lx, dz)
      class(fringeAD), intent(inout) :: this
      character(*), intent(in) :: inputfile 
      integer, intent(in) :: nx, ny, nz
      real(rkind), intent(in) :: x(nx), z(nz), dz, Lx
      integer :: k, ioUnit, ierr
      real(rkind) :: xi_st(nx), xi_en(nx), S1(nx), S2(nx), Fringe_func(nx)
      real(rkind) :: FringeAD_st = 0.875_rkind, FringeAD_en = one, FringeAD_delta_st=0.05_rkind, FringeAD_delta_en=0.075_rkind
      real(rkind) :: FringeAD_H = 10._rkind, FringeAD_deltaH=4._rkind
      logical :: use_tanh = .true.
      real(rkind) :: FringeAD_deltaH_
      real(rkind) :: xper(nx), fringe_len, delta_st_, delta_en_, sigma

      namelist /FRINGEAD/ FringeAD_st, FringeAD_en, FringeAD_delta_st, FringeAD_delta_en, FringeAD_H, FringeAD_deltaH, use_tanh

      if(allocated(this%Fringe_kernel))deallocate(this%Fringe_kernel)
      allocate(this%Fringe_kernel(nx, ny, nz))

      ioUnit = 1019
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=FRINGEAD)
      close(ioUnit)

      ! Scale up to Lx
      ! Note that FringeAD_H is already in proper units
      FringeAD_st = FringeAD_st * Lx
      FringeAD_en = FringeAD_en * Lx
      delta_st_ = max(FringeAD_delta_st, 0.005_rkind) * Lx
      delta_en_ = max(FringeAD_delta_en, 0.005_rkind) * Lx

      ! Periodic coordinate measured from fringe start, wrapped into [0,Lx)
      xper = modulo(x - FringeAD_st, Lx)

      ! Periodic forward length of the fringe region from start to end
      fringe_len = modulo(FringeAD_en - FringeAD_st, Lx)

      xi_st = xper / delta_st_
      xi_en = (xper - fringe_len) / delta_en_ + one
      
      ! FringeAD_deltaH is specified in units of vertical grid spacing.
      ! Enforce a minimum smooth transition width of 2*dz.
      FringeAD_deltaH_ = max(two, FringeAD_deltaH) * abs(dz)

      call this%S_fringe(xi_st, S1)
      call this%S_fringe(xi_en, S2)
      Fringe_func = one - (S1 - S2)

      do k = 1,nz
        if(use_tanh)then
          if (z(k) <= FringeAD_H) then
            sigma = zero
          else
            sigma = tanh((z(k) - FringeAD_H) / FringeAD_deltaH_)
            sigma = sigma*sigma ! tanh squared
          end if
          this%Fringe_kernel(:,:,k) = spread((one - sigma) + sigma*Fringe_func, dim=2, ncopies=ny)        
        else
          if(z(k) < FringeAD_H)then
            this%Fringe_kernel(:,:,k) = one
          else
            this%Fringe_kernel(:,:,k) = spread(Fringe_func, dim=2, ncopies=ny)
          end if
        end if
      end do
    end subroutine 

    subroutine S_fringe(this, x, output)
      class(fringeAD), intent(inout) :: this
      real(rkind), dimension(:), intent(in)    :: x
      real(rkind), dimension(:), intent(out)   :: output
      integer :: i
      real(rkind) :: exparg

      do i = 1,size(x)
        if (x(i) .le. zero) then
           output(i) = zero
        else if (x(i) .ge. one) then
           output(i) = one
        else
           exparg = one/(x(i) - one + 1.0D-32) + one/(x(i) + 1.0D-32)
           exparg = min(exparg,708.0d0) ! overflows if exparg > 709. need a better fix for this
           output(i) = one/(one + exp(exparg))
        end if
      end do
   end subroutine
end module fringeADMethod
