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

        if (allocated(this%Fringe_kernel)) deallocate(this%Fringe_kernel)

    end subroutine destroy


    subroutine init(this, inputfile, nx, ny, nz, x, y, z, Lx, Ly, dz, dy)
        class(fringeAD), intent(inout) :: this
        character(*),    intent(in)    :: inputfile
        integer,         intent(in)    :: nx, ny, nz
        real(rkind),     intent(in)    :: x(nx), y(ny), z(nz)
        real(rkind),     intent(in)    :: dz, dy, Lx, Ly

        integer :: j, k, ioUnit, ierr

        real(rkind) :: xi_st(nx)
        real(rkind) :: xi_en(nx)
        real(rkind) :: S1(nx)
        real(rkind) :: S2(nx)
        real(rkind) :: Fringe_func(nx)
        real(rkind) :: xper(nx)

        real(rkind) :: fringe_len
        real(rkind) :: delta_st_
        real(rkind) :: delta_en_

        real(rkind) :: FringeAD_deltaH_
        real(rkind) :: FringeAD_deltaY_

        real(rkind) :: sigma_z
        real(rkind) :: sigma_y(ny)
        real(rkind) :: alpha

        real(rkind) :: y_fringe_start
        real(rkind) :: y_fringe_end
        real(rkind) :: y_tol
        real(rkind) :: dist_to_fringe

        logical :: in_y_fringe(ny)
        logical :: use_y_taper
        logical :: y_fringe_wrap
        logical :: y_fringe_active

        ! ------------------------------------------------------------------
        ! ADFRINGE namelist variables.
        !
        ! FringeAD_st and FringeAD_en define the x-extent over which the AD
        ! fringe is active. This interval is periodic in x.
        !
        ! Example:
        !
        !     FringeAD_st = 0.9
        !     FringeAD_en = 0.1
        !
        ! activates the AD fringe from 0.9*Lx to Lx and from 0 to 0.1*Lx.
        !
        ! FringeAD_H and FringeAD_deltaH control the vertical activation.
        ! Below FringeAD_H, the AD fringe is inactive. Above FringeAD_H, it is
        ! tapered with tanh^2 if use_tanh = .true.
        !
        ! FringeAD_taper_y controls whether the AD fringe is suppressed/tapered
        ! near the main y-fringe.
        ! ------------------------------------------------------------------
        real(rkind) :: FringeAD_st       = 0.875_rkind
        real(rkind) :: FringeAD_en       = one
        real(rkind) :: FringeAD_delta_st = 0.05_rkind
        real(rkind) :: FringeAD_delta_en = 0.075_rkind
        real(rkind) :: FringeAD_H        = 10._rkind
        real(rkind) :: FringeAD_deltaH   = 4._rkind
        real(rkind) :: FringeAD_deltaY   = 4._rkind
        logical     :: use_tanh          = .true.
        logical     :: FringeAD_taper_y  = .true.

        ! ------------------------------------------------------------------
        ! FRINGE namelist variables.
        !
        ! The AD-fringe y-taper uses:
        !
        !     Apply_y_fringe
        !     Fringe_yst
        !     Fringe_yen
        !
        ! The other variables are declared so that the full FRINGE namelist can
        ! be read safely.
        ! ------------------------------------------------------------------
        logical     :: Apply_x_fringe     = .true.
        logical     :: Apply_y_fringe     = .false.
        logical     :: do_shifts          = .false.

        real(rkind) :: Fringe_xst         = 0.75_rkind
        real(rkind) :: Fringe_xen         = one
        real(rkind) :: Fringe_yst         = one
        real(rkind) :: Fringe_yen         = one

        real(rkind) :: Fringe_delta_st_x  = one
        real(rkind) :: Fringe_delta_en_x  = one
        real(rkind) :: Fringe_delta_st_y  = one
        real(rkind) :: Fringe_delta_en_y  = one

        real(rkind) :: LambdaFact         = 2.45_rkind
        real(rkind) :: LambdaFact2        = 2.45_rkind
        real(rkind) :: LambdaFactPotTemp  = 2.45_rkind

        real(rkind) :: Fringe1_xst        = 0.75_rkind
        real(rkind) :: Fringe1_xen        = one
        real(rkind) :: Fringe1_delta_st_x = one
        real(rkind) :: Fringe1_delta_en_x = one

        real(rkind) :: Fringe2_xst        = 0.75_rkind
        real(rkind) :: Fringe2_xen        = one
        real(rkind) :: Fringe2_delta_st_x = one
        real(rkind) :: Fringe2_delta_en_x = one

        real(rkind) :: xshift             = zero
        real(rkind) :: yshift             = zero

        namelist /ADFRINGE/ FringeAD_st, FringeAD_en, FringeAD_delta_st, FringeAD_delta_en, &
                            FringeAD_H, FringeAD_deltaH, FringeAD_deltaY, use_tanh,         &
                            FringeAD_taper_y

        namelist /FRINGE/ Apply_x_fringe, Apply_y_fringe, Fringe_xst, Fringe_xen,           &
                          Fringe_delta_st_x, Fringe_delta_en_x, Fringe_delta_st_y,          &
                          Fringe_delta_en_y, LambdaFact, Fringe_yen, Fringe_yst,            &
                          LambdaFactPotTemp, LambdaFact2,                                   &
                          Fringe1_delta_st_x, Fringe1_delta_en_x, Fringe1_xst, Fringe1_xen, &
                          Fringe2_delta_st_x, Fringe2_delta_en_x, Fringe2_xst, Fringe2_xen, &
                          do_shifts, xshift, yshift

        if (allocated(this%Fringe_kernel)) deallocate(this%Fringe_kernel)
        allocate(this%Fringe_kernel(nx, ny, nz))

        ioUnit = 1019

        ! ------------------------------------------------------------------
        ! Read ADFRINGE namelist.
        ! ------------------------------------------------------------------
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        if (ierr /= 0) then
            call message('fringeADMethod:init could not open input file for ADFRINGE namelist.')
            this%Fringe_kernel = one
            return
        end if

        read(unit=ioUnit, nml=ADFRINGE, iostat=ierr)
        close(ioUnit)

        if (ierr /= 0) then
            call message('fringeADMethod:init could not read ADFRINGE namelist; using defaults where needed.')
        end if

        ! ------------------------------------------------------------------
        ! Read FRINGE namelist.
        !
        ! This is needed to identify the main y-fringe location.
        ! ------------------------------------------------------------------
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        if (ierr /= 0) then
            call message('fringeADMethod:init could not open input file for FRINGE namelist.')
            Apply_y_fringe = .false.
        else
            read(unit=ioUnit, nml=FRINGE, iostat=ierr)
            close(ioUnit)

            if (ierr /= 0) then
                call message('fringeADMethod:init could not read FRINGE namelist; assuming no y-fringe.')
                Apply_y_fringe = .false.
            end if
        end if

        call message(0, 'AD fringe parameters:')
        call message(1, 'FringeAD_st', FringeAD_st)
        call message(1, 'FringeAD_en', FringeAD_en)
        call message(1, 'FringeAD_delta_st', FringeAD_delta_st)
        call message(1, 'FringeAD_delta_en', FringeAD_delta_en)
        call message(1, 'FringeAD_H', FringeAD_H)
        call message(1, 'FringeAD_deltaH', FringeAD_deltaH)
        call message(1, 'FringeAD_deltaY', FringeAD_deltaY)
        if (use_tanh) then
            call message(1, 'use_tanh = .true.')
        else
            call message(1, 'use_tanh = .false.')
        end if
        if (FringeAD_taper_y) then
            call message(1, 'FringeAD_taper_y = .true.')
        else
            call message(1, 'FringeAD_taper_y = .false.')
        end if
        if (Apply_y_fringe) then
            call message(1, 'Apply_y_fringe = .true.')
        else
            call message(1, 'Apply_y_fringe = .false.')
        end if
        call message(1, 'Fringe_yst', Fringe_yst)
        call message(1, 'Fringe_yen', Fringe_yen)

        ! ------------------------------------------------------------------
        ! x-direction: build the periodic AD-fringe function.
        !
        ! Fringe_func = 1 outside the AD-fringe x-region.
        ! Fringe_func = 0 inside the AD-fringe x-region.
        !
        ! The interval is periodic:
        !
        !     xper       = modulo(x - FringeAD_st, Lx)
        !     fringe_len = modulo(FringeAD_en - FringeAD_st, Lx)
        !
        ! so FringeAD_st > FringeAD_en naturally wraps around the x-boundary.
        ! ------------------------------------------------------------------
        FringeAD_st = FringeAD_st * Lx
        FringeAD_en = FringeAD_en * Lx

        delta_st_ = max(FringeAD_delta_st, 0.005_rkind) * Lx
        delta_en_ = max(FringeAD_delta_en, 0.005_rkind) * Lx

        xper       = modulo(x - FringeAD_st, Lx)
        fringe_len = modulo(FringeAD_en - FringeAD_st, Lx)

        xi_st = xper / delta_st_
        xi_en = (xper - fringe_len) / delta_en_ + one

        call this%S_fringe(xi_st, S1)
        call this%S_fringe(xi_en, S2)

        Fringe_func = one - (S1 - S2)

        ! ------------------------------------------------------------------
        ! z-direction: vertical activation width.
        ! ------------------------------------------------------------------
        FringeAD_deltaH_ = max(two, FringeAD_deltaH) * abs(dz)

        ! ------------------------------------------------------------------
        ! y-direction: optional generic taper near the main y-fringe.
        !
        ! The main y-fringe is identified from Fringe_yst and Fringe_yen.
        !
        ! Case 1: high-y fringe
        !
        !     Fringe_yst = 0.9
        !     Fringe_yen = 1.0
        !
        !     y >= 0.9*Ly:
        !         AD fringe inactive
        !
        !     0.9*Ly - deltaY <= y <= 0.9*Ly:
        !         AD fringe tapered
        !
        !     y < 0.9*Ly - deltaY:
        !         AD fringe fully active, subject to x and z masks
        !
        ! Case 2: low-y fringe
        !
        !     Fringe_yst = 0.0
        !     Fringe_yen = 0.1
        !
        !     y <= 0.1*Ly:
        !         AD fringe inactive
        !
        !     0.1*Ly <= y <= 0.1*Ly + deltaY:
        !         AD fringe tapered
        !
        !     y > 0.1*Ly + deltaY:
        !         AD fringe fully active, subject to x and z masks
        !
        ! Case 3: wrapped y-fringe
        !
        !     Fringe_yst = 0.9
        !     Fringe_yen = 0.1
        !
        !     y >= 0.9*Ly or y <= 0.1*Ly:
        !         AD fringe inactive
        !
        !     The AD fringe is tapered away from both exposed edges.
        !
        ! If Apply_y_fringe = .false. or FringeAD_taper_y = .false.,
        ! sigma_y = 1 everywhere and the AD-fringe kernel has no y-dependence.
        ! ------------------------------------------------------------------
        FringeAD_deltaY_ = max(two, FringeAD_deltaY) * abs(dy)

        sigma_y(:)     = one
        in_y_fringe(:) = .false.

        use_y_taper = Apply_y_fringe .and. FringeAD_taper_y

        if (use_y_taper) then

            y_fringe_start = Fringe_yst * Ly
            y_fringe_end   = Fringe_yen * Ly
            y_tol          = epsilon(one) * max(one, Ly)

            y_fringe_active = abs(y_fringe_end - y_fringe_start) > y_tol

            if (y_fringe_active) then

                y_fringe_wrap = y_fringe_start > y_fringe_end

                do j = 1, ny

                    if (.not. y_fringe_wrap) then

                        ! Non-wrapped y-fringe:
                        !
                        !     y_fringe_start <= y <= y_fringe_end
                        !
                        ! Examples:
                        !     0.9*Ly -> Ly
                        !     0      -> 0.1*Ly
                        if (y(j) >= y_fringe_start .and. y(j) <= y_fringe_end) then

                            in_y_fringe(j) = .true.
                            sigma_y(j)     = zero

                        else

                            if (y(j) < y_fringe_start) then
                                dist_to_fringe = y_fringe_start - y(j)
                            else
                                dist_to_fringe = y(j) - y_fringe_end
                            end if

                            sigma_y(j) = tanh(dist_to_fringe / FringeAD_deltaY_)
                            sigma_y(j) = sigma_y(j) * sigma_y(j)

                        end if

                    else

                        ! Wrapped y-fringe:
                        !
                        !     y >= y_fringe_start or y <= y_fringe_end
                        !
                        ! Example:
                        !     0.9*Ly -> 0.1*Ly
                        if (y(j) >= y_fringe_start .or. y(j) <= y_fringe_end) then

                            in_y_fringe(j) = .true.
                            sigma_y(j)     = zero

                        else

                            dist_to_fringe = min(y(j) - y_fringe_end, y_fringe_start - y(j))

                            sigma_y(j) = tanh(dist_to_fringe / FringeAD_deltaY_)
                            sigma_y(j) = sigma_y(j) * sigma_y(j)

                        end if

                    end if

                end do

            end if

        end if

        ! ------------------------------------------------------------------
        ! Build Fringe_kernel(nx, ny, nz).
        !
        ! The smooth form is:
        !
        !     kernel = (1 - alpha) + alpha*Fringe_func
        !
        ! where:
        !
        !     alpha = sigma_z * sigma_y
        !
        ! Therefore:
        !
        !     alpha = 0:
        !         kernel = 1
        !         AD fringe inactive
        !
        !     alpha = 1:
        !         kernel = Fringe_func
        !         AD fringe fully active in x
        !
        !     0 < alpha < 1:
        !         smooth transition
        !
        ! The sharp form use_tanh = .false. uses the same in_y_fringe array
        ! as the smooth form, so both branches are logically consistent.
        ! ------------------------------------------------------------------
        do k = 1, nz

            if (use_tanh) then

                if (z(k) <= FringeAD_H) then
                    sigma_z = zero
                else
                    sigma_z = tanh((z(k) - FringeAD_H) / FringeAD_deltaH_)
                    sigma_z = sigma_z * sigma_z
                end if

                do j = 1, ny
                    alpha = sigma_z * sigma_y(j)
                    this%Fringe_kernel(:,j,k) = (one - alpha) + alpha * Fringe_func(:)
                end do

            else

                do j = 1, ny

                    if (z(k) < FringeAD_H) then

                        ! Below the AD-fringe height:
                        ! AD fringe inactive.
                        this%Fringe_kernel(:,j,k) = one

                    else if (use_y_taper .and. in_y_fringe(j)) then

                        ! Inside the main y-fringe:
                        ! AD fringe inactive.
                        this%Fringe_kernel(:,j,k) = one

                    else

                        ! Above FringeAD_H and outside the main y-fringe:
                        ! AD fringe fully active in x.
                        !
                        ! For use_tanh = .false., the transition in y is sharp.
                        this%Fringe_kernel(:,j,k) = Fringe_func(:)

                    end if

                end do

            end if

        end do

    end subroutine init


    subroutine S_fringe(this, x, output)
        class(fringeAD),           intent(inout) :: this
        real(rkind), dimension(:), intent(in)    :: x
        real(rkind), dimension(:), intent(out)   :: output

        integer     :: i
        real(rkind) :: exparg

        do i = 1, size(x)

            if (x(i) <= zero) then

                output(i) = zero

            else if (x(i) >= one) then

                output(i) = one

            else

                exparg = one / (x(i) - one + 1.0d-32) + one / (x(i) + 1.0d-32)
                exparg = min(exparg, 708.0d0)

                output(i) = one / (one + exp(exparg))

            end if

        end do

    end subroutine S_fringe

end module fringeADMethod
