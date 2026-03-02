module derivatives_mod
    use mpi
    use exits,          only: message, gracefulExit
    use constants,      only: one, two, zero, half
    use kind_parameters,only: rkind, clen
    use PadeDerOps,     only: Pade6stagg
    use spectralMod,    only: spectral
    use decomp_2d
    use decomp_2d_io
    implicit none

    integer :: myrank, nprocs
    character(len=clen) :: inputdir, outputdir, filename
    character(len=1) :: derivative_type
    integer :: nx, ny, nz, prow=0, pcol=0
    real(rkind) :: Lx, Ly, Lz
    logical :: is_staggered = .false.
    integer :: bottom_BC=0, top_BC=0, NumericalSchemeVert=1

    type(decomp_info), target :: gpC, gpE
    type(spectral),    target :: spectC, spectE
    type(Pade6stagg)          :: Pade6opZ

    ! real buffers (X-pencil physical fields)
    real(rkind), allocatable, target :: rbuffxC(:,:,:,:), rbuffxE(:,:,:,:)

    ! complex buffers (spectral work arrays) - MUST match spect%spectdecomp
    complex(rkind), allocatable :: cbuffyC(:,:,:), cbuffyE(:,:,:), cbuffzE(:,:,:)
    complex(rkind), allocatable, target :: cbuffzC(:,:,:,:)

    abstract interface
        subroutine deriv_xy_iface(f, df)
            import rkind
            real(rkind), intent(in)  :: f(:,:,:)
            real(rkind), intent(out) :: df(:,:,:)
        end subroutine deriv_xy_iface
    end interface
    procedure(deriv_xy_iface), pointer :: ddx_ptr => null(), ddy_ptr => null()

contains

    subroutine assert_no_unit_thickness(gp, label)
        type(decomp_info), intent(in) :: gp
        character(len=*), intent(in)  :: label

        if ( gp%xsz(1)==1 .or. gp%xsz(2)==1 .or. gp%xsz(3)==1 .or. &
             gp%ysz(1)==1 .or. gp%ysz(2)==1 .or. gp%ysz(3)==1 .or. &
             gp%zsz(1)==1 .or. gp%zsz(2)==1 .or. gp%zsz(3)==1 ) then
            call message(0, "Warning: unit-thickness pencil detected in "//trim(label), 9100)
        end if
    end subroutine assert_no_unit_thickness

    !-----------------------------
    ! X-derivatives (spectral)
    !-----------------------------
    subroutine ddx_Cell(f, dfdx)
        real(rkind), intent(in)  :: f(:,:,:)
        real(rkind), intent(out) :: dfdx(:,:,:)

        call spectC%fft(f, cbuffyC)
        call spectC%mtimes_ik1_ip(cbuffyC)
        call spectC%dealias(cbuffyC)
        call spectC%ifft(cbuffyC, dfdx)
    end subroutine ddx_Cell

    subroutine ddx_Edge(f, dfdx)
        real(rkind), intent(in)  :: f(:,:,:)
        real(rkind), intent(out) :: dfdx(:,:,:)

        call spectE%fft(f, cbuffyE)
        call spectE%mtimes_ik1_ip(cbuffyE)
        call spectE%dealias(cbuffyE)
        call spectE%ifft(cbuffyE, dfdx)
    end subroutine ddx_Edge

    !-----------------------------
    ! Y-derivatives (spectral)
    !-----------------------------
    subroutine ddy_Cell(f, dfdy)
        real(rkind), intent(in)  :: f(:,:,:)
        real(rkind), intent(out) :: dfdy(:,:,:)

        call spectC%fft(f, cbuffyC)
        call spectC%mtimes_ik2_ip(cbuffyC)
        call spectC%dealias(cbuffyC)
        call spectC%ifft(cbuffyC, dfdy)
    end subroutine ddy_Cell

    subroutine ddy_Edge(f, dfdy)
        real(rkind), intent(in)  :: f(:,:,:)
        real(rkind), intent(out) :: dfdy(:,:,:)

        call spectE%fft(f, cbuffyE)
        call spectE%mtimes_ik2_ip(cbuffyE)
        call spectE%dealias(cbuffyE)
        call spectE%ifft(cbuffyE, dfdy)
    end subroutine ddy_Edge

    !-----------------------------
    ! Z-derivatives
    !-----------------------------
    subroutine ddz_Cell(f, dfdz, n1, n2)
        real(rkind), intent(in)  :: f(:,:,:)
        real(rkind), intent(out) :: dfdz(:,:,:)
        integer, intent(in) :: n1, n2

        call spectC%fft(f, cbuffyC)
        call transpose_y_to_z(cbuffyC, cbuffzC(:,:,:,1), spectC%spectdecomp)
        call Pade6opZ%ddz_C2C(cbuffzC(:,:,:,1), cbuffzC(:,:,:,2), n1, n2)
        call transpose_z_to_y(cbuffzC(:,:,:,2), cbuffyC, spectC%spectdecomp)
        call spectC%dealias(cbuffyC)
        call spectC%ifft(cbuffyC, dfdz)
    end subroutine ddz_Cell

    subroutine ddz_Edge(f, dfdz, n1, n2, n3, n4)
        real(rkind), intent(in)  :: f(:,:,:)
        real(rkind), intent(out) :: dfdz(:,:,:)
        integer, intent(in) :: n1, n2, n3, n4

        call spectE%fft(f, cbuffyE)
        call transpose_y_to_z(cbuffyE, cbuffzE, spectE%spectdecomp)
        call Pade6opZ%ddz_E2C(cbuffzE, cbuffzC(:,:,:,1), n1, n2)
        call Pade6opZ%interpz_C2E(cbuffzC(:,:,:,1), cbuffzE, n3, n4)
        call transpose_z_to_y(cbuffzE, cbuffyE, spectE%spectdecomp)
        call spectE%ifft(cbuffyE, dfdz)
    end subroutine ddz_Edge

end module derivatives_mod


program derivatives
    use derivatives_mod
    implicit none

    integer :: ierr, ioUnit
    real(rkind) :: dx, dy, dz
    character(len=clen) :: tmpname, outfile, inputfile
    character(len=3) :: tag
    logical :: exists
    real(rkind), pointer :: buffer(:,:,:), deriv(:,:,:)
    type(decomp_info), pointer :: gp => null()

    namelist /INPUT/ inputdir, outputdir, nx, ny, nz, Lx, Ly, Lz, prow, pcol, filename, derivative_type, &
                     is_staggered, bottom_BC, top_BC, NumericalSchemeVert

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call GETARG(1, inputfile)

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    if (mod(nx,2) /= 0) call gracefulExit("nx must be even.", 101)
    if (mod(ny,2) /= 0) call gracefulExit("ny must be even.", 102)

    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call decomp_info_init(nx, ny, nz,   gpC)
    call decomp_info_init(nx, ny, nz+1, gpE)

    call assert_no_unit_thickness(gpC, "gpC")
    call assert_no_unit_thickness(gpE, "gpE")

    call spectC%init("x", nx, ny, nz,   dx,dy,dz, "FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
    call spectE%init("x", nx, ny, nz+1, dx,dy,dz, "FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)

    call Pade6opZ%init(gpC, spectC%spectdecomp, gpE, spectE%spectdecomp, dz, NumericalSchemeVert, .false., spectC)

    ! Real buffers (physical, X-pencil)
    allocate(rbuffxC(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3), 2))
    allocate(rbuffxE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3), 2))

    ! Complex buffers MUST be sized from the spectral decomposition
    allocate(cbuffyC( spectC%spectdecomp%ysz(1), spectC%spectdecomp%ysz(2), spectC%spectdecomp%ysz(3) ))
    allocate(cbuffyE( spectE%spectdecomp%ysz(1), spectE%spectdecomp%ysz(2), spectE%spectdecomp%ysz(3) ))
    allocate(cbuffzC( spectC%spectdecomp%zsz(1), spectC%spectdecomp%zsz(2), spectC%spectdecomp%zsz(3), 2 ))
    allocate(cbuffzE( spectE%spectdecomp%zsz(1), spectE%spectdecomp%zsz(2), spectE%spectdecomp%zsz(3) ))

    ! Set pointers for which grid we’re operating on
    if (is_staggered) then
        buffer => rbuffxE(:,:,:,1)
        deriv  => rbuffxE(:,:,:,2)
        gp     => gpE
        ddx_ptr => ddx_Edge
        ddy_ptr => ddy_Edge
    else
        buffer => rbuffxC(:,:,:,1)
        deriv  => rbuffxC(:,:,:,2)
        gp     => gpC
        ddx_ptr => ddx_Cell
        ddy_ptr => ddy_Cell
    end if

    ! Read input
    tmpname = trim(inputdir)//"/"//trim(filename)
    inquire(file=trim(tmpname), exist=exists)
    if (.not. exists) then
        call message(1, 'Not found: '//trim(tmpname)//' ... exiting')
        call gracefulExit("Input file not found.", 2001)
    end if
    call message(1, 'Reading '//trim(tmpname))
    call decomp_2d_read_one(1, buffer, trim(tmpname), gp)

    ! Derivative selection
    select case (derivative_type)
    case ("x")
        call ddx_ptr(buffer, deriv)
        tag = "ddx"
    case ("y")
        call ddy_ptr(buffer, deriv)
        tag = "ddy"
    case ("z")
        if (is_staggered) then
            call ddz_Edge(buffer, deriv, bottom_BC, top_BC, 0, 0)
        else
            call ddz_Cell(buffer, deriv, bottom_BC, top_BC)
        end if
        tag = "ddz"
    case default
        call gracefulExit("Invalid derivative_type. Must be 'x', 'y', or 'z'.", 103)
    end select

    ! Write output
    outfile = trim(outputdir)//"/"//trim(tag)//"_"//trim(filename)
    call message(1, 'Writing '//trim(outfile))
    call decomp_2d_write_one(1, deriv, trim(outfile), gp)

    ! Cleanup
    deallocate(rbuffxC, rbuffxE, cbuffyC, cbuffyE, cbuffzC, cbuffzE)
    call spectC%destroy()
    call spectE%destroy()
    call Pade6opZ%destroy()
    call decomp_info_finalize(gpC)
    call decomp_info_finalize(gpE)
    call decomp_2d_finalize()
    call MPI_Finalize(ierr)

end program derivatives
