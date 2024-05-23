module TileFieldsMod
    use kind_parameters, only: rkind
    use exits, only: GracefulExit, message
    use constants, only : zero, one

    implicit none
    integer :: seed = 321341  ! seed for pseudo-random perturbations

contains

    ! tile arrays along axis 1
    subroutine tile_x(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, N, nx, nxf

        nx = size(arrIn, 1)
        nxf = size(arrOut, 1)

        if (.not. (modulo(nxf, nx) == 0)) then
            call GracefulExit('TILE_X: tiling must be an integer number', 999)
        endif

        N = nxf / nx  ! number of repeats

        ! tile data
        if (N == 1) then
            arrOut = arrIn
        else
            do j = 0, (N - 1) ! loop thru number of tiles
                do k = 1, nx
                    arrOut(k + nx * j, :, :) = arrIn(k, :, :)
                enddo
            enddo
        endif
    end subroutine

    ! tile arrays along axis 2
    subroutine tile_y(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, N, ny, nyf

        ny = size(arrIn, 2)
        nyf = size(arrOut, 2)

        if (.not. (modulo(nyf, ny) == 0)) then
            call GracefulExit('TILE_Y: tiling must be an integer number', 999)
        endif

        N = nyf / ny  ! number of repeats

        ! tile data
        if (N == 1) then
            arrOut = arrIn
        else
            do j = 0, (N - 1) ! loop thru number of tiles
                do k = 1, ny
                    arrOut(:, k + ny * j, :) = arrIn(:, k, :)
                enddo
            enddo
        endif
    end subroutine

    ! tile arrays along axis 3
    subroutine tile_z(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, N, nz, nzf
        logical :: is_edge

        nz = size(arrIn, 3)
        nzf = size(arrOut, 3)

        if (modulo(nzf, nz - 1) == 1) then
            is_edge = .true.  ! assume these are edge cells; handle differently
            N = nzf / (nz - 1)  ! number of repeats
        elseif (modulo(nzf, nz) == 0) then
            is_edge = .false.
            N = nzf / nz  ! number of repeats
        else
            call GracefulExit('TILE_Z: tiling must be an integer number', 999)
        endif

        ! tile data
        if (N == 1) then
            arrOut = arrIn
        elseif (.not. is_edge) then
            do j = 0, (N - 1) ! loop thru number of tiles
                do k = 1, nz
                    arrOut(:, :, k + nz * j) = arrIn(:, :, k)
                enddo
            enddo
        else
            do j = 0, (N - 1) ! loop thru number of tiles
                do k = 1, (nz - 1)
                    arrOut(:, :, k + (nz - 1) * j) = arrIn(:, :, k)
                enddo
            enddo
            arrOut(:, :, nzf) = arrIn(:, :, nz)  ! set equal top boundary
        endif
    end subroutine

    ! Add pseudo-random perturbations at the first grid level
    subroutine add_perturbations(arr, pfact, user_seed, nlevels)
        use random,             only: gaussian_random
        real(rkind), dimension(:,:,:), intent(inout) :: arr
        real(rkind), intent(in) :: pfact  ! magnitude of perturbations
        real(rkind), dimension(:,:,:), allocatable :: randArr
        integer, intent(in), optional :: nlevels
        integer, intent(in), optional :: user_seed
        integer :: N, k, local_seed

        if (present(nlevels)) then
            N = nlevels
        else
            N = 10  ! take some default value
        endif

        if (present(user_seed)) then  ! needed in parallel to break periodicity
            local_seed = user_seed + seed
        else
            local_seed = seed
        end if

        ! Add random numbers
        allocate(randArr(size(arr, 1), size(arr, 2), size(arr, 3)))
        call gaussian_random(randArr, zero, one, local_seed)
        do k = 1,min(size(arr,3), N)
            arr(:,:,k) = arr(:, :, k) + pfact*randArr(:,:,k)
        end do
        deallocate(randArr)
        call message(1, 'Added perturbations in the bottom N levels', min(size(arr,3), N))
    end subroutine
end module


program tileFields
    use kind_parameters, only: rkind, clen
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi
    use timer, only: tic, toc
    use decomp_2d_io
    use TileFieldsMod
    use basic_io

    implicit none

    character(len=clen) :: inputfile
    character(len=clen) :: outputdir, inputdir
    integer :: ioUnit, nx, ny, nz, ierr, inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID
    integer :: ntile_x=2, ntile_y=1, ntile_z=1, nxf, nyf, nzf
    integer :: i, io_status
    logical :: tileInZ = .false., isStratified = .false., periodicInZ = .false.
    type(decomp_info) :: gpC, gpE, gpC_upX, gpC_upXY, gpC_upXYZ, gpE_upX, gpE_upXY, gpE_upXYZ
    real(rkind), dimension(:,:,:), allocatable :: f, fxup_inX, fxup_inY, fxyup_inY
    real(rkind), dimension(:,:,:), allocatable :: fxyup_inZ, fxyup_inX, fxyzup_inZ, fxyzup_inY, fxyzup_inX
    real(rkind) :: tsim, frameangle=zero, pfact=one
    character(len=clen) :: tempname, fname
    character(len=clen), dimension(3) :: keys
    keys = [character(len=clen) :: "_u.", "_v.", "_T."]  !<-- cell-centered field names

    namelist /INPUT/ nx, ny, nz, ntile_x, ntile_y, ntile_z, &
        inputdir, outputdir, inputFile_TID, inputFile_RID, &
        outputFile_TID, outputFile_RID, isStratified, PeriodicInZ, pfact

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    ! Initialize decomp2d objects:
    call decomp_2d_init(nx, ny, nz, 0, 0)  !<-- original domain size, centered cells
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)  !<-- original domain size, edge cells
    nxf = nx * ntile_x
    nyf = ny * ntile_y
    nzf = nz * ntile_z
    if (ntile_z > 1) tileInZ = .true.  ! this will allocate additional arrays

    ! Print nx, ny, nz statements
    call message(0, 'Tiling up to nx', nxf)
    call message(0, 'Tiling up to ny', nyf)
    call message(0, 'Tiling up to nz', nzf)

    call decomp_info_init(nxf, ny , nz   , gpC_upX  )  !<-- up in X only
    call decomp_info_init(nxf, nyf, nz   , gpC_upXY )  !<-- up in X and Y only
    call decomp_info_init(nxf, nyf, nzf  , gpC_upXYZ)  !<-- up in X, Y, Z
    call decomp_info_init(nxf, ny , nz+1 , gpE_upX  )  !<-- up in X only (edges)
    call decomp_info_init(nxf, nyf, nz+1 , gpE_upXY )  !<-- up in X and Y only (edges)
    call decomp_info_init(nxf, nyf, nzf+1, gpE_upXYZ)  !<-- up in X, Y, Z (edges)

    !!!!!! READ HEADER !!!!!!!!!!!
    write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",inputFile_RID, "_info.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=10,file=fname,access='sequential',form='formatted')
    read (10, *)  tsim
    read (10, *, iostat=io_status)  frameangle
    close(10)
    call message(0, "Upsampling File data dumped at tSIM=", tsim)
    if (io_status == 0) then
        call message(1, 'Found frame angle to copy:', frameangle)
    end if

    !!!!!!!!!!!!! CELL FIELDS !!!!!!!!!!!!!!!
    allocate(f(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(fxup_inX(gpC_upX%xsz(1),gpC_upX%xsz(2),gpC_upX%xsz(3)))
    allocate(fxup_inY(gpC_upX%ysz(1),gpC_upX%ysz(2),gpC_upX%ysz(3)))  ! this is the same size as fxup_inX
    allocate(fxyup_inY(gpC_upXY%ysz(1),gpC_upXY%ysz(2),gpC_upXY%ysz(3)))
    allocate(fxyup_inX(gpC_upXY%xsz(1),gpC_upXY%xsz(2),gpC_upXY%xsz(3)))

    ! if tiling in Z, then make sure we have many additional arrays to allocate
    if (tileInZ) then
        allocate(fxyup_inZ(gpC_upXY%zsz(1),gpC_upXY%zsz(2),gpC_upXY%zsz(3)))
        allocate(fxyzup_inX(gpC_upXYZ%xsz(1),gpC_upXYZ%xsz(2),gpC_upXYZ%xsz(3)))
        allocate(fxyzup_inY(gpC_upXYZ%ysz(1),gpC_upXYZ%ysz(2),gpC_upXYZ%ysz(3)))
        allocate(fxyzup_inZ(gpC_upXYZ%zsz(1),gpC_upXYZ%zsz(2),gpC_upXYZ%zsz(3)))
    end if

    ! Read and Write cell fields
    do i = 1, size(keys)
        ! skip temperature if not stratified
        if ((trim(keys(i)) == "_T.") .and. (.not. isStratified)) cycle

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, trim(keys(i)),inputFile_TID
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,f,fname, gpC)  ! read original fields
        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, trim(keys(i)),outputFile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

        call tile_x(f,fxup_inX)
        call transpose_x_to_y(fxup_inX,fxup_inY,gpC_upX)
        call tile_y(fxup_inY,fxyup_inY)

        if (tileInZ) then
            if (periodicInZ) then
                ! only tile in z if periodic
                call transpose_y_to_z(fxyup_inY,fxyup_inZ,gpC_upXY)

                ! do tiling
                call tile_z(fxyup_inZ, fxyzup_inZ)

                ! add pseudo-random perturbations to break peridicity
                if ((trim(keys(i)) == "_T.") .and. (pfact > 0)) then
                    call message(0, 'Adding perturbations in the temperature field with magnitude', pfact)
                    call add_perturbations(fxyzup_inZ, pfact, nrank)
                endif

                call transpose_z_to_y(fxyzup_inZ,fxyzup_inY,gpC_upXYZ)
                call transpose_y_to_x(fxyzup_inY,fxyzup_inX,gpC_upXYZ)
                call decomp_2d_write_one(1,fxyzup_inX,fname, gpC_upXYZ)  ! write tiled fields
            else
                call GraceFulExit("ERROR: cannot tile in z if not periodic", 999)

            endif

        else
            ! still check if pfact > 0 (pseudo-random perturbations)
            if ((trim(keys(i)) == "_T.") .and. (pfact > 0)) then
                call message(0, 'Adding perturbations in the temperature field with magnitude', pfact)
                call add_perturbations(fxyup_inY, pfact, nrank)
            endif

            ! no tiling in z, write fields
            call transpose_y_to_x(fxyup_inY,fxyup_inX,gpC_upXY)
            call decomp_2d_write_one(1,fxyup_inX,fname, gpC_upXY)  ! write tiled fields
        end if

    enddo

    ! < SCALARS GO HERE > (or add to the list of keys: TODO)

    ! Deallocate variables to prep tiling edge variables
    deallocate(f, fxup_inX, fxup_inY, fxyup_inY)
    if (tileInZ) then
        deallocate(fxyup_inZ, fxyzup_inX, fxyzup_inY, fxyzup_inZ)
    else
        deallocate(fxyup_inX)
    end if

    !!!!!!!!!!!!! EDGE FIELDS !!!!!!!!!!!!!!!
    allocate(f(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(fxup_inX(gpE_upX%xsz(1),gpE_upX%xsz(2),gpE_upX%xsz(3)))
    allocate(fxup_inY(gpE_upX%ysz(1),gpE_upX%ysz(2),gpE_upX%ysz(3)))
    allocate(fxyup_inY(gpE_upXY%ysz(1),gpE_upXY%ysz(2),gpE_upXY%ysz(3)))

    if (tileInZ) then
        allocate(fxyup_inZ(gpE_upXY%zsz(1),gpE_upXY%zsz(2),gpE_upXY%zsz(3)))
        allocate(fxyzup_inX(gpE_upXYZ%xsz(1),gpE_upXYZ%xsz(2),gpE_upXYZ%xsz(3)))
        allocate(fxyzup_inY(gpE_upXYZ%ysz(1),gpE_upXYZ%ysz(2),gpE_upXYZ%ysz(3)))
        allocate(fxyzup_inZ(gpE_upXYZ%zsz(1),gpE_upXYZ%zsz(2),gpE_upXYZ%zsz(3)))
    else
        allocate(fxyup_inX(gpE_upXY%xsz(1),gpE_upXY%xsz(2),gpE_upXY%xsz(3)))
    end if

    ! Read and Write w - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_w.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpE)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_w.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

    call tile_x(f,fxup_inX)
    call transpose_x_to_y(fxup_inX,fxup_inY,gpE_upX)
    call tile_y(fxup_inY,fxyup_inY)

    if (tileInZ) then
        call transpose_y_to_z(fxyup_inY,fxyup_inZ,gpE_upXY)
        call tile_z(fxyup_inZ,fxyzup_inZ)
        call transpose_z_to_y(fxyzup_inZ,fxyzup_inY,gpE_upXYZ)
        call transpose_y_to_x(fxyzup_inY,fxyzup_inX,gpE_upXYZ)
        call decomp_2d_write_one(1,fxyzup_inX,fname, gpE_upXYZ)  ! write tiled fields
    else
        call transpose_y_to_x(fxyup_inY,fxyup_inX,gpE_upXY)
        call decomp_2d_write_one(1,fxyup_inX,fname, gpE_upXY)  ! write tiled fields
    end if

    deallocate(f, fxup_inX, fxup_inY, fxyup_inY)
    if (tileInZ) then
        deallocate(fxyup_inZ, fxyzup_inX, fxyzup_inY, fxyzup_inZ)
    else
        deallocate(fxyup_inX)
    end if

    !!!!!! WRITE HEADER !!!!!!!!!!!
    if (nrank == 0) then
        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",outputFile_RID, "_info.",outputFile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") tsim
        ! also write frame angle, if it was read in
        if (io_status == 0) then
            write(10, "(100g15.5)") frameangle
        end if
        close(10)
    end if

    call MPI_Finalize(ierr)           !<-- Terminate MPI

end program
