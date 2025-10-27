!===============================================================================
! Complete example: LES refinement with proper decomp2d integration
!
! This follows the tiling structure but uses spectral interpolation for X-Y
! refinement and physical interpolation for Z refinement.
!
! Structure:
!   1. Read field in X-pencils
!   2. Refine X-Y spectrally (combined operation in Fourier space)
!   3. If refining Z: transpose X→Y→Z, interpolate, transpose back Z→Y→X
!   4. Write refined field
!===============================================================================

module refine_fields_mod
    use mpi
    use exits, only: message, gracefulExit
    use kind_parameters,  only: rkind, clen
    use timer, only: tic, toc
    use PadeDerOps, only: Pade6stagg
    use spectralMod, only: spectral
    use decomp_2d 
    use decomp_2d_io
    use constants, only: zero
    implicit none

    type(Pade6stagg) :: Pade6opZ
  
    ! Decomposition info for intermediate grids
    ! Following the tiling pattern: original -> upX -> upXY -> upXYZ
    ! Cell-centered grids (for u, v, T)
    type(decomp_info) :: gpC        ! Original coarse grid (cell centers)
    type(decomp_info) :: gpC_XY     ! Refined in X and Y (cell centers)
    type(decomp_info) :: gpC_XYZ    ! Refined in X, Y, and Z (cell centers)

    ! Edge-based grids (for w velocity, staggered in z)
    type(decomp_info) :: gpE        ! Original coarse grid (edges, nz+1)
    type(decomp_info) :: gpE_XY     ! Refined in X and Y (edges, nz+1)
    type(decomp_info) :: gpE_XYZ    ! Refined in X, Y, and Z (edges, nz_f+1)

    type(decomp_info), pointer :: Sp_gpC_c, Sp_gpC_XY, Sp_gpE_c, Sp_gpE_XY
    type(spectral), target  :: spectE_c, spectC_c, spectE_f, spectC_f, spectC_XY, spectE_XY
    
    ! Intermediate arrays
    real(rkind), allocatable :: fxy_inX(:,:,:)   ! Refined X-Y, in X-pencils
    real(rkind), allocatable :: fxy_inY(:,:,:)   ! Refined X-Y, in Y-pencils
    real(rkind), allocatable :: fxy_inZ(:,:,:)   ! Refined X-Y, in Z-pencils
    real(rkind), allocatable :: fxyz_inY(:,:,:)   ! Refined X-Y-Z, in Y-pencils
    real(rkind), allocatable :: fxyz_inZ(:,:,:)   ! Refined X-Y-Z, in Z-pencils
    real(rkind), allocatable :: fxyE_inX(:,:,:)  ! Refined X-Y, in X-pencils
    real(rkind), allocatable :: fxyE_inY(:,:,:)  ! Refined X-Y, in Y-pencils
    real(rkind), allocatable :: fxyE_inZ(:,:,:)  ! Refined X-Y, in Z-pencils
    real(rkind), allocatable :: fxyzE_inY(:,:,:)   ! Refined X-Y-Z, in Y-pencils
    real(rkind), allocatable :: fxyzE_inZ(:,:,:)   ! Refined X-Y-Z, in Z-pencils

    complex(rkind), allocatable :: cbuffyC(:,:,:), cbuffzC1(:,:,:), cbuffzC2(:,:,:)
    complex(rkind), allocatable :: cbuffyE(:,:,:), cbuffzE1(:,:,:)
    
    ! Parity flags for ddz
    integer :: uBC_bottom, uBC_top, vBC_bottom, vBC_top, wBC_bottom, wBC_top, TBC_bottom, TBC_top, dWdzBC_bottom, dWdzBC_top

    integer :: refine_x=2, refine_y=2, refine_z=1
    logical :: isStratified=.true.

    real(rkind), allocatable :: u_c(:,:,:), v_c(:,:,:), w_c(:,:,:), T_c(:,:,:)
    real(rkind), allocatable :: u_f(:,:,:), v_f(:,:,:), w_f(:,:,:), T_f(:,:,:)
            
    contains

    subroutine write_restart_file(field, outputdir, outputFile_TID, outputFile_RID, name, gp)
      implicit none
      real(rkind), dimension(:,:,:), intent(in) :: field
      character(len=*), intent(in) :: outputdir, name
      integer, intent(in) :: outputFile_TID, outputFile_RID
      type(decomp_info), intent(in) :: gp
      character(len=clen) :: tempname, fname

      write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID,trim(name),outputFile_TID
      fname = trim(outputdir)//"/"//trim(tempname)
      call decomp_2d_write_one(1,field,fname, gp)  ! write refined fields
    end subroutine write_restart_file

    subroutine read_restart_file(field, inputdir, inputFile_TID, inputFile_RID, name, gp)
      implicit none
      real(rkind), dimension(:,:,:), intent(out) :: field
      character(len=*), intent(in) :: inputdir, name
      integer, intent(in) :: inputFile_TID, inputFile_RID
      type(decomp_info), intent(in) :: gp
      character(len=clen) :: tempname, fname

      write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID,trim(name),inputFile_TID
      fname = trim(inputdir)//"/"//trim(tempname)
      call decomp_2d_read_one(1,field,fname, gp)  ! read original fields
    end subroutine read_restart_file

    subroutine get_boundary_conditions_stencil(botWall, TopWall, botBC_Temp, topBC_Temp)
      implicit none
      integer, intent(in) :: botWall, TopWall, topBC_Temp, botBC_Temp
      
      wBC_bottom     = -1; wBC_top     = -1
      uBC_bottom     =  0; uBC_top     =  1
      vBC_bottom     =  0; vBC_top     =  1
      TBC_bottom     =  1; TBC_top     =  0
      dWdzBC_bottom  =  0; dWdzBC_top  =  0

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
         dwdzBC_bottom = -1
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
         dwdzBC_top = -1
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

      !! Temperature
      select case (topBC_Temp)
      case(0) ! Dirichlet (default)
          TBC_top = 0
      case(1)
          TBC_top = 1
        case (2) ! Inhomogeneous Neumann BC for temperature at the top
          TBC_top = 0
      case (3)
          TBC_top = 0
      end select 
      select case (botBC_Temp)
      case (0) ! Dirichlet BC for temperature at the bottom
          TBC_bottom = 0      
      case(1)  ! Homogenenous Neumann BC at the bottom
          TBC_bottom = 1
        case (2) ! Inhomogeneous Neumann BC for temperature at the bottom
          TBC_bottom = 0
      case (3) 
          TBC_bottom = 0    
      end select
      
   end subroutine get_boundary_conditions_stencil

  !-----------------------------------------------------------------------------
  ! Refine a single field
  ! 1. Refine X-Y spectrally (combined operation)
  ! 2. If refining Z: transpose, interpolate, transpose back
  ! Optional: handle staggered grids in z-direction
  !-----------------------------------------------------------------------------
  subroutine refine_single_field(field_c, field_f, dz, n1, n2)
    implicit none
    
    real(rkind), dimension(:,:,:), intent(in) :: field_c
    real(rkind), dimension(:,:,:), intent(out) :: field_f
    real(rkind), intent(in) :: dz
    integer, intent(in) :: n1, n2

    ! Step 1: Horizontal refinement (X-Y) using spectral interpolation
    call refine_horizontally(field_c, fxy_inX, spectC_c, spectC_f)

    ! Step 2: Handle Z-refinement if needed
    if (refine_z > 1) then
      
      call transpose_x_to_y(fxy_inX, fxy_inY, gpC_XY)
      call transpose_y_to_z(fxy_inY, fxy_inZ, gpC_XY)
      call refine_z_physical(fxy_inZ, fxyz_inZ, dz, .false., n1, n2)
      call transpose_z_to_y(fxyz_inZ, fxyz_inY, gpC_XYZ)
      call transpose_y_to_x(fxyz_inY, field_f, gpC_XYZ)   
      
    else
      field_f(:,:,:) = fxy_inX(:,:,:)
    end if
    
  end subroutine refine_single_field

  subroutine refine_single_fieldE(field_c, field_f, dz, n1, n2, n3, n4)
    implicit none
    
    real(rkind), dimension(:,:,:), intent(in) :: field_c
    real(rkind), dimension(:,:,:), intent(out) :: field_f
    real(rkind), intent(in) :: dz
    integer, intent(in) :: n1, n2, n3, n4

    ! Step 1: Horizontal refinement (X-Y) using spectral interpolation
    call refine_horizontally(field_c, fxyE_inX, spectE_c, spectE_f)

    ! Step 2: Handle Z-refinement if needed
    if (refine_z > 1) then      
      call transpose_x_to_y(fxyE_inX, fxyE_inY, gpE_XY)
      call transpose_y_to_z(fxyE_inY, fxyE_inZ, gpE_XY)
      call refine_z_physical(fxyE_inZ, fxyzE_inZ, dz, .true., n1, n2, n3=n3, n4=n4)
      call transpose_z_to_y(fxyzE_inZ, fxyzE_inY, gpE_XYZ)
      call transpose_y_to_x(fxyzE_inY, field_f, gpE_XYZ)       
    else
      field_f(:,:,:) = fxyE_inX(:,:,:)
    end if
    
  end subroutine refine_single_fieldE

    subroutine refine_horizontally(field_c, field_f, spect_c, spect_f)
    implicit none

    ! Arguments
    real(rkind), intent(in)  :: field_c(:,:,:)  ! Coarse Physical (X-pencil)
    real(rkind), intent(out) :: field_f(:,:,:)  ! Fine Physical (X-pencil)
    type(spectral), intent(inout) :: spect_c    
    type(spectral), intent(inout) :: spect_f    

    ! Internal Complex Buffers
    complex(rkind), allocatable :: hat_c_yp(:,:,:) ! Coarse Y-pencil
    complex(rkind), allocatable :: hat_i_yp(:,:,:) ! Intermediate Y-pencil (Fine Y, Coarse X)
    complex(rkind), allocatable :: hat_i_xp(:,:,:) ! Intermediate X-pencil (Fine Y, Coarse X)
    complex(rkind), allocatable :: hat_f_xp(:,:,:) ! Fine X-pencil (Fine Y, Fine X)
    complex(rkind), allocatable :: hat_f_yp(:,:,:) ! Fine Y-pencil (Fine Y, Fine X)

    type(decomp_info) :: decomp_inter  
    integer :: nxc_g, nyc_g, nxf_g, nyf_g, nzc_g
    integer :: nxc_hat
    !real(rkind) :: scale

    nxc_g = spect_c%nx_g ; nyc_g = spect_c%ny_g ; nzc_g = spect_c%nz_g
    nxf_g = spect_f%nx_g ; nyf_g = spect_f%ny_g
    nxc_hat = nxc_g/2 + 1

    ! Initialize Intermediate Decomposition (Coarse X_hat, Fine Y, Coarse Z)
    call decomp_info_init(nxc_hat, nyf_g, nzc_g, decomp_inter)

    !===============================================================
    ! SAFEGUARDS (single place, integer-only checks)
    !===============================================================
    ! Rationale:
    !   - This routine assumes a “slab-like” process grid where:
    !       * y-pencils cover the full global y-range locally  (ysz(2) == Ny)
    !       * x-pencils cover the full global kx-range locally (xsz(1) == Nx_hat)
    !   - It also assumes we can local-copy hat_i_xp -> hat_f_xp, which requires
    !     the y/z partitioning (in x-pencil layout) to match between decomp_inter
    !     and spect_f%spectdecomp.
    !
    ! If any of these assumptions are violated,
    ! this routine must fall back to global-index mapping instead of local slices,
    ! which is not currently implemented.
    !===============================================================

    !-----------------------------
    ! (A) Full-y slabs in y-pencils
    !-----------------------------
    if (spect_c%spectdecomp%ysz(2) /= spect_c%ny_g)&
        call gracefulExit("spect_c does not have full-y in y-pencils.", 001)
    
    if (spect_f%spectdecomp%ysz(2) /= spect_f%ny_g)&
        call gracefulExit("spect_f does not have full-y in y-pencils.", 002)

    !--------------------------------------------
    ! (B) Full-kx slabs in x-pencils for padding
    !--------------------------------------------
    ! Intermediate x-pencil should contain full kx = 1..nxc_hat locally
    ! where nxc_hat = nx_c/2 + 1
    if (decomp_inter%xsz(1) /= (spect_c%nx_g/2 + 1))&
        call gracefulExit("Intermediate decomposition does not have full-kx in x-pencils.", 003)

    ! Fine x-pencil should contain full kx = 1..nx_f/2+1 locally
    ! (not strictly required for the copy of only 1:nxc_hat, but it is the
    !  assumption behind using local indices without global mapping)
    if (spect_f%spectdecomp%xsz(1) /= (spect_f%nx_g/2 + 1))&
        call gracefulExit("Fine decomposition does not have full-kx in x-pencils.", 004)

    !---------------------------------------------------------
    ! (C) X-pencil y/z partition must match for local copy
    !     hat_f_xp(1:nxc_hat,:,:) = hat_i_xp(1:nxc_hat,:,:)
    !---------------------------------------------------------
    if (decomp_inter%xsz(2) /= spect_f%spectdecomp%xsz(2))&
        call gracefulExit("Intermediate decomposition y-size in x-pencils does not match fine decomposition.", 005)
    
    if (decomp_inter%xsz(3) /= spect_f%spectdecomp%xsz(3))&
        call gracefulExit("Intermediate decomposition z-size in x-pencils does not match fine decomposition.", 006)

    if (decomp_inter%xst(2) /= spect_f%spectdecomp%xst(2))&
        call gracefulExit("Intermediate decomposition y-start in x-pencils does not match fine decomposition.", 007)
    
    if (decomp_inter%xst(3) /= spect_f%spectdecomp%xst(3))&
        call gracefulExit("Intermediate decomposition z-start in x-pencils does not match fine decomposition.", 008)
    
    if (decomp_inter%xen(2) /= spect_f%spectdecomp%xen(2))&
        call gracefulExit("Intermediate decomposition y-end in x-pencils does not match fine decomposition.", 009)
    
    if (decomp_inter%xen(3) /= spect_f%spectdecomp%xen(3))&
        call gracefulExit("Intermediate decomposition z-end in x-pencils does not match fine decomposition.", 010)

    ! 1. Forward FFT -> Result in Y-pencil
    allocate(hat_c_yp(spect_c%spectdecomp%ysz(1), spect_c%spectdecomp%ysz(2), spect_c%spectdecomp%ysz(3)))
    if (size(hat_c_yp,2) /= nyc_g) call GracefulExit("hat_c_yp does not contain full y locally", 011)
    call spect_c%fft(field_c, hat_c_yp)

    ! 2. Pad Y-direction locally (Intermediate Y-pencil)
    allocate(hat_i_yp(decomp_inter%ysz(1), decomp_inter%ysz(2), decomp_inter%ysz(3)))
    if (size(hat_i_yp,2) /= nyf_g) call GracefulExit("hat_i_yp does not contain full y locally", 012)
    hat_i_yp = (zero, zero)

    ! Non-negative block includes Nyquist
    hat_i_yp(:, 1:nyc_g/2+1, :) = hat_c_yp(:, 1:nyc_g/2+1, :)

    ! Strictly negative modes only (length = nyc/2 - 1)
    hat_i_yp(:, nyf_g-(nyc_g/2-1)+1:nyf_g, :) = hat_c_yp(:, nyc_g/2+2:nyc_g, :)

    ! 3. Transpose to X-pencil to handle X-padding locally
    allocate(hat_i_xp(decomp_inter%xsz(1), decomp_inter%xsz(2), decomp_inter%xsz(3)))
    call transpose_y_to_x(hat_i_yp, hat_i_xp, decomp_inter)

    ! 4. Pad X-direction locally (Fine X-pencil)
    allocate(hat_f_xp(spect_f%spectdecomp%xsz(1), spect_f%spectdecomp%xsz(2), spect_f%spectdecomp%xsz(3)))
    hat_f_xp = (zero, zero)
    hat_f_xp(1:nxc_hat, :, :) = hat_i_xp(1:nxc_hat, :, :)

    ! 5. Scaling
    ! scale = (real(nxf_g, rkind)/real(nxc_g, rkind)) * (real(nyf_g, rkind)/real(nyc_g, rkind))
    ! hat_f_xp = hat_f_xp * scale

    ! 6. Transpose Fine X-pencil back to Fine Y-pencil for the IFFT
    allocate(hat_f_yp(spect_f%spectdecomp%ysz(1), spect_f%spectdecomp%ysz(2), spect_f%spectdecomp%ysz(3)))
    call transpose_x_to_y(hat_f_xp, hat_f_yp, spect_f%spectdecomp)

    ! 7. Inverse FFT (Fine Y-pencil to Fine X-physical)
    call spect_f%ifft(hat_f_yp, field_f)

    ! Cleanup
    deallocate(hat_c_yp, hat_i_yp, hat_i_xp, hat_f_xp, hat_f_yp)
    call decomp_info_finalize(decomp_inter)

  end subroutine refine_horizontally

  subroutine refine_z_physical(field_c, field_f, dz_c, staggered, bottom_flag, top_flag, n3, n4)
    !---------------------------------------------------------------------------
    ! Vertical refinement in physical z.
    !
    !  - If staggered=.true.  (edge nodes): nz_f = nz_c*refine_z, nodes are nz+1.
    !    -> pure interval refinement (no extrapolation needed); Hermite everywhere
    !       with exact top node copy.
    !
    !  - If staggered=.false. (cell centres): nz_f = nz_c*refine_z, nodes are nz.
    !    -> fine centres extend +/- (dz_c - dz_f)/2 beyond coarse-centre set.
    !       We therefore:
    !         * use one-sided Taylor extrapolation for boundary fine centres that
    !           lie outside [z_c(1), z_c(nz_c)]
    !         * use cubic Hermite in the interior
    !
    ! Uses coarse physical gradient computed by:
    !   call ddz_R2R(f, dfdz, bottom_flag, top_flag)
    ! where dfdz is in physical units (per metre).
    !---------------------------------------------------------------------------

    implicit none
    logical, intent(in) :: staggered
    real(rkind), intent(in) :: dz_c
    integer, intent(in) :: bottom_flag, top_flag
    integer, intent(in), optional :: n3, n4  ! optional arguments for staggered fields (w)

    real(rkind), intent(in)  :: field_c(:,:,:)   ! coarse (centres or edges)
    real(rkind), intent(out) :: field_f(:,:,:)   ! fine   (centres or edges)

    integer :: nx, ny
    integer :: nz3_c, nz3_f
    integer :: nz_c, nz_f
    integer :: nz_nodes_c, nz_nodes_f
    integer :: i, j, kf, kc
    integer :: q, s
    real(rkind) :: t
    real(rkind) :: f0, f1, m0, m1
    real(rkind) :: h00, h10, h01, h11
    real(rkind) :: dz_f
    real(rkind) :: zc1, zcN, zf
    real(rkind) :: z0, fbc, mbc
    real(rkind), allocatable :: dfdz_c(:,:,:)
    integer :: n3_, n4_

    nx    = size(field_c,1)
    ny    = size(field_c,2)
    nz3_c = size(field_c,3)
    nz3_f = size(field_f,3)

    if (refine_z < 1) call GracefulExit("refine_z_physical: refine_z must be >= 1", 801)
    dz_f = dz_c / real(refine_z, rkind)

    if (staggered) then
        ! edges: field has nz+1 nodes
        nz_c = nz3_c - 1
        if (nz_c < 1) call GracefulExit("refine_z_physical: staggered needs >=2 edge nodes", 802)

        nz_f = nz_c * refine_z
        if (nz3_f /= nz_f + 1) call GracefulExit("refine_z_physical: fine staggered must be nz_f+1", 803)

    else
        ! centres: field has nz nodes
        nz_c = nz3_c
        if (nz_c < 2) call GracefulExit("refine_z_physical: centred needs >=2 nodes", 804)

        nz_f = nz_c * refine_z
        if (nz3_f /= nz_f) call GracefulExit("refine_z_physical: fine centred must be nz_f", 805)
    end if

    nz_nodes_c = nz3_c
    nz_nodes_f = nz3_f

    ! Coarse physical gradient at the same nodes as field_c
    allocate(dfdz_c(nx, ny, nz_nodes_c))
    if(staggered)then
      if(present(n3) .and. present(n4)) then
        n3_ = n3
        n4_ = n4
      else
        n3_ = 0
        n4_ = 0
      end if
      call ddz_Edge(field_c, dfdz_c, bottom_flag, top_flag, n3_, n4_)
    else
      call ddz_Cell(field_c, dfdz_c, bottom_flag, top_flag)
    end if
    
    !---------------------------------------------------------------------------
    ! Precompute coarse-centre bounds only needed for cell-centred extrapolation.
    ! For centres: z_c(k) = (k-0.5)*dz_c
    ! For edges:   not used (edges nest exactly by construction)
    !---------------------------------------------------------------------------
    if (.not. staggered) then
        zc1 = 0.5_rkind * dz_c
        zcN = (real(nz_c, rkind) - 0.5_rkind) * dz_c
    end if

    do kf = 1, nz_nodes_f

        if (.not. staggered) then
        ! Fine-centre physical location: z_f = (kf-0.5)*dz_f
        zf = (real(kf, rkind) - 0.5_rkind) * dz_f

        !---------------------------
        ! Bottom one-sided extrapolation
        !---------------------------
        if (zf < zc1) then
            z0 = zc1
            do j = 1, ny
            do i = 1, nx
                fbc = field_c(i,j,1)
                mbc = dfdz_c (i,j,1)
                field_f(i,j,kf) = fbc + (zf - z0) * mbc
            end do
            end do
            cycle
        end if

        !---------------------------
        ! Top one-sided extrapolation
        !---------------------------
        if (zf > zcN) then
            z0 = zcN
            do j = 1, ny
            do i = 1, nx
                fbc = field_c(i,j,nz_c)
                mbc = dfdz_c (i,j,nz_c)
                field_f(i,j,kf) = fbc + (zf - z0) * mbc
            end do
            end do
            cycle
        end if
        end if

        !-----------------------------------------------------------------------
        ! Interior mapping (Hermite) using integer quotient+remainder:
        !   q = coarse interval index (0-based)
        !   s = sub-index within interval (0..refine_z-1)
        !   kc = left coarse node index (1-based)
        !   t  = s/refine_z in [0,1)
        !-----------------------------------------------------------------------
        q = (kf - 1) / refine_z
        s = (kf - 1) - q*refine_z

        if (staggered) then
        ! For edges, the very top fine node maps exactly to last coarse node.
        if (q >= nz_c) then
            do j = 1, ny
            do i = 1, nx
                field_f(i,j,kf) = field_c(i,j,nz_nodes_c)
            end do
            end do
            cycle
        end if
        else
        ! For centres, we are guaranteed here to be inside [z_c(1), z_c(nz_c)].
        ! Clamp q so kc+1 is safe.
        if (q > nz_c - 2) q = nz_c - 2
        end if

        kc = q + 1
        t  = real(s, rkind) / real(refine_z, rkind)

        ! Hermite basis
        h00 =  2.0_rkind*t*t*t - 3.0_rkind*t*t + 1.0_rkind
        h10 =        t*t*t - 2.0_rkind*t*t + t
        h01 = -2.0_rkind*t*t*t + 3.0_rkind*t*t
        h11 =        t*t*t -       t*t

        do j = 1, ny
        do i = 1, nx
            f0 = field_c(i,j,kc)
            f1 = field_c(i,j,kc+1)
            m0 = dfdz_c (i,j,kc)
            m1 = dfdz_c (i,j,kc+1)

            field_f(i,j,kf) = h00*f0 + h10*(dz_c*m0) + h01*f1 + h11*(dz_c*m1)
        end do
        end do

    end do

    deallocate(dfdz_c)

  end subroutine refine_z_physical

  subroutine initializeEverything(Lx, Ly, Lz, nx, ny, nz, p_row, p_col, &
                                  NumericalSchemeVert, botWall, TopWall, botBC_Temp, topBC_Temp)
    implicit none
    real(rkind), intent(in) :: Lx, Ly, Lz
    integer, intent(in) :: nx, ny, nz, p_row, p_col
    integer, intent(in) :: NumericalSchemeVert, botWall, TopWall, botBC_Temp, topBC_Temp
    integer :: nx_f, ny_f, nz_f
    real(rkind) :: dx, dy, dz
    
    ! Make sure nx, ny are even for spectral refinement
    if (mod(nx, 2) /= 0) call gracefulExit("nx must be even for spectral refinement.", 101)
    if (mod(ny, 2) /= 0) call gracefulExit("ny must be even for spectral refinement.", 102) 

    ! Calculate refined grid sizes
    nx_f = nx * refine_x
    ny_f = ny * refine_y
    nz_f = nz * refine_z

    !-----------------------------------------------------------------------------
    ! Initialize decomp2d for the original (coarse) grid
    !-----------------------------------------------------------------------------
    call decomp_2d_init(nx, ny, nz, p_row, p_col)

    ! Get local decomposition info for array allocation
    ! Cell-centered grids
    call decomp_info_init(nx, ny, nz, gpC)
    call decomp_info_init(nx_f, ny_f, nz, gpC_XY)
    call decomp_info_init(nx_f, ny_f, nz_f, gpC_XYZ)

    ! Edge grids (for staggered w)
    call decomp_info_init(nx, ny, nz+1, gpE)
    call decomp_info_init(nx_f, ny_f, nz+1, gpE_XY)
    call decomp_info_init(nx_f, ny_f, nz_f+1, gpE_XYZ)

    ! Initialize spectral
    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    call spectC_c%init("x",nx,ny,nz, dx, dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
    call spectE_c%init("x",nx,ny,nz + 1,dx,dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)

    sp_gpC_c => spectC_c%spectdecomp
    sp_gpE_c => spectE_c%spectdecomp
    
    ! Initialize spectral for fine grid
    dx = Lx/real(nx_f,rkind); dy = Ly/real(ny_f,rkind); dz = Lz/real(nz_f,rkind)
    call spectC_f%init("x",nx_f,ny_f,nz_f, dx, dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
    call spectE_f%init("x",nx_f,ny_f,nz_f + 1,dx,dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)

    ! Initialize spectral for horizontally fine grid but still coarse in z
    dx = Lx/real(nx_f,rkind); dy = Ly/real(ny_f,rkind); dz = Lz/real(nz,rkind)
    call spectC_XY%init("x",nx_f,ny_f,nz, dx, dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
    call spectE_XY%init("x",nx_f,ny_f,nz + 1,dx,dy,dz,"FOUR",'2/3rd', dimTransform=2, fixOddball=.false., init_periodicInZ=.false.)
    sp_gpC_XY => spectC_XY%spectdecomp
    sp_gpE_XY => spectE_XY%spectdecomp

    ! PadeOps
    call Pade6opz%init(gpC_XY, sp_gpC_XY, gpE_XY, sp_gpE_XY, dz, NumericalSchemeVert,.false., spectC_XY)

    allocate(cbuffyC(sp_gpC_XY%ysz(1),gpCsp_gpC_XY_XY%ysz(2),sp_gpC_XY%ysz(3)))
    allocate(cbuffyE(sp_gpE_XY%ysz(1),sp_gpE_XY%ysz(2),sp_gpE_XY%ysz(3)))
    allocate(cbuffzC1(sp_gpC_XY%zsz(1),sp_gpC_XY%zsz(2),sp_gpC_XY%zsz(3)))
    allocate(cbuffzC2(sp_gpC_XY%zsz(1),sp_gpC_XY%zsz(2),sp_gpC_XY%zsz(3)))
    allocate(cbuffzE1(sp_gpE_XY%zsz(1),sp_gpE_XY%zsz(2),sp_gpE_XY%zsz(3)))
    
    ! BC Stencils
    call get_boundary_conditions_stencil(botWall, TopWall, botBC_Temp, topBC_Temp)

    ! Allocations
    ! -------------
    ! Coarse grid arrays
    allocate(u_c(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(v_c(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))    
    allocate(w_c(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    
    ! Fine grid arrays (cell-centered)
    allocate(u_f(gpC_XYZ%xsz(1),gpC_XYZ%xsz(2),gpC_XYZ%xsz(3)))
    allocate(v_f(gpC_XYZ%xsz(1),gpC_XYZ%xsz(2),gpC_XYZ%xsz(3)))
    allocate(w_f(gpE_XYZ%xsz(1),gpE_XYZ%xsz(2),gpE_XYZ%xsz(3)))
    
    if(isStratified)then
        allocate(T_c(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
        allocate(T_f(gpC_XYZ%xsz(1),gpC_XYZ%xsz(2),gpC_XYZ%xsz(3)))
    end if
    
    ! Horizontally refined
    allocate(fxy_inX(gpC_XY%xsz(1), gpC_XY%xsz(2), gpC_XY%xsz(3)))
    allocate(fxy_inY(gpC_XY%ysz(1), gpC_XY%ysz(2), gpC_XY%ysz(3)))
    allocate(fxy_inZ(gpC_XY%zsz(1), gpC_XY%zsz(2), gpC_XY%zsz(3)))
    allocate(fxyE_inX(gpE_XY%xsz(1), gpE_XY%xsz(2), gpE_XY%xsz(3)))
    allocate(fxyE_inY(gpE_XY%ysz(1), gpE_XY%ysz(2), gpE_XY%ysz(3)))
    allocate(fxyE_inZ(gpE_XY%zsz(1), gpE_XY%zsz(2), gpE_XY%zsz(3)))

    ! Fully refined
    allocate(fxyz_inY(gpC_XYZ%ysz(1), gpC_XYZ%ysz(2), gpC_XYZ%ysz(3)))
    allocate(fxyz_inZ(gpC_XYZ%zsz(1), gpC_XYZ%zsz(2), gpC_XYZ%zsz(3)))
    allocate(fxyzE_inY(gpE_XYZ%ysz(1), gpE_XYZ%ysz(2), gpE_XYZ%ysz(3)))
    allocate(fxyzE_inZ(gpE_XYZ%zsz(1), gpE_XYZ%zsz(2), gpE_XYZ%zsz(3)))

  end subroutine initializeEverything

  subroutine cleanup()
    implicit none

    deallocate(u_c, v_c, w_c)
    deallocate(u_f, v_f, w_f)
    if (allocated(T_c)) deallocate(T_c)
    if (allocated(T_f)) deallocate(T_f)
    
    deallocate(cbuffyC, cbuffyE, cbuffzC1, cbuffzC2, cbuffzE1)
    deallocate(fxy_inX, fxy_inY, fxy_inZ)
    deallocate(fxyE_inX, fxyE_inY, fxyE_inZ)
    deallocate(fxyz_inY, fxyz_inZ)
    deallocate(fxyzE_inY, fxyzE_inZ)

    call spectC_c%destroy()
    call spectE_c%destroy()
    call spectC_f%destroy()
    call spectE_f%destroy()
    call spectC_XY%destroy()
    call spectE_XY%destroy()
    call Pade6opZ%destroy()

    ! Cell-centered grids
    call decomp_info_finalize(gpC)
    call decomp_info_finalize(gpC_XY)
    call decomp_info_finalize(gpC_XYZ)
    
    ! Edge-based grids
    call decomp_info_finalize(gpE)
    call decomp_info_finalize(gpE_XY)
    call decomp_info_finalize(gpE_XYZ)

    call decomp_2d_finalize()
  end subroutine

  subroutine ddz_Cell(f, dfdz, n1, n2)
    implicit none
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: dfdz
    integer, intent(in) :: n1, n2
    
    call spectC_XY%fft(f, cbuffyC)
    call transpose_y_to_z(cbuffyC, cbuffzC1, spectC_XY%spectdecomp)
    call Pade6opZ%ddz_C2C(cbuffzC1, cbuffzC2, n1, n2)
    call transpose_z_to_y(cbuffzC2, cbuffyC, spectC_XY%spectdecomp)
    call spectC_XY%dealias(cbuffyC)
    call spectC_XY%ifft(cbuffyC, dfdz)
  end subroutine

  subroutine ddz_Edge(f, dfdz, n1, n2, n3, n4)
    implicit none
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: dfdz
    integer, intent(in) :: n1, n2, n3, n4
    
    call spectE_XY%fft(f, cbuffyE)    
    call transpose_y_to_z(cbuffyE, cbuffzE1, spectE_XY%spectdecomp)
    call Pade6opZ%ddz_E2C(cbuffzE1, cbuffzC1, n1, n2)
    call Pade6opZ%interpz_C2E(cbuffzC1, cbuffzE1, n3, n4)
    call transpose_z_to_y(cbuffzE1, cbuffyE, spectE_XY%spectdecomp)
    call spectE_XY%ifft(cbuffyE, dfdz)
  end subroutine

end module refine_fields_mod

program refine_fields
    use refine_fields_mod
    implicit none

    ! Grid parameters
    integer :: nx, ny, nz
    integer :: ierr, ioUnit, p_row=0, p_col=0
    real(rkind) :: Lx, Ly, Lz, dz
    character(len=clen) :: inputfile
    character(len=clen) :: outputdir, inputdir
    integer :: inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID
    integer :: botWall, TopWall, botBC_Temp, topBC_Temp
    integer :: NumericalSchemeVert=1

    namelist /INPUT/ Lx, Ly, Lz, nx, ny, nz, refine_x, refine_y, refine_z, &
        inputdir, outputdir, inputFile_TID, inputFile_RID, &
        outputFile_TID, outputFile_RID, isStratified, p_row, p_col, &
        NumericalSchemeVert, botWall, TopWall, botBC_Temp, topBC_Temp

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    dz = Lz / real(nz, rkind)

    call initializeEverything(Lx, Ly, Lz, nx, ny, nz, p_row, p_col, &
                              NumericalSchemeVert, botWall, TopWall, botBC_Temp, topBC_Temp)

    !----------------------------------------------------------
    ! Read coarse fields from restart files (in X-pencils)
    !----------------------------------------------------------
    call read_restart_file(u_c, inputdir, inputFile_TID, inputFile_RID, '_u.', gpC)
    call read_restart_file(v_c, inputdir, inputFile_TID, inputFile_RID, '_v.', gpC)
    call read_restart_file(w_c, inputdir, inputFile_TID, inputFile_RID, '_w.', gpE)
    if(isStratified)then
        call read_restart_file(T_c, inputdir, inputFile_TID, inputFile_RID, '_T.', gpC)
    end if

    ! Refine cell-centered fields (u, v, T)
    call refine_single_field(u_c, u_f, dz, uBC_bottom, uBC_top)
    call refine_single_field(v_c, v_f, dz, vBC_bottom, vBC_top)
    if(isStratified)call refine_single_field(T_c, T_f, dz, TBC_bottom, TBC_top)
    
    ! Refine w velocity (staggered in z)
    call refine_single_fieldE(w_c, w_f, dz, wBC_bottom, wBC_top, dwdzBC_bottom, dwdzBC_top)
    
    ! Dump to file
    call write_restart_file(u_f, outputdir, outputFile_TID, outputFile_RID, '_u.', gpC_XYZ)
    call write_restart_file(v_f, outputdir, outputFile_TID, outputFile_RID, '_v.', gpC_XYZ)
    call write_restart_file(w_f, outputdir, outputFile_TID, outputFile_RID, '_w.', gpE_XYZ)
    if(isStratified)then
        call write_restart_file(T_f, outputdir, outputFile_TID, outputFile_RID, '_T.', gpC_XYZ)
    end if

    ! Clean up and finalize MPI
    call cleanup()
    call MPI_FINALIZE(ierr)
  
end program refine_fields