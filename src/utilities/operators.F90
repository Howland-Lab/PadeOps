module operators
    use kind_parameters, only: rkind, clen
    use FiltersMod, only: filters
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
    use exits,           only: GracefulExit
   
    implicit none

contains

    subroutine gradient(decomp, der, f, dfdx, dfdy, dfdz)
        type(decomp_info), intent(in) :: decomp
        type(derivatives) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: f
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdx
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdy
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdz

        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        if ( (size(f,1) .NE. decomp%ysz(1)) .AND. (size(f,2) .NE. decomp%ysz(2)) .AND. (size(f,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit("Either size of input array to gradient operator is inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if

        ! Get Y derivatives
        call der%ddy(f,dfdy)

        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,decomp)
        call der%ddx(xtmp,xdum)
        call transpose_x_to_y(xdum,dfdx)

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,decomp)
        call der%ddz(ztmp,zdum)
        call transpose_z_to_y(zdum,dfdz)

    end subroutine 

    subroutine curl(decomp, der, u, v, w, curlu)
        type(decomp_info), intent(in) :: decomp
        type(derivatives) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: u, v, w
        real(rkind), dimension(size(u,1), size(u,2), size(u,3),3),             intent(out) :: curlu

        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp
        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        if (    (size(u,1) .NE. decomp%ysz(1)) .AND. (size(u,2) .NE. decomp%ysz(2)) .AND. (size(u,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit("Either size of input array 'u' to the curl operator is inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if
        if (    (size(v,1) .NE. decomp%ysz(1)) .AND. (size(v,2) .NE. decomp%ysz(2)) .AND. (size(v,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit("Either size of input array 'v' to the curl operator is inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if
        if (    (size(w,1) .NE. decomp%ysz(1)) .AND. (size(w,2) .NE. decomp%ysz(2)) .AND. (size(w,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit("Either size of input array 'w' to the curl operator is inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if

        ! Get dw/dy
        call der%ddy( w, curlu(:,:,:,1) )

        ! Get dv/dz
        call transpose_y_to_z( v, ztmp, decomp)
        call der%ddz(ztmp,zdum)
        call transpose_z_to_y(zdum,ytmp)

        curlu(:,:,:,1) = curlu(:,:,:,1) - ytmp ! dw/dy - dv/dz

        ! Get du/dz
        call transpose_y_to_z( u, ztmp, decomp)
        call der%ddz(ztmp,zdum)
        call transpose_z_to_y(zdum, curlu(:,:,:,2) )

        ! Get dw/dx
        call transpose_y_to_x( w, xtmp, decomp)
        call der%ddx(xtmp,xdum)
        call transpose_x_to_y(xdum, ytmp )

        curlu(:,:,:,2) = curlu(:,:,:,2) - ytmp ! du/dz - dw/dx

        ! Get dv/dx
        call transpose_y_to_x( v, xtmp, decomp)
        call der%ddx(xtmp,xdum)
        call transpose_x_to_y(xdum, curlu(:,:,:,3) )

        ! Get du/dy
        call der%ddy( u, ytmp )

        curlu(:,:,:,3) = curlu(:,:,:,3) - ytmp ! dv/dx - du/dy

    end subroutine 

    subroutine divergence(decomp, der, u, v, w, div)
        type(decomp_info), intent(in) :: decomp
        type(derivatives) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: u, v, w
        real(rkind), dimension(size(u,1), size(u,2), size(u,3)),             intent(out) :: div

        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp
        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        if ( (size(u,1) .NE. decomp%ysz(1)) .AND. (size(u,2) .NE. decomp%ysz(2)) .AND. (size(u,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit("Either size of input array to divergence operator is inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if

        ! Get Y derivatives
        call der%ddy(v,div)

        ! Get X derivatives
        call transpose_y_to_x(u,xtmp,decomp)
        call der%ddx(xtmp,xdum)
        call transpose_x_to_y(xdum,ytmp)

        div = div + ytmp

        ! Get Z derivatives
        call transpose_y_to_z(w,ztmp,decomp)
        call der%ddz(ztmp,zdum)
        call transpose_z_to_y(zdum,ytmp)

        div = div + ytmp

    end subroutine 

end module