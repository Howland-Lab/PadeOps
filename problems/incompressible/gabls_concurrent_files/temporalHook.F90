module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max, GracefulExit
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use decomp_2d
    use reductions, only: p_sum

    implicit none 

    integer :: i, j, nt_print2screen = 1
    real(rkind) ::  maxDiv, DomMaxDiv, angle
    integer :: ierr 

contains

    subroutine doTemporalStuff(igp, simid)
        class(igrid), intent(inout) :: igp
        integer, intent(in) :: simid

        if (mod(igp%step,nt_print2screen) == 0) then
            maxDiv = maxval(igp%divergence)
            DomMaxDiv = p_maxval(maxDiv)

            select case (simid)
            case (1)
               call message(0,"Primary Simulation Info:")
            case (2)
               call message(0,"Concurrent Simulation Info:")
            end select

            call message(0,"Time",igp%tsim)
            call message(1,"u_star:",igp%sgsmodel%get_ustar())
            call message(1,"TIDX:",igp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message(1,"Inv. Ob. Len:",igp%sgsmodel%get_InvObLength())
            call message(1,"Surface Flux (K*nd velocity):",igp%wTh_surf)
            call message_min_max(1,"Bounds for u:", p_minval(minval(igp%u)), p_maxval(maxval(igp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(igp%v)), p_maxval(maxval(igp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(igp%w)), p_maxval(maxval(igp%w)))
            call message_min_max(1,"Bounds for T:", p_minval(minval(igp%T)), p_maxval(maxval(igp%T)))
            call message(1,"Control Galpha:", igp%G_alpha)
            call message(1,"frameAngle:",igp%frameAngle)
            if (igp%useCFL) then
                call message(1,"Current dt:",igp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if 

    end subroutine

end module 
