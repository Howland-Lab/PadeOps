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

    subroutine doTemporalStuff(gp, simid)
        class(igrid), intent(inout) :: gp 
        integer, intent(in) :: simid

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            select case (simid)
            case (1)
               call message(0,"Primary Simulation Info:")
            case (2)
               call message(0,"Concurrent Simulation Info:")
            end select 
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message(1,"u_star:",gp%sgsmodel%get_ustar())
            call message(1,"Inv. Ob. Len:",gp%sgsmodel%get_InvObLength())
            call message(1,"Surface Flux (K*nd velocity):",gp%wTh_surf)
            call message(1,"T_surf:",gp%sgsmodel%get_T_surf())
            call message(1,"wTh_surf:",gp%sgsmodel%get_wTh_surf())
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))

            if ((simid == 1) .and. (gp%useWindTurbines)) then
                call message(0,"Wind direction hub height", gp%WindTurbineArr%windAngle)
            end if 

            ! add controller print statements, if the controller is used
            if (gp%useControl) then
                call message(1, "Current angle controller Phi:", gp%angCont_yaw%getPhi())
                call message(1, "Frame angle:" , gp%frameAngle)
                call message(1, "Current wind angle:", gp%angCont_yaw%getPhiHub())
            end if

            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call message(0,"------------------------------------------")
            if (simid == 1) then
                if (allocated(gp%scalars)) then
                    call message_min_max(1,"Bounds for SCALAR 1:", p_minval(minval(gp%scalars(1)%F)), p_maxval(maxval(gp%scalars(1)%F)))
                    call message_min_max(1,"Bounds for SCALAR 2:", p_minval(minval(gp%scalars(2)%F)), p_maxval(maxval(gp%scalars(2)%F)))
                    call message_min_max(1,"Bounds for SCALAR 3:", p_minval(minval(gp%scalars(3)%F)), p_maxval(maxval(gp%scalars(3)%F)))
                end if
                
                if (p_maxval(maxval(gp%u))>4.) then
                    call message(1, "this step has blown up", gp%tsim)
                    call gp%dumpFullField(gp%u,"uVel")
                    call gp%dumpFullField(gp%v,"vVel")
                    call gp%dumpFullField(gp%wC,"wVel")
                    call gp%dumpFullField(gp%T, "potT")
                    call gp%dumpFullField(gp%T, "prss")
                    call GracefulExit("u-velocity has blown up",1)
                end if
            elseif (simid == 2) then
                call toc()
                call tic()
            end if 
        end if 

    end subroutine



end module 
