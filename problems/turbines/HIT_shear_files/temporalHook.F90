module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc
    use mpi

    implicit none

    integer :: nt_print2screen = 1
    real(rkind) :: maxDiv, DomMaxDiv

contains

    subroutine doTemporalStuff(gp, simid)
        class(igrid), intent(inout) :: gp
        integer, intent(in) :: simid

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            select case (simid)
              case (1)
                call message(0,"Actuator Disk Simulation Info:")
              case (2)
                call message(0,"HIT Simulation Info:")
              case (0)
                call message(0, "Empty Simulation Info:")  ! for the HIT_AD_deficit problem
            end select
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
                ! call message(1, "dt limit: "//gp%dtlimit)
            end if
            call message(0,"------------------------------------------")
            if (simid == 1) then
                if (allocated(gp%scalars)) then
                    call message_min_max(1,"Bounds for SCALAR 1:", p_minval(minval(gp%scalars(1)%F)), p_maxval(maxval(gp%scalars(1)%F)))
                    call message_min_max(1,"Bounds for SCALAR 2:", p_minval(minval(gp%scalars(2)%F)), p_maxval(maxval(gp%scalars(2)%F)))
                end if
                ! check against blowing up
                if (p_maxval(maxval(gp%u))>10.) then
                    call message(1, "this step has blown up", gp%tsim)
                    call gp%dumpFullField(gp%u,"uVel")
                    call gp%dumpFullField(gp%v,"vVel")
                    call gp%dumpFullField(gp%wC,"wVel")
                    call gp%dumpFullField(gp%T, "potT")
                    call GracefulExit("u-velocity has blown up",1)
                end if

            elseif (simid == 2) then
                call message(1,"Mean TKE for HIT:", gp%getMeanKE())
                call message(1,"Mean  uu for HIT:", gp%getMeanuu())
                call message(1,"Mean  uv for HIT:", gp%getMeanuv())
                call message(1,"Mean  uw for HIT:", gp%getMeanuw())
                call message(1,"Mean  vv for HIT:", gp%getMeanvv())
                call message(1,"Mean  vw for HIT:", gp%getMeanvw())
                call message(1,"Mean  ww for HIT:", gp%getMeanww())
                call toc()
                call tic()
            end if
        end if

    end subroutine


end module
