module frozen_igrid_mod
    use IncompressibleGrid, only: igrid
    use exits, only: message
    use kind_parameters, only: rkind
    implicit none
    public :: frozen_igrid

    ! A very simple class which overwrites the time advancement
    ! to keep velocity fields frozen.
    type, extends(igrid) :: frozen_igrid

    contains
        procedure :: timeAdvance => timeAdvanceFrozen
        procedure :: advance_SSP_RK45_Stage_1 => RK45Frozen
        procedure :: advance_SSP_RK45_Stage_2 => RK45Frozen
        procedure :: advance_SSP_RK45_Stage_3 => RK45Frozen
        procedure :: advance_SSP_RK45_Stage_4 => RK45Frozen
        procedure :: advance_SSP_RK45_Stage_5 => RK45Frozen

    end type frozen_igrid

contains

    subroutine timeAdvanceFrozen(this, dtforced)
        class(frozen_igrid), intent(inout) :: this
        real(rkind), intent(in), optional :: dtforced

        call message(0, "frozen_igrid: skipping time advance")

    end subroutine

    subroutine RK45Frozen(this)
        class(frozen_igrid), intent(inout), target :: this
        ! skip everything
    end subroutine

end module
