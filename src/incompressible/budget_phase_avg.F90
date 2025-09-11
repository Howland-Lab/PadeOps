module budgets_phase_avg_mod
    use kind_parameters, only: rkind, clen, mpirkind
   use decomp_2d
   use reductions, only: p_sum
   use incompressibleGrid, only: igrid  
   use exits, only: message, GracefulExit
   use basic_io, only: read_2d_ascii, write_2d_ascii
   use constants, only: half, zero
   use mpi 

   implicit none 

   private
   public :: budgets_phase_avg

    ! BUDGET TYPE: 
    ! BUDGET_0: 6 Reynolds stress terms + 3 temp fluxes + meanU + meanV + meanT

    ! BUDGET_0 term indices:
    ! 1:  <U> 
    ! 2:  <V>
    ! 3:  <W>
    ! 4:  <uu>
    ! 5:  <uv> 
    ! 6:  <uw>
    ! 7:  <vv>
    ! 8:  <vw>
    ! 9:  <ww>
    ! 10: <P>
    ! 11: <tau11> 
    ! 12: <tau12> 
    ! 13: <tau13>
    ! 14: <tau22> 
    ! 15: <tau23> 
    ! 16: <tau33> 
    ! 17: <p'u'>
    ! 18: <p'v'>
    ! 19: <p'w'>
    ! 20: <u'k'>
    ! 21: <v'k'>
    ! 22: <w'k'>
    ! 23: <u_j'tau_1j'>
    ! 24: <u_j'tau_2j'>
    ! 25: <u_j'tau_3j'>
    ! 26: <T>
    ! 27: <uT>
    ! 28: <vT>
    ! 29: <wT>
    ! 30: <TT>
    ! 31: <wb> buoyancy term for TKE and uiuj budgets
    ! 32: Scalar means (32-> 32+num_scalars)
    ! 32+x: Scalar variances(32+num_scalars+1:33+2*num_scalars)

    type :: budgets_phase_avg

    contains

    end type

end module 