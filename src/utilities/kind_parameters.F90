! Set the floating point precision

module kind_parameters

    use mpi
    use decomp_2d_constants, only : mytype
    implicit none
    
    private
    public :: single_kind, double_kind, rkind, mpirkind, mpickind, clen, stdin, stdout, stderr

    integer, parameter :: single_kind = kind(0.0)
    integer, parameter :: double_kind = mytype  ! kind(0.d0)
    integer, parameter :: rkind = double_kind
    integer, parameter :: mpirkind = MPI_DOUBLE_PRECISION
    integer, parameter :: mpickind = MPI_DOUBLE_COMPLEX
   
    integer, parameter :: clen = 256

    integer, parameter :: stdin  = 5
    integer, parameter :: stdout = 6
    integer, parameter :: stderr = 0

end module
