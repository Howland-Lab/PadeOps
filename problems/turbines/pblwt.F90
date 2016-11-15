! Template for PadeOps
! Grid - hit_grid
! Problem - HIT

#include "pblwt_files/initialize.F90"       
#include "pblwt_files/io.F90"            
#include "pblwt_files/temporalHook.F90"  

program pblwt
    use mpi
    use kind_parameters,  only: rkind,clen,stdout,stderr
    use IncompressibleGridWallM, only: igridWallM
    use pblwt_IO, only: start_io, finalize_io
    use constants, only: half 
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message

    implicit none

    type(igridWallM), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call start_io(igp)                !<-- Start I/O by creating a header file (see io.F90)

    call igp%printDivergence()
  
    call tic() 
    do while (igp%tsim < igp%tstop) 
       
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
    call finalize_io                  !<-- Close the header file (wrap up i/o)

    !call igp%finalize_stats()
    call igp%finalize_stats3D()
    
    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program