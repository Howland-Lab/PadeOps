! Template for PadeOps

#include "stratified_pbl_files/initialize.F90"       
#include "stratified_pbl_files/temporalHook.F90"  

program stratified_pbl
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use budgets_time_avg_mod, only: budgets_time_avg
    use stratified_pblBCsmod 
    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    type(budgets_time_avg) :: budg_tavg


    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call setup_stratified_pblBCs(inputfile) 
    
    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.true.)                !<-- Start I/O by creating a header file (see io.F90)
    
    call igp%printDivergence()
   
    call stratified_pblBCs_CorrectnessCheck(igp%tsim, igp%G_geostrophic)

    call budg_tavg%init(inputfile, igp) !<-- Budget class initialization 

    call tic() 
    do while (igp%tsim < igp%tstop) 

       call get_stratified_pblBCs(igp%tsim, igp%G_geostrophic, igp%wTh_surf, igp%G_alpha, igp%useControl)
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       call budg_tavg%doBudgets()       !<--- perform budget related operations 
       
    end do 
    
    call budg_tavg%destroy()             !<-- release memory taken by the budget class 
   
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
