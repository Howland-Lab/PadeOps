! Template for PadeOps
#include "stratified_pbl_concurrent_files/initialize.F90"       
#include "stratified_pbl_concurrent_files/temporalHook.F90"  

program stratified_pbl_concurrent
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use constants, only: one, zero
    use budgets_time_avg_mod, only: budgets_time_avg  
    use budgets_time_avg_deficit_mod, only: budgets_time_avg_deficit
    use link_turbine_to_scalar, only: setup_turb_scalar_source
    use stratified_pblBCsmod 
    implicit none

    type(igrid), allocatable, target :: primary, precursor 
    character(len=clen) :: inputfile, primary_inputFile, precursor_inputFile
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg, pre_budg_tavg
    type(budgets_time_avg_deficit) :: budg_tavg_deficit
    real(rkind) :: dt1, dt2, dt
    logical :: synchronize_RK_fringe = .true.

    namelist /concurrent/ primary_inputfile, precursor_inputfile, synchronize_RK_fringe

    call MPI_Init(ierr)                                                

    call GETARG(1,inputfile)                                            

    allocate(precursor)                                                       
    allocate(primary)                                                     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)    

    call setup_stratified_pblBCs(primary_inputFile) 

    call primary%init(primary_inputFile, .true.)                               
    call primary%start_io(.true.)                                          
    call primary%printDivergence()
    
    if (primary%useScalars) then
        call setup_turb_scalar_source(primary)
    end if 

    call precursor%init(precursor_inputFile, .false.)                                           
    precursor%Am_I_Primary = .false. 
    call precursor%start_io(.true.)                                           
   
    if (primary%usefringe) then
        call primary%fringe_x%associateFringeTargets(precursor%u, precursor%v, precursor%wC, precursor%T) 
        call primary%fringe_x%associateFringeTarget_scalar(precursor%T)
    end if 

    call budg_tavg%init(primary_inputfile, primary)   !<-- Budget class initialization 
    call pre_budg_tavg%init(precursor_inputfile, precursor)   !<-- Budget class initialization 
    call budg_tavg_deficit%init(pre_budg_tavg, primary_inputfile, budg_tavg)   !<-- Budget class initialization for the deficit 

    if (primary%useWindTurbines) then
        call primary%WindTurbineArr%link_reference_domain_for_control(primary%u, primary%v, primary%rbuffyC, primary%rbuffzC, primary%gpC)
    end if 

    call stratified_pblBCs_CorrectnessCheck(primary%tsim, primary%G_geostrophic)
    call stratified_pblBCs_CorrectnessCheck(precursor%tsim, precursor%G_geostrophic)

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic()
    do while (primary%tsim < primary%tstop) 
       dt1 = primary%get_dt(recompute=.true.)
       dt2 = precursor%get_dt(recompute=.true.)
       dt = min(dt1, dt2)
       call get_stratified_pblBCs(primary%tsim, primary%G_geostrophic, primary%wTh_surf, primary%G_alpha, primary%useControl)
       call get_stratified_pblBCs(precursor%tsim, precursor%G_geostrophic, precursor%wTh_surf, precursor%G_alpha, precursor%useControl)
       if (synchronize_RK_fringe) then
           primary%dt = dt
           precursor%dt = dt
           ! Stage 1
           call primary%advance_SSP_RK45_Stage_1()
           call precursor%advance_SSP_RK45_Stage_1()
           ! Stage 2
           call primary%advance_SSP_RK45_Stage_2()
           call precursor%advance_SSP_RK45_Stage_2()
           ! Stage 3
           call primary%advance_SSP_RK45_Stage_3()
           call precursor%advance_SSP_RK45_Stage_3()
           ! Stage 4
           call primary%advance_SSP_RK45_Stage_4()
           call precursor%advance_SSP_RK45_Stage_4()
           ! Stage 5
           call primary%advance_SSP_RK45_Stage_5()
           call precursor%advance_SSP_RK45_Stage_5()
           ! Call wrap up 
           call primary%wrapup_timestep()
           call precursor%wrapup_timestep() 
       else 
           call primary%timeAdvance(dt)
           call precursor%timeAdvance(dt)
       end if
       
       call budg_tavg%doBudgets()       
       call pre_budg_tavg%doBudgets()       
       call budg_tavg_deficit%doBudgets()       

       call doTemporalStuff(primary,   1)                                        
       call doTemporalStuff(precursor,2)                                        
     
    end do 

    call budg_tavg%destroy()           !<-- release memory taken by the budget class 
    call pre_budg_tavg%destroy()           !<-- release memory taken by the budget class 
    call budg_tavg_deficit%destroy()           !<-- release memory taken by the budget class 

    call precursor%finalize_io()          
    call primary%finalize_io()         

    call precursor%destroy()                
    call primary%destroy()               
   
    deallocate(precursor, primary)             
    
    call MPI_Finalize(ierr)        

end program
