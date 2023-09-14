module budgets_time_avg_deficit_mod
    use kind_parameters, only: rkind, clen, mpirkind
    use decomp_2d
    use reductions, only: p_sum
    use incompressibleGrid, only: igrid
    use budgets_time_avg_mod, only: budgets_time_avg
    use exits, only: message, GracefulExit
    use basic_io, only: read_2d_ascii, write_2d_ascii
    use constants, only: half, zero
    use mpi
 
    implicit none 
 
    private
    public :: budgets_time_avg_deficit
 
    ! BUDGET TYPE: 
    ! BUDGET_0: 6 Reynolds stress terms + 3 temp fluxes + meanU + meanV + meanT
    ! BUDGET_1: momentum equation terms (Budget0 also computed) 
    ! BUDGET_2: MKE budget (Budget0 and Budget1 also computed) 
    ! BUDGET 3: TKE budget (Budget 0, 1 and 2 also included)
 
 
    ! BUDGET_0 term indices:
    ! 1:  <Delta U> 
    ! 2:  <Delta V>
    ! 3:  <Delta W>
    ! 4:  <Delta P>
    ! 5:  <Delta u Delta u>
    ! 6:  <Delta u Delta v> 
    ! 7:  <Delta u Delta w>
    ! 8:  <Delta v Delta v>
    ! 9:  <Delta v Delta w>
    ! 10: <Delta w Delta w>
    ! 11: <Delta u base u>
    ! 12: <Delta u base v> 
    ! 13: <base u Delta v> 
    ! 14: <Delta u base w>
    ! 15: <base u Delta w>  
    ! 16: <Delta v base v>
    ! 17: <Delta v base w>
    ! 18: <base v Delta w>    
    ! 19: <Delta w base w>
    ! 20: delta tau11
    ! 21: delta tau12
    ! 22: delta tau13
    ! 23: delta tau22
    ! 24: delta tau23
    ! 25: delta tau33
    ! 26: <Delta T>
    ! 27: <Delta u Delta T>
    ! 28: <Delta v Delta T>
    ! 29: <Delta w Delta T>
    ! 30: <Delta T Delta T>
 
    ! BUDGET_1 term indices:  
    ! 1:  X eqn - total advection in deficit equation
    ! 2:  X eqn - mean delta advection by base flow
    ! 3:  X eqn - mean delta advection by deficit flow
    ! 4:  X eqn - mean base advection by deficit flow
    ! 5:  X eqn - fluctuating delta advection by base flow
    ! 6:  X eqn - fluctuating delta advection by deficit flow
    ! 7:  X eqn - fluctuating base advection by deficit flow
    ! 8:  X eqn - Pressure gradient term
    ! 9:  X eqn - SGS term
    ! 10: X eqn - Actuator disk/turbine term 
    
    ! 11: Y eqn - total advection in deficit equation  
    ! 12: Y eqn - mean delta advection by base flow
    ! 13: Y eqn - mean delta advection by deficit flow
    ! 14: Y eqn - mean base advection by deficit flow
    ! 15: Y eqn - fluctuating delta advection by base flow
    ! 16: Y eqn - fluctuating delta advection by deficit flow
    ! 17: Y eqn - fluctuating base advection by deficit flow
    ! 18: Y eqn - Pressure gradient term
    ! 19: Y eqn - SGS term
    ! 20: Y eqn - Actuator disk/turbine term 
 
    ! 21: Z eqn - total advection in deficit equation  
    ! 22: Z eqn - mean delta advection by base flow
    ! 23: Z eqn - mean delta advection by deficit flow
    ! 24: Z eqn - mean base advection by deficit flow
    ! 25: Z eqn - fluctuating delta advection by base flow
    ! 26: Z eqn - fluctuating delta advection by deficit flow
    ! 27: Z eqn - fluctuating base advection by deficit flow
    ! 28: Z eqn - Pressure gradient term 
    ! 29: Z eqn - SGS term 
    ! 30: Z eqn - Buoyancy term
    
    ! Coriolis terms
    ! 31:  X eqn - Coriolis Term 
    ! 32:  X eqn - Geostrophic Forcing Term 
    ! 33:  Y eqn - Coriolis Term 
    ! 34:  Y eqn - Geostrophic Forcing Term 
 
 
    ! BUDGET_2 term indices: 
    ! 1:  Loss to Resolved TKE      
    ! 2:  Advective transport       
    ! 3:  Reynolds stress transport 
    ! 4:  Pressure transport        
    ! 5:  SGS + viscous transport   
    ! 6:  Loss to SGS TKE + viscous dissipation 
    ! 7:  Actuator disk sink       
    ! 8:  Geostrophic              
    ! 9:  Coriolis   
    ! 10: Buoyancy       
    ! 11: Loss to Resolved TKE - <delta ui' delta uj'>
    ! 12: Loss to Resolved TKE - <delta ui' base uj'>
    ! 13: Loss to Resolved TKE - <base ui' delta uj'>
    ! 14: Turbulent transport  - <delta ui' delta uj'>
    ! 15: Turbulent transport  - <delta ui' base uj'>
    ! 16: Turbulent transport  - <base ui' delta uj'>
    ! 17: Mean delta MKE advection by base flow
    ! 18: Mean delta MKE advection by deficit flow
    ! 19: Mean base advection by deficit flow 
    
    ! BUDGET_3 term indices:
    ! 1. TKE production              (G)
    ! 2. advective transport        (B)
    ! 3. turbulent transport         (C)
    ! 4. Pressure transport          (D)
    ! 5. SGS + viscous transport     (E+F)
    ! 6. SGS + viscous dissipation   (H+I)
    ! 7. Actuator disk/Turbine sink  (J)
    ! 8. Buoyancy
    ! 9:  TKE production - <delta ui' delta uj'> d delta ui / dxj
    ! 10: TKE production - <delta ui' base uj'> d delta ui / dxj
    ! 11: TKE production - <delta ui' delta uj'> d base ui / dxj
    ! 12: Turbulent transport  - <delta ui' delta uj'>
    ! 13: Turbulent transport  - <delta ui' base uj'>
    ! 14: Turbulent transport  - <base ui' delta uj'>
    ! 15: Delta TKE advection by base flow
    ! 16: Delta TKE advection by deficit flow
    ! 17: Base advection by deficit flow 
    ! 18: Base advection by deficit flow - x
    ! 19: Base advection by deficit flow - y
    ! 20: Base advection by deficit flow - z

    ! BUDGET_4_ij term indices: 
    ! 1. <u'u'> : Shear Production           - N/A,      dump
    ! 2. <u'u'> : convective transport       - N/A,      dump
    ! 3. <u'u'> : Turbulent transport        - assemble, dump
    ! 4. <u'u'> : Pressure strain            - assemble, dump
    ! 5. <u'u'> : Pressure transport         - assemble, dump
    ! 6. <u'u'> : SGS + Viscous transport    - assemble, dump
    ! 7. <u'u'> : SGS + Viscous dissipation  - assemble, dump
    ! 8. <u'u'> : Buoyancy contribution      - N/A,      N/A
    ! 9. <u'u'> : Coriolis exchange          - assemble, dump
    ! 10. <u'u'> : Actuator disk sink/source - assemble, dump
 
    ! 11. <v'v'> : Shear Production          - N/A,      dump 
    ! 12. <v'v'> : convective transport      - N/A,      dump 
    ! 13. <v'v'> : Turbulent transport       - assemble, dump
    ! 14. <v'v'> : Pressure strain           - assemble, dump
    ! 15. <v'v'> : Pressure transport        - assemble, dump
    ! 16. <v'v'> : SGS + Viscous transport   - assemble, dump
    ! 17. <v'v'> : SGS + Viscous dissipation - assemble, dump
    ! 18. <v'v'> : Buoyancy contribution     - N/A,      N/A
    ! 19. <v'v'> : Coriolis exchange         - assemble, dump
    ! 20. <v'v'> : Actuator disk sink/source - assemble, dump
 
    ! 21. <u'w'> : Shear Production          - N/A,      dump
    ! 22. <u'w'> : convective transport      - N/A,      dump
    ! 23. <u'w'> : Turbulent transport       - assemble, dump
    ! 24. <u'w'> : Pressure strain           - assemble, dump 
    ! 25. <u'w'> : Pressure transport        - assemble, dump
    ! 27. <u'w'> : SGS + Viscous transport   - assemble, 
    ! 28. <u'w'> : SGS + Viscous dissipation - assemble, 
    ! 29. <u'w'> : Buoyancy contribution     - assemble,
    ! 30. <u'w'> : Coriolis exchange         - assemble,
    ! 31. <u'w'> : Actuator disk sink/source - assemble, 
 
    ! 19. <v'w'> : Shear Production          - N/A,      dump
    ! 22. <v'w'> : convective transport      - N/A,      dump
    ! 20. <v'w'> : Turbulent transport       - assemble, dump
    ! 21. <v'w'> : Pressure strain           - assemble, dump
    ! 22. <v'w'> : Pressure transport        - assemble, dump
    ! 23. <v'w'> : SGS + Viscous transport   - assemble, dump
    ! 24. <v'w'> : SGS + Viscous dissipation - assemble, dump
    ! 25. <v'w'> : Buoyancy contribution     - assemble, dump
    ! 26. <v'w'> : Coriolis exchange         - assemble, dump
    ! 27. <v'w'> : Actuator disk sink/source - assemble, dump
 
    ! 28. <w'w'> : Shear Production          - N/A,      dump
    ! 22. <w'w'> : convective transport      - N/A,      dump
    ! 29. <w'w'> : Turbulent transport       - assemble, dump
    ! 30. <w'w'> : Pressure strain           - assemble, dump
    ! 31. <w'w'> : Pressure transport        - assemble, dump
    ! 32. <w'w'> : SGS + Viscous transport   - assemble, dump
    ! 33. <w'w'> : SGS + Viscous dissipation - assemble, dump
    ! 34. <w'w'> : Buoyancy contribution     - assemble, dump
    ! 35. <w'w'> : Coriolis exchange         - assemble, dump
    ! 36. <w'w'> : Actuator disk sink/source - N/A,      dump
 
    type :: budgets_time_avg_deficit
         private
         integer :: budgetType = 1, run_id, nz
 
         type(budgets_time_avg), pointer :: pre_budget, prim_budget
         
         real(rkind), dimension(:,:,:,:), allocatable :: budget_0, budget_1, budget_2, budget_3, budget_4_11, budget_4_22, budget_4_33, budget_4_13, budget_4_23
         integer :: counter
         character(len=clen) :: budgets_dir
 
         logical :: useWindTurbines, isStratified, useCoriolis
         real(rkind), allocatable, dimension(:) :: runningSum_sc, runningSum_sc_turb, runningSum_turb
         logical :: HaveScalars
         integer :: tidx_dump 
         integer :: tidx_compute
         integer :: tidx_budget_start 
         real(rkind) :: time_budget_start 
         logical :: do_budgets
         logical :: forceDump
         logical :: splitPressureDNS
 
     contains
         procedure           :: init        
         procedure           :: destroy
         procedure           :: ResetBudget
         procedure           :: RestartBudget
         procedure, private  :: restart_budget_field
         procedure           :: DoBudgets
         
         procedure, private  :: updateBudget
         procedure, private  :: DumpBudget
         procedure, private  :: dump_budget_field 
         
         procedure, private  :: AssembleBudget0
         procedure, private  :: DumpBudget0 
         
         procedure, private  :: AssembleBudget1
         procedure, private  :: DumpBudget1 
 
         procedure, private  :: AssembleBudget2
         procedure, private  :: DumpBudget2 
        
         procedure, private  :: AssembleBudget3
         procedure, private  :: DumpBudget3 
 
         procedure, private :: ddx_R2R
         procedure, private :: ddy_R2R
         procedure, private :: ddz_R2R
         procedure, private :: ddx_C2R
         procedure, private :: ddy_C2R
         procedure, private :: ddz_C2R
         procedure, private :: interp_Edge2Cell
         procedure, private :: interp_Cell2Edge
         procedure, private :: multiply_CellFieldsOnEdges
     end type 
 
 
 contains 
 
     subroutine init(this, pre_budget, primary_inputfile, prim_budget) 
         class(budgets_time_avg_deficit), intent(inout) :: this
         character(len=*), intent(in) :: primary_inputfile 
         type(budgets_time_avg), intent(inout), target :: pre_budget, prim_budget
         
         character(len=clen) :: budgets_dir = "NULL"
         character(len=clen) :: restart_dir = "NULL"
         integer :: ioUnit, ierr,  budgetType = 1, restart_tid = 0, restart_rid = 0, restart_counter = 0
         logical :: restart_budgets = .false. 
         integer :: tidx_compute = 1000000, tidx_dump = 1000000, tidx_budget_start = -100
         real(rkind) :: time_budget_start = -1.0d0
         logical :: do_budgets = .false. 
         namelist /BUDGET_TIME_AVG_DEFICIT/ budgetType, budgets_dir, restart_budgets, restart_dir, restart_rid, restart_tid, restart_counter, tidx_dump, tidx_compute, do_budgets, tidx_budget_start, time_budget_start
         
         restart_dir = "NULL"

         ! STEP 1: Read in inputs, link pointers and allocate budget vectors
         ioUnit = 534
         open(unit=ioUnit, file=trim(primary_inputfile), form='FORMATTED', iostat=ierr)
         read(unit=ioUnit, NML=BUDGET_TIME_AVG_DEFICIT)
         close(ioUnit)
 
         this%pre_budget => pre_budget 
         this%prim_budget => prim_budget
         this%run_id = this%prim_budget%igrid_sim%runid
         this%nz = this%prim_budget%igrid_sim%nz
         this%do_budgets = do_budgets
         this%tidx_dump = tidx_dump
         this%tidx_compute = tidx_compute
         this%tidx_budget_start = tidx_budget_start  
         this%time_budget_start = time_budget_start  
         this%useWindTurbines = this%prim_budget%igrid_sim%useWindTurbines
         this%isStratified    = this%prim_budget%igrid_sim%isStratified
         this%useCoriolis    = this%prim_budget%igrid_sim%useCoriolis
         this%forceDump = .false.
 
         this%budgets_dir = budgets_dir
         this%budgetType = budgetType
 
         this%splitPressureDNS = this%prim_budget%igrid_sim%computeDNSPressure
 
         this%HaveScalars = this%prim_budget%igrid_sim%useScalars
 
         if((this%tidx_budget_start > 0) .and. (this%time_budget_start > 0.0d0)) then
             call GracefulExit("Both tidx_budget_start and time_budget_start in budget_time_avg are positive. Turn one negative", 100)
         endif
 
         if(this%do_budgets) then 
             !if (this%isStratified) then
             ! Always assume that you are stratified
 
                 if (this%HaveScalars) then
                     allocate(this%budget_0(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),30+2*this%prim_budget%igrid_sim%n_scalars))
                 else
                     allocate(this%budget_0(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),30))
                 end if
                 allocate(this%budget_2(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),19))
                 allocate(this%budget_1(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),34))
             !else
             !    allocate(this%budget_0(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),25))
             !    allocate(this%budget_2(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),07))
             !    allocate(this%budget_1(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),10))
             !end if
             allocate(this%budget_3(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),20))
             allocate(this%budget_4_11(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),10))
             allocate(this%budget_4_22(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),10))
             allocate(this%budget_4_13(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),10))
             allocate(this%budget_4_23(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),10))
             allocate(this%budget_4_33(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3),10))
 
             if ((trim(budgets_dir) .eq. "null") .or.(trim(budgets_dir) .eq. "NULL")) then 
                this%budgets_dir = this%prim_budget%igrid_sim%outputDir
             end if 


            if ((trim(restart_dir) .eq. "null") .or.(trim(restart_dir) .eq. "NULL")) then
                restart_dir = this%budgets_dir
            end if 

 
             if (restart_budgets) then
                 call message(0,"Budget deficit  restart")
                 call this%RestartBudget(restart_dir, restart_rid, restart_tid, restart_counter)
             else
                 call this%resetBudget()
             end if
 
             ! ! STEP 4: For horizontally-averaged surface quantities (called Scalar here), and turbine statistics
             ! allocate(this%inst_horz_avg(5)) ! [ustar, uw, vw, Linv, wT]
             ! allocate(this%runningSum_sc(5))
             ! this%runningSum_sc = zero
             ! if(this%useWindTurbines) then
             !     allocate(this%runningSum_sc_turb(8*this%prim_budget%igrid_sim%WindTurbineArr%nTurbines))
             !     allocate(this%runningSum_turb   (8*this%prim_budget%igrid_sim%WindTurbineArr%nTurbines))
             !     this%runningSum_sc_turb = zero
             !     this%runningSum_turb = zero
             ! endif
 
         end if
 
     end subroutine 
 
 
     subroutine doBudgets(this, forceDump)
         class(budgets_time_avg_deficit), intent(inout) :: this
         logical, intent(in), optional :: forceDump
 
         if(present(forceDump)) then
             this%forceDump = forceDump
         endif
 
         if (this%do_budgets)  then
             if( ( (this%tidx_budget_start>0) .and. (this%prim_budget%igrid_sim%step>this%tidx_budget_start) ) .or. &
                 ( (this%time_budget_start>0) .and. (this%prim_budget%igrid_sim%tsim>this%time_budget_start) ) ) then
         
                 if (mod(this%prim_budget%igrid_sim%step,this%tidx_compute) .eq. 0) then
                     call this%updateBudget()
                 end if
 
                 if ((mod(this%prim_budget%igrid_sim%step,this%tidx_dump) .eq. 0) .or. this%forceDump) then
                     call this%dumpBudget()
                     call message(0,"Dumped a deficit budget .stt file")
                 end if 
             end if 
         end if 
 
         this%forceDump = .false. ! reset to default value
 
     end subroutine 
 
     subroutine updateBudget(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
    
         call this%prim_budget%igrid_sim%getMomentumTerms()  
         call this%pre_budget%igrid_sim%getMomentumTerms()  
 
         select case (this%budgetType)
         case(0)
             call this%AssembleBudget0()
         case(1)
             call this%AssembleBudget0()
             call this%AssembleBudget1()
         case(2)
             call this%AssembleBudget0()
             call this%AssembleBudget1()
         case(3)
            call this%AssembleBudget0()
            call this%AssembleBudget1()
            call this%AssembleBudget3()
         end select
 
         this%counter = this%counter + 1
 
     end subroutine 
 
 
     subroutine DumpBudget(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
         
         ! MKE budget is only assembled before dumping
        if (this%budgetType>1) call this%AssembleBudget2() 
        
         ! Budget 0: 
         call this%dumpbudget0()
 
         ! Budget 1: 
         if (this%budgetType>0) then
             call this%dumpbudget1()
         end if
 
         if (this%budgetType>1) then
             call this%dumpbudget2()
         end if 

         if (this%budgetType>2) then
            call this%dumpbudget3()
        end if 
 
     end subroutine 
 
     ! ---------------------- Budget 0 ------------------------
     subroutine DumpBudget0(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
         integer :: idx 
         
         ! Step 1: Get the average from sum
         this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
         this%pre_budget%budget_0 = this%pre_budget%budget_0/(real(this%counter,rkind) + 1.d-18) 
         
         ! Step 2: Get the <Rij> from <ui uj>
         this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
         this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
         this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
         this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
         this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
         this%budget_0(:,:,:,10)  = this%budget_0(:,:,:,10)  - this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
 
         ! Step 3: Get the mixed <Rij> from mixed <ui uj>
         this%budget_0(:,:,:,11)  = this%budget_0(:,:,:,11)  - this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1) ! R11
         this%budget_0(:,:,:,12)  = this%budget_0(:,:,:,12)  - this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2) ! R12
         this%budget_0(:,:,:,13)  = this%budget_0(:,:,:,13)  - this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1) ! R12
         this%budget_0(:,:,:,14)  = this%budget_0(:,:,:,14)  - this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3) ! R13
         this%budget_0(:,:,:,15)  = this%budget_0(:,:,:,15)  - this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1) ! R13
         this%budget_0(:,:,:,16)  = this%budget_0(:,:,:,16)  - this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2) ! R22
         this%budget_0(:,:,:,17)  = this%budget_0(:,:,:,17)  - this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3) ! R23
         this%budget_0(:,:,:,18)  = this%budget_0(:,:,:,18)  - this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2) ! R23
         this%budget_0(:,:,:,19)  = this%budget_0(:,:,:,19)  - this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3) ! R33
        
        ! STEP 4: Compute Delta heat fluxes
        if (this%isStratified) then
            this%budget_0(:,:,:,27) = this%budget_0(:,:,:,27) - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,28) = this%budget_0(:,:,:,28) - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,29) = this%budget_0(:,:,:,29) - this%budget_0(:,:,:,3)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,30) = this%budget_0(:,:,:,30) - this%budget_0(:,:,:,26)*this%budget_0(:,:,:,26)
        end if 


         ! Step 7: Dump the full budget 
         do idx = 1,size(this%budget_0,4)
             call this%dump_budget_field(this%budget_0(:,:,:,idx),idx,0)
         end do 
         
        ! STEP 4: Compute Delta heat fluxes
         if (this%isStratified) then
            this%budget_0(:,:,:,27) = this%budget_0(:,:,:,27) + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,28) = this%budget_0(:,:,:,28) + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,29) = this%budget_0(:,:,:,29) + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,30) = this%budget_0(:,:,:,30) + this%budget_0(:,:,:,26)*this%budget_0(:,:,:,26)
        end if 

         ! Step 9: Go back to <ui uj> from <Rij>
         this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
         this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
         this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
         this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
         this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
         this%budget_0(:,:,:,10)  = this%budget_0(:,:,:,10)  + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
         
         this%budget_0(:,:,:,11)  = this%budget_0(:,:,:,11)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1) ! R11
         this%budget_0(:,:,:,12)  = this%budget_0(:,:,:,12)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2) ! R12
         this%budget_0(:,:,:,13)  = this%budget_0(:,:,:,13)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1) ! R12
         this%budget_0(:,:,:,14)  = this%budget_0(:,:,:,14)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3) ! R13
         this%budget_0(:,:,:,15)  = this%budget_0(:,:,:,15)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1) ! R13
         this%budget_0(:,:,:,16)  = this%budget_0(:,:,:,16)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2) ! R22
         this%budget_0(:,:,:,17)  = this%budget_0(:,:,:,17)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3) ! R23
         this%budget_0(:,:,:,18)  = this%budget_0(:,:,:,18)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2) ! R23
         this%budget_0(:,:,:,19)  = this%budget_0(:,:,:,19)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3) ! R33
        
 
         ! Step 11: Go back to summing instead of averaging
         this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)
         this%pre_budget%budget_0 = this%pre_budget%budget_0*(real(this%counter,rkind) + 1.d-18) 
 
     end subroutine 
 
     subroutine AssembleBudget0(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
 
         ! STEP 1: Compute mean Delta U, Delta V, and Delta W
         this%budget_0(:,:,:,1) = this%budget_0(:,:,:,1) + this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u
         this%budget_0(:,:,:,2) = this%budget_0(:,:,:,2) + this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v
         this%budget_0(:,:,:,3) = this%budget_0(:,:,:,3) + this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC
 
         ! STEP 2: Pressure
         this%budget_0(:,:,:,4) = this%budget_0(:,:,:,4) + this%prim_budget%igrid_sim%pressure - this%pre_budget%igrid_sim%pressure
 
         ! STEP 3: Get Reynolds stresses (IMPORTANT: need to correct for fluctuation before dumping)
         this%budget_0(:,:,:,5) = this%budget_0(:,:,:,5) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) * (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
         this%budget_0(:,:,:,6) = this%budget_0(:,:,:,6) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) * (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
         this%budget_0(:,:,:,7) = this%budget_0(:,:,:,7) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) * (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
         this%budget_0(:,:,:,8) = this%budget_0(:,:,:,8) + (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v) * (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
         this%budget_0(:,:,:,9) = this%budget_0(:,:,:,9) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) * (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
         this%budget_0(:,:,:,10) = this%budget_0(:,:,:,10) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) * (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
         
         ! STEP 4: Get mixed Reynolds stresses
         this%budget_0(:,:,:,11) = this%budget_0(:,:,:,11) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) * (this%pre_budget%igrid_sim%u) 
         this%budget_0(:,:,:,12) = this%budget_0(:,:,:,12) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) * (this%pre_budget%igrid_sim%v)
         this%budget_0(:,:,:,13) = this%budget_0(:,:,:,13) + (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v) * (this%pre_budget%igrid_sim%u)
         this%budget_0(:,:,:,14) = this%budget_0(:,:,:,14) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) * (this%pre_budget%igrid_sim%wC)
         this%budget_0(:,:,:,15) = this%budget_0(:,:,:,15) + (this%pre_budget%igrid_sim%u) * (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
         this%budget_0(:,:,:,16) = this%budget_0(:,:,:,16) + (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v) * (this%pre_budget%igrid_sim%v)
         this%budget_0(:,:,:,17) = this%budget_0(:,:,:,17) + (this%pre_budget%igrid_sim%wC) * (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
         this%budget_0(:,:,:,18) = this%budget_0(:,:,:,18) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) * (this%pre_budget%igrid_sim%v)
         this%budget_0(:,:,:,19) = this%budget_0(:,:,:,19) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) * (this%pre_budget%igrid_sim%wC)
         
         ! STEP 5: SGS stresses (also viscous stress if finite reynolds number is being used)
         call this%pre_budget%igrid_sim%sgsmodel%populate_tauij_E_to_C()
         call this%prim_budget%igrid_sim%sgsmodel%populate_tauij_E_to_C()
         this%budget_0(:,:,:,20:25) = this%budget_0(:,:,:,20:25) + this%prim_budget%igrid_sim%tauSGS_ij - this%pre_budget%igrid_sim%tauSGS_ij 
 
         ! STEP 6: Compute mean T and heat fluxes
         if (this%isStratified) then
            this%budget_0(:,:,:,26) = this%budget_0(:,:,:,26) + this%prim_budget%igrid_sim%T - this%pre_budget%igrid_sim%T
            this%budget_0(:,:,:,27) = this%budget_0(:,:,:,27) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) &
                                        *(this%prim_budget%igrid_sim%T - this%pre_budget%igrid_sim%T)
            this%budget_0(:,:,:,28) = this%budget_0(:,:,:,28) + (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v) &
                                        *(this%prim_budget%igrid_sim%T - this%pre_budget%igrid_sim%T)
            this%budget_0(:,:,:,29) = this%budget_0(:,:,:,29) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) &
                                        *(this%prim_budget%igrid_sim%T - this%pre_budget%igrid_sim%T)
            this%budget_0(:,:,:,30) = this%budget_0(:,:,:,30) + (this%prim_budget%igrid_sim%T - this%pre_budget%igrid_sim%T) &
                                        *(this%prim_budget%igrid_sim%T - this%pre_budget%igrid_sim%T)
        end if 
 
     end subroutine
 
 
 
     ! ---------------------- Budget 1 ------------------------
     subroutine AssembleBudget1(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind) :: tmp1, tmp2, tmp3, tmp4
 
         ! STEP 1: Get terms from u-equation
         ! Total advection in deficit equation
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%uc - this%pre_budget%uc,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,1) = this%budget_1(:,:,:,1) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
 
         ! Fluctuating advection of deficit by base flow
         call this%ddx_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%u
 
         call this%ddy_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%v
 
         call this%ddz_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%wC
 
         ! Fluctuating advection of deficit by deficit flow
         call this%ddx_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
 
         call this%ddy_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
 
         call this%ddz_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
         
         ! Fluctuating advection of base by deficit flow                           
         call this%ddx_R2R(this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
 
         call this%ddy_R2R(this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
 
         call this%ddz_R2R(this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
 
         ! Pressure gradient dP/dx
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%px - this%pre_budget%px,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,8) = this%budget_1(:,:,:,8) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
         
         ! SGS stresses
         call this%prim_budget%igrid_sim%spectC%ifft((this%prim_budget%usgs - this%pre_budget%usgs),this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,9) = this%budget_1(:,:,:,9) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
         
         ! Acuator Disk
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%uturb,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,10) = this%budget_1(:,:,:,10) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
 
         ! STEP 2: Get terms from v-equation 
         ! Total advection in deficit equation
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%vc - this%pre_budget%vc,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,11) = this%budget_1(:,:,:,11) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
 
         ! Fluctuating advection of deficit by base flow
         call this%ddx_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%u
 
         call this%ddy_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%v
 
         call this%ddz_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%wC
 
         ! Fluctuating advection of deficit by deficit flow
         call this%ddx_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
 
         call this%ddy_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
 
         call this%ddz_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
                             
         ! Fluctuating advection of base by deficit flow
         call this%ddx_R2R(this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
 
         call this%ddy_R2R(this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
 
         call this%ddz_R2R(this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
 
         ! Pressure gradient dP/dy
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%py - this%pre_budget%py,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,18) = this%budget_1(:,:,:,18) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        
         ! SGS stresses
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%vsgs - this%pre_budget%vsgs,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,19) = this%budget_1(:,:,:,19) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
         
         ! Actuator Disk
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%vturb,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,20) = this%budget_1(:,:,:,20) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
         
         ! STEP 2: Get terms from w-equation
         ! Total advection in deficit equation
         call this%prim_budget%igrid_sim%spectE%ifft(this%prim_budget%wc - this%pre_budget%wc,this%prim_budget%igrid_sim%rbuffxE(:,:,:,1))
         call this%interp_Edge2Cell(this%prim_budget%igrid_sim%rbuffxE(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,21) = this%budget_1(:,:,:,21) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
 
         ! Fluctuating advection of deficit by base flow
         call this%ddx_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%u
 
         call this%ddy_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%v
 
         call this%ddz_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%igrid_sim%wC
 
         ! Fluctuating advection of deficit by deficit flow
         call this%ddx_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
 
         call this%ddy_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
 
         call this%ddz_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
          
         ! Fluctuating advection of base by deficit flow
         call this%ddx_R2R(this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
 
         call this%ddy_R2R(this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)
 
         call this%ddz_R2R(this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                 *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
 
         ! Pressure gradient dP/dz
         call this%prim_budget%igrid_sim%spectE%ifft(this%prim_budget%pz - this%pre_budget%pz,this%prim_budget%igrid_sim%rbuffxE(:,:,:,1))
         call this%interp_Edge2Cell(this%prim_budget%igrid_sim%rbuffxE(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,28) = this%budget_1(:,:,:,28) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
         
         ! SGS stresses
         call this%prim_budget%igrid_sim%spectE%ifft(this%prim_budget%wsgs - this%pre_budget%wsgs,this%prim_budget%igrid_sim%rbuffxE(:,:,:,1))
         call this%interp_Edge2Cell(this%prim_budget%igrid_sim%rbuffxE(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,29) = this%budget_1(:,:,:,29) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
 
         ! Buoyancy
         call this%prim_budget%igrid_sim%spectE%ifft(this%prim_budget%wb - this%pre_budget%wb,this%prim_budget%igrid_sim%rbuffxE(:,:,:,1))
         call this%interp_Edge2Cell(this%prim_budget%igrid_sim%rbuffxE(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
         this%budget_1(:,:,:,30) = this%budget_1(:,:,:,30) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        
         if (this%useCoriolis) then
             ! Get the geostrophic forcing 
             call this%pre_budget%igrid_sim%get_geostrophic_forcing(tmp1, tmp2)         ! Forcing in x and y directions respectively 
             call this%prim_budget%igrid_sim%get_geostrophic_forcing(tmp3, tmp4)         ! Forcing in x and y directions respectively 
             ! Coriolis term, X       
             call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%ucor,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
             call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%ucor,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
             this%budget_1(:,:,:,31) = this%budget_1(:,:,:,31) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - tmp3 &
                                          - (this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) - tmp1) ! Remove the geostrophic forcing term
             ! Geostrophic term, X
             this%budget_1(:,:,:,32) = this%budget_1(:,:,:,32) + tmp3 - tmp1
             ! Coriolis term, Y       
             call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%vcor,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
             call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%vcor,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
             this%budget_1(:,:,:,33) = this%budget_1(:,:,:,33) + this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - tmp4 & 
                                         - (this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) - tmp2) ! Remove the geostrophic forcing term
             ! Geostrophic term, Y
             this%budget_1(:,:,:,34) = this%budget_1(:,:,:,34) + tmp4 - tmp2
         end if 
         
     end subroutine
 
     subroutine DumpBudget1(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
         integer :: idx 
 
         ! Step 1: Get the average from sum
         this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
         this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)
         this%pre_budget%budget_0 = this%pre_budget%budget_0/(real(this%counter,rkind) + 1.d-18) ! this assumes the counter is the same for both
         
         !!!! Advection terms !!!!
         ! X - eqn
         ! mean advection of deficit by base flow
         call this%ddx_R2R(this%budget_0(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,2) = -this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%budget_0(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,2) = this%budget_1(:,:,:,2) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%budget_0(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,2) = this%budget_1(:,:,:,2) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)
 
         ! fluctuating advection of deficit by base flow
         this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) - this%budget_1(:,:,:,2)
 
         ! mean advection of deficit by deficit flow
         call this%ddx_R2R(this%budget_0(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,3) = -this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%budget_0(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,3) = this%budget_1(:,:,:,3) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%budget_0(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,3) = this%budget_1(:,:,:,3) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,3)
 
         ! fluctuating advection of deficit by deficit flow
         this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) - this%budget_1(:,:,:,3)
 
         ! mean advection of base by deficit flow
         call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,4) = -this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,4) = this%budget_1(:,:,:,4) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,4) = this%budget_1(:,:,:,4) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,3)
 
         ! fluctuating advection of base by deficit flow
         this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) - this%budget_1(:,:,:,4)
 
         ! Y - eqn
         ! mean advection of deficit by base flow
         call this%ddx_R2R(this%budget_0(:,:,:,2), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,12) = -this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%budget_0(:,:,:,2), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,12) = this%budget_1(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%budget_0(:,:,:,2), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,12) = this%budget_1(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)
 
         ! fluctuating advection of deficit by base flow
         this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) - this%budget_1(:,:,:,12)
 
         ! mean advection of deficit by deficit flow
         call this%ddx_R2R(this%budget_0(:,:,:,2), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,13) = -this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%budget_0(:,:,:,2), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,13) = this%budget_1(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%budget_0(:,:,:,2), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,13) = this%budget_1(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,3)
 
         ! fluctuating advection of deficit by deficit flow
         this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) - this%budget_1(:,:,:,13)
 
         ! mean advection of base by deficit flow
         call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,14) = -this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,14) = this%budget_1(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,14) = this%budget_1(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,3)
 
         ! fluctuating advection of deficit by deficit flow
         this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) - this%budget_1(:,:,:,14)
 
         ! Z - eqn
         ! mean advection of deficit by base flow
         call this%ddx_R2R(this%budget_0(:,:,:,3), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,22) = -this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%budget_0(:,:,:,3), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,22) = this%budget_1(:,:,:,22) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%budget_0(:,:,:,3), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,22) = this%budget_1(:,:,:,22) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)   
     
         ! fluctuating advection of deficit by base flow
         this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) - this%budget_1(:,:,:,22)
 
         ! mean advection of deficit by deficit flow
         call this%ddx_R2R(this%budget_0(:,:,:,3), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,23) = -this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%budget_0(:,:,:,3), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,23) = this%budget_1(:,:,:,23) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%budget_0(:,:,:,3), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,23) = this%budget_1(:,:,:,23) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,3)
 
         ! fluctuating advection of deficit by deficit flow
         this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) - this%budget_1(:,:,:,23)
 
         ! mean advection of base by deficit flow
         call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,24) = -this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,1)
 
         call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,24) = this%budget_1(:,:,:,24) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,2)
 
         call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
         this%budget_1(:,:,:,24) = this%budget_1(:,:,:,24) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)*this%budget_0(:,:,:,3)
 
         ! fluctuating advection of base by deficit flow
         this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) - this%budget_1(:,:,:,24)
 
         ! Step 2: Dump the full budget 
         do idx = 1,size(this%budget_1,4)
             call this%dump_budget_field(this%budget_1(:,:,:,idx),idx,1)
         end do 
         
         ! Step 3: Go back to summing instead of averaging
         this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) + this%budget_1(:,:,:,2)
         this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) + this%budget_1(:,:,:,3)
         this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) + this%budget_1(:,:,:,4)
         this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) + this%budget_1(:,:,:,12)
         this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) + this%budget_1(:,:,:,13)
         this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) + this%budget_1(:,:,:,14)
         this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) + this%budget_1(:,:,:,22)
         this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) + this%budget_1(:,:,:,23)
         this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) + this%budget_1(:,:,:,24)
 
         this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)
         this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)
         this%pre_budget%budget_0 = this%pre_budget%budget_0*(real(this%counter,rkind) + 1.d-18) ! this assumes that the counter is the same for both and it should be
 
     end subroutine 
 
     ! ---------------------- Budget 2 ------------------------
     subroutine AssembleBudget2(this)
         class(budgets_time_avg_deficit), intent(inout), target :: this
         real(rkind), dimension(:,:,:), pointer :: delta_Umn, delta_Vmn, delta_Wmn, base_Umn, base_Vmn, base_Wmn 
         real(rkind), dimension(:,:,:), pointer :: delta_R11, delta_R12, delta_R13, delta_R22, delta_R23, delta_R33
         real(rkind), dimension(:,:,:), pointer :: delta_base_R11, delta_base_R12, delta_base_R13, delta_base_R22 
         real(rkind), dimension(:,:,:), pointer :: delta_base_R23, delta_base_R33, delta_base_R21, delta_base_R31, delta_base_R32
         real(rkind), dimension(:,:,:), pointer :: delta_tau11, delta_tau12, delta_tau13, delta_tau22, delta_tau23, delta_tau33
         real(rkind), dimension(:,:,:), pointer :: buff, buff2
 
         if (this%counter > 0) then
 
             ! Calculate terms
             ! Get the average from sum
             this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
             this%pre_budget%budget_0 = this%pre_budget%budget_0/(real(this%counter,rkind) + 1.d-18)
        
             ! Step 2: Get the <Rij> from <ui uj>
             this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
             this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
             this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
             this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
             this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
             this%budget_0(:,:,:,10) = this%budget_0(:,:,:,10) - this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
     
             ! Step 2: Get the mixed <Rij> from mixed <ui uj>
             this%budget_0(:,:,:,11)  = this%budget_0(:,:,:,11)  - this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1) ! R11
             this%budget_0(:,:,:,12)  = this%budget_0(:,:,:,12)  - this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2) ! R12
             this%budget_0(:,:,:,13)  = this%budget_0(:,:,:,13)  - this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1) ! R21
             this%budget_0(:,:,:,14)  = this%budget_0(:,:,:,14)  - this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3) ! R13
             this%budget_0(:,:,:,15)  = this%budget_0(:,:,:,15)  - this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1) ! R31
             this%budget_0(:,:,:,16)  = this%budget_0(:,:,:,16)  - this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2) ! R22
             this%budget_0(:,:,:,17)  = this%budget_0(:,:,:,17)  - this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3) ! R23
             this%budget_0(:,:,:,18)  = this%budget_0(:,:,:,18)  - this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2) ! R32
             this%budget_0(:,:,:,19)  = this%budget_0(:,:,:,19)  - this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3) ! R33
             
             ! Budget 1 stress divergence terms
             this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)

             this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) - this%budget_1(:,:,:,2)
             this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) - this%budget_1(:,:,:,3)
             this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) - this%budget_1(:,:,:,4)
             this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) - this%budget_1(:,:,:,12)
             this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) - this%budget_1(:,:,:,13)
             this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) - this%budget_1(:,:,:,14)
             this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) - this%budget_1(:,:,:,22)
             this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) - this%budget_1(:,:,:,23)
             this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) - this%budget_1(:,:,:,24)
 
             delta_Umn => this%budget_0(:,:,:,1);    delta_Vmn => this%budget_0(:,:,:,2);      delta_Wmn => this%budget_0(:,:,:,3);
             base_Umn => this%pre_budget%budget_0(:,:,:,1);    base_Vmn => this%pre_budget%budget_0(:,:,:,2);      base_Wmn => this%pre_budget%budget_0(:,:,:,3);
             delta_R11 => this%budget_0(:,:,:,5);    delta_R12 => this%budget_0(:,:,:,6);      delta_R13 => this%budget_0(:,:,:,7)
             delta_R22 => this%budget_0(:,:,:,8);    delta_R23 => this%budget_0(:,:,:,9);      delta_R33 => this%budget_0(:,:,:,10)
             delta_base_R11 => this%budget_0(:,:,:,11);    delta_base_R12 => this%budget_0(:,:,:,12);      delta_base_R21 => this%budget_0(:,:,:,13)
             delta_base_R13 => this%budget_0(:,:,:,14);    delta_base_R31 => this%budget_0(:,:,:,15);      delta_base_R22 => this%budget_0(:,:,:,16)    
             delta_base_R23 => this%budget_0(:,:,:,17);    delta_base_R32 => this%budget_0(:,:,:,18);      delta_base_R33 => this%budget_0(:,:,:,19)
             delta_tau11 => this%budget_0(:,:,:,20);   delta_tau12 => this%budget_0(:,:,:,21); delta_tau13 => this%budget_0(:,:,:,22) 
             delta_tau22 => this%budget_0(:,:,:,23);   delta_tau23 => this%budget_0(:,:,:,24); delta_tau33 => this%budget_0(:,:,:,25)
             buff => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1); buff2 => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
 
             ! 1:  Loss to Resolved TKE components, dissipation,
             ! and mean advection components
             call this%ddx_R2R(delta_Umn, buff); 
             this%budget_2(:,:,:,11) = delta_R11*buff
             this%budget_2(:,:,:,12) = delta_base_R11*buff
             this%budget_2(:,:,:,13) = delta_base_R11*buff
 
             this%budget_2(:,:,:,6) = delta_tau11*buff
 
             this%budget_2(:,:,:,17) = - base_Umn*delta_Umn*buff
             this%budget_2(:,:,:,18) = - delta_Umn*delta_Umn*buff
             
             call this%ddx_R2R(delta_Vmn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R12*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R21*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R12*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau12*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Umn*delta_Vmn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Umn*delta_Vmn*buff
 
             call this%ddx_R2R(delta_Wmn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R13*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R31*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R13*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau13*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Umn*delta_Wmn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Umn*delta_Wmn*buff
 
             call this%ddy_R2R(delta_Umn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R12*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R12*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R21*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau12*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Vmn*delta_Umn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Vmn*delta_Umn*buff
 
             call this%ddy_R2R(delta_Vmn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R22*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R22*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R22*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau22*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Vmn*delta_Vmn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Vmn*delta_Vmn*buff
 
             call this%ddy_R2R(delta_Wmn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R23*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R32*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R23*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau23*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Vmn*delta_Wmn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Vmn*delta_Wmn*buff
 
             call this%ddz_R2R(delta_Umn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R13*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R13*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R31*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau13*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Wmn*delta_Umn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Wmn*delta_Umn*buff
 
             call this%ddz_R2R(delta_Vmn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R23*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R23*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R32*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau23*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Wmn*delta_Vmn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Wmn*delta_Vmn*buff
 
             call this%ddz_R2R(delta_Wmn, buff);
             this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + delta_R33*buff
             this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + delta_base_R33*buff
             this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + delta_base_R33*buff
 
             this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + delta_tau33*buff
 
             this%budget_2(:,:,:,17) = this%budget_2(:,:,:,17) - base_Wmn*delta_Wmn*buff
             this%budget_2(:,:,:,18) = this%budget_2(:,:,:,18) - delta_Wmn*delta_Wmn*buff
 
             this%budget_2(:,:,:,1) = this%budget_2(:,:,:,11) + this%budget_2(:,:,:,12) + this%budget_2(:,:,:,13)  
 
             ! 2: Advection - total
             buff2 = half*(delta_Umn*delta_Umn + delta_Vmn*delta_Vmn + delta_Wmn*delta_Wmn)
             call this%ddx_R2R(buff2,buff); this%budget_2(:,:,:,2) = -(base_Umn + delta_Umn)*buff
             call this%ddy_R2R(buff2,buff); this%budget_2(:,:,:,2) = this%budget_2(:,:,:,2) - (base_Vmn + delta_Vmn)*buff
             call this%ddz_R2R(buff2,buff); this%budget_2(:,:,:,2) = this%budget_2(:,:,:,2) - (base_Wmn + delta_Wmn)*buff
 
             ! 4:  Pressure transport        (C)
             this%budget_2(:,:,:,4) = delta_Umn*this%budget_1(:,:,:,8) + delta_Vmn*this%budget_1(:,:,:,18) + delta_Wmn*this%budget_1(:,:,:,28)
 
             ! 5:  SGS + viscous transport   (D+F)
             this%budget_2(:,:,:,5) = delta_Umn*this%budget_1(:,:,:,9) + delta_Vmn*this%budget_1(:,:,:,19) + delta_Wmn*this%budget_1(:,:,:,29)
             this%budget_2(:,:,:,5) = this%budget_2(:,:,:,5) - this%budget_2(:,:,:,6)
 
             ! 7:  Actuator disk sink        (J)
             !!!!!!!!!!!!!!!!!!!
             ! this assumes that the turbine forcing is only in the x direction
             ! according to the u velocity
             !!!!!!!!!!!!!!!!!!!
             this%budget_2(:,:,:,7) = delta_Umn*this%budget_1(:,:,:,10)
 
             ! 8: Coriolis terms
             if (this%useCoriolis) then
                 ! Geostrophic forcing term
                 this%budget_2(:,:,:,8) = delta_Umn*this%budget_1(:,:,:,32) &
                                + delta_Vmn*this%budget_1(:,:,:,34)
             
                 ! Coriolis forcing term (should be 0)
                 this%budget_2(:,:,:,9) = delta_Umn*this%budget_1(:,:,:,31) &
                                + delta_Vmn*this%budget_1(:,:,:,33)
             end if
 
             ! 10: Buoyancy
             this%budget_2(:,:,:,10) = delta_Wmn*this%budget_1(:,:,:,30)
    
            ! 3:  Turbulent transport
            this%budget_2(:,:,:,14) = this%budget_0(:,:,:,1)*this%budget_1(:,:,:,6) + this%budget_0(:,:,:,2)*this%budget_1(:,:,:,16) &
                                        + this%budget_0(:,:,:,3)*this%budget_1(:,:,:,26) - this%budget_2(:,:,:,11)  
            this%budget_2(:,:,:,15) = this%budget_0(:,:,:,1)*this%budget_1(:,:,:,5) + this%budget_0(:,:,:,2)*this%budget_1(:,:,:,15) &
                                        + this%budget_0(:,:,:,3)*this%budget_1(:,:,:,25) - this%budget_2(:,:,:,12)  
            this%budget_2(:,:,:,16) = this%budget_0(:,:,:,1)*this%budget_1(:,:,:,7) + this%budget_0(:,:,:,2)*this%budget_1(:,:,:,17) &
                                        + this%budget_0(:,:,:,3)*this%budget_1(:,:,:,27) - this%budget_2(:,:,:,13) 
    
            this%budget_2(:,:,:,3) = this%budget_2(:,:,:,14) + this%budget_2(:,:,:,15) + this%budget_2(:,:,:,16)      
    
             ! Mixed mean advection
             ! Mean advection components
             call this%ddx_R2R(base_Umn, buff); 
             this%budget_2(:,:,:,19) = - delta_Umn*delta_Umn*buff
             
             call this%ddx_R2R(base_Vmn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Umn*delta_Vmn*buff
 
             call this%ddx_R2R(base_Wmn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Umn*delta_Wmn*buff
 
             call this%ddy_R2R(base_Umn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Vmn*delta_Umn*buff
 
             call this%ddy_R2R(base_Vmn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Vmn*delta_Vmn*buff
 
             call this%ddy_R2R(base_Wmn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Vmn*delta_Wmn*buff
 
             call this%ddz_R2R(base_Umn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Wmn*delta_Umn*buff
 
             call this%ddz_R2R(base_Vmn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Wmn*delta_Vmn*buff
 
             call this%ddz_R2R(base_Wmn, buff);
             this%budget_2(:,:,:,19) = this%budget_2(:,:,:,19) - delta_Wmn*delta_Wmn*buff            
 
             
             nullify(delta_Umn,delta_Vmn,delta_Wmn,base_Umn,base_Vmn,base_Wmn, &
                     delta_R11,delta_R12,delta_R13, &
                     delta_R22,delta_R23,delta_R33, &
                     delta_base_R11, delta_base_R12, delta_base_R21, &
                     delta_base_R22, delta_base_R13, delta_base_R31, &
                     delta_base_R23, delta_base_R32, &
                     delta_base_R33, delta_tau11, delta_tau12, &
                     delta_tau22, delta_tau23, delta_tau13, delta_tau33, buff,buff2)
 
             ! Go back to sum
            this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) + this%budget_1(:,:,:,2)
            this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) + this%budget_1(:,:,:,3)
            this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) + this%budget_1(:,:,:,4)
            this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) + this%budget_1(:,:,:,12)
            this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) + this%budget_1(:,:,:,13)
            this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) + this%budget_1(:,:,:,14)
            this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) + this%budget_1(:,:,:,22)
            this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) + this%budget_1(:,:,:,23)
            this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) + this%budget_1(:,:,:,24)
             this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)   
 
             this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
             this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
             this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
             this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
             this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
             this%budget_0(:,:,:,10) = this%budget_0(:,:,:,10) + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
             
             this%budget_0(:,:,:,11)  = this%budget_0(:,:,:,11)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1) ! R11
             this%budget_0(:,:,:,12)  = this%budget_0(:,:,:,12)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2) ! R12
             this%budget_0(:,:,:,13)  = this%budget_0(:,:,:,13)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1) ! R12
             this%budget_0(:,:,:,14)  = this%budget_0(:,:,:,14)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3) ! R13
             this%budget_0(:,:,:,15)  = this%budget_0(:,:,:,15)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1) ! R31
             this%budget_0(:,:,:,16)  = this%budget_0(:,:,:,16)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2) ! R22
             this%budget_0(:,:,:,17)  = this%budget_0(:,:,:,17)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3) ! R23
             this%budget_0(:,:,:,18)  = this%budget_0(:,:,:,18)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2) ! R32
             this%budget_0(:,:,:,19)  = this%budget_0(:,:,:,19)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3) ! R33
 
             this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)
             this%pre_budget%budget_0 = this%pre_budget%budget_0*(real(this%counter,rkind) + 1.d-18)
 
         end if 
 
     end subroutine 
     
     subroutine DumpBudget2(this)
        class(budgets_time_avg_deficit), intent(inout) :: this
        integer :: idx    

        ! Dump the full budget 
        do idx = 1,size(this%budget_2,4)
            call this%dump_budget_field(this%budget_2(:,:,:,idx),idx,2)
        end do 
 
    end subroutine 
     
    ! ---------------------- Budget 3 ------------------------
    subroutine AssembleBudget3(this)
        class(budgets_time_avg_deficit), intent(inout) :: this

        call this%pre_budget%igrid_sim%sgsmodel%populate_tauij_E_to_C()
        call this%prim_budget%igrid_sim%sgsmodel%populate_tauij_E_to_C()

        !!! 3. turbulent transport        
        ! < delta ui' delta uj' d delta ui'/dxj >
        call this%ddx_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddy_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddz_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddx_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddy_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddz_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddx_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        call this%ddy_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        call this%ddz_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        ! < delta ui' base uj' d delta ui'/dxj >
        call this%ddx_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddy_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddz_R2R(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddx_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddy_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddz_R2R(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddx_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        call this%ddy_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        call this%ddz_R2R(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)


        ! < delta ui' delta uj' d base ui'/dxj >
        call this%ddx_R2R(this%pre_budget%igrid_sim%u, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        ! base advection term in x
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) 

        call this%ddy_R2R(this%pre_budget%igrid_sim%u, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        ! base advection term in y
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddz_R2R(this%pre_budget%igrid_sim%u, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)
                      
        ! base advection term in z
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)

        call this%ddx_R2R(this%pre_budget%igrid_sim%v, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        ! base advection term in x
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddy_R2R(this%pre_budget%igrid_sim%v, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        ! base advection term in y
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddz_R2R(this%pre_budget%igrid_sim%v, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        ! base advection term in z
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)

        call this%ddx_R2R(this%pre_budget%igrid_sim%wC, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        ! base advection term in x
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        call this%ddy_R2R(this%pre_budget%igrid_sim%wC, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        ! base advection term in y
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        call this%ddz_R2R(this%pre_budget%igrid_sim%wC, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)*(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)
                
        ! base advection term in z
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1) &
                                *(this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        ! 4. Pressure transport      
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%px,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%px,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) &
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%py,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%py,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v) &
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))

        call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%pz,this%pre_budget%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%pre_budget%igrid_sim%rbuffxE(:,:,:,1), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectE%ifft(this%prim_budget%pz,this%prim_budget%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%prim_budget%igrid_sim%rbuffxE(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) &
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))

        ! 5. SGS + viscous transport     (E+F)
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%usgs,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%usgs,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) &
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%vsgs,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%vsgs,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v) &
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))

        call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%wsgs,this%pre_budget%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%pre_budget%igrid_sim%rbuffxE(:,:,:,1), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectE%ifft(this%prim_budget%wsgs,this%prim_budget%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%prim_budget%igrid_sim%rbuffxE(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))

        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5)  + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) &
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        

        ! 6. SGS + viscous dissipation 
        call this%ddx_R2R(this%pre_budget%igrid_sim%u, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddx_R2R(this%prim_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    *(this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,1) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,1))

        call this%ddx_R2R(this%pre_budget%igrid_sim%v, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddx_R2R(this%prim_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    * (this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,2) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,2))

        call this%ddx_R2R(this%pre_budget%igrid_sim%wC, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddx_R2R(this%prim_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    * (this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,3) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,3))

        call this%ddy_R2R(this%pre_budget%igrid_sim%u, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddy_R2R(this%prim_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    *(this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,2) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,2))

        call this%ddy_R2R(this%pre_budget%igrid_sim%v, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddy_R2R(this%prim_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    * (this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,4) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,4))

        call this%ddy_R2R(this%pre_budget%igrid_sim%wC, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddy_R2R(this%prim_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    * (this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,5) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,5))

        call this%ddz_R2R(this%pre_budget%igrid_sim%u, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddz_R2R(this%prim_budget%igrid_sim%u, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    *(this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,3) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,3))

        call this%ddz_R2R(this%pre_budget%igrid_sim%v, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddz_R2R(this%prim_budget%igrid_sim%v, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    * (this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,5) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,5))

        call this%ddz_R2R(this%pre_budget%igrid_sim%wC, this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        call this%ddz_R2R(this%prim_budget%igrid_sim%wC, this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)) &
                                    * (this%prim_budget%igrid_sim%tauSGS_ij(:,:,:,6) - this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,6))

        ! 7. Actuator disk/Turbine sink  (J)
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%vturb,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%vturb,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + (this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v) & 
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))

        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%uturb,this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%uturb,this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + (this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u) & 
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))

 
        ! 8. Buoyancy transfer: 
        call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%wb,this%pre_budget%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%pre_budget%igrid_sim%rbuffxE(:,:,:,1), this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectE%ifft(this%prim_budget%wb,this%prim_budget%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%prim_budget%igrid_sim%rbuffxE(:,:,:,1), this%prim_budget%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC) &
                                    * (this%prim_budget%igrid_sim%rbuffxC(:,:,:,1) - this%pre_budget%igrid_sim%rbuffxC(:,:,:,1))

    end subroutine 
    
    subroutine DumpBudget3(this)
        class(budgets_time_avg_deficit), intent(inout), target :: this
        integer :: idx
        real(rkind), dimension(:,:,:), pointer :: delta_Umn, delta_Vmn, delta_Wmn, base_Umn, base_Vmn, base_Wmn 
        real(rkind), dimension(:,:,:), pointer :: delta_R11, delta_R12, delta_R13, delta_R22, delta_R23, delta_R33
        real(rkind), dimension(:,:,:), pointer :: buff, buff2

        ! Calculate terms
        ! Get the average from sum
        this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
        this%pre_budget%budget_0 = this%pre_budget%budget_0/(real(this%counter,rkind) + 1.d-18)
    
        ! Step 2: Get the <Rij> from <ui uj>
        this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
        this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
        this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
        this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
        this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
        this%budget_0(:,:,:,10) = this%budget_0(:,:,:,10) - this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33

        
        ! Budget 1 stress divergence terms
        this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)

        this%budget_3 = this%budget_3/(real(this%counter,rkind) + 1.d-18)

        delta_Umn => this%budget_0(:,:,:,1);    delta_Vmn => this%budget_0(:,:,:,2);      delta_Wmn => this%budget_0(:,:,:,3);
        base_Umn => this%pre_budget%budget_0(:,:,:,1);    base_Vmn => this%pre_budget%budget_0(:,:,:,2);      base_Wmn => this%pre_budget%budget_0(:,:,:,3);
        delta_R11 => this%budget_0(:,:,:,5);    delta_R12 => this%budget_0(:,:,:,6);      delta_R13 => this%budget_0(:,:,:,7)
        delta_R22 => this%budget_0(:,:,:,8);    delta_R23 => this%budget_0(:,:,:,9);      delta_R33 => this%budget_0(:,:,:,10);
        buff => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1); buff2 => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)

        !!! 1. TKE production   
        ! TKE production - <delta ui' delta uj'> d delta ui / dxj
        this%budget_3(:,:,:,9) = - this%budget_2(:,:,:,11)
   
        ! TKE production - <delta ui' base uj'> d delta ui / dxj
        this%budget_3(:,:,:,10) = - this%budget_2(:,:,:,12)

        ! Total TKE production           
        ! this does not include production from interaction of base flow gradient with the deficit stresses
        this%budget_3(:,:,:,1) = this%budget_3(:,:,:,9) + this%budget_3(:,:,:,10)

        ! 2. advective transport     
        buff2 = half*(delta_R11 + delta_R22 + delta_R33)
        call this%ddx_R2R(buff2,buff); this%budget_3(:,:,:,16) = -delta_Umn*buff
        this%budget_3(:,:,:,15) = -base_Umn*buff
        call this%ddy_R2R(buff2,buff); this%budget_3(:,:,:,16) = this%budget_3(:,:,:,16) - delta_Vmn*buff
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) - base_Vmn*buff
        call this%ddz_R2R(buff2,buff); this%budget_3(:,:,:,16) = this%budget_3(:,:,:,16) - delta_Wmn*buff
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) - base_Wmn*buff

        ! total advective transport
        this%budget_3(:,:,:,2) = this%budget_3(:,:,:,15) + this%budget_3(:,:,:,16)

        ! 17. Base advective term by delta
        buff2 = base_Umn
        call this%ddx_R2R(buff2,buff)
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + buff*delta_Umn

        ! 11. TKE production - <delta uj' delta ui'> d base ui / dxj
        this%budget_3(:,:,:,11) = - buff*delta_R11

        call this%ddy_R2R(buff2,buff);
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + buff*delta_Umn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R12

        call this%ddz_R2R(buff2,buff);
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) +buff*delta_Umn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R13

        buff2 = base_Vmn
        call this%ddx_R2R(buff2,buff)
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + buff*delta_Vmn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R12

        call this%ddy_R2R(buff2,buff);
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + buff*delta_Vmn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R22

        call this%ddz_R2R(buff2,buff);
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) + buff*delta_Vmn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R23

        buff2 = base_Wmn
        call this%ddx_R2R(buff2,buff)
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + buff*delta_Wmn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R13

        call this%ddy_R2R(buff2,buff);
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + buff*delta_Wmn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R23

        call this%ddz_R2R(buff2,buff);
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) + buff*delta_Wmn
        this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) - buff*delta_R33

        this%budget_3(:,:,:,17) = delta_Umn*this%budget_3(:,:,:,18) + delta_Vmn*this%budget_3(:,:,:,19) + delta_Wmn*this%budget_3(:,:,:,20)
                                    
        !!! 3. turbulent transport 
        ! < delta ui' delta uj' d delta ui'/dxj >        
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) - this%budget_3(:,:,:,16) - this%budget_2(:,:,:,18) - this%budget_2(:,:,:,14)

        ! < delta ui' base uj' d delta ui'/dxj >        
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) - this%budget_3(:,:,:,15) - this%budget_2(:,:,:,17) - this%budget_2(:,:,:,15)

        ! < delta ui' delta uj' d base ui'/dxj >        
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) - this%budget_3(:,:,:,17) - this%budget_2(:,:,:,19) - this%budget_2(:,:,:,16) &
             - this%budget_2(:,:,:,13) -  this%budget_3(:,:,:,11)

        ! total (choosing to only include true transport terms for now)
        this%budget_3(:,:,:,3) = this%budget_3(:,:,:,12) + this%budget_3(:,:,:,13)

        ! 4. Pressure transport         
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) - this%budget_2(:,:,:,4)

        ! 5. SGS + viscous transport     
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - this%budget_3(:,:,:,6) - this%budget_2(:,:,:,5) 

        !- this%budget_2(:,:,:,6) kktodo

        ! 6. SGS + viscous dissipation  
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) - this%budget_2(:,:,:,6)
        
        ! 7. Actuator disk/Turbine sink 
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) - delta_Umn*this%budget_1(:,:,:,10) - delta_Vmn*this%budget_1(:,:,:,20)

        ! 8. Buoyancy
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) - delta_Wmn*this%budget_1(:,:,:,30)

        ! Dump the full budget 
        do idx = 1,size(this%budget_3,4)
            call this%dump_budget_field(this%budget_3(:,:,:,idx),idx,3)
        end do 

        ! Revert arrays to the correct state for Assemble (Order is very
        ! important throughout this subroutine, particularly indices 5 and 6)       
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + this%budget_3(:,:,:,17) + this%budget_2(:,:,:,19) + this%budget_2(:,:,:,16) &
             + this%budget_2(:,:,:,13) + this%budget_3(:,:,:,11)
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + this%budget_3(:,:,:,15) + this%budget_2(:,:,:,17) + this%budget_2(:,:,:,15)
        this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) + this%budget_3(:,:,:,16) + this%budget_2(:,:,:,18) + this%budget_2(:,:,:,14)
        
        buff2 = base_Umn
        call this%ddx_R2R(buff2,buff)
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - buff*delta_Umn

        call this%ddy_R2R(buff2,buff);
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - buff*delta_Umn

        call this%ddz_R2R(buff2,buff);
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - buff*delta_Umn

        buff2 = base_Vmn
        call this%ddx_R2R(buff2,buff)
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - buff*delta_Vmn

        call this%ddy_R2R(buff2,buff);
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - buff*delta_Vmn

        call this%ddz_R2R(buff2,buff);
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - buff*delta_Vmn

        buff2 = base_Wmn
        call this%ddx_R2R(buff2,buff)
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - buff*delta_Wmn

        call this%ddy_R2R(buff2,buff);
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - buff*delta_Wmn

        call this%ddz_R2R(buff2,buff);
        this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - buff*delta_Wmn

        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + delta_Wmn*this%budget_1(:,:,:,30)
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + delta_Umn*this%budget_1(:,:,:,10) + delta_Vmn*this%budget_1(:,:,:,20)
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%budget_2(:,:,:,6)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + this%budget_3(:,:,:,6) + this%budget_2(:,:,:,5) !+ this%budget_2(:,:,:,6) kktodo
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + this%budget_2(:,:,:,4)

        nullify(delta_Umn,delta_Vmn,delta_Wmn,base_Umn,base_Vmn,base_Wmn, &
        delta_R11,delta_R12,delta_R13, &
        delta_R22,delta_R23,delta_R33, buff,buff2)

        ! Go back to sum
        this%budget_3 = this%budget_3*(real(this%counter,rkind) + 1.d-18)
        this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)

        this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
        this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
        this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
        this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
        this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
        this%budget_0(:,:,:,10) = this%budget_0(:,:,:,10) + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33

        this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)
        this%pre_budget%budget_0 = this%pre_budget%budget_0*(real(this%counter,rkind) + 1.d-18)

    end subroutine 

 

     ! ----------------------supporting subroutines ------------------------
     subroutine dump_budget_field(this, field, fieldID, BudgetID)
         use decomp_2d_io
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(in) :: field
         integer, intent(in) :: fieldID, BudgetID
         character(len=clen) :: fname, tempname 
 
         write(tempname,"(A3,I2.2,A15,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_deficit_budget",BudgetID,"_term",fieldID,"_t",this%prim_budget%igrid_sim%step,"_n",this%counter,".s3D"
         fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
 
         call decomp_2d_write_one(1,field,fname, this%prim_budget%igrid_sim%gpC)
 
     end subroutine 
 
     subroutine RestartBudget(this, dir, rid, tid, cid)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(:,:,:), pointer :: buff, buff2
         integer, intent(in) :: rid, cid, tid
         character(len=clen) :: dir
         integer :: idx
 
         buff => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
         buff2 => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
 
         this%counter = cid
 
         ! Budget 0: 
         do idx = 1,size(this%budget_0,4)
            !          if (allocated(this%budget_0)) deallocate(this%budget_0)
            call this%restart_budget_field(this%budget_0(:,:,:,idx), dir, rid, tid, cid, 0, idx)
         end do
          
         this%pre_budget%budget_0 = this%pre_budget%budget_0/(real(this%counter,rkind) + 1.d-18)

         ! Step 1: Go back to <ui uj> from <Rij>
         this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
         this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
         this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
         this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
         this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
         this%budget_0(:,:,:,10)  = this%budget_0(:,:,:,10)  + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
         
         this%budget_0(:,:,:,11)  = this%budget_0(:,:,:,11)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1) ! R11
         this%budget_0(:,:,:,12)  = this%budget_0(:,:,:,12)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2) ! R12
         this%budget_0(:,:,:,13)  = this%budget_0(:,:,:,13)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1) ! R12
         this%budget_0(:,:,:,14)  = this%budget_0(:,:,:,14)  + this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3) ! R13
         this%budget_0(:,:,:,15)  = this%budget_0(:,:,:,15)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1) ! R13
         this%budget_0(:,:,:,16)  = this%budget_0(:,:,:,16)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2) ! R22
         this%budget_0(:,:,:,17)  = this%budget_0(:,:,:,17)  + this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3) ! R23
         this%budget_0(:,:,:,18)  = this%budget_0(:,:,:,18)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2) ! R23
         this%budget_0(:,:,:,19)  = this%budget_0(:,:,:,19)  + this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3) ! R33
         
        ! STEP 4: Compute Delta heat fluxes
         if (this%isStratified) then
            this%budget_0(:,:,:,27) = this%budget_0(:,:,:,27) + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,28) = this%budget_0(:,:,:,28) + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,29) = this%budget_0(:,:,:,29) + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,26)
            this%budget_0(:,:,:,30) = this%budget_0(:,:,:,30) + this%budget_0(:,:,:,26)*this%budget_0(:,:,:,26)
        end if 

        ! Step 11: Go back to summing instead of averaging
        this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)
        this%pre_budget%budget_0 = this%pre_budget%budget_0*(real(this%counter,rkind) + 1.d-18)
 
         ! Budget 1: 
        if (this%budgetType>0) then
            do idx = 1,size(this%budget_1,4)
            !          if (allocated(this%budget_1)) deallocate(this%budget_1)
            call this%restart_budget_field(this%budget_1(:,:,:,idx), dir, rid, tid, cid, 1, idx)
            end do

            ! fluctuating advection of deficit by base flow
            this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) + this%budget_1(:,:,:,2)
            this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) + this%budget_1(:,:,:,12)
            this%budget_1(:,:,:,25) = this%budget_1(:,:,:,25) + this%budget_1(:,:,:,22)
    
            ! fluctuating advection of deficit by deficit flow
            this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) + this%budget_1(:,:,:,3)
            this%budget_1(:,:,:,16) = this%budget_1(:,:,:,16) + this%budget_1(:,:,:,13)
            this%budget_1(:,:,:,26) = this%budget_1(:,:,:,26) + this%budget_1(:,:,:,23)
    
            ! fluctuating advection of deficit by deficit flow
            this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) + this%budget_1(:,:,:,4)
            this%budget_1(:,:,:,17) = this%budget_1(:,:,:,17) + this%budget_1(:,:,:,14)
            this%budget_1(:,:,:,27) = this%budget_1(:,:,:,27) + this%budget_1(:,:,:,24)
    
            this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)

        end if 



        ! Budget 2
         if (this%budgetType>1) then
            do idx = 1,size(this%budget_2,4)
               call this%restart_budget_field(this%budget_2(:,:,:,idx), dir, rid, tid, cid, 2, idx)   
            end do
         end if

        ! Budget 3
         if (this%budgetType>2) then
            do idx = 1,size(this%budget_3,4)
                call this%restart_budget_field(this%budget_3(:,:,:,idx), dir, rid, tid, cid, 3, idx)   
             end do

            ! Get the average from sum
            this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
            this%pre_budget%budget_0 = this%pre_budget%budget_0/(real(this%counter,rkind) + 1.d-18)
            this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)

            !  Revert arrays to the correct state for Assemble (Order is very
            ! important throughout this subroutine, particularly indices 5 and 6)       
            this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + this%budget_3(:,:,:,17) + this%budget_2(:,:,:,19) + this%budget_2(:,:,:,16) &
                 + this%budget_2(:,:,:,13) + this%budget_3(:,:,:,11)
            this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + this%budget_3(:,:,:,15) + this%budget_2(:,:,:,17) + this%budget_2(:,:,:,15)
            this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) + this%budget_3(:,:,:,16) + this%budget_2(:,:,:,18) + this%budget_2(:,:,:,14)
            
            buff2 = this%pre_budget%budget_0(:,:,:,1)
            call this%ddx_R2R(buff2,buff)
            this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - buff*this%budget_0(:,:,:,1)
    
            call this%ddy_R2R(buff2,buff)
            this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - buff*this%budget_0(:,:,:,1)
    
            call this%ddz_R2R(buff2,buff)
            this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - buff*this%budget_0(:,:,:,1)
    
            buff2 = this%pre_budget%budget_0(:,:,:,2)
            call this%ddx_R2R(buff2,buff)
            this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - buff*this%budget_0(:,:,:,2)
    
            call this%ddy_R2R(buff2,buff)
            this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - buff*this%budget_0(:,:,:,2)
    
            call this%ddz_R2R(buff2,buff)
            this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - buff*this%budget_0(:,:,:,2)
    
            buff2 = this%pre_budget%budget_0(:,:,:,3)
            call this%ddx_R2R(buff2,buff)
            this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) - buff*this%budget_0(:,:,:,3)
    
            call this%ddy_R2R(buff2,buff)
            this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) - buff*this%budget_0(:,:,:,3)
    
            call this%ddz_R2R(buff2,buff)
            this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) - buff*this%budget_0(:,:,:,3)
    
            this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + this%budget_0(:,:,:,3)*this%budget_1(:,:,:,30)
            this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + this%budget_0(:,:,:,1)*this%budget_1(:,:,:,10) + this%budget_0(:,:,:,2)*this%budget_1(:,:,:,20)
            this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%budget_2(:,:,:,6)
            this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + this%budget_3(:,:,:,6) + this%budget_2(:,:,:,5)
            this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + this%budget_2(:,:,:,4)
            
            ! Budget 1 stress divergence terms
            this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)
            this%pre_budget%budget_0 = this%pre_budget%budget_0*(real(this%counter,rkind) + 1.d-18)
            this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)
            this%budget_3 = this%budget_3*(real(this%counter,rkind) + 1.d-18)
         end if

         nullify(buff, buff2)
         
     end subroutine     
 
     subroutine restart_budget_field(this, field, dir, runID, timeID, counterID, budgetID, fieldID)
         use decomp_2d_io
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: field
         integer, intent(in) :: runID, counterID, timeID, budgetID, fieldID
         character(len=clen) :: fname, tempname
         character(len=clen), intent(in) :: dir
 
         write(tempname,"(A3,I2.2,A15,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A4)") "Run",runID,"_deficit_budget",budgetID,"_term",fieldID,"_t",timeID,"_n",counterID,".s3D"
         fname = dir(:len_trim(dir))//"/"//trim(tempname)
         call decomp_2d_read_one(1,field,fname, this%prim_budget%igrid_sim%gpC)           
     end subroutine 
 
     subroutine ResetBudget(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
         
         this%counter = 0
         this%budget_0 = 0.d0 
         this%budget_1 = 0.d0 
         this%budget_2 = 0.d0 
         this%budget_3 = 0.d0 
         this%budget_4_11 = 0.d0 
         this%budget_4_22 = 0.d0 
         this%budget_4_33 = 0.d0 
         this%budget_4_13 = 0.d0 
         this%budget_4_23 = 0.d0 
         
     end subroutine 
     
     subroutine destroy(this)
         class(budgets_time_avg_deficit), intent(inout) :: this
 
         nullify(this%prim_budget%igrid_sim)
         if(this%do_budgets) then
      !       deallocate(this%uc, this%vc, this%wc, this%usgs, this%vsgs, this%wsgs, this%px, this%py, this%pz, this%uturb)  
             deallocate(this%budget_0, this%budget_1)
             deallocate(this%runningSum_sc)
         end if
         if(this%useWindTurbines) then
             deallocate(this%runningSum_sc_turb)
             deallocate(this%runningSum_turb)
         endif
 
     end subroutine 
 
     ! ----------------------private derivative operators ------------------------
     subroutine ddx_R2R(this, f, dfdx)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(in) :: f
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: dfdx
         
         call this%prim_budget%igrid_sim%spectC%fft(f,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%mtimes_ik1_ip(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdx)
     end subroutine 
 
     subroutine ddx_C2R(this, fhat, dfdx)
         class(budgets_time_avg_deficit), intent(inout) :: this
         complex(rkind), dimension(this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(1),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(2),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: dfdx
         
         call this%prim_budget%igrid_sim%spectC%mtimes_ik1_oop(fhat,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdx)
     end subroutine 
     
     subroutine ddy_R2R(this, f, dfdy)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(in) :: f
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: dfdy
         
         call this%prim_budget%igrid_sim%spectC%fft(f,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%mtimes_ik2_ip(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdy)
     end subroutine 
 
     subroutine ddy_C2R(this, fhat, dfdy)
         class(budgets_time_avg_deficit), intent(inout) :: this
         complex(rkind), dimension(this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(1),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(2),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: dfdy
         
         call this%prim_budget%igrid_sim%spectC%mtimes_ik2_oop(fhat,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdy)
     end subroutine 
     
     subroutine ddz_R2R(this, f, dfdz)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(in) :: f
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: dfdz
         
         !call transpose_x_to_y(f,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         !call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         !call this%prim_budget%igrid_sim%Pade6opZ%ddz_C2C(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,2),0,0)
         !call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzC(:,:,:,2),this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         !call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),dfdz,this%prim_budget%igrid_sim%gpC)
 
         call this%prim_budget%igrid_sim%spectC%fft(f,this%prim_budget%igrid_sim%cbuffyC(:,:,:,2))
         call this%ddz_C2R(this%prim_budget%igrid_sim%cbuffyC(:,:,:,2), dfdz)
 
     end subroutine 
     
     subroutine ddz_C2R(this, fhat, dfdz)
         class(budgets_time_avg_deficit), intent(inout) :: this
         complex(rkind), dimension(this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(1),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(2),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: dfdz
         
         call transpose_y_to_z(fhat,this%prim_budget%igrid_sim%cbuffzC(:,:,:,1),this%prim_budget%igrid_sim%sp_gpC)
         call this%prim_budget%igrid_sim%Pade6opZ%ddz_C2C(this%prim_budget%igrid_sim%cbuffzC(:,:,:,1),this%prim_budget%igrid_sim%cbuffzC(:,:,:,2),0,0)
         call transpose_z_to_y(this%prim_budget%igrid_sim%cbuffzC(:,:,:,2),this%prim_budget%igrid_sim%cbuffyC(:,:,:,1),this%prim_budget%igrid_sim%sp_gpC)
         call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
         call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdz)
         
     end subroutine 
 
     subroutine interp_Edge2Cell(this, fE, fC)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(in) :: fE
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: fC
 
         call transpose_x_to_y(fE,this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),this%prim_budget%igrid_sim%gpE)
         call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%gpE)
         call this%prim_budget%igrid_sim%Pade6opZ%interpz_E2C(this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,2),0,0)
         call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzC(:,:,:,2),this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),fC,this%prim_budget%igrid_sim%gpC)
         
     end subroutine 
 
     subroutine interp_Cell2Edge(this, fC, fE)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(in) :: fC
         real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(out) :: fE
 
         call transpose_x_to_y(fC,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call this%prim_budget%igrid_sim%Pade6opZ%interpz_C2E(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),0,0)
         call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),this%prim_budget%igrid_sim%gpE)
         call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),fE,this%prim_budget%igrid_sim%gpE)
         
     end subroutine 
         
     subroutine multiply_CellFieldsOnEdges(this, f1C, f2C, fmultC)
         class(budgets_time_avg_deficit), intent(inout) :: this
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(in) :: f1C,f2C
         real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)), intent(out) :: fmultC
 
         ! interpolate 1st Cell field
         call transpose_x_to_y(f1C,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call this%prim_budget%igrid_sim%Pade6opZ%interpz_C2E(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),0,0)
 
         ! interpolate 2nd Cell field
         call transpose_x_to_y(f2C,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call this%prim_budget%igrid_sim%Pade6opZ%interpz_C2E(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,2),0,0)
 
         ! multiply on Edges and interpolate back to Cells
         this%prim_budget%igrid_sim%rbuffzE(:,:,:,1) = this%prim_budget%igrid_sim%rbuffzE(:,:,:,1) * this%prim_budget%igrid_sim%rbuffzE(:,:,:,2)
         call this%prim_budget%igrid_sim%Pade6opZ%interpz_E2C(this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),0,0)
         call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
         call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),fmultC,this%prim_budget%igrid_sim%gpC)
 
     end subroutine 
 end module 
