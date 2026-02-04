module budgets_time_avg_deficit_compact_mod
    use kind_parameters, only: rkind, clen
    use decomp_2d
    use budgets_time_avg_mod, only: budgets_time_avg
    use exits, only: message, GracefulExit
    use constants, only: zero, half, two
    use mpi
    use incompressibleGrid, only : uBC_bottom, uBC_top, vBC_bottom, vBC_top, wBC_bottom, wBC_top, &
                      TBC_bottom, TBC_top, UWBC_bottom, UWBC_top, VWBC_bottom, VWBC_top, &
                      WTBC_bottom, WTBC_top
 
    implicit none 
 
    private
    public :: budgets_time_avg_deficit_compact

    ! Comments here

    type :: budgets_time_avg_deficit_compact
        private
        integer :: run_id, nx, ny, nz
        logical :: do_budget0=.false., do_budget1=.false., do_budget2=.false., do_budget3=.false.
        
        type(budgets_time_avg), pointer :: pre_budget, prim_budget
        
        real(rkind), dimension(:,:,:,:), allocatable :: budget_0, budget_1, budget_2, budget_3
        integer :: size_budget_0, size_budget_1, size_budget_2, size_budget_3
        real(rkind), dimension(:,:,:,:), allocatable :: MCG
        logical :: doMCG = .false.
        integer :: counter
        real(rkind) :: timeSum, weight
        character(len=clen) :: budgets_dir
        logical :: time_weighted_average=.false.

        logical :: useWindTurbines=.false., isStratified=.true., useCoriolis=.false.
        integer :: tidx_dump 
        integer :: tidx_compute
        integer :: tidx_budget_start 
        real(rkind) :: time_budget_start 
        logical :: do_budgets
        logical :: forceDump

        ! Avoid allocating a new holder of delta_tauij with every call to AssembleBudget3
        real(rkind), dimension(:,:,:,:), allocatable :: delta_tauij
 
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
        procedure, private  :: AssembleBudget1
        procedure, private  :: AssembleBudget2
        procedure, private  :: AssembleBudget3
        procedure, private  :: AssembleMCG
        procedure, private  :: restartMCG
   
        procedure, private  :: getProductOfMeans
        ! procedure, private  :: writeTimeSum
        ! procedure, private  :: readTimeSum

        procedure, private :: ddx_R2R
        procedure, private :: ddy_R2R
        procedure, private :: ddz_R2R
        !procedure, private :: ddz_C2R
        procedure, private :: dealias
        procedure, private :: interp_Edge2Cell
        ! procedure, private :: interp_Cell2Edge
        ! procedure, private :: multiply_CellFieldsOnEdges
        ! procedure, private :: multiply_edges_interp_cell
     end type

    contains

    subroutine init(this, pre_budget, primary_inputfile, prim_budget) 
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        character(len=*), intent(in) :: primary_inputfile 
        type(budgets_time_avg), intent(inout), target :: pre_budget, prim_budget        
        character(len=clen) :: budgets_dir = "NULL"
        character(len=clen) :: restart_dir = "NULL"
        integer :: ioUnit, ierr,  restart_tid = 0, restart_rid = 0, restart_counter = 0
        logical :: restart_budgets = .false. 
        integer :: tidx_compute = 10000, tidx_dump = 10000, tidx_budget_start = -100
        real(rkind) :: time_budget_start = -1.0d0
        logical :: use_time_weighted_average=.false.
        logical :: do_budgets = .false. 
        logical :: do_budget0=.false., do_budget1=.false., do_budget2=.false., do_budget3=.false.
        namelist /BUDGET_TIME_AVG_DEFICIT_COMPACT/ budgets_dir, restart_budgets, restart_dir, &
            restart_rid, restart_tid, restart_counter, tidx_dump, tidx_compute, do_budgets, &
            use_time_weighted_average, tidx_budget_start, time_budget_start, &
            do_budget0, do_budget1, do_budget2, do_budget3 

        ! STEP 1: Read in inputs, link pointers and allocate budget vectors
        ioUnit = 534
        open(unit=ioUnit, file=trim(primary_inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_TIME_AVG_DEFICIT_COMPACT)
        close(ioUnit)

        this%pre_budget => pre_budget 
        this%prim_budget => prim_budget
        this%run_id = this%prim_budget%igrid_sim%runid
        this%nx = this%prim_budget%igrid_sim%gpC%xsz(1)
        this%ny = this%prim_budget%igrid_sim%gpC%xsz(2)
        this%nz = this%prim_budget%igrid_sim%gpC%xsz(3)  ! centered grid x, y, z
        ! this%nxE = this%prim_budget%igrid_sim%gpE%xsz(1)
        ! this%nyE = this%prim_budget%igrid_sim%gpE%xsz(2)
        ! this%nzE = this%prim_budget%igrid_sim%gpE%xsz(3) 
        this%do_budgets = do_budgets
        this%tidx_dump = tidx_dump
        this%tidx_compute = tidx_compute
        this%tidx_budget_start = tidx_budget_start  
        this%time_budget_start = time_budget_start  
        this%useWindTurbines = this%prim_budget%igrid_sim%useWindTurbines
        this%isStratified    = this%prim_budget%igrid_sim%isStratified
        this%useCoriolis    = this%prim_budget%igrid_sim%useCoriolis
        ! Deactivate time-weighted sum till time-averaged budgets are weighted similarily
        !this%time_weighted_average = use_time_weighted_average
        this%time_weighted_average = .False.
        this%forceDump = .false.
        this%do_budget0 = do_budget0
        this%do_budget1 = do_budget1
        this%do_budget2 = do_budget2
        this%do_budget3 = do_budget3
        
        if(this%do_budget1)this%do_budget0=.true.
        if(this%do_budget2)this%do_budget0=.true.
        if(this%do_budget3)then
            this%do_budget0=.true.
            this%do_budget1=.true.
            this%do_budget2=.true.
        end if
        if(this%do_budget2) this%doMCG = .true.
        this%budgets_dir = budgets_dir

        if(this%do_budgets) then 
            if((this%tidx_budget_start > 0) .and. (this%time_budget_start > zero)) then
                call GracefulExit("Both tidx_budget_start and time_budget_start in budget_time_avg are positive. Turn one negative", 100)
            endif

            if(this%do_budget0)then
                if(this%useWindTurbines)then
                    this%size_budget_0 = 22
                else
                    this%size_budget_0 = 20
                end if
                allocate(this%budget_0(this%nx,this%ny,this%nz,this%size_budget_0))
            end if

            if(this%do_budget1)then
                this%size_budget_1 = 15
                allocate(this%budget_1(this%nx,this%ny,this%nz,this%size_budget_1))
            end if

            if(this%do_budget2)then
                this%size_budget_2 = 15
                allocate(this%budget_2(this%nx,this%ny,this%nz,this%size_budget_2))
            end if

            if(this%do_budget3)then
                if(this%useWindTurbines)then
                    this%size_budget_3 = 21
                else
                    this%size_budget_3 = 19
                end if
                allocate(this%budget_3(this%nx,this%ny,this%nz,this%size_budget_3))
                allocate(this%delta_tauij(this%nx,this%ny,this%nz,6))
            end if

            if(this%doMCG)allocate(this%MCG(this%nx,this%ny,this%nz,18))

            if ((trim(budgets_dir) .eq. "null") .or.(trim(budgets_dir) .eq. "NULL")) then 
                this%budgets_dir = this%prim_budget%igrid_sim%outputDir
            end if 

            if ((trim(restart_dir) .eq. "null") .or.(trim(restart_dir) .eq. "NULL")) then
                restart_dir = this%budgets_dir
            end if 

            if (restart_budgets) then
                call message(0,"Budget deficit restart")
                call this%RestartBudget(restart_dir, restart_rid, restart_tid, restart_counter)
            else
                call this%resetBudget()
            end if 
        end if 
     end subroutine

     subroutine doBudgets(this, forceDump)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        logical, intent(in), optional :: forceDump

        if(present(forceDump)) then
            this%forceDump = forceDump
        endif

        if(this%prim_budget%igrid_sim%tsim > this%prim_budget%igrid_sim%tstop) then
            this%forceDump = .TRUE.
        endif

        if (this%do_budgets)  then
            if( ( (this%tidx_budget_start>0) .and. (this%prim_budget%igrid_sim%step>this%tidx_budget_start) ) .or. &
                ( (this%time_budget_start>0) .and. (this%prim_budget%igrid_sim%tsim>this%time_budget_start) ) ) then
        
                if (mod(this%prim_budget%igrid_sim%step,this%tidx_compute) .eq. 0) then
                    call this%updateBudget()
                end if

                if ((mod(this%prim_budget%igrid_sim%step,this%tidx_dump) .eq. 0) .or. this%forceDump) then
                    call this%dumpBudget()
                    call message(0,"Dumped a compact deficit budget file")
                end if 
            end if 
        end if 

        this%forceDump = .false. ! reset to default value
    end subroutine

    subroutine updateBudget(this)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this

        ! This step computes the pressure field of the primary and precursor simulations.
        call this%prim_budget%igrid_sim%getMomentumTerms()  
        call this%pre_budget%igrid_sim%getMomentumTerms()  

        ! Interpolate SGS stresses to cells
        call this%pre_budget%igrid_sim%sgsmodel%populate_tauij_E_to_C()
        call this%prim_budget%igrid_sim%sgsmodel%populate_tauij_E_to_C()
        this%delta_tauij = this%prim_budget%igrid_sim%tauSGS_ij - this%pre_budget%igrid_sim%tauSGS_ij

        if(this%doMCG) call this%AssembleMCG()
        if(this%do_budget0) call this%AssembleBudget0()
        if(this%do_budget1) call this%AssembleBudget1()
        if(this%do_budget2) call this%AssembleBudget2()
        if(this%do_budget3) call this%AssembleBudget3()

        this%counter = this%counter + 1
    end subroutine 

    subroutine DumpBudget(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind) :: totalWeight
        integer :: idx, budgetid, budgetsize
        real(rkind), dimension(:,:,:), pointer :: buffer
        real(rkind), dimension(:,:,:,:), pointer :: budget
        logical :: doBudget

        totalWeight = real(this%counter,rkind) + 1.d-18

        ! Cell x-pencil buffers 
        ! Buffers 1 and 2 are used locally inside getProductOfMeans
        buffer => this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)

        ! Convert assembled budgets to mean instead of sum
        if(this%do_budget0) this%budget_0 = this%budget_0/totalWeight
        if(this%do_budget1) this%budget_1 = this%budget_1/totalWeight
        if(this%do_budget2) this%budget_2 = this%budget_2/totalWeight
        if(this%do_budget3) this%budget_3 = this%budget_3/totalWeight
        if(this%doMCG) this%MCG = this%MCG/totalWeight
        this%pre_budget%budget_0 = this%pre_budget%budget_0/totalWeight
        this%pre_budget%budget_1 = this%pre_budget%budget_1/totalWeight

        ! Budget 0
        if(this%do_budget0)then
            budgetid = 0
            do idx = 1, this%size_budget_0
                if((idx.eq.15).or.(idx.eq.16))then
                    if(.not. this%useCoriolis)cycle
                end if
                if((idx.eq.5).or.(idx.eq.17))then
                    if(.not. this%isStratified)cycle
                end if
                call this%dump_budget_field(this%budget_0(:,:,:,idx), idx, budgetid)
            end do
        end if

        ! Dealias budgets 1-3 as they hold product of multiple fields
        if(this%do_budget1)then
            do idx = 1, this%size_budget_1
                call this%dealias(this%budget_1(:,:,:,idx))
            end do
        end if
        if(this%do_budget2)then
            do idx = 1, this%size_budget_2
                call this%dealias(this%budget_2(:,:,:,idx))
            end do
        end if
        if(this%do_budget3)then
            do idx = 1, this%size_budget_3
                call this%dealias(this%budget_3(:,:,:,idx))
            end do
        end if

        do budgetid=1,3
            select case(budgetid)
            case(1)
                budget => this%budget_1
                budgetsize = this%size_budget_1
                doBudget = this%do_budget1
            case(2)
                budget => this%budget_2
                budgetsize = this%size_budget_2
                doBudget = this%do_budget2
            case(3)
                budget => this%budget_3
                budgetsize = this%size_budget_3
                doBudget = this%do_budget3
            end select

            if(doBudget)then
                do idx = 1,budgetsize

                    ! Skip Buoyancy covariance in TKE budget
                    if(budgetid.eq.3)then
                        if((idx.eq.10).or.(idx.eq.11).or.(idx.eq.12))then
                            if(.not. this%isStratified) cycle
                        end if
                    end if

                    ! Get the product of means. buffer is dealiased inside getProductOfMeans
                    call this%getProductOfMeans(budgetid, idx, buffer)

                    ! Remove product of means. The original budget is not impacted
                    buffer = budget(:,:,:,idx) - buffer

                    ! Dump
                    call this%dump_budget_field(buffer, idx, budgetid)
                end do
            end if
        end do

        ! Return to summing
        if(this%do_budget0) this%budget_0 = this%budget_0*totalWeight
        if(this%do_budget1) this%budget_1 = this%budget_1*totalWeight
        if(this%do_budget2) this%budget_2 = this%budget_2*totalWeight
        if(this%do_budget3) this%budget_3 = this%budget_3*totalWeight
        if(this%doMCG) this%MCG = this%MCG*totalWeight
        this%pre_budget%budget_0 = this%pre_budget%budget_0*totalWeight
        this%pre_budget%budget_1 = this%pre_budget%budget_1*totalWeight
    end subroutine

    ! ---------------------- Mean Cell Gradients (MCG) ------------------------
    subroutine AssembleMCG(this)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this  
        this%MCG(:,:,:,1:9) = this%MCG(:,:,:,1:9) + this%prim_budget%igrid_sim%duidxjC(:,:,:,1:9) - this%pre_budget%igrid_sim%duidxjC(:,:,:,1:9)
        this%MCG(:,:,:,10:18) = this%MCG(:,:,:,10:18) + this%pre_budget%igrid_sim%duidxjC(:,:,:,1:9)
    end subroutine    

    ! ---------------------- Budget 0 ------------------------
    subroutine AssembleBudget0(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rbuffxE1, rbuffxC1, rbuffxC2 
        complex(rkind), dimension(:,:,:), pointer :: cbuffyE1, cbuffyC1
        
        ! Link pointers
        cbuffyE1 => this%prim_budget%igrid_sim%cbuffyE(:,:,:,1)
        cbuffyC1 => this%prim_budget%igrid_sim%cbuffyC(:,:,:,2) ! 1 is used in ddx, ddy, ddz routines        
        rbuffxE1 => this%prim_budget%igrid_sim%rbuffxE(:,:,:,1)
        rbuffxC1 => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        rbuffxC2 => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)        
        
        ! STEP 1: Compute mean Delta U, Delta V, and Delta W
        this%budget_0(:,:,:,1) = this%budget_0(:,:,:,1) + (this%prim_budget%igrid_sim%u  - this%pre_budget%igrid_sim%u)
        this%budget_0(:,:,:,2) = this%budget_0(:,:,:,2) + (this%prim_budget%igrid_sim%v  - this%pre_budget%igrid_sim%v)
        this%budget_0(:,:,:,3) = this%budget_0(:,:,:,3) + (this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC)

        ! STEP 2: Pressure
        this%budget_0(:,:,:,4) = this%budget_0(:,:,:,4) + (this%prim_budget%igrid_sim%pressure - this%pre_budget%igrid_sim%pressure)

        ! STEP 3: Potential temperature
        if (this%isStratified)then 
            this%budget_0(:,:,:,5) = this%budget_0(:,:,:,5) + (this%prim_budget%igrid_sim%T - this%pre_budget%igrid_sim%T)
            
            cbuffyE1 = this%prim_budget%wb - this%pre_budget%wb
            call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
            call this%interp_Edge2Cell(rbuffxE1, rbuffxC1, TBC_bottom, TBC_top)
            this%budget_0(:,:,:,17) = this%budget_0(:,:,:,17) + rbuffxC1
        end if

        ! Step 4: SGS stresses (also viscous stress if finite reynolds number is being used)
        this%budget_0(:,:,:,6:11) = this%budget_0(:,:,:,6:11) + this%delta_tauij

        ! Step 5: SGS stress gradients
        ! Reverse signs of usgs, vsgs, wsgs
        cbuffyC1 = this%pre_budget%usgs - this%prim_budget%usgs
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)
        this%budget_0(:,:,:,12) = this%budget_0(:,:,:,12) + rbuffxC1

        cbuffyC1 = this%pre_budget%vsgs - this%prim_budget%vsgs
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)
        this%budget_0(:,:,:,13) = this%budget_0(:,:,:,13) + rbuffxC1

        ! wsgs is odd
        cbuffyE1 = this%pre_budget%wsgs - this%prim_budget%wsgs
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, rbuffxC1, -1, -1)
        this%budget_0(:,:,:,14) = this%budget_0(:,:,:,14) + rbuffxC1
        
        ! Step 6: Coriolis
        if(this%useCoriolis) then
            ! Remove the geostrophic forcing term from exported Coriolis force  
            call this%pre_budget%igrid_sim%get_geostrophic_forcing(rbuffxC1, rbuffxC2)  
            this%budget_0(:,:,:,15) = this%budget_0(:,:,:,15) + rbuffxC1
            this%budget_0(:,:,:,16) = this%budget_0(:,:,:,16) + rbuffxC2      

            call this%prim_budget%igrid_sim%get_geostrophic_forcing(rbuffxC1, rbuffxC2)
            this%budget_0(:,:,:,15) = this%budget_0(:,:,:,15) - rbuffxC1
            this%budget_0(:,:,:,16) = this%budget_0(:,:,:,16) - rbuffxC2              
            
            ! Coriolis term, X 
            cbuffyC1 = this%prim_budget%ucor - this%pre_budget%ucor      
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)
            this%budget_0(:,:,:,15) = this%budget_0(:,:,:,15) + rbuffxC1

            ! Coriolis term, Y       
            cbuffyC1 = this%prim_budget%vcor - this%pre_budget%vcor
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)            
            this%budget_0(:,:,:,16) = this%budget_0(:,:,:,16) + rbuffxC1
        end if

        ! Step 7: Pressure gradient force
        ! px sign is reversed
        cbuffyC1 = this%pre_budget%px - this%prim_budget%px
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)
        this%budget_0(:,:,:,18) = this%budget_0(:,:,:,18) + rbuffxC1

        ! py sign is reversed
        cbuffyC1 = this%pre_budget%py - this%prim_budget%py
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)
        this%budget_0(:,:,:,19) = this%budget_0(:,:,:,19) + rbuffxC1

        ! pz sign is reversed
        ! pz is odd
        cbuffyE1 = this%pre_budget%pz - this%prim_budget%pz
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, rbuffxC1, -1, -1)
        this%budget_0(:,:,:,20) = this%budget_0(:,:,:,20) + rbuffxC1         

        ! Step 8: turbine forcing
        if(this%useWindTurbines)then        
            cbuffyC1 = this%prim_budget%uturb - this%pre_budget%uturb
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)
            this%budget_0(:,:,:,21) = this%budget_0(:,:,:,21) + rbuffxC1

            cbuffyC1 = this%prim_budget%vturb - this%pre_budget%vturb
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, rbuffxC1)
            this%budget_0(:,:,:,22) = this%budget_0(:,:,:,22) + rbuffxC1
        end if

        nullify(rbuffxE1, rbuffxC1, rbuffxC2, cbuffyC1, cbuffyE1)
    end subroutine

    ! ---------------------- Budget 1 ------------------------
    subroutine AssembleBudget1(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: du, dv, dw, duE, dvE, dwE, buffer, buffE

        ! Cell x-pencil buffers 
        du =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        dv =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
        dw =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,3)
        buffer =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)
        
        ! Edge x-pencil buffers (only 2 are allocated in igrid.F90)
        duE => this%prim_budget%igrid_sim%rbuffxE(:,:,:,1)
        dvE => this%prim_budget%igrid_sim%rbuffxE(:,:,:,2)
        dwE => this%pre_budget%igrid_sim%rbuffxE(:,:,:,1)
        buffE => this%pre_budget%igrid_sim%rbuffxE(:,:,:,2)

        ! Perturbation fields
        du = this%prim_budget%igrid_sim%u  - this%pre_budget%igrid_sim%u
        dv = this%prim_budget%igrid_sim%v  - this%pre_budget%igrid_sim%v
        dw = this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC
        
        duE = this%prim_budget%igrid_sim%uE - this%pre_budget%igrid_sim%uE
        dvE = this%prim_budget%igrid_sim%vE - this%pre_budget%igrid_sim%vE
        dwE = this%prim_budget%igrid_sim%w  - this%pre_budget%igrid_sim%w

        ! Reynolds stresses
        this%budget_1(:,:,:,1) = this%budget_1(:,:,:,1) + du * du
        this%budget_1(:,:,:,2) = this%budget_1(:,:,:,2) + du * dv   
        buffE = duE * dwE; call this%interp_Edge2Cell(buffE, buffer, UWBC_bottom, UWBC_top)  
        this%budget_1(:,:,:,3) = this%budget_1(:,:,:,3) + buffer
        this%budget_1(:,:,:,4) = this%budget_1(:,:,:,4) + dv * dv
        buffE = dvE * dwE; call this%interp_Edge2Cell(buffE, buffer, VWBC_bottom, VWBC_top)
        this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) + buffer
        this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) + dw * dw
         
        ! Mixed Reynolds stresses
        this%budget_1(:,:,:,7)  = this%budget_1(:,:,:,7) + du * this%pre_budget%igrid_sim%u
        this%budget_1(:,:,:,8)  = this%budget_1(:,:,:,8) + du * this%pre_budget%igrid_sim%v
        this%budget_1(:,:,:,9)  = this%budget_1(:,:,:,9) + dv * this%pre_budget%igrid_sim%u

        buffE = duE * this%pre_budget%igrid_sim%w; ; call this%interp_Edge2Cell(buffE, buffer, UWBC_bottom, UWBC_top)
        this%budget_1(:,:,:,10) = this%budget_1(:,:,:,10) + buffer
        buffE = dwE * this%pre_budget%igrid_sim%uE; call this%interp_Edge2Cell(buffE, buffer, UWBC_bottom, UWBC_top)
        this%budget_1(:,:,:,11) = this%budget_1(:,:,:,11) + buffer
        this%budget_1(:,:,:,12) = this%budget_1(:,:,:,12) + dv * this%pre_budget%igrid_sim%v
        buffE = dvE * this%pre_budget%igrid_sim%w; call this%interp_Edge2Cell(buffE, buffer, VWBC_bottom, VWBC_top)
        this%budget_1(:,:,:,13) = this%budget_1(:,:,:,13) + buffer
        buffE = dwE * this%pre_budget%igrid_sim%vE; call this%interp_Edge2Cell(buffE, buffer, VWBC_bottom, VWBC_top)
        this%budget_1(:,:,:,14) = this%budget_1(:,:,:,14) + buffer
        this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) + dw * this%pre_budget%igrid_sim%wC

        nullify(du, dv, dw, duE, dvE, dwE, buffer, buffE)
    end subroutine

    ! ---------------------- Budget 2 ------------------------
    subroutine AssembleBudget2(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: du, dv, buffC, dw
        real(rkind), dimension(:,:,:), pointer :: ubase, vbase, wbase
        real(rkind), dimension(:,:,:), pointer :: dudxC_prim, dudyC_prim, dudzC_prim, dudxC_pre, dudyC_pre, dudzC_pre
        real(rkind), dimension(:,:,:), pointer :: dvdxC_prim, dvdyC_prim, dvdzC_prim, dvdxC_pre, dvdyC_pre, dvdzC_pre
        real(rkind), dimension(:,:,:), pointer :: dwdxC_prim, dwdyC_prim, dwdzC_prim, dwdxC_pre, dwdyC_pre, dwdzC_pre

        ! Cell x-pencil buffers 
        du => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        dv => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)        
        dw => this%prim_budget%igrid_sim%rbuffxC(:,:,:,3)
        buffC => this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)

        ! Perturbation fields
        du = this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u
        dv = this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v
        dw = this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC

        ! Base-flow fields
        ubase => this%pre_budget%igrid_sim%u
        vbase => this%pre_budget%igrid_sim%v
        wbase => this%pre_budget%igrid_sim%wC

        ! Primary simulation:
        dudxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,1)
        dudyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,2)
        dudzC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,3)
        dvdxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,4)
        dvdyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,5)
        dvdzC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,6)
        dwdxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,7)
        dwdyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,8)
        dwdzC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,9)

        ! Precursor simulation:
        dudxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,1)
        dudyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,2)
        dudzC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,3)
        dvdxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,4)
        dvdyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,5)
        dvdzC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,6) 
        dwdxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,7)
        dwdyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,8)
        dwdzC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,9)       

        ! delta u_j d_j(delta u)
        this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + du * (dudxC_prim - dudxC_pre) + dv * (dudyC_prim - dudyC_pre) + dw * (dudzC_prim - dudzC_pre)
        
        ! delta u_j d_j(delta v)
        this%budget_2(:,:,:,2) = this%budget_2(:,:,:,2) + du * (dvdxC_prim - dvdxC_pre) + dv * (dvdyC_prim - dvdyC_pre) + dw * (dvdzC_prim - dvdzC_pre)
        
        ! delta u_j d_j(delta w)
        this%budget_2(:,:,:,3) = this%budget_2(:,:,:,3) + du * (dwdxC_prim - dwdxC_pre) + dv * (dwdyC_prim - dwdyC_pre) + dw * (dwdzC_prim - dwdzC_pre)
        
        ! delta u_j d_j(base u)
        this%budget_2(:,:,:,4) = this%budget_2(:,:,:,4) + du * dudxC_pre + dv * dudyC_pre + dw * dudzC_pre
        
        ! delta u_j d_j(base v)
        this%budget_2(:,:,:,5) = this%budget_2(:,:,:,5) + du * dvdxC_pre + dv * dvdyC_pre + dw * dvdzC_pre
        
        ! delta u_j d_j(base w)
        this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + du * dwdxC_pre + dv * dwdyC_pre + dw * dwdzC_pre
        
        ! base u_j d_j(delta u)
        this%budget_2(:,:,:,7) = this%budget_2(:,:,:,7) + ubase * (dudxC_prim - dudxC_pre) + vbase * (dudyC_prim - dudyC_pre) + wbase * (dudzC_prim - dudzC_pre)
        
        ! base u_j d_j(delta v)
        this%budget_2(:,:,:,8) = this%budget_2(:,:,:,8) + ubase * (dvdxC_prim - dvdxC_pre) + vbase * (dvdyC_prim - dvdyC_pre) + wbase * (dvdzC_prim - dvdzC_pre)
        
        ! base u_j d_j(delta w)
        this%budget_2(:,:,:,9) = this%budget_2(:,:,:,9) + ubase * (dwdxC_prim - dwdxC_pre) + vbase * (dwdyC_prim - dwdyC_pre) + wbase * (dwdzC_prim - dwdzC_pre)
        
        ! base u_j d_j(base u)
        this%budget_2(:,:,:,10) = this%budget_2(:,:,:,10) + ubase * dudxC_pre + vbase * dudyC_pre + wbase * dudzC_pre
        
        ! base u_j d_j(base v)
        this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + ubase * dvdxC_pre + vbase * dvdyC_pre + wbase * dvdzC_pre
       
        ! base u_j d_j(base w)
        this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + ubase * dwdxC_pre + vbase * dwdyC_pre + wbase * dwdzC_pre  
        
        ! base u_i d_1(delta u_i)
        this%budget_2(:,:,:,13) = this%budget_2(:,:,:,13) + ubase * (dudxC_prim - dudxC_pre) + vbase * (dvdxC_prim - dvdxC_pre) + wbase * (dwdxC_prim - dwdxC_pre)
        
        ! base u_i d_2(delta u_i)
        this%budget_2(:,:,:,14) = this%budget_2(:,:,:,14) + ubase * (dudyC_prim - dudyC_pre) + vbase * (dvdyC_prim - dvdyC_pre) + wbase * (dwdyC_prim - dwdyC_pre)
        
        ! base u_i d_3(delta u_i)
        this%budget_2(:,:,:,15) = this%budget_2(:,:,:,15) + ubase * (dudzC_prim - dudzC_pre) + vbase * (dvdzC_prim - dvdzC_pre) + wbase * (dwdzC_prim - dwdzC_pre)    

        ! Release memory        
        nullify(du, dv, dw, buffC)
        nullify(ubase, vbase, wbase)
        nullify(dudxC_prim, dudyC_prim, dudzC_prim, dudxC_pre, dudyC_pre, dudzC_pre)
        nullify(dvdxC_prim, dvdyC_prim, dvdzC_prim, dvdxC_pre, dvdyC_pre, dvdzC_pre)
        nullify(dwdxC_prim, dwdyC_prim, dwdzC_prim, dwdxC_pre, dwdyC_pre, dwdzC_pre)
    end subroutine

    ! ---------------------- Budget 3 ------------------------
    subroutine AssembleBudget3(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: du, dv, dw
        real(rkind), dimension(:,:,:), pointer :: ubase, vbase, wbase
        real(rkind), dimension(:,:,:), pointer :: rbuffxE1, rbuffxE2, buffer
        complex(rkind), dimension(:,:,:), pointer :: cbuffyE1, cbuffyC1
        real(rkind), dimension(:,:,:), pointer :: dudxC_prim, dudyC_prim, dudzC_prim, dudxC_pre, dudyC_pre, dudzC_pre
        real(rkind), dimension(:,:,:), pointer :: dvdxC_prim, dvdyC_prim, dvdzC_prim, dvdxC_pre, dvdyC_pre, dvdzC_pre
        real(rkind), dimension(:,:,:), pointer :: dwdxC_prim, dwdyC_prim, dwdzC_prim, dwdxC_pre, dwdyC_pre, dwdzC_pre
        real(rkind), dimension(:,:,:,:), pointer :: base_tauij

        ! Cell x-pencil buffers 
        du => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        dv => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
        dw => this%prim_budget%igrid_sim%rbuffxC(:,:,:,3)
        buffer => this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)
        
        ! Cell y-pencil buffer 
        cbuffyC1 => this%prim_budget%igrid_sim%cbuffyC(:,:,:,2) ! 1 is used in ddx, ddy, ddz routines 

        ! Edge x-pencil buffer
        rbuffxE1 => this%prim_budget%igrid_sim%rbuffxE(:,:,:,1)
        rbuffxE2 => this%prim_budget%igrid_sim%rbuffxE(:,:,:,2)

        ! Edge y-pencil buffer
        cbuffyE1 => this%prim_budget%igrid_sim%cbuffyE(:,:,:,1)

        ! Perturbation fields
        du = this%prim_budget%igrid_sim%u  - this%pre_budget%igrid_sim%u
        dv = this%prim_budget%igrid_sim%v  - this%pre_budget%igrid_sim%v
        dw = this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC

        ubase => this%pre_budget%igrid_sim%u
        vbase => this%pre_budget%igrid_sim%v
        wbase => this%pre_budget%igrid_sim%wC

        ! Primary simulation gradients:
        dudxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,1)
        dudyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,2)
        dudzC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,3)
        dvdxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,4)
        dvdyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,5)
        dvdzC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,6)
        dwdxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,7)
        dwdyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,8)
        dwdzC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,9)
        
        ! Precursor simulation gradients:
        dudxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,1)
        dudyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,2)
        dudzC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,3)
        dvdxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,4)
        dvdyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,5)
        dvdzC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,6) 
        dwdxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,7)
        dwdyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,8)
        dwdzC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,9)
        base_tauij => this%pre_budget%igrid_sim%tauSGS_ij
        
        ! Term 1: delta u_j' d_j(delta p')
        ! Term 2: base  u_j' d_j(delta p')        
        ! px, py, pz signs are reversed
        cbuffyC1 = this%pre_budget%px - this%prim_budget%px
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer)
        this%budget_3(:,:,:,1)=this%budget_3(:,:,:,1)+ buffer * du
        this%budget_3(:,:,:,2)=this%budget_3(:,:,:,2)+ buffer * ubase

        cbuffyC1 = this%pre_budget%py - this%prim_budget%py
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer)
        this%budget_3(:,:,:,1)=this%budget_3(:,:,:,1)+ buffer * dv
        this%budget_3(:,:,:,2)=this%budget_3(:,:,:,2)+ buffer * vbase

        ! pz is odd
        cbuffyE1 = this%pre_budget%pz - this%prim_budget%pz
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, buffer, -1, -1)
        this%budget_3(:,:,:,1)=this%budget_3(:,:,:,1)+ buffer * dw        
        this%budget_3(:,:,:,2)=this%budget_3(:,:,:,2)+ buffer * wbase

        ! Term 3: delta u_j' d_j(base p')
        ! px, py, pz signs are reversed
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%px, buffer)
        this%budget_3(:,:,:,3)=this%budget_3(:,:,:,3)- buffer * du

        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%py, buffer)
        this%budget_3(:,:,:,3)=this%budget_3(:,:,:,3)- buffer * dv

        ! pz is odd
        call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%pz, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, buffer, -1, -1)
        this%budget_3(:,:,:,3)=this%budget_3(:,:,:,3)- buffer * dw

        ! Term 4: d_j(base  u_i' * delta tau_ij') [SGS transport] 
        ! Term 6: d_j(delta u_i' * delta tau_ij')  [SGS transport]
        ! sign of usgs, vsgs, and wsgs are reversed.
        cbuffyC1 = this%pre_budget%usgs - this%prim_budget%usgs
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer)
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer * ubase 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer * du

        cbuffyC1 = this%pre_budget%vsgs - this%prim_budget%vsgs
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer) 
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer * vbase  
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer * dv

        ! wsgs is odd
        cbuffyE1 = this%pre_budget%wsgs - this%prim_budget%wsgs
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, buffer, -1, -1)
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer * wbase
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer * dw

        ! Term 5: d_j(delta  u_i' base tau_ij') [SGS transport]         
        ! sign of usgs, vsgs, and wsgs are reversed. 
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%usgs, buffer)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - buffer * du

        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%vsgs, buffer)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - buffer * dv 

        ! wsgs is odd
        call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%wsgs, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, buffer, -1, -1)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - buffer * dw

        ! The remaining of B3(4) is exactly B3(7). Calculation is done once
        ! Term 4: d_j(base  u_i' * delta tau_ij')  [SGS transport]
        ! Term 7: delta tau_ij' d_j(base u_i')     [SGS dissipation]
        buffer = dudxC_pre*this%delta_tauij(:,:,:,1) + dudyC_pre*this%delta_tauij(:,:,:,2) + dudzC_pre*this%delta_tauij(:,:,:,3)+&
                 dvdxC_pre*this%delta_tauij(:,:,:,2) + dvdyC_pre*this%delta_tauij(:,:,:,4) + dvdzC_pre*this%delta_tauij(:,:,:,5)+&
                 dwdxC_pre*this%delta_tauij(:,:,:,3) + dwdyC_pre*this%delta_tauij(:,:,:,5) + dwdzC_pre*this%delta_tauij(:,:,:,6)
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer        

        ! The remaining of B3(5) is exactly B3(8). Calculation is done once
        ! Term 5: d_j(delta  u_i' base tau_ij')    [SGS transport]  
        ! Term 8: base  tau_ij' * d_j(delta u_i')  [SGS dissipation]
        buffer = (dudxC_prim-dudxC_pre)*base_tauij(:,:,:,1)+(dudyC_prim-dudyC_pre)*base_tauij(:,:,:,2)+(dudzC_prim-dudzC_pre)*base_tauij(:,:,:,3)+&
                 (dvdxC_prim-dvdxC_pre)*base_tauij(:,:,:,2)+(dvdyC_prim-dvdyC_pre)*base_tauij(:,:,:,4)+(dvdzC_prim-dvdzC_pre)*base_tauij(:,:,:,5)+&
                 (dwdxC_prim-dwdxC_pre)*base_tauij(:,:,:,3)+(dwdyC_prim-dwdyC_pre)*base_tauij(:,:,:,5)+(dwdzC_prim-dwdzC_pre)*base_tauij(:,:,:,6)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer

        ! The remaining of B3(6) is exactly B3(9). Calculation is done once
        ! Term 6: d_j(delta u_i' * delta tau_ij')  [SGS transport]
        ! Term 9: delta tau_ij' * d_j(delta u_i')  [SGS dissipation]
        buffer = (dudxC_prim-dudxC_pre)*this%delta_tauij(:,:,:,1)+(dudyC_prim-dudyC_pre)*this%delta_tauij(:,:,:,2)+(dudzC_prim-dudzC_pre)*this%delta_tauij(:,:,:,3)+&
                 (dvdxC_prim-dvdxC_pre)*this%delta_tauij(:,:,:,2)+(dvdyC_prim-dvdyC_pre)*this%delta_tauij(:,:,:,4)+(dvdzC_prim-dvdzC_pre)*this%delta_tauij(:,:,:,5)+&
                 (dwdxC_prim-dwdxC_pre)*this%delta_tauij(:,:,:,3)+(dwdyC_prim-dwdyC_pre)*this%delta_tauij(:,:,:,5)+(dwdzC_prim-dwdzC_pre)*this%delta_tauij(:,:,:,6)
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer

        ! Term 10: delta u_3' delta wb'
        ! Term 11: delta u_3' base wb'
        ! Term 12: base u_3' delta wb'
        ! Multiply on edges
        if(this%isStratified)then
            cbuffyE1 = this%prim_budget%wb - this%pre_budget%wb 
            call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
            
            rbuffxE2 = rbuffxE1 * (this%prim_budget%igrid_sim%w - this%pre_budget%igrid_sim%w)
            call this%interp_Edge2Cell(rbuffxE2, buffer, WTBC_bottom, WTBC_top)        
            this%budget_3(:,:,:,10) = this%budget_3(:,:,:,10) + buffer

            rbuffxE2 = rbuffxE1 * this%pre_budget%igrid_sim%w
            call this%interp_Edge2Cell(rbuffxE2, buffer, WTBC_bottom, WTBC_top)
            this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) + buffer
            
            call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%wb, rbuffxE1)
            rbuffxE2 = (this%prim_budget%igrid_sim%w - this%pre_budget%igrid_sim%w) * rbuffxE1
            call this%interp_Edge2Cell(rbuffxE2, buffer, WTBC_bottom, WTBC_top)
            this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) + buffer    
        end if  

        ! Term 13: base  u_i' delta u_j' d_j(base u_i')  [Turbulent transport of TKE]
        ! Term 17: delta u_i' delta u_j' d_j(base u_i')  [Turbulent transport of TKE]
        buffer = du * dudxC_pre + dv * dudyC_pre + dw * dudzC_pre
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + ubase * buffer
        this%budget_3(:,:,:,17) = this%budget_3(:,:,:,17) + du * buffer

        buffer = du * dvdxC_pre + dv * dvdyC_pre + dw * dvdzC_pre
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + vbase * buffer
        this%budget_3(:,:,:,17) = this%budget_3(:,:,:,17) + dv * buffer

        buffer = du * dwdxC_pre + dv * dwdyC_pre + dw * dwdzC_pre
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + wbase * buffer
        this%budget_3(:,:,:,17) = this%budget_3(:,:,:,17) + dw * buffer

        ! Term 14: base  u_i' base u_j' d_j(delta u_i')  [Turbulent transport of TKE]
        ! Term 18: delta u_i' base u_j' d_j(delta u_i')  [Turbulent transport of TKE]
        buffer = ubase*(dudxC_prim-dudxC_pre) + vbase*(dudyC_prim-dudyC_pre) + wbase*(dudzC_prim-dudzC_pre)
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + ubase * buffer
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + du * buffer

        buffer = ubase*(dvdxC_prim-dvdxC_pre) + vbase*(dvdyC_prim-dvdyC_pre) + wbase*(dvdzC_prim-dvdzC_pre)
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + vbase * buffer
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + dv * buffer 

        buffer = ubase*(dwdxC_prim-dwdxC_pre) + vbase*(dwdyC_prim-dwdyC_pre) + wbase*(dwdzC_prim-dwdzC_pre)
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + wbase * buffer
        this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + dw * buffer   

        ! Term 15: delta u_i' base u_j' d_j(base u_i')  [Turbulent transport of TKE]
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + &
            du*(ubase * dudxC_pre + vbase * dudyC_pre + wbase * dudzC_pre) + &
            dv*(ubase * dvdxC_pre + vbase * dvdyC_pre + wbase * dvdzC_pre) + &
            dw*(ubase * dwdxC_pre + vbase * dwdyC_pre + wbase * dwdzC_pre)

        ! Term 16: base  u_i' delta u_j' d_j(delta u_i')  [Turbulent transport of TKE]
        ! Term 19: delta u_i' delta u_j' d_j(delta u_i')  [Turbulent transport of TKE]
        buffer = du*(dudxC_prim-dudxC_pre) + dv*(dudyC_prim-dudyC_pre) + dw*(dudzC_prim-dudzC_pre)
        this%budget_3(:,:,:,16) = this%budget_3(:,:,:,16) + ubase * buffer
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + du * buffer

        buffer = du*(dvdxC_prim-dvdxC_pre) + dv*(dvdyC_prim-dvdyC_pre) + dw*(dvdzC_prim-dvdzC_pre)
        this%budget_3(:,:,:,16) = this%budget_3(:,:,:,16) + vbase * buffer
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + dv * buffer

        buffer = du*(dwdxC_prim-dwdxC_pre) + dv*(dwdyC_prim-dwdyC_pre) + dw*(dwdzC_prim-dwdzC_pre)
        this%budget_3(:,:,:,16) = this%budget_3(:,:,:,16) + wbase * buffer
        this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + dw * buffer

        if (this%useWindTurbines)then
            cbuffyC1 = this%prim_budget%uturb - this%pre_budget%uturb
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer)
            this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) + du * buffer 
            this%budget_3(:,:,:,21) = this%budget_3(:,:,:,21) + ubase * buffer

            cbuffyC1 = this%prim_budget%vturb - this%pre_budget%vturb
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer)
            this%budget_3(:,:,:,20) = this%budget_3(:,:,:,20) + dv * buffer 
            this%budget_3(:,:,:,21) = this%budget_3(:,:,:,21) + vbase * buffer
        end if 

        nullify(du, dv, dw, rbuffxE1, rbuffxE2, buffer, buffer, cbuffyE1, cbuffyC1, ubase, vbase, wbase)  
        nullify(dudxC_prim, dudyC_prim, dudzC_prim, dudxC_pre, dudyC_pre, dudzC_pre)
        nullify(dvdxC_prim, dvdyC_prim, dvdzC_prim, dvdxC_pre, dvdyC_pre, dvdzC_pre)
        nullify(dwdxC_prim, dwdyC_prim, dwdzC_prim, dwdxC_pre, dwdyC_pre, dwdzC_pre)      
    end subroutine

    subroutine getProductOfMeans(this, budgetid, idx, buffer)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        integer, intent(in) :: idx, budgetid
        real(rkind), dimension(:,:,:), intent(out) :: buffer
        real(rkind), dimension(:,:,:), pointer :: bf, bf2
        
        ! Cell x-pencil buffers 
        bf => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        bf2 => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
        buffer = 0.d0

        if(budgetid.eq.1)then
            select case(idx)
            case(1)
                buffer = this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1)
            case(2)
                buffer = this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2)
            case(3)
                buffer = this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3)
            case(4)
                buffer = this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2)
            case(5)
                buffer = this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3)
            case(6)
                buffer = this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3)
            case(7)
                buffer = this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)
            case(8)
                buffer = this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)
            case(9)
                buffer = this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1)
            case(10)
                buffer = this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)
            case(11)
                buffer = this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1)
            case(12)
                buffer = this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2)
            case(13)
                buffer = this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3)
            case(14)
                buffer = this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2)
            case(15)
                buffer = this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3)
            end select
        
        else if(budgetid.eq.2)then
            select case(idx)
            case(1)
                buffer = this%budget_0(:,:,:,1)*this%MCG(:,:,:,1) + &
                         this%budget_0(:,:,:,2)*this%MCG(:,:,:,2) + &
                         this%budget_0(:,:,:,3)*this%MCG(:,:,:,3)
            case(2)
                buffer = this%budget_0(:,:,:,1)*this%MCG(:,:,:,4) + &
                         this%budget_0(:,:,:,2)*this%MCG(:,:,:,5) + &
                         this%budget_0(:,:,:,3)*this%MCG(:,:,:,6)
            case(3)
                buffer = this%budget_0(:,:,:,1)*this%MCG(:,:,:,7) + &
                         this%budget_0(:,:,:,2)*this%MCG(:,:,:,8) + &
                         this%budget_0(:,:,:,3)*this%MCG(:,:,:,9)
            case(4)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,1) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,2) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,3)
            case(5)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,4) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,5) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,6)
            case(6)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,7) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,8) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,9)
            case(7) 
                buffer = this%budget_0(:,:,:,1)*this%MCG(:,:,:,10) + &
                         this%budget_0(:,:,:,2)*this%MCG(:,:,:,11) + &
                         this%budget_0(:,:,:,3)*this%MCG(:,:,:,12)
            case(8)
                buffer = this%budget_0(:,:,:,1)*this%MCG(:,:,:,13) + &
                         this%budget_0(:,:,:,2)*this%MCG(:,:,:,14) + &
                         this%budget_0(:,:,:,3)*this%MCG(:,:,:,15)
            case(9)
                buffer = this%budget_0(:,:,:,1)*this%MCG(:,:,:,16) + &
                         this%budget_0(:,:,:,2)*this%MCG(:,:,:,17) + &
                         this%budget_0(:,:,:,3)*this%MCG(:,:,:,18)
            case(10)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,10) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,11) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,12)
            case(11)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,13) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,14) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,15)
            case(12)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,16) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,17) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,18)
            case(13)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,1) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,4) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,7)
            case(14)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,2) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,5) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,8)
            case(15)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%MCG(:,:,:,3) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%MCG(:,:,:,6) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%MCG(:,:,:,9)
            end select

        else if(budgetid.eq.3)then
            select case(idx)
            case(1) ! delta u_j' d_j(delta p')
                buffer = this%budget_0(:,:,:,1)*this%budget_0(:,:,:,18) + &
                         this%budget_0(:,:,:,2)*this%budget_0(:,:,:,19) + &
                         this%budget_0(:,:,:,3)*this%budget_0(:,:,:,20)
            
            case(2) ! base  u_j' d_j(delta p')
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,18) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,19) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,20)

            case(3) ! delta u_j' d_j(base p')
                ! px, py, pz signs are reversed in base-flow budget
                buffer = - this%budget_0(:,:,:,1)*this%pre_budget%budget_1(:,:,:,2) &
                         - this%budget_0(:,:,:,2)*this%pre_budget%budget_1(:,:,:,6) &
                         - this%budget_0(:,:,:,3)*this%pre_budget%budget_1(:,:,:,9)

            case(4) ! d_j(base u_i' delta tau_ij') [SGS transport]
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,12) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,13) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,14) + &
                         this%MCG(:,:,:,10) * this%budget_0(:,:,:,6)               + &
                         this%MCG(:,:,:,11) * this%budget_0(:,:,:,7)               + &
                         this%MCG(:,:,:,12) * this%budget_0(:,:,:,8)               + &
                         this%MCG(:,:,:,13) * this%budget_0(:,:,:,7)               + &
                         this%MCG(:,:,:,14) * this%budget_0(:,:,:,9)               + &
                         this%MCG(:,:,:,15) * this%budget_0(:,:,:,10)              + &
                         this%MCG(:,:,:,16) * this%budget_0(:,:,:,8)               + &
                         this%MCG(:,:,:,17) * this%budget_0(:,:,:,10)              + &
                         this%MCG(:,:,:,18) * this%budget_0(:,:,:,11)                

            case(5) ! d_j(delta u_i' base tau_ij') [SGS transport]                
                ! The sign of ui_sgs in this%pre_budget%budget_1 is reversed
                buffer = - this%budget_0(:,:,:,1)*this%pre_budget%budget_1(:,:,:,3)    &
                         - this%budget_0(:,:,:,2)*this%pre_budget%budget_1(:,:,:,7)    &
                         - this%budget_0(:,:,:,3)*this%pre_budget%budget_1(:,:,:,10) + &
                         this%MCG(:,:,:,1) * this%pre_budget%budget_0(:,:,:,11)      + &
                         this%MCG(:,:,:,2) * this%pre_budget%budget_0(:,:,:,12)      + &
                         this%MCG(:,:,:,3) * this%pre_budget%budget_0(:,:,:,13)      + &
                         this%MCG(:,:,:,4) * this%pre_budget%budget_0(:,:,:,12)      + &
                         this%MCG(:,:,:,5) * this%pre_budget%budget_0(:,:,:,14)      + &
                         this%MCG(:,:,:,6) * this%pre_budget%budget_0(:,:,:,15)      + &
                         this%MCG(:,:,:,7) * this%pre_budget%budget_0(:,:,:,13)      + &
                         this%MCG(:,:,:,8) * this%pre_budget%budget_0(:,:,:,15)      + &
                         this%MCG(:,:,:,9) * this%pre_budget%budget_0(:,:,:,16)
     
            case(6) ! d_j(delta u_i' * delta tau_ij')  [SGS transport]
                buffer = this%budget_0(:,:,:,1)*this%budget_0(:,:,:,12)  + &
                         this%budget_0(:,:,:,2)*this%budget_0(:,:,:,13)  + &
                         this%budget_0(:,:,:,3)*this%budget_0(:,:,:,14)  + &
                         this%MCG(:,:,:,1) * this%budget_0(:,:,:,6)      + &
                         this%MCG(:,:,:,2) * this%budget_0(:,:,:,7)      + &
                         this%MCG(:,:,:,3) * this%budget_0(:,:,:,8)      + &
                         this%MCG(:,:,:,4) * this%budget_0(:,:,:,7)      + &
                         this%MCG(:,:,:,5) * this%budget_0(:,:,:,9)      + &
                         this%MCG(:,:,:,6) * this%budget_0(:,:,:,10)     + &
                         this%MCG(:,:,:,7) * this%budget_0(:,:,:,8)      + &
                         this%MCG(:,:,:,8) * this%budget_0(:,:,:,10)     + &
                         this%MCG(:,:,:,9) * this%budget_0(:,:,:,11) 

            case(7) ! delta tau_ij' * d_j(base u_i')     [SGS dissipation]
                buffer = this%MCG(:,:,:,10) * this%budget_0(:,:,:,6)                + &
                         this%MCG(:,:,:,11) * this%budget_0(:,:,:,7)                + &
                         this%MCG(:,:,:,12) * this%budget_0(:,:,:,8)                + &
                         this%MCG(:,:,:,13) * this%budget_0(:,:,:,7)                + &
                         this%MCG(:,:,:,14) * this%budget_0(:,:,:,9)                + &
                         this%MCG(:,:,:,15) * this%budget_0(:,:,:,10)               + &
                         this%MCG(:,:,:,16) * this%budget_0(:,:,:,8)                + &
                         this%MCG(:,:,:,17) * this%budget_0(:,:,:,10)               + &
                         this%MCG(:,:,:,18) * this%budget_0(:,:,:,11)

            case(8) ! base  tau_ij' * d_j(delta u_i')     [SGS dissipation]                
                buffer = this%MCG(:,:,:,1) * this%pre_budget%budget_0(:,:,:,11)      + &
                         this%MCG(:,:,:,2) * this%pre_budget%budget_0(:,:,:,12)      + &
                         this%MCG(:,:,:,3) * this%pre_budget%budget_0(:,:,:,13)      + &
                         this%MCG(:,:,:,4) * this%pre_budget%budget_0(:,:,:,12)      + &
                         this%MCG(:,:,:,5) * this%pre_budget%budget_0(:,:,:,14)      + &
                         this%MCG(:,:,:,6) * this%pre_budget%budget_0(:,:,:,15)      + &
                         this%MCG(:,:,:,7) * this%pre_budget%budget_0(:,:,:,13)      + &
                         this%MCG(:,:,:,8) * this%pre_budget%budget_0(:,:,:,15)      + &
                         this%MCG(:,:,:,9) * this%pre_budget%budget_0(:,:,:,16)

            case(9) ! delta tau_ij' * d_j(delta u_i')     [SGS dissipation]
                buffer = this%MCG(:,:,:,1) * this%budget_0(:,:,:,6)      + &
                         this%MCG(:,:,:,2) * this%budget_0(:,:,:,7)      + &
                         this%MCG(:,:,:,3) * this%budget_0(:,:,:,8)      + &
                         this%MCG(:,:,:,4) * this%budget_0(:,:,:,7)      + &
                         this%MCG(:,:,:,5) * this%budget_0(:,:,:,9)      + &
                         this%MCG(:,:,:,6) * this%budget_0(:,:,:,10)     + &
                         this%MCG(:,:,:,7) * this%budget_0(:,:,:,8)      + &
                         this%MCG(:,:,:,8) * this%budget_0(:,:,:,10)     + &
                         this%MCG(:,:,:,9) * this%budget_0(:,:,:,11) 

            case(10) ! delta u_3' delta wb'
                buffer = this%budget_0(:,:,:,3)*this%budget_0(:,:,:,17)

            case(11) ! delta u_3' base wb'
                buffer = this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,31)

            case(12) ! base u_3' delta wb'
                buffer = this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,17)

            case(13) ! base  u_i' delta u_j' d_j(base u_i')  [Turbulent transport of TKE]
                ! Differentiate mean(base u_i base u_i) numerically
                ! (base u_i * base u_i) is even at the boundaries, so use a flag of 1 at bottom and top
                bf = half*(this%pre_budget%budget_0(:,:,:,4) + this%pre_budget%budget_0(:,:,:,7) + this%pre_budget%budget_0(:,:,:,9))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,2)*bf2
                call this%ddz_R2R(bf, bf2, 1, 1); buffer = buffer + this%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*this%budget_2(:,:,:,4) + &
                                  this%pre_budget%budget_0(:,:,:,2)*this%budget_2(:,:,:,5) + &
                                  this%pre_budget%budget_0(:,:,:,3)*this%budget_2(:,:,:,6) + &
                    this%MCG(:,:,:,10) * (this%budget_1(:,:,:,7)  - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,11) * (this%budget_1(:,:,:,9)  - two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,12) * (this%budget_1(:,:,:,11) - two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,13) * (this%budget_1(:,:,:,8)  - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,14) * (this%budget_1(:,:,:,12) - two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,15) * (this%budget_1(:,:,:,14) - two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,16) * (this%budget_1(:,:,:,10) - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,17) * (this%budget_1(:,:,:,13) - two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,18) * (this%budget_1(:,:,:,15) - two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3))
                
            case(14) ! base  u_i' base u_j' d_j(delta u_i')  [Turbulent transport of TKE]
                buffer = this%pre_budget%budget_0(:,:,:,1)*(this%budget_2(:,:,:,13) + this%budget_2(:,:,:,7)) + &
                         this%pre_budget%budget_0(:,:,:,2)*(this%budget_2(:,:,:,14) + this%budget_2(:,:,:,8)) + &
                         this%pre_budget%budget_0(:,:,:,3)*(this%budget_2(:,:,:,15) + this%budget_2(:,:,:,9)) + &
                         this%MCG(:,:,:,1)*(this%pre_budget%budget_0(:,:,:,4) - two*this%pre_budget%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)) + &
                         this%MCG(:,:,:,2)*(this%pre_budget%budget_0(:,:,:,5) - two*this%pre_budget%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)) + &
                         this%MCG(:,:,:,3)*(this%pre_budget%budget_0(:,:,:,6) - two*this%pre_budget%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)) + &
                         this%MCG(:,:,:,4)*(this%pre_budget%budget_0(:,:,:,5) - two*this%pre_budget%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1)) + &
                         this%MCG(:,:,:,5)*(this%pre_budget%budget_0(:,:,:,7) - two*this%pre_budget%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2)) + &
                         this%MCG(:,:,:,6)*(this%pre_budget%budget_0(:,:,:,8) - two*this%pre_budget%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3)) + &
                         this%MCG(:,:,:,7)*(this%pre_budget%budget_0(:,:,:,6) - two*this%pre_budget%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1)) + &
                         this%MCG(:,:,:,8)*(this%pre_budget%budget_0(:,:,:,8) - two*this%pre_budget%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2)) + &
                         this%MCG(:,:,:,9)*(this%pre_budget%budget_0(:,:,:,9) - two*this%pre_budget%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3))

            case(15) ! delta u_i' base u_j' d_j(base u_i')  [Turbulent transport of TKE]
                bf = this%budget_1(:,:,:,7) + this%budget_1(:,:,:,12) + this%budget_1(:,:,:,15)
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*(bf2 - this%budget_2(:,:,:,13))
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,2)*(bf2 - this%budget_2(:,:,:,14))
                ! bf is an even function. Use a flag of 1 for ddz at both top and bottom
                call this%ddz_R2R(bf, bf2, 1, 1); buffer = buffer + this%pre_budget%budget_0(:,:,:,3)*(bf2 - this%budget_2(:,:,:,15))

                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_2(:,:,:,10) + &
                                  this%budget_0(:,:,:,2)*this%budget_2(:,:,:,11) + &
                                  this%budget_0(:,:,:,3)*this%budget_2(:,:,:,12) + &
                        this%MCG(:,:,:,10) * (this%budget_1(:,:,:,7) - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)) + &
                        this%MCG(:,:,:,11) * (this%budget_1(:,:,:,8) - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)) + &
                        this%MCG(:,:,:,12) * (this%budget_1(:,:,:,10)- two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)) + &
                        this%MCG(:,:,:,13) * (this%budget_1(:,:,:,9) - two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1)) + &
                        this%MCG(:,:,:,14) * (this%budget_1(:,:,:,12)- two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2)) + &
                        this%MCG(:,:,:,15) * (this%budget_1(:,:,:,13)- two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3)) + &
                        this%MCG(:,:,:,16) * (this%budget_1(:,:,:,11)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1)) + &
                        this%MCG(:,:,:,17) * (this%budget_1(:,:,:,14)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2)) + &
                        this%MCG(:,:,:,18) * (this%budget_1(:,:,:,15)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3))
            
            case(16) ! base  u_i' delta u_j' d_j(delta u_i')  [Turbulent transport of TKE]
                buffer = this%budget_0(:,:,:,1)*this%budget_2(:,:,:,13)  + &
                         this%budget_0(:,:,:,2)*this%budget_2(:,:,:,14)  + &
                         this%budget_0(:,:,:,3)*this%budget_2(:,:,:,15)  + &
                         this%pre_budget%budget_0(:,:,:,1)*this%budget_2(:,:,:,1) + &
                         this%pre_budget%budget_0(:,:,:,2)*this%budget_2(:,:,:,2) + &
                         this%pre_budget%budget_0(:,:,:,3)*this%budget_2(:,:,:,3) + &
                    this%MCG(:,:,:,1)*(this%budget_1(:,:,:,7) - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,2)*(this%budget_1(:,:,:,9) - two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,3)*(this%budget_1(:,:,:,11)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,4)*(this%budget_1(:,:,:,8) - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,5)*(this%budget_1(:,:,:,12)- two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,6)*(this%budget_1(:,:,:,14)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,7)*(this%budget_1(:,:,:,10)- two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,8)*(this%budget_1(:,:,:,13)- two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,9)*(this%budget_1(:,:,:,15)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3))
                            
            case(17) ! delta u_i' delta u_j' d_j(base u_i')  [Turbulent transport of TKE]
                ! Differentiate mean(base u_i delta u_i) numerically
                bf = this%budget_1(:,:,:,7) + this%budget_1(:,:,:,12) + this%budget_1(:,:,:,15)
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,1)*(bf2 - this%budget_2(:,:,:,13) + this%budget_2(:,:,:,4))
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,2)*(bf2 - this%budget_2(:,:,:,14) + this%budget_2(:,:,:,5))
                ! bf is an even function. Use a flag of 1 for ddz at both top and bottom
                call this%ddz_R2R(bf, bf2, 1, 1); buffer = buffer + this%budget_0(:,:,:,3)*(bf2 - this%budget_2(:,:,:,15) + this%budget_2(:,:,:,6))
                buffer = buffer                                                                                          + &
                         this%MCG(:,:,:,10)*(this%budget_1(:,:,:,1) - two*this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1)) + &
                         this%MCG(:,:,:,11)*(this%budget_1(:,:,:,2) - two*this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2)) + &
                         this%MCG(:,:,:,12)*(this%budget_1(:,:,:,3) - two*this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3)) + &
                         this%MCG(:,:,:,13)*(this%budget_1(:,:,:,2) - two*this%budget_0(:,:,:,2)*this%budget_0(:,:,:,1)) + &
                         this%MCG(:,:,:,14)*(this%budget_1(:,:,:,4) - two*this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2)) + &
                         this%MCG(:,:,:,15)*(this%budget_1(:,:,:,5) - two*this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3)) + &
                         this%MCG(:,:,:,16)*(this%budget_1(:,:,:,3) - two*this%budget_0(:,:,:,3)*this%budget_0(:,:,:,1)) + &
                         this%MCG(:,:,:,17)*(this%budget_1(:,:,:,5) - two*this%budget_0(:,:,:,3)*this%budget_0(:,:,:,2)) + &
                         this%MCG(:,:,:,18)*(this%budget_1(:,:,:,6) - two*this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3))
            case(18) ! delta u_i' base u_j' d_j(delta u_i')  [Turbulent transport of TKE]
                ! Differentiate mean(delta u_i delta u_i) numerically
                bf = half*(this%budget_1(:,:,:,1) + this%budget_1(:,:,:,4) + this%budget_1(:,:,:,6))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,2)*bf2
                ! bf is an even function. Use a flag of 1 for ddz at both top and bottom
                call this%ddz_R2R(bf, bf2, 1, 1); buffer = buffer + this%pre_budget%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_2(:,:,:,7) + &
                                  this%budget_0(:,:,:,2)*this%budget_2(:,:,:,8) + &
                                  this%budget_0(:,:,:,3)*this%budget_2(:,:,:,9) + &
                    this%MCG(:,:,:,1) * (this%budget_1(:,:,:,7) - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,2) * (this%budget_1(:,:,:,8) - two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,3) * (this%budget_1(:,:,:,10)- two*this%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,4) * (this%budget_1(:,:,:,9) - two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,5) * (this%budget_1(:,:,:,12)- two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,6) * (this%budget_1(:,:,:,13)- two*this%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,7) * (this%budget_1(:,:,:,11)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,8) * (this%budget_1(:,:,:,14)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,9) * (this%budget_1(:,:,:,15)- two*this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3))

            case(19) ! delta u_i' delta u_j' d_j(delta u_i')  [Turbulent transport of TKE]
                ! Differentiate mean(delta u_i delta u_i) numerically
                bf = half*(this%budget_1(:,:,:,1) + this%budget_1(:,:,:,4) + this%budget_1(:,:,:,6))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,2)*bf2
                ! bf is an even function. Use a flag of 1 for ddz at both top and bottom
                call this%ddz_R2R(bf, bf2, 1, 1); buffer = buffer + this%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_2(:,:,:,1) + &
                                  this%budget_0(:,:,:,2)*this%budget_2(:,:,:,2) + &
                                  this%budget_0(:,:,:,3)*this%budget_2(:,:,:,3) + &
                    this%MCG(:,:,:,1) * (this%budget_1(:,:,:,1) - two*this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,2) * (this%budget_1(:,:,:,2) - two*this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,3) * (this%budget_1(:,:,:,3) - two*this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,4) * (this%budget_1(:,:,:,2) - two*this%budget_0(:,:,:,2)*this%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,5) * (this%budget_1(:,:,:,4) - two*this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,6) * (this%budget_1(:,:,:,5) - two*this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3)) + &
                    this%MCG(:,:,:,7) * (this%budget_1(:,:,:,3) - two*this%budget_0(:,:,:,3)*this%budget_0(:,:,:,1)) + &
                    this%MCG(:,:,:,8) * (this%budget_1(:,:,:,5) - two*this%budget_0(:,:,:,3)*this%budget_0(:,:,:,2)) + &
                    this%MCG(:,:,:,9) * (this%budget_1(:,:,:,6) - two*this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3))

            case(20)
                buffer = this%budget_0(:,:,:,1)*this%budget_0(:,:,:,21) + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,22)
            case(21)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,21) + this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,22)
            end select
        end if

        ! Dealias the product of means
        call this%dealias(buffer)
        
        ! Nullify pointers
        nullify(bf, bf2)
    end subroutine

    ! ----------------------supporting subroutines ------------------------
    ! subroutine writeTimeSum(this)
    !     class(budgets_time_avg_deficit_compact), intent(inout), target :: this
    !     character(len=clen) :: fname, tempname 
    !     integer :: ios

    !     write(tempname,"(A3,I2.2,A14,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_time_weight_t",this%prim_budget%igrid_sim%step,"_n",this%counter,".txt"
    !     fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
    !     open(unit=10, file=trim(fname), status='replace', action='write', form='formatted', iostat=ios)
    !     write(10,'(ES23.15)') this%timeSum
    !     close(10)        
    ! end subroutine

    ! subroutine readTimeSum(this, dir, rid, tid, cid)
    !     class(budgets_time_avg_deficit_compact), intent(inout), target :: this
    !     integer, intent(in) :: rid, cid, tid
    !     character(len=clen) :: dir
    !     character(len=clen) :: fname, tempname 
    !     integer :: ios

    !     write(tempname,"(A3,I2.2,A14,I6.6,A2,I6.6,A4)") "Run",rid,"_time_weight_t",tid,"_n",cid,".txt"
    !     fname = trim(dir)//"/"//trim(tempname)
    !     open(unit=10, file=trim(fname), status='old', action='read', form='formatted', iostat=ios)
    !     read(10,'(ES23.15)') this%timeSum
    !     close(10)        
    ! end subroutine

    subroutine dump_budget_field(this, field, fieldID, BudgetID)
        use decomp_2d_io
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: field
        integer, intent(in) :: fieldID, BudgetID
        character(len=clen) :: fname, tempname 

        write(tempname,"(A3,I2.2,A20,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_comp_deficit_budget",BudgetID,"_term",fieldID,"_t",this%prim_budget%igrid_sim%step,"_n",this%counter,".s3D"
        fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)

        call decomp_2d_write_one(1,field,fname, this%prim_budget%igrid_sim%gpC)
    end subroutine 

    subroutine restart_budget_field(this, field, dir, runID, timeID, counterID, budgetID, fieldID)
        use decomp_2d_io
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: field
        integer, intent(in) :: runID, counterID, timeID, budgetID, fieldID
        character(len=clen) :: fname, tempname
        character(len=clen), intent(in) :: dir

        write(tempname,"(A3,I2.2,A20,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A4)") "Run",runID,"_comp_deficit_budget",budgetID,"_term",fieldID,"_t",timeID,"_n",counterID,".s3D"
        fname = dir(:len_trim(dir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,field,fname, this%prim_budget%igrid_sim%gpC)           
     end subroutine 

     subroutine RestartBudget(this, dir, rid, tid, cid)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        integer, intent(in) :: rid, cid, tid
        character(len=clen) :: dir
        integer :: idx
        real(rkind), dimension(:,:,:), pointer :: buffer
        real(rkind) :: totalWeight

        ! Cell x-pencil buffers 
        buffer => this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)
        this%counter = cid

        ! if(this%time_weighted_average)then
        !     ! If this is time-weighted averaging, we should read the sum of times
        !     call this%readTimeSum(trim(dir),rid,tid,cid)
        !     totalWeight = this%timeSum + 1.d-18
        ! else
        !     totalWeight = real(this%counter,rkind) + 1.d-18
        ! end if        
        totalWeight = real(this%counter,rkind) + 1.d-18

        ! I assume here that this%pre_budget%budget_0 and 
        ! this%pre_budget%budget_1 are already restarted 
        ! and are in summing mode
        this%pre_budget%budget_0 = this%pre_budget%budget_0/totalWeight
        this%pre_budget%budget_1 = this%pre_budget%budget_1/totalWeight

        ! Budget 0 
        if(this%do_budget0)then
            do idx = 1, this%size_budget_0
                if((idx.eq.15).or.(idx.eq.16))then
                    if(.not. this%useCoriolis)cycle
                end if
                if((idx.eq.5).or.(idx.eq.17))then
                    if(.not. this%isStratified)cycle
                end if

                call this%restart_budget_field(this%budget_0(:,:,:,idx), dir, rid, tid, cid, 0, idx)
            end do
        end if

        ! Budget 1
        if(this%do_budget1)then
            do idx = 1, this%size_budget_1
                call this%restart_budget_field(this%budget_1(:,:,:,idx), dir, rid, tid, cid, 1, idx)
                call this%getProductOfMeans(1, idx, buffer)
                this%budget_1(:,:,:,idx) = this%budget_1(:,:,:,idx) + buffer
            end do
        end if

        ! Budget 2
        if(this%do_budget2)then
            do idx = 1, this%size_budget_2
                call this%restart_budget_field(this%budget_2(:,:,:,idx), dir, rid, tid, cid, 2, idx)
                call this%getProductOfMeans(2, idx, buffer)
                this%budget_2(:,:,:,idx) = this%budget_2(:,:,:,idx) + buffer
            end do
        end if

        ! Budget 3
        if(this%do_budget3)then
            do idx = 1, this%size_budget_3
                if((idx.eq.10).or.(idx.eq.11).or.(idx.eq.12))then
                    if(.not. this%isStratified) cycle
                end if

                call this%restart_budget_field(this%budget_3(:,:,:,idx), dir, rid, tid, cid, 3, idx)
                call this%getProductOfMeans(3, idx, buffer)
                this%budget_3(:,:,:,idx) = this%budget_3(:,:,:,idx) + buffer
            end do
        end if

        ! Return to summing
        if(this%do_budget0) this%budget_0 = this%budget_0*totalWeight
        if(this%do_budget1) this%budget_1 = this%budget_1*totalWeight
        if(this%do_budget2) this%budget_2 = this%budget_2*totalWeight
        if(this%do_budget3) this%budget_3 = this%budget_3*totalWeight
        this%pre_budget%budget_0 = this%pre_budget%budget_0*totalWeight
        this%pre_budget%budget_1 = this%pre_budget%budget_1*totalWeight

        ! To save time and storage, MCG were not written to file.
        ! We restart MCG by numerically differentiating the mean flow
        ! MCG is automatically in the summing mode because we
        ! differentiate budget 0 in the summing mode
        if(this%doMCG) call this%restartMCG()

        nullify(buffer)         
    end subroutine

    subroutine restartMCG(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: dudx_def, dudy_def, dudz_def, dudx_pre, dudy_pre, dudz_pre
        real(rkind), dimension(:,:,:), pointer :: dvdx_def, dvdy_def, dvdz_def, dvdx_pre, dvdy_pre, dvdz_pre
        real(rkind), dimension(:,:,:), pointer :: dwdx_def, dwdy_def, dwdz_def, dwdx_pre, dwdy_pre, dwdz_pre

        dudx_def => this%MCG(:,:,:,1)
        dudy_def => this%MCG(:,:,:,2)
        dudz_def => this%MCG(:,:,:,3)
        dvdx_def => this%MCG(:,:,:,4)
        dvdy_def => this%MCG(:,:,:,5)
        dvdz_def => this%MCG(:,:,:,6)
        dwdx_def => this%MCG(:,:,:,7)
        dwdy_def => this%MCG(:,:,:,8)
        dwdz_def => this%MCG(:,:,:,9)
        dudx_pre => this%MCG(:,:,:,10)
        dudy_pre => this%MCG(:,:,:,11)
        dudz_pre => this%MCG(:,:,:,12)
        dvdx_pre => this%MCG(:,:,:,13)
        dvdy_pre => this%MCG(:,:,:,14)
        dvdz_pre => this%MCG(:,:,:,15)
        dwdx_pre => this%MCG(:,:,:,16)
        dwdy_pre => this%MCG(:,:,:,17)
        dwdz_pre => this%MCG(:,:,:,18)

        call this%ddx_R2R(this%budget_0(:,:,:,1), dudx_def)
        call this%ddy_R2R(this%budget_0(:,:,:,1), dudy_def)
        call this%ddz_R2R(this%budget_0(:,:,:,1), dudz_def, uBC_bottom, uBC_top)
        call this%ddx_R2R(this%budget_0(:,:,:,2), dvdx_def)
        call this%ddy_R2R(this%budget_0(:,:,:,2), dvdy_def)
        call this%ddz_R2R(this%budget_0(:,:,:,2), dvdz_def, vBC_bottom, vBC_top)
        call this%ddx_R2R(this%budget_0(:,:,:,3), dwdx_def)
        call this%ddy_R2R(this%budget_0(:,:,:,3), dwdy_def)
        call this%ddz_R2R(this%budget_0(:,:,:,3), dwdz_def, wBC_bottom, wBC_top)
        call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1), dudx_pre)
        call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1), dudy_pre)
        call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1), dudz_pre, uBC_bottom, uBC_top)
        call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2), dvdx_pre)
        call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2), dvdy_pre)
        call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2), dvdz_pre, vBC_bottom, vBC_top)
        call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3), dwdx_pre)
        call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3), dwdy_pre)
        call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3), dwdz_pre, wBC_bottom, wBC_top)

        nullify(dudx_def, dudy_def, dudz_def, dudx_pre, dudy_pre, dudz_pre)
        nullify(dvdx_def, dvdy_def, dvdz_def, dvdx_pre, dvdy_pre, dvdz_pre)
        nullify(dwdx_def, dwdy_def, dwdz_def, dwdx_pre, dwdy_pre, dwdz_pre)
    end subroutine

    subroutine ResetBudget(this)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        
        this%counter = 0
        this%timeSum = zero
        if(allocated(this%budget_0)) this%budget_0 = zero
        if(allocated(this%budget_1)) this%budget_1 = zero 
        if(allocated(this%budget_2)) this%budget_2 = zero 
        if(allocated(this%budget_3)) this%budget_3 = zero
        if(allocated(this%delta_tauij)) this%delta_tauij = zero 
        if(allocated(this%MCG)) this%MCG = zero   
    end subroutine 

    subroutine destroy(this)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this

        nullify(this%prim_budget, this%pre_budget)
        if(this%do_budgets) then
            if(allocated(this%budget_0)) deallocate(this%budget_0)
            if(allocated(this%budget_1)) deallocate(this%budget_1)
            if(allocated(this%budget_2)) deallocate(this%budget_2)
            if(allocated(this%budget_3)) deallocate(this%budget_3)
            if(allocated(this%delta_tauij)) deallocate(this%delta_tauij)
            if(allocated(this%MCG)) deallocate(this%MCG)
        end if
    end subroutine 

    ! ----------------------private derivative operators ------------------------
    subroutine dealias(this, f)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(inout) :: f
        complex(rkind), dimension(:,:,:), pointer :: cbuffyC

        cbuffyC => this%prim_budget%igrid_sim%cbuffyC(:,:,:,1)
        
        call this%prim_budget%igrid_sim%spectC%fft(f, cbuffyC)
        call this%prim_budget%igrid_sim%spectC%dealias(cbuffyC)
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC, f)
    end subroutine

    subroutine ddx_R2R(this, f, dfdx)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdx
        complex(rkind), dimension(:,:,:), pointer :: cbuffyC

        cbuffyC => this%prim_budget%igrid_sim%cbuffyC(:,:,:,1)
        
        call this%prim_budget%igrid_sim%spectC%fft(f, cbuffyC)
        call this%prim_budget%igrid_sim%spectC%mtimes_ik1_ip(cbuffyC)
        call this%prim_budget%igrid_sim%spectC%dealias(cbuffyC)
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC, dfdx)

        nullify(cbuffyC)
    end subroutine 

    subroutine ddy_R2R(this, f, dfdy)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdy
        complex(rkind), dimension(:,:,:), pointer :: cbuffyC

        cbuffyC => this%prim_budget%igrid_sim%cbuffyC(:,:,:,1)
        
        call this%prim_budget%igrid_sim%spectC%fft(f, cbuffyC)
        call this%prim_budget%igrid_sim%spectC%mtimes_ik2_ip(cbuffyC)
        call this%prim_budget%igrid_sim%spectC%dealias(cbuffyC)
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC, dfdy)

        nullify(cbuffyC)
    end subroutine 
     
    subroutine ddz_R2R(this, f, dfdz, n1, n2)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdz
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(:,:,:), pointer :: cbuffyC, cbuffzC1, cbuffzC2

        cbuffyC => this%prim_budget%igrid_sim%cbuffyC(:,:,:,1)
        cbuffzC1 => this%prim_budget%igrid_sim%cbuffzC(:,:,:,1)
        cbuffzC2 => this%prim_budget%igrid_sim%cbuffzC(:,:,:,2)

        call this%prim_budget%igrid_sim%spectC%fft(f, cbuffyC)
        call transpose_y_to_z(cbuffyC, cbuffzC1, this%prim_budget%igrid_sim%sp_gpC)
        call this%prim_budget%igrid_sim%Pade6opZ%ddz_C2C(cbuffzC1, cbuffzC2, n1, n2)
        call transpose_z_to_y(cbuffzC2, cbuffyC, this%prim_budget%igrid_sim%sp_gpC)
        call this%prim_budget%igrid_sim%spectC%dealias(cbuffyC)
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC, dfdz)

        nullify(cbuffyC, cbuffzC1, cbuffzC2)
    end subroutine
     
    ! subroutine ddz_C2R(this, fhat, dfdz, n1, n2)
    !     class(budgets_time_avg_deficit_compact), intent(inout) :: this
    !     complex(rkind), dimension(this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(1),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(2),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
    !     real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdz
    !     integer, intent(in) :: n1, n2
        
    !     call transpose_y_to_z(fhat,this%prim_budget%igrid_sim%cbuffzC(:,:,:,1),this%prim_budget%igrid_sim%sp_gpC)
    !     call this%prim_budget%igrid_sim%Pade6opZ%ddz_C2C(this%prim_budget%igrid_sim%cbuffzC(:,:,:,1),this%prim_budget%igrid_sim%cbuffzC(:,:,:,2),n1,n2)
    !     call transpose_z_to_y(this%prim_budget%igrid_sim%cbuffzC(:,:,:,2),this%prim_budget%igrid_sim%cbuffyC(:,:,:,1),this%prim_budget%igrid_sim%sp_gpC)
    !     call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
    !     call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdz)
    ! end subroutine 
 
    subroutine interp_Edge2Cell(this, fE, fC, n1, n2)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(in) :: fE
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: fC
        integer, intent(in) :: n1, n2
        real(rkind), dimension(:,:,:), pointer :: rbuffyE, rbuffzE, rbuffzC, rbuffyC

        rbuffyE => this%prim_budget%igrid_sim%rbuffyE(:,:,:,1)
        rbuffzE => this%prim_budget%igrid_sim%rbuffzE(:,:,:,1)
        rbuffzC => this%prim_budget%igrid_sim%rbuffzC(:,:,:,2)
        rbuffyC => this%prim_budget%igrid_sim%rbuffyC(:,:,:,1)

        call transpose_x_to_y(fE, rbuffyE, this%prim_budget%igrid_sim%gpE)
        call transpose_y_to_z(rbuffyE, rbuffzE, this%prim_budget%igrid_sim%gpE)
        call this%prim_budget%igrid_sim%Pade6opZ%interpz_E2C(rbuffzE, rbuffzC, n1, n2)
        call transpose_z_to_y(rbuffzC, rbuffyC, this%prim_budget%igrid_sim%gpC)
        call transpose_y_to_x(rbuffyC, fC, this%prim_budget%igrid_sim%gpC)

        nullify(rbuffyE, rbuffzE, rbuffzC, rbuffyC)
    end subroutine 
 
    ! subroutine interp_Cell2Edge(this, fC, fE, n1, n2)
    !     class(budgets_time_avg_deficit_compact), intent(inout) :: this
    !     real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: fC
    !     real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(out) :: fE

    !     call transpose_x_to_y(fC,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
    !     call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
    !     call this%prim_budget%igrid_sim%Pade6opZ%interpz_C2E(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),n1,n2)
    !     call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),this%prim_budget%igrid_sim%gpE)
    !     call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),fE,this%prim_budget%igrid_sim%gpE)
    ! end subroutine 
         
    ! subroutine multiply_CellFieldsOnEdges(this, f1C, f2C, fmultC, n1, n2)
    !     class(budgets_time_avg_deficit_compact), intent(inout) :: this
    !     real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f1C,f2C
    !     real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: fmultC
    !     integer, intent(in) :: n1, n2

    !     ! interpolate 1st Cell field
    !     call transpose_x_to_y(f1C,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
    !     call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
    !     call this%prim_budget%igrid_sim%Pade6opZ%interpz_C2E(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),n1,n2)

    !     ! interpolate 2nd Cell field
    !     call transpose_x_to_y(f2C,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
    !     call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
    !     call this%prim_budget%igrid_sim%Pade6opZ%interpz_C2E(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,2),n1,n2)

    !     ! multiply on Edges and interpolate back to Cells
    !     this%prim_budget%igrid_sim%rbuffzE(:,:,:,1) = this%prim_budget%igrid_sim%rbuffzE(:,:,:,1) * this%prim_budget%igrid_sim%rbuffzE(:,:,:,2)
    !     call this%prim_budget%igrid_sim%Pade6opZ%interpz_E2C(this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),n1,n2)
    !     call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
    !     call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),fmultC,this%prim_budget%igrid_sim%gpC)
    ! end subroutine 

    ! multiply on edge cells and interpolate to cell centers to reduce aliasing issues
    ! function multiply_Edges_interp_cell(this, f1E, f2E, n1, n2) result(fmultC)
    !     class(budgets_time_avg_deficit_compact), intent(inout) :: this
    !     real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(in) :: f1E,f2E
    !     real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)) :: fmultC
    !     integer, intent(in) :: n1, n2

    !     call this%interp_Edge2Cell(f1E * f2E, fmultC, n1, n2)
    ! end function
end module
