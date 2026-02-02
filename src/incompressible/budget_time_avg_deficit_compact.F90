module budgets_time_avg_deficit_compact_mod
    use kind_parameters, only: rkind, clen
    use decomp_2d
    use budgets_time_avg_mod, only: budgets_time_avg
    use exits, only: message, GracefulExit
    use constants, only: zero
    use mpi
 
    implicit none 
 
    private
    public :: budgets_time_avg_deficit_compact

    ! Comments here

    type :: budgets_time_avg_deficit_compact
        private
        integer :: run_id, nx, ny, nz
        logical :: do_budget0=.false., do_budget1=.false., do_budget2=.false., do_budget3=.false.
        logical :: write_budget0=.false., write_budget1=.false., write_budget2=.false., write_budget3=.false.
        
        type(budgets_time_avg), pointer :: pre_budget, prim_budget
        
        real(rkind), dimension(:,:,:,:), allocatable :: budget_0, budget_1, budget_2, budget_3
        integer :: size_budget_0, size_budget_1, size_budget_2, size_budget_3
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
   
        procedure, private  :: getProductOfMeans
        ! procedure, private  :: writeTimeSum
        ! procedure, private  :: readTimeSum

        procedure, private :: ddx_R2R
        procedure, private :: ddy_R2R
        procedure, private :: ddz_R2R
        procedure, private :: ddz_C2R
        procedure, private :: dealias
        procedure, private :: interp_Edge2Cell
        procedure, private :: interp_Cell2Edge
        procedure, private :: multiply_CellFieldsOnEdges
        procedure, private :: multiply_edges_interp_cell
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
        logical :: write_budget0=.false., write_budget1=.false., write_budget2=.false., write_budget3=.false.
        namelist /BUDGET_TIME_AVG_DEFICIT_COMPACT/ budgets_dir, restart_budgets, restart_dir, &
            restart_rid, restart_tid, restart_counter, tidx_dump, tidx_compute, do_budgets, &
            use_time_weighted_average, tidx_budget_start, time_budget_start, &
            write_budget0, write_budget1, write_budget2, write_budget3 

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
        this%write_budget0 = write_budget0
        this%write_budget1 = write_budget1
        this%write_budget2 = write_budget2
        this%write_budget3 = write_budget3

        if(write_budget0)this%do_budget0=.true.
        if(write_budget1)this%do_budget1=.true.
        if(write_budget2)this%do_budget2=.true.
        if(write_budget3)this%do_budget3=.true.
        
        if(this%do_budget1)this%do_budget0=.true.
        if(this%do_budget2)this%do_budget0=.true.
        if(this%do_budget3)then
            this%do_budget0=.true.
            this%do_budget1=.true.
            this%do_budget2=.true.
        end if
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
                this%size_budget_2 = 12
                allocate(this%budget_2(this%nx,this%ny,this%nz,this%size_budget_2))
            end if

            if(this%do_budget3)then
                if(this%useWindTurbines)then
                    this%size_budget_3 = 19
                else
                    this%size_budget_3 = 17
                end if
                allocate(this%budget_3(this%nx,this%ny,this%nz,this%size_budget_3))
                allocate(this%delta_tauij(this%nx,this%ny,this%nz,6))
            end if

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

        ! To be multiplied by every term added to the sum
        ! if(this%time_weighted_average)then
        !     this%weight = this%prim_budget%igrid_sim%dt
        ! else
        !     this%weight = real(1., rkind)
        ! end if

        if(this%do_budget0) call this%AssembleBudget0()
        if(this%do_budget1) call this%AssembleBudget1()
        if(this%do_budget2) call this%AssembleBudget2()
        if(this%do_budget3) call this%AssembleBudget3()

        this%counter = this%counter + 1
        ! this%timeSum = this%timeSum + this%prim_budget%igrid_sim%dt
    end subroutine 

    subroutine DumpBudget(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind) :: totalWeight
        integer :: idx, budgetid, budgetsize
        real(rkind), dimension(:,:,:), pointer :: buffer
        real(rkind), dimension(:,:,:,:), pointer :: budget
        logical :: writeBudget

        ! if(this%time_weighted_average)then
        !     totalWeight = this%timeSum + 1.d-18
        !     call this%writeTimeSum()
        ! else
        !     totalWeight = real(this%counter,rkind) + 1.d-18
        ! end if
        totalWeight = real(this%counter,rkind) + 1.d-18

        ! Cell x-pencil buffers 
        ! Buffers 1 and 2 are used locally inside getProductOfMeans
        buffer => this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)

        ! Convert assembled budgets to mean instead of sum
        if(this%do_budget0) this%budget_0 = this%budget_0/totalWeight
        if(this%do_budget1) this%budget_1 = this%budget_1/totalWeight
        if(this%do_budget2) this%budget_2 = this%budget_2/totalWeight
        if(this%do_budget3) this%budget_3 = this%budget_3/totalWeight
        this%pre_budget%budget_0 = this%pre_budget%budget_0/totalWeight
        this%pre_budget%budget_1 = this%pre_budget%budget_1/totalWeight

        ! Budget 0
        if(this%write_budget0)then
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

        do budgetid=1,3
            select case(budgetid)
            case(1)
                budget => this%budget_1
                budgetsize = this%size_budget_1
                writeBudget = this%write_budget1
            case(2)
                budget => this%budget_2
                budgetsize = this%size_budget_2
                writeBudget = this%write_budget2
            case(3)
                budget => this%budget_3
                budgetsize = this%size_budget_3
                writeBudget = this%write_budget3
            end select

            if(writeBudget)then
                do idx = 1,budgetsize

                    ! Skip Buoyancy covariance in TKE budget
                    if(budgetid.eq.3)then
                        if((idx.eq.10).or.(idx.eq.11).or.(idx.eq.12))then
                            if(.not. this%isStratified) cycle
                        end if
                    end if

                    ! Get the product of means
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
        this%pre_budget%budget_0 = this%pre_budget%budget_0*totalWeight
        this%pre_budget%budget_1 = this%pre_budget%budget_1*totalWeight
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
            call this%interp_Edge2Cell(rbuffxE1, rbuffxC1)
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

        cbuffyE1 = this%pre_budget%wsgs - this%prim_budget%wsgs
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, rbuffxC1)
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
        cbuffyE1 = this%pre_budget%pz - this%prim_budget%pz
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, rbuffxC1)
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
        real(rkind), dimension(:,:,:), pointer :: du, dv, dw, duE, dvE, dwE, buffer

        ! Cell x-pencil buffers 
        du =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        dv =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
        dw =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,3)
        buffer =>  this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)
        
        ! Edge x-pencil buffers (only 2 are allocated in igrid.F90)
        duE => this%prim_budget%igrid_sim%rbuffxE(:,:,:,1)
        dvE => this%prim_budget%igrid_sim%rbuffxE(:,:,:,2)
        dwE => this%pre_budget%igrid_sim%rbuffxE(:,:,:,1)

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
        buffer = this%multiply_Edges_interp_cell(duE, dwE)
        this%budget_1(:,:,:,3) = this%budget_1(:,:,:,3) + buffer
        this%budget_1(:,:,:,4) = this%budget_1(:,:,:,4) + dv * dv
        buffer = this%multiply_Edges_interp_cell(dvE, dwE)
        this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) + buffer
        this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) + dw * dw
         
        ! Mixed Reynolds stresses
        this%budget_1(:,:,:,7)  = this%budget_1(:,:,:,7) + du * this%pre_budget%igrid_sim%u
        this%budget_1(:,:,:,8)  = this%budget_1(:,:,:,8) + du * this%pre_budget%igrid_sim%v
        this%budget_1(:,:,:,9)  = this%budget_1(:,:,:,9) + dv * this%pre_budget%igrid_sim%u
        buffer = this%multiply_Edges_interp_cell(duE, this%pre_budget%igrid_sim%w)
        this%budget_1(:,:,:,10) = this%budget_1(:,:,:,10) + buffer
        buffer = this%multiply_Edges_interp_cell(dwE, this%pre_budget%igrid_sim%uE)
        this%budget_1(:,:,:,11) = this%budget_1(:,:,:,11) + buffer
        this%budget_1(:,:,:,12) = this%budget_1(:,:,:,12) + dv * this%pre_budget%igrid_sim%v
        buffer = this%multiply_Edges_interp_cell(dvE, this%pre_budget%igrid_sim%w)
        this%budget_1(:,:,:,13) = this%budget_1(:,:,:,13) + buffer
        buffer = this%multiply_Edges_interp_cell(dwE, this%pre_budget%igrid_sim%vE)
        this%budget_1(:,:,:,14) = this%budget_1(:,:,:,14) + buffer
        this%budget_1(:,:,:,15) = this%budget_1(:,:,:,15) + dw * this%pre_budget%igrid_sim%wC

        nullify(du, dv, dw, duE, dvE, dwE, buffer)
    end subroutine

    ! ---------------------- Budget 2 ------------------------
    subroutine AssembleBudget2(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: du, dv, buffC, dw, buffer
        real(rkind), dimension(:,:,:), pointer :: dwE, buffE, duE, dvE
        real(rkind), dimension(:,:,:), pointer :: ubase, vbase, wbaseE, ubaseE, vbaseE
        real(rkind), dimension(:,:,:), pointer :: dudxC_prim, dudyC_prim, dudzE_prim, dudxC_pre, dudyC_pre, dudzE_pre
        real(rkind), dimension(:,:,:), pointer :: dudzC_prim, dvdzC_prim, dudzC_pre, dvdzC_pre 
        real(rkind), dimension(:,:,:), pointer :: dvdxC_prim, dvdyC_prim, dvdzE_prim, dvdxC_pre, dvdyC_pre, dvdzE_pre
        real(rkind), dimension(:,:,:), pointer :: dwdxE_prim, dwdyE_prim, dwdzE_prim, dwdxE_pre, dwdyE_pre, dwdzE_pre
    
        ! Cell x-pencil buffers 
        du => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        dv => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
        buffC => this%prim_budget%igrid_sim%rbuffxC(:,:,:,3)
        dw => this%prim_budget%igrid_sim%rbuffxC(:,:,:,4)
        buffer => this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)
        dwE => this%prim_budget%igrid_sim%rbuffxE(:,:,:,1)         
        buffE => this%prim_budget%igrid_sim%rbuffxE(:,:,:,2) 
        duE => this%pre_budget%igrid_sim%rbuffxE(:,:,:,1)         
        dvE => this%pre_budget%igrid_sim%rbuffxE(:,:,:,2)     
        
        ! Perturbation fields
        du = this%prim_budget%igrid_sim%u - this%pre_budget%igrid_sim%u
        dv = this%prim_budget%igrid_sim%v - this%pre_budget%igrid_sim%v
        dw = this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC
        duE = this%prim_budget%igrid_sim%uE - this%pre_budget%igrid_sim%uE
        dvE = this%prim_budget%igrid_sim%vE - this%pre_budget%igrid_sim%vE
        dwE = this%prim_budget%igrid_sim%w - this%pre_budget%igrid_sim%w

        ! Base-flow fields
        ubase => this%pre_budget%igrid_sim%u
        vbase => this%pre_budget%igrid_sim%v
        ubaseE => this%pre_budget%igrid_sim%uE
        vbaseE => this%pre_budget%igrid_sim%vE
        wbaseE=> this%pre_budget%igrid_sim%w
        ! -----------------------------------------------------------
        !
        ! Primary simulation:
        ! Cell gradients
        dudxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,1)
        dudyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,2)
        dudzC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,3)
        dvdxC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,4)
        dvdyC_prim => this%prim_budget%igrid_sim%duidxjC(:,:,:,5)

        ! Edge gradients
        dudzE_prim => this%prim_budget%igrid_sim%duidxjE(:,:,:,3)
        dvdzE_prim => this%prim_budget%igrid_sim%duidxjE(:,:,:,6)
        dwdxE_prim => this%prim_budget%igrid_sim%duidxjE(:,:,:,7)
        dwdyE_prim => this%prim_budget%igrid_sim%duidxjE(:,:,:,8)
        dwdzE_prim => this%prim_budget%igrid_sim%duidxjE(:,:,:,9)        
        ! -----------------------------------------------------------
        !
        ! Precursor simulation:
        ! Cell gradients
        dudxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,1)
        dudyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,2)
        dudzC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,3)
        dvdxC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,4)
        dvdyC_pre => this%pre_budget%igrid_sim%duidxjC(:,:,:,5)       

        ! Edge gradients
        dudzE_pre => this%pre_budget%igrid_sim%duidxjE(:,:,:,3)
        dvdzE_pre => this%pre_budget%igrid_sim%duidxjE(:,:,:,6)
        dwdxE_pre => this%pre_budget%igrid_sim%duidxjE(:,:,:,7)
        dwdyE_pre => this%pre_budget%igrid_sim%duidxjE(:,:,:,8)
        dwdzE_pre => this%pre_budget%igrid_sim%duidxjE(:,:,:,9)
        ! -----------------------------------------------------------

        ! Term 1: delta u_j d_j(delta u)
        ! buffE = dwE * (dudzE_prim - dudzE_pre)
        ! call this%interp_Edge2Cell(buffE, buffC)
        buffC = du * (dudxC_prim - dudxC_pre) + dv * (dudyC_prim - dudyC_pre) + dw * (dudzC_prim - dudzC_pre)
        call this%dealias(buffC)
        this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + buffC
        
        ! this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + & 
        !         du * (dudxC_prim - dudxC_pre) + dv * (dudyC_prim - dudyC_pre) + dw * (dudzC_prim - dudzC_pre)
        !buffE = dwE * (dudzE_prim - dudzE_pre)
        ! buffE = dudzE_prim - dudzE_pre
        ! call this%interp_Edge2Cell(buffE, buffC)
        ! this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + buffC * dw

        ! Term 2: delta u_j d_j(delta v)
        this%budget_2(:,:,:,2) = this%budget_2(:,:,:,2) +  du * (dvdxC_prim - dvdxC_pre) + dv * (dvdyC_prim - dvdyC_pre)
        buffE = dwE * (dvdzE_prim - dvdzE_pre)
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,2) = this%budget_2(:,:,:,2) + buffC

        ! Term 3: delta u_j d_j(delta w)
        buffE = duE * (dwdxE_prim - dwdxE_pre) + dvE * (dwdyE_prim - dwdyE_pre) + dwE * (dwdzE_prim - dwdzE_pre)
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,3) = this%budget_2(:,:,:,3) + buffC
        
        ! Term 4: delta u_j d_j(base u)
        this%budget_2(:,:,:,4) = this%budget_2(:,:,:,4) + du * dudxC_pre +  dv * dudyC_pre
        buffE = dwE * dudzE_pre
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,4) = this%budget_2(:,:,:,4) + buffC

        ! Term 5: delta u_j d_j(base v)
        this%budget_2(:,:,:,5) = this%budget_2(:,:,:,5) +  du * dvdxC_pre + dv * dvdyC_pre
        buffE = dwE * dvdzE_pre
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,5) = this%budget_2(:,:,:,5) + buffC

        ! Term 6: delta u_j d_j(base w)
        buffE = duE * dwdxE_pre + dvE * dwdyE_pre + dwE * dwdzE_pre
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + buffC
        
        ! Term 7: base u_j d_j(delta u)
        this%budget_2(:,:,:,7) = this%budget_2(:,:,:,7) + ubase * (dudxC_prim - dudxC_pre) + vbase * (dudyC_prim - dudyC_pre)
        buffE = wbaseE * (dudzE_prim - dudzE_pre)
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,7) = this%budget_2(:,:,:,7) + buffC

        ! Term 8: base u_j d_j(delta v)
        this%budget_2(:,:,:,8) = this%budget_2(:,:,:,8) + ubase * (dvdxC_prim - dvdxC_pre) + vbase * (dvdyC_prim - dvdyC_pre)
        buffE = wbaseE * (dvdzE_prim - dvdzE_pre)
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,8) = this%budget_2(:,:,:,8) + buffC

        ! Term 9: base u_j d_j(delta w)
        buffE = ubaseE * (dwdxE_prim - dwdxE_pre) + vbaseE * (dwdyE_prim-dwdyE_pre) + wbaseE * (dwdzE_prim-dwdzE_pre)
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,9) = this%budget_2(:,:,:,9) + buffC

        ! Term 10: base u_j d_j(base u)
        this%budget_2(:,:,:,10) = this%budget_2(:,:,:,10) + ubase * dudxC_pre + vbase * dudyC_pre
        buffE = wbaseE * dudzE_pre
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,10) = this%budget_2(:,:,:,10) + buffC

        ! Term 11: base u_j d_j(base v)
        this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + ubase * dvdxC_pre +  vbase * dvdyC_pre
        buffE = wbaseE * dvdzE_pre
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,11) = this%budget_2(:,:,:,11) + buffC

        ! Term 12: base u_j d_j(base w)
        buffE=ubaseE * dwdxE_pre + vbaseE * dwdyE_pre + wbaseE * dwdzE_pre        
        call this%interp_Edge2Cell(buffE, buffC)
        this%budget_2(:,:,:,12) = this%budget_2(:,:,:,12) + buffC

        ! Release memory        
        nullify(du, dv, dw, buffC)
        nullify(dwE, buffE, duE, dvE)
        nullify(ubase, vbase, wbaseE, ubaseE, vbaseE)
        nullify(dudxC_prim, dudyC_prim, dudzE_prim, dudxC_pre, dudyC_pre, dudzE_pre)
        nullify(dvdxC_prim, dvdyC_prim, dvdzE_prim, dvdxC_pre, dvdyC_pre, dvdzE_pre)
        nullify(dwdxE_prim, dwdyE_prim, dwdzE_prim, dwdxE_pre, dwdyE_pre, dwdzE_pre)
    end subroutine

    ! ---------------------- Budget 3 ------------------------
    subroutine AssembleBudget3(this)
        class(budgets_time_avg_deficit_compact), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: du, dv, dw
        real(rkind), dimension(:,:,:), pointer :: ubase, vbase, wcbase
        real(rkind), dimension(:,:,:), pointer :: rbuffxE1, buffer, bf
        complex(rkind), dimension(:,:,:), pointer :: cbuffyE1, cbuffyC1

        ! Cell x-pencil buffers 
        du => this%prim_budget%igrid_sim%rbuffxC(:,:,:,1)
        dv => this%prim_budget%igrid_sim%rbuffxC(:,:,:,2)
        dw => this%prim_budget%igrid_sim%rbuffxC(:,:,:,3)

        buffer => this%pre_budget%igrid_sim%rbuffxC(:,:,:,1)
        bf     => this%pre_budget%igrid_sim%rbuffxC(:,:,:,2)

        ! Cell y-pencil buffer 
        cbuffyC1 => this%prim_budget%igrid_sim%cbuffyC(:,:,:,2)

        ! Edge x-pencil buffer
        rbuffxE1 => this%prim_budget%igrid_sim%rbuffxE(:,:,:,1)

        ! Edge y-pencil buffer
        cbuffyE1 => this%prim_budget%igrid_sim%cbuffyE(:,:,:,1)

        ! Perturbation fields
        du = this%prim_budget%igrid_sim%u  - this%pre_budget%igrid_sim%u
        dv = this%prim_budget%igrid_sim%v  - this%pre_budget%igrid_sim%v
        dw = this%prim_budget%igrid_sim%wC - this%pre_budget%igrid_sim%wC

        ubase => this%pre_budget%igrid_sim%u
        vbase => this%pre_budget%igrid_sim%v
        wcbase => this%pre_budget%igrid_sim%wC
        
        ! Term 1: delta u_j' d_j(delta p')
        ! Term 2: base  u_j' d_j(delta p')        
        ! px, py, pz signs are reversed
        cbuffyC1 = this%pre_budget%px - this%prim_budget%px
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, bf)
        this%budget_3(:,:,:,1)=this%budget_3(:,:,:,1)+ bf * du
        this%budget_3(:,:,:,2)=this%budget_3(:,:,:,2)+ bf * ubase

        cbuffyC1 = this%pre_budget%py - this%prim_budget%py
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, bf)
        this%budget_3(:,:,:,1)=this%budget_3(:,:,:,1)+ bf * dv
        this%budget_3(:,:,:,2)=this%budget_3(:,:,:,2)+ bf * vbase

        cbuffyE1 = this%pre_budget%pz - this%prim_budget%pz
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, bf)
        this%budget_3(:,:,:,1)=this%budget_3(:,:,:,1)+ bf * dw        
        this%budget_3(:,:,:,2)=this%budget_3(:,:,:,2)+ bf * wcbase

        ! Term 3: delta u_j' d_j(base p')
        ! px, py, pz signs are reversed
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%px, bf)
        this%budget_3(:,:,:,3)=this%budget_3(:,:,:,3)- bf * du

        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%py, bf)
        this%budget_3(:,:,:,3)=this%budget_3(:,:,:,3)- bf * dv

        call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%pz, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, bf)
        this%budget_3(:,:,:,3)=this%budget_3(:,:,:,3)- bf * dw

        ! Term 4: d_j(base  u_i' * delta tau_ij') [SGS transport] 
        ! Term 6: d_j(delta u_i' * delta tau_ij')  [SGS transport]
        ! sign of usgs, vsgs, and wsgs are reversed.
        cbuffyC1 = this%pre_budget%usgs - this%prim_budget%usgs
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, bf)
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + bf * ubase 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + bf * du

        cbuffyC1 = this%pre_budget%vsgs - this%prim_budget%vsgs
        call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, bf) 
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + bf * vbase  
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + bf * dv

        cbuffyE1 = this%pre_budget%wsgs - this%prim_budget%wsgs
        call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, bf)
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + bf * wcbase
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + bf * dw
        
        ! The remaining of B3(4) is exactly B3(7). Calculation is done once
        ! Term 7: delta tau_ij' d_j(base u_i')     [SGS dissipation]  
        ! Term 13: d_j(delta u_j' base u_i' base u_i')/2  [Turbulent transport of TKE]
        ! Term 14: d_j(base  u_j' base u_i' delta u_i')   [Turbulent transport of TKE] 
        ! Term 15: d_j(delta u_j' base u_i' delta u_i')   [Turbulent transport of TKE]
        call this%ddx_R2R(ubase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,1) ! i=1, j=1
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        buffer =   bf * du * ubase
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + buffer
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + buffer
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * du * du
        
        call this%ddy_R2R(ubase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,2) ! i=1, j=2
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + bf * dv * ubase
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + bf * vbase * du
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * dv * du
        
        call this%ddz_R2R(ubase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,3) ! i=1, j=3
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + bf * dw * ubase
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + bf * wcbase * du
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * dw * du

        call this%ddx_R2R(vbase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,2) ! i=2, j=1
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + bf * du * vbase
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + bf * ubase * dv
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * du * dv
        
        call this%ddy_R2R(vbase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,4) ! i=2, j=2
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        buffer =   bf * dv * vbase
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + buffer
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + buffer
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * dv * dv
        
        call this%ddz_R2R(vbase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,5) ! i=2, j=3
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + bf * dw * vbase
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + bf * wcbase * dv
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * dw * dv

        call this%ddx_R2R(wcbase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,3) ! i=3, j=1
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + bf * du * wcbase
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + bf * ubase * dw
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * du * dw
        
        call this%ddy_R2R(wcbase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,5) ! i=3, j=2
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + bf * dv * wcbase
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + bf * vbase * dw
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * dv * dw
        
        call this%ddz_R2R(wcbase,bf)
        buffer =   bf * this%delta_tauij(:,:,:,6) ! i=3, j=3
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + buffer
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + buffer
        buffer =   bf * dw * wcbase
        this%budget_3(:,:,:,13) = this%budget_3(:,:,:,13) + buffer
        this%budget_3(:,:,:,14) = this%budget_3(:,:,:,14) + buffer
        this%budget_3(:,:,:,15) = this%budget_3(:,:,:,15) + bf * dw * dw

        ! Term 5: d_j(delta  u_i' base tau_ij') [SGS transport]         
        ! sign of usgs, vsgs, and wsgs are reversed. 
        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%usgs, bf)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - bf * du

        call this%pre_budget%igrid_sim%spectC%ifft(this%pre_budget%vsgs, bf)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - bf * dv 

        call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%wsgs, rbuffxE1)
        call this%interp_Edge2Cell(rbuffxE1, bf)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - bf * dw

        ! The remaining of B3(5) is the exactly as B3(8)
        ! Term 8: base  tau_ij' * d_j(delta u_i')     [SGS dissipation]
        ! Do the rest of B3(6): d_j(delta u_i' * delta tau_ij')  [SGS transport]
        ! Term 9: delta tau_ij' * d_j(delta u_i')     [SGS dissipation]
        ! Term 14: d_j(base  u_j' base u_i' delta u_i')  [Turbulent transport of TKE] 
        ! Term 15: d_j(delta u_j' base u_i' delta u_i')  [Turbulent transport of TKE]
        ! Term 16: d_j(base  u_j' delta u_i' delta u_i')/2  [Turbulent transport of TKE]
        ! Term 17: d_j(delta u_j' delta u_i' delta u_i')/2 [Turbulent transport of TKE] 
        
        call this%ddx_R2R(du, bf)! i=1, j=1
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,1) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer  
        buffer =   bf * this%delta_tauij(:,:,:,1)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer  
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * ubase * ubase
        buffer =   bf * du * ubase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ buffer
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ buffer
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * du * du         
        
        call this%ddy_R2R(du, bf)! i=1, j=2
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,2) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer  
        buffer =   bf * this%delta_tauij(:,:,:,2)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer  
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer  
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * ubase * vbase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ bf * dv * ubase
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ bf * vbase * du
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * dv * du  
        
        call this%ddz_R2R(du, bf)! i=1, j=3
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,3) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer
        buffer =   bf * this%delta_tauij(:,:,:,3)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * ubase * wcbase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ bf * dw * ubase
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ bf * wcbase * du
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * dw * du
        
        call this%ddx_R2R(dv, bf)! i=2, j=1
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,2) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer
        buffer =   bf * this%delta_tauij(:,:,:,2)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * vbase * ubase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ bf * du * vbase
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ bf * ubase * dv
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * du * dv
        
        call this%ddy_R2R(dv, bf)! i=2, j=2
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,4) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer
        buffer =   bf * this%delta_tauij(:,:,:,4)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * vbase * vbase
        buffer =   bf * dv * vbase 
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ buffer
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ buffer
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * dv * dv
        
        call this%ddz_R2R(dv, bf)! i=2, j=3
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,5) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer
        buffer =   bf * this%delta_tauij(:,:,:,5)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * vbase * wcbase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ bf * dw * vbase
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ bf * wcbase * dv
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * dw * dv
        
        call this%ddx_R2R(dw, bf)! i=3, j=1
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,3) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer
        buffer =   bf * this%delta_tauij(:,:,:,3)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * wcbase * ubase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ bf * du * wcbase
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ bf * ubase * dw
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * du * dw
        
        call this%ddy_R2R(dw, bf)! i=3, j=2
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,5) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer
        buffer =   bf * this%delta_tauij(:,:,:,5)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * wcbase * vbase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ bf * dv * wcbase
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ bf * vbase * dw
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * dv * dw
        
        call this%ddz_R2R(dw, bf)! i=3, j=3
        buffer =   bf * this%pre_budget%igrid_sim%tauSGS_ij(:,:,:,6) 
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + buffer
        this%budget_3(:,:,:,8) = this%budget_3(:,:,:,8) + buffer    
        buffer =   bf * this%delta_tauij(:,:,:,6)   
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + buffer 
        this%budget_3(:,:,:,9) = this%budget_3(:,:,:,9) + buffer 
        this%budget_3(:,:,:,14)= this%budget_3(:,:,:,14)+ bf * wcbase * wcbase
        buffer =   bf * dw * wcbase
        this%budget_3(:,:,:,15)= this%budget_3(:,:,:,15)+ buffer
        this%budget_3(:,:,:,16)= this%budget_3(:,:,:,16)+ buffer
        this%budget_3(:,:,:,17)= this%budget_3(:,:,:,17)+ bf * dw * dw 

        ! Term 10: delta u_3' delta wb'
        ! Term 11: delta u_3' base wb'
        ! Term 12: base u_3' delta wb'
        if(this%isStratified)then
            cbuffyE1 = this%prim_budget%wb - this%pre_budget%wb 
            call this%prim_budget%igrid_sim%spectE%ifft(cbuffyE1, rbuffxE1)
            call this%interp_Edge2Cell(rbuffxE1, buffer)        
            this%budget_3(:,:,:,10) = this%budget_3(:,:,:,10) + dw * buffer
            this%budget_3(:,:,:,12) = this%budget_3(:,:,:,12) + wcbase * buffer
            
            call this%pre_budget%igrid_sim%spectE%ifft(this%pre_budget%wb, rbuffxE1)
            call this%interp_Edge2Cell(rbuffxE1, buffer)
            this%budget_3(:,:,:,11) = this%budget_3(:,:,:,11) + dw * buffer    
        end if  

        if (this%useWindTurbines)then
            cbuffyC1 = this%prim_budget%uturb - this%pre_budget%uturb
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer)
            this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + du * buffer 
            this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + ubase * buffer

            cbuffyC1 = this%prim_budget%vturb - this%pre_budget%vturb
            call this%prim_budget%igrid_sim%spectC%ifft(cbuffyC1, buffer)
            this%budget_3(:,:,:,18) = this%budget_3(:,:,:,18) + dv * buffer 
            this%budget_3(:,:,:,19) = this%budget_3(:,:,:,19) + vbase * buffer
        end if 

        nullify(du, dv, dw, rbuffxE1, buffer, bf, cbuffyE1, cbuffyC1, ubase, vbase, wcbase)        
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
                call this%ddx_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_0(:,:,:,1)
                call this%ddy_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_0(:,:,:,2)
                call this%ddz_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_0(:,:,:,3)
            case(2)
                call this%ddx_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_0(:,:,:,1)
                call this%ddy_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_0(:,:,:,2)
                call this%ddz_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_0(:,:,:,3)
            case(3)
                call this%ddx_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_0(:,:,:,1)
                call this%ddy_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_0(:,:,:,2)
                call this%ddz_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_0(:,:,:,3)
            case(4)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_0(:,:,:,1)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_0(:,:,:,2)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_0(:,:,:,3)
            case(5)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_0(:,:,:,1)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_0(:,:,:,2)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_0(:,:,:,3)
            case(6)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_0(:,:,:,1)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_0(:,:,:,2)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_0(:,:,:,3)
            case(7)
                call this%ddx_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,1)
                call this%ddy_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,2)
                call this%ddz_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,3)
            case(8)
                call this%ddx_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,1)
                call this%ddy_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,2)
                call this%ddz_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,3)
            case(9)
                call this%ddx_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,1)
                call this%ddy_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,2)
                call this%ddz_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,3)
            case(10)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,1)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,2)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,3)
            case(11)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,1)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,2)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,3)
            case(12)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,1)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,2)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,3)
            end select
            call this%dealias(buffer)

        else if(budgetid.eq.3)then
            select case(idx)
            case(1) ! d_j(delta u_j' delta p')
                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,18)
                buffer = buffer + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,19)
                buffer = buffer + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,20)
            case(2) ! d_j(base  u_j' delta p')
                buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,18)
                buffer = buffer + this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,19)
                buffer = buffer + this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,20)

            case(3) ! d_j(delta u_j' base p')
                ! px, py, pz signs are reversed in base-flow budget
                buffer = buffer - this%budget_0(:,:,:,1)*this%pre_budget%budget_1(:,:,:,2)
                buffer = buffer - this%budget_0(:,:,:,2)*this%pre_budget%budget_1(:,:,:,6)
                buffer = buffer - this%budget_0(:,:,:,3)*this%pre_budget%budget_1(:,:,:,9)

            case(4) ! d_j(base u_i' delta tau_ij') [SGS transport]
                buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,12)
                buffer = buffer + this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,13)
                buffer = buffer + this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,14)

                ! The rest of the term is the same as that of B3(7)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,6) ! i=1, j=1
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,7) ! i=1, j=2
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,8) ! i=1, j=3

                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,7) ! i=2, j=1
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,9) ! i=2, j=2
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,10) ! i=2, j=3

                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,8) ! i=3, j=1
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,10) ! i=3, j=2
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,11) ! i=3, j=3

            case(5) ! d_j(delta u_i' base tau_ij') [SGS transport]                
                ! The sign of ui_sgs in this%pre_budget%budget_1 is reversed
                buffer = buffer - this%budget_0(:,:,:,1)*this%pre_budget%budget_1(:,:,:,3)
                buffer = buffer - this%budget_0(:,:,:,2)*this%pre_budget%budget_1(:,:,:,7)
                buffer = buffer - this%budget_0(:,:,:,3)*this%pre_budget%budget_1(:,:,:,10)

                ! The rest of this term is the same as B3(8)
                call this%ddx_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,11) ! i=1, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,12) ! i=1, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,13) ! i=1, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,12) ! i=2, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,14) ! i=2, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,15) ! i=2, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,13) ! i=3, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,15) ! i=3, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,16) ! i=3, j=3

            case(6) ! d_j(delta u_i' * delta tau_ij')  [SGS transport]
                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,12)
                buffer = buffer + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,13)
                buffer = buffer + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,14)

                ! The rest of this term is the same as B3(9)
                call this%ddx_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,6) ! i=1, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,7) ! i=1, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,8) ! i=1, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,7) ! i=2, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,9) ! i=2, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,10) ! i=2, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,8) ! i=3, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,10) ! i=3, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,11) ! i=3, j=3

            case(7) ! delta tau_ij' * d_j(base u_i')     [SGS dissipation]
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,6) ! i=1, j=1
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,7) ! i=1, j=2
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,8) ! i=1, j=3

                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,7) ! i=2, j=1
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,9) ! i=2, j=2
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,10) ! i=2, j=3

                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,8) ! i=3, j=1
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,10) ! i=3, j=2
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3),bf)
                buffer=buffer + bf * this%budget_0(:,:,:,11) ! i=3, j=3

            case(8) ! base  tau_ij' * d_j(delta u_i')     [SGS dissipation]                
                call this%ddx_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,11) ! i=1, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,12) ! i=1, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,13) ! i=1, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,12) ! i=2, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,14) ! i=2, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,15) ! i=2, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,13) ! i=3, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,15) ! i=3, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%pre_budget%budget_0(:,:,:,16) ! i=3, j=3

            case(9) ! delta tau_ij' * d_j(delta u_i')     [SGS dissipation]
                call this%ddx_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,6) ! i=1, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,7) ! i=1, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,1),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,8) ! i=1, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,7) ! i=2, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,9) ! i=2, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,2),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,10) ! i=2, j=3

                call this%ddx_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,8) ! i=3, j=1
                call this%ddy_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,10) ! i=3, j=2
                call this%ddz_R2R(this%budget_0(:,:,:,3),bf)
                buffer=buffer+bf*this%budget_0(:,:,:,11) ! i=3, j=3

            case(10) ! delta u_3' delta wb'
                buffer = this%budget_0(:,:,:,3)*this%budget_0(:,:,:,17)

            case(11) ! delta u_3' base wb'
                buffer = this%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,31)

            case(12) ! base u_3' delta wb'
                buffer = this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,17)

            case(13) ! d_j(delta u_j' base u_i' base u_i')/2  [Turbulent transport of TKE]
                bf = 0.5d0*(this%pre_budget%budget_0(:,:,:,4) + this%pre_budget%budget_0(:,:,:,7) + this%pre_budget%budget_0(:,:,:,9)) &
                     - (this%pre_budget%budget_0(:,:,:,1)*this%pre_budget%budget_0(:,:,:,1) + &
                        this%pre_budget%budget_0(:,:,:,2)*this%pre_budget%budget_0(:,:,:,2) + &
                        this%pre_budget%budget_0(:,:,:,3)*this%pre_budget%budget_0(:,:,:,3))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,2)*bf2
                call this%ddz_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*this%budget_2(:,:,:,4) + &
                                  this%pre_budget%budget_0(:,:,:,2)*this%budget_2(:,:,:,5) + &
                                  this%pre_budget%budget_0(:,:,:,3)*this%budget_2(:,:,:,6)
                
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,7)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,9)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,11)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,8)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,12)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,14)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,10)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,13)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,15)
            
            case(14) ! d_j(base  u_j' base u_i' delta u_i')   [Turbulent transport of TKE]
                bf = this%budget_1(:,:,:,7) + this%budget_1(:,:,:,12) + this%budget_1(:,:,:,15) - &
                     2.d0 * (this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) + &
                            this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) + &
                            this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,3))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,2)*bf2
                call this%ddz_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*this%budget_2(:,:,:,7) + &
                                  this%pre_budget%budget_0(:,:,:,2)*this%budget_2(:,:,:,8) + &
                                  this%pre_budget%budget_0(:,:,:,3)*this%budget_2(:,:,:,9)

                call this%ddx_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,4)
                call this%ddy_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,5)
                call this%ddz_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,6)
                call this%ddx_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,5)
                call this%ddy_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,7)
                call this%ddz_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,8)
                call this%ddx_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,6)
                call this%ddy_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,8)
                call this%ddz_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%pre_budget%budget_0(:,:,:,9)

                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_2(:,:,:,10) + &
                                  this%budget_0(:,:,:,2)*this%budget_2(:,:,:,11) + &
                                  this%budget_0(:,:,:,3)*this%budget_2(:,:,:,12)

                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,7)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,8)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,10)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,9)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,12)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,13)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,11)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,14)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,15)
            
            case(15) ! d_j(delta u_j' base u_i' delta u_i')  [Turbulent transport of TKE]
                bf = this%budget_1(:,:,:,7) + this%budget_1(:,:,:,12) + this%budget_1(:,:,:,15) - &
                     2.d0 * (this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) + &
                            this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) + &
                            this%pre_budget%budget_0(:,:,:,3)*this%budget_0(:,:,:,3))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,2)*bf2
                call this%ddz_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*this%budget_2(:,:,:,1) + &
                                  this%pre_budget%budget_0(:,:,:,2)*this%budget_2(:,:,:,2) + &
                                  this%pre_budget%budget_0(:,:,:,3)*this%budget_2(:,:,:,3)
                call this%ddx_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,7)
                call this%ddy_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,9)
                call this%ddz_R2R(this%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,11)
                call this%ddx_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,8)
                call this%ddy_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,12)
                call this%ddz_R2R(this%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,14)
                call this%ddx_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,10)
                call this%ddy_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,13)
                call this%ddz_R2R(this%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,15)

                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_2(:,:,:,4) + &
                                  this%budget_0(:,:,:,2)*this%budget_2(:,:,:,5) + &
                                  this%budget_0(:,:,:,3)*this%budget_2(:,:,:,6)

                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,1)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,2)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,1), bf); buffer = buffer + bf*this%budget_1(:,:,:,3)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,2)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,4)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,2), bf); buffer = buffer + bf*this%budget_1(:,:,:,5)
                call this%ddx_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,3)
                call this%ddy_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,5)
                call this%ddz_R2R(this%pre_budget%budget_0(:,:,:,3), bf); buffer = buffer + bf*this%budget_1(:,:,:,6)
            
            case(16) ! d_j(base  u_j' delta u_i' delta u_i')/2  [Turbulent transport of TKE]
                bf = 0.5d0*(this%budget_1(:,:,:,1) + this%budget_1(:,:,:,4) + this%budget_1(:,:,:,6)) &
                     - (this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) + &
                        this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) + &
                        this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,2)*bf2
                call this%ddz_R2R(bf, bf2); buffer = buffer + this%pre_budget%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_2(:,:,:,7) + &
                                  this%budget_0(:,:,:,2)*this%budget_2(:,:,:,8) + &
                                  this%budget_0(:,:,:,3)*this%budget_2(:,:,:,9)
                
                call this%ddx_R2R(this%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,7)
                call this%ddy_R2R(this%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,8)
                call this%ddz_R2R(this%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,10)
                call this%ddx_R2R(this%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,9)
                call this%ddy_R2R(this%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,12)
                call this%ddz_R2R(this%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,13)
                call this%ddx_R2R(this%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,11)
                call this%ddy_R2R(this%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,14)
                call this%ddz_R2R(this%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,15)
            
            case(17) ! d_j(delta u_j' delta u_i' delta u_i')/2 [Turbulent transport of TKE]
                bf = 0.5d0*(this%budget_1(:,:,:,1) + this%budget_1(:,:,:,4) + this%budget_1(:,:,:,6)) &
                     - (this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) + &
                        this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) + &
                        this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3))
                call this%ddx_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,1)*bf2
                call this%ddy_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,2)*bf2
                call this%ddz_R2R(bf, bf2); buffer = buffer + this%budget_0(:,:,:,3)*bf2

                buffer = buffer + this%budget_0(:,:,:,1)*this%budget_2(:,:,:,1) + &
                                  this%budget_0(:,:,:,2)*this%budget_2(:,:,:,2) + &
                                  this%budget_0(:,:,:,3)*this%budget_2(:,:,:,3)
                
                call this%ddx_R2R(this%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,1)
                call this%ddy_R2R(this%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,2)
                call this%ddz_R2R(this%budget_0(:,:,:,1),bf); buffer=buffer+bf*this%budget_1(:,:,:,3)
                call this%ddx_R2R(this%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,2)
                call this%ddy_R2R(this%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,4)
                call this%ddz_R2R(this%budget_0(:,:,:,2),bf); buffer=buffer+bf*this%budget_1(:,:,:,5)
                call this%ddx_R2R(this%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,3)
                call this%ddy_R2R(this%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,5)
                call this%ddz_R2R(this%budget_0(:,:,:,3),bf); buffer=buffer+bf*this%budget_1(:,:,:,6)

            case(18)
                buffer = this%budget_0(:,:,:,1)*this%budget_0(:,:,:,21) + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,22)
            case(19)
                buffer = this%pre_budget%budget_0(:,:,:,1)*this%budget_0(:,:,:,21) + this%pre_budget%budget_0(:,:,:,2)*this%budget_0(:,:,:,22)
            end select
        end if
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

        nullify(buffer)         
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
        end if
    end subroutine 

    ! ----------------------private derivative operators ------------------------
    subroutine dealias(this, f)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(inout) :: f
        
        call this%prim_budget%igrid_sim%spectC%fft(f,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), f)
    end subroutine

    subroutine ddx_R2R(this, f, dfdx)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdx
        
        call this%prim_budget%igrid_sim%spectC%fft(f,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%mtimes_ik1_ip(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdx)
    end subroutine 

    subroutine ddy_R2R(this, f, dfdy)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdy
        
        call this%prim_budget%igrid_sim%spectC%fft(f,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%mtimes_ik2_ip(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdy)
    end subroutine 
     
    subroutine ddz_R2R(this, f, dfdz)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdz

        call this%prim_budget%igrid_sim%spectC%fft(f,this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%ddz_C2R(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdz)
    end subroutine
     
    subroutine ddz_C2R(this, fhat, dfdz)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        complex(rkind), dimension(this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(1),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(2),this%prim_budget%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: dfdz
        
        call transpose_y_to_z(fhat,this%prim_budget%igrid_sim%cbuffzC(:,:,:,1),this%prim_budget%igrid_sim%sp_gpC)
        call this%prim_budget%igrid_sim%Pade6opZ%ddz_C2C(this%prim_budget%igrid_sim%cbuffzC(:,:,:,1),this%prim_budget%igrid_sim%cbuffzC(:,:,:,2),0,0)
        call transpose_z_to_y(this%prim_budget%igrid_sim%cbuffzC(:,:,:,2),this%prim_budget%igrid_sim%cbuffyC(:,:,:,1),this%prim_budget%igrid_sim%sp_gpC)
        call this%prim_budget%igrid_sim%spectC%dealias(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1))
        call this%prim_budget%igrid_sim%spectC%ifft(this%prim_budget%igrid_sim%cbuffyC(:,:,:,1), dfdz)
    end subroutine 
 
    subroutine interp_Edge2Cell(this, fE, fC)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(in) :: fE
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: fC

        call transpose_x_to_y(fE,this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),this%prim_budget%igrid_sim%gpE)
        call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%gpE)
        call this%prim_budget%igrid_sim%Pade6opZ%interpz_E2C(this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,2),0,0)
        call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzC(:,:,:,2),this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
        call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),fC,this%prim_budget%igrid_sim%gpC)
    end subroutine 
 
    subroutine interp_Cell2Edge(this, fC, fE)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: fC
        real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(out) :: fE

        call transpose_x_to_y(fC,this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
        call transpose_y_to_z(this%prim_budget%igrid_sim%rbuffyC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%gpC)
        call this%prim_budget%igrid_sim%Pade6opZ%interpz_C2E(this%prim_budget%igrid_sim%rbuffzC(:,:,:,1),this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),0,0)
        call transpose_z_to_y(this%prim_budget%igrid_sim%rbuffzE(:,:,:,1),this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),this%prim_budget%igrid_sim%gpE)
        call transpose_y_to_x(this%prim_budget%igrid_sim%rbuffyE(:,:,:,1),fE,this%prim_budget%igrid_sim%gpE)
    end subroutine 
         
    subroutine multiply_CellFieldsOnEdges(this, f1C, f2C, fmultC)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(in) :: f1C,f2C
        real(rkind), dimension(this%nx,this%ny,this%nz), intent(out) :: fmultC

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

    ! multiply on edge cells and interpolate to cell centers to reduce aliasing issues
    function multiply_Edges_interp_cell(this, f1E, f2E) result(fmultC)
        class(budgets_time_avg_deficit_compact), intent(inout) :: this
        real(rkind), dimension(this%prim_budget%igrid_sim%gpE%xsz(1),this%prim_budget%igrid_sim%gpE%xsz(2),this%prim_budget%igrid_sim%gpE%xsz(3)), intent(in) :: f1E,f2E
        real(rkind), dimension(this%prim_budget%igrid_sim%gpC%xsz(1),this%prim_budget%igrid_sim%gpC%xsz(2),this%prim_budget%igrid_sim%gpC%xsz(3)) :: fmultC

        call this%interp_Edge2Cell(f1E * f2E, fmultC)
    end function
end module