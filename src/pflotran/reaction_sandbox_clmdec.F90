module Reaction_Sandbox_CLMDec_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

! ------------------------------------------------------------------------------
! Description
! to be used to implement CLM-CN, and CLM-Microbe decomposition reactions
! extended from reaction_sandbox_clm_cn  
! 1) pools can be either immobile or aqueous species (e.g., DOM, acetate-, )
! 2) separate N into NH3 (or NH4+) and NO3-; Must have NH3 or NH4+, NO3- is used 
!    if it is specified in the input file
! 3) include flexibilities to have multiple downstream pools, and variable 
!    respiration fraction as in CLM-Microbe; 
! 4) add residual concentrations for upstream pools, NH3, and NO3- to 
!    keep reactant concentrations above 0 (used if > 0); 
! 5) add shut off down regulation for NH3 and NO3- (used when the first > 0)
! 6) include NH3 oxidation in decomposition using Parton et al. 2001 (used when 
!    N2O(aq) is specified in the input file) 
! 7) add optional immobile species to track respiration, N mineralization, and 
!    immobilization
! Author: Guoping Tang
! Date:   07/08/14 
! -----------------------------------------------------------------------------

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  PetscInt, parameter :: LITTER_DECOMP_CLMCN = 1 
  PetscInt, parameter :: LITTER_DECOMP_CLMMICROBE = 2 

                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter :: CN_ratio_mass_to_mol = 1.16616d0

  ! Sinsabaugh et al. 2013 Ecology Letters, 16, 930-939
  PetscReal, parameter :: CN_ratio_microbe = 9.32928d0   ! 8.0d0 
  PetscReal, parameter :: CUE_max = 0.6d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clmdec_type

    PetscInt  :: temperature_response_function
    PetscReal :: Q10
    PetscInt  :: moisture_response_function

    PetscInt  :: litter_decomp_type          ! CLM-CN or CLM-Microbe

    PetscReal :: half_saturation_nh3
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh3_no3
    PetscReal :: n2o_frac_mineralization     ! fraction of n2o from net N mineralization

    PetscReal :: residual_conc_cpool
    PetscReal :: residual_conc_nh3
    PetscReal :: residual_conc_no3

    PetscReal :: downreg_no3_0               ! shut off
    PetscReal :: downreg_no3_1               ! start to decrease from 1
    PetscReal :: downreg_nh3_0               ! shut off
    PetscReal :: downreg_nh3_1               ! start to decrease from 1

    PetscReal :: net_n_min_rate_smooth_0     ! start from 0
    PetscReal :: net_n_min_rate_smooth_1     ! rise to 1

    PetscReal :: nc_bacteria
    PetscReal :: nc_fungi
    PetscReal :: fraction_bacteria

    PetscInt :: npool                        ! litter or variable CN ration pools
    PetscReal, pointer :: pool_nc_ratio(:)   ! NC ratio in mole  npool

    PetscInt :: nrxn
    PetscReal, pointer :: rate_constant(:)           ! nrxn

    PetscBool, pointer :: is_litter_decomp(:)        ! nrxn
    PetscInt,  pointer :: upstream_c_id(:)           ! nrxn
    PetscInt,  pointer :: upstream_n_id(:)           ! nrxn
    PetscReal, pointer :: upstream_nc(:)             ! nrxn
    PetscBool, pointer :: upstream_is_aqueous(:)     ! nrxn

    PetscInt,  pointer :: n_downstream_pools(:)      ! maximum # of downstream pools
    PetscInt,  pointer :: downstream_id(:,:)         ! nrxn by maximum # of downstream pools
    PetscBool, pointer :: downstream_is_aqueous(:,:) ! nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_stoich(:,:)     ! nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_nc(:,:)         ! nrxn by maximum # of downstream pools
    PetscReal, pointer :: mineral_c_stoich(:)        ! nrxn
    PetscReal, pointer :: mineral_n_stoich(:)        ! nrxn

    PetscInt :: species_id_co2
    PetscInt :: species_id_nh3
    PetscInt :: species_id_no3
    PetscInt :: species_id_n2o
    PetscInt :: species_id_dom
    PetscInt :: species_id_bacteria
    PetscInt :: species_id_fungi

    PetscInt :: species_id_hrimm
    PetscInt :: species_id_nmin
    PetscInt :: species_id_nimm
    PetscInt :: species_id_ngasmin
    PetscInt :: species_id_proton

    type(pool_type), pointer :: pools
    type(clmdec_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => CLMDec_Read
    procedure, public :: Setup => CLMDec_Setup
    procedure, public :: Evaluate => CLMDec_React
    procedure, public :: Destroy => CLMDec_Destroy
  end type reaction_sandbox_clmdec_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: stoich
    PetscReal :: nc_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: clmdec_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    type(pool_type), pointer :: downstream_pools
    PetscReal :: rate_constant
    type(clmdec_reaction_type), pointer :: next
  end type clmdec_reaction_type
  
  public :: CLMDec_Create

contains

! **************************************************************************** !

function CLMDec_Create()
  ! Allocates CLMDec reaction sandbox object.

  implicit none
  
  type(reaction_sandbox_clmdec_type), pointer :: CLMDec_Create
  
  allocate(CLMDec_Create)

  CLMDec_Create%Q10 = 1.5d0
  CLMDec_Create%litter_decomp_type=LITTER_DECOMP_CLMCN
  CLMDec_Create%half_saturation_nh3 = 1.0d-15
  CLMDec_Create%half_saturation_no3 = 1.0d-15
  CLMDec_Create%inhibition_nh3_no3 = 1.0d-15
  CLMDec_Create%n2o_frac_mineralization = 0.02d0  ! Parton et al. 2001

  CLMDec_Create%residual_conc_cpool = 1.0d-20
  CLMDec_Create%residual_conc_nh3 = 1.0d-20
  CLMDec_Create%residual_conc_no3 = 1.0d-20

  CLMDec_Create%downreg_no3_0 = -1.0d-9 
  CLMDec_Create%downreg_no3_1 = 1.0d-7
  CLMDec_Create%downreg_nh3_0 = -1.0d-9 
  CLMDec_Create%downreg_nh3_1 = 1.0d-7

  CLMDec_Create%net_n_min_rate_smooth_0 = 0.0d0 
  CLMDec_Create%net_n_min_rate_smooth_1 = 1.0d-20

  CLMDec_Create%nc_bacteria = 0.17150d0

  ! CN_ratio_fungi = 17.4924d0     !15.0d0 ! or 10.0
  CLMDec_Create%nc_fungi = 0.05717d0

  CLMDec_Create%fraction_bacteria = 0.340927d0

  CLMDec_Create%npool = 0

  nullify(CLMDec_Create%pool_nc_ratio)

  CLMDec_Create%nrxn = 0
  nullify(CLMDec_Create%rate_constant)
  nullify(CLMDec_Create%is_litter_decomp)
  nullify(CLMDec_Create%upstream_c_id)
  nullify(CLMDec_Create%upstream_n_id)
  nullify(CLMDec_Create%upstream_nc)
  nullify(CLMDec_Create%upstream_is_aqueous)
  
  nullify(CLMDec_Create%n_downstream_pools)
  nullify(CLMDec_Create%downstream_id)
  nullify(CLMDec_Create%downstream_is_aqueous)
  nullify(CLMDec_Create%downstream_stoich)
  nullify(CLMDec_Create%mineral_c_stoich)
  nullify(CLMDec_Create%mineral_n_stoich)

  CLMDec_Create%species_id_co2 = 0
  CLMDec_Create%species_id_nh3 = 0
  CLMDec_Create%species_id_no3 = 0
  CLMDec_Create%species_id_n2o = 0
  CLMDec_Create%species_id_dom = 0
  CLMDec_Create%species_id_proton = 0
  CLMDec_Create%species_id_bacteria = 0
  CLMDec_Create%species_id_fungi = 0
  CLMDec_Create%species_id_hrimm = 0
  CLMDec_Create%species_id_nmin = 0
  CLMDec_Create%species_id_nimm = 0
  CLMDec_Create%species_id_ngasmin = 0

  nullify(CLMDec_Create%next)
  nullify(CLMDec_Create%pools)
  nullify(CLMDec_Create%reactions)

end function CLMDec_Create

! **************************************************************************** !

subroutine CLMDec_Read(this,input,option)
  ! 
  ! Reads input deck for reaction sandbox parameters
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal
 
  implicit none
  
  class(reaction_sandbox_clmdec_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(pool_type), pointer :: new_pool_rxn, prev_pool_rxn
  type(clmdec_reaction_type), pointer :: new_reaction, prev_reaction
  
  PetscReal :: rate_constant, turnover_time
  PetscReal :: temp_real
  
  nullify(new_pool)
  nullify(prev_pool)

  nullify(new_pool_rxn)
  nullify(prev_pool_rxn)

  nullify(new_reaction)
  nullify(prev_reaction)
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
      'CHEMISTRY,REACTION_SANDBOX,CLMDec')
    call StringToUpper(word)   
    select case(trim(word))

      case('CLM-MICROBE-LITTER-DECOMPOSITION')
        this%litter_decomp_type = LITTER_DECOMP_CLMMICROBE    

      case('RESIDUAL_CONC_CPOOL')
        call InputReadDouble(input,option,this%residual_conc_cpool)
        call InputErrorMsg(input,option,'residual_con_cpool', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

      case('RESIDUAL_CONC_NH3')
        call InputReadDouble(input,option,this%residual_conc_nh3)
        call InputErrorMsg(input,option,'residual_con_nh3', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

      case('RESIDUAL_CONC_NO3')
        call InputReadDouble(input,option,this%residual_conc_no3)
        call InputErrorMsg(input,option,'residual_con_no3', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

      case('HALF_SATURATION_NH3')
        call InputReadDouble(input,option,this%half_saturation_nh3)
        call InputErrorMsg(input,option,'NH3 half saturation', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

      case('HALF_SATURATION_NO3')
        call InputReadDouble(input,option,this%half_saturation_no3)
        call InputErrorMsg(input,option,'NO3 half saturation', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

      case('DOWNREGULATE_NH3')
        call InputReadDouble(input,option,this%downreg_nh3_0)
        call InputErrorMsg(input,option,'downreg_nh3_0', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
        call InputReadDouble(input,option,this%downreg_nh3_1)
        call InputErrorMsg(input,option,'downreg_nh3_1', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
        if (this%downreg_nh3_0 > this%downreg_nh3_1) then
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec,' // &
            'NH4+ down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('DOWNREGULATE_NO3')
        call InputReadDouble(input,option,this%downreg_no3_0)
        call InputErrorMsg(input,option,'downreg_no3_0', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
        call InputReadDouble(input,option,this%downreg_no3_1)
        call InputErrorMsg(input,option,'downreg_no3_1', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
        if (this%downreg_no3_0 > this%downreg_no3_1) then
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec,' // &
            'NO3- down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif

      case('SMOOTH_NET_N_MINERALIZATION')
        call InputReadDouble(input,option,this%net_n_min_rate_smooth_0)
        call InputErrorMsg(input,option,'net_n_min_rate_smooth_0', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
        call InputReadDouble(input,option,this%net_n_min_rate_smooth_1)
        call InputErrorMsg(input,option,'net_n_min_rate_smooth_1', &
          'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
        if (this%net_n_min_rate_smooth_0 > this%net_n_min_rate_smooth_1) then
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec,' // &
            'Net N mineralization smooth 0 concentration > 1 concentration.'
          call printErrMsg(option)
        endif

     case('NH3_INHIBITION_NO3')
       call InputReadDouble(input,option,this%inhibition_nh3_no3)
       call InputErrorMsg(input,option,'NH3 inhibition coefficient', &
         'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

     case('N2O_FRAC_MINERALIZATION')
       call InputReadDouble(input,option,this%n2o_frac_mineralization)
       call InputErrorMsg(input,option,'n2o fraction from mineralization', &
         'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

     case('POOLS')
       do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit   

         allocate(new_pool)
         new_pool%name = ''
         new_pool%nc_ratio = UNINITIALIZED_DOUBLE
         nullify(new_pool%next)

         call InputReadWord(input,option,new_pool%name,PETSC_TRUE)
         call InputErrorMsg(input,option,'pool name', &
           'CHEMISTRY,REACTION_SANDBOX,CLMDec,POOLS')
         call InputReadDouble(input,option,temp_real)
         if (InputError(input)) then
           new_pool%nc_ratio = UNINITIALIZED_DOUBLE
         else
           ! convert CN ratio from mass C/mass N to mol C/mol N
           if (temp_real > 0.0d0 ) then
             new_pool%nc_ratio = 1.0d0/temp_real/CN_ratio_mass_to_mol
           endif
         endif

         if (associated(this%pools)) then
           prev_pool%next => new_pool
         else
           this%pools => new_pool
         endif
         prev_pool => new_pool
         nullify(new_pool)
       enddo

      case('REACTION')
      
        allocate(new_reaction)
        new_reaction%upstream_pool_name = ''
        new_reaction%rate_constant = UNINITIALIZED_DOUBLE
        nullify(new_reaction%downstream_pools)
        nullify(new_reaction%next)
        
        ! need to set these temporarily in order to check that they
        ! are not both set.
        turnover_time = 0.d0
        rate_constant = 0.d0
        
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
            case('DOWNSTREAM_POOL')
              allocate(new_pool_rxn)
              new_pool_rxn%name = ''
              new_pool_rxn%stoich = 0.d0
              nullify(new_pool_rxn%next)

              call InputReadWord(input,option, &
                 new_pool_rxn%name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
              call InputReadDouble(input,option,new_pool_rxn%stoich)
              call InputErrorMsg(input,option,'Downstream pool stoich', &
                'CHEMISTRY,REACTION_SANDBOX_CLMDec ' // &
                'TEMPERATURE RESPONSE FUNCTION')

              if (associated(new_reaction%downstream_pools)) then
                prev_pool_rxn%next => new_pool_rxn
              else
                new_reaction%downstream_pools => new_pool_rxn
              endif
              prev_pool_rxn => new_pool_rxn
              nullify(new_pool_rxn)

            case('RATE_CONSTANT')
              call InputReadDouble(input,option,rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMDec RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,option)
              endif
            case('TURNOVER_TIME')
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMDec TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,option)
              endif
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
              call printErrMsg(option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
            'be included in a CLMDec reaction definition, but not both. ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        else if (turnover_time > 0.d0) then
          new_reaction%rate_constant = 1.d0 / turnover_time
        else
          new_reaction%rate_constant = rate_constant
        endif
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine CLMDec_Read

! **************************************************************************** !

subroutine CLMDec_Setup(this,reaction,option)
  ! 
  ! Sets up CLMDec reaction after it has been read from input
  ! 

  use Reaction_Aux_module
  use Option_module
  use String_module
  use Immobile_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(reaction_sandbox_clmdec_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH), allocatable :: pool_names(:)
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt, pointer :: species_id_pool_c(:)
  PetscInt, pointer :: species_id_pool_n(:)
  PetscBool, pointer :: pool_is_aqueous(:)

  PetscInt :: icount, jcount, max_downstream_pools, ipool
  PetscReal :: stoich_c, stoich_n

  type(pool_type), pointer :: cur_pool
  type(clmdec_reaction_type), pointer :: cur_rxn
  
  ! count # pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    cur_pool => cur_pool%next
  enddo
  this%npool = icount
  
  ! count # reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    cur_rxn => cur_rxn%next
  enddo
  this%nrxn = icount
 
  allocate(this%n_downstream_pools(this%nrxn))
 
  ! count # downstream pools in each reaction
  max_downstream_pools = -1
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1

    jcount = 0
    cur_pool => cur_rxn%downstream_pools

    do
      if (.not.associated(cur_pool)) exit
      jcount = jcount + 1
      cur_pool => cur_pool%next
    enddo

    this%n_downstream_pools(icount) = jcount

    if (max_downstream_pools < jcount) then
      max_downstream_pools = jcount
    endif 

    cur_rxn => cur_rxn%next
  enddo

  ! allocate and initialize arrays
  allocate(this%pool_nc_ratio(this%npool))

  allocate(this%rate_constant(this%nrxn))
  allocate(this%is_litter_decomp(this%nrxn))
  allocate(this%upstream_c_id(this%nrxn))
  allocate(this%upstream_n_id(this%nrxn))
  allocate(this%upstream_nc(this%nrxn))
  allocate(this%upstream_is_aqueous(this%nrxn))
  
  allocate(this%downstream_id(this%nrxn,max_downstream_pools))
  allocate(this%downstream_stoich(this%nrxn,max_downstream_pools))
  allocate(this%downstream_nc(this%nrxn,max_downstream_pools))
  allocate(this%downstream_is_aqueous(this%nrxn,max_downstream_pools))
  allocate(this%mineral_c_stoich(this%nrxn))
  allocate(this%mineral_n_stoich(this%nrxn))

  this%pool_nc_ratio = 0.d0
  this%rate_constant = 0.d0
  this%is_litter_decomp = PETSC_FALSE
  this%upstream_c_id = 0
  this%upstream_n_id = 0
  this%upstream_nc = UNINITIALIZED_DOUBLE
  this%upstream_is_aqueous = PETSC_FALSE

  this%downstream_id = 0
  this%downstream_is_aqueous = PETSC_FALSE
  this%downstream_stoich = 0.d0
  this%mineral_c_stoich = 0.d0
  this%mineral_n_stoich = 0.d0
  
  ! temporary array for mapping pools in reactions
  allocate(pool_names(this%npool))
  allocate(pool_is_aqueous(this%npool))
  allocate(species_id_pool_c(this%npool))
  allocate(species_id_pool_n(this%npool))

  pool_names = ''
  pool_is_aqueous = PETSC_FALSE
  species_id_pool_c = UNINITIALIZED_INTEGER 
  species_id_pool_n = UNINITIALIZED_INTEGER 

  ! pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    this%pool_nc_ratio(icount) = cur_pool%nc_ratio
    pool_names(icount) = cur_pool%name

    if (cur_pool%nc_ratio < 0.d0) then
      ! Since no CN ratio provided, must provide two species with the
      ! same name as the pool with C or N appended.
      word = trim(cur_pool%name) // 'C'
      species_id_pool_c(icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      word = trim(cur_pool%name) // 'N'
      species_id_pool_n(icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      if (species_id_pool_c(icount)<=0 .or. species_id_pool_n(icount)<=0) then
        option%io_buffer = 'For CLMDec pools with no CN ratio defined, ' // &
          'the user must define two immobile species with the same root ' // &
          'name as the pool with "C" or "N" appended, respectively.'
        call printErrMsg(option)
      endif
    else ! only one species (e.g. SOMX)
      species_id_pool_c(icount) = &
        GetImmobileSpeciesIDFromName(cur_pool%name,reaction%immobile, &
                                     PETSC_FALSE,option)
      if (species_id_pool_c(icount) <= 0) then
        species_id_pool_c(icount) = GetPrimarySpeciesIDFromName( &
          cur_pool%name, reaction, PETSC_FALSE,option)
        if (species_id_pool_c(icount) <= 0) then
          option%io_buffer = 'CLMDec pool: ' // cur_pool%name // 'is not ' // &
            'specified either in the IMMOBILE_SPECIES or PRIMARY_SPECIES!'
          call printErrMsg(option)
        else
          pool_is_aqueous(icount) = PETSC_TRUE
        endif
      endif
      
      if (StringCompare(cur_pool%name, 'Bacteria')) then
        this%nc_bacteria = cur_pool%nc_ratio
      endif 

      if (StringCompare(cur_pool%name, 'Fungi')) then
        this%nc_fungi = cur_pool%nc_ratio
      endif 
    endif
    cur_pool => cur_pool%next
  enddo
 
  ! reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    ! upstream pools
    icount = icount + 1
    ipool = StringFindEntryInList(cur_rxn%upstream_pool_name,pool_names)
    if (ipool == 0) then
      option%io_buffer = 'Upstream pool ' // &
        trim(cur_rxn%upstream_pool_name) // &
        'in reaction not found in list of pools.'
      call printErrMsg(option)
    else
      this%upstream_c_id(icount) = species_id_pool_c(ipool)
      this%upstream_n_id(icount) = species_id_pool_n(ipool)
      this%upstream_nc(icount) = this%pool_nc_ratio(ipool) 
      this%upstream_is_aqueous(icount) = pool_is_aqueous(ipool) 
      if (this%upstream_n_id(icount) > 0) then
        this%is_litter_decomp(icount) = PETSC_TRUE
      else
        if (this%upstream_nc(icount) < 0.0d0) then
          option%io_buffer = 'SOM decomp. reaction with upstream pool ' // &
            trim(cur_rxn%upstream_pool_name) // &
            'has negative C:N ratio in upstream pool.'
          call printErrMsg(option)
        endif
      endif
    endif

    ! downstream pools
    jcount = 0
    cur_pool => cur_rxn%downstream_pools

    do
      if (.not.associated(cur_pool)) exit
      jcount = jcount + 1

      if (len_trim(cur_pool%name) > 0) then
        ipool = StringFindEntryInList(cur_pool%name,pool_names)
        if (ipool == 0) then
          option%io_buffer = 'Downstream pool "' // trim(cur_pool%name) // &
            '" in reaction with upstream pool "' // &
            trim(cur_rxn%upstream_pool_name) // '" not found in list of pools.'
          call printErrMsg(option)
        else
          this%downstream_id(icount, jcount) = species_id_pool_c(ipool)
          this%downstream_stoich(icount, jcount) = cur_pool%stoich 
          this%downstream_nc(icount, jcount) = this%pool_nc_ratio(ipool) 
          this%downstream_is_aqueous(icount, jcount) = pool_is_aqueous(ipool) 

          if (this%downstream_nc(icount,jcount) < 0.d0) then
            option%io_buffer = 'For CLMDec reactions, downstream pools ' // &
              'must have a constant C:N ratio (i.e. C and N are not ' // &
              'tracked individually).  Therefore, pool "' // &
              trim(cur_pool%name) // &
             '" may not be used as a downstream pool.'
            call printErrMsg(option)
          endif
        endif
      endif

      cur_pool => cur_pool%next

    enddo

    this%rate_constant(icount) = cur_rxn%rate_constant
    cur_rxn => cur_rxn%next
  enddo 
  
  deallocate(pool_names)
  call DeallocateArray(pool_is_aqueous)
  call DeallocateArray(species_id_pool_c)
  call DeallocateArray(species_id_pool_n)

  ! set stoichiometric coefficients for som decomposition reactions  
  ! as they are constant due to fixed CN ratio
  do icount = 1, this%nrxn
    if (this%is_litter_decomp(icount)) then
      cycle
    else
      ! calculate respiration factor
      stoich_c = 1.0d0
      stoich_n = this%upstream_nc(icount)

      do jcount = 1, this%n_downstream_pools(icount)
        stoich_c = stoich_c - this%downstream_stoich(icount, jcount)
        stoich_n = stoich_n - this%downstream_stoich(icount, jcount) * &
                              this%downstream_nc(icount, jcount)
      enddo

      if (stoich_c < -1.0d-10) then
        option%io_buffer = 'CLMDec SOM decomposition reaction has negative' // &
          ' respiration fraction!'
        call printErrMsg(option)
      endif

      this%mineral_c_stoich(icount) = stoich_c
      this%mineral_n_stoich(icount) = stoich_n

     endif
  enddo

  word = 'HCO3-'
  this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%species_id_co2 < 0) then
     word = 'CO2(aq)'
     this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  if (this%species_id_co2 <= 0) then
    option%io_buffer = 'Neither HCO3- nor CO2(aq) is specified in the ' // &
      'input file for CLMDec!'
    call printErrMsg(option)
  endif

  word = 'NH4+'
  this%species_id_nh3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%species_id_nh3 < 0) then
    word = 'NH3(aq)'
    this%species_id_nh3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  if (this%species_id_nh3 <= 0) then
    option%io_buffer = 'Neither NH4+ nor NH3(aq) is specified in the input' // &
      'file for CLMDec!'
    call printErrMsg(option)
  endif

  word = 'NO3-'
  this%species_id_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  
  word = 'N2O(aq)'
  this%species_id_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'H+'
  this%species_id_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'Bacteria'
  this%species_id_bacteria = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  word = 'Fungi'
  this%species_id_fungi = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'HRimm'
  this%species_id_hrimm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  word = 'Nmin'
  this%species_id_nmin = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'Nimm'
  this%species_id_nimm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'NGASmin'
  this%species_id_ngasmin = GetImmobileSpeciesIDFromName( &
    word,reaction%immobile,PETSC_FALSE,option)

  if (this%species_id_bacteria > 0 .and. this%species_id_fungi > 0 .and. & 
    this%nc_bacteria > 0.0d0 .and. this%nc_fungi > 0.0d0 ) then
    this%fraction_bacteria = (1.0d0/this%nc_bacteria) ** 0.6d0 / & 
      ((1.0d0/this%nc_bacteria) ** 0.6d0 + (1.0d0/this%nc_fungi) ** 0.6d0) 
  endif 

end subroutine CLMDec_Setup

! ************************************************************************** !
subroutine CLMDec_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                       global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type

  implicit none

  class(reaction_sandbox_clmdec_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscReal :: saturation
  PetscErrorCode :: ierr

  PetscInt :: local_id
  PetscReal :: theta
  PetscReal :: psi

  PetscReal :: c_nh3              ! concentration (mole/L)
  PetscReal :: ac_nh3             ! activity coefficient
  PetscReal :: f_nh3              ! nh3 / (half_saturation + nh3)
  PetscReal :: d_nh3              ! half_saturation / (half_saturation + nh3)^2
  PetscReal :: f_nh3_inhibit      ! inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: d_nh3_inhibit_dnh3 ! -inhibition_coef/(inhibition_coef + nh3)^2

  PetscReal :: c_no3              ! concentration (mole/L)
  PetscReal :: ac_no3             ! activity coefficient 
  PetscReal :: f_no3              ! no3 / (half_saturation + no3)
  PetscReal :: d_no3              ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: temp_real

  PetscInt :: irxn
  PetscInt :: ipool_up, ipool_down
  PetscReal :: CN_ratio_up, CN_ratio_down
  PetscBool :: constant_CN_ratio_up
  PetscReal :: resp_frac

  PetscReal :: c_uc, c_un         ! upstream c pool, n pool concentration
  PetscReal :: ac_uc, ac_un       ! activity coefficient, if aqueous

  PetscInt :: ispec_uc, ispec_un, ispec_d   ! species id for upstream C, N, and downstream
  PetscInt :: ires_uc, ires_un, ires_d      ! id used for residual and Jacobian
  PetscInt :: ires_co2, ires_nh3, ires_n2o, ires_no3
  PetscInt :: ires_hrimm, ires_nmin, ires_nimm, ires_ngasmin
  PetscReal :: stoich_c, stoich_n

  PetscReal :: scaled_rate_const

  PetscReal :: rate_nh3        ! mole/s 
  PetscReal :: drate_nh3_duc   ! d Rate / d upstream c
  PetscReal :: drate_nh3_dnh3  ! d Rate / d nh3 ammonia limitation

  PetscReal :: Rdu_duc, Rdn_duc, Rdc_duc, Rdb_duc, Rdf_duc  ! u = Lit1N/Lit1C, c, b, f for CLM-Microbe
  PetscReal :: Rdu_dun, Rdn_dun, Rdc_dun, Rdb_dun, Rdf_dun
  PetscReal :: Rno3du_duc, Rno3dn_duc, Rno3dc_duc, Rno3db_duc, Rno3df_duc
  PetscReal :: Rno3du_dun, Rno3dn_dun, Rno3dc_dun, Rno3db_dun, Rno3df_dun

  ! for N immobilization reactions with NO3 as N source 
  PetscReal :: rate_no3       ! mole/s
  PetscReal :: drate_no3_dno3 ! d Rate_no3 / d no3 
  PetscReal :: drate_no3_duc  ! d Rate_no3 / d uc 
  PetscReal :: drate_no3_dnh3 ! d Rate_no3 / d nh3 

  PetscInt :: i, j
  PetscReal :: tc     ! temperature in C
  PetscReal :: f_t    ! temperature response function
  PetscReal :: f_w    ! moisture response function

  ! save mineral N fraction and decomposition rate for net N mineralization and N2O calculation 
  PetscReal :: net_n_mineralization_rate
  PetscReal :: dnet_n_mineralization_rate_dnh3
  PetscReal :: dnet_n_mineralization_rate_dno3
  PetscReal :: dnet_n_mineralization_rate_duc(this%nrxn)
  PetscReal :: ph, f_ph
  PetscReal :: rate_n2o, drate_n2o_dnh3, drate_n2o_dno3, drate_n2o_duc
  PetscReal :: f_rate_n2o, df_rate_n2o

  PetscInt :: ires_b, ires_f
  PetscReal :: xxx, delta, regulator, dregulator

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume
    
  c_nh3     = rt_auxvar%pri_molal(this%species_id_nh3)
  ac_nh3    = rt_auxvar%pri_act_coef(this%species_id_nh3)
  temp_real = (c_nh3 - this%residual_conc_nh3) * ac_nh3 &
            + this%half_saturation_nh3
  f_nh3     = (c_nh3 - this%residual_conc_nh3) * ac_nh3 / temp_real 
  if (compute_derivative) then
    d_nh3     = ac_nh3 * this%half_saturation_nh3 / temp_real / temp_real
  endif

  if (this%downreg_nh3_0 > 0.0d0) then
    ! additional down regulation for NH4+ immobilization 
    if (c_nh3 <= this%downreg_nh3_0) then
      regulator = 0.0d0
      dregulator = 0.0d0
    elseif (c_nh3 >= this%downreg_nh3_1 .or. &
            this%downreg_nh3_1 - this%downreg_nh3_0 <= 1.0d-20) then
      regulator = 1.0d0
      dregulator = 0.0d0
    else
      xxx = c_nh3 - this%downreg_nh3_0
      delta = this%downreg_nh3_1 - this%downreg_nh3_0
      regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
      dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                 / delta / delta
    endif
    
    ! rate = rate_orginal * regulator
    ! drate = drate_original * regulator + rate_orginal * dregulator
    if (compute_derivative) then
      d_nh3 = d_nh3 * regulator + f_nh3 * dregulator
    endif

    f_nh3 = f_nh3 * regulator

  endif

  if (this%species_id_no3 > 0) then
    c_no3     = rt_auxvar%pri_molal(this%species_id_no3)
    ac_no3    = rt_auxvar%pri_act_coef(this%species_id_no3)
    temp_real = (c_no3 - this%residual_conc_no3) * ac_no3 &
              + this%half_saturation_no3
    f_no3 = (c_no3 - this%residual_conc_no3) * ac_no3 / temp_real

    if (compute_derivative) then
      d_no3 = ac_no3 * this%half_saturation_no3 / temp_real / temp_real 
    endif

    if (this%downreg_no3_0 > 0.0d0) then
      ! additional down regulation for NO3- immobilization
      if (c_no3 <= this%downreg_no3_0) then
        regulator = 0.0d0
        dregulator = 0.0d0
      elseif (c_no3 >= this%downreg_no3_1 .or. &
              this%downreg_no3_1 - this%downreg_no3_0 <= 1.0d-20) then
        regulator = 1.0d0
        dregulator = 0.0d0
      else
        xxx = c_no3 - this%downreg_no3_0
        delta = this%downreg_no3_1 - this%downreg_no3_0
        regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                   / delta / delta
      endif

      ! rate = rate_orginal * regulator
      ! drate = drate_original * regulator + rate_orginal * dregulator
      if (compute_derivative) then
        d_no3 = d_no3 * regulator + f_no3 * dregulator
      endif
      f_no3 = f_no3 * regulator

    endif

    temp_real = this%inhibition_nh3_no3 + c_nh3 * ac_nh3
    f_nh3_inhibit = this%inhibition_nh3_no3/temp_real
    if (compute_derivative) then
      d_nh3_inhibit_dnh3 = -1.0d0 * this%inhibition_nh3_no3 * ac_nh3 &
                         / temp_real / temp_real
    endif
  endif 

  ires_co2 = this%species_id_co2
  ires_nh3 = this%species_id_nh3
  ires_no3 = this%species_id_no3
  ires_n2o = this%species_id_n2o

  if (this%species_id_hrimm > 0) then
    ires_hrimm = this%species_id_hrimm + reaction%offset_immobile
  endif

  if (this%species_id_nmin > 0) then
    ires_nmin = this%species_id_nmin + reaction%offset_immobile
  endif

  if (this%species_id_nimm > 0) then
    ires_nimm = this%species_id_nimm + reaction%offset_immobile
  endif

  if (this%species_id_ngasmin > 0) then
    ires_ngasmin = this%species_id_ngasmin + reaction%offset_immobile
  endif

  ! temperature response function
  tc = global_auxvar%temp

  f_t = 1.0d0

  saturation = global_auxvar%sat(1)
  theta = saturation * porosity 
  ! if positive, saturated soil's psi is nearly zero
  psi = min(global_auxvar%pres(1) - option%reference_pressure, -1.d-20)   

  ! moisture response function 
  f_w = 1.0d0

  if(f_t < 1.0d-20 .or. f_w < 1.0d-20) then
     return
  endif

  if (this%species_id_n2o > 0) then
    net_n_mineralization_rate = 0.0d0
    dnet_n_mineralization_rate_dnh3 = 0.0d0
    dnet_n_mineralization_rate_dno3 = 0.0d0
    do irxn = 1, this%nrxn
      dnet_n_mineralization_rate_duc(irxn) = 0.0d0
    enddo
  endif

  do irxn = 1, this%nrxn
  
    ! upstream pool
    ispec_uc = this%upstream_c_id(irxn)

    if (this%upstream_is_aqueous(irxn)) then
      c_uc    = rt_auxvar%pri_molal(ispec_uc)
      ac_uc   = rt_auxvar%pri_act_coef(ispec_uc)
      ires_uc = ispec_uc
    else
      c_uc    = rt_auxvar%immobile(ispec_uc)
      ac_uc   = 1.0d0
      ires_uc = reaction%offset_immobile + ispec_uc
    endif

    ! for litter decomposition reactions, stoich needs to be calculated on the fly
    if (this%is_litter_decomp(irxn)) then

      ispec_un = this%upstream_n_id(irxn)
      if (this%upstream_is_aqueous(irxn)) then
        c_un    = rt_auxvar%pri_molal(ispec_un)
        ac_un   = rt_auxvar%pri_act_coef(ispec_un)
        ires_un = ispec_un
      else
        c_un    = rt_auxvar%immobile(ispec_un)
        ac_un   = 1.0d0
        ires_un = ispec_un + reaction%offset_immobile
      endif
      this%upstream_nc(irxn) = c_un / c_uc

      if (this%litter_decomp_type == LITTER_DECOMP_CLMCN) then

        ! calculate respiration factor (CO2 stoichiometry)
        stoich_c = 1.0d0

        do j = 1, this%n_downstream_pools(irxn)
          stoich_c = stoich_c - this%downstream_stoich(irxn, j)
        enddo

        if (stoich_c < 0.0d0) then
          option%io_buffer = 'CLMDec litter decomposition reaction has' // &
                             'negative respiration fraction!'
          call printErrMsg(option)
        endif

        this%mineral_c_stoich(irxn) = stoich_c

        ! calculate N stoichiometry
        stoich_n = this%upstream_nc(irxn)

        do j = 1, this%n_downstream_pools(irxn)
          stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                              this%downstream_nc(irxn, j)
        enddo

        this%mineral_n_stoich(irxn) = stoich_n

      elseif(this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then

        ! Sinsabaugh et al. 2013 Ecology Letters, 16, 930-939
        resp_frac = CN_ratio_microbe * this%upstream_nc(irxn) !c_un/c_uc

        if(resp_frac > CUE_max) then
          resp_frac = CUE_max
        endif

        ! c pools
        this%mineral_c_stoich(irxn) = resp_frac
       
        if (this%n_downstream_pools(irxn) .ne. 2) then
          option%io_buffer = 'CLM_Microbe litter decomposition reaction ' // &
                              'more than 2 (bacteria and fungi pools)!'
          call printErrMsg(option)
        endif

        do i = 1, this%n_downstream_pools(irxn)
          if (this%downstream_id(irxn, i) == this%species_id_bacteria) then
            this%downstream_stoich(irxn, i) = this%fraction_bacteria * &
              (1.0d0  - resp_frac)
          else
            this%downstream_stoich(irxn, i) = (1.0d0 - this%fraction_bacteria) &
                                            * (1.0d0  - resp_frac)
          endif
        enddo

        stoich_n = this%upstream_nc(irxn)

        do j = 1, this%n_downstream_pools(irxn)
          stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                                this%downstream_nc(irxn, j)
        enddo

        this%mineral_n_stoich(irxn) = stoich_n

      endif

    endif

    if (this%upstream_is_aqueous(irxn)) then
      ! scaled_rate_const units: (kg water / s) = (1/s) * (kg water)
      scaled_rate_const = this%rate_constant(irxn) * volume * f_t * f_w &
                        * 1000.0d0 * theta
      ! will need to replace 1000.0d0 with water density
    else
      ! scaled_rate_const units: (m^3 bulk / s) = (1/s) * (m^3 bulk)
      scaled_rate_const = this%rate_constant(irxn)*volume*f_t*f_w
    endif

    ! residual units: (mol/sec) = (kg water/s) * (mol/kg water) or
    ! residual units: (mol/sec) = (m^3 bulk/s) * (mol/m^3 bulk)
    rate_nh3 = scaled_rate_const * (c_uc - this%residual_conc_cpool) * ac_uc

    if (compute_derivative) then
      drate_nh3_duc = scaled_rate_const * ac_uc
    endif

    ! NH3 limiting
    if (this%mineral_n_stoich(irxn) < 0.0d0) then
      if (compute_derivative) then
        drate_nh3_dnh3 = rate_nh3 * d_nh3 
      endif
      rate_nh3       = rate_nh3 * f_nh3
      if (compute_derivative) then
        drate_nh3_duc  = drate_nh3_duc * f_nh3
      endif
    else
      drate_nh3_dnh3 = 0.d0
    endif 

    ! CO2
    Residual(ires_co2) = Residual(ires_co2) - &
                         this%mineral_c_stoich(irxn) * rate_nh3
    if (this%species_id_hrimm > 0) then
      Residual(ires_hrimm) = Residual(ires_hrimm) - &
                             this%mineral_c_stoich(irxn) * rate_nh3
    endif
    
    ! NH3
    Residual(ires_nh3) = Residual(ires_nh3) - &
                         this%mineral_n_stoich(irxn) * rate_nh3

    if (this%species_id_nimm > 0 .and. this%mineral_n_stoich(irxn) < 0.0d0) then
      Residual(ires_nimm) = Residual(ires_nimm) + &
                            this%mineral_n_stoich(irxn) * rate_nh3
    endif

    if(this%species_id_nmin > 0 .and. this%mineral_n_stoich(irxn) > 0.0d0) then
       Residual(ires_nmin) = Residual(ires_nmin) - &
                             this%mineral_n_stoich(irxn) * rate_nh3
    endif

    ! upstream c
    Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * rate_nh3

    ! upstream n
    if (this%is_litter_decomp(irxn)) then
      Residual(ires_un) = Residual(ires_un) - &
                          (-1.d0) * this%upstream_nc(irxn) * rate_nh3
    endif
    
    ! downstream pools
    do j = 1, this%n_downstream_pools(irxn)
      ispec_d = this%downstream_id(irxn, j)
      if (this%downstream_is_aqueous(irxn, j)) then
        ires_d = ispec_d
      else
        ires_d = reaction%offset_immobile + ispec_d
      endif
      if (ispec_d > 0) then
        Residual(ires_d) = Residual(ires_d) - &
                           this%downstream_stoich(irxn, j) * rate_nh3
      endif
    enddo

    if (this%species_id_n2o > 0) then
      net_n_mineralization_rate = net_n_mineralization_rate + &
        this%mineral_n_stoich(irxn) * rate_nh3

      if (compute_derivative) then
        dnet_n_mineralization_rate_dnh3 = dnet_n_mineralization_rate_dnh3 + &
          this%mineral_n_stoich(irxn) * drate_nh3_dnh3
        dnet_n_mineralization_rate_duc(irxn) = &
          dnet_n_mineralization_rate_duc(irxn) + &
          this%mineral_n_stoich(irxn) * drate_nh3_duc
      endif
    endif

    ! start residual calculation for N immobilization reaction with NO3 uptake
    ! if nitrate is available, N immobilization decomposition reactions occurs
    ! with rate depending on NH3, with reduced rate if NH3 is abundent  
    if (this%species_id_no3 > 0 .and. this%mineral_n_stoich(irxn) < 0.d0) then

      rate_no3       = scaled_rate_const * (c_uc - this%residual_conc_cpool) &
                     * ac_uc * f_no3 * f_nh3_inhibit

      if (compute_derivative) then
        drate_no3_duc  = scaled_rate_const * ac_uc * f_no3 * f_nh3_inhibit
        drate_no3_dno3 = scaled_rate_const * (c_uc - this%residual_conc_cpool) &
                       * ac_uc * d_no3 * f_nh3_inhibit
        drate_no3_dnh3 = scaled_rate_const * (c_uc - this%residual_conc_cpool) &
                       * ac_uc * f_no3 * d_nh3_inhibit_dnh3
      endif

      ! carbon
      Residual(ires_co2) = Residual(ires_co2) - &
        this%mineral_c_stoich(irxn) * rate_no3
      if (this%species_id_hrimm > 0) then
        Residual(ires_hrimm) = Residual(ires_hrimm) - &
        this%mineral_c_stoich(irxn) * rate_no3
      endif
    
      ! NO3
      Residual(ires_no3) = Residual(ires_no3) - &
        this%mineral_n_stoich(irxn) * rate_no3

      if (this%species_id_nimm > 0) then
        Residual(ires_nimm) = Residual(ires_nimm) + &
          this%mineral_n_stoich(irxn) * rate_no3
      endif

      ! upstream c
      Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * rate_no3
    
      ! upstream n
      if (this%is_litter_decomp(irxn)) then
        Residual(ires_un) = Residual(ires_un) + this%upstream_nc(irxn) *rate_no3
      endif
    
      ! downstream pools
      do j = 1, this%n_downstream_pools(irxn)
        ispec_d = this%downstream_id(irxn, j)
        if (this%downstream_is_aqueous(irxn, j)) then
          ires_d = ispec_d
        else
          ires_d = reaction%offset_immobile + ispec_d
        endif
        if (ispec_d > 0) then
          Residual(ires_d) = Residual(ires_d) - &
            this%downstream_stoich(irxn, j) * rate_no3
         endif
      enddo

      if (this%species_id_n2o > 0) then
        net_n_mineralization_rate = net_n_mineralization_rate + &
          this%mineral_n_stoich(irxn)*rate_no3
        if (compute_derivative) then
          dnet_n_mineralization_rate_dnh3 = dnet_n_mineralization_rate_dnh3 + &
            this%mineral_n_stoich(irxn) * drate_no3_dnh3
          dnet_n_mineralization_rate_dno3 = dnet_n_mineralization_rate_dno3 + &
            this%mineral_n_stoich(irxn) * drate_no3_dno3
          dnet_n_mineralization_rate_duc(irxn) = &
            dnet_n_mineralization_rate_duc(irxn) + &
            this%mineral_n_stoich(irxn) * drate_no3_duc
        endif
      endif 
    endif

    if (compute_derivative) then

      if (this%is_litter_decomp(irxn)) then
        if (this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
          ! LitC + u LitN -> di SOMi + (1 - di) CO2 + n N
          ! Rdu/duc = R (-1) LitN/LitC^2 = - u R / LitC 
          Rdu_duc = -1.0d0 * this%upstream_nc(irxn) * rate_nh3 / c_uc
         
          ! n = u - (1 - di) ni
          ! dn/dLitC = du/dLitC
          Rdn_duc = Rdu_duc

          ! Rdu/dun = R /LitC 
          Rdu_dun = rate_nh3 / c_uc 

          ! Rdn/dun = Rdu/dLitN = Rdu/dun
          Rdn_dun = Rdu_dun

        elseif (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! Lit1C + u Lit1N -> b Bacteria + f Fungi + c CO2 + n N
          ! c = min(CUEmax, Lit1N/Lit1C*CN_ratio_microbe)
          ! g = CNbacteria^0.6/(CNbacterial^0.6 + CNfungi^0.6)
          ! n = u - b nb - f nf
          ! b = g (1 - c)
          ! f = (1 - g) (1 - c)

          Rdu_duc = -1.0d0 * this%upstream_nc(irxn) * rate_nh3 / c_uc
   
          if (resp_frac < CUE_max) then
            ! Rdc/dLit1C = -RLit1N/Lit1C^2*CN_ratio_microbe           
            Rdc_duc = -1.0d0 * this%upstream_nc(irxn) * rate_nh3 / c_uc  &
                    * CN_ratio_microbe
          else
            Rdc_duc = 0.0d0 
          endif

          ! Rdb/dLitC = -g Rdc/dLitC
          Rdb_duc = -1.0d0 * this%fraction_bacteria * Rdc_duc
  
          ! Rdf/dLitC = -(1 - g) Rdc/dLitC
          Rdf_duc = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rdc_duc

          ! Rdn/dLitC = Rdu/dLitC - nb Rdb/dLitC - nf Rdf/dLitC
          Rdn_duc = Rdu_duc - this%nc_bacteria * Rdb_duc &
                  - this%nc_fungi * Rdf_duc

          ! Rdu/dun = R/LitC = dR/duc
          Rdu_dun = rate_nh3 / c_uc 

          if (resp_frac < CUE_max) then
            ! Rdc/dLitN = R/LitC*CN_ratio_microbe
            Rdc_dun = rate_nh3 / c_uc * CN_ratio_microbe
          else
            Rdc_dun = 0.0d0 
          endif

          ! Rdb/dLitN = -g Rdc/dLitN
          Rdb_dun = -1.0d0 * this%fraction_bacteria * Rdc_dun

          ! Rdf/dLitN = -(1 - g) Rdc/dLitN
          Rdf_dun = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rdc_dun

          ! Rdn/dLitN = Rdu/dLitN - nb Rdb/dLitN - nf Rdf/dLitN
          Rdn_dun = Rdu_dun - this%nc_bacteria * Rdb_dun &
                  - this%nc_fungi * Rdf_dun

          ires_b = reaction%offset_immobile + this%species_id_bacteria 
          ires_f = reaction%offset_immobile + this%species_id_fungi 
            
        endif

      endif   

      ! with respect to upstream C
      ! CO2
      Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
        this%mineral_c_stoich(irxn) * drate_nh3_duc

      if (this%is_litter_decomp(irxn) .and. &
        this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
        ! dRco2/dLitC = dcR/dLitC = cdR/dLitC + R dc/dLitC
        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - Rdc_duc
      endif

      if (this%species_id_hrimm > 0) then
        Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_nh3_duc

        if (this%is_litter_decomp(irxn) .and. &
           this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
           Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - Rdc_duc
        endif
      endif

      ! N
      ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
      ! first, n dR/dC_u
      Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
        this%mineral_n_stoich(irxn) * drate_nh3_duc

      if (this%species_id_nimm > 0 .and. &
          this%mineral_n_stoich(irxn) < 0.0d0) then
        Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + &
          this%mineral_n_stoich(irxn) * drate_nh3_duc
      endif

      if (this%species_id_nmin > 0 .and. &
          this%mineral_n_stoich(irxn) > 0.0d0) then
        Jacobian(ires_nmin,ires_uc) = Jacobian(ires_nmin,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_nh3_duc
      endif

      if (this%is_litter_decomp(irxn)) then
        ! litter pool is immobile
        ! second, Rdn/dC_u
        Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - Rdn_duc

        if (this%species_id_nimm > 0 .and. &
            this%mineral_n_stoich(irxn) < 0.0d0) then
          Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + Rdn_duc
        endif

        if (this%species_id_nmin > 0 .and. &
            this%mineral_n_stoich(irxn) > 0.0d0) then
          Jacobian(ires_nmin,ires_uc) = Jacobian(ires_nmin,ires_uc) - Rdn_duc
        endif

      endif

      ! upstream C pool
      Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) + drate_nh3_duc

      ! upstream N pool
      if (this%is_litter_decomp(irxn)) then
        ! litter pools are immobile
        ! R_Nu = Nu/Cu * R_Cu
        ! dR_Nu/dCu = Nu/Cu dR_Cu/dCu + R du/dCu
        ! = 0 only when residual_conc_cpool = 0
        Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) + &
          this%upstream_nc(irxn) * drate_nh3_duc + Rdu_duc
        
      endif

      ! downstream pools
      do j = 1, this%n_downstream_pools(irxn)
        ispec_d = this%downstream_id(irxn, j)
        if (ispec_d < 0) then
          option%io_buffer = 'Downstream pool species not specified!'
          call printErrMsg(option)
        endif

        if (this%downstream_is_aqueous(irxn, j)) then
          ires_d = ispec_d
        else
          ires_d = reaction%offset_immobile + ispec_d
        endif
         
        Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
          this%downstream_stoich(irxn, j) * drate_nh3_duc

        ! additional term if downstream stoich is variable
        if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          if (ispec_d == this%species_id_bacteria) then
            ! dRbacteria/dLit1C = b dR/dLit1C + R db/dLit1C
            Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rdb_duc
          elseif(ispec_d == this%species_id_fungi) then
            Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rdf_duc
          else
            option%io_buffer = 'Downstream pool for CLM-Microbe should be' // &
                               'either bacteria or fungi!'
            call printErrMsg(option)
          endif
        endif
      enddo

      ! with respect to upstream n (due to variable CN ratio)
      if (this%is_litter_decomp(irxn)) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
        Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) + Rdu_dun

        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu
        Jacobian(ires_nh3,ires_un) = Jacobian(ires_nh3,ires_un) - Rdn_dun

        if (this%species_id_nimm > 0 .and. &
            this%mineral_n_stoich(irxn) < 0.0d0) then
          Jacobian(ires_nimm,ires_un) = Jacobian(ires_nimm,ires_un) + Rdn_dun
        endif

        if (this%species_id_nmin > 0 .and. &
            this%mineral_n_stoich(irxn) > 0.0d0) then
          Jacobian(ires_nmin,ires_un) = Jacobian(ires_nmin,ires_un) - Rdn_dun
        endif

        if (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! CO2 dR_co2/dNu = d cR/dNu = Rdc/Nu
          Jacobian(ires_co2,ires_un) = Jacobian(ires_co2,ires_un) - Rdc_dun
       
          if (this%species_id_hrimm > 0) then
            Jacobian(ires_hrimm,ires_un) = Jacobian(ires_hrimm,ires_un) -Rdc_dun
          endif

          ! bacteria dRbacteria/dNu = Rdb/dNu
          Jacobian(ires_b,ires_un) = Jacobian(ires_b,ires_un) - Rdb_dun

          ! fungi dRfungi/dNu = Rdf/dNu
          Jacobian(ires_f,ires_un) = Jacobian(ires_f,ires_un) - Rdf_dun
        endif 
      endif

      ! with respect to nh3
      if (this%mineral_n_stoich(irxn) < 0.0d0) then
        ! CO2
        Jacobian(ires_co2,ires_nh3) = Jacobian(ires_co2,ires_nh3) - &
          this%mineral_c_stoich(irxn) * drate_nh3_dnh3

        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_nh3) = Jacobian(ires_hrimm,ires_nh3) - &
            this%mineral_c_stoich(irxn) * drate_nh3_dnh3
        endif

        ! N
        Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) - &
          this%mineral_n_stoich(irxn) * drate_nh3_dnh3
  
        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_nh3) = Jacobian(ires_nimm,ires_nh3) + &
            this%mineral_n_stoich(irxn) * drate_nh3_dnh3
        endif

        ! upstream C
        Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
          (-1.d0) * drate_nh3_dnh3
 
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_nh3_dnh3 
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif
          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
            this%downstream_stoich(irxn, j) * drate_nh3_dnh3
        enddo

      endif

      if (this%species_id_no3 > 0 .and. this%mineral_n_stoich(irxn) < 0.d0) then

        if (this%is_litter_decomp(irxn)) then
          if (this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
            ! Lit1C + u Lit1N -> di SOMi + (1 - di) CO2 + n N
            ! Rdu/duc = R (-1) Lit1N/Lit1C^2 
            Rno3du_duc = -1.0d0 * this%upstream_nc(irxn) * rate_no3 / c_uc

            ! n = u - (1 - di) ni
            ! dn/dLit1C = du/dLit1C
            Rno3dn_duc = Rno3du_duc

            ! Rdu/dun = R /Lit1C 
            Rno3du_dun = rate_no3 / c_uc

            ! Rdn/dun = du/dLit1C = dR/duc 
            Rno3dn_dun = Rno3du_dun

          elseif (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            ! Lit1C + u Lit1N -> b Bacteria + f Fungi + c CO2 + n N
            ! c = min(CUEmax, Lit1N/Lit1C*CN_ratio_microbe)
            ! g = CNbacteria^0.6/(CNbacterial^0.6 + CNfungi^0.6)
            ! n = u - b nb - f nf
            ! b = g (1 - c)
            ! f = (1 - g) (1 - c)

            Rno3du_duc = -1.0d0 * this%upstream_nc(irxn) * rate_no3 / c_uc

            if (resp_frac < CUE_max) then
              ! Rdc/dLit1C = -RLit1N/Lit1C^2*CN_ratio_microbe           
              Rno3dc_duc = -1.0d0 * this%upstream_nc(irxn) * rate_no3 / c_uc &
                         * CN_ratio_microbe
            else
              Rno3dc_duc = 0.0d0
            endif

            ! Rdb/dLit1C = -g Rdc/dLit1C  
            Rno3db_duc = -1.0d0 * this%fraction_bacteria * Rno3dc_duc

            ! Rdf/dLit1C = -(1 - g) Rdc/dLit1C  
            Rno3df_duc = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rno3dc_duc

            ! Rdn/dLit1C = Rdu/dLit1C - nb Rdb/dLit1C - nf Rdf/dLit1C  
            Rno3dn_duc = Rno3du_duc - this%nc_bacteria * Rno3db_duc &
                       - this%nc_fungi * Rno3df_duc

            ! Rdu/dun = R /Lit1N 
            Rno3du_dun = rate_no3 / c_uc

            if (resp_frac < CUE_max) then
              ! Rdc/dLit1N = R/Lit1C*CN_ratio_microbe = dR/dLit1C*CN_ratio_microbe           
              Rno3dc_dun = rate_no3 / c_uc * CN_ratio_microbe
            else
              Rno3dc_dun = 0.0d0
            endif

            ! Rdb/dLit1N = -g Rdc/dLit1N 
            Rno3db_dun = -1.0d0 * this%fraction_bacteria * Rno3dc_dun

            ! Rdf/dLit1N = -(1 - g) Rdc/dLit1N  
            Rno3df_dun = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rno3dc_dun

            ! Rdn/dLit1N = Rdu/dLit1N - nb Rdb/dLit1N - nf Rdf/dLit1N  
            Rno3dn_dun = Rno3du_dun - this%nc_bacteria * Rno3db_dun &
                       - this%nc_fungi * Rno3df_dun

          endif
        endif

        ! with respect to upstream
        ! CO2

        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_no3_duc

        if (this%is_litter_decomp(irxn) .and. &
          this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! dRco2/dLitC = dcR/dLitC = cdR/dLitC + R dc/dLitC
          Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - Rno3dc_duc
        endif

        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - &
            this%mineral_c_stoich(irxn) * drate_no3_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE ) then 
            Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) &
                                         - Rno3dc_duc
          endif
        endif

        ! NO3
        ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
        ! first term: n dR/dC_u
        Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_no3_duc

        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + &
            this%mineral_n_stoich(irxn) * drate_no3_duc
        endif

        ! second term: R dn/dC_u
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - Rno3dn_duc

          if (this%species_id_nimm > 0) then
            Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) &
                                        + Rno3dn_duc
          endif

        endif

        ! upstream C pool
        Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) + drate_no3_duc

        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) + &
            this%upstream_nc(irxn) * drate_no3_duc + Rno3du_duc
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
         
          Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
            this%downstream_stoich(irxn, j) * drate_no3_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            if (ispec_d == this%species_id_bacteria) then
              Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rno3db_duc
            elseif (ispec_d == this%species_id_fungi) then
              Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rno3df_duc
            else
              option%io_buffer = 'Downstream pool for CLM-Microbe should ' // &
                                 'be either bacteria or fungi!'
            call printErrMsg(option)
            endif
          endif
        enddo

        ! with respect to upstream n (due to variable CN ratio)
        if (this%is_litter_decomp(irxn)) then

          Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) + Rno3du_dun

          Jacobian(ires_no3,ires_un) = Jacobian(ires_no3,ires_un) - Rno3dn_dun

          if (this%species_id_nimm > 0) then
            Jacobian(ires_nimm,ires_un) = Jacobian(ires_nimm,ires_un) &
                                        + Rno3dn_dun
          endif

          if (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            Jacobian(ires_co2,ires_un) = Jacobian(ires_co2,ires_un) - Rno3dc_dun

            if (this%species_id_hrimm > 0) then
              Jacobian(ires_hrimm,ires_un) = Jacobian(ires_hrimm,ires_un) &
                                           - Rno3dc_dun
            endif

            Jacobian(ires_b,ires_un) = Jacobian(ires_b,ires_un) - Rno3db_dun

            Jacobian(ires_f,ires_un) = Jacobian(ires_f,ires_un) - Rno3df_dun

          endif 
        endif

        ! with respect to no3
        ! CO2
        Jacobian(ires_co2,ires_no3) = Jacobian(ires_co2,ires_no3) - &
          this%mineral_c_stoich(irxn) * drate_no3_dno3
  
        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_no3) = Jacobian(ires_hrimm,ires_no3) - &
            this%mineral_c_stoich(irxn) * drate_no3_dno3
        endif

        ! N
        Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) - &
          this%mineral_n_stoich(irxn) * drate_no3_dno3
  
        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_no3) = Jacobian(ires_nimm,ires_no3) + &
            this%mineral_n_stoich(irxn) * drate_no3_dno3
        endif

        ! upstream C pool
        Jacobian(ires_uc,ires_no3) = Jacobian(ires_uc,ires_no3) - &
          (-1.d0) * drate_no3_dno3
 
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_no3) = Jacobian(ires_un,ires_no3) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_no3_dno3
        endif
  
        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)

          ispec_d = this%downstream_id(irxn, j)

          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          Jacobian(ires_d,ires_no3) = Jacobian(ires_d,ires_no3) - &
            this%downstream_stoich(irxn, j) * drate_no3_dno3
        enddo

        ! with respect to nh3 (due to nh3 inhibition on no3 immobilization)
        ! CO2
        Jacobian(ires_co2,ires_nh3) = Jacobian(ires_co2,ires_nh3) - &
          this%mineral_c_stoich(irxn) * drate_no3_dnh3

        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_nh3) = Jacobian(ires_hrimm,ires_nh3) - &
            this%mineral_c_stoich(irxn) * drate_no3_dnh3
        endif

        ! N
        Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - &
          this%mineral_n_stoich(irxn) * drate_no3_dnh3
  
        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_nh3) = Jacobian(ires_nimm,ires_nh3) + &
            this%mineral_n_stoich(irxn) * drate_no3_dnh3
        endif

        ! upstream C pool
        Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
          (-1.d0) * drate_no3_dnh3
  
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_no3_dnh3 
        endif
  
        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)

          ispec_d = this%downstream_id(irxn, j)

          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
            this%downstream_stoich(irxn, j) * drate_no3_dnh3
        enddo

      endif

    endif
 
  enddo

  if (this%species_id_n2o > 0) then

    f_t = 1.0d0
    f_w = 1.0d0
    f_ph = 1.0d0

    if (f_t > 1.0d-20 .and. f_w > 1.0d-20 .and. f_ph > 1.0d-20) then
      temp_real = f_t * f_w * f_ph

      if (temp_real > 1.0d0) then
        temp_real = 1.0d0
      endif

      temp_real = temp_real * this%n2o_frac_mineralization 
      
      if (net_n_mineralization_rate <= this%net_n_min_rate_smooth_0) then
        f_rate_n2o = 0.0d0
        df_rate_n2o = 0.0d0
      elseif (net_n_mineralization_rate >= this%net_n_min_rate_smooth_1 .or. &
        this%net_n_min_rate_smooth_1-this%net_n_min_rate_smooth_1 > 1.d-20) then
        f_rate_n2o = 1.0d0
        df_rate_n2o = 0.0d0
      else
        xxx = net_n_mineralization_rate - this%net_n_min_rate_smooth_0
        delta = this%net_n_min_rate_smooth_1 - this%net_n_min_rate_smooth_0
        f_rate_n2o  = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        df_rate_n2o = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                   / delta / delta
      endif

      ! residuals 
      rate_n2o = temp_real * net_n_mineralization_rate * f_nh3 * f_rate_n2o
 
      Residual(ires_nh3) = Residual(ires_nh3) + rate_n2o 

      Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o

      if (this%species_id_ngasmin > 0) then
         Residual(ires_ngasmin) = Residual(ires_ngasmin) - 0.5d0 * rate_n2o
      endif

      if (compute_derivative) then
        drate_n2o_dnh3 = temp_real * dnet_n_mineralization_rate_dnh3 * f_nh3 &
                       + temp_real * net_n_mineralization_rate * d_nh3
   
        drate_n2o_dnh3 = drate_n2o_dnh3 * f_rate_n2o + rate_n2o * df_rate_n2o &
                  * dnet_n_mineralization_rate_dnh3

        Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_n2o_dnh3
        Jacobian(ires_n2o,ires_nh3) = Jacobian(ires_n2o,ires_nh3) &
                                    - 0.5d0 * drate_n2o_dnh3
        if (this%species_id_ngasmin > 0) then
           Jacobian(ires_ngasmin,ires_nh3) = Jacobian(ires_ngasmin,ires_nh3) &
                                           - 0.5d0 * drate_n2o_dnh3
        endif

        if (this%species_id_no3 > 0) then
          drate_n2o_dno3 = temp_real * dnet_n_mineralization_rate_dno3 * f_nh3
   
          drate_n2o_dno3 = drate_n2o_dno3 * f_rate_n2o + rate_n2o * df_rate_n2o &
                  * dnet_n_mineralization_rate_dno3

          Jacobian(ires_n2o,ires_no3) = Jacobian(ires_n2o,ires_no3) &
                                    - 0.5d0 * drate_n2o_dno3

          if (this%species_id_ngasmin > 0) then
             Jacobian(ires_ngasmin,ires_no3) = Jacobian(ires_ngasmin,ires_no3) &
                                           - 0.5d0 * drate_n2o_dno3
          endif
        endif       

        do irxn = 1, this%nrxn
          ispec_uc = this%upstream_c_id(irxn)

          if (this%upstream_is_aqueous(irxn)) then
            ires_uc = ispec_uc
          else
            ires_uc = reaction%offset_immobile + ispec_uc
          endif
      
          drate_n2o_duc = temp_real * dnet_n_mineralization_rate_duc(irxn) * f_nh3
   
          drate_n2o_duc = drate_n2o_duc * f_rate_n2o + rate_n2o * df_rate_n2o &
                  * dnet_n_mineralization_rate_duc(irxn)

          Jacobian(ires_n2o,ires_uc) = Jacobian(ires_n2o,ires_uc) &
                                    - 0.5d0 * drate_n2o_duc

          if (this%species_id_ngasmin > 0) then
            Jacobian(ires_ngasmin,ires_uc) = Jacobian(ires_ngasmin,ires_uc) &
                                           - 0.5d0 * drate_n2o_duc
          endif
       
        enddo

      endif

    endif

  endif

end subroutine CLMDec_React

! **************************************************************************** !
!
! CLMDecDestroy: Destroys allocatable or pointer objects created in this module
!
! **************************************************************************** !
subroutine CLMDec_Destroy(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_clmdec_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(clmdec_reaction_type), pointer :: cur_reaction, prev_reaction
  
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    prev_pool => cur_pool
    cur_pool => cur_pool%next
    deallocate(prev_pool)
    nullify(prev_pool)
  enddo
  
  cur_reaction => this%reactions
  do
    if (.not.associated(cur_reaction)) exit

    cur_pool => cur_reaction%downstream_pools  
    do
      if (.not.associated(cur_pool)) exit
      prev_pool => cur_pool
      cur_pool => cur_pool%next
      deallocate(prev_pool)
      nullify(prev_pool)
    enddo

    prev_reaction => cur_reaction
    cur_reaction => cur_reaction%next

    deallocate(prev_reaction)
    nullify(prev_reaction)
  enddo
  
  call DeallocateArray(this%pool_nc_ratio)
  call DeallocateArray(this%rate_constant)
  call DeallocateArray(this%is_litter_decomp)
  call DeallocateArray(this%upstream_c_id)
  call DeallocateArray(this%upstream_n_id)
  call DeallocateArray(this%upstream_nc)
  call DeallocateArray(this%upstream_is_aqueous)
  call DeallocateArray(this%downstream_id)
  call DeallocateArray(this%downstream_stoich)
  call DeallocateArray(this%downstream_is_aqueous)
  call DeallocateArray(this%mineral_c_stoich) 
  call DeallocateArray(this%mineral_n_stoich) 
 
end subroutine CLMDec_Destroy

end module Reaction_Sandbox_CLMDec_class
