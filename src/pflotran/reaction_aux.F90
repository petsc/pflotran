module Reaction_Aux_module
  
  use Database_Aux_module
  use Mineral_Aux_module
  use Microbial_Aux_module
  use Immobile_Aux_module
  use Surface_Complexation_Aux_module
  
#ifdef SOLID_SOLUTION  
  use Solid_Solution_Aux_module
#endif

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  ! activity coefficients
  PetscInt, parameter, public :: ACT_COEF_FREQUENCY_OFF = 0
  PetscInt, parameter, public :: ACT_COEF_FREQUENCY_TIMESTEP = 1
  PetscInt, parameter, public :: ACT_COEF_FREQUENCY_NEWTON_ITER = 2
  PetscInt, parameter, public :: ACT_COEF_ALGORITHM_LAG = 3
  PetscInt, parameter, public :: ACT_COEF_ALGORITHM_NEWTON = 4
  PetscInt, parameter, public :: NO_BDOT = 5

  type, public :: species_idx_type
    PetscInt :: h2o_aq_id
    PetscInt :: h_ion_id
    PetscInt :: na_ion_id
    PetscInt :: cl_ion_id
    PetscInt :: co2_aq_id
    PetscInt :: tracer_aq_id
    PetscInt :: co2_gas_id
    PetscInt :: o2_gas_id
    PetscInt :: water_age_id
    PetscInt :: tracer_age_id
  end type species_idx_type
  
  type, public :: aq_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: a0
    PetscReal :: molar_weight
    PetscReal :: Z
    PetscBool :: print_me
    PetscBool :: is_redox
    type(database_rxn_type), pointer :: dbaserxn
    type(aq_species_type), pointer :: next
  end type aq_species_type

  type, public :: gas_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(gas_species_type), pointer :: next    
  end type gas_species_type

  type, public :: colloid_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: mobile_fraction
    PetscReal :: forward_rate
    PetscReal :: backward_rate
    PetscReal :: surface_area
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(colloid_type), pointer :: next
  end type colloid_type

  type, public :: ion_exchange_rxn_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: mineral_name
    type(ion_exchange_cation_type), pointer :: cation_list
    PetscReal :: CEC
    type (ion_exchange_rxn_type), pointer :: next
  end type ion_exchange_rxn_type

  type, public :: ion_exchange_cation_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: k
    type (ion_exchange_cation_type), pointer :: next
  end type ion_exchange_cation_type

  type, public :: kd_rxn_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: Kd
    PetscReal :: Langmuir_B
    PetscReal :: Freundlich_n
    type (kd_rxn_type), pointer :: next
  end type kd_rxn_type    

  type, public :: radioactive_decay_rxn_type
    PetscInt :: id
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: rate_constant
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(radioactive_decay_rxn_type), pointer :: next
  end type radioactive_decay_rxn_type

  type, public :: general_rxn_type
    PetscInt :: id
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: forward_rate
    PetscReal :: backward_rate
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(general_rxn_type), pointer :: next
  end type general_rxn_type

  type, public :: aq_species_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscReal, pointer :: basis_molarity(:)
    PetscInt, pointer :: constraint_type(:)
    PetscInt, pointer :: constraint_spec_id(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type aq_species_constraint_type

  type, public :: colloid_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc_mob(:)
    PetscReal, pointer :: constraint_conc_imb(:)
    PetscReal, pointer :: basis_conc_mob(:)
    PetscReal, pointer :: basis_conc_imb(:)
  end type colloid_constraint_type

  type, public :: reaction_type
    character(len=MAXSTRINGLENGTH) :: database_filename
    PetscBool :: use_full_geochemistry
    PetscBool :: use_log_formulation ! flag for solving for the change in the log of the concentration
    PetscBool :: check_update
    PetscBool :: print_all_species
    PetscBool :: print_all_primary_species
    PetscBool :: print_all_secondary_species
    PetscBool :: print_all_gas_species
    PetscBool :: print_pH
    PetscBool :: print_Eh
    PetscBool :: print_pe
    PetscBool :: print_O2
    PetscBool :: print_kd
    PetscBool :: print_total_sorb
    PetscBool :: print_total_sorb_mobile
    PetscBool :: print_colloid
    PetscBool :: print_act_coefs
    PetscBool :: print_total_component
    PetscBool :: print_free_ion
    PetscBool :: print_total_bulk ! total in aq and sorbed phases
    PetscBool :: initialize_with_molality
    PetscBool :: print_age
    PetscBool :: use_geothermal_hpt
    PetscInt :: print_free_conc_type
    PetscInt :: print_tot_conc_type
    PetscInt :: print_secondary_conc_type
    PetscInt :: num_dbase_temperatures
    PetscInt :: num_dbase_parameters
    PetscReal, pointer :: dbase_temperatures(:)
    type(species_idx_type), pointer :: species_idx

    type(aq_species_type), pointer :: primary_species_list
    type(aq_species_type), pointer :: secondary_species_list
    type(gas_species_type), pointer :: gas_species_list
    type(colloid_type), pointer :: colloid_list
    type(ion_exchange_rxn_type), pointer :: ion_exchange_rxn_list
    type(general_rxn_type), pointer :: general_rxn_list
    type(radioactive_decay_rxn_type), pointer :: radioactive_decay_rxn_list
    type(kd_rxn_type), pointer :: kd_rxn_list
    type(aq_species_type), pointer :: redox_species_list
    PetscInt :: act_coef_update_frequency
    PetscInt :: act_coef_update_algorithm
    PetscBool :: checkpoint_activity_coefs
    PetscBool :: act_coef_use_bdot
    PetscBool :: use_activity_h2o
    PetscBool :: calculate_water_age
    PetscBool :: calculate_tracer_age
    
    ! new reaction objects
    type(surface_complexation_type), pointer :: surface_complexation
    type(mineral_type), pointer :: mineral
    type(microbial_type), pointer :: microbial
    type(immobile_type), pointer :: immobile
    
#ifdef SOLID_SOLUTION    
    type(solid_solution_type), pointer :: solid_solution_list
#endif    
    
    ! compressed arrays for efficient computation
    ! primary aqueous complexes
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: ncollcomp
    PetscInt :: nimcomp
    
    ! offsets
    PetscInt :: offset_aqueous
    PetscInt :: offset_colloid
    PetscInt :: offset_collcomp
    PetscInt :: offset_immobile
    
    character(len=MAXWORDLENGTH), pointer :: primary_species_names(:)
    PetscBool, pointer :: primary_species_print(:)
    PetscReal, pointer :: primary_spec_a0(:)
    PetscReal, pointer :: primary_spec_Z(:)
    PetscReal, pointer :: primary_spec_molar_wt(:)
    
    ! aqueous complexes
    PetscInt :: neqcplx
    character(len=MAXWORDLENGTH), pointer :: secondary_species_names(:)
    PetscBool, pointer :: secondary_species_print(:)
    character(len=MAXWORDLENGTH), pointer :: eqcplx_basis_names(:,:)
    PetscBool, pointer :: eqcplx_basis_print(:)
    PetscInt, pointer :: eqcplxspecid(:,:)   ! (0:ncomp in rxn)
    PetscReal, pointer :: eqcplxstoich(:,:)
    PetscInt, pointer :: eqcplxh2oid(:)       ! id of water, if present
    PetscReal, pointer :: eqcplxh2ostoich(:)  ! stoichiometry of water, if present
    PetscReal, pointer :: eqcplx_a0(:)  ! Debye-Huckel constant
    PetscReal, pointer :: eqcplx_Z(:)
    PetscReal, pointer :: eqcplx_molar_wt(:)
    PetscReal, pointer :: eqcplx_logK(:)
    PetscReal, pointer :: eqcplx_logKcoef(:,:)
    ! Debye-Huckel
    PetscReal :: debyeA  ! Debye-Huckel A coefficient
    PetscReal :: debyeB  ! Debye-Huckel B coefficient
    PetscReal :: debyeBdot  ! Debye-Huckel Bdot coefficient
    
    ! gas species
    PetscInt :: ngas
    character(len=MAXWORDLENGTH), pointer :: gas_species_names(:)
    PetscBool, pointer :: gas_species_print(:)
    PetscInt, pointer :: eqgasspecid(:,:)   ! (0:ncomp in rxn)
    PetscReal, pointer :: eqgasstoich(:,:)
    PetscInt, pointer :: eqgash2oid(:)       ! id of water, if present
    PetscReal, pointer :: eqgash2ostoich(:)  ! stoichiometry of water, if present
    PetscReal, pointer :: eqgas_logK(:)
    PetscReal, pointer :: eqgas_logKcoef(:,:)
!#ifdef CHUAN_CO2
!   PetscReal :: scco2_eq_logK ! SC CO2 
!#endif

    PetscInt :: nsorb
    PetscInt :: neqsorb
    PetscBool, pointer :: kd_print(:)
    PetscBool, pointer :: total_sorb_print(:)

    ! ionx exchange reactions
    PetscInt :: neqionxrxn
    PetscInt :: neqionxcation 
    PetscBool, pointer :: eqionx_rxn_Z_flag(:)
    PetscInt, pointer :: eqionx_rxn_cation_X_offset(:)
    PetscReal, pointer :: eqionx_rxn_CEC(:)
    PetscInt, pointer :: eqionx_rxn_to_surf(:)
    PetscReal, pointer :: eqionx_rxn_k(:,:)
    PetscInt, pointer :: eqionx_rxn_cationid(:,:)
#if 0    
    PetscReal, pointer :: kinionx_rxn_CEC(:)
    PetscReal, pointer :: kinionx_rxn_k(:,:)
    PetscInt, pointer :: kinionx_rxn_cationid(:)
#endif    

    ! colloids
    PetscInt :: ncoll
    character(len=MAXWORDLENGTH), pointer :: colloid_names(:)
    character(len=MAXWORDLENGTH), pointer :: colloid_species_names(:)
    PetscReal, pointer :: colloid_mobile_fraction(:)
    PetscInt, pointer :: pri_spec_to_coll_spec(:)
    PetscInt, pointer :: coll_spec_to_pri_spec(:)
    PetscBool, pointer :: total_sorb_mobile_print(:)
    PetscBool, pointer :: colloid_print(:)
    
    ! radioactive decay rxn
    PetscInt :: nradiodecay_rxn
    ! ids and stoichiometries for species involved in reaction
    PetscInt, pointer :: radiodecayspecid(:,:)
    PetscReal, pointer :: radiodecaystoich(:,:)
    ! index of radiodecayspecid for species in forward
    ! reaction equation 
    PetscInt, pointer :: radiodecayforwardspecid(:)
    PetscReal, pointer :: radiodecay_kf(:)

    ! general rxn
    PetscInt :: ngeneral_rxn
    ! ids and stoichiometries for species involved in reaction
    PetscInt, pointer :: generalspecid(:,:)
    PetscReal, pointer :: generalstoich(:,:)
    ! index of generalspecid & generalstoich for species in forward
    ! reaction equation 
    PetscInt, pointer :: generalforwardspecid(:,:)
    PetscReal, pointer :: generalforwardstoich(:,:)
    ! index of generalspecid & generalstoich for species in backward
    ! reaction equation 
    PetscInt, pointer :: generalbackwardspecid(:,:)
    PetscReal, pointer :: generalbackwardstoich(:,:)
    PetscInt, pointer :: generalh2oid(:)
    PetscReal, pointer :: generalh2ostoich(:)
    PetscReal, pointer :: general_kf(:)
    PetscReal, pointer :: general_kr(:)  
    
    ! kd rxn
    PetscInt :: neqkdrxn
    PetscInt, pointer :: eqkdspecid(:)
    PetscInt, pointer :: eqkdtype(:)
    PetscReal, pointer :: eqkddistcoef(:)
    PetscReal, pointer :: eqkdlangmuirb(:)
    PetscReal, pointer :: eqkdfreundlichn(:)
    
    PetscReal :: max_dlnC
    PetscReal :: max_relative_change_tolerance
    PetscReal :: max_residual_tolerance
    
    PetscBool :: update_permeability
    PetscBool :: update_tortuosity
    PetscBool :: update_porosity
    PetscReal :: minimum_porosity
    PetscBool :: update_mineral_surface_area
    PetscBool :: update_mnrl_surf_with_porosity

    PetscBool :: update_armor_mineral_surface
    PetscInt :: update_armor_mineral_surface_flag


    PetscBool :: use_sandbox
    
  end type reaction_type

  public :: ReactionCreate, &
            SpeciesIndexCreate, &
            GasSpeciesCreate, &
            GetPrimarySpeciesCount, &
            GetPrimarySpeciesNames, &
            GetPrimarySpeciesIDFromName, &
            GetSecondarySpeciesCount, &
            GetSecondarySpeciesNames, &
            GetGasCount, &
            GetGasNames, &
            GetGasIDFromName, &
            GetMineralCount, &
            GetMineralNames, &
            GetMineralIDFromName, &
            GetColloidCount, &
            GetColloidNames, &
            GetColloidIDFromName, &
            GetImmobileCount, &
            ReactionFitLogKCoef, &
            ReactionInitializeLogK, &
            ReactionInterpolateLogK, &
            ReactionInitializeLogK_hpt, &
            ReactionInterpolateLogK_hpt, &
            TransitionStateTheoryRxnCreate, &
            TransitionStatePrefactorCreate, &
            TSPrefactorSpeciesCreate, &
            TransitionStateTheoryRxnDestroy, &
            AqueousSpeciesCreate, &
            AqueousSpeciesDestroy, &
            AqueousSpeciesConstraintCreate, &
            AqueousSpeciesConstraintDestroy, &
            MineralConstraintCreate, &
            MineralConstraintDestroy, &
            RadioactiveDecayRxnCreate, &
            RadioactiveDecayRxnDestroy, &
            GeneralRxnCreate, &
            GeneralRxnDestroy, &
            KDRxnCreate, &
            KDRxnDestroy, &
            ColloidCreate, &
            ColloidDestroy, &
            ColloidConstraintCreate, &
            ColloidConstraintDestroy, &
            IonExchangeRxnCreate, &
            IonExchangeCationCreate, &
            ReactionDestroy, &
            LogKeh
             
contains

! ************************************************************************** !
!
! ReactionCreate: Allocate and initialize reaction object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function ReactionCreate()

  use Option_module

  implicit none
  
  type(reaction_type), pointer :: ReactionCreate
  
  type(reaction_type), pointer :: reaction

  allocate(reaction)  

  reaction%database_filename = ''
  reaction%num_dbase_temperatures = 0
  nullify(reaction%dbase_temperatures)

  reaction%act_coef_use_bdot = PETSC_TRUE
  reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
  reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_LAG
  reaction%checkpoint_activity_coefs = PETSC_TRUE
  reaction%print_all_species = PETSC_FALSE
  reaction%print_all_primary_species = PETSC_FALSE
  reaction%print_all_secondary_species = PETSC_FALSE
  reaction%print_all_gas_species = PETSC_FALSE
  reaction%print_pH = PETSC_FALSE
  reaction%print_Eh = PETSC_FALSE
  reaction%print_pe = PETSC_FALSE
  reaction%print_O2 = PETSC_FALSE
  reaction%print_kd = PETSC_FALSE
  reaction%print_total_sorb = PETSC_FALSE
  reaction%print_total_sorb_mobile = PETSC_FALSE
  reaction%print_colloid = PETSC_FALSE
  reaction%print_act_coefs = PETSC_FALSE
  reaction%use_log_formulation = PETSC_FALSE
  reaction%check_update = PETSC_TRUE
  reaction%use_full_geochemistry = PETSC_FALSE
  reaction%use_activity_h2o = PETSC_FALSE
  reaction%calculate_tracer_age = PETSC_FALSE
  reaction%calculate_water_age = PETSC_FALSE
  reaction%print_age = PETSC_FALSE
  reaction%print_total_component = PETSC_TRUE
  reaction%print_free_ion = PETSC_FALSE
  reaction%print_total_bulk = PETSC_FALSE
  reaction%use_geothermal_hpt = PETSC_FALSE

  reaction%initialize_with_molality = PETSC_FALSE
  reaction%print_free_conc_type = 0
  reaction%print_tot_conc_type = 0
  reaction%print_secondary_conc_type = 0
  
  nullify(reaction%species_idx)

  nullify(reaction%primary_species_list)
  nullify(reaction%secondary_species_list)
  nullify(reaction%gas_species_list)
  nullify(reaction%colloid_list)
  nullify(reaction%ion_exchange_rxn_list)
  nullify(reaction%radioactive_decay_rxn_list)
  nullify(reaction%general_rxn_list)
  nullify(reaction%kd_rxn_list)
  nullify(reaction%redox_species_list)
  
  ! new reaction objects
  reaction%surface_complexation => SurfaceComplexationCreate()
  reaction%mineral => MineralCreate()
  reaction%microbial => MicrobialCreate()
  reaction%immobile => ImmobileCreate()
#ifdef SOLID_SOLUTION  
  nullify(reaction%solid_solution_list)
#endif
  
  nullify(reaction%primary_species_names)
  nullify(reaction%secondary_species_names)
  nullify(reaction%eqcplx_basis_names)
  nullify(reaction%gas_species_names)
  nullify(reaction%colloid_names)
  nullify(reaction%colloid_species_names)

  nullify(reaction%primary_species_print)
  nullify(reaction%secondary_species_print)
  nullify(reaction%eqcplx_basis_print)
  nullify(reaction%gas_species_print)
  nullify(reaction%kd_print)
  nullify(reaction%total_sorb_print)
  nullify(reaction%total_sorb_mobile_print)
  nullify(reaction%colloid_print)
  
  reaction%ncomp = 0
  reaction%naqcomp = 0
  reaction%ncoll = 0
  reaction%ncollcomp = 0
  reaction%nimcomp = 0
  reaction%offset_aqueous = 0
  reaction%offset_colloid = 0
  reaction%offset_collcomp = 0
  reaction%offset_immobile = 0
  nullify(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_Z)
  nullify(reaction%primary_spec_molar_wt)

  reaction%ngas = 0
  nullify(reaction%eqgasspecid)
  nullify(reaction%eqgasstoich)
  nullify(reaction%eqgash2oid)
  nullify(reaction%eqgash2ostoich)
  nullify(reaction%eqgas_logK)
  nullify(reaction%eqgas_logKcoef)
!#ifdef CHUAN_CO2
! reaction%scco2_eq_logK = 0.d0
!#endif
  
  reaction%neqcplx = 0
  nullify(reaction%eqcplxspecid)
  nullify(reaction%eqcplxstoich)
  nullify(reaction%eqcplxh2oid)
  nullify(reaction%eqcplxh2ostoich)
  nullify(reaction%eqcplx_a0)
  nullify(reaction%eqcplx_Z)
  nullify(reaction%eqcplx_molar_wt)
  nullify(reaction%eqcplx_logK)
  nullify(reaction%eqcplx_logKcoef)
  
  reaction%debyeA = 0.5114d0 
  reaction%debyeB = 0.3288d0 
  reaction%debyeBdot = 0.0410d0 

  reaction%nsorb = 0
  reaction%neqsorb = 0

  reaction%neqionxrxn = 0
  reaction%neqionxcation = 0
  nullify(reaction%eqionx_rxn_Z_flag)
  nullify(reaction%eqionx_rxn_cation_X_offset)
  nullify(reaction%eqionx_rxn_CEC)
  nullify(reaction%eqionx_rxn_to_surf)
  nullify(reaction%eqionx_rxn_k)
  nullify(reaction%eqionx_rxn_cationid)
#if 0  
  nullify(reaction%kinionx_CEC)
  nullify(reaction%kinionx_k)
  nullify(reaction%kinionx_cationid)
#endif

  reaction%ncoll = 0
  nullify(reaction%pri_spec_to_coll_spec)
  nullify(reaction%coll_spec_to_pri_spec)
  nullify(reaction%colloid_mobile_fraction)
  
  reaction%ngeneral_rxn = 0
  nullify(reaction%generalspecid)
  nullify(reaction%generalstoich)
  nullify(reaction%generalforwardspecid)
  nullify(reaction%generalforwardstoich)
  nullify(reaction%generalbackwardspecid)
  nullify(reaction%generalbackwardstoich)
  nullify(reaction%generalh2oid)
  nullify(reaction%generalh2ostoich)
  nullify(reaction%general_kf)
  nullify(reaction%general_kr)
  
  reaction%nradiodecay_rxn = 0
  nullify(reaction%radiodecayspecid)
  nullify(reaction%radiodecaystoich)
  nullify(reaction%radiodecayforwardspecid)
  nullify(reaction%radiodecay_kf)

  reaction%neqkdrxn = 0
  nullify(reaction%eqkdspecid)
  nullify(reaction%eqkdtype)
  nullify(reaction%eqkddistcoef)
  nullify(reaction%eqkdlangmuirb)
  nullify(reaction%eqkdfreundlichn)
      
  reaction%max_dlnC = 5.d0
  reaction%max_relative_change_tolerance = 1.d-6
  reaction%max_residual_tolerance = 1.d-12

  reaction%update_permeability = PETSC_FALSE
  reaction%update_tortuosity = PETSC_FALSE
  reaction%update_porosity = PETSC_FALSE
  reaction%minimum_porosity = 0.d0
  reaction%update_mineral_surface_area = PETSC_FALSE
  reaction%update_mnrl_surf_with_porosity = PETSC_FALSE

  reaction%update_armor_mineral_surface = PETSC_FALSE
  reaction%update_armor_mineral_surface_flag = 0

  reaction%use_sandbox = PETSC_FALSE

  ReactionCreate => reaction
  
end function ReactionCreate

! ************************************************************************** !
!
! SpeciesIndexCreate: Allocate and initialize a species index object
! author: Peter Lichtner
! date: 01/29/10
!
! ************************************************************************** !
function SpeciesIndexCreate()

  use Option_module

  implicit none
  
  type(species_idx_type), pointer :: SpeciesIndexCreate
  
  type(species_idx_type), pointer :: species_idx

  allocate(species_idx) 

  species_idx%h2o_aq_id = 0
  species_idx%h_ion_id = 0
  species_idx%na_ion_id = 0
  species_idx%cl_ion_id = 0
  species_idx%co2_aq_id = 0
  species_idx%tracer_aq_id = 0
  species_idx%co2_gas_id = 0
  species_idx%o2_gas_id = 0
  species_idx%tracer_age_id = 0
  species_idx%water_age_id = 0

  SpeciesIndexCreate => species_idx
  
end function SpeciesIndexCreate

! ************************************************************************** !
!
! AqueousSpeciesCreate: Allocate and initialize an aqueous species object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function AqueousSpeciesCreate()

  use Option_module

  implicit none
  
  type(aq_species_type), pointer :: AqueousSpeciesCreate
  
  type(aq_species_type), pointer :: species

  allocate(species) 
  species%id = 0 
  species%name = ''
  species%a0 = 0.d0
  species%molar_weight = 0.d0
  species%Z = 0.d0
  species%print_me = PETSC_FALSE
  species%is_redox = PETSC_FALSE
  nullify(species%dbaserxn)
  nullify(species%next)

  AqueousSpeciesCreate => species
  
end function AqueousSpeciesCreate

! ************************************************************************** !
!
! GasSpeciesCreate: Allocate and initialize a gas species object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function GasSpeciesCreate()

  use Option_module

  implicit none
  
  type(gas_species_type), pointer :: GasSpeciesCreate
  
  type(gas_species_type), pointer :: species

  allocate(species)  
  species%id = 0
  species%name = ''
  species%molar_volume = 0.d0
  species%molar_weight = 0.d0
  species%print_me = PETSC_FALSE
  nullify(species%dbaserxn)
  nullify(species%next)

  GasSpeciesCreate => species
  
end function GasSpeciesCreate

! ************************************************************************** !
!
! ColloidCreate: Allocate and initialize a colloid object
! author: Glenn Hammond
! date: 02/24/10
!
! ************************************************************************** !
function ColloidCreate()

  use Option_module

  implicit none
  
  type(colloid_type), pointer :: ColloidCreate
  
  type(colloid_type), pointer :: colloid

  allocate(colloid)  
  colloid%id = 0
  colloid%itype = 0
  colloid%name = ''
  colloid%mobile_fraction = 0.5d0
  colloid%forward_rate = 0.d0
  colloid%backward_rate = 0.d0
  colloid%surface_area = 1.d0
  colloid%molar_weight = 0.d0
  colloid%print_me = PETSC_FALSE
  nullify(colloid%next)
  
  ColloidCreate => colloid
  
end function ColloidCreate

! ************************************************************************** !
!
! IonExchangeRxnCreate: Allocate and initialize an ion exchange reaction
! author: Peter Lichtner
! date: 10/24/08
!
! ************************************************************************** !
function IonExchangeRxnCreate()

  implicit none
    
  type(ion_exchange_rxn_type), pointer :: IonExchangeRxnCreate

  type(ion_exchange_rxn_type), pointer :: ionxrxn
  
  allocate(ionxrxn)
  ionxrxn%id = 0
  ionxrxn%mineral_name = ''
  ionxrxn%CEC = 0.d0
  nullify(ionxrxn%cation_list)
  nullify(ionxrxn%next)
  
  IonExchangeRxnCreate => ionxrxn
  
end function IonExchangeRxnCreate

! ************************************************************************** !
!
! IonExchangeCationCreate: Allocate and initialize a cation associated with
!                          an ion exchange reaction
! author: Peter Lichtner
! date: 10/24/08
!
! ************************************************************************** !
function IonExchangeCationCreate()

  implicit none
    
  type(ion_exchange_cation_type), pointer :: IonExchangeCationCreate

  type(ion_exchange_cation_type), pointer :: cation
  
  allocate(cation)
  cation%name = ''
  cation%k = 0.d0
  nullify(cation%next)
  
  IonExchangeCationCreate => cation
  
end function IonExchangeCationCreate

! ************************************************************************** !
!
! RadioactiveDecayRxnCreate: Allocate and initialize a radioactive decay
!                            reaction
! author: Glenn Hammond
! date: 01/07/14
!
! ************************************************************************** !
function RadioactiveDecayRxnCreate()

  implicit none
    
  type(radioactive_decay_rxn_type), pointer :: RadioactiveDecayRxnCreate

  type(radioactive_decay_rxn_type), pointer :: rxn
  
  allocate(rxn)
  rxn%id = 0
  rxn%reaction = ''
  rxn%rate_constant = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%dbaserxn)
  nullify(rxn%next)
  
  RadioactiveDecayRxnCreate => rxn
  
end function RadioactiveDecayRxnCreate

! ************************************************************************** !
!
! GeneralRxnCreate: Allocate and initialize a general reaction
! author: Glenn Hammond
! date: 09/03/10
!
! ************************************************************************** !
function GeneralRxnCreate()

  implicit none
    
  type(general_rxn_type), pointer :: GeneralRxnCreate

  type(general_rxn_type), pointer :: rxn
  
  allocate(rxn)
  rxn%id = 0
  rxn%reaction = ''
  rxn%forward_rate = 0.d0
  rxn%backward_rate = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%dbaserxn)
  nullify(rxn%next)
  
  GeneralRxnCreate => rxn
  
end function GeneralRxnCreate

! ************************************************************************** !
!
! KDRxnCreate: Allocate and initialize a KD sorption reaction
! author: Glenn Hammond
! date: 09/32/10
!
! ************************************************************************** !
function KDRxnCreate()

  implicit none
    
  type(kd_rxn_type), pointer :: KDRxnCreate

  type(kd_rxn_type), pointer :: rxn
  
  allocate(rxn)
  rxn%id = 0
  rxn%itype = 0
  rxn%species_name = ''
  rxn%Kd = 0.d0
  rxn%Langmuir_B = 0.d0
  rxn%Freundlich_n = 0.d0
  nullify(rxn%next)
  
  KDRxnCreate => rxn
  
end function KDRxnCreate

! ************************************************************************** !
!
! AqueousSpeciesConstraintCreate: Creates an aqueous species constraint 
!                                 object
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function AqueousSpeciesConstraintCreate(reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  type(aq_species_constraint_type), pointer :: AqueousSpeciesConstraintCreate

  type(aq_species_constraint_type), pointer :: constraint
  
  allocate(constraint)
  allocate(constraint%names(reaction%naqcomp))
  constraint%names = ''
  allocate(constraint%constraint_conc(reaction%naqcomp))
  constraint%constraint_conc = 0.d0
  allocate(constraint%basis_molarity(reaction%naqcomp))
  constraint%basis_molarity = 0.d0
  allocate(constraint%constraint_spec_id(reaction%naqcomp))
  constraint%constraint_spec_id = 0
  allocate(constraint%constraint_type(reaction%naqcomp))
  constraint%constraint_type = 0
  allocate(constraint%constraint_aux_string(reaction%naqcomp))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(reaction%naqcomp))
  constraint%external_dataset = PETSC_FALSE

  AqueousSpeciesConstraintCreate => constraint

end function AqueousSpeciesConstraintCreate
 
! ************************************************************************** !
!
! ColloidConstraintCreate: Creates a colloid constraint object
! author: Glenn Hammond
! date: 03/12/10
!
! ************************************************************************** !
function ColloidConstraintCreate(reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  type(colloid_constraint_type), pointer :: ColloidConstraintCreate

  type(colloid_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(reaction%ncoll))
  constraint%names = ''
  allocate(constraint%constraint_conc_mob(reaction%ncoll))
  constraint%constraint_conc_mob = 0.d0
  allocate(constraint%constraint_conc_imb(reaction%ncoll))
  constraint%constraint_conc_imb = 0.d0
  allocate(constraint%basis_conc_mob(reaction%ncoll))
  constraint%basis_conc_mob = 0.d0
  allocate(constraint%basis_conc_imb(reaction%ncoll))
  constraint%basis_conc_imb = 0.d0

  ColloidConstraintCreate => constraint

end function ColloidConstraintCreate
  
! ************************************************************************** !
!
! GetPrimarySpeciesNames: Returns the names of primary species in an array
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetPrimarySpeciesNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetPrimarySpeciesNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = GetPrimarySpeciesCount(reaction)
  allocate(names(count))
  
  count = 1
  species => reaction%primary_species_list
  do
    if (.not.associated(species)) exit
    names(count) = species%name
    count = count + 1
    species => species%next
  enddo

  GetPrimarySpeciesNames => names
  
end function GetPrimarySpeciesNames

! ************************************************************************** !
!
! GetPrimarySpeciesCount: Returns the number of primary species
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetPrimarySpeciesCount(reaction)

  implicit none
  
  PetscInt :: GetPrimarySpeciesCount
  type(reaction_type) :: reaction

  type(aq_species_type), pointer :: species

  GetPrimarySpeciesCount = 0
  species => reaction%primary_species_list
  do
    if (.not.associated(species)) exit
    GetPrimarySpeciesCount = GetPrimarySpeciesCount + 1
    species => species%next
  enddo

end function GetPrimarySpeciesCount

! ************************************************************************** !
!
! GetPrimarySpeciesIDFromName: Returns the id of named primary species
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
function GetPrimarySpeciesIDFromName(name,reaction,option)

  use Option_module
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(reaction_type) :: reaction
  type(option_type) :: option

  PetscInt :: GetPrimarySpeciesIDFromName

  type(aq_species_type), pointer :: species
  PetscInt :: i

  GetPrimarySpeciesIDFromName = -999
  
  ! if the primary species name list exists
  if (associated(reaction%primary_species_names)) then
    do i = 1, size(reaction%primary_species_names)
      if (StringCompare(name,reaction%primary_species_names(i), &
                        MAXWORDLENGTH)) then
        GetPrimarySpeciesIDFromName = i
        exit
      endif
    enddo
  else
    species => reaction%primary_species_list
    i = 0
    do
      if (.not.associated(species)) exit
      i = i + 1
      if (StringCompare(name,species%name,MAXWORDLENGTH)) then
        GetPrimarySpeciesIDFromName = i
        exit
      endif
      species => species%next
    enddo
  endif

  if (GetPrimarySpeciesIDFromName <= 0) then
    option%io_buffer = 'Species "' // trim(name) // &
      '" not founds among primary species in GetPrimarySpeciesIDFromName().'
    call printErrMsg(option)
  endif
  
end function GetPrimarySpeciesIDFromName

! ************************************************************************** !
!
! GetSecondarySpeciesNames: Returns the names of secondary species in an array
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetSecondarySpeciesNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetSecondarySpeciesNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = GetSecondarySpeciesCount(reaction)
  allocate(names(count))
  
  count = 1
  species => reaction%secondary_species_list
  do
    if (.not.associated(species)) exit
    names(count) = species%name
    count = count + 1
    species => species%next
  enddo

  GetSecondarySpeciesNames => names
  
end function GetSecondarySpeciesNames

! ************************************************************************** !
!
! GetSecondarySpeciesCount: Returns the number of secondary species
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetSecondarySpeciesCount(reaction)

  implicit none
  
  PetscInt :: GetSecondarySpeciesCount
  type(reaction_type) :: reaction

  type(aq_species_type), pointer :: species

  GetSecondarySpeciesCount = 0
  species => reaction%secondary_species_list
  do
    if (.not.associated(species)) exit
    GetSecondarySpeciesCount = GetSecondarySpeciesCount + 1
    species => species%next
  enddo

end function GetSecondarySpeciesCount

! ************************************************************************** !
!
! GetGasNames: Returns the names of gases in an array
! author: Glenn Hammond
! date: 10/21/08
!
! ************************************************************************** !
function GetGasNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetGasNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(gas_species_type), pointer :: gas

  count = GetGasCount(reaction)
  allocate(names(count))
  
  count = 1
  gas => reaction%gas_species_list
  do
    if (.not.associated(gas)) exit
    names(count) = gas%name
    count = count + 1
    gas => gas%next
  enddo

  GetGasNames => names
  
end function GetGasNames

! ************************************************************************** !
!
! GetGasCount: Returns the number of primary species
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetGasCount(reaction)

  implicit none
  
  PetscInt :: GetGasCount
  type(reaction_type) :: reaction

  type(gas_species_type), pointer :: gas

  GetGasCount = 0
  gas => reaction%gas_species_list
  do
    if (.not.associated(gas)) exit
    GetGasCount = GetGasCount + 1
    gas => gas%next
  enddo

end function GetGasCount

! ************************************************************************** !
!
! GetGasIDFromName: Returns the id of gas with the corresponding name
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
function GetGasIDFromName(reaction,name)

  use String_module
  
  implicit none
  
  type(reaction_type) :: reaction
  character(len=MAXWORDLENGTH) :: name

  PetscInt :: GetGasIDFromName
  type(gas_species_type), pointer :: gas

  GetGasIDFromName = -1
 
  gas => reaction%gas_species_list
  do
    if (.not.associated(gas)) exit
    if (StringCompare(name,gas%name,MAXWORDLENGTH)) then
      GetGasIDFromName = gas%id
      exit
    endif
    gas => gas%next
  enddo

end function GetGasIDFromName

! ************************************************************************** !
!
! GetColloidIDFromName: Returns the id of colloid with the corresponding name
! author: Glenn Hammond
! date: 02/24/10
!
! ************************************************************************** !
function GetColloidIDFromName(reaction,name)

  use String_module
  
  implicit none
  
  type(reaction_type) :: reaction
  character(len=MAXWORDLENGTH) :: name

  PetscInt :: GetColloidIDFromName
  type(colloid_type), pointer :: colloid

  GetColloidIDFromName = -1
 
  colloid => reaction%colloid_list
  do
    if (.not.associated(colloid)) exit
    if (StringCompare(name,colloid%name,MAXWORDLENGTH)) then
      GetColloidIDFromName = colloid%id
      exit
    endif
    colloid => colloid%next
  enddo

end function GetColloidIDFromName

! ************************************************************************** !
!
! GetColloidNames: Returns the names of colloids in an array
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
function GetColloidNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetColloidNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(colloid_type), pointer :: colloid

  count = GetColloidCount(reaction)
  allocate(names(count))
  
  count = 1
  colloid => reaction%colloid_list
  do
    if (.not.associated(colloid)) exit
    names(count) = colloid%name
    count = count + 1
    colloid => colloid%next
  enddo

  GetColloidNames => names
  
end function GetColloidNames

! ************************************************************************** !
!
! GetColloidCount: Returns the number of colloids
! author: Glenn Hammond
! date: 02/24/10
!
! ************************************************************************** !
function GetColloidCount(reaction)

  implicit none
  
  PetscInt :: GetColloidCount
  type(reaction_type) :: reaction

  type(colloid_type), pointer :: colloid

  GetColloidCount = 0
  colloid => reaction%colloid_list
  do
    if (.not.associated(colloid)) exit
    GetColloidCount = GetColloidCount + 1
    colloid => colloid%next
  enddo

end function GetColloidCount


! ************************************************************************** !
!
! GetImmobileCount: Returns the number of immobile species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function GetImmobileCount(reaction)

  implicit none
  
  PetscInt :: GetImmobileCount
  type(reaction_type) :: reaction

  GetImmobileCount = ImmobileGetCount(reaction%immobile)
  
end function GetImmobileCount

! ************************************************************************** !
!
! ReactionFitLogKCoef: Least squares fit to log K over database temperature 
!                      range
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionFitLogKCoef(coefs,logK,name,option,reaction)

  use Option_module
  use Utility_module

  implicit none
  
  type(reaction_type) :: reaction
  PetscReal :: coefs(FIVE_INTEGER)
  character(len=MAXWORDLENGTH) :: name 
  PetscReal :: logK(reaction%num_dbase_temperatures)
  type(option_type) :: option

  PetscInt :: temp_int(reaction%num_dbase_temperatures), &
              indx(reaction%num_dbase_temperatures)
  PetscReal :: a(FIVE_INTEGER,FIVE_INTEGER), &
               vec(FIVE_INTEGER,reaction%num_dbase_temperatures), temperature_kelvin

  PetscInt :: i, j, k, iflag
  
  ! need to fill in vec with equations for temperatures vs coefs.
  
  do i = 1, reaction%num_dbase_temperatures
    temperature_kelvin = reaction%dbase_temperatures(i) + 273.15d0
    vec(1,i) = log(temperature_kelvin)
    vec(2,i) = 1.d0
    vec(3,i) = temperature_kelvin
    vec(4,i) = 1.d0/temperature_kelvin
    vec(5,i) = 1.d0/(temperature_kelvin*temperature_kelvin)
  enddo
  
  iflag = 0
  do j = 1, FIVE_INTEGER
    coefs(j) = 0.d0
    do i = 1, reaction%num_dbase_temperatures
      if (dabs(logK(i) - 500.) < 1.d-10) then
        iflag = 1
        temp_int(i) = ZERO_INTEGER
      else
        coefs(j) = coefs(j) + vec(j,i)*logK(i)
        temp_int(i) = ONE_INTEGER
      endif
    enddo
  enddo

  if (iflag == 1) then
    option%io_buffer = 'In ReactionFitLogKCoef: log K = 500 for ' // trim(name)
    call printWrnMsg(option)
  endif

  do j = 1, FIVE_INTEGER
    do k = j, FIVE_INTEGER
      a(j,k) = 0.d0
      do i = 1, reaction%num_dbase_temperatures
        if (temp_int(i) == 1) then
          a(j,k) = a(j,k) + vec(j,i)*vec(k,i)
        endif
      enddo
      if (j .ne. k) a(k,j) = a(j,k)
    enddo
  enddo

  call ludcmp(a,FIVE_INTEGER,indx,i)
  call lubksb(a,FIVE_INTEGER,indx,coefs)

end subroutine ReactionFitLogKCoef

! ************************************************************************** !
!
! ReactionInitializeLogK: Least squares fit to log K over database temperature range
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionInitializeLogK(logKcoef,logKs,logK,option,reaction)

  use Option_module

  implicit none
  
  type(reaction_type) :: reaction
  PetscReal :: logKcoef(FIVE_INTEGER)
  PetscReal :: logKs(reaction%num_dbase_temperatures)
  PetscReal :: logK, logK_1D_Array(ONE_INTEGER)
  type(option_type) :: option
  
  PetscReal :: coefs(FIVE_INTEGER,ONE_INTEGER)
  PetscReal :: temperature
  PetscInt :: itemperature
  PetscInt :: i
  
  ! we always initialize on reference temperature
  temperature = option%reference_temperature
  
  itemperature = 0
  if (option%use_isothermal) then ! find database temperature if relevant
    do i = 1, reaction%num_dbase_temperatures
      if (dabs(option%reference_temperature - &
               reaction%dbase_temperatures(i)) < 1.d-10) then
        itemperature = i
        exit
      endif
    enddo
  endif
  
  if (itemperature > 0) then ! use database temperature
    logK = logKs(itemperature)
  else                       ! interpolate
    coefs(:,ONE_INTEGER) = logKcoef(:)
    call ReactionInterpolateLogK(coefs,logK_1D_Array,temperature,ONE_INTEGER)
    logK = logK_1D_Array(ONE_INTEGER)
  endif

end subroutine ReactionInitializeLogK

! ************************************************************************** !
!
! ReactionInterpolateLogK: Interpolation log K function: temp - temperature [C]
!                             b - fit coefficients determined from fit(...)
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionInterpolateLogK(coefs,logKs,temp,n)

  implicit none
  
  PetscInt :: n
  PetscReal :: coefs(5,n), logKs(n), temp

  PetscInt :: i
  PetscReal :: temp_kelvin
  
  temp_kelvin = temp + 273.15d0
  
  do i = 1, n
    logKs(i) = coefs(1,i)*log(temp_kelvin) &
             + coefs(2,i)           &
             + coefs(3,i)*temp_kelvin      &
             + coefs(4,i)/temp_kelvin      &
             + coefs(5,i)/(temp_kelvin*temp_kelvin)
  enddo
  
end subroutine ReactionInterpolateLogK

! ************************************************************************** !
!
! ReactionInitializeLogK: Least squares fit to log K over database temperature range
! author: Chuan Lu
! date: 12/29/11
!
! ************************************************************************** !
subroutine ReactionInitializeLogK_hpt(logKcoef,logK,option,reaction)

  use Option_module

  implicit none
  
  type(reaction_type) :: reaction
  PetscReal :: logKcoef(17)
  PetscReal :: logK, logK_1D_Array(ONE_INTEGER)
  type(option_type) :: option
  
  PetscReal :: coefs(17,ONE_INTEGER)
  PetscReal :: temperature, pressure
  PetscInt :: itemperature
  PetscInt :: i
  
  ! we always initialize on reference temperature
  temperature = option%reference_temperature
  pressure = option%reference_pressure 
  
  
  coefs(:,ONE_INTEGER) = logKcoef(:)
  call ReactionInterpolateLogK_hpt(coefs,logK_1D_Array,temperature,pressure, &
                                   ONE_INTEGER)
  logK = logK_1D_Array(ONE_INTEGER)
!   print *,'ReactionInitializeLogK_hpt: ', pressure,temperature, logK

end subroutine ReactionInitializeLogK_hpt

! ************************************************************************** !
!
! ReactionInterpolateLogK: Interpolation log K function: temp - temperature [C]
!                             b - fit coefficients determined from fit(...)
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionInterpolateLogK_hpt(coefs,logKs,temp,pres,n)

  implicit none
  
  PetscInt :: n
  PetscReal :: coefs(17,n), logKs(n), temp, pres

  PetscInt :: i
  PetscReal :: temp_kelvin, tr, pr, logtr
  
  temp_kelvin = temp + 273.15d0
  tr = temp_kelvin/273.15d0
  pr = pres/1.d7
  logtr = log(tr)/log(10.d0) 
  
  do i = 1, n
    logKs(i) = coefs(1,i)                 &
             + coefs(2,i) * tr            &
             + coefs(3,i) / tr            &
             + coefs(4,i) * logtr         &
             + coefs(5,i) * tr * tr       &
             + coefs(6,i) / tr / tr       &
             + coefs(7,i) * sqrt(tr)      &
             + coefs(8,i) * pr            &
             + coefs(9,i) * pr * tr       &
             + coefs(10,i) * pr / tr      &
             + coefs(11,i) * pr * logtr   &
             + coefs(12,i) / pr           &
             + coefs(13,i) / pr * tr      &
             + coefs(14,i) / pr / tr      &
             + coefs(15,i) * pr * pr      &
             + coefs(16,i) * pr * pr * tr &
             + coefs(17,i) * pr * pr / tr 
  enddo
 ! print *,'ReactionInterpolateLogK_hpt: ', pres,temp, logKs, coefs
end subroutine ReactionInterpolateLogK_hpt


! ************************************************************************** !
!
! Function logkeh: Maier-Kelly fit to equilibrium constant half-cell reaction
! 2 H2O - 4 H+ - 4 e- = O2, to compute Eh and pe.
! author: Peter Lichtner
! date: 04/27/13
!
! ************************************************************************** !
PetscReal function logkeh(tk)

  implicit none

  PetscReal, intent(in) :: tk

  PetscReal, parameter :: cm1 = 6.745529048112373d0
  PetscReal, parameter :: c0 = -48.295936593543715d0
  PetscReal, parameter :: c1 = 0.0005578156078778505d0
  PetscReal, parameter :: c2 = 27780.749538022003d0
  PetscReal, parameter :: c3 = 4027.3376948579394d0

  logkeh = cm1 * log(tk) + c0 + c1 * tk + c2 / tk + c3 / (tk * tk)

end function logkeh

! ************************************************************************** !
!
! SpeciesIndexDestroy: Deallocates a species index object
! author: Glenn Hammond
! date: 01/29/10
!
! ************************************************************************** !
subroutine SpeciesIndexDestroy(species_idx)

  implicit none
    
  type(species_idx_type), pointer :: species_idx

  if (associated(species_idx)) deallocate(species_idx)  
  nullify(species_idx)

end subroutine SpeciesIndexDestroy

! ************************************************************************** !
!
! AqueousSpeciesDestroy: Deallocates an aqueous species
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine AqueousSpeciesDestroy(species)

  implicit none
    
  type(aq_species_type), pointer :: species

  if (associated(species%dbaserxn)) call DatabaseRxnDestroy(species%dbaserxn)
  deallocate(species)  
  nullify(species)

end subroutine AqueousSpeciesDestroy

! ************************************************************************** !
!
! AqueousSpeciesListDestroy: Deallocates an aqueous species
! author: Glenn Hammond
! date: 09/03/10
!
! ************************************************************************** !
subroutine AqueousSpeciesListDestroy(aq_species_list)

  implicit none
    
  type(aq_species_type), pointer :: aq_species_list  
    
  type(aq_species_type), pointer :: species, prev_species

  species => aq_species_list
  do
    if (.not.associated(species)) exit
    prev_species => species
    species => species%next
    call AqueousSpeciesDestroy(prev_species)
  enddo  
  nullify(aq_species_list)

end subroutine AqueousSpeciesListDestroy

! ************************************************************************** !
!
! GasSpeciesDestroy: Deallocates a gas species
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine GasSpeciesDestroy(species)

  implicit none
    
  type(gas_species_type), pointer :: species

  if (associated(species%dbaserxn)) call DatabaseRxnDestroy(species%dbaserxn)
  deallocate(species)  
  nullify(species)

end subroutine GasSpeciesDestroy

! ************************************************************************** !
!
! ColloidDestroy: Deallocates a colloid
! author: Glenn Hammond
! date: 02/24/10
!
! ************************************************************************** !
subroutine ColloidDestroy(colloid)

  implicit none
    
  type(colloid_type), pointer :: colloid

  deallocate(colloid)  
  nullify(colloid)

end subroutine ColloidDestroy

! ************************************************************************** !
!
! IonExchangeRxnDestroy: Deallocates an ion exchange reaction
! author: Glenn Hammond
! date: 10/24/08
!
! ************************************************************************** !
subroutine IonExchangeRxnDestroy(ionxrxn)

  implicit none
    
  type(ion_exchange_rxn_type), pointer :: ionxrxn
  
  type(ion_exchange_cation_type), pointer :: cur_cation, prev_cation

  if (.not.associated(ionxrxn)) return
  
  cur_cation => ionxrxn%cation_list
  do
    if (.not.associated(cur_cation)) exit
    prev_cation => cur_cation
    cur_cation => cur_cation%next
    deallocate(prev_cation)
    nullify(prev_cation)
  enddo
  
  nullify(ionxrxn%next)

  deallocate(ionxrxn)  
  nullify(ionxrxn)

end subroutine IonExchangeRxnDestroy

! ************************************************************************** !
!
! RadioactiveDecayRxnDestroy: Deallocates a general reaction
! author: Glenn Hammond
! date: 01/07/14
!
! ************************************************************************** !
subroutine RadioactiveDecayRxnDestroy(rxn)

  implicit none
    
  type(radioactive_decay_rxn_type), pointer :: rxn

  if (.not.associated(rxn)) return
  
  if (associated(rxn%dbaserxn)) &
    call DatabaseRxnDestroy(rxn%dbaserxn)
  nullify(rxn%dbaserxn)
  nullify(rxn%next)

  deallocate(rxn)  
  nullify(rxn)

end subroutine RadioactiveDecayRxnDestroy

! ************************************************************************** !
!
! GeneralRxnDestroy: Deallocates a general reaction
! author: Glenn Hammond
! date: 09/03/10
!
! ************************************************************************** !
subroutine GeneralRxnDestroy(rxn)

  implicit none
    
  type(general_rxn_type), pointer :: rxn

  if (.not.associated(rxn)) return
  
  if (associated(rxn%dbaserxn)) &
    call DatabaseRxnDestroy(rxn%dbaserxn)
  nullify(rxn%dbaserxn)
  nullify(rxn%next)

  deallocate(rxn)  
  nullify(rxn)

end subroutine GeneralRxnDestroy

! ************************************************************************** !
!
! KDRxnDestroy: Deallocates a KD reaction
! author: Glenn Hammond
! date: 09/30/10
!
! ************************************************************************** !
subroutine KDRxnDestroy(rxn)

  implicit none
    
  type(kd_rxn_type), pointer :: rxn

  if (.not.associated(rxn)) return
  
  deallocate(rxn)  
  nullify(rxn)

end subroutine KDRxnDestroy

! ************************************************************************** !
!
! AqueousSpeciesConstraintDestroy: Destroys an aqueous species constraint 
!                                  object
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine AqueousSpeciesConstraintDestroy(constraint)

  use Utility_module, only: DeallocateArray
  
  implicit none
  
  type(aq_species_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%basis_molarity)
  call DeallocateArray(constraint%constraint_type)
  call DeallocateArray(constraint%constraint_spec_id)
  call DeallocateArray(constraint%constraint_aux_string)
  call DeallocateArray(constraint%external_dataset)

  deallocate(constraint)
  nullify(constraint)

end subroutine AqueousSpeciesConstraintDestroy

! ************************************************************************** !
!
! ColloidConstraintDestroy: Destroys a colloid constraint object
! author: Glenn Hammond
! date: 03/12/10
!
! ************************************************************************** !
subroutine ColloidConstraintDestroy(constraint)

  use Utility_module, only: DeallocateArray

  implicit none
  
  type(colloid_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc_mob)
  call DeallocateArray(constraint%constraint_conc_imb)
  call DeallocateArray(constraint%basis_conc_mob)
  call DeallocateArray(constraint%basis_conc_imb)

  deallocate(constraint)
  nullify(constraint)

end subroutine ColloidConstraintDestroy

! ************************************************************************** !
!
! ReactionDestroy: Deallocates a reaction object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine ReactionDestroy(reaction)

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(reaction_type), pointer :: reaction
  
  type(aq_species_type), pointer :: aq_species, prev_aq_species
  type(gas_species_type), pointer :: gas_species, prev_gas_species
  type(mineral_rxn_type), pointer :: mineral, prev_mineral
  type(colloid_type), pointer :: colloid, prev_colloid
  type(ion_exchange_rxn_type), pointer :: ionxrxn, prev_ionxrxn
  type(surface_complexation_rxn_type), pointer :: srfcplxrxn, prev_srfcplxrxn
  type(general_rxn_type), pointer :: general_rxn, prev_general_rxn
  type(radioactive_decay_rxn_type), pointer :: radioactive_decay_rxn, &
                                               prev_radioactive_decay_rxn
  type(kd_rxn_type), pointer :: kd_rxn, prev_kd_rxn

  if (.not.associated(reaction)) return
  
  !species index
  call SpeciesIndexDestroy(reaction%species_idx)

  ! primary species
  if (associated(reaction%primary_species_list)) &
    call AqueousSpeciesListDestroy(reaction%primary_species_list)
  nullify(reaction%primary_species_list)

  ! secondary species
  if (associated(reaction%secondary_species_list)) &
    call AqueousSpeciesListDestroy(reaction%secondary_species_list)
  nullify(reaction%secondary_species_list)

  ! gas species
  gas_species => reaction%gas_species_list
  do
    if (.not.associated(gas_species)) exit
    prev_gas_species => gas_species
    gas_species => gas_species%next
    call GasSpeciesDestroy(prev_gas_species)
  enddo  
  nullify(reaction%gas_species_list)
  
  ! colloid species
  colloid => reaction%colloid_list
  do
    if (.not.associated(colloid)) exit
    prev_colloid => colloid
    colloid => colloid%next
    call ColloidDestroy(prev_colloid)
  enddo    
  nullify(reaction%colloid_list)
  
  ! ionx exchange reactions
  ionxrxn => reaction%ion_exchange_rxn_list
  do
    if (.not.associated(ionxrxn)) exit
    prev_ionxrxn => ionxrxn
    ionxrxn => ionxrxn%next
    call IonExchangeRxnDestroy(prev_ionxrxn)
  enddo    
  nullify(reaction%ion_exchange_rxn_list)

  ! general reactions
  radioactive_decay_rxn => reaction%radioactive_decay_rxn_list
  do
    if (.not.associated(radioactive_decay_rxn)) exit
    prev_radioactive_decay_rxn => radioactive_decay_rxn
    radioactive_decay_rxn => radioactive_decay_rxn%next
    call GeneralRxnDestroy(prev_general_rxn)
  enddo    
  nullify(reaction%radioactive_decay_rxn_list)
  
  ! general reactions
  general_rxn => reaction%general_rxn_list
  do
    if (.not.associated(general_rxn)) exit
    prev_general_rxn => general_rxn
    general_rxn => general_rxn%next
    call GeneralRxnDestroy(prev_general_rxn)
  enddo    
  nullify(reaction%general_rxn_list)
  
  ! kd reactions
  kd_rxn => reaction%kd_rxn_list
  do
    if (.not.associated(kd_rxn)) exit
    prev_kd_rxn => kd_rxn
    kd_rxn => kd_rxn%next
    call KDRxnDestroy(prev_kd_rxn)
  enddo    
  nullify(reaction%kd_rxn_list)
  
  call SurfaceComplexationDestroy(reaction%surface_complexation)
  call MineralDestroy(reaction%mineral)
  call MicrobialDestroy(reaction%microbial)
  call ImmobileDestroy(reaction%immobile)
#ifdef SOLID_SOLUTION  
  call SolidSolutionDestroy(reaction%solid_solution_list)
#endif  

  if (associated(reaction%dbase_temperatures)) &
    deallocate(reaction%dbase_temperatures)
  nullify(reaction%dbase_temperatures)  
  
  ! redox species
  if (associated(reaction%redox_species_list)) &
    call AqueousSpeciesListDestroy(reaction%redox_species_list)
  nullify(reaction%redox_species_list)
  
  call DeallocateArray(reaction%primary_species_names)
  call DeallocateArray(reaction%secondary_species_names)
  call DeallocateArray(reaction%gas_species_names)
  call DeallocateArray(reaction%eqcplx_basis_names)
  call DeallocateArray(reaction%colloid_names)
  call DeallocateArray(reaction%colloid_species_names)  
  
  call DeallocateArray(reaction%primary_species_print)
  call DeallocateArray(reaction%secondary_species_print)
  call DeallocateArray(reaction%gas_species_print)
  call DeallocateArray(reaction%eqcplx_basis_print)
  call DeallocateArray(reaction%kd_print)
  call DeallocateArray(reaction%total_sorb_print)
  call DeallocateArray(reaction%total_sorb_mobile_print)
  call DeallocateArray(reaction%colloid_print)
  
  call DeallocateArray(reaction%primary_spec_a0)
  call DeallocateArray(reaction%primary_spec_Z)
  call DeallocateArray(reaction%primary_spec_molar_wt)
  
  call DeallocateArray(reaction%eqcplxspecid)
  call DeallocateArray(reaction%eqcplxstoich)
  call DeallocateArray(reaction%eqcplxh2oid)
  call DeallocateArray(reaction%eqcplxh2ostoich)
  call DeallocateArray(reaction%eqcplx_a0)
  call DeallocateArray(reaction%eqcplx_Z)
  call DeallocateArray(reaction%eqcplx_molar_wt)
  call DeallocateArray(reaction%eqcplx_logK)
  call DeallocateArray(reaction%eqcplx_logKcoef)
  
  call DeallocateArray(reaction%eqgasspecid)
  call DeallocateArray(reaction%eqgasstoich)
  call DeallocateArray(reaction%eqgash2oid)
  call DeallocateArray(reaction%eqgash2ostoich)
  call DeallocateArray(reaction%eqgas_logK)
  call DeallocateArray(reaction%eqgas_logKcoef)
  
  call DeallocateArray(reaction%eqionx_rxn_Z_flag)
  call DeallocateArray(reaction%eqionx_rxn_cation_X_offset)
  call DeallocateArray(reaction%eqionx_rxn_to_surf)
  call DeallocateArray(reaction%eqionx_rxn_CEC)
  call DeallocateArray(reaction%eqionx_rxn_k)
  call DeallocateArray(reaction%eqionx_rxn_cationid)
  
#if 0  
  call DeallocateArray(reaction%kinionx_CEC)
  call DeallocateArray(reaction%kinionx_k)
  call DeallocateArray(reaction%kinionx_cationid)
#endif
  
  call DeallocateArray(reaction%pri_spec_to_coll_spec)
  call DeallocateArray(reaction%coll_spec_to_pri_spec)
  call DeallocateArray(reaction%colloid_mobile_fraction)
  
  call DeallocateArray(reaction%radiodecayspecid)
  call DeallocateArray(reaction%radiodecaystoich)
  call DeallocateArray(reaction%radiodecayforwardspecid)
  call DeallocateArray(reaction%radiodecay_kf)
  
  call DeallocateArray(reaction%generalspecid)
  call DeallocateArray(reaction%generalstoich)
  call DeallocateArray(reaction%generalforwardspecid)
  call DeallocateArray(reaction%generalforwardstoich)
  call DeallocateArray(reaction%generalbackwardspecid)
  call DeallocateArray(reaction%generalbackwardstoich)
  call DeallocateArray(reaction%generalh2oid)
  call DeallocateArray(reaction%generalh2ostoich)
  call DeallocateArray(reaction%general_kf)
  call DeallocateArray(reaction%general_kr)
  
  call DeallocateArray(reaction%eqkdspecid)
  call DeallocateArray(reaction%eqkdtype)
  call DeallocateArray(reaction%eqkdspecid)
  call DeallocateArray(reaction%eqkddistcoef)
  call DeallocateArray(reaction%eqkdlangmuirb)
  call DeallocateArray(reaction%eqkdfreundlichn)
  
  deallocate(reaction)
  nullify(reaction)

end subroutine ReactionDestroy

end module Reaction_Aux_module
