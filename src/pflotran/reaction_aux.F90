module Reaction_Aux_module
  
  implicit none

  private 

#include "definitions.h"
  
  PetscInt, parameter, public :: SRFCMPLX_RXN_NULL = 0
  PetscInt, parameter, public :: SRFCMPLX_RXN_EQUILIBRIUM = 1
  PetscInt, parameter, public :: SRFCMPLX_RXN_MULTIRATE_KINETIC = 2
  PetscInt, parameter, public :: SRFCMPLX_RXN_KINETIC = 3
  
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

  type, public :: database_rxn_type
    PetscInt :: nspec
    character(len=MAXWORDLENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK(:)
  end type database_rxn_type

  type, public :: mineral_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(transition_state_rxn_type), pointer :: tstrxn
    type(mineral_type), pointer :: next
  end type mineral_type

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

  type, public :: transition_state_rxn_type
    PetscReal :: affinity_factor_sigma
    PetscReal :: affinity_factor_beta
    PetscReal :: affinity_threshold
    PetscReal :: rate_limiter
    PetscInt :: irreversible
    PetscReal :: rate
    PetscReal :: activation_energy
    type(transition_state_prefactor_type), pointer :: prefactor
    type(transition_state_rxn_type), pointer :: next
  end type transition_state_rxn_type
  
  type, public :: transition_state_prefactor_type
    type(ts_prefactor_species_type), pointer :: species
    ! these supercede the those above in transition_state_rxn_type
    PetscReal :: rate
    PetscReal :: activation_energy
    type(transition_state_prefactor_type), pointer :: next
  end type transition_state_prefactor_type

  type, public :: ts_prefactor_species_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    PetscReal :: alpha
    PetscReal :: beta
    PetscReal :: attenuation_coef
    type(ts_prefactor_species_type), pointer :: next
  end type ts_prefactor_species_type

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

  type, public :: surface_complex_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: free_site_stoich
    PetscReal :: Z
    PetscReal :: forward_rate
    PetscReal :: backward_rate
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(surface_complex_type), pointer :: next
  end type surface_complex_type

  type, public :: surface_complexation_rxn_type
    PetscInt :: id
    PetscInt :: itype
    PetscInt :: free_site_id
    character(len=MAXWORDLENGTH) :: free_site_name
    PetscBool :: free_site_print_me
    PetscInt :: mineral_id
    character(len=MAXWORDLENGTH) :: mineral_name
    character(len=MAXWORDLENGTH) :: colloid_name
    PetscReal :: site_density ! site density in mol/m^3 bulk
    type(surface_complex_type), pointer :: complex_list
    type (surface_complexation_rxn_type), pointer :: next
  end type surface_complexation_rxn_type    

  type, public :: kd_rxn_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: Kd
    PetscReal :: Langmuir_B
    PetscReal :: Freundlich_n
    type (kd_rxn_type), pointer :: next
  end type kd_rxn_type    

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

  type, public :: mineral_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_vol_frac(:)
    PetscReal, pointer :: constraint_area(:)
    PetscReal, pointer :: basis_vol_frac(:)
    PetscReal, pointer :: basis_area(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type mineral_constraint_type

  type, public :: srfcplx_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscReal, pointer :: basis_conc(:)
    PetscReal, pointer :: constraint_free_site_conc(:)
    PetscReal, pointer :: basis_free_site_conc(:)
  end type srfcplx_constraint_type

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
    PetscBool :: print_all_mineral_species
    PetscBool :: print_pH
    PetscBool :: print_kd
    PetscBool :: print_total_sorb
    PetscBool :: print_total_sorb_mobile
    PetscBool :: print_colloid
    PetscBool :: print_act_coefs
    PetscBool :: print_total_component
    PetscBool :: print_free_ion
    PetscBool :: initialize_with_molality
    PetscBool :: print_age
    PetscInt :: print_free_conc_type
    PetscInt :: print_tot_conc_type
    PetscInt :: print_secondary_conc_type
    PetscInt :: num_dbase_temperatures
    PetscReal, pointer :: dbase_temperatures(:)
    
    type(species_idx_type), pointer :: species_idx

    type(aq_species_type), pointer :: primary_species_list
    type(aq_species_type), pointer :: secondary_species_list
    type(gas_species_type), pointer :: gas_species_list
    type(mineral_type), pointer :: mineral_list
    type(colloid_type), pointer :: colloid_list
    type(ion_exchange_rxn_type), pointer :: ion_exchange_rxn_list
    type(surface_complexation_rxn_type), pointer :: surface_complexation_rxn_list
    type(general_rxn_type), pointer :: general_rxn_list
    type(kd_rxn_type), pointer :: kd_rxn_list
    type(aq_species_type), pointer :: redox_species_list
    PetscInt :: act_coef_update_frequency
    PetscInt :: act_coef_update_algorithm
    PetscBool :: checkpoint_activity_coefs
    PetscBool :: act_coef_use_bdot
    PetscBool :: use_activity_h2o
    PetscBool :: calculate_water_age
    PetscBool :: calculate_tracer_age
    
    ! compressed arrays for efficient computation
    ! primary aqueous complexes
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: ncollcomp
    
    ! offsets
    PetscInt :: offset_aq
    PetscInt :: offset_coll
    PetscInt :: offset_collcomp
    
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
    
    PetscInt :: neqsorb
    PetscBool, pointer :: kd_print(:)
    PetscBool, pointer :: total_sorb_print(:)

    ! ionx exchange reactions
    PetscInt :: neqionxrxn
    PetscInt :: neqionxcation 
    PetscBool, pointer :: eqionx_rxn_Z_flag(:)
    PetscInt, pointer :: eqionx_rxn_cation_X_offset(:)
    PetscReal, pointer :: eqionx_rxn_CEC(:)
    PetscReal, pointer :: eqionx_rxn_k(:,:)
    PetscInt, pointer :: eqionx_rxn_cationid(:,:)
#if 0    
    PetscReal, pointer :: kinionx_rxn_CEC(:)
    PetscReal, pointer :: kinionx_rxn_k(:,:)
    PetscInt, pointer :: kinionx_rxn_cationid(:)
#endif    

    ! surface complexation reactions
    PetscInt :: neqsrfcplx
    PetscInt :: neqsrfcplxrxn 
    PetscInt, pointer :: eqsrfcplx_rxn_to_surf(:)
    PetscInt, pointer :: eqsrfcplx_rxn_surf_type(:)
    PetscInt, pointer :: eqsrfcplx_rxn_to_complex(:,:)
    PetscReal, pointer :: eqsrfcplx_rxn_site_density(:) ! site density in mol/m^3 bulk
    PetscBool, pointer :: eqsrfcplx_rxn_stoich_flag(:)
    character(len=MAXWORDLENGTH), pointer :: eqsrfcplx_site_names(:)
    PetscBool, pointer :: eqsrfcplx_site_print(:)
    character(len=MAXWORDLENGTH), pointer :: eqsrfcplx_names(:)
    PetscBool, pointer :: eqsrfcplx_print(:)
    PetscInt, pointer :: eqsrfcplxspecid(:,:)
    PetscReal, pointer :: eqsrfcplxstoich(:,:)
    PetscInt, pointer :: eqsrfcplxh2oid(:)
    PetscReal, pointer :: eqsrfcplxh2ostoich(:)
    PetscInt, pointer :: eqsrfcplx_free_site_id(:)
    PetscReal, pointer :: eqsrfcplx_free_site_stoich(:)
!    PetscInt, pointer :: eqsrfcplx_mineral_id(:)
    PetscReal, pointer :: eqsrfcplx_logK(:)
    PetscReal, pointer :: eqsrfcplx_logKcoef(:,:)
    PetscReal, pointer :: eqsrfcplx_Z(:)  ! valence

    PetscInt :: nkinsrfcplx
    PetscInt :: nkinsrfcplxrxn
    PetscInt, pointer :: kinsrfcplx_rxn_to_surf(:)
    PetscInt, pointer :: kinsrfcplx_rxn_surf_type(:)
    PetscInt, pointer :: kinsrfcplx_rxn_to_complex(:,:) 
    PetscInt, pointer :: kinsrfcplx_rxn_to_site(:)
    PetscReal, pointer :: kinsrfcplx_rxn_site_density(:)
    PetscBool, pointer :: kinsrfcplx_rxn_stoich_flag(:)
    character(len=MAXWORDLENGTH), pointer :: kinsrfcplx_site_names(:)
    PetscBool, pointer :: kinsrfcplx_site_print(:)
    character(len=MAXWORDLENGTH), pointer :: kinsrfcplx_names(:)
    PetscBool, pointer :: kinsrfcplx_print(:)
    PetscInt, pointer :: kinsrfcplxspecid(:,:)
    PetscReal, pointer :: kinsrfcplxstoich(:,:)
    PetscInt, pointer :: kinsrfcplxh2oid(:)
    PetscReal, pointer :: kinsrfcplxh2ostoich(:)
    PetscInt, pointer :: kinsrfcplx_free_site_id(:)
    PetscReal, pointer :: kinsrfcplx_free_site_stoich(:)
    PetscReal, pointer :: kinsrfcplx_forward_rate(:)
    PetscReal, pointer :: kinsrfcplx_backward_rate(:)  
!    PetscReal, pointer :: kinsrfcplx_logK(:)
!    PetscReal, pointer :: kinsrfcplx_logKcoef(:,:)
    PetscReal, pointer :: kinsrfcplx_Z(:)  ! valence

    ! multirate kinetic surface complexation
    PetscInt :: kinmr_nrate 
    PetscReal, pointer :: kinmr_rate(:)
    PetscReal, pointer :: kinmr_frac(:)
    PetscReal :: kinmr_scale_factor
    
    ! mineral reactions
    PetscInt :: nmnrl
    character(len=MAXWORDLENGTH), pointer :: mineral_names(:)
    
    ! colloids
    PetscInt :: ncoll
    character(len=MAXWORDLENGTH), pointer :: colloid_names(:)
    character(len=MAXWORDLENGTH), pointer :: colloid_species_names(:)
    PetscReal, pointer :: colloid_mobile_fraction(:)
    PetscInt, pointer :: pri_spec_to_coll_spec(:)
    PetscInt, pointer :: coll_spec_to_pri_spec(:)
    PetscBool, pointer :: total_sorb_mobile_print(:)
    PetscBool, pointer :: colloid_print(:)
    
      ! for saturation states
    PetscInt, pointer :: mnrlspecid(:,:)
    PetscReal, pointer :: mnrlstoich(:,:)
    PetscInt, pointer :: mnrlh2oid(:)
    PetscReal, pointer :: mnrlh2ostoich(:)
    PetscReal, pointer :: mnrl_logK(:)
    PetscReal, pointer :: mnrl_logKcoef(:,:)
    PetscBool, pointer :: mnrl_print(:)
    
      ! for kinetic reactions
    PetscInt :: nkinmnrl
    character(len=MAXWORDLENGTH), pointer :: kinmnrl_names(:)
    PetscBool, pointer :: kinmnrl_print(:)
    PetscInt, pointer :: kinmnrlspecid(:,:)
    PetscReal, pointer :: kinmnrlstoich(:,:)
    PetscInt, pointer :: kinmnrlh2oid(:)
    PetscReal, pointer :: kinmnrlh2ostoich(:)
    PetscReal, pointer :: kinmnrl_logK(:)
    PetscReal, pointer :: kinmnrl_logKcoef(:,:)
    PetscReal, pointer :: kinmnrl_rate(:)
    PetscReal, pointer :: kinmnrl_activation_energy(:)
    PetscReal, pointer :: kinmnrl_molar_vol(:)
    PetscReal, pointer :: kinmnrl_molar_wt(:)
    PetscInt, pointer :: kinmnrl_num_prefactors(:)
    PetscInt, pointer :: kinmnrl_prefactor_id(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_alpha(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_beta(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_atten_coef(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_rate(:,:)
    PetscReal, pointer :: kinmnrl_pref_activation_energy(:,:)
    PetscReal, pointer :: kinmnrl_Tempkin_const(:)
    PetscReal, pointer :: kinmnrl_affinity_power(:)
    PetscReal, pointer :: kinmnrl_affinity_threshold(:)
    PetscReal, pointer :: kinmnrl_rate_limiter(:)
    PetscInt, pointer :: kinmnrl_irreversible(:)
    
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

  end type reaction_type

  public :: ReactionCreate, &
            SpeciesIndexCreate, &
            GasSpeciesCreate, &
            GetPrimarySpeciesCount, &
            GetPrimarySpeciesNames, &
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
            DatabaseRxnCreate, &
            DatabaseRxnDestroy, &
            TransitionStateTheoryRxnCreate, &
            TransitionStatePrefactorCreate, &
            TSPrefactorSpeciesCreate, &
            TransitionStateTheoryRxnDestroy, &
            SurfaceComplexationRxnCreate, &
            AqueousSpeciesCreate, &
            AqueousSpeciesDestroy, &
            AqueousSpeciesConstraintCreate, &
            AqueousSpeciesConstraintDestroy, &
            MineralCreate, &
            MineralDestroy, &
            MineralConstraintCreate, &
            MineralConstraintDestroy, &
            SurfaceComplexCreate, &
            SurfaceComplexDestroy, &
            SurfaceComplexConstraintCreate, &
            SurfaceComplexConstraintDestroy, &
            GeneralRxnCreate, &
            GeneralRxnDestroy, &
            KDRxnCreate, &
            KDRxnDestroy, &
            ColloidCreate, &
            ColloidDestroy, &
            ColloidConstraintCreate, &
            ColloidConstraintDestroy, &
            IonExchangeRxnCreate, &
            IonExchangeCationCreate
             
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
  reaction%print_all_mineral_species = PETSC_FALSE
  reaction%print_pH = PETSC_FALSE
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

  reaction%initialize_with_molality = PETSC_FALSE
  reaction%print_free_conc_type = 0
  reaction%print_tot_conc_type = 0
  reaction%print_secondary_conc_type = 0
  
  nullify(reaction%species_idx)

  nullify(reaction%primary_species_list)
  nullify(reaction%secondary_species_list)
  nullify(reaction%gas_species_list)
  nullify(reaction%mineral_list)
  nullify(reaction%colloid_list)
  nullify(reaction%ion_exchange_rxn_list)
  nullify(reaction%surface_complexation_rxn_list)
  nullify(reaction%general_rxn_list)
  nullify(reaction%kd_rxn_list)
  nullify(reaction%redox_species_list)
  
  nullify(reaction%primary_species_names)
  nullify(reaction%secondary_species_names)
  nullify(reaction%eqcplx_basis_names)
  nullify(reaction%gas_species_names)
  nullify(reaction%eqsrfcplx_site_names)
  nullify(reaction%eqsrfcplx_names)
  nullify(reaction%kinsrfcplx_site_names)
  nullify(reaction%kinsrfcplx_names)
  nullify(reaction%mineral_names)
  nullify(reaction%colloid_names)
  nullify(reaction%colloid_species_names)
  nullify(reaction%kinmnrl_names)

  nullify(reaction%primary_species_print)
  nullify(reaction%secondary_species_print)
  nullify(reaction%eqcplx_basis_print)
  nullify(reaction%gas_species_print)
  nullify(reaction%eqsrfcplx_site_print)
  nullify(reaction%eqsrfcplx_print)
  nullify(reaction%mnrl_print)
  nullify(reaction%kinsrfcplx_site_print)
  nullify(reaction%kinsrfcplx_print)
  nullify(reaction%kinmnrl_print)
  nullify(reaction%kd_print)
  nullify(reaction%total_sorb_print)
  nullify(reaction%total_sorb_mobile_print)
  nullify(reaction%colloid_print)
  
  reaction%ncomp = 0
  reaction%naqcomp = 0
  reaction%ncoll = 0
  reaction%ncollcomp = 0
  reaction%offset_aq = 0
  reaction%offset_coll = 0
  reaction%offset_collcomp = 0
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

  reaction%neqionxrxn = 0
  reaction%neqionxcation = 0
  nullify(reaction%eqionx_rxn_Z_flag)
  nullify(reaction%eqionx_rxn_cation_X_offset)
  nullify(reaction%eqionx_rxn_CEC)
  nullify(reaction%eqionx_rxn_k)
  nullify(reaction%eqionx_rxn_cationid)
#if 0  
  nullify(reaction%kinionx_CEC)
  nullify(reaction%kinionx_k)
  nullify(reaction%kinionx_cationid)
#endif

  reaction%neqsorb = 0
  reaction%neqsrfcplx = 0
  reaction%neqsrfcplxrxn = 0
  nullify(reaction%eqsrfcplx_rxn_to_surf)
  nullify(reaction%eqsrfcplx_rxn_surf_type)
  nullify(reaction%eqsrfcplx_rxn_to_complex)
  nullify(reaction%eqsrfcplx_rxn_site_density)
  nullify(reaction%eqsrfcplx_rxn_stoich_flag) 
  nullify(reaction%eqsrfcplx_site_names)
  nullify(reaction%eqsrfcplx_names)
  nullify(reaction%eqsrfcplxspecid)
  nullify(reaction%eqsrfcplxstoich)
  nullify(reaction%eqsrfcplxh2oid)
  nullify(reaction%eqsrfcplxh2ostoich)
!  nullify(reaction%eqsrfcplx_mineral_id)
  nullify(reaction%eqsrfcplx_free_site_id)
  nullify(reaction%eqsrfcplx_free_site_stoich)
  nullify(reaction%eqsrfcplx_logK)
  nullify(reaction%eqsrfcplx_logKcoef)
  nullify(reaction%eqsrfcplx_Z)

  reaction%nkinsrfcplx = 0
  reaction%nkinsrfcplxrxn = 0
  nullify(reaction%kinsrfcplx_rxn_to_surf)
  nullify(reaction%kinsrfcplx_rxn_to_complex)
  nullify(reaction%kinsrfcplx_rxn_to_site)
  nullify(reaction%kinsrfcplx_rxn_site_density)
  nullify(reaction%kinsrfcplx_rxn_stoich_flag) 
  nullify(reaction%kinsrfcplx_site_names)
  nullify(reaction%kinsrfcplx_names)
  nullify(reaction%kinsrfcplxspecid)
  nullify(reaction%kinsrfcplxstoich)
  nullify(reaction%kinsrfcplxh2oid)
  nullify(reaction%kinsrfcplxh2ostoich)
  nullify(reaction%kinsrfcplx_free_site_id)
  nullify(reaction%kinsrfcplx_free_site_stoich)
!  nullify(reaction%kinsrfcplx_logK)
!  nullify(reaction%kinsrfcplx_logKcoef)
  nullify(reaction%kinsrfcplx_Z)

  reaction%kinmr_nrate = 0
  nullify(reaction%kinmr_rate)
  nullify(reaction%kinmr_frac)
  reaction%kinmr_scale_factor = 1.d0

  reaction%nmnrl = 0  
  nullify(reaction%mnrlspecid)
  nullify(reaction%mnrlstoich)
  nullify(reaction%mnrlh2oid)
  nullify(reaction%mnrlstoich)
  nullify(reaction%mnrlh2ostoich)
  nullify(reaction%mnrl_logKcoef)

  reaction%ncoll = 0
  nullify(reaction%pri_spec_to_coll_spec)
  nullify(reaction%coll_spec_to_pri_spec)
  nullify(reaction%colloid_mobile_fraction)
  
  reaction%nkinmnrl = 0  
  nullify(reaction%kinmnrlspecid)
  nullify(reaction%kinmnrlstoich)
  nullify(reaction%kinmnrlh2oid)
  nullify(reaction%kinmnrlh2ostoich)
  nullify(reaction%kinmnrl_logK)
  nullify(reaction%kinmnrl_logKcoef)
  nullify(reaction%kinmnrl_rate)
  nullify(reaction%kinmnrl_activation_energy)
  nullify(reaction%kinmnrl_molar_vol)
  nullify(reaction%kinmnrl_molar_wt)

  nullify(reaction%kinmnrl_num_prefactors)
  nullify(reaction%kinmnrl_prefactor_id)
  nullify(reaction%kinmnrl_pref_alpha)
  nullify(reaction%kinmnrl_pref_beta)
  nullify(reaction%kinmnrl_pref_atten_coef)
  nullify(reaction%kinmnrl_pref_rate)
  nullify(reaction%kinmnrl_pref_activation_energy)

  nullify(reaction%kinmnrl_Tempkin_const)
  nullify(reaction%kinmnrl_affinity_power)
  nullify(reaction%kinmnrl_affinity_threshold)
  nullify(reaction%kinmnrl_irreversible)
  nullify(reaction%kinmnrl_rate_limiter)
  
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

  reaction%neqkdrxn = 0
  nullify(reaction%eqkdspecid)
  nullify(reaction%eqkdtype)
  nullify(reaction%eqkddistcoef)
  nullify(reaction%eqkdlangmuirb)
  nullify(reaction%eqkdfreundlichn)
      
  reaction%max_dlnC = 5.d0
  reaction%max_relative_change_tolerance = 1.d-6
  reaction%max_residual_tolerance = 1.d-12

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
! MineralCreate: Allocate and initialize a mineral object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function MineralCreate()

  use Option_module

  implicit none
  
  type(mineral_type), pointer :: MineralCreate
  
  type(mineral_type), pointer :: mineral

  allocate(mineral)  
  mineral%id = 0
  mineral%itype = 0
  mineral%name = ''
  mineral%molar_volume = 0.d0
  mineral%molar_weight = 0.d0
  mineral%print_me = PETSC_FALSE
  nullify(mineral%tstrxn)
  nullify(mineral%next)
  
  MineralCreate => mineral
  
end function MineralCreate

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
! DatabaseRxnCreate: Allocate and initialize an equilibrium reaction
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
function DatabaseRxnCreate()

  implicit none
    
  type(database_rxn_type), pointer :: DatabaseRxnCreate

  type(database_rxn_type), pointer :: dbaserxn

  allocate(dbaserxn)
  dbaserxn%nspec = 0
  nullify(dbaserxn%spec_name)
  nullify(dbaserxn%stoich)
  nullify(dbaserxn%spec_ids)
  nullify(dbaserxn%logK)
  
  DatabaseRxnCreate => dbaserxn
  
end function DatabaseRxnCreate

! ************************************************************************** !
!
! TransitionStateTheoryRxnCreate: Allocate and initialize a transition state
!                                 theory reaction
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
function TransitionStateTheoryRxnCreate()

  implicit none
    
  type(transition_state_rxn_type), pointer :: TransitionStateTheoryRxnCreate

  type(transition_state_rxn_type), pointer :: tstrxn

  allocate(tstrxn)
  tstrxn%affinity_factor_sigma = 0.d0
  tstrxn%affinity_factor_beta = 0.d0
  tstrxn%affinity_threshold = 0.d0
  tstrxn%rate_limiter = 0.d0
  tstrxn%irreversible = 0
  tstrxn%activation_energy = 0.d0
  tstrxn%rate = 0.d0
  nullify(tstrxn%prefactor)
  nullify(tstrxn%next)
  
  TransitionStateTheoryRxnCreate => tstrxn
  
end function TransitionStateTheoryRxnCreate

! ************************************************************************** !
!
! TransitionStatePrefactorCreate: Allocate and initialize a transition state
!                                 theory prefactor
! author: Glenn Hammond
! date: 07/29/11
!
! ************************************************************************** !
function TransitionStatePrefactorCreate()

  implicit none
    
  type(transition_state_prefactor_type), pointer :: TransitionStatePrefactorCreate

  type(transition_state_prefactor_type), pointer :: prefactor

  allocate(prefactor)
  prefactor%rate = 0.d0
  prefactor%activation_energy = 0.d0
  nullify(prefactor%species)
  nullify(prefactor%next)
  
  TransitionStatePrefactorCreate => prefactor
  
end function TransitionStatePrefactorCreate

! ************************************************************************** !
!
! TSPrefactorSpeciesCreate: Allocate and initialize a transition state
!                           theory prefactor species
! author: Glenn Hammond
! date: 08/01/11
!
! ************************************************************************** !
function TSPrefactorSpeciesCreate()

  implicit none
    
  type(ts_prefactor_species_type), pointer :: TSPrefactorSpeciesCreate

  type(ts_prefactor_species_type), pointer :: species

  allocate(species)
  species%name = ''
  species%id = 0
  species%alpha = 0.d0
  species%beta = 0.d0
  species%attenuation_coef = 0.d0
  nullify(species%next)
  
  TSPrefactorSpeciesCreate => species
  
end function TSPrefactorSpeciesCreate

! ************************************************************************** !
!
! SurfaceComplexationRxnCreate: Allocate and initialize a surface complexation
!                               reaction
! author: Peter Lichtner
! date: 10/21/08
!
! ************************************************************************** !
function SurfaceComplexationRxnCreate()

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: SurfaceComplexationRxnCreate

  type(surface_complexation_rxn_type), pointer :: srfcplxrxn
  
  allocate(srfcplxrxn)
  srfcplxrxn%free_site_id = 0
  srfcplxrxn%itype = SRFCMPLX_RXN_NULL
  srfcplxrxn%free_site_name = ''
  srfcplxrxn%free_site_print_me = PETSC_FALSE

  srfcplxrxn%mineral_id = 0
  srfcplxrxn%mineral_name = ''
  srfcplxrxn%colloid_name = ''
  srfcplxrxn%site_density = 0.d0
  
  nullify(srfcplxrxn%complex_list)
  nullify(srfcplxrxn%next)
  
  SurfaceComplexationRxnCreate => srfcplxrxn
  
end function SurfaceComplexationRxnCreate

! ************************************************************************** !
!
! SurfaceComplexCreate: Allocate and initialize a surface complex reaction
! author: Peter Lichtner
! date: 10/21/08
!
! ************************************************************************** !
function SurfaceComplexCreate()

  implicit none
    
  type(surface_complex_type), pointer :: SurfaceComplexCreate

  type(surface_complex_type), pointer :: srfcplx
  
  allocate(srfcplx)
  srfcplx%id = 0
  srfcplx%name = ''
  srfcplx%Z = 0.d0
  srfcplx%free_site_stoich = 0.d0
  srfcplx%forward_rate = 0.d0
  srfcplx%backward_rate = -999.d0
  srfcplx%print_me = PETSC_FALSE
  nullify(srfcplx%dbaserxn)
  nullify(srfcplx%next)
  
  SurfaceComplexCreate => srfcplx
  
end function SurfaceComplexCreate

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
! MineralConstraintCreate: Creates a mineral constraint object
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function MineralConstraintCreate(reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  type(mineral_constraint_type), pointer :: MineralConstraintCreate

  type(mineral_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(reaction%nkinmnrl))
  constraint%names = ''
  allocate(constraint%constraint_vol_frac(reaction%nkinmnrl))
  constraint%constraint_vol_frac = 0.d0
  allocate(constraint%constraint_area(reaction%nkinmnrl))
  constraint%constraint_area = 0.d0
  allocate(constraint%basis_vol_frac(reaction%nkinmnrl))
  constraint%basis_vol_frac = 0.d0
  allocate(constraint%basis_area(reaction%nkinmnrl))
  constraint%basis_area = 0.d0
  allocate(constraint%constraint_aux_string(reaction%nkinmnrl))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(reaction%nkinmnrl))
  constraint%external_dataset = PETSC_FALSE

  MineralConstraintCreate => constraint

end function MineralConstraintCreate
  
! ************************************************************************** !
!
! SurfaceComplexConstraintCreate: Creates a surface complex constraint object
! author: Glenn Hammond
! date: 12/21/09
!
! ************************************************************************** !
function SurfaceComplexConstraintCreate(reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  type(srfcplx_constraint_type), pointer :: SurfaceComplexConstraintCreate

  type(srfcplx_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(reaction%nkinsrfcplx))
  constraint%names = ''
  allocate(constraint%constraint_conc(reaction%nkinsrfcplx))
  constraint%constraint_conc = 0.d0
  allocate(constraint%basis_conc(reaction%nkinsrfcplx))
  constraint%basis_conc = 0.d0
  allocate(constraint%constraint_free_site_conc(reaction%nkinsrfcplxrxn))
  constraint%constraint_free_site_conc = 0.d0
  allocate(constraint%basis_free_site_conc(reaction%nkinsrfcplxrxn))
  constraint%basis_free_site_conc = 0.d0

  SurfaceComplexConstraintCreate => constraint

end function SurfaceComplexConstraintCreate
  
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
! GetMineralNames: Returns the names of minerals in an array
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
function GetMineralNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetMineralNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(mineral_type), pointer :: mineral

  count = GetMineralCount(reaction)
  allocate(names(count))
  
  count = 1
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    names(count) = mineral%name
    count = count + 1
    mineral => mineral%next
  enddo

  GetMineralNames => names
  
end function GetMineralNames

! ************************************************************************** !
!
! GetMineralCount: Returns the number of minerals
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetMineralCount(reaction)

  implicit none
  
  PetscInt :: GetMineralCount
  type(reaction_type) :: reaction

  type(mineral_type), pointer :: mineral

  GetMineralCount = 0
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    GetMineralCount = GetMineralCount + 1
    mineral => mineral%next
  enddo

end function GetMineralCount

! ************************************************************************** !
!
! GetMineralIDFromName: Returns the id of mineral with the corresponding name
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
function GetMineralIDFromName(reaction,name)

  use String_module
  
  implicit none
  
  type(reaction_type) :: reaction
  character(len=MAXWORDLENGTH) :: name

  PetscInt :: GetMineralIDFromName
  type(mineral_type), pointer :: mineral

  GetMineralIDFromName = -1
 
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    if (StringCompare(name,mineral%name,MAXWORDLENGTH)) then
      GetMineralIDFromName = mineral%id
      exit
    endif
    mineral => mineral%next
  enddo

end function GetMineralIDFromName

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
! SpeciesIndexDestroy: Deallocates a species index object
! author: Glenn Hammond
! date: 01/29/10
!
! ************************************************************************** !
subroutine SpeciesIndexDestroy(species_idx)

  implicit none
    
  type(species_idx_type), pointer :: species_idx

  deallocate(species_idx)  
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
! MineralDestroy: Deallocates a mineral
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine MineralDestroy(mineral)

  implicit none
    
  type(mineral_type), pointer :: mineral

  if (associated(mineral%dbaserxn)) &
    call DatabaseRxnDestroy(mineral%dbaserxn)
  if (associated(mineral%tstrxn)) &
    call TransitionStateTheoryRxnDestroy(mineral%tstrxn)

  deallocate(mineral)  
  nullify(mineral)

end subroutine MineralDestroy

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
! DatabaseRxnDestroy: Deallocates a database reaction
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine DatabaseRxnDestroy(dbaserxn)

  implicit none
    
  type(database_rxn_type), pointer :: dbaserxn

  if (.not.associated(dbaserxn)) return
  
  if (associated(dbaserxn%spec_name)) deallocate(dbaserxn%spec_name)
  nullify(dbaserxn%spec_name)
  if (associated(dbaserxn%spec_ids)) deallocate(dbaserxn%spec_ids)
  nullify(dbaserxn%spec_ids)
  if (associated(dbaserxn%stoich)) deallocate(dbaserxn%stoich)
  nullify(dbaserxn%stoich)
  if (associated(dbaserxn%logK)) deallocate(dbaserxn%logK)
  nullify(dbaserxn%logK)

  deallocate(dbaserxn)  
  nullify(dbaserxn)

end subroutine DatabaseRxnDestroy

! ************************************************************************** !
!
! TransitionStateTheoryRxnDestroy: Deallocates a transition state reaction
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
recursive subroutine TransitionStateTheoryRxnDestroy(tstrxn)

  implicit none
    
  type(transition_state_rxn_type), pointer :: tstrxn

  if (.not.associated(tstrxn)) return
  
  call TransitionStateTheoryRxnDestroy(tstrxn%next)
  call TransitionStatePrefactorDestroy(tstrxn%prefactor)

  deallocate(tstrxn)  
  nullify(tstrxn)

end subroutine TransitionStateTheoryRxnDestroy

! ************************************************************************** !
!
! TransitionStatePrefactorDestroy: Deallocates a transition state prefactor
! author: Glenn Hammond
! date: 07/29/11
!
! ************************************************************************** !
recursive subroutine TransitionStatePrefactorDestroy(prefactor)

  implicit none
    
  type(transition_state_prefactor_type), pointer :: prefactor

  if (.not.associated(prefactor)) return
  
  call TransitionStatePrefactorDestroy(prefactor%next)
  call TSPrefactorSpeciesDestroy(prefactor%species)

  deallocate(prefactor)  
  nullify(prefactor)

end subroutine TransitionStatePrefactorDestroy

! ************************************************************************** !
!
! TSPrefactorSpeciesDestroy: Deallocates a transition state prefactor
! author: Glenn Hammond
! date: 08/01/11
!
! ************************************************************************** !
recursive subroutine TSPrefactorSpeciesDestroy(species)

  implicit none
    
  type(ts_prefactor_species_type), pointer :: species

  if (.not.associated(species)) return
  
  call TSPrefactorSpeciesDestroy(species%next)

  deallocate(species)  
  nullify(species)

end subroutine TSPrefactorSpeciesDestroy

! ************************************************************************** !
!
! SurfaceComplexationRxnDestroy: Deallocates a surface complexation reaction
! author: Glenn Hammond
! date: 10/21/08
!
! ************************************************************************** !
subroutine SurfaceComplexationRxnDestroy(srfcplxrxn)

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: srfcplxrxn

  type(surface_complex_type), pointer :: cur_srfcplx, prev_srfcplx
  
  if (.not.associated(srfcplxrxn)) return
  
  cur_srfcplx => srfcplxrxn%complex_list
  do
    if (.not.associated(cur_srfcplx)) exit
    prev_srfcplx => cur_srfcplx
    cur_srfcplx => cur_srfcplx%next
    call SurfaceComplexDestroy(prev_srfcplx)
    nullify(prev_srfcplx)
  enddo
  
  deallocate(srfcplxrxn)  
  nullify(srfcplxrxn)

end subroutine SurfaceComplexationRxnDestroy

! ************************************************************************** !
!
! SurfaceComplexDestroy: Deallocates a surface complex
! author: Glenn Hammond
! date: 10/21/08
!
! ************************************************************************** !
subroutine SurfaceComplexDestroy(srfcplx)

  implicit none
    
  type(surface_complex_type), pointer :: srfcplx

  if (.not.associated(srfcplx)) return
  
  if (associated(srfcplx%dbaserxn)) &
    call DatabaseRxnDestroy(srfcplx%dbaserxn)
  nullify(srfcplx%dbaserxn)
  nullify(srfcplx%next)

  deallocate(srfcplx)  
  nullify(srfcplx)

end subroutine SurfaceComplexDestroy

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

  implicit none
  
  type(aq_species_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  if (associated(constraint%names)) &
    deallocate(constraint%names)
  nullify(constraint%names)
  if (associated(constraint%constraint_conc)) &
    deallocate(constraint%constraint_conc)
  nullify(constraint%constraint_conc)
  if (associated(constraint%basis_molarity)) &
    deallocate(constraint%basis_molarity)
  nullify(constraint%basis_molarity)
  if (associated(constraint%constraint_type)) &
    deallocate(constraint%constraint_type)
  nullify(constraint%constraint_type)
  if (associated(constraint%constraint_spec_id)) &
    deallocate(constraint%constraint_spec_id)
  nullify(constraint%constraint_spec_id)
  if (associated(constraint%constraint_aux_string)) &
    deallocate(constraint%constraint_aux_string)
  nullify(constraint%constraint_aux_string)
  if (associated(constraint%constraint_aux_string)) &
    deallocate(constraint%constraint_aux_string)
  nullify(constraint%constraint_aux_string)

  deallocate(constraint)
  nullify(constraint)

end subroutine AqueousSpeciesConstraintDestroy

! ************************************************************************** !
!
! MineralConstraintDestroy: Destroys a mineral constraint object
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine MineralConstraintDestroy(constraint)

  implicit none
  
  type(mineral_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  if (associated(constraint%names)) &
    deallocate(constraint%names)
  nullify(constraint%names)
  if (associated(constraint%constraint_vol_frac)) &
    deallocate(constraint%constraint_vol_frac)
  nullify(constraint%constraint_vol_frac)
  if (associated(constraint%constraint_area)) &
    deallocate(constraint%constraint_area)
  nullify(constraint%constraint_area)
  if (associated(constraint%basis_vol_frac)) &
    deallocate(constraint%basis_vol_frac)
  nullify(constraint%basis_vol_frac)
  if (associated(constraint%basis_area)) &
    deallocate(constraint%basis_area)
  nullify(constraint%basis_area)
  if (associated(constraint%constraint_aux_string)) &
    deallocate(constraint%constraint_aux_string)
  nullify(constraint%constraint_aux_string)
  if (associated(constraint%constraint_aux_string)) &
    deallocate(constraint%constraint_aux_string)
  nullify(constraint%constraint_aux_string)

  deallocate(constraint)
  nullify(constraint)

end subroutine MineralConstraintDestroy

! ************************************************************************** !
!
! SurfaceComplexConstraintDestroy: Destroys a surface complex constraint 
!                                  object
! author: Glenn Hammond
! date: 12/21/09
!
! ************************************************************************** !
subroutine SurfaceComplexConstraintDestroy(constraint)

  implicit none
  
  type(srfcplx_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  if (associated(constraint%names)) &
    deallocate(constraint%names)
  nullify(constraint%names)
  if (associated(constraint%constraint_conc)) &
    deallocate(constraint%constraint_conc)
  nullify(constraint%constraint_conc)
  if (associated(constraint%basis_conc)) &
    deallocate(constraint%basis_conc)
  nullify(constraint%basis_conc)
  if (associated(constraint%constraint_free_site_conc)) &
    deallocate(constraint%constraint_free_site_conc)
  nullify(constraint%constraint_free_site_conc)
  if (associated(constraint%basis_free_site_conc)) &
    deallocate(constraint%basis_free_site_conc)
  nullify(constraint%basis_free_site_conc)

  deallocate(constraint)
  nullify(constraint)

end subroutine SurfaceComplexConstraintDestroy

! ************************************************************************** !
!
! ColloidConstraintDestroy: Destroys a colloid constraint object
! author: Glenn Hammond
! date: 03/12/10
!
! ************************************************************************** !
subroutine ColloidConstraintDestroy(constraint)

  implicit none
  
  type(colloid_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  if (associated(constraint%names)) &
    deallocate(constraint%names)
  nullify(constraint%names)
  if (associated(constraint%constraint_conc_mob)) &
    deallocate(constraint%constraint_conc_mob)
  nullify(constraint%constraint_conc_mob)
  if (associated(constraint%constraint_conc_imb)) &
    deallocate(constraint%constraint_conc_imb)
  nullify(constraint%constraint_conc_imb)
  if (associated(constraint%basis_conc_mob)) &
    deallocate(constraint%basis_conc_mob)
  nullify(constraint%basis_conc_mob)
  if (associated(constraint%basis_conc_imb)) &
    deallocate(constraint%basis_conc_imb)
  nullify(constraint%basis_conc_imb)

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

  implicit none

  type(reaction_type), pointer :: reaction
  
  type(aq_species_type), pointer :: aq_species, prev_aq_species
  type(gas_species_type), pointer :: gas_species, prev_gas_species
  type(mineral_type), pointer :: mineral, prev_mineral
  type(colloid_type), pointer :: colloid, prev_colloid
  type(ion_exchange_rxn_type), pointer :: ionxrxn, prev_ionxrxn
  type(surface_complexation_rxn_type), pointer :: srfcplxrxn, prev_srfcplxrxn
  type(general_rxn_type), pointer :: general_rxn, prev_general_rxn
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
  
  ! mineral species
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    prev_mineral => mineral
    mineral => mineral%next
    call MineralDestroy(prev_mineral)
  enddo    
  nullify(reaction%mineral_list)
  
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

  ! surface complexation reactions
  srfcplxrxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(srfcplxrxn)) exit
    prev_srfcplxrxn => srfcplxrxn
    srfcplxrxn => srfcplxrxn%next
    call SurfaceComplexationRxnDestroy(prev_srfcplxrxn)
  enddo    
  nullify(reaction%surface_complexation_rxn_list)
  
  ! general reactions
  general_rxn => reaction%general_rxn_list
  do
    if (.not.associated(general_rxn)) exit
    prev_general_rxn => general_rxn
    general_rxn => general_rxn%next
    call GeneralRxnDestroy(prev_general_rxn)
  enddo    
  nullify(reaction%surface_complexation_rxn_list)
  
  ! kd reactions
  kd_rxn => reaction%kd_rxn_list
  do
    if (.not.associated(kd_rxn)) exit
    prev_kd_rxn => kd_rxn
    kd_rxn => kd_rxn%next
    call KDRxnDestroy(prev_kd_rxn)
  enddo    
  nullify(reaction%kd_rxn_list)

  ! redox species
  if (associated(reaction%redox_species_list)) &
    call AqueousSpeciesListDestroy(reaction%redox_species_list)
  nullify(reaction%redox_species_list)

  if (associated(reaction%primary_species_names)) &
    deallocate(reaction%primary_species_names)
  nullify(reaction%primary_species_names)
  if (associated(reaction%secondary_species_names)) &
    deallocate(reaction%secondary_species_names)
  nullify(reaction%secondary_species_names)
  if (associated(reaction%gas_species_names)) &
    deallocate(reaction%gas_species_names)
  nullify(reaction%gas_species_names)
  if (associated(reaction%eqsrfcplx_site_names)) &
    deallocate(reaction%eqsrfcplx_site_names)
  nullify(reaction%eqsrfcplx_site_names)
  if (associated(reaction%eqsrfcplx_names)) &
    deallocate(reaction%eqsrfcplx_names)
  nullify(reaction%eqsrfcplx_names)
  if (associated(reaction%kinsrfcplx_site_names)) &
    deallocate(reaction%kinsrfcplx_site_names)
  nullify(reaction%kinsrfcplx_site_names)
  if (associated(reaction%kinsrfcplx_names)) &
    deallocate(reaction%kinsrfcplx_names)
  nullify(reaction%kinsrfcplx_names)
  if (associated(reaction%mineral_names)) &
    deallocate(reaction%mineral_names)
  nullify(reaction%mineral_names)
  if (associated(reaction%colloid_names)) &
    deallocate(reaction%colloid_names)
  nullify(reaction%colloid_names)
  if (associated(reaction%colloid_species_names)) &
    deallocate(reaction%colloid_species_names)
  nullify(reaction%colloid_species_names)
  if (associated(reaction%kinmnrl_names)) &
    deallocate(reaction%kinmnrl_names)
  nullify(reaction%kinmnrl_names)

  if (associated(reaction%primary_species_print)) &
    deallocate(reaction%primary_species_print)
  nullify(reaction%primary_species_print)
  if (associated(reaction%secondary_species_print)) &
    deallocate(reaction%secondary_species_print)
  nullify(reaction%secondary_species_print)
  if (associated(reaction%eqcplx_basis_print)) &
    deallocate(reaction%eqcplx_basis_print)
  nullify(reaction%eqcplx_basis_print)
  if (associated(reaction%gas_species_print)) &
    deallocate(reaction%gas_species_print)
  nullify(reaction%gas_species_print)
  if (associated(reaction%eqsrfcplx_site_print)) &
    deallocate(reaction%eqsrfcplx_site_print)
  nullify(reaction%eqsrfcplx_site_print)
  if (associated(reaction%eqsrfcplx_print)) &
    deallocate(reaction%eqsrfcplx_print)
  nullify(reaction%eqsrfcplx_print)
  if (associated(reaction%mnrl_print)) &
    deallocate(reaction%mnrl_print)
  nullify(reaction%mnrl_print)
  if (associated(reaction%kinsrfcplx_site_print)) &
    deallocate(reaction%kinsrfcplx_site_print)
  nullify(reaction%kinsrfcplx_site_print)
  if (associated(reaction%kinsrfcplx_print)) &
    deallocate(reaction%kinsrfcplx_print)
  nullify(reaction%kinsrfcplx_print)
  if (associated(reaction%kinmnrl_print)) &
    deallocate(reaction%kinmnrl_print)
  nullify(reaction%kinmnrl_print)
  if (associated(reaction%kd_print)) &
    deallocate(reaction%kd_print)
  nullify(reaction%kd_print)
  if (associated(reaction%total_sorb_print)) &
    deallocate(reaction%total_sorb_print)
  nullify(reaction%total_sorb_print)
  if (associated(reaction%total_sorb_mobile_print)) &
    deallocate(reaction%total_sorb_mobile_print)
  nullify(reaction%total_sorb_mobile_print)
  if (associated(reaction%colloid_print)) &
    deallocate(reaction%colloid_print)
  nullify(reaction%colloid_print)
    
    
  if (associated(reaction%primary_spec_a0)) &
    deallocate(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_a0)
  if (associated(reaction%primary_spec_Z)) &
    deallocate(reaction%primary_spec_Z)
  nullify(reaction%primary_spec_Z)
  if (associated(reaction%primary_spec_molar_wt)) &
    deallocate(reaction%primary_spec_molar_wt)
  nullify(reaction%primary_spec_molar_wt)
  
  if (associated(reaction%eqcplxspecid)) &
    deallocate(reaction%eqcplxspecid)
  nullify(reaction%eqcplxspecid)
  if (associated(reaction%eqcplxstoich)) &
    deallocate(reaction%eqcplxstoich)
  nullify(reaction%eqcplxstoich)
  if (associated(reaction%eqcplxh2oid)) &
    deallocate(reaction%eqcplxh2oid)
  nullify(reaction%eqcplxh2oid)
  if (associated(reaction%eqcplxh2ostoich)) &
    deallocate(reaction%eqcplxh2ostoich)
  nullify(reaction%eqcplxh2ostoich)
  if (associated(reaction%eqcplx_a0)) &
    deallocate(reaction%eqcplx_a0)
  nullify(reaction%eqcplx_a0)
  if (associated(reaction%eqcplx_Z)) &
    deallocate(reaction%eqcplx_Z)
  nullify(reaction%eqcplx_Z)
  if (associated(reaction%eqcplx_molar_wt)) &
    deallocate(reaction%eqcplx_molar_wt)
  nullify(reaction%eqcplx_molar_wt)
  if (associated(reaction%eqcplx_logK)) &
    deallocate(reaction%eqcplx_logK)
  nullify(reaction%eqcplx_logK)
  if (associated(reaction%eqcplx_logKcoef)) &
    deallocate(reaction%eqcplx_logKcoef)
  nullify(reaction%eqcplx_logKcoef)
  
  if (associated(reaction%eqgasspecid)) &
    deallocate(reaction%eqgasspecid)
  nullify(reaction%eqgasspecid)
  if (associated(reaction%eqgasstoich)) &
    deallocate(reaction%eqgasstoich)
  nullify(reaction%eqgasstoich)
  if (associated(reaction%eqgash2oid)) &
    deallocate(reaction%eqgash2oid)
  nullify(reaction%eqgash2oid)
  if (associated(reaction%eqgash2ostoich)) &
    deallocate(reaction%eqgash2ostoich)
  nullify(reaction%eqgash2ostoich)
  if (associated(reaction%eqgas_logK)) &
    deallocate(reaction%eqgas_logK)
  nullify(reaction%eqgas_logK)
  if (associated(reaction%eqgas_logKcoef)) &
    deallocate(reaction%eqgas_logKcoef)
  nullify(reaction%eqgas_logKcoef)
  
  if (associated(reaction%eqionx_rxn_Z_flag)) &
    deallocate(reaction%eqionx_rxn_Z_flag)
  nullify(reaction%eqionx_rxn_Z_flag)
  if (associated(reaction%eqionx_rxn_cation_X_offset)) &
    deallocate(reaction%eqionx_rxn_cation_X_offset)
  nullify(reaction%eqionx_rxn_cation_X_offset)
  if (associated(reaction%eqionx_rxn_CEC)) &
    deallocate(reaction%eqionx_rxn_CEC)
  nullify(reaction%eqionx_rxn_CEC)
  if (associated(reaction%eqionx_rxn_k)) &
    deallocate(reaction%eqionx_rxn_k)
  nullify(reaction%eqionx_rxn_k)
  if (associated(reaction%eqionx_rxn_cationid)) &
    deallocate(reaction%eqionx_rxn_cationid)
  nullify(reaction%eqionx_rxn_cationid)

#if 0  
  if (associated(reaction%kinionx_CEC)) &
    deallocate(reaction%kinionx_CEC)
  nullify(reaction%kinionx_CEC)
  if (associated(reaction%kinionx_k)) &
    deallocate(reaction%kinionx_k)
  nullify(reaction%kinionx_k)
  if (associated(reaction%kinionx_cationid)) &
    deallocate(reaction%kinionx_cationid)
  nullify(reaction%kinionx_cationid)
#endif
  
  if (associated(reaction%eqsrfcplx_rxn_to_surf)) &
    deallocate(reaction%eqsrfcplx_rxn_to_surf)
  nullify(reaction%eqsrfcplx_rxn_to_surf)
  if (associated(reaction%eqsrfcplx_rxn_surf_type)) &
    deallocate(reaction%eqsrfcplx_rxn_surf_type)
  nullify(reaction%eqsrfcplx_rxn_surf_type)
  if (associated(reaction%eqsrfcplx_rxn_to_complex)) &
    deallocate(reaction%eqsrfcplx_rxn_to_complex)
  nullify(reaction%eqsrfcplx_rxn_to_complex)

  if (associated(reaction%eqsrfcplx_rxn_site_density)) &
    deallocate(reaction%eqsrfcplx_rxn_site_density)
  nullify(reaction%eqsrfcplx_rxn_site_density)
  if (associated(reaction%eqsrfcplx_rxn_stoich_flag)) &
    deallocate(reaction%eqsrfcplx_rxn_stoich_flag)
  nullify(reaction%eqsrfcplx_rxn_stoich_flag) 
! these deallocates above
!  if (associated(reaction%eqsrfcplx_site_names)) &
!    deallocate(reaction%eqsrfcplx_site_names)
!  nullify(reaction%eqsrfcplx_site_names)
!  if (associated(reaction%eqsrfcplx_names)) &
!    deallocate(reaction%eqsrfcplx_names)
!  nullify(reaction%eqsrfcplx_names)
  if (associated(reaction%eqsrfcplxspecid)) &
    deallocate(reaction%eqsrfcplxspecid)
  nullify(reaction%eqsrfcplxspecid)
  if (associated(reaction%eqsrfcplxstoich)) &
    deallocate(reaction%eqsrfcplxstoich)
  nullify(reaction%eqsrfcplxstoich)
  if (associated(reaction%eqsrfcplxh2oid)) &
    deallocate(reaction%eqsrfcplxh2oid)
  nullify(reaction%eqsrfcplxh2oid)
  if (associated(reaction%eqsrfcplxh2ostoich)) &
    deallocate(reaction%eqsrfcplxh2ostoich)
  nullify(reaction%eqsrfcplxh2ostoich)
!  if (associated(reaction%eqsrfcplx_mineral_id)) &
!    deallocate(reaction%eqsrfcplx_mineral_id)
!  nullify(reaction%eqsrfcplx_mineral_id)
  if (associated(reaction%eqsrfcplx_free_site_id)) &
    deallocate(reaction%eqsrfcplx_free_site_id)
  nullify(reaction%eqsrfcplx_free_site_id)
  if (associated(reaction%eqsrfcplx_free_site_stoich)) &
    deallocate(reaction%eqsrfcplx_free_site_stoich)
  nullify(reaction%eqsrfcplx_free_site_stoich)
  if (associated(reaction%eqsrfcplx_logK)) &
    deallocate(reaction%eqsrfcplx_logK)
  nullify(reaction%eqsrfcplx_logK)
  if (associated(reaction%eqsrfcplx_logKcoef)) &
    deallocate(reaction%eqsrfcplx_logKcoef)
  nullify(reaction%eqsrfcplx_logKcoef)
  if (associated(reaction%eqsrfcplx_Z)) &
    deallocate(reaction%eqsrfcplx_Z)
  nullify(reaction%eqsrfcplx_Z)

  if (associated(reaction%kinsrfcplx_rxn_to_surf)) &
    deallocate(reaction%kinsrfcplx_rxn_to_surf)
  nullify(reaction%kinsrfcplx_rxn_to_surf)
  if (associated(reaction%kinsrfcplx_rxn_surf_type)) &
    deallocate(reaction%kinsrfcplx_rxn_surf_type)
  nullify(reaction%kinsrfcplx_rxn_surf_type)
  if (associated(reaction%kinsrfcplx_rxn_to_complex)) &
    deallocate(reaction%kinsrfcplx_rxn_to_complex)
  nullify(reaction%kinsrfcplx_rxn_to_complex)
  if (associated(reaction%kinsrfcplx_rxn_to_site)) &
    deallocate(reaction%kinsrfcplx_rxn_to_site)
  nullify(reaction%kinsrfcplx_rxn_to_site)
  if (associated(reaction%kinsrfcplx_rxn_site_density)) &
    deallocate(reaction%kinsrfcplx_rxn_site_density)
  nullify(reaction%kinsrfcplx_rxn_site_density)
  if (associated(reaction%kinsrfcplx_rxn_stoich_flag)) &
    deallocate(reaction%kinsrfcplx_rxn_stoich_flag)
  nullify(reaction%kinsrfcplx_rxn_stoich_flag)
! these deallocates above
!  if (associated(reaction%kinsrfcplx_site_names)) &
!    deallocate(reaction%kinsrfcplx_site_names)
!  nullify(reaction%kinsrfcplx_site_names)
!  if (associated(reaction%kinsrfcplx_names)) &
!    deallocate(reaction%kinsrfcplx_names)
!  nullify(reaction%kinsrfcplx_names)
  if (associated(reaction%kinsrfcplxspecid)) &
    deallocate(reaction%kinsrfcplxspecid)
  nullify(reaction%kinsrfcplxspecid)
  if (associated(reaction%kinsrfcplxstoich)) &
    deallocate(reaction%kinsrfcplxstoich)
  nullify(reaction%kinsrfcplxstoich)
  if (associated(reaction%kinsrfcplx_free_site_stoich)) &
    deallocate(reaction%kinsrfcplx_free_site_stoich)
  nullify(reaction%kinsrfcplx_free_site_stoich)
  if (associated(reaction%kinsrfcplx_forward_rate)) &
    deallocate(reaction%kinsrfcplx_forward_rate)
  nullify(reaction%kinsrfcplx_forward_rate)
  if (associated(reaction%kinsrfcplx_backward_rate)) &
    deallocate(reaction%kinsrfcplx_backward_rate)
  nullify(reaction%kinsrfcplx_backward_rate)
!  if (associated(reaction%kinsrfcplx_logK)) &
!    deallocate(reaction%kinsrfcplx_logK)
!  nullify(reaction%kinsrfcplx_logK)
!  if (associated(reaction%kinsrfcplx_logKcoef)) &
!    deallocate(reaction%kinsrfcplx_logKcoef)
!  nullify(reaction%kinsrfcplx_logKcoef)
  if (associated(reaction%kinsrfcplx_Z)) &
    deallocate(reaction%kinsrfcplx_Z)
  nullify(reaction%kinsrfcplx_Z)
  
  if (associated(reaction%mnrlspecid)) deallocate(reaction%mnrlspecid)
  nullify(reaction%mnrlspecid)
  if (associated(reaction%mnrlstoich)) deallocate(reaction%mnrlstoich)
  nullify(reaction%mnrlstoich)
  if (associated(reaction%mnrlh2oid)) deallocate(reaction%mnrlh2oid)
  nullify(reaction%mnrlh2oid)
  if (associated(reaction%mnrlh2ostoich)) deallocate(reaction%mnrlh2ostoich)
  nullify(reaction%mnrlh2ostoich)
  if (associated(reaction%mnrl_logK)) deallocate(reaction%mnrl_logK)
  nullify(reaction%mnrl_logK)
  if (associated(reaction%mnrl_logKcoef)) deallocate(reaction%mnrl_logKcoef)
  nullify(reaction%mnrl_logKcoef)
  
  if (associated(reaction%kinmnrlspecid)) deallocate(reaction%kinmnrlspecid)
  nullify(reaction%kinmnrlspecid)
  if (associated(reaction%kinmnrlstoich)) deallocate(reaction%kinmnrlstoich)
  nullify(reaction%kinmnrlstoich)
  if (associated(reaction%kinmnrlh2oid)) deallocate(reaction%kinmnrlh2oid)
  nullify(reaction%kinmnrlh2oid)
  if (associated(reaction%kinmnrlh2ostoich)) &
    deallocate(reaction%kinmnrlh2ostoich)
  nullify(reaction%kinmnrlh2ostoich)
  if (associated(reaction%kinmnrl_logK)) deallocate(reaction%kinmnrl_logK)
  nullify(reaction%kinmnrl_logK)
  if (associated(reaction%kinmnrl_logKcoef)) &
    deallocate(reaction%kinmnrl_logKcoef)
  nullify(reaction%kinmnrl_logKcoef)
  if (associated(reaction%kinmnrl_rate)) deallocate(reaction%kinmnrl_rate)
  nullify(reaction%kinmnrl_rate)
  if (associated(reaction%kinmnrl_molar_vol)) &
    deallocate(reaction%kinmnrl_molar_vol)
  nullify(reaction%kinmnrl_molar_vol)  
   if (associated(reaction%kinmnrl_molar_wt)) &
    deallocate(reaction%kinmnrl_molar_wt)
  nullify(reaction%kinmnrl_molar_wt)  

  if (associated(reaction%kinmnrl_num_prefactors)) &
    deallocate(reaction%kinmnrl_num_prefactors)
  nullify(reaction%kinmnrl_num_prefactors)
  if (associated(reaction%kinmnrl_prefactor_id)) &
    deallocate(reaction%kinmnrl_prefactor_id)
  nullify(reaction%kinmnrl_prefactor_id)
  if (associated(reaction%kinmnrl_pref_alpha)) &
    deallocate(reaction%kinmnrl_pref_alpha)
  nullify(reaction%kinmnrl_pref_alpha)
  if (associated(reaction%kinmnrl_pref_beta)) &
    deallocate(reaction%kinmnrl_pref_beta)
  nullify(reaction%kinmnrl_pref_beta)
  if (associated(reaction%kinmnrl_pref_atten_coef)) &
    deallocate(reaction%kinmnrl_pref_atten_coef)
  nullify(reaction%kinmnrl_pref_atten_coef)
  if (associated(reaction%kinmnrl_pref_rate)) &
    deallocate(reaction%kinmnrl_pref_rate)
  nullify(reaction%kinmnrl_pref_rate)
  if (associated(reaction%kinmnrl_pref_activation_energy)) &
    deallocate(reaction%kinmnrl_pref_activation_energy)
  nullify(reaction%kinmnrl_pref_activation_energy)

  if (associated(reaction%kinmnrl_Tempkin_const)) &
    deallocate(reaction%kinmnrl_Tempkin_const)
  nullify(reaction%kinmnrl_Tempkin_const)
  if (associated(reaction%kinmnrl_affinity_power)) &
    deallocate(reaction%kinmnrl_affinity_power)
  nullify(reaction%kinmnrl_affinity_power)
  if (associated(reaction%kinmnrl_affinity_threshold)) &
    deallocate(reaction%kinmnrl_affinity_threshold)
  nullify(reaction%kinmnrl_affinity_threshold)
  if (associated(reaction%kinmnrl_activation_energy)) &
    deallocate(reaction%kinmnrl_activation_energy)
  nullify(reaction%kinmnrl_activation_energy)
  if (associated(reaction%kinmnrl_rate_limiter)) &
    deallocate(reaction%kinmnrl_rate_limiter)
  nullify(reaction%kinmnrl_rate_limiter)
  if (associated(reaction%kinmnrl_irreversible)) &
    deallocate(reaction%kinmnrl_irreversible)
  nullify(reaction%kinmnrl_irreversible)

  if (associated(reaction%kinmr_rate)) deallocate(reaction%kinmr_rate)
  nullify(reaction%kinmr_rate)

  if (associated(reaction%kinmr_frac)) deallocate(reaction%kinmr_frac)
  nullify(reaction%kinmr_frac)
  
  if (associated(reaction%pri_spec_to_coll_spec)) &
    deallocate(reaction%pri_spec_to_coll_spec)
  nullify(reaction%pri_spec_to_coll_spec)
  if (associated(reaction%coll_spec_to_pri_spec)) &
    deallocate(reaction%coll_spec_to_pri_spec)
  nullify(reaction%coll_spec_to_pri_spec)
  if (associated(reaction%colloid_mobile_fraction)) &
    deallocate(reaction%colloid_mobile_fraction)
  nullify(reaction%colloid_mobile_fraction)
  
  if (associated(reaction%generalspecid)) &
    deallocate(reaction%generalspecid)
  nullify(reaction%generalspecid)
  if (associated(reaction%generalstoich)) &
    deallocate(reaction%generalstoich)
  nullify(reaction%generalstoich)
  if (associated(reaction%generalforwardspecid)) &
    deallocate(reaction%generalforwardspecid)
  nullify(reaction%generalforwardspecid)
  if (associated(reaction%generalforwardstoich)) &
    deallocate(reaction%generalforwardstoich)
  nullify(reaction%generalforwardstoich)
  if (associated(reaction%generalbackwardspecid)) &
    deallocate(reaction%generalbackwardspecid)
  nullify(reaction%generalbackwardspecid)
  if (associated(reaction%generalbackwardstoich)) &
    deallocate(reaction%generalbackwardstoich)
  nullify(reaction%generalbackwardstoich)
  if (associated(reaction%generalh2oid)) &
    deallocate(reaction%generalh2oid)
  nullify(reaction%generalh2oid)
  if (associated(reaction%generalh2ostoich)) &
    deallocate(reaction%generalh2ostoich)
  nullify(reaction%generalh2ostoich)
  if (associated(reaction%general_kf)) &
    deallocate(reaction%general_kf)
  nullify(reaction%general_kf)
  if (associated(reaction%general_kr)) &
    deallocate(reaction%general_kr)
  nullify(reaction%general_kr)
  
  if (associated(reaction%eqkdspecid)) &
    deallocate(reaction%eqkdspecid)
  nullify(reaction%eqkdspecid)
  if (associated(reaction%eqkdtype)) &
    deallocate(reaction%eqkdtype)
  nullify(reaction%eqkdtype)
  if (associated(reaction%eqkdspecid)) &
    deallocate(reaction%eqkdspecid)
  nullify(reaction%eqkdspecid)
  if (associated(reaction%eqkddistcoef)) &
    deallocate(reaction%eqkddistcoef)
  nullify(reaction%eqkddistcoef)
  if (associated(reaction%eqkdlangmuirb)) &
    deallocate(reaction%eqkdlangmuirb)
  nullify(reaction%eqkdlangmuirb)
  if (associated(reaction%eqkdfreundlichn)) &
    deallocate(reaction%eqkdfreundlichn)
  nullify(reaction%eqkdfreundlichn)
  
  deallocate(reaction)
  nullify(reaction)

end subroutine ReactionDestroy

end module Reaction_Aux_module
