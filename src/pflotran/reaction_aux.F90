module Reaction_Aux_module
  
  implicit none

  private 

#include "definitions.h"
  
  type, public :: aq_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: a0
    PetscReal :: molar_weight
    PetscReal :: Z
    PetscTruth :: print_me
    type(equilibrium_rxn_type), pointer :: eqrxn
    type(aq_species_type), pointer :: next
  end type aq_species_type

  type, public :: gas_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    PetscTruth :: print_me
    type(equilibrium_rxn_type), pointer :: eqrxn
    type(gas_species_type), pointer :: next    
  end type gas_species_type

  type, public :: equilibrium_rxn_type
    PetscInt :: nspec
    character(len=MAXWORDLENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK(:)
  end type equilibrium_rxn_type

  type, public :: kinetic_rxn_type
    PetscInt :: nspec
    character(len=MAXWORDLENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK(:)
    PetscReal :: rate_forward
    PetscReal :: rate_reverse
  end type kinetic_rxn_type

  type, public :: mineral_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    PetscTruth :: print_me
    type(transition_state_rxn_type), pointer :: tstrxn
    type(mineral_type), pointer :: next
  end type mineral_type

  type, public :: transition_state_rxn_type
    PetscInt :: nspec
    character(len=MAXWORDLENGTH), pointer :: spec_name(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: stoich(:)
    PetscReal, pointer :: logK(:)
    PetscInt :: nspec_primary_prefactor
    character(len=MAXWORDLENGTH), pointer :: spec_name_primary_prefactor(:)
    PetscInt, pointer :: spec_ids_primary_prefactor(:)
    PetscReal, pointer :: stoich_primary_prefactor(:)
    PetscInt :: nspec_secondary_prefactor
    character(len=MAXWORDLENGTH), pointer :: spec_name_secondary_prefactor(:)
    PetscInt, pointer :: spec_ids_secondary_prefactor(:)
    PetscReal, pointer :: stoich_secondary_prefactor(:)
    PetscReal :: affinity_factor_sigma
    PetscReal :: affinity_factor_beta
    PetscReal :: rate
  end type transition_state_rxn_type
  
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
    PetscTruth :: print_me
    type(equilibrium_rxn_type), pointer :: eqrxn
    type(surface_complex_type), pointer :: next
  end type surface_complex_type

  type, public :: surface_complexation_rxn_type
    PetscInt :: id
    PetscInt :: free_site_id
    character(len=MAXWORDLENGTH) :: free_site_name
    PetscTruth :: free_site_print_me
    PetscInt :: mineral_id
    character(len=MAXWORDLENGTH) :: mineral_name
    PetscReal :: site_density
    type(surface_complex_type), pointer :: complex_list
    type (surface_complexation_rxn_type), pointer :: next
  end type surface_complexation_rxn_type    

  type, public :: aq_species_constraint_type
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscReal, pointer :: basis_molarity(:)
    PetscInt, pointer :: constraint_type(:)
    PetscInt, pointer :: constraint_spec_id(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_spec_name(:)
  end type aq_species_constraint_type

  type, public :: mineral_constraint_type
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_vol_frac(:)
    PetscReal, pointer :: constraint_area(:)
    PetscReal, pointer :: basis_vol_frac(:)
    PetscReal, pointer :: basis_area(:)
  end type mineral_constraint_type

  type, public :: reaction_type
    character(len=MAXSTRINGLENGTH) :: database_filename
    PetscTruth :: use_full_geochemistry
    PetscTruth :: use_log_formulation ! flag for solving for the change in the log of the concentration
    PetscTruth :: print_all_species
    PetscTruth :: print_pH
    PetscTruth :: print_kd
    PetscTruth :: print_total_sorb
    PetscTruth :: print_act_coefs
    PetscTruth :: print_total_component
    PetscTruth :: print_free_ion
    PetscInt :: num_dbase_temperatures
    PetscInt :: h_ion_id
    PetscInt :: na_ion_id
    PetscInt :: cl_ion_id
    PetscInt :: o2_gas_id
    PetscInt :: co2_gas_id
    PetscReal, pointer :: dbase_temperatures(:)
    type(aq_species_type), pointer :: primary_species_list
    type(aq_species_type), pointer :: secondary_species_list
    type(gas_species_type), pointer :: gas_species_list
    type(mineral_type), pointer :: mineral_list
    type(ion_exchange_rxn_type), pointer :: ion_exchange_rxn_list
    type(surface_complexation_rxn_type), pointer :: surface_complexation_rxn_list
    PetscInt :: act_coef_update_frequency
    PetscInt :: act_coef_update_algorithm
    PetscTruth :: checkpoint_activity_coefs
    PetscTruth :: act_coef_use_bdot
    PetscTruth :: use_activity_h2o
    
    ! compressed arrays for efficient computation
    ! primary aqueous complexes
    PetscInt :: ncomp
    character(len=MAXWORDLENGTH), pointer :: primary_species_names(:)
    PetscTruth, pointer :: primary_species_print(:)
    PetscReal, pointer :: primary_spec_a0(:)
    PetscReal, pointer :: primary_spec_Z(:)
    
    ! aqueous complexes
    PetscInt :: neqcmplx
    character(len=MAXWORDLENGTH), pointer :: secondary_species_names(:)
    PetscTruth, pointer :: secondary_species_print(:)
    character(len=MAXWORDLENGTH), pointer :: eqcmplx_basis_names(:,:)
    PetscTruth, pointer :: eqcmplx_basis_print(:)
    PetscInt, pointer :: eqcmplxspecid(:,:)   ! (0:ncomp in rxn)
    PetscReal, pointer :: eqcmplxstoich(:,:)
    PetscInt, pointer :: eqcmplxh2oid(:)       ! id of water, if present
    PetscReal, pointer :: eqcmplxh2ostoich(:)  ! stoichiometry of water, if present
    PetscReal, pointer :: eqcmplx_a0(:)  ! Debye-Huckel constant
    PetscReal, pointer :: eqcmplx_Z(:)
    PetscReal, pointer :: eqcmplx_logK(:)
    PetscReal, pointer :: eqcmplx_logKcoef(:,:)
    ! Debye-Huckel
    PetscReal :: debyeA  ! Debye-Huckel A coefficient
    PetscReal :: debyeB  ! Debye-Huckel B coefficient
    PetscReal :: debyeBdot  ! Debye-Huckel Bdot coefficient
    
    ! gas species
    PetscInt :: ngas
    character(len=MAXWORDLENGTH), pointer :: gas_species_names(:)
    PetscTruth, pointer :: gas_species_print(:)
    PetscInt, pointer :: eqgasspecid(:,:)   ! (0:ncomp in rxn)
    PetscReal, pointer :: eqgasstoich(:,:)
    PetscInt, pointer :: eqgash2oid(:)       ! id of water, if present
    PetscReal, pointer :: eqgash2ostoich(:)  ! stoichiometry of water, if present
    PetscReal, pointer :: eqgas_logK(:)
    PetscReal, pointer :: eqgas_logKcoef(:,:)
    
    PetscInt :: neqsorb
    PetscTruth, pointer :: kd_print(:)
    PetscTruth, pointer :: total_sorb_print(:)

    ! ionx exchange reactions
    PetscInt :: neqionxrxn
    PetscInt :: neqionxcation 
    PetscTruth, pointer :: eqionx_rxn_Z_flag(:)
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
    PetscInt :: neqsurfcmplx
    PetscInt :: neqsurfcmplxrxn 
    PetscInt, pointer :: eqsurfcmplx_rxn_to_mineral(:)
    PetscInt, pointer :: eqsurfcmplx_rxn_to_complex(:,:)
    PetscReal, pointer :: eqsurfcmplx_rxn_site_density(:)
    PetscTruth, pointer :: eqsurfcmplx_rxn_stoich_flag(:)
    character(len=MAXWORDLENGTH), pointer :: surface_site_names(:)
    PetscTruth, pointer :: surface_site_print(:)
    character(len=MAXWORDLENGTH), pointer :: surface_complex_names(:)
    PetscTruth, pointer :: surface_complex_print(:)
    PetscInt, pointer :: eqsurfcmplxspecid(:,:)
    PetscReal, pointer :: eqsurfcmplxstoich(:,:)
    PetscInt, pointer :: eqsurfcmplxh2oid(:)
    PetscReal, pointer :: eqsurfcmplxh2ostoich(:)
    PetscInt, pointer :: eqsurfcmplx_free_site_id(:)
    PetscReal, pointer :: eqsurfcmplx_free_site_stoich(:)
    PetscInt, pointer :: eqsurfcmplx_mineral_id(:)
    PetscReal, pointer :: eqsurfcmplx_logK(:)
    PetscReal, pointer :: eqsurfcmplx_logKcoef(:,:)
    PetscReal, pointer :: eqsurfcmplx_Z(:)  ! valence

#if 0    
    PetscInt, pointer :: kinsurfcmplxspecid(:,:)
    PetscReal, pointer :: kinsurfcmplxstoich(:,:)
    PetscInt, pointer :: kinsurfcmplxh2oid(:)
    PetscReal, pointer :: kinsurfcmplxh2ostoich(:)
    PetscReal, pointer :: kinsurfcmplx_freesite_stoich(:)
    PetscReal, pointer :: kinsurfcmplx_logK(:)
    PetscReal, pointer :: kinsurfcmplx_logKcoef(:,:)
    PetscReal, pointer :: kinsurfcmplx_Z(:)  ! valence
#endif

    ! multirate reaction rates
    PetscInt :: kinmr_nrate 
    PetscReal, pointer :: kinmr_rate(:)
    
    ! mineral reactions
    PetscInt :: nmnrl
    character(len=MAXWORDLENGTH), pointer :: mineral_names(:)
    
      ! for saturation states
    PetscInt, pointer :: mnrlspecid(:,:)
    PetscReal, pointer :: mnrlstoich(:,:)
    PetscInt, pointer :: mnrlh2oid(:)
    PetscReal, pointer :: mnrlh2ostoich(:)
    PetscReal, pointer :: mnrl_logK(:)
    PetscReal, pointer :: mnrl_logKcoef(:,:)
    
      ! for kinetic reactions
    PetscInt :: nkinmnrl
    character(len=MAXWORDLENGTH), pointer :: kinmnrl_names(:)
    PetscTruth, pointer :: kinmnrl_print(:)
    PetscInt, pointer :: kinmnrlspecid(:,:)
    PetscReal, pointer :: kinmnrlstoich(:,:)
    PetscInt, pointer :: kinmnrlh2oid(:)
    PetscReal, pointer :: kinmnrlh2ostoich(:)
    PetscReal, pointer :: kinmnrl_logK(:)
    PetscReal, pointer :: kinmnrl_logKcoef(:,:)
    PetscReal, pointer :: kinmnrl_rate(:,:)
    PetscReal, pointer :: kinmnrl_molar_vol(:)
    PetscInt, pointer :: kinmnrl_num_prefactors(:)
    PetscInt, pointer :: kinmnrl_pri_prefactor_id(:,:,:)
    PetscReal, pointer :: kinmnrl_pri_pref_alpha_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_pri_pref_beta_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_pri_pref_atten_coef(:,:,:)
    PetscInt, pointer :: kinmnrl_sec_prefactor_id(:,:,:)
    PetscReal, pointer :: kinmnrl_sec_pref_alpha_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_sec_pref_beta_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_sec_pref_atten_coef(:,:,:)
    PetscReal, pointer :: kinmnrl_Tempkin_const(:)
    PetscReal, pointer :: kinmnrl_affinity_power(:)
    
    PetscReal :: max_dlnC
    
  end type reaction_type

  public :: ReactionCreate, &
            AqueousSpeciesCreate, &
            GasSpeciesCreate, &
            MineralCreate, &
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
            GetKineticMineralCount, &
            EquilibriumRxnCreate, &
            EquilibriumRxnDestroy, &
            TransitionStateTheoryRxnCreate, &
            TransitionStateTheoryRxnDestroy, &
            SurfaceComplexationRxnCreate, &
            AqueousSpeciesConstraintDestroy, &
            MineralConstraintDestroy, &
            AqueousSpeciesConstraintCreate, &
            MineralConstraintCreate, &
            SurfaceComplexCreate, &
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
  reaction%print_all_species = PETSC_TRUE
  reaction%print_pH = PETSC_FALSE
  reaction%print_kd = PETSC_FALSE
  reaction%print_total_sorb = PETSC_FALSE
  reaction%print_act_coefs = PETSC_FALSE
  reaction%use_log_formulation = PETSC_FALSE
  reaction%use_full_geochemistry = PETSC_FALSE
  reaction%use_activity_h2o = PETSC_FALSE
  reaction%print_total_component = PETSC_TRUE
  reaction%print_free_ion = PETSC_FALSE
  
  reaction%h_ion_id = 0
  reaction%na_ion_id = 0
  reaction%cl_ion_id = 0
  reaction%o2_gas_id = 0
  reaction%co2_gas_id = 0

  nullify(reaction%primary_species_list)
  nullify(reaction%secondary_species_list)
  nullify(reaction%gas_species_list)
  nullify(reaction%mineral_list)
  nullify(reaction%ion_exchange_rxn_list)
  nullify(reaction%surface_complexation_rxn_list)
  
  nullify(reaction%primary_species_names)
  nullify(reaction%secondary_species_names)
  nullify(reaction%eqcmplx_basis_names)
  nullify(reaction%gas_species_names)
  nullify(reaction%surface_site_names)
  nullify(reaction%surface_complex_names)
  nullify(reaction%mineral_names)
  nullify(reaction%kinmnrl_names)

  nullify(reaction%primary_species_print)
  nullify(reaction%secondary_species_print)
  nullify(reaction%eqcmplx_basis_print)
  nullify(reaction%gas_species_print)
  nullify(reaction%surface_site_print)
  nullify(reaction%surface_complex_print)
  nullify(reaction%kinmnrl_print)
  nullify(reaction%kd_print)
  nullify(reaction%total_sorb_print)
  
  reaction%ncomp = 0
  nullify(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_Z)

  reaction%ngas = 0
  nullify(reaction%eqgasspecid)
  nullify(reaction%eqgasstoich)
  nullify(reaction%eqgash2oid)
  nullify(reaction%eqgash2ostoich)
  nullify(reaction%eqgas_logK)
  nullify(reaction%eqgas_logKcoef)
  
  reaction%neqcmplx = 0
  nullify(reaction%eqcmplxspecid)
  nullify(reaction%eqcmplxstoich)
  nullify(reaction%eqcmplxh2oid)
  nullify(reaction%eqcmplxh2ostoich)
  nullify(reaction%eqcmplx_a0)
  nullify(reaction%eqcmplx_Z)
  nullify(reaction%eqcmplx_logK)
  nullify(reaction%eqcmplx_logKcoef)
  
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
  reaction%neqsurfcmplx = 0
  reaction%neqsurfcmplxrxn = 0
  nullify(reaction%eqsurfcmplx_rxn_to_mineral)
  nullify(reaction%eqsurfcmplx_rxn_to_complex)
  nullify(reaction%eqsurfcmplx_rxn_site_density)
  nullify(reaction%eqsurfcmplx_rxn_stoich_flag) 
  nullify(reaction%surface_site_names)
  nullify(reaction%surface_complex_names)
  nullify(reaction%eqsurfcmplxspecid)
  nullify(reaction%eqsurfcmplxstoich)
  nullify(reaction%eqsurfcmplxh2oid)
  nullify(reaction%eqsurfcmplxh2ostoich)
  nullify(reaction%eqsurfcmplx_mineral_id)
  nullify(reaction%eqsurfcmplx_free_site_id)
  nullify(reaction%eqsurfcmplx_free_site_stoich)
  nullify(reaction%eqsurfcmplx_logK)
  nullify(reaction%eqsurfcmplx_logKcoef)
  nullify(reaction%eqsurfcmplx_Z)

#if 0  
  nullify(reaction%kinsurfcmplxspecid)
  nullify(reaction%kinsurfcmplxstoich)
  nullify(reaction%kinsurfcmplxh2oid)
  nullify(reaction%kinsurfcmplxh2ostoich)
  nullify(reaction%kinsurfcmplx_freesite_stoich)
  nullify(reaction%kinsurfcmplx_logK)
  nullify(reaction%kinsurfcmplx_logKcoef)
  nullify(reaction%kinsurfcmplx_Z)
#endif

  reaction%kinmr_nrate = 0
  nullify(reaction%kinmr_rate)

  reaction%nmnrl = 0  
  nullify(reaction%mnrlspecid)
  nullify(reaction%mnrlstoich)
  nullify(reaction%mnrlh2oid)
  nullify(reaction%mnrlstoich)
  nullify(reaction%mnrlh2ostoich)
  nullify(reaction%mnrl_logKcoef)
  
  reaction%nkinmnrl = 0  
  nullify(reaction%kinmnrlspecid)
  nullify(reaction%kinmnrlstoich)
  nullify(reaction%kinmnrlh2oid)
  nullify(reaction%kinmnrlh2ostoich)
  nullify(reaction%kinmnrl_logK)
  nullify(reaction%kinmnrl_logKcoef)
  nullify(reaction%kinmnrl_rate)
  nullify(reaction%kinmnrl_molar_vol)
  nullify(reaction%kinmnrl_num_prefactors)
  nullify(reaction%kinmnrl_pri_prefactor_id)
  nullify(reaction%kinmnrl_pri_pref_alpha_stoich)
  nullify(reaction%kinmnrl_pri_pref_beta_stoich)
  nullify(reaction%kinmnrl_pri_pref_atten_coef)
  nullify(reaction%kinmnrl_sec_prefactor_id)
  nullify(reaction%kinmnrl_sec_pref_alpha_stoich)
  nullify(reaction%kinmnrl_sec_pref_beta_stoich)
  nullify(reaction%kinmnrl_sec_pref_atten_coef)
  nullify(reaction%kinmnrl_Tempkin_const)
  nullify(reaction%kinmnrl_affinity_power)
  
  reaction%max_dlnC = 5.d0

  ReactionCreate => reaction
  
end function ReactionCreate

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
  nullify(species%eqrxn)
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
  nullify(species%eqrxn)
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
! EquilibriumRxnCreate: Allocate and initialize an equilibrium reaction
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
function EquilibriumRxnCreate()

  implicit none
    
  type(equilibrium_rxn_type), pointer :: EquilibriumRxnCreate

  type(equilibrium_rxn_type), pointer :: eqrxn

  allocate(eqrxn)
  eqrxn%nspec = 0
  nullify(eqrxn%spec_name)
  nullify(eqrxn%stoich)
  nullify(eqrxn%spec_ids)
  nullify(eqrxn%logK)
  
  EquilibriumRxnCreate => eqrxn
  
end function EquilibriumRxnCreate

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
  tstrxn%nspec = 0
  nullify(tstrxn%spec_name)
  nullify(tstrxn%stoich)
  nullify(tstrxn%spec_ids)
  nullify(tstrxn%logK)
  tstrxn%nspec_primary_prefactor = 0
  nullify(tstrxn%spec_name_primary_prefactor)
  nullify(tstrxn%spec_ids_primary_prefactor)
  nullify(tstrxn%stoich_primary_prefactor)
  tstrxn%nspec_secondary_prefactor = 0
  nullify(tstrxn%spec_name_secondary_prefactor)
  nullify(tstrxn%spec_ids_secondary_prefactor)
  nullify(tstrxn%stoich_secondary_prefactor)
  tstrxn%affinity_factor_sigma = 0.d0
  tstrxn%affinity_factor_beta = 0.d0
  tstrxn%rate = 0.d0
  
  TransitionStateTheoryRxnCreate => tstrxn
  
end function TransitionStateTheoryRxnCreate

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

  type(surface_complexation_rxn_type), pointer :: surfcplxrxn
  
  allocate(surfcplxrxn)
  surfcplxrxn%free_site_id = 0
  surfcplxrxn%free_site_name = ''
  surfcplxrxn%free_site_print_me = PETSC_FALSE

  surfcplxrxn%mineral_id = 0
  surfcplxrxn%mineral_name = ''
  surfcplxrxn%site_density = 0.d0
  
  nullify(surfcplxrxn%complex_list)
  nullify(surfcplxrxn%next)
  
  SurfaceComplexationRxnCreate => surfcplxrxn
  
end function SurfaceComplexationRxnCreate

! ************************************************************************** !
!
! SurfaceComplexCreate: Allocate and initialize a surface complexreaction
! author: Peter Lichtner
! date: 10/21/08
!
! ************************************************************************** !
function SurfaceComplexCreate()

  implicit none
    
  type(surface_complex_type), pointer :: SurfaceComplexCreate

  type(surface_complex_type), pointer :: srfcmplx
  
  allocate(srfcmplx)
  srfcmplx%id = 0
  srfcmplx%name = ''
  srfcmplx%Z = 0.d0
  srfcmplx%free_site_stoich = 0.d0
  srfcmplx%print_me = PETSC_FALSE
  nullify(srfcmplx%eqrxn)
  nullify(srfcmplx%next)
  
  SurfaceComplexCreate => srfcmplx
  
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
  allocate(constraint%names(reaction%ncomp))
  constraint%names = ''
  allocate(constraint%constraint_conc(reaction%ncomp))
  constraint%constraint_conc = 0.d0
  allocate(constraint%basis_molarity(reaction%ncomp))
  constraint%basis_molarity = 0.d0
  allocate(constraint%constraint_spec_id(reaction%ncomp))
  constraint%constraint_spec_id = 0
  allocate(constraint%constraint_type(reaction%ncomp))
  constraint%constraint_type = 0
  allocate(constraint%constraint_spec_name(reaction%ncomp))
  constraint%constraint_spec_name = ''

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

  MineralConstraintCreate => constraint

end function MineralConstraintCreate
  
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
  
  integer :: GetPrimarySpeciesCount
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
  
  integer :: GetSecondarySpeciesCount
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
  
  integer :: GetGasCount
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
  
  integer :: GetMineralCount
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
! GetKineticMineralCount: Returns the number of kinetic minerals
! author: Glenn Hammond
! date: 11/04/08
!
! ************************************************************************** !
function GetKineticMineralCount(reaction)

  implicit none
  
  integer :: GetKineticMineralCount
  type(reaction_type) :: reaction

  type(mineral_type), pointer :: mineral

  GetKineticMineralCount = 0
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    if (mineral%itype == MINERAL_KINETIC) &
      GetKineticMineralCount = GetKineticMineralCount + 1
    mineral => mineral%next
  enddo

end function GetKineticMineralCount

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
! AqueousSpeciesDestroy: Deallocates an aqueous species
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine AqueousSpeciesDestroy(species)

  implicit none
    
  type(aq_species_type), pointer :: species

  if (associated(species%eqrxn)) call EquilibriumRxnDestroy(species%eqrxn)
  deallocate(species)  
  nullify(species)

end subroutine AqueousSpeciesDestroy

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

  if (associated(species%eqrxn)) call EquilibriumRxnDestroy(species%eqrxn)
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

  if (associated(mineral%tstrxn)) &
    call TransitionStateTheoryRxnDestroy(mineral%tstrxn)
  deallocate(mineral)  
  nullify(mineral)

end subroutine MineralDestroy

! ************************************************************************** !
!
! EquilibriumRxnDestroy: Deallocates an equilibrium reaction
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine EquilibriumRxnDestroy(eqrxn)

  implicit none
    
  type(equilibrium_rxn_type), pointer :: eqrxn

  if (.not.associated(eqrxn)) return
  
  if (associated(eqrxn%spec_name)) deallocate(eqrxn%spec_name)
  nullify(eqrxn%spec_name)
  if (associated(eqrxn%spec_ids)) deallocate(eqrxn%spec_ids)
  nullify(eqrxn%spec_ids)
  if (associated(eqrxn%stoich)) deallocate(eqrxn%stoich)
  nullify(eqrxn%stoich)
  if (associated(eqrxn%logK)) deallocate(eqrxn%logK)
  nullify(eqrxn%logK)

  deallocate(eqrxn)  
  nullify(eqrxn)

end subroutine EquilibriumRxnDestroy

! ************************************************************************** !
!
! TransitionStateTheoryRxnDestroy: Deallocates a transition state reaction
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine TransitionStateTheoryRxnDestroy(tstrxn)

  implicit none
    
  type(transition_state_rxn_type), pointer :: tstrxn

  if (.not.associated(tstrxn)) return
  
  if (associated(tstrxn%spec_name)) deallocate(tstrxn%spec_name)
  nullify(tstrxn%spec_name)
  if (associated(tstrxn%spec_ids)) deallocate(tstrxn%spec_ids)
  nullify(tstrxn%spec_ids)
  if (associated(tstrxn%stoich)) deallocate(tstrxn%stoich)
  nullify(tstrxn%stoich)
  if (associated(tstrxn%spec_name_primary_prefactor)) deallocate(tstrxn%spec_name_primary_prefactor)
  nullify(tstrxn%spec_name_primary_prefactor)
  if (associated(tstrxn%spec_ids_primary_prefactor)) deallocate(tstrxn%spec_ids_primary_prefactor)
  nullify(tstrxn%spec_ids_primary_prefactor)
  if (associated(tstrxn%stoich_primary_prefactor)) deallocate(tstrxn%stoich_primary_prefactor)
  nullify(tstrxn%stoich_primary_prefactor)
  if (associated(tstrxn%spec_name_secondary_prefactor)) deallocate(tstrxn%spec_name_secondary_prefactor)
  nullify(tstrxn%spec_name_secondary_prefactor)
  if (associated(tstrxn%spec_ids_secondary_prefactor)) deallocate(tstrxn%spec_ids_secondary_prefactor)
  nullify(tstrxn%spec_ids_secondary_prefactor)
  if (associated(tstrxn%stoich_secondary_prefactor)) deallocate(tstrxn%stoich_secondary_prefactor)
  nullify(tstrxn%stoich_secondary_prefactor)
  if (associated(tstrxn%logK)) deallocate(tstrxn%logK)
  nullify(tstrxn%logK)

  deallocate(tstrxn)  
  nullify(tstrxn)

end subroutine TransitionStateTheoryRxnDestroy

! ************************************************************************** !
!
! SurfaceComplexationRxnDestroy: Deallocates a surface complexation reaction
! author: Glenn Hammond
! date: 10/21/08
!
! ************************************************************************** !
subroutine SurfaceComplexationRxnDestroy(surfcplxrxn)

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: surfcplxrxn

  type(surface_complex_type), pointer :: cur_surfcplx, prev_surfcplx
  
  if (.not.associated(surfcplxrxn)) return
  
  cur_surfcplx => surfcplxrxn%complex_list
  do
    if (.not.associated(cur_surfcplx)) exit
    prev_surfcplx => cur_surfcplx
    cur_surfcplx => cur_surfcplx%next
    call SurfaceComplexDestroy(prev_surfcplx)
    nullify(prev_surfcplx)
  enddo
  
  deallocate(surfcplxrxn)  
  nullify(surfcplxrxn)

end subroutine SurfaceComplexationRxnDestroy

! ************************************************************************** !
!
! SurfaceComplexDestroy: Deallocates a surface complex
! author: Glenn Hammond
! date: 10/21/08
!
! ************************************************************************** !
subroutine SurfaceComplexDestroy(surfcplx)

  implicit none
    
  type(surface_complex_type), pointer :: surfcplx

  if (.not.associated(surfcplx)) return
  
  if (associated(surfcplx%eqrxn)) &
    call EquilibriumRxnDestroy(surfcplx%eqrxn)
  nullify(surfcplx%eqrxn)
  nullify(surfcplx%next)

  deallocate(surfcplx)  
  nullify(surfcplx)

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
  if (associated(constraint%constraint_spec_name)) &
    deallocate(constraint%constraint_spec_name)
  nullify(constraint%constraint_spec_name)

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

  deallocate(constraint)
  nullify(constraint)

end subroutine MineralConstraintDestroy

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
  type(ion_exchange_rxn_type), pointer :: ionxrxn, prev_ionxrxn
  type(surface_complexation_rxn_type), pointer :: surfcplxrxn, prev_surfcplxrxn

  if (.not.associated(reaction)) return 

  ! primary species
  aq_species => reaction%primary_species_list
  do
    if (.not.associated(aq_species)) exit
    prev_aq_species => aq_species
    aq_species => aq_species%next
    call AqueousSpeciesDestroy(prev_aq_species)
  enddo  
  nullify(reaction%primary_species_list)

  ! secondary species
  aq_species => reaction%secondary_species_list
  do
    if (.not.associated(aq_species)) exit
    prev_aq_species => aq_species
    aq_species => aq_species%next
    call AqueousSpeciesDestroy(prev_aq_species)
  enddo  
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
  surfcplxrxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(surfcplxrxn)) exit
    prev_surfcplxrxn => surfcplxrxn
    surfcplxrxn => surfcplxrxn%next
    call SurfaceComplexationRxnDestroy(prev_surfcplxrxn)
  enddo    
  nullify(reaction%surface_complexation_rxn_list)
  
  if (associated(reaction%primary_species_names)) deallocate(reaction%primary_species_names)
  nullify(reaction%primary_species_names)
  if (associated(reaction%secondary_species_names)) deallocate(reaction%secondary_species_names)
  nullify(reaction%secondary_species_names)
  if (associated(reaction%gas_species_names)) deallocate(reaction%gas_species_names)
  nullify(reaction%gas_species_names)
  if (associated(reaction%surface_site_names)) deallocate(reaction%surface_site_names)
  nullify(reaction%surface_site_names)
  if (associated(reaction%surface_complex_names)) deallocate(reaction%surface_complex_names)
  nullify(reaction%surface_complex_names)
  if (associated(reaction%mineral_names)) deallocate(reaction%mineral_names)
  nullify(reaction%mineral_names)
  if (associated(reaction%kinmnrl_names)) deallocate(reaction%kinmnrl_names)
  nullify(reaction%kinmnrl_names)

  if (associated(reaction%primary_species_print)) deallocate(reaction%primary_species_print)
  nullify(reaction%primary_species_print)
  if (associated(reaction%secondary_species_print)) deallocate(reaction%secondary_species_print)
  nullify(reaction%primary_species_print)
  if (associated(reaction%eqcmplx_basis_print)) deallocate(reaction%eqcmplx_basis_print)
  nullify(reaction%eqcmplx_basis_print)
  if (associated(reaction%gas_species_print)) deallocate(reaction%gas_species_print)
  nullify(reaction%gas_species_print)
  if (associated(reaction%surface_site_print)) deallocate(reaction%surface_site_print)
  nullify(reaction%primary_species_print)
  if (associated(reaction%surface_complex_print)) deallocate(reaction%surface_complex_print)
  nullify(reaction%surface_complex_print)
  if (associated(reaction%kinmnrl_print)) deallocate(reaction%kinmnrl_print)
  nullify(reaction%kinmnrl_print)
  if (associated(reaction%kd_print)) deallocate(reaction%kd_print)
  nullify(reaction%kd_print)
  if (associated(reaction%total_sorb_print)) deallocate(reaction%total_sorb_print)
  nullify(reaction%total_sorb_print)
    
  if (associated(reaction%primary_spec_a0)) deallocate(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_a0)
  if (associated(reaction%primary_spec_Z)) deallocate(reaction%primary_spec_Z)
  nullify(reaction%primary_spec_Z)
  
  if (associated(reaction%eqcmplxspecid)) deallocate(reaction%eqcmplxspecid)
  nullify(reaction%eqcmplxspecid)
  if (associated(reaction%eqcmplxstoich)) deallocate(reaction%eqcmplxstoich)
  nullify(reaction%eqcmplxstoich)
  if (associated(reaction%eqcmplxh2oid)) deallocate(reaction%eqcmplxh2oid)
  nullify(reaction%eqcmplxh2oid)
  if (associated(reaction%eqcmplxh2ostoich)) deallocate(reaction%eqcmplxh2ostoich)
  nullify(reaction%eqcmplxh2ostoich)
  if (associated(reaction%eqcmplx_a0)) deallocate(reaction%eqcmplx_a0)
  nullify(reaction%eqcmplx_a0)
  if (associated(reaction%eqcmplx_Z)) deallocate(reaction%eqcmplx_Z)
  nullify(reaction%eqcmplx_Z)
  if (associated(reaction%eqcmplx_logK)) deallocate(reaction%eqcmplx_logK)
  nullify(reaction%eqcmplx_logK)
  if (associated(reaction%eqcmplx_logKcoef)) deallocate(reaction%eqcmplx_logKcoef)
  nullify(reaction%eqcmplx_logKcoef)
  
  if (associated(reaction%eqgasspecid)) deallocate(reaction%eqgasspecid)
  nullify(reaction%eqgasspecid)
  if (associated(reaction%eqgasstoich)) deallocate(reaction%eqgasstoich)
  nullify(reaction%eqgasstoich)
  if (associated(reaction%eqgash2oid)) deallocate(reaction%eqgash2oid)
  nullify(reaction%eqgash2oid)
  if (associated(reaction%eqgash2ostoich)) deallocate(reaction%eqgash2ostoich)
  nullify(reaction%eqgash2ostoich)
  if (associated(reaction%eqgas_logK)) deallocate(reaction%eqgas_logK)
  nullify(reaction%eqgas_logK)
  if (associated(reaction%eqgas_logKcoef)) deallocate(reaction%eqgas_logKcoef)
  nullify(reaction%eqgas_logKcoef)
  
  if (associated(reaction%eqionx_rxn_Z_flag)) deallocate(reaction%eqionx_rxn_Z_flag)
  nullify(reaction%eqionx_rxn_Z_flag)
  if (associated(reaction%eqionx_rxn_cation_X_offset)) deallocate(reaction%eqionx_rxn_cation_X_offset)
  nullify(reaction%eqionx_rxn_cation_X_offset)
  if (associated(reaction%eqionx_rxn_CEC)) deallocate(reaction%eqionx_rxn_CEC)
  nullify(reaction%eqionx_rxn_CEC)
  if (associated(reaction%eqionx_rxn_k)) deallocate(reaction%eqionx_rxn_k)
  nullify(reaction%eqionx_rxn_k)
  if (associated(reaction%eqionx_rxn_cationid)) deallocate(reaction%eqionx_rxn_cationid)
  nullify(reaction%eqionx_rxn_cationid)

#if 0  
  if (associated(reaction%kinionx_CEC)) deallocate(reaction%kinionx_CEC)
  nullify(reaction%kinionx_CEC)
  if (associated(reaction%kinionx_k)) deallocate(reaction%kinionx_k)
  nullify(reaction%kinionx_k)
  if (associated(reaction%kinionx_cationid)) deallocate(reaction%kinionx_cationid)
  nullify(reaction%kinionx_cationid)
#endif
  
  if (associated(reaction%eqsurfcmplx_rxn_to_mineral)) deallocate(reaction%eqsurfcmplx_rxn_to_mineral)
  nullify(reaction%eqsurfcmplx_rxn_to_mineral)
  if (associated(reaction%eqsurfcmplx_rxn_to_complex)) deallocate(reaction%eqsurfcmplx_rxn_to_complex)
  nullify(reaction%eqsurfcmplx_rxn_to_complex)

  if (associated(reaction%eqsurfcmplx_rxn_site_density)) deallocate(reaction%eqsurfcmplx_rxn_site_density)
  nullify(reaction%eqsurfcmplx_rxn_site_density)
  if (associated(reaction%eqsurfcmplx_rxn_stoich_flag)) deallocate(reaction%eqsurfcmplx_rxn_stoich_flag)
  nullify(reaction%eqsurfcmplx_rxn_stoich_flag) 

  if (associated(reaction%surface_site_names)) deallocate(reaction%surface_site_names)
  nullify(reaction%surface_site_names)
  if (associated(reaction%surface_complex_names)) deallocate(reaction%surface_complex_names)
  nullify(reaction%surface_complex_names)
  if (associated(reaction%eqsurfcmplxspecid)) deallocate(reaction%eqsurfcmplxspecid)
  nullify(reaction%eqsurfcmplxspecid)
  if (associated(reaction%eqsurfcmplxstoich)) deallocate(reaction%eqsurfcmplxstoich)
  nullify(reaction%eqsurfcmplxstoich)
  if (associated(reaction%eqsurfcmplxh2oid)) deallocate(reaction%eqsurfcmplxh2oid)
  nullify(reaction%eqsurfcmplxh2oid)
  if (associated(reaction%eqsurfcmplxh2ostoich)) deallocate(reaction%eqsurfcmplxh2ostoich)
  nullify(reaction%eqsurfcmplxh2ostoich)
  if (associated(reaction%eqsurfcmplx_mineral_id)) deallocate(reaction%eqsurfcmplx_mineral_id)
  nullify(reaction%eqsurfcmplx_mineral_id)
  if (associated(reaction%eqsurfcmplx_free_site_id)) deallocate(reaction%eqsurfcmplx_free_site_id)
  nullify(reaction%eqsurfcmplx_free_site_id)
  if (associated(reaction%eqsurfcmplx_free_site_stoich)) deallocate(reaction%eqsurfcmplx_free_site_stoich)
  nullify(reaction%eqsurfcmplx_free_site_stoich)
  if (associated(reaction%eqsurfcmplx_logK)) deallocate(reaction%eqsurfcmplx_logK)
  nullify(reaction%eqsurfcmplx_logK)
  if (associated(reaction%eqsurfcmplx_logKcoef)) deallocate(reaction%eqsurfcmplx_logKcoef)
  nullify(reaction%eqsurfcmplx_logKcoef)
  if (associated(reaction%eqsurfcmplx_Z)) deallocate(reaction%eqsurfcmplx_Z)
  nullify(reaction%eqsurfcmplx_Z)

#if 0  
  if (associated(reaction%kinsurfcmplxspecid)) deallocate(reaction%kinsurfcmplxspecid)
  nullify(reaction%kinsurfcmplxspecid)
  if (associated(reaction%kinsurfcmplxstoich)) deallocate(reaction%kinsurfcmplxstoich)
  nullify(reaction%kinsurfcmplxstoich)
  if (associated(reaction%kinsurfcmplx_freesite_stoich)) deallocate(reaction%kinsurfcmplx_freesite_stoich)
  nullify(reaction%kinsurfcmplx_freesite_stoich)
  if (associated(reaction%kinsurfcmplx_logK)) deallocate(reaction%kinsurfcmplx_logK)
  nullify(reaction%kinsurfcmplx_logK)
  if (associated(reaction%kinsurfcmplx_logKcoef)) deallocate(reaction%kinsurfcmplx_logKcoef)
  nullify(reaction%kinsurfcmplx_logKcoef)
  if (associated(reaction%kinsurfcmplx_Z)) deallocate(reaction%kinsurfcmplx_Z)
  nullify(reaction%kinsurfcmplx_Z)
#endif
  
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
  if (associated(reaction%kinmnrlh2ostoich)) deallocate(reaction%kinmnrlh2ostoich)
  nullify(reaction%kinmnrlh2ostoich)
  if (associated(reaction%kinmnrl_logK)) deallocate(reaction%kinmnrl_logK)
  nullify(reaction%kinmnrl_logK)
  if (associated(reaction%kinmnrl_logKcoef)) deallocate(reaction%kinmnrl_logKcoef)
  nullify(reaction%kinmnrl_logKcoef)
  if (associated(reaction%kinmnrl_rate)) deallocate(reaction%kinmnrl_rate)
  nullify(reaction%kinmnrl_rate)
  if (associated(reaction%kinmnrl_molar_vol)) deallocate(reaction%kinmnrl_molar_vol)
  nullify(reaction%kinmnrl_molar_vol)  
  if (associated(reaction%kinmnrl_num_prefactors)) deallocate(reaction%kinmnrl_num_prefactors)
  nullify(reaction%kinmnrl_num_prefactors)
  if (associated(reaction%kinmnrl_pri_prefactor_id)) deallocate(reaction%kinmnrl_pri_prefactor_id)
  nullify(reaction%kinmnrl_pri_prefactor_id)
  if (associated(reaction%kinmnrl_pri_pref_alpha_stoich)) deallocate(reaction%kinmnrl_pri_pref_alpha_stoich)
  nullify(reaction%kinmnrl_pri_pref_alpha_stoich)
  if (associated(reaction%kinmnrl_pri_pref_beta_stoich)) deallocate(reaction%kinmnrl_pri_pref_beta_stoich)
  nullify(reaction%kinmnrl_pri_pref_beta_stoich)
  if (associated(reaction%kinmnrl_pri_pref_atten_coef)) deallocate(reaction%kinmnrl_pri_pref_atten_coef)
  nullify(reaction%kinmnrl_pri_pref_atten_coef)
  if (associated(reaction%kinmnrl_sec_prefactor_id)) deallocate(reaction%kinmnrl_sec_prefactor_id)
  nullify(reaction%kinmnrl_sec_prefactor_id)
  if (associated(reaction%kinmnrl_sec_pref_alpha_stoich)) deallocate(reaction%kinmnrl_sec_pref_alpha_stoich)
  nullify(reaction%kinmnrl_sec_pref_alpha_stoich)
  if (associated(reaction%kinmnrl_sec_pref_beta_stoich)) deallocate(reaction%kinmnrl_sec_pref_beta_stoich)
  nullify(reaction%kinmnrl_sec_pref_beta_stoich)
  if (associated(reaction%kinmnrl_sec_pref_atten_coef)) deallocate(reaction%kinmnrl_sec_pref_atten_coef)
  nullify(reaction%kinmnrl_sec_pref_atten_coef)
  if (associated(reaction%kinmnrl_Tempkin_const)) deallocate(reaction%kinmnrl_Tempkin_const)
  nullify(reaction%kinmnrl_Tempkin_const)
  if (associated(reaction%kinmnrl_affinity_power)) deallocate(reaction%kinmnrl_affinity_power)
  nullify(reaction%kinmnrl_affinity_power)

  if (associated(reaction%kinmr_rate)) deallocate(reaction%kinmr_rate)
  nullify(reaction%kinmr_rate)
  
  deallocate(reaction)
  nullify(reaction)

end subroutine ReactionDestroy

end module Reaction_Aux_module
