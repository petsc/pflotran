module Reaction_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"
  
  type, public :: aq_species_type
    PetscInt :: id
    character(len=MAXNAMELENGTH) :: name
    PetscReal :: a0
    PetscReal :: molar_weight
    PetscReal :: Z
    type(equilibrium_rxn_type), pointer :: eqrxn
    type(aq_species_type), pointer :: next
  end type aq_species_type

  type, public :: gas_species_type
    PetscInt :: id
    character(len=MAXNAMELENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    type(equilibrium_rxn_type), pointer :: eqrxn
    type(gas_species_type), pointer :: next    
  end type gas_species_type

  type, public :: equilibrium_rxn_type
    PetscInt :: nspec
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK(:)
  end type equilibrium_rxn_type

  type, public :: kinetic_rxn_type
    PetscInt :: nspec
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK(:)
    PetscReal :: rate_forward
    PetscReal :: rate_reverse
  end type kinetic_rxn_type

  type, public :: mineral_type
    PetscInt :: id
    character(len=MAXNAMELENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    type(transition_state_rxn_type), pointer :: tstrxn
    type(mineral_type), pointer :: next
  end type mineral_type

  type, public :: transition_state_rxn_type
    PetscInt :: nspec
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: stoich(:)
    PetscReal, pointer :: logK(:)
    PetscInt :: nspec_primary_prefactor
    character(len=MAXNAMELENGTH), pointer :: spec_name_primary_prefactor(:)
    PetscInt, pointer :: spec_ids_primary_prefactor(:)
    PetscReal, pointer :: stoich_primary_prefactor(:)
    PetscInt :: nspec_secondary_prefactor
    character(len=MAXNAMELENGTH), pointer :: spec_name_secondary_prefactor(:)
    PetscInt, pointer :: spec_ids_secondary_prefactor(:)
    PetscReal, pointer :: stoich_secondary_prefactor(:)
    PetscReal :: affinity_factor_sigma
    PetscReal :: affinity_factor_beta
    PetscReal :: rate
    PetscReal :: area0
  end type transition_state_rxn_type
  
  type, public :: ion_exchange_rxn_type
    character(len=MAXNAMELENGTH) :: reference_cation_name
    character(len=MAXNAMELENGTH) :: mnrl_name
    PetscInt :: ncation
    character(len=MAXNAMELENGTH), pointer :: cation_name(:)
    PetscInt, pointer :: cation_ids(:)
    PetscReal, pointer :: k(:)  ! selectivity coefficient
    PetscReal :: CEC
    type (ion_exchange_rxn_type), pointer :: next
  end type ion_exchange_rxn_type

  type, public :: surface_complexation_rxn_type
    PetscInt :: id
    character(len=MAXNAMELENGTH) :: name
    PetscInt :: mineral_id
    character(len=MAXNAMELENGTH) :: mineral_name
    PetscInt :: nspec
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: stoich(:)
    PetscReal :: free_site_stoich
    PetscReal, pointer :: logK(:)
    PetscReal :: Z
    type (surface_complexation_rxn_type), pointer :: next
  end type surface_complexation_rxn_type

  type, public :: aq_species_constraint_type
    character(len=MAXNAMELENGTH), pointer :: names(:)
    PetscReal, pointer :: conc(:)
    PetscReal, pointer :: basis_conc(:)
    PetscInt, pointer :: constraint_type(:)
    character(len=MAXNAMELENGTH), pointer :: constraint_spec_name(:)
  end type aq_species_constraint_type

  type, public :: mineral_constraint_type
    character(len=MAXNAMELENGTH), pointer :: names(:)
    PetscReal, pointer :: conc(:)
    PetscReal, pointer :: basis_conc(:)
  end type mineral_constraint_type

  type, public :: reaction_type
    character(len=MAXSTRINGLENGTH) :: database_filename
    PetscInt :: num_dbase_temperatures
    PetscReal, pointer :: dbase_temperatures(:)
    type(aq_species_type), pointer :: primary_species_list
    type(aq_species_type), pointer :: secondary_species_list
    type(gas_species_type), pointer :: gas_species_list
    type(mineral_type), pointer :: mineral_list
    type(ion_exchange_rxn_type), pointer :: ion_exchange_list
    type(surface_complexation_rxn_type), pointer :: surface_complex_list
    ! compressed arrays for efficient computation
    ! primary aqueous complexes
    PetscInt :: ncomp
    character(len=MAXNAMELENGTH), pointer :: primary_species_names(:)
    PetscReal, pointer :: primary_spec_a0(:)
    PetscReal, pointer :: primary_spec_Z(:)
    ! aqueous complexes
    PetscInt :: neqcmplx
    character(len=MAXNAMELENGTH), pointer :: secondary_species_names(:)
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
    ! ionx exchange reactions
    character(len=MAXNAMELENGTH), pointer :: ion_exchange_names(:)
    PetscInt, pointer :: eqionx_ncation(:)
    PetscReal, pointer :: eqionx_CEC(:)
    PetscReal, pointer :: eqionx_k(:,:)
    PetscInt, pointer :: eqionx_cationid(:)
    PetscInt, pointer :: eqionx_rxn_offset(:)
    PetscInt, pointer :: kinionx_ncation(:)
    PetscReal, pointer :: kinionx_CEC(:)
    PetscReal, pointer :: kinionx_k(:,:)
    PetscInt, pointer :: kinionx_cationid(:)
    PetscInt, pointer :: kinionx_rxn_offset(:)
    ! surface complexation reactions
    character(len=MAXNAMELENGTH), pointer :: surface_complex_names(:)
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
    ! mineral reactions
    PetscInt :: nmnrl
    character(len=MAXNAMELENGTH), pointer :: mineral_names(:)
      ! for saturation states
    PetscInt, pointer :: mnrlspecid(:,:)
    PetscReal, pointer :: mnrlstoich(:,:)
    PetscReal, pointer :: mnrl_logK(:)
    PetscReal, pointer :: mnrl_molar_vol(:)
      ! for kinetic reactions
    PetscInt :: nkinmnrl
    character(len=MAXNAMELENGTH), pointer :: kinmnrl_names(:)
    PetscInt, pointer :: kinmnrlspecid(:,:)
    PetscReal, pointer :: kinmnrlstoich(:,:)
    PetscInt, pointer :: kinmnrlh2oid(:)
    PetscReal, pointer :: kinmnrlh2ostoich(:)
    PetscReal, pointer :: kinmnrl_logK(:)
    PetscReal, pointer :: kinmnrl_logKcoef(:,:)
    PetscReal, pointer :: kinmnrl_rate(:,:)
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
  end type reaction_type

  public :: ReactionCreate, &
            AqueousSpeciesCreate, &
            GasSpeciesCreate, &
            MineralCreate, &
            GetPrimarySpeciesCount, &
            GetPrimarySpeciesNames, &
            GetSecondarySpeciesCount, &
            GetSecondarySpeciesNames, &
            GetMineralCount, &
            GetMineralNames, &
            EquilibriumRxnCreate, &
            EquilibriumRxnDestroy, &
            TransitionStateTheoryRxnCreate, &
            TransitionStateTheoryRxnDestroy, &
            SurfaceComplexationRxnCreate, &
            AqueousSpeciesConstraintDestroy, &
            MineralConstraintDestroy, &
            AqueousSpeciesConstraintCreate, &
            MineralConstraintCreate
             
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

  nullify(reaction%primary_species_list)
  nullify(reaction%secondary_species_list)
  nullify(reaction%gas_species_list)
  nullify(reaction%mineral_list)
  nullify(reaction%ion_exchange_list)
  nullify(reaction%surface_complex_list)
  
  nullify(reaction%primary_species_names)
  nullify(reaction%secondary_species_names)
  nullify(reaction%surface_complex_names)
  nullify(reaction%ion_exchange_names)
  nullify(reaction%mineral_names)
  nullify(reaction%kinmnrl_names)
  
  reaction%ncomp = 0
  nullify(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_Z)
  
  reaction%neqcmplx = 0
  nullify(reaction%eqcmplxspecid)
  nullify(reaction%eqcmplxstoich)
  nullify(reaction%eqcmplxh2oid)
  nullify(reaction%eqcmplxh2ostoich)
  nullify(reaction%eqcmplx_a0)
  nullify(reaction%eqcmplx_Z)
  nullify(reaction%eqcmplx_logK)
  nullify(reaction%eqcmplx_logKcoef)
  
  reaction%debyeA = 0.d0
  reaction%debyeB = 0.d0
  reaction%debyeBdot = 0.d0
  
  nullify(reaction%eqionx_ncation)
  nullify(reaction%eqionx_CEC)
  nullify(reaction%eqionx_k)
  nullify(reaction%eqionx_cationid)
  nullify(reaction%eqionx_rxn_offset)
  
  nullify(reaction%kinionx_ncation)
  nullify(reaction%kinionx_CEC)
  nullify(reaction%kinionx_k)
  nullify(reaction%kinionx_cationid)
  nullify(reaction%kinionx_rxn_offset)
  
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

  reaction%nmnrl = 0  
  nullify(reaction%mnrlspecid)
  nullify(reaction%mnrlstoich)
  nullify(reaction%mnrl_logK)
  nullify(reaction%mnrl_molar_vol)
  
  reaction%nkinmnrl = 0  
  nullify(reaction%kinmnrlspecid)
  nullify(reaction%kinmnrlstoich)
  nullify(reaction%kinmnrlh2oid)
  nullify(reaction%kinmnrlh2ostoich)
  nullify(reaction%kinmnrl_logK)
  nullify(reaction%kinmnrl_logKcoef)
  nullify(reaction%kinmnrl_rate)
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
  mineral%name = ''
  mineral%molar_volume = 0.d0
  mineral%molar_weight = 0.d0
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
  tstrxn%area0 = 0.d0
  
  TransitionStateTheoryRxnCreate => tstrxn
  
end function TransitionStateTheoryRxnCreate

! ************************************************************************** !
!
! SurfaceComplexationRxnCreate: Allocate and initialize a surface complexation
!                               reaction
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
function SurfaceComplexationRxnCreate()

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: SurfaceComplexationRxnCreate

  type(surface_complexation_rxn_type), pointer :: surfcplxrxn
  
  allocate(surfcplxrxn)
  surfcplxrxn%id = 0
  surfcplxrxn%name = ''
  surfcplxrxn%mineral_id = 0
  surfcplxrxn%mineral_name = ''
  surfcplxrxn%free_site_stoich = 0
  surfcplxrxn%nspec = 0
  surfcplxrxn%Z = 0.d0
  nullify(surfcplxrxn%spec_name)
  nullify(surfcplxrxn%stoich)
  nullify(surfcplxrxn%spec_ids)
  nullify(surfcplxrxn%logK)
  
  SurfaceComplexationRxnCreate => surfcplxrxn
  
end function SurfaceComplexationRxnCreate

! ************************************************************************** !
!
! AqueousSpeciesConstraintCreate: Creates an aqueous species constraint 
!                                 object
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function AqueousSpeciesConstraintCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(aq_species_constraint_type), pointer :: AqueousSpeciesConstraintCreate

  type(aq_species_constraint_type), pointer :: constraint
  
  allocate(constraint)
  allocate(constraint%names(option%ncomp))
  constraint%names = ''
  allocate(constraint%conc(option%ncomp))
  constraint%conc = 0.d0
  allocate(constraint%basis_conc(option%ncomp))
  constraint%basis_conc = 0.d0
  allocate(constraint%constraint_type(option%ncomp))
  constraint%constraint_type = 0
  allocate(constraint%constraint_spec_name(option%ncomp))
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
function MineralConstraintCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(mineral_constraint_type), pointer :: MineralConstraintCreate

  type(mineral_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint)
  allocate(constraint%names(option%nmnrl))
  constraint%names = ''
  allocate(constraint%conc(option%nmnrl))
  constraint%conc = 0.d0
  allocate(constraint%basis_conc(option%nmnrl))
  constraint%basis_conc = 0.d0

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
! GetMineralCount: Returns the number of primary species
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
! IonExchangeRxnDestroy: Deallocates an ion exchange reaction
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine IonExchangeRxnDestroy(ionxrxn)

  implicit none
    
  type(ion_exchange_rxn_type), pointer :: ionxrxn

  if (.not.associated(ionxrxn)) return
  
  if (associated(ionxrxn%cation_name)) deallocate(ionxrxn%cation_name)
  nullify(ionxrxn%cation_name)
  if (associated(ionxrxn%cation_ids)) deallocate(ionxrxn%cation_ids)
  nullify(ionxrxn%cation_ids)
  if (associated(ionxrxn%k)) deallocate(ionxrxn%k)
  nullify(ionxrxn%k)

  deallocate(ionxrxn)  
  nullify(ionxrxn)

end subroutine IonExchangeRxnDestroy

! ************************************************************************** !
!
! SurfaceComplexationRxnDestroy: Deallocates a surface complexation reaction
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine SurfaceComplexationRxnDestroy(surfcplxrxn)

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: surfcplxrxn

  if (.not.associated(surfcplxrxn)) return
  
  if (associated(surfcplxrxn%spec_name)) deallocate(surfcplxrxn%spec_name)
  nullify(surfcplxrxn%spec_name)
  if (associated(surfcplxrxn%spec_ids)) deallocate(surfcplxrxn%spec_ids)
  nullify(surfcplxrxn%spec_ids)
  if (associated(surfcplxrxn%stoich)) deallocate(surfcplxrxn%stoich)
  nullify(surfcplxrxn%stoich)
  if (associated(surfcplxrxn%logK)) deallocate(surfcplxrxn%logK)
  nullify(surfcplxrxn%logK)

  deallocate(surfcplxrxn)  
  nullify(surfcplxrxn)

end subroutine SurfaceComplexationRxnDestroy

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
  if (associated(constraint%conc)) &
    deallocate(constraint%conc)
  nullify(constraint%conc)
  if (associated(constraint%basis_conc)) &
    deallocate(constraint%basis_conc)
  nullify(constraint%basis_conc)
  if (associated(constraint%constraint_type)) &
    deallocate(constraint%constraint_type)
  nullify(constraint%constraint_type)
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
  if (associated(constraint%conc)) &
    deallocate(constraint%conc)
  nullify(constraint%conc)
  if (associated(constraint%basis_conc)) &
    deallocate(constraint%basis_conc)
  nullify(constraint%basis_conc)

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
  ionxrxn => reaction%ion_exchange_list
  do
    if (.not.associated(ionxrxn)) exit
    prev_ionxrxn => ionxrxn
    ionxrxn => ionxrxn%next
    call IonExchangeRxnDestroy(prev_ionxrxn)
  enddo    
  nullify(reaction%ion_exchange_list)

  ! surface complexation reactions
  surfcplxrxn => reaction%surface_complex_list
  do
    if (.not.associated(surfcplxrxn)) exit
    prev_surfcplxrxn => surfcplxrxn
    surfcplxrxn => surfcplxrxn%next
    call SurfaceComplexationRxnDestroy(prev_surfcplxrxn)
  enddo    
  nullify(reaction%surface_complex_list)
  
  if (associated(reaction%primary_species_names)) deallocate(reaction%primary_species_names)
  nullify(reaction%primary_species_names)
  if (associated(reaction%secondary_species_names)) deallocate(reaction%secondary_species_names)
  nullify(reaction%secondary_species_names)
  if (associated(reaction%ion_exchange_names)) deallocate(reaction%ion_exchange_names)
  nullify(reaction%ion_exchange_names)
  if (associated(reaction%surface_complex_names)) deallocate(reaction%surface_complex_names)
  nullify(reaction%surface_complex_names)
  if (associated(reaction%mineral_names)) deallocate(reaction%mineral_names)
  nullify(reaction%mineral_names)
  if (associated(reaction%kinmnrl_names)) deallocate(reaction%kinmnrl_names)
  nullify(reaction%kinmnrl_names)
  
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
  
  if (associated(reaction%eqionx_ncation)) deallocate(reaction%eqionx_ncation)
  nullify(reaction%eqionx_ncation)
  if (associated(reaction%eqionx_CEC)) deallocate(reaction%eqionx_CEC)
  nullify(reaction%eqionx_CEC)
  if (associated(reaction%eqionx_k)) deallocate(reaction%eqionx_k)
  nullify(reaction%eqionx_k)
  if (associated(reaction%eqionx_cationid)) deallocate(reaction%eqionx_cationid)
  nullify(reaction%eqionx_cationid)
  if (associated(reaction%eqionx_cationid)) deallocate(reaction%eqionx_cationid)
  nullify(reaction%eqionx_rxn_offset)
  
  if (associated(reaction%kinionx_ncation)) deallocate(reaction%kinionx_ncation)
  nullify(reaction%kinionx_ncation)
  if (associated(reaction%kinionx_CEC)) deallocate(reaction%kinionx_CEC)
  nullify(reaction%kinionx_CEC)
  if (associated(reaction%kinionx_k)) deallocate(reaction%kinionx_k)
  nullify(reaction%kinionx_k)
  if (associated(reaction%kinionx_cationid)) deallocate(reaction%kinionx_cationid)
  nullify(reaction%kinionx_cationid)
  if (associated(reaction%kinionx_rxn_offset)) deallocate(reaction%kinionx_rxn_offset)
  nullify(reaction%kinionx_rxn_offset)
  
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
  if (associated(reaction%mnrl_logK)) deallocate(reaction%mnrl_logK)
  nullify(reaction%mnrl_logK)
  if (associated(reaction%mnrl_molar_vol)) deallocate(reaction%mnrl_molar_vol)
  nullify(reaction%mnrl_molar_vol)
  
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

end subroutine ReactionDestroy

end module Reaction_Aux_module
