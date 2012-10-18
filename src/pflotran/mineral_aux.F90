module Mineral_Aux_module
  
  use Database_Aux_module

  implicit none
  
  private 

#include "definitions.h"
  
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

  type, public :: transition_state_rxn_type
    PetscReal :: affinity_factor_sigma
    PetscReal :: affinity_factor_beta
    PetscReal :: affinity_threshold
    PetscReal :: rate_limiter
    PetscReal :: surf_area_vol_frac_pwr
    PetscReal :: surf_area_porosity_pwr
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
  
  type, public :: mineral_rxn_type

    PetscInt :: nmnrl
    character(len=MAXWORDLENGTH), pointer :: mineral_names(:)
    
    type(mineral_type), pointer :: mineral_list

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
    PetscReal, pointer :: kinmnrl_surf_area_vol_frac_pwr(:)
    PetscReal, pointer :: kinmnrl_surf_area_porosity_pwr(:)
    PetscInt, pointer :: kinmnrl_irreversible(:)
   
  end type mineral_rxn_type

  public :: MineralReactionCreate, &
            GetMineralCount, &
            GetMineralNames, &
            GetMineralIDFromName, &
            TransitionStateTheoryRxnCreate, &
            TransitionStatePrefactorCreate, &
            TSPrefactorSpeciesCreate, &
            TransitionStateTheoryRxnDestroy, &
            MineralCreate, &
            MineralDestroy, &
            MineralConstraintCreate, &
            MineralConstraintDestroy, &
            MineralReactionDestroy
             
contains

! ************************************************************************** !
!
! MineralReactionCreate: Allocate and initialize mineral reaction object
! author: Glenn Hammond
! date: 08/16/12
!
! ************************************************************************** !
function MineralReactionCreate()

  implicit none
  
  type(mineral_rxn_type), pointer :: MineralReactionCreate
  
  type(mineral_rxn_type), pointer :: mineral_reaction

  allocate(mineral_reaction)  
    
  nullify(mineral_reaction%mineral_list)
  
  ! for saturation states
  mineral_reaction%nmnrl = 0  
  nullify(mineral_reaction%mineral_names)
  nullify(mineral_reaction%mnrl_print)
  nullify(mineral_reaction%mnrlspecid)
  nullify(mineral_reaction%mnrlh2oid)
  nullify(mineral_reaction%mnrlstoich)
  nullify(mineral_reaction%mnrlh2ostoich)
  nullify(mineral_reaction%mnrl_logK)
  nullify(mineral_reaction%mnrl_logKcoef)
  
  ! for kinetic mineral reactions
  mineral_reaction%nkinmnrl = 0  
  nullify(mineral_reaction%kinmnrl_names)
  nullify(mineral_reaction%kinmnrl_print)
  nullify(mineral_reaction%kinmnrlspecid)
  nullify(mineral_reaction%kinmnrlstoich)
  nullify(mineral_reaction%kinmnrlh2oid)
  nullify(mineral_reaction%kinmnrlh2ostoich)
  nullify(mineral_reaction%kinmnrl_logK)
  nullify(mineral_reaction%kinmnrl_logKcoef)
  nullify(mineral_reaction%kinmnrl_rate)
  nullify(mineral_reaction%kinmnrl_activation_energy)
  nullify(mineral_reaction%kinmnrl_molar_vol)
  nullify(mineral_reaction%kinmnrl_molar_wt)

  nullify(mineral_reaction%kinmnrl_num_prefactors)
  nullify(mineral_reaction%kinmnrl_prefactor_id)
  nullify(mineral_reaction%kinmnrl_pref_alpha)
  nullify(mineral_reaction%kinmnrl_pref_beta)
  nullify(mineral_reaction%kinmnrl_pref_atten_coef)
  nullify(mineral_reaction%kinmnrl_pref_rate)
  nullify(mineral_reaction%kinmnrl_pref_activation_energy)

  nullify(mineral_reaction%kinmnrl_Tempkin_const)
  nullify(mineral_reaction%kinmnrl_affinity_power)
  nullify(mineral_reaction%kinmnrl_affinity_threshold)
  nullify(mineral_reaction%kinmnrl_irreversible)
  nullify(mineral_reaction%kinmnrl_rate_limiter)
  nullify(mineral_reaction%kinmnrl_surf_area_vol_frac_pwr)
  nullify(mineral_reaction%kinmnrl_surf_area_porosity_pwr)

  MineralReactionCreate => mineral_reaction
  
end function MineralReactionCreate

! ************************************************************************** !
!
! MineralCreate: Allocate and initialize a mineral object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function MineralCreate()

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
  tstrxn%affinity_factor_sigma = -999.d0
  tstrxn%affinity_factor_beta = -999.d0
  tstrxn%affinity_threshold = 0.d0
  tstrxn%surf_area_vol_frac_pwr = 0.d0
  tstrxn%surf_area_porosity_pwr = 0.d0
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
! MineralConstraintCreate: Creates a mineral constraint object
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function MineralConstraintCreate(mineral_reaction,option)

  use Option_module
  
  implicit none
  
  type(mineral_rxn_type) :: mineral_reaction
  type(option_type) :: option
  type(mineral_constraint_type), pointer :: MineralConstraintCreate

  type(mineral_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(mineral_reaction%nkinmnrl))
  constraint%names = ''
  allocate(constraint%constraint_vol_frac(mineral_reaction%nkinmnrl))
  constraint%constraint_vol_frac = 0.d0
  allocate(constraint%constraint_area(mineral_reaction%nkinmnrl))
  constraint%constraint_area = 0.d0
  allocate(constraint%basis_vol_frac(mineral_reaction%nkinmnrl))
  constraint%basis_vol_frac = 0.d0
  allocate(constraint%basis_area(mineral_reaction%nkinmnrl))
  constraint%basis_area = 0.d0
  allocate(constraint%constraint_aux_string(mineral_reaction%nkinmnrl))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(mineral_reaction%nkinmnrl))
  constraint%external_dataset = PETSC_FALSE

  MineralConstraintCreate => constraint

end function MineralConstraintCreate
  
! ************************************************************************** !
!
! GetMineralNames: Returns the names of minerals in an array
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
function GetMineralNames(mineral_reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetMineralNames(:)
  type(mineral_rxn_type) :: mineral_reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(mineral_type), pointer :: mineral

  count = GetMineralCount(mineral_reaction)
  allocate(names(count))
  
  count = 1
  mineral => mineral_reaction%mineral_list
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
function GetMineralCount(mineral_reaction)

  implicit none
  
  PetscInt :: GetMineralCount
  type(mineral_rxn_type) :: mineral_reaction

  type(mineral_type), pointer :: mineral

  GetMineralCount = 0
  mineral => mineral_reaction%mineral_list
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
function GetMineralIDFromName(mineral_reaction,name)

  use String_module
  
  implicit none
  
  type(mineral_rxn_type) :: mineral_reaction
  character(len=MAXWORDLENGTH) :: name

  PetscInt :: GetMineralIDFromName
  type(mineral_type), pointer :: mineral

  GetMineralIDFromName = -1
 
  mineral => mineral_reaction%mineral_list
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
  if (associated(constraint%external_dataset)) &
    deallocate(constraint%external_dataset)
  nullify(constraint%external_dataset)

  deallocate(constraint)
  nullify(constraint)

end subroutine MineralConstraintDestroy

! ************************************************************************** !
!
! MineralReactionDestroy: Deallocates a mineral reaction object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine MineralReactionDestroy(mineral_reaction)

  implicit none

  type(mineral_rxn_type), pointer :: mineral_reaction
  
  type(mineral_type), pointer :: mineral, prev_mineral

  if (.not.associated(mineral_reaction)) return
  
  ! mineral species
  mineral => mineral_reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    prev_mineral => mineral
    mineral => mineral%next
    call MineralDestroy(prev_mineral)
  enddo    
  nullify(mineral_reaction%mineral_list)
  
  if (associated(mineral_reaction%mineral_names)) &
    deallocate(mineral_reaction%mineral_names)
  nullify(mineral_reaction%mineral_names)
  if (associated(mineral_reaction%kinmnrl_names)) &
    deallocate(mineral_reaction%kinmnrl_names)
  nullify(mineral_reaction%kinmnrl_names)

  if (associated(mineral_reaction%mnrl_print)) &
    deallocate(mineral_reaction%mnrl_print)
  nullify(mineral_reaction%mnrl_print)
  if (associated(mineral_reaction%kinmnrl_print)) &
    deallocate(mineral_reaction%kinmnrl_print)
  nullify(mineral_reaction%kinmnrl_print)

  if (associated(mineral_reaction%mnrlspecid)) &
    deallocate(mineral_reaction%mnrlspecid)
  nullify(mineral_reaction%mnrlspecid)
  if (associated(mineral_reaction%mnrlstoich)) &
    deallocate(mineral_reaction%mnrlstoich)
  nullify(mineral_reaction%mnrlstoich)
  if (associated(mineral_reaction%mnrlh2oid)) &
    deallocate(mineral_reaction%mnrlh2oid)
  nullify(mineral_reaction%mnrlh2oid)
  if (associated(mineral_reaction%mnrlh2ostoich)) &
    deallocate(mineral_reaction%mnrlh2ostoich)
  nullify(mineral_reaction%mnrlh2ostoich)
  if (associated(mineral_reaction%mnrl_logK)) &
    deallocate(mineral_reaction%mnrl_logK)
  nullify(mineral_reaction%mnrl_logK)
  if (associated(mineral_reaction%mnrl_logKcoef)) &
    deallocate(mineral_reaction%mnrl_logKcoef)
  nullify(mineral_reaction%mnrl_logKcoef)
  
  if (associated(mineral_reaction%kinmnrlspecid)) &
    deallocate(mineral_reaction%kinmnrlspecid)
  nullify(mineral_reaction%kinmnrlspecid)
  if (associated(mineral_reaction%kinmnrlstoich)) &
    deallocate(mineral_reaction%kinmnrlstoich)
  nullify(mineral_reaction%kinmnrlstoich)
  if (associated(mineral_reaction%kinmnrlh2oid)) &
    deallocate(mineral_reaction%kinmnrlh2oid)
  nullify(mineral_reaction%kinmnrlh2oid)
  if (associated(mineral_reaction%kinmnrlh2ostoich)) &
    deallocate(mineral_reaction%kinmnrlh2ostoich)
  nullify(mineral_reaction%kinmnrlh2ostoich)
  if (associated(mineral_reaction%kinmnrl_logK)) &
    deallocate(mineral_reaction%kinmnrl_logK)
  nullify(mineral_reaction%kinmnrl_logK)
  if (associated(mineral_reaction%kinmnrl_logKcoef)) &
    deallocate(mineral_reaction%kinmnrl_logKcoef)
  nullify(mineral_reaction%kinmnrl_logKcoef)
  if (associated(mineral_reaction%kinmnrl_rate)) &
    deallocate(mineral_reaction%kinmnrl_rate)
  nullify(mineral_reaction%kinmnrl_rate)
  if (associated(mineral_reaction%kinmnrl_molar_vol)) &
    deallocate(mineral_reaction%kinmnrl_molar_vol)
  nullify(mineral_reaction%kinmnrl_molar_vol)  
   if (associated(mineral_reaction%kinmnrl_molar_wt)) &
    deallocate(mineral_reaction%kinmnrl_molar_wt)
  nullify(mineral_reaction%kinmnrl_molar_wt)  

  if (associated(mineral_reaction%kinmnrl_num_prefactors)) &
    deallocate(mineral_reaction%kinmnrl_num_prefactors)
  nullify(mineral_reaction%kinmnrl_num_prefactors)
  if (associated(mineral_reaction%kinmnrl_prefactor_id)) &
    deallocate(mineral_reaction%kinmnrl_prefactor_id)
  nullify(mineral_reaction%kinmnrl_prefactor_id)
  if (associated(mineral_reaction%kinmnrl_pref_alpha)) &
    deallocate(mineral_reaction%kinmnrl_pref_alpha)
  nullify(mineral_reaction%kinmnrl_pref_alpha)
  if (associated(mineral_reaction%kinmnrl_pref_beta)) &
    deallocate(mineral_reaction%kinmnrl_pref_beta)
  nullify(mineral_reaction%kinmnrl_pref_beta)
  if (associated(mineral_reaction%kinmnrl_pref_atten_coef)) &
    deallocate(mineral_reaction%kinmnrl_pref_atten_coef)
  nullify(mineral_reaction%kinmnrl_pref_atten_coef)
  if (associated(mineral_reaction%kinmnrl_pref_rate)) &
    deallocate(mineral_reaction%kinmnrl_pref_rate)
  nullify(mineral_reaction%kinmnrl_pref_rate)
  if (associated(mineral_reaction%kinmnrl_pref_activation_energy)) &
    deallocate(mineral_reaction%kinmnrl_pref_activation_energy)
  nullify(mineral_reaction%kinmnrl_pref_activation_energy)

  if (associated(mineral_reaction%kinmnrl_Tempkin_const)) &
    deallocate(mineral_reaction%kinmnrl_Tempkin_const)
  nullify(mineral_reaction%kinmnrl_Tempkin_const)
  if (associated(mineral_reaction%kinmnrl_affinity_power)) &
    deallocate(mineral_reaction%kinmnrl_affinity_power)
  nullify(mineral_reaction%kinmnrl_affinity_power)
  if (associated(mineral_reaction%kinmnrl_affinity_threshold)) &
    deallocate(mineral_reaction%kinmnrl_affinity_threshold)
  nullify(mineral_reaction%kinmnrl_affinity_threshold)
  if (associated(mineral_reaction%kinmnrl_activation_energy)) &
    deallocate(mineral_reaction%kinmnrl_activation_energy)
  nullify(mineral_reaction%kinmnrl_activation_energy)
  if (associated(mineral_reaction%kinmnrl_rate_limiter)) &
    deallocate(mineral_reaction%kinmnrl_rate_limiter)
  nullify(mineral_reaction%kinmnrl_rate_limiter)
  if (associated(mineral_reaction%kinmnrl_surf_area_vol_frac_pwr)) &
    deallocate(mineral_reaction%kinmnrl_surf_area_vol_frac_pwr)
  nullify(mineral_reaction%kinmnrl_surf_area_vol_frac_pwr)
  if (associated(mineral_reaction%kinmnrl_surf_area_porosity_pwr)) &
    deallocate(mineral_reaction%kinmnrl_surf_area_porosity_pwr)
  nullify(mineral_reaction%kinmnrl_surf_area_porosity_pwr)
  if (associated(mineral_reaction%kinmnrl_irreversible)) &
    deallocate(mineral_reaction%kinmnrl_irreversible)
  nullify(mineral_reaction%kinmnrl_irreversible)

  deallocate(mineral_reaction)
  nullify(mineral_reaction)

end subroutine MineralReactionDestroy

end module Mineral_Aux_module
