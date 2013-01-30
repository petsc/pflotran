module Mineral_Aux_module
  
  use Database_Aux_module

  implicit none
  
  private 

#include "definitions.h"
  
  type, public :: mineral_rxn_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(transition_state_rxn_type), pointer :: tstrxn
    type(mineral_rxn_type), pointer :: next
  end type mineral_rxn_type

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
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type mineral_constraint_type
  
  type, public :: mineral_type
  
    PetscInt :: nmnrl
    PetscBool :: print_all
    character(len=MAXWORDLENGTH), pointer :: mineral_names(:)
    
    type(mineral_rxn_type), pointer :: mineral_list

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
   
  end type mineral_type

  public :: MineralCreate, &
            GetMineralCount, &
            GetMineralNames, &
            GetMineralIDFromName, &
            TransitionStateTheoryRxnCreate, &
            TransitionStatePrefactorCreate, &
            TSPrefactorSpeciesCreate, &
            TransitionStateTheoryRxnDestroy, &
            MineralRxnCreate, &
            MineralRxnDestroy, &
            MineralConstraintCreate, &
            MineralConstraintDestroy, &
            MineralDestroy
             
contains

! ************************************************************************** !
!
! MineralCreate: Allocate and initialize mineral reaction object
! author: Glenn Hammond
! date: 08/16/12
!
! ************************************************************************** !
function MineralCreate()

  implicit none
  
  type(mineral_type), pointer :: MineralCreate
  
  type(mineral_type), pointer :: mineral

  allocate(mineral)  
    
  nullify(mineral%mineral_list)
  mineral%print_all = PETSC_FALSE
  
  ! for saturation states
  mineral%nmnrl = 0  
  nullify(mineral%mineral_names)
  nullify(mineral%mnrl_print)
  nullify(mineral%mnrlspecid)
  nullify(mineral%mnrlh2oid)
  nullify(mineral%mnrlstoich)
  nullify(mineral%mnrlh2ostoich)
  nullify(mineral%mnrl_logK)
  nullify(mineral%mnrl_logKcoef)
  
  ! for kinetic mineral reactions
  mineral%nkinmnrl = 0  
  nullify(mineral%kinmnrl_names)
  nullify(mineral%kinmnrl_print)
  nullify(mineral%kinmnrlspecid)
  nullify(mineral%kinmnrlstoich)
  nullify(mineral%kinmnrlh2oid)
  nullify(mineral%kinmnrlh2ostoich)
  nullify(mineral%kinmnrl_logK)
  nullify(mineral%kinmnrl_logKcoef)
  nullify(mineral%kinmnrl_rate)
  nullify(mineral%kinmnrl_activation_energy)
  nullify(mineral%kinmnrl_molar_vol)
  nullify(mineral%kinmnrl_molar_wt)

  nullify(mineral%kinmnrl_num_prefactors)
  nullify(mineral%kinmnrl_prefactor_id)
  nullify(mineral%kinmnrl_pref_alpha)
  nullify(mineral%kinmnrl_pref_beta)
  nullify(mineral%kinmnrl_pref_atten_coef)
  nullify(mineral%kinmnrl_pref_rate)
  nullify(mineral%kinmnrl_pref_activation_energy)

  nullify(mineral%kinmnrl_Tempkin_const)
  nullify(mineral%kinmnrl_affinity_power)
  nullify(mineral%kinmnrl_affinity_threshold)
  nullify(mineral%kinmnrl_irreversible)
  nullify(mineral%kinmnrl_rate_limiter)
  nullify(mineral%kinmnrl_surf_area_vol_frac_pwr)
  nullify(mineral%kinmnrl_surf_area_porosity_pwr)

  MineralCreate => mineral
  
end function MineralCreate

! ************************************************************************** !
!
! MineralRxnCreate: Allocate and initialize a mineral object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function MineralRxnCreate()

  implicit none
  
  type(mineral_rxn_type), pointer :: MineralRxnCreate
  
  type(mineral_rxn_type), pointer :: mineral

  allocate(mineral)  
  mineral%id = 0
  mineral%itype = 0
  mineral%name = ''
  mineral%molar_volume = 0.d0
  mineral%molar_weight = 0.d0
  mineral%print_me = PETSC_FALSE
  nullify(mineral%tstrxn)
  nullify(mineral%next)
  
  MineralRxnCreate => mineral
  
end function MineralRxnCreate

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
function MineralConstraintCreate(mineral,option)

  use Option_module
  
  implicit none
  
  type(mineral_type) :: mineral
  type(option_type) :: option
  type(mineral_constraint_type), pointer :: MineralConstraintCreate

  type(mineral_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(mineral%nkinmnrl))
  constraint%names = ''
  allocate(constraint%constraint_vol_frac(mineral%nkinmnrl))
  constraint%constraint_vol_frac = 0.d0
  allocate(constraint%constraint_area(mineral%nkinmnrl))
  constraint%constraint_area = 0.d0
  allocate(constraint%constraint_aux_string(mineral%nkinmnrl))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(mineral%nkinmnrl))
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
function GetMineralNames(mineral)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetMineralNames(:)
  type(mineral_type) :: mineral

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(mineral_rxn_type), pointer :: cur_mineral

  count = GetMineralCount(mineral)
  allocate(names(count))
  
  count = 1
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    names(count) = cur_mineral%name
    count = count + 1
    cur_mineral => cur_mineral%next
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
function GetMineralCount(mineral)

  implicit none
  
  PetscInt :: GetMineralCount
  type(mineral_type) :: mineral

  type(mineral_rxn_type), pointer :: cur_mineral

  GetMineralCount = 0
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    GetMineralCount = GetMineralCount + 1
    cur_mineral => cur_mineral%next
  enddo

end function GetMineralCount

! ************************************************************************** !
!
! GetMineralIDFromName: Returns the id of mineral with the corresponding name
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
function GetMineralIDFromName(mineral,name)

  use String_module
  
  implicit none
  
  type(mineral_type) :: mineral
  character(len=MAXWORDLENGTH) :: name

  PetscInt :: GetMineralIDFromName
  type(mineral_rxn_type), pointer :: cur_mineral

  GetMineralIDFromName = -1
 
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    if (StringCompare(name,cur_mineral%name,MAXWORDLENGTH)) then
      GetMineralIDFromName = cur_mineral%id
      exit
    endif
    cur_mineral => cur_mineral%next
  enddo

end function GetMineralIDFromName

! ************************************************************************** !
!
! MineralDestroy: Deallocates a mineral rxn object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine MineralRxnDestroy(mineral)

  implicit none
    
  type(mineral_rxn_type), pointer :: mineral

  if (associated(mineral%dbaserxn)) &
    call DatabaseRxnDestroy(mineral%dbaserxn)
  if (associated(mineral%tstrxn)) &
    call TransitionStateTheoryRxnDestroy(mineral%tstrxn)

  deallocate(mineral)  
  nullify(mineral)

end subroutine MineralRxnDestroy

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

  use Utility_module, only: DeallocateArray
  
  implicit none
  
  type(mineral_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_vol_frac)
  call DeallocateArray(constraint%constraint_area)
  call DeallocateArray(constraint%constraint_aux_string)
  call DeallocateArray(constraint%external_dataset)

  deallocate(constraint)
  nullify(constraint)

end subroutine MineralConstraintDestroy

! ************************************************************************** !
!
! MineralDestroy: Deallocates a mineral object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine MineralDestroy(mineral)

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(mineral_type), pointer :: mineral
  
  type(mineral_rxn_type), pointer :: cur_mineral, prev_mineral

  if (.not.associated(mineral)) return
  
  ! mineral species
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    prev_mineral => cur_mineral
    cur_mineral => cur_mineral%next
    call MineralRxnDestroy(prev_mineral)
  enddo    
  nullify(mineral%mineral_list)
  
  call DeallocateArray(mineral%mineral_names)
  call DeallocateArray(mineral%kinmnrl_names)
  call DeallocateArray(mineral%mnrl_print)
  call DeallocateArray(mineral%mnrlspecid)
  call DeallocateArray(mineral%mnrlstoich)
  call DeallocateArray(mineral%mnrlh2oid)
  call DeallocateArray(mineral%mnrlh2ostoich)
  call DeallocateArray(mineral%mnrl_logK)
  call DeallocateArray(mineral%mnrl_logKcoef)
  
  call DeallocateArray(mineral%kinmnrlspecid)
  call DeallocateArray(mineral%kinmnrlstoich)
  call DeallocateArray(mineral%kinmnrlh2oid)
  call DeallocateArray(mineral%kinmnrlh2ostoich)
  call DeallocateArray(mineral%kinmnrl_logK)
  call DeallocateArray(mineral%kinmnrl_logKcoef)
  call DeallocateArray(mineral%kinmnrl_rate)
  call DeallocateArray(mineral%kinmnrl_molar_vol)
  call DeallocateArray(mineral%kinmnrl_molar_wt)

  call DeallocateArray(mineral%kinmnrl_num_prefactors)
  call DeallocateArray(mineral%kinmnrl_prefactor_id)
  call DeallocateArray(mineral%kinmnrl_pref_alpha)
  call DeallocateArray(mineral%kinmnrl_pref_beta)
  call DeallocateArray(mineral%kinmnrl_pref_atten_coef)
  call DeallocateArray(mineral%kinmnrl_pref_rate)
  call DeallocateArray(mineral%kinmnrl_pref_activation_energy)
  
  call DeallocateArray(mineral%kinmnrl_Tempkin_const)
  call DeallocateArray(mineral%kinmnrl_affinity_power)
  call DeallocateArray(mineral%kinmnrl_affinity_threshold)
  call DeallocateArray(mineral%kinmnrl_activation_energy)
  call DeallocateArray(mineral%kinmnrl_rate_limiter)
  call DeallocateArray(mineral%kinmnrl_surf_area_vol_frac_pwr)
  call DeallocateArray(mineral%kinmnrl_surf_area_porosity_pwr)
  call DeallocateArray(mineral%kinmnrl_irreversible)
  
  deallocate(mineral)
  nullify(mineral)

end subroutine MineralDestroy

end module Mineral_Aux_module
