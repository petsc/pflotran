module Immobile_Aux_module
  
  use Database_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: immobile_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(immobile_species_type), pointer :: next    
  end type immobile_species_type
  
  type, public :: immobile_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type immobile_constraint_type
  
  type, public :: immobile_type

    PetscInt :: nimmobile
    PetscBool :: print_all
    
    type(immobile_species_type), pointer :: list

    ! immobile species
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscBool, pointer :: print_me(:)    

  end type immobile_type

  public :: ImmobileCreate, &
            ImmobileSpeciesCreate, &
            ImmobileConstraintCreate, &
            ImmobileGetCount, &
            ImmobileConstraintDestroy, &
            GetImmobileSpeciesIDFromName, &
            ImmobileDestroy
             
contains

! ************************************************************************** !
!
! ImmobileCreate: Allocate and initialize immobile object
! author: Glenn Hammond
! date: 01/11/13
!
! ************************************************************************** !
function ImmobileCreate()

  implicit none
  
  type(immobile_type), pointer :: ImmobileCreate
  
  type(immobile_type), pointer :: immobile

  allocate(immobile)  
  nullify(immobile%list)
  immobile%nimmobile = 0
  immobile%print_all = PETSC_FALSE
  nullify(immobile%names)
  nullify(immobile%print_me)

  ImmobileCreate => immobile
  
end function ImmobileCreate

! ************************************************************************** !
!
! ImmobileSpeciesCreate: Allocate and initialize a immobile species object
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function ImmobileSpeciesCreate()

  implicit none
  
  type(immobile_species_type), pointer :: ImmobileSpeciesCreate
  
  type(immobile_species_type), pointer :: species

  allocate(species)  
  species%id = 0
  species%name = ''
  species%molar_weight = 0.d0
  species%print_me = PETSC_FALSE
  nullify(species%next)

  ImmobileSpeciesCreate => species
  
end function ImmobileSpeciesCreate

! ************************************************************************** !
!
! ImmobileConstraintCreate: Creates a immobile constraint object
! author: Glenn Hammond
! date: 01/07/13
!
! ************************************************************************** !
function ImmobileConstraintCreate(immobile,option)

  use Option_module
  
  implicit none
  
  type(immobile_type) :: immobile
  type(option_type) :: option
  type(immobile_constraint_type), pointer :: ImmobileConstraintCreate

  type(immobile_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(immobile%nimmobile))
  constraint%names = ''
  allocate(constraint%constraint_conc(immobile%nimmobile))
  constraint%constraint_conc = 0.d0
  allocate(constraint%constraint_aux_string(immobile%nimmobile))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(immobile%nimmobile))
  constraint%external_dataset = PETSC_FALSE

  ImmobileConstraintCreate => constraint

end function ImmobileConstraintCreate
  
! ************************************************************************** !
!
! ImmobileGetCount: Returns the number of immobile species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function ImmobileGetCount(immobile)

  implicit none
  
  PetscInt :: ImmobileGetCount
  type(immobile_type) :: immobile

  type(immobile_species_type), pointer :: immobile_species

  ImmobileGetCount = 0
  immobile_species => immobile%list
  do
    if (.not.associated(immobile_species)) exit
    ImmobileGetCount = ImmobileGetCount + 1
    immobile_species => immobile_species%next
  enddo

end function ImmobileGetCount

! ************************************************************************** !
!
! GetImmobileSpeciesIDFromName: Returns the id of named immobile species
! author: Glenn Hammond
! date: 01/28/13
!
! ************************************************************************** !
function GetImmobileSpeciesIDFromName(name,immobile,option)

  use Option_module
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(immobile_type) :: immobile
  type(option_type) :: option

  PetscInt :: GetImmobileSpeciesIDFromName

  type(immobile_species_type), pointer :: species
  PetscInt :: i

  GetImmobileSpeciesIDFromName = -999
  
  ! if the primary species name list exists
  if (associated(immobile%names)) then
    do i = 1, size(immobile%names)
      if (StringCompare(name,immobile%names(i), &
                        MAXWORDLENGTH)) then
        GetImmobileSpeciesIDFromName = i
        exit
      endif
    enddo
  else
    species => immobile%list
    i = 0
    do
      if (.not.associated(species)) exit
      i = i + 1
      if (StringCompare(name,species%name,MAXWORDLENGTH)) then
        GetImmobileSpeciesIDFromName = i
        exit
      endif
      species => species%next
    enddo
  endif

  if (GetImmobileSpeciesIDFromName <= 0) then
    option%io_buffer = 'Species "' // trim(name) // &
      '" not founds among immobile species in GetImmobileSpeciesIDFromName().'
    call printErrMsg(option)
  endif
  
end function GetImmobileSpeciesIDFromName

! ************************************************************************** !
!
! ImmobileSpeciesDestroy: Deallocates a immobile species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
subroutine ImmobileSpeciesDestroy(species)

  implicit none
    
  type(immobile_species_type), pointer :: species

  if (.not.associated(species)) return
  
  deallocate(species)  
  nullify(species)

end subroutine ImmobileSpeciesDestroy

! ************************************************************************** !
!
! ImmobileConstraintDestroy: Destroys a colloid constraint object
! author: Glenn Hammond
! date: 03/12/10
!
! ************************************************************************** !
subroutine ImmobileConstraintDestroy(constraint)

  use Utility_module, only: DeallocateArray

  implicit none
  
  type(immobile_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%constraint_aux_string)
  call DeallocateArray(constraint%external_dataset)
  
  deallocate(constraint)
  nullify(constraint)

end subroutine ImmobileConstraintDestroy

! ************************************************************************** !
!
! ImmobileDestroy: Deallocates a immobile object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine ImmobileDestroy(immobile)

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(immobile_type), pointer :: immobile
  
  type(immobile_species_type), pointer :: cur_immobile_species, &
                                         prev_immobile_species

  if (.not.associated(immobile)) return

  ! immobile species
  cur_immobile_species => immobile%list
  do
    if (.not.associated(cur_immobile_species)) exit
    prev_immobile_species => cur_immobile_species
    cur_immobile_species => cur_immobile_species%next
    call ImmobileSpeciesDestroy(prev_immobile_species)
  enddo    
  nullify(immobile%list)
  
  call DeallocateArray(immobile%names)
  call DeallocateArray(immobile%print_me)
  
  deallocate(immobile)
  nullify(immobile)

end subroutine ImmobileDestroy

end module Immobile_Aux_module
