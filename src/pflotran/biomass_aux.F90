module Biomass_Aux_module
  
  use Database_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: biomass_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(biomass_species_type), pointer :: next    
  end type biomass_species_type
  
  type, public :: biomass_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type biomass_constraint_type
  
  type, public :: biomass_type

    PetscInt :: nbiomass
    PetscBool :: print_all
    
    type(biomass_species_type), pointer :: list

    ! biomass species
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscBool, pointer :: print_me(:)    
    PetscInt, pointer :: immobile_id(:)    

  end type biomass_type

  public :: BiomassCreate, &
            BiomassSpeciesCreate, &
            BiomassConstraintCreate, &
            BiomassGetCount, &
            BiomassConstraintDestroy, &
            BiomassDestroy
             
contains

! ************************************************************************** !
!
! BiomassCreate: Allocate and initialize biomass object
! author: Glenn Hammond
! date: 01/11/13
!
! ************************************************************************** !
function BiomassCreate()

  implicit none
  
  type(biomass_type), pointer :: BiomassCreate
  
  type(biomass_type), pointer :: biomass

  allocate(biomass)  
  nullify(biomass%list)
  biomass%nbiomass = 0
  biomass%print_all = PETSC_FALSE
  nullify(biomass%names)
  nullify(biomass%print_me)
  nullify(biomass%immobile_id)

  BiomassCreate => biomass
  
end function BiomassCreate

! ************************************************************************** !
!
! BiomassSpeciesCreate: Allocate and initialize a biomass species object
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function BiomassSpeciesCreate()

  implicit none
  
  type(biomass_species_type), pointer :: BiomassSpeciesCreate
  
  type(biomass_species_type), pointer :: species

  allocate(species)  
  species%id = 0
  species%name = ''
  species%molar_weight = 0.d0
  species%print_me = PETSC_FALSE
  nullify(species%next)

  BiomassSpeciesCreate => species
  
end function BiomassSpeciesCreate

! ************************************************************************** !
!
! BiomassConstraintCreate: Creates a biomass constraint object
! author: Glenn Hammond
! date: 01/07/13
!
! ************************************************************************** !
function BiomassConstraintCreate(biomass,option)

  use Option_module
  
  implicit none
  
  type(biomass_type) :: biomass
  type(option_type) :: option
  type(biomass_constraint_type), pointer :: BiomassConstraintCreate

  type(biomass_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(biomass%nbiomass))
  constraint%names = ''
  allocate(constraint%constraint_conc(biomass%nbiomass))
  constraint%constraint_conc = 0.d0
  allocate(constraint%constraint_aux_string(biomass%nbiomass))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(biomass%nbiomass))
  constraint%external_dataset = PETSC_FALSE

  BiomassConstraintCreate => constraint

end function BiomassConstraintCreate
  
! ************************************************************************** !
!
! BiomassGetCount: Returns the number of biomass species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function BiomassGetCount(biomass)

  implicit none
  
  PetscInt :: BiomassGetCount
  type(biomass_type) :: biomass

  type(biomass_species_type), pointer :: biomass_species

  BiomassGetCount = 0
  biomass_species => biomass%list
  do
    if (.not.associated(biomass_species)) exit
    BiomassGetCount = BiomassGetCount + 1
    biomass_species => biomass_species%next
  enddo

end function BiomassGetCount

! ************************************************************************** !
!
! BiomassSpeciesDestroy: Deallocates a biomass species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
subroutine BiomassSpeciesDestroy(species)

  implicit none
    
  type(biomass_species_type), pointer :: species

  if (.not.associated(species)) return
  
  deallocate(species)  
  nullify(species)

end subroutine BiomassSpeciesDestroy

! ************************************************************************** !
!
! BiomassConstraintDestroy: Destroys a colloid constraint object
! author: Glenn Hammond
! date: 03/12/10
!
! ************************************************************************** !
subroutine BiomassConstraintDestroy(constraint)

  use Utility_module, only: DeallocateArray

  implicit none
  
  type(biomass_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%constraint_aux_string)
  call DeallocateArray(constraint%external_dataset)
  
  deallocate(constraint)
  nullify(constraint)

end subroutine BiomassConstraintDestroy

! ************************************************************************** !
!
! BiomassDestroy: Deallocates a biomass object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine BiomassDestroy(biomass)

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(biomass_type), pointer :: biomass
  
  type(biomass_species_type), pointer :: cur_biomass_species, &
                                         prev_biomass_species

  if (.not.associated(biomass)) return

  ! biomass species
  cur_biomass_species => biomass%list
  do
    if (.not.associated(cur_biomass_species)) exit
    prev_biomass_species => cur_biomass_species
    cur_biomass_species => cur_biomass_species%next
    call BiomassSpeciesDestroy(prev_biomass_species)
  enddo    
  nullify(biomass%list)
  
  call DeallocateArray(biomass%names)
  call DeallocateArray(biomass%print_me)
  call DeallocateArray(biomass%immobile_id)
  
  deallocate(biomass)
  nullify(biomass)

end subroutine BiomassDestroy

end module Biomass_Aux_module
