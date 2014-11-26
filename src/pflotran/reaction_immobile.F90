module Immobile_module

  use Immobile_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  public :: ImmobileRead, &
            ImmobileProcessConstraint

contains

! ************************************************************************** !

subroutine ImmobileRead(immobile,input,option)
  ! 
  ! Reads immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(immobile_type) :: immobile
  type(input_type) :: input
  type(option_type) :: option
  
  type(immobile_species_type), pointer :: new_immobile_species, &
                                         prev_immobile_species
           
  nullify(prev_immobile_species)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
          
    immobile%nimmobile = immobile%nimmobile + 1
          
    new_immobile_species => ImmobileSpeciesCreate()
    call InputReadWord(input,option,new_immobile_species%name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,IMMOBILE_SPECIES')    
    if (.not.associated(immobile%list)) then
      immobile%list => new_immobile_species
      new_immobile_species%id = 1
    endif
    if (associated(prev_immobile_species)) then
      prev_immobile_species%next => new_immobile_species
      new_immobile_species%id = prev_immobile_species%id + 1
    endif
    prev_immobile_species => new_immobile_species
    nullify(new_immobile_species)
  enddo

end subroutine ImmobileRead

! ************************************************************************** !

subroutine ImmobileProcessConstraint(immobile,constraint_name, &
                                    constraint,option)
  ! 
  ! Initializes constraints based on immobile
  ! species in system
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/07/13
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(immobile_type), pointer :: immobile
  character(len=MAXWORDLENGTH) :: constraint_name
  type(immobile_constraint_type), pointer :: constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: iimmobile, jimmobile
  
  character(len=MAXWORDLENGTH) :: immobile_name(immobile%nimmobile)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(immobile%nimmobile)
  PetscReal :: constraint_conc(immobile%nimmobile)
  PetscBool :: external_dataset(immobile%nimmobile)
  
  if (.not.associated(constraint)) return
  
  immobile_name = ''
  constraint_aux_string = ''
  external_dataset = PETSC_FALSE
  do iimmobile = 1, immobile%nimmobile
    found = PETSC_FALSE
    do jimmobile = 1, immobile%nimmobile
      if (StringCompare(constraint%names(iimmobile), &
                        immobile%names(jimmobile), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
                'Immobile species "' // trim(constraint%names(iimmobile)) // &
                '" from CONSTRAINT "' // trim(constraint_name) // &
                '" not found among immobile species.'
      call printErrMsg(option)
    else
      immobile_name(iimmobile) = constraint%names(iimmobile)
      constraint_conc(iimmobile) = &
        constraint%constraint_conc(iimmobile)
      constraint_aux_string(iimmobile) = &
        constraint%constraint_aux_string(iimmobile)
      external_dataset(iimmobile) = constraint%external_dataset(iimmobile)
    endif  
  enddo
  constraint%names = immobile_name
  constraint%constraint_conc = constraint_conc
  constraint%constraint_aux_string = constraint_aux_string
  constraint%external_dataset = external_dataset

end subroutine ImmobileProcessConstraint


end module Immobile_module
