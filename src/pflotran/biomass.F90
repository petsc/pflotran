module Biomass_module

  use Biomass_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  public :: BiomassRead, &
            BiomassProcessConstraint

contains

! ************************************************************************** !
!
! BiomassRead: Reads biomass species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
subroutine BiomassRead(biomass,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(biomass_type) :: biomass
  type(input_type) :: input
  type(option_type) :: option
  
  type(biomass_species_type), pointer :: new_biomass_species, &
                                         prev_biomass_species
           
  nullify(prev_biomass_species)
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
          
    biomass%nbiomass = biomass%nbiomass + 1
          
    new_biomass_species => BiomassSpeciesCreate()
    call InputReadWord(input,option,new_biomass_species%name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,BIOMASS')    
    if (.not.associated(biomass%list)) then
      biomass%list => new_biomass_species
      new_biomass_species%id = 1
    endif
    if (associated(prev_biomass_species)) then
      prev_biomass_species%next => new_biomass_species
      new_biomass_species%id = prev_biomass_species%id + 1
    endif
    prev_biomass_species => new_biomass_species
    nullify(new_biomass_species)
  enddo

end subroutine BiomassRead

! ************************************************************************** !
!
! BiomassProcessConstraint: Initializes constraints based on biomass
!                           species in system
! author: Glenn Hammond
! date: 01/07/13
!
! ************************************************************************** !
subroutine BiomassProcessConstraint(biomass,constraint_name, &
                                    constraint,option)
  use Option_module
  use Input_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(biomass_type), pointer :: biomass
  character(len=MAXWORDLENGTH) :: constraint_name
  type(biomass_constraint_type), pointer :: constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: ibiomass, jbiomass
  
  character(len=MAXWORDLENGTH) :: biomass_name(biomass%nbiomass)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(biomass%nbiomass)
  PetscReal :: constraint_conc(biomass%nbiomass)
  PetscBool :: external_dataset(biomass%nbiomass)
  
  if (.not.associated(constraint)) return
  
  biomass_name = ''
  constraint_aux_string = ''
  external_dataset = PETSC_FALSE
  do ibiomass = 1, biomass%nbiomass
    found = PETSC_FALSE
    do jbiomass = 1, biomass%nbiomass
      if (StringCompare(constraint%names(ibiomass), &
                        biomass%names(jbiomass), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
                'Biomass species "' // trim(constraint%names(ibiomass)) // &
                '" from CONSTRAINT "' // trim(constraint_name) // &
                '" not found among biomass species.'
      call printErrMsg(option)
    else
      biomass_name(ibiomass) = constraint%names(ibiomass)
      constraint_conc(ibiomass) = &
        constraint%constraint_conc(ibiomass)
      constraint_aux_string(ibiomass) = &
        constraint%constraint_aux_string(ibiomass)
      external_dataset(ibiomass) = constraint%external_dataset(ibiomass)
    endif  
  enddo
  constraint%names = biomass_name
  constraint%constraint_conc = constraint_conc
  constraint%constraint_aux_string = constraint_aux_string
  constraint%external_dataset = external_dataset

end subroutine BiomassProcessConstraint


end module Biomass_module
