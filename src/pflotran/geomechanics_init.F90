#ifdef GEOMECH

module Geomechanics_Init_module


  implicit none

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"
#include "finclude/petscts.h"


  public :: GeomechicsInitReadRequiredCards, &
            GeomechanicsInitReadInput

contains

! ************************************************************************** !
!
! GeomechicsInitReadRequiredCards: Reads the required input file cards
! related to geomechanics
! author: Satish Karra
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechicsInitReadRequiredCards(geomech_realization)

  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Geomech_Grid_module
  use Input_module
  use String_module
  use Patch_module
  use Option_module
  
  implicit none
  
  type(geomech_realization_type)             :: geomech_realization
  type(geomech_discretization_type), pointer :: discretization
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(geomech_patch_type), pointer   :: patch
  type(option_type), pointer          :: option
  type(input_type), pointer           :: input
  
  patch          => geomech_realization%geomech_patch
  option         => geomech_realization%option
  discretization => geomech_realization%discretization
  
  input => geomech_realization%input
  
! Read in select required cards
!.........................................................................

  ! GEOMECHICS information
  string = "GEOMECHANICS"
  call InputFindStringInFile(input,option,string)
  if(InputError(input)) return
  option%ngeomechdof = 3  ! displacements in x, y, z directions
  
  call GeomechanicsInit(geomech_realization,input,option)  

  discretization%itype = UNSTRUCTURED_GRID ! Assuming only unstructured for now
  select case(discretization%itype)
    case(UNSTRUCTURED_GRID)
      patch => GeomechanicsPatchCreate()
      patch%geomech_grid => discretization%grid
      geomech_realization%geomech_patch => patch
  end select


end subroutine GeomechicsInitReadRequiredCards

! ************************************************************************** !
!
! GeomechanicsInit: Reads the required geomechanics data from input file
! author: Satish Karra
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsInit(geomech_realization,input,option)

  use Option_module
  use Input_module
  use String_module
  use Geomech_Grid_module
  use Geomech_Grid_Aux_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  
  implicit none
  
  type(geomech_realization_type)             :: geomech_realization
  type(geomech_discretization_type), pointer :: discretization
  type(input_type)                           :: input
  type(option_type)                          :: option
  type(geomech_grid_type)                    :: geomech_grid
  character(len=MAXWORDLENGTH)               :: word
  
  discretization => geomech_realization%discretization
  
  discretization%grid => GMGridCreate()
  
  if (option%myrank == option%io_rank) then
    write(*,*),'GEOMECHANICS: Subsurface unstructured grid will be used for '// &
               'geomechanics.'
  endif
     
end subroutine GeomechanicsInit

! ************************************************************************** !
!
! GeomechanicsInitReadInput: Reads the geomechanics input data 
! author: Satish Karra
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsInitReadInput(geomech_realization,input,option)

  use Option_module
  use Input_module
  use String_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomech_Grid_module
  use Geomech_Grid_Aux_module
  use Dataset_Aux_module
  use Waypoint_module
  use Geomechanics_Material_module
  ! Still need to add other geomech modules for output, etc once created
  
  implicit none
  
  type(geomech_realization_type)               :: geomech_realization
  type(input_type)                             :: input
  type(option_type)                            :: option
  
  type(geomech_discretization_type), pointer   :: discretization
  type(dataset_type), pointer                  :: dataset
  type(geomech_material_property_type),pointer :: geomech_material_property
  type(waypoint_type), pointer                 :: waypoint
  
  character(len=MAXWORDLENGTH)                 :: word
  character(len=MAXWORDLENGTH)                 :: card
  character(len=1) :: backslash

  PetscReal :: temp_real, temp_real2
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''
  
  discretization => geomech_realization%discretization
  
  
  do
    call InputReadFlotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)
    write(*,*),'word :: ',trim(word)

    select case(trim(word))
      !.........................................................................
      ! Read geomechanics grid information
      case ('GEOMECHANICS_GRID')
        call InputSkipToEND(input,option,trim(word))

       !.........................................................................
      ! Read geomechanics material information
      case ('GEOMECHANICS_MATERIAL_PROPERTY')
        geomech_material_property => GeomechanicsMaterialPropertyCreate()

        call InputReadWord(input,option,geomech_material_property%name, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'name','GEOMECHANICS_MATERIAL_PROPERTY')
        call GeomechanicsMaterialPropertyRead(geomech_material_property,input, &
                                              option)
        call GeomechanicsMaterialPropertyAddToList(geomech_material_property, &
                                geomech_realization%geomech_material_properties)
        nullify(geomech_material_property)
             
      !.........................................................................
      case default
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select
  enddo
  
end subroutine GeomechanicsInitReadInput
 
end module Geomechanics_Init_module
#endif
! GEOMECH
