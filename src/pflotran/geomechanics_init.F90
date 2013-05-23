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


  public :: GeomechicsInitReadRequiredCards

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
  use Geomech_Grid_module
  use Input_module
  use String_module
  use Patch_module
  use Level_module
  use Option_module
  
  implicit none
  
  type(geomech_realization_type)             :: geomech_realization
  type(geomech_discretization_type), pointer :: discretization
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(patch_type), pointer   :: patch
  type(level_type), pointer   :: level
  type(option_type), pointer  :: option
  type(input_type), pointer   :: input
  
  patch          => geomech_realization%patch
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
  
  string = "GEOMECHANICS_GRID"
  call InputFindStringInFile(input,option,string)
!  call SurfaceFlowReadRequiredCardsFromInput(surf_realization,input,option)
  call GeomechanicsInit(geomech_realization,input,option)  

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%geomech_grid => discretization%grid
      if (.not.associated(geomech_realization%level_list)) then
        geomech_realization%level_list => LevelCreateList()
      endif
      level => LevelCreate()
      call LevelAddToList(level,geomech_realization%level_list)
      call PatchAddToList(patch,level%patch_list)
      geomech_realization%patch => patch
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
  
  input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  call InputReadFlotranString(input,option)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
  call StringToUpper(word)
    
  select case(trim(word))
    case ('TYPE')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'keyword','TYPE')
      call StringToUpper(word)

      select case(trim(word))
        case ('UNSTRUCTURED')
          discretization%itype = UNSTRUCTURED_GRID
          discretization%ctype = 'unstructured'
          call InputReadNChars(input,option, &
                               discretization%filename, &
                               MAXSTRINGLENGTH, &
                               PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','filename')
          discretization%grid => GMGridCreate()
          call GeomechGridRead(discretization%grid,discretization%filename, &
                               option)
        case default
          option%io_buffer = 'Geomechanics supports only unstructured grid'
          call printErrMsg(option)          
      end select
  end select       
     
  
  
  
  
  
  
  
  

end subroutine GeomechanicsInit

end module Geomechanics_Init_module
#endif
! GEOMECH
