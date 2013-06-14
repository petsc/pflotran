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
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechicsInitReadRequiredCards(geomech_realization)

  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Input_module
  use String_module
  use Patch_module
  use Option_module
  
  implicit none
  
  type(geomech_realization_type)             :: geomech_realization
  type(geomech_discretization_type), pointer :: discretization
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(option_type), pointer          :: option
  type(input_type), pointer           :: input
  
  option         => geomech_realization%option  
  input          => geomech_realization%input
  
! Read in select required cards
!.........................................................................

  ! GEOMECHANICS information
  string = "GEOMECHANICS"
  call InputFindStringInFile(input,option,string)
  if(InputError(input)) return
  option%ngeomechdof = 3  ! displacements in x, y, z directions

end subroutine GeomechicsInitReadRequiredCards

! ************************************************************************** !
!
! GeomechanicsInit: Reads the required geomechanics data from input file
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsInit(geomech_realization,input,option)

  use Option_module
  use Input_module
  use String_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  
  implicit none
  
  type(geomech_realization_type)             :: geomech_realization
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_patch_type), pointer          :: patch
  type(input_type)                           :: input
  type(option_type)                          :: option
  character(len=MAXWORDLENGTH)               :: word
  
  discretization       => geomech_realization%discretization
  discretization%grid  => GMGridCreate()
  discretization%itype = UNSTRUCTURED_GRID ! Assuming only unstructured for now
  
  select case(discretization%itype)
    case(UNSTRUCTURED_GRID)
      patch => GeomechanicsPatchCreate()
      patch%geomech_grid => discretization%grid
      geomech_realization%geomech_patch => patch
  end select  
     
end subroutine GeomechanicsInit

! ************************************************************************** !
!
! GeomechanicsInitReadInput: Reads the geomechanics input data 
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsInitReadInput(geomech_realization,geomech_solver, &
                                     input,option)

  use Option_module
  use Input_module
  use String_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Material_module
  use Geomechanics_Region_module
  use Geomechanics_Debug_module
  use Geomechanics_Strata_module
  use Geomechanics_Condition_module
  use Geomechanics_Coupler_module
  use Solver_module
  use Waypoint_module

  ! Still need to add other geomech modules for output, etc once created
  
  implicit none
  
  type(geomech_realization_type)               :: geomech_realization
  type(solver_type)                            :: geomech_solver
  type(input_type)                             :: input
  type(option_type)                            :: option
  
  type(geomech_discretization_type), pointer   :: discretization
  type(geomech_material_property_type),pointer :: geomech_material_property
  type(waypoint_type), pointer                 :: waypoint
  type(geomech_grid_type), pointer             :: grid
  type(gm_region_type), pointer                :: region
  type(geomech_debug_type), pointer            :: debug
  type(geomech_strata_type), pointer           :: strata
  type(geomech_condition_type), pointer        :: condition
  type(geomech_coupler_type), pointer          :: coupler

  
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
  
  call GeomechanicsInit(geomech_realization,input,option)  
  
  if (associated(geomech_realization%geomech_patch)) grid => &
    geomech_realization%geomech_patch%geomech_grid
    
  do
    call InputReadFlotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)
    option%io_buffer = 'word :: ' // trim(word)
    call printMsg(option)   

    select case(trim(word))
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
      case ('GEOMECHANICS_REGION')
        region => GeomechRegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','GEOMECHANICS_REGION')
        call printMsg(option,region%name)
        call GeomechRegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call GeomechRegionAddToList(region,geomech_realization%geomech_regions)
        nullify(region)
 
      !.........................................................................
      case ('GEOMECHANICS_CONDITION')
        condition => GeomechConditionCreate(option)
        call InputReadWord(input,option,condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'SURF_FLOW_CONDITION','name')
        call printMsg(option,condition%name)
        call GeomechConditionRead(condition,input,option)
        call GeomechConditionAddToList(condition,geomech_realization%geomech_conditions)
        nullify(condition)
        
     !.........................................................................
      case ('GEOMECHANICS_BOUNDARY_CONDITION')
        coupler =>  GeomechCouplerCreate(GM_BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Geomech Boundary Condition name')
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)
        
      !.........................................................................
      case ('GEOMECHANICS_SRC_SINK')
        coupler => GeomechCouplerCreate(GM_SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)
                 
      !.........................................................................
      case('NEWTON_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('GEOMECHANICS')
            call SolverReadNewton(geomech_solver,input,option)
        end select     
        
      !.........................................................................
      case ('GEOMECHANICS_DEBUG')
        call GeomechDebugRead(geomech_realization%debug,input,option)       
        
       !.........................................................................
      case ('GEOMECHANICS_STRATIGRAPHY','GEOMECHANICS_STRATA')
        strata => GeomechStrataCreate()
        call GeomechStrataRead(strata,input,option)
        call GeomechRealizAddStrata(geomech_realization,strata)
        nullify(strata)       
        
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
