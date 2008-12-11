module Debug_module

  implicit none
  
  private
  
#include "definitions.h"

  type, public :: pflow_debug_type
    PetscTruth :: vecview_residual
    PetscTruth :: vecview_solution
    PetscTruth :: matview_Jacobian
    PetscTruth :: matview_Jacobian_detailed
    PetscTruth :: norm_Jacobian

    PetscTruth :: print_numerical_derivatives

    PetscTruth :: print_couplers
  end type pflow_debug_type
  
  type, public :: ptran_debug_type
    PetscTruth :: vecview_residual
    PetscTruth :: vecview_solution
    PetscTruth :: matview_Jacobian
    PetscTruth :: matview_Jacobian_detailed
    PetscTruth :: norm_Jacobian
    PetscTruth :: print_couplers    
  end type ptran_debug_type

  interface DebugRead
    module procedure DebugReadPflow
  end interface DebugRead
  
  
  public :: DebugCreatePflow, DebugCreatePtran, DebugRead
  
contains

! ************************************************************************** !
!
! DebugCreatePflow: Create object that stores debugging options for PFLOW
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
function DebugCreatePflow()

  implicit none
  
  type(pflow_debug_type), pointer :: DebugCreatePflow
  
  type(pflow_debug_type), pointer :: debug
  
  allocate(debug)
  
  debug%vecview_residual = PETSC_FALSE
  debug%vecview_solution = PETSC_FALSE
  debug%matview_Jacobian = PETSC_FALSE
  debug%matview_Jacobian_detailed = PETSC_FALSE
  debug%norm_Jacobian = PETSC_FALSE
  
  debug%print_numerical_derivatives = PETSC_FALSE
  
  debug%print_couplers = PETSC_FALSE

  DebugCreatePflow => debug

end function DebugCreatePflow

! ************************************************************************** !
!
! DebugCreatePtran: Create object that stores debugging options for PFLOW
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
function DebugCreatePtran()

  implicit none
  
  type(ptran_debug_type), pointer :: DebugCreatePtran
  
  type(ptran_debug_type), pointer :: debug
  
  allocate(debug)
  debug%vecview_residual = PETSC_FALSE
  debug%vecview_solution = PETSC_FALSE
  debug%matview_Jacobian = PETSC_FALSE
  debug%matview_Jacobian_detailed = PETSC_FALSE
  debug%norm_Jacobian = PETSC_FALSE
  
  debug%print_couplers = PETSC_FALSE

  
  DebugCreatePtran => debug

end function DebugCreatePtran

! ************************************************************************** !
!
! DebugReadPflow: Reads debugging data from the input file
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
subroutine DebugReadPflow(debug,input,option)

  use Option_module
  use Input_module
  
  implicit none
    
  type(pflow_debug_type) :: debug
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','DEBUG')   
      
    select case(trim(keyword))
    
      case('PRINT_SOLUTION','VECVIEW_SOLUTION','VIEW_SOLUTION')
        debug%vecview_solution = PETSC_TRUE
      case('PRINT_RESIDUAL','VECVIEW_RESIDUAL','VIEW_RESIDUAL')
        debug%vecview_residual = PETSC_TRUE
      case('PRINT_JACOBIAN','MATVIEW_JACOBIAN','VIEW_JACOBIAN')
        debug%matview_Jacobian = PETSC_TRUE
      case('PRINT_JACOBIAN_NORM','NORM_JACOBIAN')
        debug%norm_Jacobian = PETSC_TRUE
      case('PRINT_COUPLERS','PRINT_COUPLER')
        debug%print_couplers = PETSC_TRUE
      case('PRINT_JACOBIAN_DETAILED','MATVIEW_JACOBIAN_DETAILED','VIEW_JACOBIAN_DETAILED')
        debug%matview_Jacobian_detailed = PETSC_TRUE
      case('PRINT_NUMERICAL_DERIVATIVES','VIEW_NUMERICAL_DERIVATIVES')
        debug%print_numerical_derivatives = PETSC_TRUE

    end select 
  
  enddo  

end subroutine DebugReadPflow

end module Debug_module
