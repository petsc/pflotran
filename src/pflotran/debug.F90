module Debug_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public :: debug_type
    PetscBool :: vecview_residual
    PetscBool :: vecview_solution
    PetscBool :: matview_Jacobian
    PetscBool :: matview_Jacobian_detailed
    PetscBool :: norm_Jacobian

    PetscBool :: print_numerical_derivatives

    PetscBool :: print_couplers
    character(len=MAXSTRINGLENGTH) :: coupler_string
    PetscBool :: print_waypoints
  end type debug_type

  public :: DebugCreate, DebugRead
  
contains

! ************************************************************************** !
!
! DebugCreate: Create object that stores debugging options for PFLOW
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
function DebugCreate()

  implicit none
  
  type(debug_type), pointer :: DebugCreate
  
  type(debug_type), pointer :: debug
  
  allocate(debug)
  
  debug%vecview_residual = PETSC_FALSE
  debug%vecview_solution = PETSC_FALSE
  debug%matview_Jacobian = PETSC_FALSE
  debug%matview_Jacobian_detailed = PETSC_FALSE
  debug%norm_Jacobian = PETSC_FALSE
  
  debug%print_numerical_derivatives = PETSC_FALSE
  
  debug%print_couplers = PETSC_FALSE
  debug%coupler_string = ''
  debug%print_waypoints = PETSC_FALSE

  DebugCreate => debug

end function DebugCreate

! ************************************************************************** !
!
! DebugRead: Reads debugging data from the input file
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
subroutine DebugRead(debug,input,option)

  use Option_module
  use Input_Aux_module
  
  implicit none
    
  type(debug_type) :: debug
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
        debug%coupler_string = trim(adjustl(input%buf))
      case('PRINT_JACOBIAN_DETAILED','MATVIEW_JACOBIAN_DETAILED','VIEW_JACOBIAN_DETAILED')
        debug%matview_Jacobian_detailed = PETSC_TRUE
      case('PRINT_NUMERICAL_DERIVATIVES','VIEW_NUMERICAL_DERIVATIVES')
        debug%print_numerical_derivatives = PETSC_TRUE
      case('WAYPOINTS')
        debug%print_waypoints = PETSC_TRUE
      case default
        option%io_buffer = 'Option "' // trim(keyword) // &
          '" not recognized under DEBUG.'
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine DebugRead

end module Debug_module
