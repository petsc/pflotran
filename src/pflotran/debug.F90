module Debug_module

  implicit none
  
  private
  
  type, public :: pflow_debug_type
    logical :: vecview_residual
    logical :: vecview_solution
    logical :: matview_Jacobian
    logical :: matview_Jacobian_detailed
    logical :: norm_Jacobian

    logical :: print_numerical_derivatives

    logical :: print_couplers
  end type pflow_debug_type
  
  type, public :: ptran_debug_type
    logical :: placeholder
  end type ptran_debug_type

  interface DebugRead
    module procedure DebugReadPflow
  end interface DebugRead
  
  
  public :: DebugCreatePflow, DebugRead
  
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
  
  debug%vecview_residual = .false.
  debug%vecview_solution = .false.
  debug%matview_Jacobian = .false.
  debug%matview_Jacobian_detailed = .false.
  debug%norm_Jacobian = .false.
  
  debug%print_numerical_derivatives = .false.
  
  debug%print_couplers = .false.

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
  
  DebugCreatePtran => debug

end function DebugCreatePtran

! ************************************************************************** !
!
! DebugReadPflow: Reads debugging data from the input file
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
subroutine DebugReadPflow(debug,fid,myrank)

  use Fileio_module
  
  implicit none

#include "definitions.h"
    
  type(pflow_debug_type) :: debug
  integer :: fid
  integer :: myrank
  
  character(len=MAXSTRINGLENGTH) :: string, error_string
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  integer :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',3)) exit  

    call fiReadWord(string,keyword,.true.,ierr)
    call fiErrorMsg(myrank,'keyword','DEBUG', ierr)   
      
    select case(trim(keyword))
    
      case('PRINT_SOLUTION','VECVIEW_SOLUTION','VIEW_SOLUTION')
        debug%vecview_solution = .true.
      case('PRINT_RESIDUAL','VECVIEW_RESIDUAL','VIEW_RESIDUAL')
        debug%vecview_residual = .true.
      case('PRINT_JACOBIAN','MATVIEW_JACOBIAN','VIEW_JACOBIAN')
        debug%matview_Jacobian = .true.
      case('PRINT_JACOBIAN_NORM','NORM_JACOBIAN')
        debug%norm_Jacobian = .true.
      case('PRINT_COUPLERS','PRINT_COUPLER')
        debug%print_couplers = .true.
      case('PRINT_JACOBIAN_DETAILED','MATVIEW_JACOBIAN_DETAILED','VIEW_JACOBIAN_DETAILED')
        debug%matview_Jacobian_detailed = .true.
      case('PRINT_NUMERICAL_DERIVATIVES','VIEW_NUMERICAL_DERIVATIVES')
        debug%print_numerical_derivatives = .true.

    end select 
  
  enddo  

end subroutine DebugReadPflow

end module Debug_module
