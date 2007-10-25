module Coupler_module
 
  use Condition_module
  use Connection_module
  use Region_module
 
  implicit none

  private
 
#include "definitions.h"
 
  type, public :: coupler_type
    integer :: id                                       ! id of coupler
    character(len=MAXWORDLENGTH) :: condition_name      ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: region_name         ! character string defining name of region to be applied
    integer :: icondition                               ! id of condition in condition array/list
    integer :: iregion                                  ! id of region in region array/list
    type(condition_type), pointer :: flow_condition     ! pointer to condition in condition array/list
    type(condition_type), pointer :: transport_condition ! pointer to condition in condition array/list
    type(region_type), pointer :: region                ! pointer to region in region array/list
    type(connection_type), pointer :: connection        ! pointer to an array/list of connections
    type(coupler_type), pointer :: next                 ! pointer to next coupler
  end type coupler_type
  
  integer, save :: num_couplers = 0
  
contains

! ************************************************************************** !
!
! createCoupler: Creates a coupler
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function createCoupler()

  implicit none

  type(coupler_type), pointer :: createCoupler
  
  allocate(createCoupler)
  createCoupler%id = 0
  createCoupler%condition_name = ""
  createCoupler%region_name = ""
  createCoupler%icondition = 0
  createCoupler%iregion = 0
  nullify(createCoupler%flow_condition)
  nullify(createCoupler%transport_condition)
  nullify(createCoupler%region)
  nullify(createCoupler%connection)
  nullify(createCoupler%next)
  
  num_couplers = num_couplers + 1
  createCoupler%id = num_couplers

end function createCoupler

! ************************************************************************** !
!
! destroyCoupler: Destroys a coupler
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine destroyCoupler(coupler)

  implicit none
  
  type(coupler_type), pointer :: coupler
  
  if (associated(coupler%flow_condition)) then
    call destroyCondition(coupler%flow_condition)
    nullify(coupler%flow_condition)
  endif
  
  if (associated(coupler%transport_condition)) then
    call destroyCondition(coupler%transport_condition)
    nullify(coupler%transport_condition)
  endif

  if (associated(coupler%region)) then
    call destroyRegion(coupler%region)
    nullify(coupler%region)
  endif
  
  if (associated(coupler%connection)) then
    call destroyConnection(coupler%connection)
    nullify(coupler%connection)
  endif

end subroutine destroyCoupler

end module Coupler_module
