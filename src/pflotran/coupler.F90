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
    type(condition_type), pointer :: flow_condition_ptr ! pointer to condition in condition array/list
    type(condition_type), pointer :: transport_condition_ptr ! pointer to condition in condition array/list
    type(region_type), pointer :: region_ptr            ! pointer to region in region array/list
    type(connection_type), pointer :: connection_ptr    ! pointer to an array/list of connections
    type(coupler_type), pointer :: next                 ! pointer to next coupler
  end type coupler_type
  
end module Coupler_module
