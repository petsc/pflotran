module Condition_module
 
  implicit none

  private
  
#include "definitions.h"
 
  type, public :: condition_type
    integer :: id                                 ! id from which condition can be referenced
    integer :: itype                              ! integer describing type of condition
    character(len=MAXWORDLENGTH) :: ctype         ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
    integer :: num_values                         ! number of entries in the arrays of values
    real*8, pointer :: times(:)                   ! array of times between which linear interpolation of values occurs
    real*8, pointer :: values(:,:)                ! array of condition values, size(ndof,max_time_index)
    real*8, pointer :: cur_value(:)               ! current value of condition a time t, size(ndof)
    real*8 :: datum(3)                            ! location of reference value(s) in domain
    real*8 :: gradient(3)                         ! rate at which reference value(s) change(s) over 3D space
    integer :: cur_index, max_index               ! current and maximum index in arrays
    type(condition_type), pointer :: next         ! pointer to next condition_type for linked-lists
  end type condition_type
  
end module Condition_module
