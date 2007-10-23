module Region_module
 
  implicit none

  private

#include "definitions.h"
 
  type, public :: region_type
    integer :: id
    character(len=MAXWORDLENGTH) :: name
    integer :: i1,i2,j1,j2,k1,k2
    integer :: num_cells
    integer, pointer :: cell_ids(:)
    type(region_type), pointer :: next
  end type region_type
  
  type, public :: region_ptr_type
    type(region_type), pointer :: ptr
  end type region_ptr_type
  
end module Region_module
