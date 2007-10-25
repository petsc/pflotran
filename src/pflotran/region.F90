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
  
  integer, save :: num_regions = 0
  
  interface createRegion
    module procedure createRegionWithBlock
    module procedure createRegionWithList
    module procedure createRegionWithNothing
  end interface createRegion
  
  public :: createRegion, destroyRegion
  
contains

! ************************************************************************** !
!
! createRegionWithNothing: Creates a region with no arguments
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function createRegionWithNothing()

  implicit none
  
  type(region_type), pointer :: createRegionWithNothing
  
  allocate(createRegionWithNothing)
  createRegionWithNothing%id = 0
  createRegionWithNothing%name = ""
  createRegionWithNothing%i1 = 0
  createRegionWithNothing%i2 = 0
  createRegionWithNothing%j1 = 0
  createRegionWithNothing%j2 = 0
  createRegionWithNothing%k1 = 0
  createRegionWithNothing%k2 = 0
  createRegionWithNothing%num_cells = 0
  nullify(createRegionWithNothing%cell_ids)
  nullify(createRegionWithNothing%next)
  
  num_regions = num_regions + 1
  
  createRegionWithNothing%id = num_regions

end function createRegionWithNothing

! ************************************************************************** !
!
! createRegionWithBlock: Creates a region with i,j,k indices for arguments
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function createRegionWithBlock(i1,i2,j1,j2,k1,k2)

  implicit none
  
  integer :: i1, i2, j1, j2, k1, k2
  
  type(region_type), pointer :: createRegionWithBlock
  
  createRegionWithBlock => createRegionWithNothing()
  createRegionWithBlock%i1 = i1
  createRegionWithBlock%i2 = i2
  createRegionWithBlock%j1 = j2
  createRegionWithBlock%j2 = j2
  createRegionWithBlock%k1 = k1
  createRegionWithBlock%k2 = k2
  createRegionWithBlock%num_cells = (abs(i2-i1)+1)*(abs(j2-j1)+1)* &
                                    (abs(k2-k1)+1)

end function createRegionWithBlock

! ************************************************************************** !
!
! createRegion: Creates a region from a list of cells
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function createRegionWithList(list)

  implicit none
  
  integer :: list(:)
  
  type(region_type), pointer :: createRegionWithList
  
  createRegionWithList => createRegionWithNothing()
  createRegionWithList%num_cells = size(list)
  allocate(createRegionWithList%cell_ids(createRegionWithList%num_cells))
  createRegionWithList%cell_ids = list

end function createRegionWithList

! ************************************************************************** !
!
! destroyRegion: Deallocates a region
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine destroyRegion(region)

  implicit none
  
  type(region_type), pointer :: region
  
  if (associated(region%cell_ids)) deallocate(region%cell_ids)
  nullify(region%cell_ids)
  nullify(region%next)
  deallocate(region)
  nullify(region)

end subroutine destroyRegion

end module Region_module
