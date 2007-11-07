module Region_module
 
  implicit none

  private

#include "definitions.h"
 
  type, public :: block_type
    integer :: i1,i2,j1,j2,k1,k2    
    type(block_type), pointer :: next
  end type block_type
 
  type, public :: region_type
    integer :: id
    character(len=MAXWORDLENGTH) :: name
    type(block_type), pointer :: block_list
    integer :: i1,i2,j1,j2,k1,k2
    integer :: num_cells
    integer, pointer :: cell_ids(:)
    type(region_type), pointer :: next
  end type region_type
  
  type, public :: region_ptr_type
    type(region_type), pointer :: ptr
  end type region_ptr_type
  
  type, public :: region_list_type
    integer :: num_regions
    type(region_type), pointer :: first
    type(region_type), pointer :: last
    type(region_type), pointer :: array(:)
  end type region_list_type
  
  integer, save :: num_regions = 0
  
  interface RegionCreate
    module procedure RegionCreateWithBlock
    module procedure RegionCreateWithList
    module procedure RegionCreateWithNothing
  end interface RegionCreate
  
  interface RegionReadFromFile
    module procedure RegionReadFromInputFile
    module procedure RegionReadFromExternalFile
  end interface RegionReadFromFile
  
  public :: RegionCreate, RegionDestroy, RegionAddToList, RegionReadFromFile, &
            RegionInitList, RegionDestroyList, RegionGetPtrFromList
  
contains

! ************************************************************************** !
!
! RegionCreateWithNothing: Creates a region with no arguments
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function RegionCreateWithNothing()

  implicit none
  
  type(region_type), pointer :: RegionCreateWithNothing
  
  allocate(RegionCreateWithNothing)
  RegionCreateWithNothing%id = 0
  RegionCreateWithNothing%name = ""
  RegionCreateWithNothing%i1 = 0
  RegionCreateWithNothing%i2 = 0
  RegionCreateWithNothing%j1 = 0
  RegionCreateWithNothing%j2 = 0
  RegionCreateWithNothing%k1 = 0
  RegionCreateWithNothing%k2 = 0
  RegionCreateWithNothing%num_cells = 0
  nullify(RegionCreateWithNothing%cell_ids)
  nullify(RegionCreateWithNothing%next)
  
  num_regions = num_regions + 1
  
  RegionCreateWithNothing%id = num_regions

end function RegionCreateWithNothing

! ************************************************************************** !
!
! RegionCreateWithBlock: Creates a region with i,j,k indices for arguments
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function RegionCreateWithBlock(i1,i2,j1,j2,k1,k2)

  implicit none
  
  integer :: i1, i2, j1, j2, k1, k2
  
  type(region_type), pointer :: RegionCreateWithBlock
  
  RegionCreateWithBlock => RegionCreateWithNothing()
  RegionCreateWithBlock%i1 = i1
  RegionCreateWithBlock%i2 = i2
  RegionCreateWithBlock%j1 = j2
  RegionCreateWithBlock%j2 = j2
  RegionCreateWithBlock%k1 = k1
  RegionCreateWithBlock%k2 = k2
  RegionCreateWithBlock%num_cells = (abs(i2-i1)+1)*(abs(j2-j1)+1)* &
                                    (abs(k2-k1)+1)

end function RegionCreateWithBlock

! ************************************************************************** !
!
! RegionCreate: Creates a region from a list of cells
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function RegionCreateWithList(list)

  implicit none
  
  integer :: list(:)
  
  type(region_type), pointer :: RegionCreateWithList
  
  RegionCreateWithList => RegionCreateWithNothing()
  RegionCreateWithList%num_cells = size(list)
  allocate(RegionCreateWithList%cell_ids(RegionCreateWithList%num_cells))
  RegionCreateWithList%cell_ids = list

end function RegionCreateWithList

! ************************************************************************** !
!
! RegionInitList: Initializes a region list
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine RegionInitList(list)

  implicit none

  type(region_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_regions = 0

end subroutine RegionInitList

! ************************************************************************** !
!
! RegionAddToList: Adds a new region to a region list
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine RegionAddToList(new_region,list)

  implicit none
  
  type(region_type), pointer :: new_region
  type(region_list_type) :: list
  
  list%num_regions = list%num_regions + 1
  new_region%id = list%num_regions
  if (.not.associated(list%first)) list%first => new_region
  if (associated(list%last)) list%last%next => new_region
  list%last => new_region
  
end subroutine RegionAddToList

! ************************************************************************** !
!
! RegionReadFromExternalFile: Reads a list of cells from an external file
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine RegionReadFromExternalFile(region,filename)

  use Fileio_module
  use Utility_module
  
  implicit none
  
  type(region_type), pointer :: region
  character(len=MAXNAMELENGTH) :: filename
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  integer :: fid
  integer, pointer :: temp_int_array(:)
  integer :: max_size = 1000
  integer :: count, temp_int, ierr
  
  fid = 86
  open(unit=fid,file=filename)
  
  allocate(temp_int_array(max_size))
  temp_int_array = 0
  
  count = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    if (ierr /= 0) exit
    do
      if (ierr /= 0) exit
      ierr = 0
      call fiReadInt(string,temp_int,ierr)
      if (ierr == 0) then
        count = count + 1
        temp_int_array(count) = temp_int
      endif
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(temp_int_array,max_size) 
      endif
    enddo
  enddo
  if (count > 0) then
    region%num_cells = count
    allocate(region%cell_ids(count))
    region%cell_ids(1:count) = temp_int_array(1:count)
  else 
    region%num_cells = 0
    nullify(region%cell_ids)
  endif
  
  deallocate(temp_int_array)
          
  close(fid)          

end subroutine RegionReadFromExternalFile

! ************************************************************************** !
!
! RegionReadFromInputFile: Reads a list of cells from the pflotran input file
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine RegionReadFromInputFile(region,fid)

  use Fileio_module
  use Utility_module
  
  implicit none
  
  type(region_type), pointer :: region
  integer :: fid
  
  logical :: continuation_flag
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  integer, pointer :: temp_int_array(:)
  integer :: max_size = 1000
  integer :: count, temp_int, ierr

  allocate(temp_int_array(max_size))
  temp_int_array = 0
  
  count = 0
  continuation_flag = .true.
  do
    if (.not.continuation_flag) exit
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit
    continuation_flag = .false.
    if (index(string,"\") > 0) continuation_flag = .true.
    ierr = 0
    do
      if (ierr /= 0) exit
      call fiReadInt(string,temp_int,ierr)
      if (ierr == 0) then
        count = count + 1
        temp_int_array(count) = temp_int
      endif
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(temp_int_array,max_size) 
      endif
    enddo
  enddo
  if (count > 0) then
    region%num_cells = count
    allocate(region%cell_ids(count))
    region%cell_ids(1:count) = temp_int_array(1:count)
  else
    region%num_cells = 0
    nullify(region%cell_ids)
  endif

  deallocate(temp_int_array) 

end subroutine RegionReadFromInputFile

! ************************************************************************** !
!
! RegionGetPtrFromList: Returns a pointer to the region matching region_name
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
function RegionGetPtrFromList(region_name,region_list)

  use Fileio_module

  implicit none
  
  type(region_type), pointer :: RegionGetPtrFromList
  character(len=MAXNAMELENGTH) :: region_name
  type(region_list_type) :: region_list

  type(region_type), pointer :: region
    
  nullify(RegionGetPtrFromList)
  region => region_list%first
  
  do 
    if (.not.associated(region)) exit
    if (fiStringCompare(region%name,region_name,len_trim(region_name))) then
      RegionGetPtrFromList => region
      return
    endif
    region => region%next
  enddo
  
end function RegionGetPtrFromList

! ************************************************************************** !
!
! RegionDestroyList: Deallocates a list of regions
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine RegionDestroyList(region_list)

  implicit none
  
  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: region, prev_region
  
  
  region => region_list%first
  do 
    if (.not.associated(region)) exit
    prev_region => region
    region => region%next
    call RegionDestroy(prev_region)
  enddo
  
  region_list%num_regions = 0
  nullify(region_list%first)
  nullify(region_list%last)
  deallocate(region_list%array)
  nullify(region_list%array)
  
  deallocate(region_list)
  nullify(region_list)

end subroutine RegionDestroyList

! ************************************************************************** !
!
! RegionDestroy: Deallocates a region
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine RegionDestroy(region)

  implicit none
  
  type(region_type), pointer :: region
  
  if (.not.associated(region)) return
  
  if (associated(region%cell_ids)) deallocate(region%cell_ids)
  nullify(region%cell_ids)
  nullify(region%next)
  
  deallocate(region)
  nullify(region)

end subroutine RegionDestroy

end module Region_module
