module Region_module

  implicit none

#include "definitions.h"

  private

  type, public :: region
  
    integer :: id          ! id of region
    character(len=MAXWORDLENGTH) :: name
    
    integer :: num_cells   ! for structure/unstructured
    integer, pointer :: cell_ids(:)
    
    integer :: i1,i2,j1,j2,k1,k2  ! for structured
    
    type(region), pointer :: next ! pointer for linked list
  
  end type

  ! pointer data structure required for making an array of region pointers in F90
  type :: region_ptr
    type(region), pointer :: ptr
  end type region_ptr

  ! linked list of regions
  type(region), pointer, save :: region_list
  ! pointer to last region; facilitates appendage of regions
  type(region), pointer, save :: last_region
  ! array of pointers to regions in list
  type(region_ptr), pointer, save :: region_array(:)
  ! number of regions in list
  integer, save :: num_regions
  
  public :: initRegionModule, readRegion, getRegionFromName
  
contains

! ************************************************************************** !
!
! InitRegionModule: Initializes module variables, lists, arrays.
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine initRegionModule()

  implicit none

  nullify(region_list)
  nullify(last_region)
  nullify(region_array)
  num_regions = 0

end subroutine InitRegionModule

! ************************************************************************** !
!
! readRegion: Reads a region from a string
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
integer function readRegion(string)

  use fileio_module

  implicit none
  
  integer :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  word = adjustl(string(1:min(31,len_trim(string))))  ! remove leading spaces
  call fiCharsToLower(string,4)
  if (fiStringCompare(string,'file',4)) then ! unstructured
    call fiReadWord(string,word,.true.,ierr) ! strip 'file' from string
    call fiReadWord(string,word,.true.,ierr) ! filename
    call fiDefaultMsg('file:name',ierr)
    readRegion = readUnstructuredRegion(word)
  else
    readRegion = readStructuredRegion(string)
  endif

end function readRegion

! ************************************************************************** !
!
! readUnstructuredRegion: Reads a generalized region from a file
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
integer function readUnstructuredRegion(filename)

  use fileio_module
  
  implicit none
  integer :: fid
  character(len=MAXWORDLENGTH) :: filename

  fid = 86
!  open(unit=fid,file=filename,action='read',status='old')
!  close(fid)
  
  print *, 'Code for region.F90:readUnstructuredRegion still to be written'
  stop
  
  readUnstructuredRegion = 0
  

end function readUnstructuredRegion

! ************************************************************************** !
!
! readStructuredRegion: Reads a structured region from i,j,k block
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
integer function readStructuredRegion(string)

  use fileio_module
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: string

  type(region), pointer :: new_region
  integer :: ierr
  
  new_region = createRegion()
  
  call fiReadInt(string,new_region%i1,ierr) 
  call fiDefaultMsg('i1',ierr)
  call fiReadInt(string,new_region%i2,ierr) 
  call fiDefaultMsg('i2',ierr)
  call fiReadInt(string,new_region%j1,ierr) 
  call fiDefaultMsg('j1',ierr)
  call fiReadInt(string,new_region%j2,ierr) 
  call fiDefaultMsg('j2',ierr)
  call fiReadInt(string,new_region%k1,ierr) 
  call fiDefaultMsg('k1',ierr)
  call fiReadInt(string,new_region%k2,ierr) 
  call fiDefaultMsg('k2',ierr)
  
  call addRegionToList(new_region)
  
  readStructuredRegion = new_region%id

end function readStructuredRegion

! ************************************************************************** !
!
! getRegionFromName: Returns the pointer of the first region in the  
!                    region list with a matching name
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
function getRegionFromName(name)

  use fileio_module
  
  implicit none
  
  type(region), pointer :: getRegionFromName
  character(len=MAXWORDLENGTH) ::  name
  
  type(region), pointer :: cur_region

  nullify(getRegionFromName)
  
  cur_region => region_list
  do 
    if (.not.associated(cur_region)) exit
    if (fiStringCompare(cur_region%name,name,len_trim(name))) then
      getRegionFromName => cur_region 
      return
    endif
    cur_region => cur_region%next
  enddo

end function getRegionFromName

! ************************************************************************** !
!
! createRegion: Allocates and initializes a new region
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function createRegion()

  implicit none
  
  type(region), pointer :: createRegion
  allocate(createRegion)
  createRegion%id = 0
  createRegion%num_cells = 0
  nullify(createRegion%cell_ids)
  createRegion%i1 = 0
  createRegion%i2 = 0
  createRegion%j1 = 0
  createRegion%j2 = 0
  createRegion%k1 = 0
  createRegion%k2 = 0

end function createRegion

! ************************************************************************** !
!
! addRegionToList: Adds a new region of the module global list of regions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine addRegionToList(new_region)

  implicit none
  
  type(region), pointer :: new_region
  
  num_regions = num_regions + 1
  new_region%id = num_regions
  if (.not.associated(region_list)) region_list => new_region
  if (associated(last_region)) last_region%next => new_region
  last_region => new_region
  
end subroutine addRegionToList

! ************************************************************************** !
!
! convertRegionListToArray: Creates an array of pointers to the regions
!                           in the region list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine convertRegionListToArray()

  implicit none
  
  integer :: count
  type(region), pointer :: cur_region
  
  allocate(region_array(num_regions))
  
  cur_region => region_list
  do 
    if (.not.associated(cur_region)) exit
    region_array(cur_region%id)%ptr => cur_region
    cur_region => cur_region%next
  enddo

end subroutine convertRegionListToArray

! ************************************************************************** !
!
! destroyRegionList: Deallocates the module global list and array of regions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine destroyRegionList

  implicit none
  
  type(region), pointer :: cur_region, prev_region
  
  deallocate(region_array)
  nullify(region_array)
  
  cur_region => region_list
  do 
    if (.not.associated(cur_region)) exit
    prev_region => cur_region
    cur_region => cur_region%next
    deallocate(prev_region)
  enddo
  
  nullify(region_list)
  nullify(last_region)
  num_regions = 0

end subroutine destroyRegionList

end module Region_module