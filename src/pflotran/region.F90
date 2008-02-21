module Region_module
 
  implicit none

  private

#include "definitions.h"
 
  type, public :: block_type
    PetscInt :: i1,i2,j1,j2,k1,k2    
    type(block_type), pointer :: next
  end type block_type
 
  type, public :: region_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: filename
    type(block_type), pointer :: block_list
    PetscInt :: i1,i2,j1,j2,k1,k2
    PetscInt :: iface
    PetscInt :: num_cells
    PetscInt, pointer :: cell_ids(:)
    PetscInt, pointer :: faces(:)
    type(region_type), pointer :: next
  end type region_type
  
  type, public :: region_ptr_type
    type(region_type), pointer :: ptr
  end type region_ptr_type
  
  type, public :: region_list_type
    PetscInt :: num_regions
    type(region_type), pointer :: first
    type(region_type), pointer :: last
    type(region_type), pointer :: array(:)
  end type region_list_type
  
  PetscInt, save :: num_regions = 0
  
  interface RegionCreate
    module procedure RegionCreateWithBlock
    module procedure RegionCreateWithList
    module procedure RegionCreateWithNothing
  end interface RegionCreate
  
  interface RegionReadFromFile
    module procedure RegionReadFromFileId
    module procedure RegionReadFromFilename
  end interface RegionReadFromFile
  
  public :: RegionCreate, RegionDestroy, RegionAddToList, RegionReadFromFile, &
            RegionInitList, RegionDestroyList, RegionGetPtrFromList, RegionRead
  
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
  
  type(region_type), pointer :: region
  
  allocate(region)
  region%id = 0
  region%name = ""
  region%filename = ""
  region%i1 = 0
  region%i2 = 0
  region%j1 = 0
  region%j2 = 0
  region%k1 = 0
  region%k2 = 0
  region%iface = 0
  region%num_cells = 0
  nullify(region%cell_ids)
  nullify(region%faces)
  nullify(region%next)
  
  num_regions = num_regions + 1
  
  region%id = num_regions
  
  RegionCreateWithNothing => region

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
  
  PetscInt :: i1, i2, j1, j2, k1, k2
  
  type(region_type), pointer :: RegionCreateWithBlock

  type(region_type), pointer :: region
  
  region => RegionCreateWithNothing()
  region%i1 = i1
  region%i2 = i2
  region%j1 = j2
  region%j2 = j2
  region%k1 = k1
  region%k2 = k2
  region%num_cells = (abs(i2-i1)+1)*(abs(j2-j1)+1)* &
                                    (abs(k2-k1)+1)
                                    
  RegionCreateWithBlock => region                                    

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
  
  PetscInt :: list(:)
  
  type(region_type), pointer :: RegionCreateWithList
  
  type(region_type), pointer :: region

  region => RegionCreateWithNothing()
  region%num_cells = size(list)
  allocate(region%cell_ids(region%num_cells))
  region%cell_ids = list
  
  RegionCreateWithList => region

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
! RegionRead: Reads a region from the input file
! author: Glenn Hammond
! date: 02/20/07
!
! ************************************************************************** !
subroutine RegionRead(region,string,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(region_type) :: region
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: length
  PetscErrorCode :: ierr

  call fiReadWord(string,region%name,.true.,ierr)
  call fiErrorMsg(option%myrank,'regn','name',ierr) 
  call printMsg(option,region%name)

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','REGION', ierr)
    call fiWordToUpper(word)   

    select case(trim(word))
    
      case('BLOCK')
        call fiReadInt(string,region%i1,ierr) 
        if (ierr /= 0) then
          ierr = 0
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'REGN',ierr)
          call fiReadInt(string,region%i1,ierr) 
        endif
        call fiErrorMsg(option%myrank,'i1','REGN', ierr)
        call fiReadInt(string,region%i2,ierr)
        call fiErrorMsg(option%myrank,'i2','REGN', ierr)
        call fiReadInt(string,region%j1,ierr)
        call fiErrorMsg(option%myrank,'j1','REGN', ierr)
        call fiReadInt(string,region%j2,ierr)
        call fiErrorMsg(option%myrank,'j2','REGN', ierr)
        call fiReadInt(string,region%k1,ierr)
        call fiErrorMsg(option%myrank,'k1','REGN', ierr)
        call fiReadInt(string,region%k2,ierr)
        call fiErrorMsg(option%myrank,'k2','REGN', ierr)
      case('FILE')
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(option%myrank,'filename','REGN', ierr)
        region%filename = word
      case('FACE')
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(option%myrank,'face','REGN', ierr)
        call fiWordToUpper(word)
        select case(word)
          case('WEST')
            region%iface = WEST_FACE
          case('EAST')
            region%iface = EAST_FACE
          case('NORTH')
            region%iface = NORTH_FACE
          case('SOUTH')
            region%iface = SOUTH_FACE
          case('BOTTOM')
            region%iface = BOTTOM_FACE
          case('TOP')
            region%iface = TOP_FACE
        end select
      case('END')
        exit        
      case default
        call printErrMsg(option,"REGION keyword not recognized")
    end select
  enddo
 
end subroutine RegionRead

! ************************************************************************** !
!
! RegionReadFromFilename: Reads a list of cells from a file named filename
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine RegionReadFromFilename(region,filename)

  use Fileio_module
  use Utility_module
  
  implicit none
  
  type(region_type), pointer :: region
  character(len=MAXNAMELENGTH) :: filename
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt :: fid
  PetscInt, pointer :: temp_int_array(:)
  PetscInt :: max_size = 1000
  PetscInt :: count, temp_int
  PetscErrorCode :: ierr
  
  fid = 86
  open(unit=fid,file=filename)
  call RegionReadFromFileId(region,fid)          
  close(fid)          

end subroutine RegionReadFromFilename

! ************************************************************************** !
!
! RegionReadFromFileId: Reads a list of cells from an open file
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine RegionReadFromFileId(region,fid)

  use Fileio_module
  use Utility_module
  
  implicit none
  
  type(region_type), pointer :: region
  PetscInt :: fid
  
  logical :: continuation_flag
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=1) :: backslash

  PetscInt, pointer :: temp_int_array(:)
  PetscInt :: max_size = 1000
  PetscInt :: count, temp_int
  PetscErrorCode :: ierr

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  allocate(temp_int_array(max_size))
  temp_int_array = 0
  
  count = 0
  continuation_flag = .true.
  do
    if (.not.continuation_flag) exit
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit
    continuation_flag = .false.
    if (index(string,backslash) > 0) continuation_flag = .true.
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

end subroutine RegionReadFromFileId

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
  PetscInt :: length
  type(region_list_type) :: region_list

  type(region_type), pointer :: region
    
  nullify(RegionGetPtrFromList)
  region => region_list%first
  
  do 
    if (.not.associated(region)) exit
    length = len_trim(region_name)
    if (length == len_trim(region%name) .and. &
        fiStringCompare(region%name,region_name,length)) then
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
  
  if (.not.associated(region_list)) return
  
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
  if (associated(region_list%array)) deallocate(region_list%array)
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
  if (associated(region%faces)) deallocate(region%faces)
  nullify(region%faces)
  nullify(region%next)
  
  deallocate(region)
  nullify(region)

end subroutine RegionDestroy

end module Region_module
