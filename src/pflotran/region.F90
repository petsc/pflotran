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
  
  type, public :: region_list_type
    integer :: num_regions
    type(region_type), pointer :: first
    type(region_type), pointer :: last
    type(region_type), pointer :: array(:)
  end type region_list_type
  
  integer, save :: num_regions = 0
  
  interface createRegion
    module procedure createRegionWithBlock
    module procedure createRegionWithList
    module procedure createRegionWithNothing
  end interface createRegion
  
  interface readRegionFromFile
    module procedure readRegionFromInputFile
    module procedure readRegionFromExternalFile
  end interface readRegionFromFile
  
  public :: createRegion, destroyRegion, addRegionToList, readRegionFromFile, &
            initRegionList
  
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
! initRegionList: Initializes a region list
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine initRegionList(list)

  implicit none

  type(region_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_regions = 0

end subroutine initRegionList

! ************************************************************************** !
!
! addRegionToList: Adds a new region to a region list
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine addRegionToList(new_region,list)

  implicit none
  
  type(region_type), pointer :: new_region
  type(region_list_type) :: list
  
  list%num_regions = list%num_regions + 1
  new_region%id = list%num_regions
  if (.not.associated(list%first)) list%first => new_region
  if (associated(list%last)) list%last%next => new_region
  list%last => new_region
  
end subroutine addRegionToList

! ************************************************************************** !
!
! readRegionFromExternalFile: Reads a list of cells from an external file
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine readRegionFromExternalFile(region,filename)

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
  integer :: count1, ierr
  
  fid = 86
  open(unit=fid,file=filename)
  
  count1 = 0
  do
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit
    do
      if (ierr /= 0) exit
      ierr = 0
      call fiReadInt(string,temp_int_array(count1),ierr)
      if (ierr == 0) count1 = count1 + 1
      if (count1 > max_size) then ! resize temporary array
        call reallocateIntArray(temp_int_array,max_size) 
      endif
    enddo
  enddo
  if (count1 > 0) then
    region%num_cells = count1
    allocate(region%cell_ids(count1))
    region%cell_ids(1:count1) = temp_int_array(1:count1)
  endif
          
  close(fid)          

end subroutine readRegionFromExternalFile

! ************************************************************************** !
!
! readRegionFromInputFile: Reads a list of cells from the pflotran input file
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine readRegionFromInputFile(region,fid)

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
  integer :: count1, ierr

  
  count1 = 0
  continuation_flag = .true.
  do
    if (.not.continuation_flag) exit
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit
    continuation_flag = .false.
    if (index("\",string) > 0) continuation_flag = .true.
    ierr = 0
    do
      if (ierr /= 0) exit
      call fiReadInt(string,temp_int_array(count1),ierr)
      if (ierr == 0) count1 = count1 + 1
      if (count1 > max_size) then ! resize temporary array
        call reallocateIntArray(temp_int_array,max_size) 
      endif
    enddo
  enddo
  if (count1 > 0) then
    region%num_cells = count1
    allocate(region%cell_ids(count1))
    region%cell_ids(1:count1) = temp_int_array(1:count1)
  endif
        

end subroutine readRegionFromInputFile

! ************************************************************************** !
!
! localizeRegions: Resticts regions to cells local to processor
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine localizeRegions(region_list,grid,option)

  use Option_module
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module

  implicit none
  
  type(region_list_type), pointer :: region_list
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  type(region_type), pointer :: region
  integer, allocatable :: temp_int_array(:)
  integer :: i, j, k, count, local_count, local_ghosted_id, local_id
  
  region => region_list%first
  do
  
    if (.not.associated(region)) exit
    
    if (.not.associated(region%cell_ids) .and. &
        region%i1 > 0 .and. region%i2 > 0 .and. &
        region%j1 > 0 .and. region%j2 > 0 .and. &
        region%k1 > 0 .and. region%k2 > 0) then
        
      region%num_cells = (min(region%i2,grid%structured_grid%nxe) - &
                          max(region%i1,grid%structured_grid%nxs+1)) * &
                         (min(region%j2,grid%structured_grid%nye) - &
                          max(region%j1,grid%structured_grid%nys+1)) * &
                         (min(region%k2,grid%structured_grid%nze) - &
                          max(region%k1,grid%structured_grid%nzs+1))
      allocate(region%cell_ids(region%num_cells))
        
      count = 0  
      do k=max(region%k1,grid%structured_grid%nzs+1), &
           min(region%k2,grid%structured_grid%nze)
        do j=max(region%j1,grid%structured_grid%nys+1), &
             min(region%j2,grid%structured_grid%nye)
          do i=max(region%i1,grid%structured_grid%nxs+1), &
               min(region%i2,grid%structured_grid%nxe)
            count = count + 1
            region%cell_ids(count) = &
                   i + j*grid%structured_grid%nlx + k*grid%structured_grid%nlxy
          enddo
        enddo
      enddo
      if (count /= region%num_cells) &
        call printErrMsg(option,"Mismatch in number of cells in block region")
    else
      allocate(temp_int_array(region%num_cells))
      if (grid%is_structured) then
        do count=1,region%num_cells
          i = mod(region%cell_ids(count),grid%structured_grid%nx) - &
                grid%structured_grid%nxs
          j = mod((region%cell_ids(count)-1)/grid%structured_grid%nx, &
                  grid%structured_grid%ny)+1 - &
                grid%structured_grid%nys
          k = ((region%cell_ids(count)-1)/grid%structured_grid%nxy)+1 - &
                grid%structured_grid%nzs
          if (i > 0 .and. i <= grid%structured_grid%nlx .and. &
              j > 0 .and. j <= grid%structured_grid%nly .and. &
              k > 0 .and. k <= grid%structured_grid%nlz) then
            temp_int_array(local_count) = &
                i + j*grid%structured_grid%nlx + k*grid%structured_grid%nlxy
            local_count = local_count + 1
          endif
        enddo
      else
        do count=1,region%num_cells
          local_ghosted_id = GetLocalGhostedIdFromHash(grid%unstructured_grid, &
                                                       region%cell_ids(count))
          if (local_ghosted_id > -1) then
            local_id = grid%nG2L(local_ghosted_id)
            if (local_id > -1) then
              temp_int_array(local_count) = local_id
              local_count = local_count + 1
            endif
          endif
        enddo
      endif
      if (local_count /= region%num_cells) then
        deallocate(region%cell_ids)
        allocate(region%cell_ids(local_count))
        region%cell_ids(1:local_count) = temp_int_array(1:local_count)
      endif
      deallocate(temp_int_array)
    endif
    
    if (region%num_cells == 0) deallocate(region%cell_ids)
    region => region%next
    
  enddo

end subroutine localizeRegions

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
