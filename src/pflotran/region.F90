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
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: i1,i2,j1,j2,k1,k2
    type(point3d_type), pointer :: coordinates(:)
    PetscInt :: iface
    PetscInt :: num_cells
    PetscInt, pointer :: cell_ids(:)
    PetscInt, pointer :: faces(:)
    PetscInt, pointer :: vert_ids(:,:) ! For Unstructured mesh
    PetscInt :: num_verts              ! For Unstructured mesh
    PetscInt :: grid_type  ! To identify whether region is applicable to a Structured or Unstructred mesh
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
  
  type, public :: point3d_type
    PetscReal :: x
    PetscReal :: y
    PetscReal :: z
  end type point3d_type
  
  interface RegionCreate
    module procedure RegionCreateWithBlock
    module procedure RegionCreateWithList
    module procedure RegionCreateWithNothing
    module procedure RegionCreateWithRegion    
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
  region%grid_type = STRUCTURED_GRID ! By default it is assumed that the region 
                   ! is applicable to strucutred grid, unless
                   ! explicitly stated in pflotran input file
  region%num_verts = 0
  nullify(region%coordinates)
  nullify(region%cell_ids)
  nullify(region%faces)
  nullify(region%vert_ids)
  nullify(region%next)
  
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
! RegionCreateWithRegion: Creates a copy of a region
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
function RegionCreateWithRegion(region)

  implicit none
  
  type(region_type), pointer :: RegionCreateWithRegion
  type(region_type), pointer :: region
  
  type(region_type), pointer :: new_region
  PetscInt :: icount
  
  new_region => RegionCreateWithNothing()
  
  new_region%id = region%id
  new_region%name = region%name
  new_region%filename = region%filename
  new_region%i1 = region%i1
  new_region%i2 = region%i2
  new_region%j1 = region%j1
  new_region%j2 = region%j2
  new_region%k1 = region%k1
  new_region%k2 = region%k2
  new_region%iface = region%iface
  new_region%num_cells = region%num_cells
  new_region%num_verts = region%num_verts
  new_region%grid_type = region%grid_type
  if (associated(region%coordinates)) then
    allocate(new_region%coordinates(size(region%coordinates)))
    do icount = 1, size(new_region%coordinates)
      new_region%coordinates(icount)%x = region%coordinates(icount)%x
      new_region%coordinates(icount)%y = region%coordinates(icount)%y
      new_region%coordinates(icount)%z = region%coordinates(icount)%z
    enddo
  endif
  if (associated(region%cell_ids)) then
    allocate(new_region%cell_ids(new_region%num_cells))
    new_region%cell_ids(1:new_region%num_cells) = region%cell_ids(1:region%num_cells)
  endif
  if (associated(region%faces)) then
    allocate(new_region%faces(new_region%num_cells))
    new_region%faces(1:new_region%num_cells) = region%faces(1:region%num_cells)
  endif
  if (associated(region%vert_ids)) then
    allocate(new_region%vert_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts))
    new_region%vert_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts) = &
    region%vert_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts)
  endif
  
  RegionCreateWithRegion => new_region
  
end function RegionCreateWithRegion
  
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
! date: 02/20/08
!
! ************************************************************************** !
subroutine RegionRead(region,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(region_type) :: region
  type(input_type) :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  PetscInt :: icount
  PetscInt, parameter :: max_num_coordinates = 30
  type(point3d_type) :: coordinates(max_num_coordinates)

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','REGION')
    call StringToUpper(keyword)   

    select case(trim(keyword))
    
      case('BLOCK')
        call InputReadInt(input,option,region%i1) 
        if (InputError(input)) then
          input%ierr = 0
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'REGION')
          call InputReadInt(input,option,region%i1) 
        endif
        call InputErrorMsg(input,option,'i1','REGION')
        call InputReadInt(input,option,region%i2)
        call InputErrorMsg(input,option,'i2','REGION')
        call InputReadInt(input,option,region%j1)
        call InputErrorMsg(input,option,'j1','REGION')
        call InputReadInt(input,option,region%j2)
        call InputErrorMsg(input,option,'j2','REGION')
        call InputReadInt(input,option,region%k1)
        call InputErrorMsg(input,option,'k1','REGION')
        call InputReadInt(input,option,region%k2)
        call InputErrorMsg(input,option,'k2','REGION')
      case('COORDINATE')
        allocate(region%coordinates(1))
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x) 
        if (InputError(input)) then
          input%ierr = 0
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'REGION')
          call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x)
        endif
        call InputErrorMsg(input,option,'x-coordinate','REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%y)
        call InputErrorMsg(input,option,'y-coordinate','REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%z)
        call InputErrorMsg(input,option,'z-coordinate','REGION')
      case('COORDINATES')
        icount = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'REGION')
          if (InputCheckExit(input,option)) exit              
          icount = icount + 1
          if (icount > max_num_coordinates) then
            write(option%io_buffer, &
                  '(''Number of coordinates in region '',a, &
                &'' exceeds limit of '',i3)') region%name, max_num_coordinates
            call printErrMsg(option)
          endif
          call InputReadDouble(input,option,coordinates(icount)%x) 
          call InputErrorMsg(input,option,'x-coordinate','REGION')
          call InputReadDouble(input,option,coordinates(icount)%y)
          call InputErrorMsg(input,option,'y-coordinate','REGION')
          call InputReadDouble(input,option,coordinates(icount)%z)
          call InputErrorMsg(input,option,'z-coordinate','REGION')
        enddo
        allocate(region%coordinates(icount))
        do icount = 1, size(region%coordinates)
          region%coordinates(icount)%x = coordinates(icount)%x
          region%coordinates(icount)%y = coordinates(icount)%y
          region%coordinates(icount)%z = coordinates(icount)%z
        enddo
      case('FILE')
        call InputReadNChars(input,option,region%filename,MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','REGION')
      case('LIST')
        option%io_buffer = 'REGION LIST currently not implemented'
        call printErrMsg(option)
      case('FACE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'face','REGION')
        call StringToUpper(word)
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
    case('GRID')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'GRID','REGION')
        call StringToUpper(word)
        select case(trim(word))
          case('STRUCTURED')
            region%grid_type = STRUCTURED_GRID
          case('UNSTRUCTURED')
            region%grid_type = UNSTRUCTURED_GRID
          case default
            option%io_buffer = 'REGION keyword: GRID = '//trim(word)//'not supported yet'
          call printErrMsg(option)
      end select
      case default
        option%io_buffer = 'REGION keyword: '//trim(keyword)//' not recognized'
        call printErrMsg(option)
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
subroutine RegionReadFromFilename(region,option,filename)

  use Input_module
  use Option_module
  use Utility_module
  
  implicit none
  
  type(region_type) :: region
  type(option_type) :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: filename
  
  input => InputCreate(IUNIT_TEMP,filename)
  call RegionReadFromFileId(region,input,option)          
  call InputDestroy(input)         

end subroutine RegionReadFromFilename

! ************************************************************************** !
!
! RegionReadFromFileId: Reads a list of cells from an open file
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine RegionReadFromFileId(region,input,option)

  use Input_module
  use Option_module
  use Utility_module
  use Logging_module
  
  implicit none
  
  type(region_type) :: region
  type(option_type) :: option
  type(input_type) :: input
  
  PetscBool :: continuation_flag
  character(len=MAXWORDLENGTH) :: word
  character(len=1) :: backslash

  PetscInt, pointer :: temp_int_array(:)
  PetscInt, pointer :: cell_ids_p(:)
  PetscInt, pointer :: face_ids_p(:)
  PetscInt, pointer :: vert_id_0_p(:)
  PetscInt, pointer :: vert_id_1_p(:)
  PetscInt, pointer :: vert_id_2_p(:)
  PetscInt, pointer :: vert_id_3_p(:)
  PetscInt, pointer :: vert_id_4_p(:)
  PetscInt :: max_size
  PetscInt :: count
  PetscInt :: temp_int
  PetscInt :: input_data_type
  PetscInt :: ii
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: remainder
  PetscErrorCode :: ierr

  PetscInt, parameter :: CELL_IDS_ONLY = 1
  PetscInt, parameter :: CELL_IDS_WITH_FACE_IDS = 2
  PetscInt, parameter :: VERTEX_IDS = 3

  call PetscLogEventBegin(logging%event_region_read_ascii,ierr)

  max_size = 1000
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  allocate(temp_int_array(max_size))
  allocate(cell_ids_p(max_size))
  allocate(face_ids_p(max_size))
  allocate(vert_id_0_p(max_size))
  allocate(vert_id_0_p(max_size))
  allocate(vert_id_1_p(max_size))
  allocate(vert_id_2_p(max_size))
  allocate(vert_id_3_p(max_size))
  allocate(vert_id_4_p(max_size))
  
  temp_int_array = 0
  cell_ids_p = 0
  face_ids_p = 0
  vert_id_0_p = 0
  vert_id_1_p = -1
  vert_id_2_p = -1
  vert_id_3_p = -1
  vert_id_4_p = -1
  
  
  count = 0
#if 0
  continuation_flag = PETSC_TRUE
  do
    if (.not.continuation_flag) exit
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    continuation_flag = PETSC_FALSE
    if (index(input%buf,backslash) > 0) &
      continuation_flag = PETSC_TRUE
    input%ierr = 0
    do
      if (InputError(input)) exit
      call InputReadInt(input,option,temp_int)
      if (.not.InputError(input)) then
        count = count + 1
        temp_int_array(count) = temp_int
      endif
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(temp_int_array,max_size) 
      endif
    enddo
  enddo
#endif

  ! Determine if region definition in the input data is one of the following:
  !  1) Contains cell ids only : Only ONE entry per line
  !  2) Contains cell ids and face ids: TWO entries per line
  !  3) Contains vertex ids that make up the face: MORE than two entries per
  !     line
  count = 0
  call InputReadFlotranString(input, option)
  do 
    call InputReadInt(input, option, temp_int)
    if(InputError(input)) exit
    count = count + 1
    temp_int_array(count) = temp_int
  enddo

  if(count == 1) then
    !
    ! Input data contains only cell ids
    !
    input_data_type = CELL_IDS_ONLY
    cell_ids_p(1) = temp_int_array(1)
    count = 1

    ! Read the data
    do
      call InputReadFlotranString(input, option)
      if (InputError(input)) exit
      call InputReadInt(input, option, temp_int)
      if (.not.InputError(input)) then
        count = count + 1
        cell_ids_p(count) = temp_int
      endif
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(cell_ids_p, max_size)
      endif
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_cells = count/option%mycommsize
      remainder = count - region%num_cells*option%mycommsize
    if(option%myrank < remainder) region%num_cells = region%num_cells + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_cells, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                    MPI_SUM, option%mycomm, ierr)
    call MPI_Scan(region%num_cells, iend, ONE_INTEGER_MPI, MPIU_INTEGER, &
                   MPI_SUM, option%mycomm, ierr)

    ! Allocate memory and save the data
    region%num_cells = iend - istart
    allocate(region%cell_ids(region%num_cells))
    region%cell_ids(1:region%num_cells) = cell_ids_p(istart+1:iend)
    deallocate(cell_ids_p)

  else if(count == 2) then
    !
    ! Input data contains cell ids + face ids
    !
    input_data_type = CELL_IDS_WITH_FACE_IDS
    cell_ids_p(1) = temp_int_array(1)
    face_ids_p(1) = temp_int_array(2)
    count = 1 ! reset the counter to represent the num of rows read

    ! Read the data
    do
      call InputReadFlotranString(input, option)
      if (InputError(input)) exit
      call InputReadInt(input, option, temp_int)
      if(InputError(input)) exit
      count = count + 1
      cell_ids_p(count) = temp_int

      call InputReadInt(input,option,temp_int)
      if(InputError(input)) then
        option%io_buffer = 'ERROR while reading the region from file'
        call printErrMsg(option)
      endif
      face_ids_p(count) = temp_int
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(cell_ids_p, max_size)
        call reallocateIntArray(face_ids_p, max_size)
      endif
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_cells = count/option%mycommsize
      remainder = count - region%num_cells*option%mycommsize
    if(option%myrank < remainder) region%num_cells = region%num_cells + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_cells,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(region%num_cells,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                   MPI_SUM,option%mycomm,ierr)

    ! Allocate memory and save the data
    allocate(region%cell_ids(region%num_cells))
    allocate(region%faces(region%num_cells))
    region%cell_ids(1:region%num_cells) = cell_ids_p(istart + 1:iend)
    region%faces(1:region%num_cells) = face_ids_p(istart + 1:iend)
    deallocate(cell_ids_p)
    deallocate(face_ids_p)

  else
    !
    ! Input data contains vertices
    !
    input_data_type = VERTEX_IDS
    vert_id_0_p(1) = temp_int_array(1)
    vert_id_1_p(1) = temp_int_array(2)
    vert_id_2_p(1) = temp_int_array(3)
    vert_id_3_p(1) = temp_int_array(4)
    if(vert_id_0_p(1) == 4 ) vert_id_4_p(1) = temp_int_array(5)
    count = 1 ! reset the counter to represent the num of rows read

    ! Read the data
    do
      call InputReadFlotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,temp_int)
      if(InputError(input)) exit
      count = count + 1
      vert_id_0_p(count) = temp_int

      vert_id_4_p(count) = -999
      do ii = 1, vert_id_0_p(count)
        call InputReadInt(input,option,temp_int)
        if(InputError(input)) then
          option%io_buffer = 'ERROR while reading the region from file'
          call printErrMsg(option)
        endif

        select case(ii)
          case(1)
            vert_id_1_p(count) = temp_int
          case(2)
            vert_id_2_p(count) = temp_int
          case(3)
            vert_id_3_p(count) = temp_int
          case(4)
            vert_id_4_p(count) = temp_int
        end select

        if (count+1 > max_size) then ! resize temporary array
          call reallocateIntArray(vert_id_0_p,max_size)
          call reallocateIntArray(vert_id_1_p,max_size)
          call reallocateIntArray(vert_id_2_p,max_size)
          call reallocateIntArray(vert_id_3_p,max_size)
          call reallocateIntArray(vert_id_4_p,max_size)
        endif
      enddo
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_verts = count/option%mycommsize
      remainder = count - region%num_verts*option%mycommsize
    if(option%myrank < remainder) region%num_verts = region%num_verts + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_verts,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(region%num_verts,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                   MPI_SUM,option%mycomm,ierr)

    ! Allocate memory and save the data
    region%num_verts = iend - istart
    allocate(region%vert_ids(0:MAX_VERT_PER_FACE,1:region%num_verts))
    region%vert_ids(0,1:region%num_verts) = vert_id_0_p(istart + 1: iend)
    region%vert_ids(1,1:region%num_verts) = vert_id_1_p(istart + 1: iend)
    region%vert_ids(2,1:region%num_verts) = vert_id_2_p(istart + 1: iend)
    region%vert_ids(3,1:region%num_verts) = vert_id_3_p(istart + 1: iend)
    region%vert_ids(4,1:region%num_verts) = vert_id_4_p(istart + 1: iend)
    deallocate(vert_id_0_p)
    deallocate(vert_id_1_p)
    deallocate(vert_id_2_p)
    deallocate(vert_id_3_p)
    deallocate(vert_id_4_p)

  endif
  
#if 0  
  count = 1
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    call InputReadInt(input,option,temp_int)
    if (.not.InputError(input)) then
      count = count + 1
      temp_int_array(count) = temp_int
      write(*,*),count,temp_int
    endif
    if (count+1 > max_size) then ! resize temporary array
      call reallocateIntArray(temp_int_array,max_size)
    endif
  enddo

  if (count > 0) then
    region%num_cells = count
    allocate(region%cell_ids(count))
    region%cell_ids(1:count) = temp_int_array(1:count)
  else
    region%num_cells = 0
    nullify(region%cell_ids)
  endif
#endif
  deallocate(temp_int_array) 

  call PetscLogEventEnd(logging%event_region_read_ascii,ierr)

end subroutine RegionReadFromFileId

! ************************************************************************** !
!
! RegionGetPtrFromList: Returns a pointer to the region matching region_name
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
function RegionGetPtrFromList(region_name,region_list)

  use String_module

  implicit none
  
  type(region_type), pointer :: RegionGetPtrFromList
  character(len=MAXWORDLENGTH) :: region_name
  PetscInt :: length
  type(region_list_type) :: region_list

  type(region_type), pointer :: region
    
  nullify(RegionGetPtrFromList)
  region => region_list%first
  
  do 
    if (.not.associated(region)) exit
    length = len_trim(region_name)
    if (length == len_trim(region%name) .and. &
        StringCompare(region%name,region_name,length)) then
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
  if (associated(region%vert_ids)) deallocate(region%vert_ids)
  
  deallocate(region)
  nullify(region)

end subroutine RegionDestroy

end module Region_module
