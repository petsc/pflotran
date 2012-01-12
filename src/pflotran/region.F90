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
    !TODO(geh): Tear anything to do with structured/unstructured grids other
    !           than cell id ane face id out of region.
    PetscInt, pointer :: vertex_ids(:,:) ! For Unstructured mesh
    PetscInt :: num_verts              ! For Unstructured mesh
    PetscInt :: grid_type  ! To identify whether region is applicable to a Structured or Unstructred mesh
    type(region_sideset_type), pointer :: sideset
    type(region_polygonal_volume_type), pointer :: polygonal_volume
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

  type, public :: region_sideset_type
    PetscInt :: nfaces
    PetscInt, pointer :: face_vertices(:,:)
  end type region_sideset_type
  
  type, public :: region_polygonal_volume_type
    type(point3d_type), pointer :: xy_coordinates(:)
    type(point3d_type), pointer :: xz_coordinates(:)
    type(point3d_type), pointer :: yz_coordinates(:)
  end type region_polygonal_volume_type

  interface RegionPointInPolygon
    module procedure RegionPointInPolygon1
    module procedure RegionPointInPolygon2
  end interface RegionPointInPolygon
    
  interface RegionCreate
    module procedure RegionCreateWithBlock
    module procedure RegionCreateWithList
    module procedure RegionCreateWithNothing
    module procedure RegionCreateWithRegion    
  end interface RegionCreate
  
  interface RegionReadFromFile
    module procedure RegionReadFromFileId
    module procedure RegionReadFromFilename
    module procedure RegionReadSideset
  end interface RegionReadFromFile
  
  public :: RegionCreate, RegionDestroy, RegionAddToList, RegionReadFromFile, &
            RegionInitList, RegionDestroyList, RegionGetPtrFromList, & 
            RegionRead, RegionReadSideset, RegionCreateSideset
  
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
  nullify(region%vertex_ids)
  nullify(region%sideset)
  nullify(region%polygonal_volume)
  nullify(region%next)
  
  RegionCreateWithNothing => region

end function RegionCreateWithNothing

! ************************************************************************** !
!
! RegionCreateSideset: Creates a sideset
! author: Glenn Hammond
! date: 12/19/11
!
! ************************************************************************** !
function RegionCreateSideset()

  implicit none
  
  type(region_sideset_type), pointer :: RegionCreateSideset
  
  type(region_sideset_type), pointer :: sideset
  
  allocate(sideset)
  sideset%nfaces = 0
  nullify(sideset%face_vertices)
  
  RegionCreateSideset => sideset

end function RegionCreateSideset

! ************************************************************************** !
!
! RegionCreatePolygonalVolume: Creates a polygonal volume.  I.e. a volume
!                              defined by 3 polygons, on each plane of the
!                              principle coordinate system x,y,z
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function RegionCreatePolygonalVolume()

  implicit none
  
  type(region_polygonal_volume_type), pointer :: RegionCreatePolygonalVolume
  
  type(region_polygonal_volume_type), pointer :: polygonal_volume
  
  allocate(polygonal_volume)
  nullify(polygonal_volume%xy_coordinates)
  nullify(polygonal_volume%xz_coordinates)
  nullify(polygonal_volume%yz_coordinates)
  
  RegionCreatePolygonalVolume => polygonal_volume

end function RegionCreatePolygonalVolume

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
  PetscInt :: icount, temp_int
  
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
    call RegionCopyCoordinates(region%coordinates, &
                                new_region%coordinates)
  endif
  if (associated(region%cell_ids)) then
    allocate(new_region%cell_ids(new_region%num_cells))
    new_region%cell_ids(1:new_region%num_cells) = &
      region%cell_ids(1:region%num_cells)
  endif
  if (associated(region%faces)) then
    allocate(new_region%faces(new_region%num_cells))
    new_region%faces(1:new_region%num_cells) = &
      region%faces(1:region%num_cells)
  endif
  if (associated(region%vertex_ids)) then
    allocate(new_region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts))
    new_region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts) = &
    region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts)
  endif
  if (associated(region%sideset)) then
    new_region%sideset => RegionCreateSideSet()
    new_region%sideset%nfaces = region%sideset%nfaces
    allocate(new_region%sideset%face_vertices( &
               size(region%sideset%face_vertices,1), &
               size(region%sideset%face_vertices,2)))
    new_region%sideset%face_vertices = region%sideset%face_vertices
  endif
  if (associated(region%polygonal_volume)) then
    new_region%polygonal_volume => RegionCreatePolygonalVolume()
    if (associated(region%polygonal_volume%xy_coordinates)) then
      call RegionCopyCoordinates(region%polygonal_volume%xy_coordinates, &
                                 new_region%polygonal_volume%xy_coordinates)
    endif
    if (associated(region%polygonal_volume%xz_coordinates)) then
      call RegionCopyCoordinates(region%polygonal_volume%xz_coordinates, &
                                 new_region%polygonal_volume%xz_coordinates)
    endif
    if (associated(region%polygonal_volume%yz_coordinates)) then
      call RegionCopyCoordinates(region%polygonal_volume%yz_coordinates, &
                                 new_region%polygonal_volume%yz_coordinates)
    endif
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
        call RegionReadCoordinates(input,option,region%name,region%coordinates)
      case('POLYGON')
        if (.not.associated(region%polygonal_volume)) then
          region%polygonal_volume => RegionCreatePolygonalVolume()
        endif
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'plane','REGION')
        call StringToUpper(word)
        select case(word)
          case('XY')
            call RegionReadCoordinates(input,option,region%name, &
                                       region%polygonal_volume%xy_coordinates)
          case('XZ')
            call RegionReadCoordinates(input,option,region%name, &
                                       region%polygonal_volume%xz_coordinates)
          case('YZ')
            call RegionReadCoordinates(input,option,region%name, &
                                       region%polygonal_volume%yz_coordinates)
          case default
            option%io_buffer = 'PLANE not recognized for REGION POLYGON.  ' // &
              'Use either XY, XZ or YZ.'
            call printErrMsg(option)
        end select
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
! RegionReadCoordinates: Reads a list of coordinates
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine RegionReadCoordinates(input,option,region_name,coordinates)

  use Input_module
  use Option_module

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: region_name
  type(point3d_type), pointer :: coordinates(:)
  
  PetscInt :: icount
  PetscInt, parameter :: max_num_coordinates = 100
  type(point3d_type) :: temp_coordinates(max_num_coordinates)

  icount = 0
  do
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'REGION')
    if (InputCheckExit(input,option)) exit              
    icount = icount + 1
    if (icount > max_num_coordinates) then
      write(option%io_buffer, &
            '(''Number of coordinates in region '',a, &
          &'' exceeds limit of '',i3)') region_name, max_num_coordinates
      option%io_buffer = trim(option%io_buffer) // &
        ' Increase size of PetscInt, parameter :: max_num_coordinates ' // &
        ' in RegionReadCoordinates()'
      call printErrMsg(option)
    endif
    call InputReadDouble(input,option,temp_coordinates(icount)%x) 
    call InputErrorMsg(input,option,'x-coordinate','REGION')
    call InputReadDouble(input,option,temp_coordinates(icount)%y)
    call InputErrorMsg(input,option,'y-coordinate','REGION')
    call InputReadDouble(input,option,temp_coordinates(icount)%z)
    call InputErrorMsg(input,option,'z-coordinate','REGION')
  enddo
  allocate(coordinates(icount))
  do icount = 1, size(coordinates)
    coordinates(icount)%x = temp_coordinates(icount)%x
    coordinates(icount)%y = temp_coordinates(icount)%y
    coordinates(icount)%z = temp_coordinates(icount)%z
  enddo

end subroutine RegionReadCoordinates

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
  
  !TODO(geh): clean and optimize this subroutine
  
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
    if (InputError(input)) exit
    count = count + 1
    temp_int_array(count) = temp_int
  enddo

  if (count == 1) then
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
    if (option%myrank < remainder) region%num_cells = region%num_cells + 1
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

  else if (count == 2) then
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
      if (InputError(input)) exit
      count = count + 1
      cell_ids_p(count) = temp_int

      call InputReadInt(input,option,temp_int)
      if (InputError(input)) then
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
    if (option%myrank < remainder) region%num_cells = region%num_cells + 1
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
    if (vert_id_0_p(1) == 4 ) vert_id_4_p(1) = temp_int_array(5)
    count = 1 ! reset the counter to represent the num of rows read

    ! Read the data
    do
      call InputReadFlotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,temp_int)
      if (InputError(input)) exit
      count = count + 1
      vert_id_0_p(count) = temp_int

      vert_id_4_p(count) = -999
      do ii = 1, vert_id_0_p(count)
        call InputReadInt(input,option,temp_int)
        if (InputError(input)) then
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
    if (option%myrank < remainder) region%num_verts = region%num_verts + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_verts,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(region%num_verts,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                   MPI_SUM,option%mycomm,ierr)

    ! Allocate memory and save the data
    region%num_verts = iend - istart
    allocate(region%vertex_ids(0:MAX_VERT_PER_FACE,1:region%num_verts))
    region%vertex_ids(0,1:region%num_verts) = vert_id_0_p(istart + 1: iend)
    region%vertex_ids(1,1:region%num_verts) = vert_id_1_p(istart + 1: iend)
    region%vertex_ids(2,1:region%num_verts) = vert_id_2_p(istart + 1: iend)
    region%vertex_ids(3,1:region%num_verts) = vert_id_3_p(istart + 1: iend)
    region%vertex_ids(4,1:region%num_verts) = vert_id_4_p(istart + 1: iend)
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
! RegionReadSideSet: Reads an unstructured grid sideset
! author: Glenn Hammond
! date: 12/19/11
!
! ************************************************************************** !
subroutine RegionReadSideSet(sideset,filename,option)

  use Input_module
  use Option_module
  use String_module
  
  implicit none
  
  type(region_sideset_type) :: sideset
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card, word
  PetscInt :: num_faces_local_save
  PetscInt :: num_faces_local
  PetscInt :: num_to_read
  PetscInt, parameter :: max_nvert_per_face = 4
  PetscInt, allocatable :: temp_int_array(:,:)

  PetscInt :: iface, ivertex, irank, num_vertices
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename)

! Format of sideset file
! type: T=triangle, Q=quadrilateral
! vertn(Q) = 4
! vertn(T) = 3
! -----------------------------------------------------------------
! num_faces  (integer)
! type vert1 vert2 ... vertn  ! for face 1 (integers)
! type vert1 vert2 ... vertn  ! for face 2
! ...
! ...
! type vert1 vert2 ... vertn  ! for face num_faces
! -----------------------------------------------------------------

  card = 'Unstructured Sideset'

  call InputReadFlotranString(input,option)
  string = 'unstructured sideset'
  call InputReadStringErrorMsg(input,option,card)  

  ! read num_faces
  call InputReadInt(input,option,sideset%nfaces)
  call InputErrorMsg(input,option,'number of faces',card)

#ifndef PARALLEL_SIDESET
  num_to_read = sideset%nfaces
  allocate(sideset%face_vertices(max_nvert_per_face,num_to_read))
  sideset%face_vertices = -999
#else
  ! divide faces across ranks
  num_faces_local = sideset%nfaces/option%mycommsize 
  num_faces_local_save = num_faces_local
  remainder = sideset%nfaces - num_faces_local*option%mycommsize
  if (option%myrank < remainder) num_faces_local = &
                                 num_faces_local + 1

  ! allocate array to store vertices for each faces
  allocate(sideset%face_vertices(max_nvert_per_face, &
                                 num_faces_local))
  sideset%face_vertices = -999

  ! for now, read all faces from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(max_nvert_per_face, &
                            num_faces_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_int_array = -999
      num_to_read = num_faces_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
#endif      
      do iface = 1, num_to_read
        ! read in the vertices defining the cell face
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,card)  
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'face type',card)
        call StringToUpper(word)
        select case(word)
          case('Q')
            num_vertices = 4
          case('T')
            num_vertices = 3
        end select
        do ivertex = 1, num_vertices
#ifdef PARALLEL_SIDESET
          call InputReadInt(input,option,temp_int_array(ivertex,iface))
#else
          call InputReadInt(input,option,sideset%face_vertices(ivertex,iface))
#endif
          call InputErrorMsg(input,option,'vertex id',card)
        enddo
      enddo
#ifdef PARALLEL_SIDESET      
      ! if the faces reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_faces_local
        string = trim(adjustl(string)) // ' faces stored on p0'
        print *, trim(string)
#endif
        sideset%face_vertices(:,1:num_faces_local) = &
          temp_int_array(:,1:num_faces_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' faces sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*max_nvert_per_face
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_faces_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' faces received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    int_mpi = num_faces_local*max_nvert_per_face
    call MPI_Recv(sideset%face_vertices,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif

!  unstructured_grid%nlmax = num_faces_local
!  unstructured_grid%num_vertices_local = num_vertices_local
#endif

  call InputDestroy(input)

end subroutine RegionReadSideSet

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
! RegionDestroySideset: Deallocates a unstructured grid side set
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine RegionDestroySideset(sideset)

  implicit none
  
  type(region_sideset_type), pointer :: sideset
  
  if (.not.associated(sideset)) return
  
  if (associated(sideset%face_vertices)) deallocate(sideset%face_vertices)
  nullify(sideset%face_vertices)
  
  deallocate(sideset)
  nullify(sideset)
  
end subroutine RegionDestroySideset

! ************************************************************************** !
!
! RegionPointInPolygonalVolume: Determines whether a point in xyz space is 
!                               within a polygonal volume defined by polygons
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function RegionPointInPolygonalVolume(x,y,z,region,option)
 
  use Option_module

  implicit none
  
  PetscReal :: x, y, z
  type(region_type) :: region
  type(option_type) :: option
  
  PetscBool :: xy, xz, yz
  PetscBool :: RegionPointInPolygonalVolume
    
  if (.not.associated(region%polygonal_volume)) then
    option%io_buffer = 'No polygonal volume defined in REGION: ' // &
      trim(region%name)
  endif
  
  RegionPointInPolygonalVolume = PETSC_FALSE
  
  ! if a polygon is not defined in a particular direction, it is assumed that
  ! the point is within the "infinite" polygon
  
  ! XY plane
  if (associated(region%polygonal_volume%xy_coordinates)) then
    if (.not.RegionPointInPolygon(x,y,z,Z_DIRECTION, &
                                  region%polygonal_volume%xy_coordinates)) then
      return
    endif
  endif

  ! XZ plane
  if (associated(region%polygonal_volume%xz_coordinates)) then
    if (.not.RegionPointInPolygon(x,y,z,Y_DIRECTION, &
                                  region%polygonal_volume%xz_coordinates)) then
      return
    endif
  endif

  ! YZ plane
  if (associated(region%polygonal_volume%yz_coordinates)) then
    if (.not.RegionPointInPolygon(x,y,z,X_DIRECTION, &
                                  region%polygonal_volume%yz_coordinates)) then
      return
    endif
  endif
  
  ! if point is not within any of polygons above, the function will return
  ! prior to this point.
  RegionPointInPolygonalVolume = PETSC_TRUE
  
end function RegionPointInPolygonalVolume

! ************************************************************************** !
!
! RegionPointInPolygon1: Determines whether a point in xyz space is within
!                        a polygon based on coordinate object
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function RegionPointInPolygon1(x,y,z,axis,coordinates)

  implicit none
  
  PetscReal :: x, y, z
  PetscInt :: axis
  type(point3d_type) :: coordinates(:)
  
  PetscInt, parameter :: max_num_coordinates = 100
  PetscReal :: xx, yy
  PetscInt :: i, num_coordinates
  PetscReal :: xx_array(max_num_coordinates), yy_array(max_num_coordinates)
  
  PetscBool :: RegionPointInPolygon1
  
  RegionPointInPolygon1 = PETSC_FALSE
  
  num_coordinates = size(coordinates)
  select case(axis)
    case(Z_DIRECTION)
      xx = x
      yy = y
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%x
        yy_array(i) = coordinates(i)%y
      enddo
    case(Y_DIRECTION)
      xx = x
      yy = z
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%x
        yy_array(i) = coordinates(i)%z
      enddo
    case(X_DIRECTION)
      xx = y
      yy = z
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%y
        yy_array(i) = coordinates(i)%z
      enddo
  end select
  RegionPointInPolygon1 = RegionPointInPolygon2(xx,yy,xx_array,yy_array)
  
end function RegionPointInPolygon1

! ************************************************************************** !
!
! RegionPointInPolygon2: Determines whether a point in xyz space is within
!                       a polygon
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function RegionPointInPolygon2(x,y,x_array,y_array)

  implicit none
  
  PetscReal :: x, y
  PetscReal :: x_array(:), y_array(:)
  
  PetscBool :: RegionPointInPolygon2
  PetscInt :: num_coordinates
  PetscInt :: i, j
  
  RegionPointInPolygon2 = PETSC_FALSE

  num_coordinates = size(x_array)
  j = 1
  do i = 1, num_coordinates
    j = i + 1
    if (j == num_coordinates) j = 1
    if ((y_array(i) < y .and. y_array(j) >= y) .or. &
        (y_array(j) < y .and. y_array(i) >= y)) then
      if ((x_array(i) + &
           (y-y_array(i))/(y_array(j)-y_array(i))*(x_array(j)-x_array(i))) &
          < x) then
        RegionPointInPolygon2 = .not.RegionPointInPolygon2
      endif
    endif    
  enddo

end function RegionPointInPolygon2

! ************************************************************************** !
!
! RegionCopyCoordinates: Deallocates a polygonal volume object
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine RegionCopyCoordinates(coordinates_in,coordinates_out)

  implicit none
  
  type(point3d_type) :: coordinates_in(:)
  type(point3d_type), pointer :: coordinates_out(:)
  
  PetscInt :: num_coordinates
  PetscInt :: i
  
  num_coordinates = size(coordinates_in)
  allocate(coordinates_out(num_coordinates))
  do i = 1, num_coordinates
   coordinates_out(i)%x = coordinates_in(i)%x 
   coordinates_out(i)%y = coordinates_in(i)%y
   coordinates_out(i)%z = coordinates_in(i)%z 
  enddo
 
end subroutine RegionCopyCoordinates

! ************************************************************************** !
!
! RegionDestroyPolygonalVolume: Deallocates a polygonal volume object
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine RegionDestroyPolygonalVolume(polygonal_volume)

  implicit none
  
  type(region_polygonal_volume_type), pointer :: polygonal_volume
  
  if (.not.associated(polygonal_volume)) return
  
  if (associated(polygonal_volume%xy_coordinates)) &
    deallocate(polygonal_volume%xy_coordinates)
  nullify(polygonal_volume%xy_coordinates)
  if (associated(polygonal_volume%xz_coordinates)) &
    deallocate(polygonal_volume%xz_coordinates)
  nullify(polygonal_volume%xz_coordinates)
  if (associated(polygonal_volume%yz_coordinates)) &
    deallocate(polygonal_volume%yz_coordinates)
  nullify(polygonal_volume%yz_coordinates)
  
  deallocate(polygonal_volume)
  nullify(polygonal_volume)
  
end subroutine RegionDestroyPolygonalVolume

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
  if (associated(region%coordinates)) deallocate(region%coordinates)
  nullify(region%coordinates)
  call RegionDestroySideset(region%sideset)
  call RegionDestroyPolygonalVolume(region%polygonal_volume)
  
  if (associated(region%vertex_ids)) deallocate(region%vertex_ids)
  nullify(region%vertex_ids)
  
  nullify(region%next)

  deallocate(region)
  nullify(region)

end subroutine RegionDestroy

end module Region_module
