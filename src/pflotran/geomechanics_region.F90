module Geomechanics_Region_module
 
  use Geometry_module
  
  implicit none

  private

#include "definitions.h"
 
  type, public :: gm_region_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXSTRINGLENGTH) :: filename
    type(point3d_type), pointer :: coordinates(:)
    PetscInt :: num_verts 
    PetscInt, pointer :: vertex_ids(:) 
    type(gm_region_type), pointer :: next
  end type gm_region_type
  
  type, public :: gm_region_ptr_type
    type(gm_region_type), pointer :: ptr
  end type gm_region_ptr_type
  
  type, public :: gm_region_list_type
    PetscInt :: num_regions
    type(gm_region_type), pointer :: first
    type(gm_region_type), pointer :: last
    type(gm_region_type), pointer :: array(:)
  end type gm_region_list_type

  interface GeomechRegionCreate
    module procedure GeomechRegionCreateWithList
    module procedure GeomechRegionCreateWithNothing
    module procedure GeomechRegionCreateWithGeomechRegion
  end interface GeomechRegionCreate
  
  interface GeomechRegionReadFromFile
    module procedure GeomechRegionReadFromFileId
    module procedure GeomechRegionReadFromFilename
  end interface GeomechRegionReadFromFile
  
   public :: GeomechRegionCreate, &
             GeomechRegionDestroy, &
             GeomechRegionAddToList, &
             GeomechRegionReadFromFile, &
             GeomechRegionDestroyList, &
             GeomechRegionRead, &
             GeomechRegionInitList, &
             GeomechRegionGetPtrFromList
             
 contains
  
! ************************************************************************** !
!
! GeomechRegionCreateWithNothing: Creates a region with no arguments for 
!                                 geomechanics
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
function GeomechRegionCreateWithNothing()

  implicit none
  
  type(gm_region_type), pointer        :: GeomechRegionCreateWithNothing
  
  type(gm_region_type), pointer        :: region
  
  allocate(region)
  region%id = 0
  region%name = ""
  region%filename = ""
  region%num_verts = 0
  nullify(region%coordinates)
  nullify(region%vertex_ids)
  nullify(region%next)
  
  GeomechRegionCreateWithNothing => region

end function GeomechRegionCreateWithNothing

! ************************************************************************** !
!
! GeomechRegionCreateWithList: Creates a region from a list of vertices
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
function GeomechRegionCreateWithList(list)

  implicit none
  
  PetscInt :: list(:)
  
  type(gm_region_type), pointer     :: GeomechRegionCreateWithList
  
  type(gm_region_type), pointer     :: region

  region => GeomechRegionCreateWithNothing()
  region%num_verts = size(list)
  allocate(region%vertex_ids(region%num_verts))
  region%vertex_ids = list
  
  GeomechRegionCreateWithList => region

end function GeomechRegionCreateWithList

! ************************************************************************** !
!
! GeomechRegionCreateWithGeomechRegion: Creates a copy of a region
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
function GeomechRegionCreateWithGeomechRegion(region)

  use Unstructured_Cell_module

  implicit none
  
  type(gm_region_type), pointer     :: GeomechRegionCreateWithGeomechRegion
  type(gm_region_type), pointer     :: region
  
  type(gm_region_type), pointer     :: new_region
  PetscInt                          :: icount, temp_int
  
  new_region => GeomechRegionCreateWithNothing()
  
  new_region%id = region%id
  new_region%name = region%name
  new_region%filename = region%filename
  new_region%num_verts = region%num_verts
  if (associated(region%coordinates)) then
    call GeometryCopyCoordinates(region%coordinates, &
                                 new_region%coordinates)
  endif
  if (associated(region%vertex_ids)) then
    allocate(new_region%vertex_ids(new_region%num_verts))
    new_region%vertex_ids(1:new_region%num_verts) = &
    region%vertex_ids(1:new_region%num_verts)
  endif
 
  GeomechRegionCreateWithGeomechRegion => new_region
  
end function GeomechRegionCreateWithGeomechRegion

! ************************************************************************** !
!
! GeomechRegionInitList: Initializes a region list
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechRegionInitList(list)

  implicit none

  type(gm_region_list_type)     :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_regions = 0

end subroutine GeomechRegionInitList

! ************************************************************************** !
!
! GeomechRegionAddToList: Adds a new region to a region list
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechRegionAddToList(new_region,list)

  implicit none
  
  type(gm_region_type), pointer     :: new_region
  type(gm_region_list_type)         :: list
  
  list%num_regions = list%num_regions + 1
  new_region%id = list%num_regions
  if (.not.associated(list%first)) list%first => new_region
  if (associated(list%last)) list%last%next => new_region
  list%last => new_region
  
end subroutine GeomechRegionAddToList

! ************************************************************************** !
!
! GeomechRegionRead: Reads a region from the input file
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechRegionRead(region,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type)             :: option
  type(gm_region_type)          :: region
  type(input_type)              :: input
  
  character(len=MAXWORDLENGTH)  :: keyword, word
 
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS_REGION')
    call StringToUpper(keyword)   

    select case(trim(keyword))
      case('COORDINATE')
        allocate(region%coordinates(1))
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x) 
        if (InputError(input)) then
          input%ierr = 0
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'GEOMECHANICS_REGION')
          call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x)
        endif
        call InputErrorMsg(input,option,'x-coordinate','GEOMECHANICS_REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%y)
        call InputErrorMsg(input,option,'y-coordinate','GEOMECHANICS_REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%z)
        call InputErrorMsg(input,option,'z-coordinate','GEOMECHANICS_REGION')
      case('COORDINATES')
        call GeometryReadCoordinates(input,option,region%name, &
                                     region%coordinates)
      case('FILE')
        call InputReadNChars(input,option,region%filename,MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','GEOMECHANICS_REGION')
        call GeomechRegionReadFromFilename(region,option,region%filename)
      case('LIST')
        option%io_buffer = 'GEOMECHANICS_REGION LIST currently not implemented'
        call printErrMsg(option)
      case default
        option%io_buffer = 'REGION keyword: '//trim(keyword)//' not recognized'
        call printErrMsg(option)
    end select
  enddo
 
end subroutine GeomechRegionRead

! ************************************************************************** !
!
! GeomechRegionReadFromFilename: Reads a list of vertex ids from a file named 
!                                filename
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechRegionReadFromFilename(region,option,filename)

  use Input_module
  use Option_module
  use Utility_module
  
  implicit none
  
  type(gm_region_type)               :: region
  type(option_type)                  :: option
  type(input_type), pointer          :: input
  character(len=MAXSTRINGLENGTH)     :: filename
  
  input => InputCreate(IUNIT_TEMP,filename,option)
  call GeomechRegionReadFromFileId(region,input,option)          
  call InputDestroy(input)         

end subroutine GeomechRegionReadFromFilename

! ************************************************************************** !
!
! GeomechRegionReadFromFileId: Reads a list of vertex ids from an open file
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechRegionReadFromFileId(region,input,option)

  use Input_module
  use Option_module
  use Utility_module
  use Logging_module
  use Unstructured_Cell_module
  
  implicit none
  
  type(gm_region_type)              :: region
  type(option_type)                 :: option
  type(input_type)                  :: input
  
  PetscBool                         :: continuation_flag
  character(len=MAXWORDLENGTH)      :: word
  character(len=1)                  :: backslash
  character(len=MAXSTRINGLENGTH)    :: string, string1

  PetscInt, pointer :: temp_int_array(:)
  PetscInt, pointer :: vertex_ids(:)
  PetscInt :: max_size
  PetscInt :: count
  PetscInt :: temp_int
  PetscInt :: input_data_type
  PetscInt :: ii
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: remainder
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_region_read_ascii,ierr)
    
  max_size = 1000
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  allocate(temp_int_array(max_size))
  allocate(vertex_ids(max_size))
  
  temp_int_array = 0
  vertex_ids = 0
  
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
    vertex_ids(1) = temp_int_array(1)
    count = 1

    ! Read the data
    do
      call InputReadFlotranString(input, option)
      if (InputError(input)) exit
      call InputReadInt(input, option, temp_int)
      if (.not.InputError(input)) then
        count = count + 1
        vertex_ids(count) = temp_int
      endif
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(vertex_ids, max_size)
      endif
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_verts = count/option%mycommsize
      remainder = count - region%num_verts*option%mycommsize
    if (option%myrank < remainder) region%num_verts = region%num_verts + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_verts, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                    MPI_SUM, option%mycomm, ierr)
    call MPI_Scan(region%num_verts, iend, ONE_INTEGER_MPI, MPIU_INTEGER, &
                   MPI_SUM, option%mycomm, ierr)

    ! Allocate memory and save the data
    region%num_verts = iend - istart
    allocate(region%vertex_ids(region%num_verts))
    region%vertex_ids(1:region%num_verts) = vertex_ids(istart+1:iend)
    deallocate(vertex_ids)
  else
   option%io_buffer = 'Provide one vertex_id per line, GEOMECHANICS_REGION.'
   call printErrMsg(option) 
  endif
  
  deallocate(temp_int_array)
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  write(string1,*) region%name
  string = 'geomech_region_' // trim(adjustl(string1)) // '_vertex_ids' &
    // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ii = 1, region%num_verts
    write(86,'(i5)') region%vertex_ids(ii)
  enddo  
  close(86)
#endif    
    

  call PetscLogEventEnd(logging%event_region_read_ascii,ierr)

end subroutine GeomechRegionReadFromFileId

! ************************************************************************** !
!
! GeomechRegionGetPtrFromList: Returns a pointer to the region matching region_name
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
function GeomechRegionGetPtrFromList(region_name,region_list)

  use String_module

  implicit none
  
  type(gm_region_type), pointer          :: GeomechRegionGetPtrFromList
  character(len=MAXWORDLENGTH)           :: region_name
  PetscInt :: length
  type(gm_region_list_type)              :: region_list

  type(gm_region_type), pointer          :: region
    
  nullify(GeomechRegionGetPtrFromList)
  region => region_list%first
  
  do 
    if (.not.associated(region)) exit
    length = len_trim(region_name)
    if (length == len_trim(region%name) .and. &
        StringCompare(region%name,region_name,length)) then
      GeomechRegionGetPtrFromList => region
      return
    endif
    region => region%next
  enddo
  
end function GeomechRegionGetPtrFromList

! ************************************************************************** !
!
! GeomechRegionDestroyList: Deallocates a list of regions
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechRegionDestroyList(region_list)

  implicit none
  
  type(gm_region_list_type), pointer        :: region_list
  
  type(gm_region_type), pointer             :: region, prev_region
  
  if (.not.associated(region_list)) return
  
  region => region_list%first
  do 
    if (.not.associated(region)) exit
    prev_region => region
    region => region%next
    call GeomechRegionDestroy(prev_region)
  enddo
  
  region_list%num_regions = 0
  nullify(region_list%first)
  nullify(region_list%last)
  if (associated(region_list%array)) deallocate(region_list%array)
  nullify(region_list%array)
  
  deallocate(region_list)
  nullify(region_list)

end subroutine GeomechRegionDestroyList

! ************************************************************************** !
!
! GeomechRegionDestroy: Deallocates a region
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechRegionDestroy(region)

  implicit none
  
  type(gm_region_type), pointer :: region
  
  if (.not.associated(region)) return
  
  if (associated(region%vertex_ids)) deallocate(region%vertex_ids)
  nullify(region%vertex_ids)
  if (associated(region%coordinates)) deallocate(region%coordinates)
  nullify(region%coordinates)
  
  nullify(region%next)

  deallocate(region)
  nullify(region)
  
end subroutine GeomechRegionDestroy

end module Geomechanics_Region_module
