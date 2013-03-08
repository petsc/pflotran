module Dataset_Aux_module
 
  implicit none

  private

#include "definitions.h"

  type, public :: dataset_map_type
    character(len=MAXSTRINGLENGTH) :: h5_dataset_map_name
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt,pointer :: map(:,:)
    PetscInt         :: map_dims_global(2)
    PetscInt         :: map_dims_local(2)
    PetscInt,pointer :: datatocell_ids(:)
    PetscInt,pointer :: cell_ids_local(:)
    PetscBool        :: first_time
  end type dataset_map_type
 
  type, public :: dataset_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: h5_dataset_name
    character(len=MAXSTRINGLENGTH) :: filename
    PetscBool :: is_transient
    PetscBool :: is_cell_indexed
    PetscInt :: max_buffer_size
    PetscInt :: data_type
    PetscInt :: data_dim ! dimensions of data: XY, X, YXZ, etc.
    PetscBool :: realization_dependent
    PetscBool :: cell_centered
    ! all data stored internally as a 1D array
    PetscInt, pointer :: iarray(:)
    PetscReal, pointer :: rarray(:)
    PetscInt :: array_size
    PetscReal :: rmax ! maximum rarray value in dataset
    PetscReal :: rmin ! maximum rarray value in dataset
    PetscInt :: ndims ! excludes the time dimension
    PetscInt, pointer :: dims(:) ! excludes the time dimension
    PetscReal, pointer :: origin(:)
    PetscReal, pointer :: discretization(:)
    type(dataset_buffer_type), pointer :: buffer
    type(dataset_type), pointer :: next
    type(dataset_map_type), pointer :: dataset_map
  end type dataset_type
  
  ! the dataset buffer is designed to hold multiple times
  type, public :: dataset_buffer_type
    PetscInt :: num_times_total
    PetscInt :: num_times_in_buffer
    PetscReal, pointer :: time_array(:)
    PetscReal :: cur_time
    PetscInt :: cur_time_index  
    PetscInt :: time_offset     ! zero-based offset of buffer start from first time
    PetscInt :: array_size      ! size of iarray/rarrays
    PetscInt, pointer :: iarray(:)
    PetscReal, pointer :: rarray(:)
  end type dataset_buffer_type

  ! dataset types
  PetscInt, parameter :: DATASET_INTEGER = 1
  PetscInt, parameter :: DATASET_REAL = 2
  
  PetscInt, parameter, public :: DIM_NULL = 0
  PetscInt, parameter, public :: DIM_X = 1
  PetscInt, parameter, public :: DIM_Y = 2
  PetscInt, parameter, public :: DIM_Z = 3
  PetscInt, parameter, public :: DIM_XY = 4
  PetscInt, parameter, public :: DIM_XZ = 5
  PetscInt, parameter, public :: DIM_YZ = 6
  PetscInt, parameter, public :: DIM_XYZ = 7
  PetscInt, parameter, public :: DIM_CELL = 8
    
  public :: DatasetCreate, &
            DatasetBufferCreate, &
            DatasetRead, &
            DatasetAddToList, &
            DatasetGetPointer, &
            DatasetInterpolateReal, &
            DatasetInterpolateBetweenTimes, &
            DatasetSetDimension, &
            DatasetGetNDimensions, &
            DatasetReorder, &
            DatasetPrint, &
            DatasetGetTimes, &
            DatasetDestroy

contains

! ************************************************************************** !
!
! DatasetCreate: Creates a dataset object
! author: Glenn Hammond
! date: 01/12/11, 10/25,11
!
! ************************************************************************** !
function DatasetCreate()
  
  implicit none

  type(dataset_type), pointer :: DatasetCreate
  
  type(dataset_type), pointer :: dataset
  
  allocate(dataset)
  dataset%name = ''
  dataset%h5_dataset_name = ''
  dataset%filename = ''
  dataset%realization_dependent = PETSC_FALSE
  dataset%is_transient = PETSC_FALSE
  dataset%is_cell_indexed = PETSC_FALSE
  dataset%data_type = 0
  dataset%max_buffer_size = 0

  dataset%cell_centered = PETSC_FALSE
  dataset%ndims = 0 ! ndims should not include time dimension
  dataset%rmax = -1.d20
  dataset%rmin = 1.d20
  dataset%data_dim = DIM_NULL
  dataset%array_size = 0
  nullify(dataset%iarray)
  nullify(dataset%rarray)
  nullify(dataset%dims)
  nullify(dataset%discretization)
  nullify(dataset%origin)
  nullify(dataset%buffer)
  nullify(dataset%next)
  nullify(dataset%dataset_map)

  DatasetCreate => dataset

end function DatasetCreate

! ************************************************************************** !
!
! DatasetBufferCreate: Creates a dataset buffer object
! author: Glenn Hammond
! date: 10/28/11
!
! ************************************************************************** !
function DatasetBufferCreate()
  
  implicit none

  type(dataset_buffer_type), pointer :: DatasetBufferCreate
  
  type(dataset_buffer_type), pointer :: dataset_buffer

  allocate(dataset_buffer)
  dataset_buffer%num_times_total = 0
  dataset_buffer%num_times_in_buffer = 0
  dataset_buffer%cur_time_index = 0
  dataset_buffer%time_offset = 0
  dataset_buffer%array_size = 0
  dataset_buffer%cur_time = 0.d0
  nullify(dataset_buffer%time_array)
  nullify(dataset_buffer%iarray)
  nullify(dataset_buffer%rarray)
  
  DatasetBufferCreate => dataset_buffer

end function DatasetBufferCreate

! ************************************************************************** !
!
! DatasetRead: Reads in contents of a dataset card
! author: Glenn Hammond
! date: 01/12/11
! 
! ************************************************************************** !
subroutine DatasetRead(dataset,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  type(dataset_type) :: dataset
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','DATASET')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('NAME') 
        call InputReadWord(input,option,dataset%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','DATASET')
      case('HDF5_DATASET_NAME') 
        call InputReadWord(input,option,dataset%h5_dataset_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'hdf5_dataset_name','DATASET')
      case('FILENAME') 
        call InputReadNChars(input,option,dataset%filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','DATASET')
      case('REALIZATION_DEPENDENT')
        dataset%realization_dependent = PETSC_TRUE
      case('MAX_BUFFER_SIZE') 
        call InputReadInt(input,option,dataset%max_buffer_size)
        call InputErrorMsg(input,option,'max_buffer_size','DATASET')
      case('MAP_HDF5_DATASET_NAME')
        if(.not.associated(dataset%dataset_map)) then
          allocate(dataset%dataset_map)
          dataset%dataset_map%h5_dataset_map_name=''
          dataset%dataset_map%filename=''
          dataset%dataset_map%first_time=PETSC_TRUE
          nullify(dataset%dataset_map%map)
          nullify(dataset%dataset_map%datatocell_ids)
          nullify(dataset%dataset_map%cell_ids_local)
        endif
        call InputReadWord(input,option,dataset%dataset_map%h5_dataset_map_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'MAP_HDF5_DATASET_NAME','DATASET')
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in dataset'    
        call printErrMsg(option)
    end select  
  
  enddo
  
  if (len_trim(dataset%h5_dataset_name) < 1) then
    dataset%h5_dataset_name = dataset%name
  endif

  if (associated(dataset%dataset_map)) then
    if(len_trim(dataset%dataset_map%filename)<1) then
      dataset%dataset_map%filename=dataset%filename
    endif
  endif

end subroutine DatasetRead

! ************************************************************************** !
!
! DatasetAddToList: Adds a dataset to linked list
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine DatasetAddToList(dataset,list)

  implicit none
  
  type(dataset_type), pointer :: dataset
  type(dataset_type), pointer :: list

  type(dataset_type), pointer :: cur_dataset
  
  if (associated(list)) then
    cur_dataset => list
    ! loop to end of list
    do
      if (.not.associated(cur_dataset%next)) exit
      cur_dataset => cur_dataset%next
    enddo
    cur_dataset%next => dataset
  else
    list => dataset
  endif
  
end subroutine DatasetAddToList

! ************************************************************************** !
!
! DatasetGetPointer: Returns the pointer to the dataset named "name"
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
function DatasetGetPointer(dataset_list, dataset_name, debug_string, option)

  use Option_module
  use String_module
  
  type(dataset_type), pointer :: dataset_list
  character(len=MAXWORDLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: debug_string
  type(option_type) :: option

  type(dataset_type), pointer :: DatasetGetPointer
  PetscBool :: found
  type(dataset_type), pointer :: cur_dataset

  found = PETSC_FALSE
  cur_dataset => dataset_list
  do 
    if (.not.associated(cur_dataset)) exit
    if (StringCompare(dataset_name, &
                      cur_dataset%name,MAXWORDLENGTH)) then
      found = PETSC_TRUE
      DatasetGetPointer => cur_dataset
      return
    endif
    cur_dataset => cur_dataset%next
  enddo
  if (.not.found) then
    option%io_buffer = 'Dataset "' // trim(dataset_name) // '" in "' // &
             trim(debug_string) // '" not found among available datasets.'
    call printErrMsgByRank(option)    
  endif

end function DatasetGetPointer

! ************************************************************************** !
!
! DatasetSetDimension: Sets the dimension of the dataset
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
subroutine DatasetSetDimension(dataset,word)

  use String_module

  implicit none
  
  type(dataset_type) :: dataset
  character(len=MAXWORDLENGTH) :: word

  call StringToUpper(word)
  select case(word)
    case('X')
      dataset%data_dim = DIM_X
    case('Y')
      dataset%data_dim = DIM_Y
    case('Z')
      dataset%data_dim = DIM_Z
    case('XY')
      dataset%data_dim = DIM_XY
    case('XZ')
      dataset%data_dim = DIM_XZ
    case('YZ')
      dataset%data_dim = DIM_YZ
    case('XYZ')
      dataset%data_dim = DIM_XYZ
    case('CELL')
      dataset%data_dim = DIM_CELL
  end select
      
end subroutine DatasetSetDimension

! ************************************************************************** !
!
! DatasetGetNDimensions: Returns the number of dimensions
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
function DatasetGetNDimensions(dataset)

  implicit none
  
  type(dataset_type) :: dataset
  
  PetscInt :: DatasetGetNDimensions

  select case(dataset%data_dim)
    case(DIM_X,DIM_Y,DIM_Z,DIM_CELL)
      DatasetGetNDimensions = ONE_INTEGER
    case(DIM_XY,DIM_XZ,DIM_YZ)
      DatasetGetNDimensions = TWO_INTEGER
    case(DIM_XYZ)
      DatasetGetNDimensions = THREE_INTEGER
    case default
      DatasetGetNDimensions = ZERO_INTEGER
  end select
      
end function DatasetGetNDimensions
      
! ************************************************************************** !
!
! DatasetInterpolateBetweenTimes: Interpolates dataset between two buffer times
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetInterpolateBetweenTimes(dataset,option)

  use Option_module

  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  PetscInt :: array_size,i
  PetscInt :: time_interpolation_method
  PetscReal :: weight2
  PetscInt :: time1_start, time1_end, time2_start, time2_end
  
  time_interpolation_method = INTERPOLATION_STEP
  
  if (dataset%buffer%cur_time_index >= &
      dataset%buffer%num_times_total) then
    ! dataset has reached the end of the time array. no interpolation
    ! needed
    array_size = size(dataset%rarray)
    time1_end = array_size*(dataset%buffer%cur_time_index - &
                            dataset%buffer%time_offset)
    time1_start = time1_end - array_size + 1
    dataset%rarray = dataset%buffer%rarray(time1_start:time1_end)
    dataset%buffer%cur_time_index = -1
    return
  endif
  
  if (associated(dataset%rarray)) then
    select case(time_interpolation_method)
      case(INTERPOLATION_STEP)
        array_size = size(dataset%rarray)
        time1_end = array_size*(dataset%buffer%cur_time_index - &
                                dataset%buffer%time_offset)
        time1_start = time1_end - array_size + 1
        dataset%rarray = dataset%buffer%rarray(time1_start:time1_end)
      case(INTERPOLATION_LINEAR)
        ! weight2 = (t-t1)/(t2-t1)
        weight2 = &
          (option%time - &
            dataset%buffer%time_array(dataset%buffer%cur_time_index)) / &
          (dataset%buffer%time_array(dataset%buffer%cur_time_index+1) - &
            dataset%buffer%time_array(dataset%buffer%cur_time_index))
        array_size = size(dataset%rarray)
        time1_end = array_size*(dataset%buffer%cur_time_index - &
                                dataset%buffer%time_offset)
        time1_start = time1_end - array_size + 1
        time2_end = time1_end + array_size
        time2_start = time1_start + array_size
        dataset%rarray = (1.d0-weight2) * &
                          dataset%buffer%rarray(time1_start:time1_end) + &
                          weight2 * &
                          dataset%buffer%rarray(time2_start:time2_end)
    end select
  endif

end subroutine DatasetInterpolateBetweenTimes

! ************************************************************************** !
!
! DatasetInterpolate: Interpolates data from the dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetInterpolateReal(dataset,xx,yy,zz,time,real_value,option)

  use Utility_module, only : InterpolateBilinear
  use Option_module
  
  implicit none
  
  type(dataset_type) :: dataset
  PetscReal, intent(in) :: xx, yy, zz
  PetscReal :: time
  PetscReal :: real_value
  type(option_type) :: option
  
  PetscInt :: spatial_interpolation_method
  PetscInt :: i, j, k
  PetscReal :: x, y, z
  PetscReal :: x1, x2, y1, y2
  PetscReal :: v1, v2, v3, v4
  PetscInt :: index
  PetscReal :: dx, dy
  PetscInt :: nx
  character(len=MAXWORDLENGTH) :: word
  
  call DatasetGetIndices(dataset,xx,yy,zz,i,j,k,x,y,z)
  
  spatial_interpolation_method = INTERPOLATION_LINEAR
  
  ! in the below, i,j,k,xx,yy,zz to not reflect the 
  ! coordinates of the problem domain in 3D.  They
  ! are transfored to the dimensions of the dataset
  select case(spatial_interpolation_method)
    case(INTERPOLATION_STEP)
    case(INTERPOLATION_LINEAR)
      select case(dataset%data_dim)
        case(DIM_X,DIM_Y,DIM_Z)
          if (i < 1 .or. i+1 > dataset%dims(1)) then 
            write(word,*) i
            word = adjustl(word)
            select case(dataset%data_dim)
              case(DIM_X)
                option%io_buffer = 'Out of x bounds, i = ' // trim(word)
              case(DIM_Y)
                option%io_buffer = 'Out of y bounds, j = ' // trim(word)
              case(DIM_Z)
                option%io_buffer = 'Out of z bounds, k = ' // trim(word)
            end select
            call printErrMsgByRank(option)
          endif
          dx = dataset%discretization(1)
          x1 = dataset%origin(1) + (i-1)*dx
          if (dataset%cell_centered) x1 = x1 + 0.5d0*dx
          v1 = dataset%rarray(i)
          v2 = dataset%rarray(i+1)
          real_value = v1 + (x-x1)/dx*(v2-v1)
        case(DIM_XY,DIM_XZ,DIM_YZ)
          if (i < 1 .or. i+1 > dataset%dims(1)) then
            write(word,*) i
            word = adjustl(word)
            select case(dataset%data_dim)
              case(DIM_XY,DIM_XZ)
                option%io_buffer = 'Out of x bounds, i = ' // trim(word)
              case(DIM_YZ)
                option%io_buffer = 'Out of y bounds, j = ' // trim(word)
            end select
            call printErrMsgByRank(option)
          endif
          if (j < 1 .or. j+1 > dataset%dims(2)) then
            write(word,*) j
            word = adjustl(word)
            select case(dataset%data_dim)
              case(DIM_XY)
                option%io_buffer = 'Out of y bounds, j = ' // trim(word)
              case(DIM_YZ,DIM_XZ)
                option%io_buffer = 'Out of z bounds, k = ' // trim(word)
            end select
            call printErrMsgByRank(option)
          endif
          dx = dataset%discretization(1)
          dy = dataset%discretization(2)
          nx = dataset%dims(1)

          x1 = dataset%origin(1) + (i-1)*dx
          if (dataset%cell_centered) x1 = x1 + 0.5d0*dx
          x2 = x1 + dx
          
          index = i + (j-1)*nx
          v1 = dataset%rarray(index)
          v2 = dataset%rarray(index+1)
          
          y1 = dataset%origin(2) + (j-1)*dy
          if (dataset%cell_centered) y1 = y1 + 0.5d0*dy
          y2 = y1 + dy
          
           ! really (j1-1+1)
          index = i + j*nx
          v3 = dataset%rarray(index)
          v4 = dataset%rarray(index+1)
          
          real_value = InterpolateBilinear(x,y,x1,x2,y1,y2,v1,v2,v3,v4)
        case(DIM_XYZ)
          option%io_buffer = 'Trilinear interpolation not yet supported'
          call printErrMsgByRank(option)
      end select
  end select
  
end subroutine DatasetInterpolateReal

! ************************************************************************** !
!
! DatasetGetIndices: Returns bounding indices for point in dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetGetIndices(dataset,xx,yy,zz,i,j,k,x,y,z)

  implicit none

  type(dataset_type) :: dataset
  PetscReal, intent(in)  :: xx, yy, zz
  PetscInt :: i, j, k
  PetscReal :: x, y, z
  
  select case(dataset%data_dim)
    ! since these are 1D array, always use first dimension
    case(DIM_X)
      x = xx
    case(DIM_Y)
      x = yy
    case(DIM_XY)
      x = xx
      y = yy
    case(DIM_XYZ)
      x = xx
      y = yy
      z = zz
    case(DIM_Z)
      x = zz
    case(DIM_XZ)
      x = xx
      y = zz
    case(DIM_YZ)
      x = yy
      y = zz
  end select

  if (dataset%cell_centered) then
    i = int((x - dataset%origin(1))/ &
            dataset%discretization(1) + 0.5d0)
    i = max(1,min(i,dataset%dims(1)-1))
    if (dataset%data_dim > DIM_Z) then ! at least 2D
      j = int((y - dataset%origin(2))/ &
              dataset%discretization(2) + 0.5d0)
      j = max(1,min(j,dataset%dims(2)-1))
    endif
    if (dataset%data_dim > DIM_YZ) then ! at least 3D
      k = int((z - dataset%origin(3))/ &
              dataset%discretization(3) + 0.5d0)
      k = max(1,min(k,dataset%dims(3)-1))
    endif
  else
    i = int((x - dataset%origin(1))/ &
            dataset%discretization(1) + 1.d0)
    if (dataset%data_dim > DIM_Z) then ! at least 2D
      j = int((y - dataset%origin(2))/ &
              dataset%discretization(2) + 1.d0)
    endif
    if (dataset%data_dim > DIM_YZ) then ! at least 3D
      k = int((z - dataset%origin(3))/ &
              dataset%discretization(3) + 1.d0)
    endif
  endif
  
end subroutine DatasetGetIndices

! ************************************************************************** !
!
! DatasetReorder: If a dataset is loaded from an HDF5 file, and it was
!                 multidimensional in the HDF5 file, the array needs to be
!                 reordered fro Fortran -> C indexing.  This subroutine
!                 takes care of the reordering.
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetReorder(dataset,option)

  use Option_module
  
  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  PetscInt, allocatable :: dims(:)
  PetscReal, allocatable :: temp_real(:)
  PetscInt :: i, j, k, t
  PetscInt :: nx, ny, nz, nt, nxXny, nyXnz, nxXnyXnz
  PetscInt :: count, index, idim
  PetscReal, pointer :: rarray(:)
  
  ! Not necessary for 1D arrays
  if (dataset%ndims < 2) return

  if (associated(dataset%buffer)) then
    rarray => dataset%buffer%rarray
  else
    rarray => dataset%rarray
  endif
  
  select case(dataset%ndims)
    case(TWO_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      if (associated(dataset%buffer)) then
        ny = dataset%buffer%num_times_in_buffer
      else
        ny = dataset%dims(TWO_INTEGER)
      endif
      count = 0
      allocate(temp_real(nx*ny))
      do i = 1, nx
        do j = 0, ny-1
          index = j*nx+i
          count = count+1
          temp_real(index) = rarray(count)
        enddo
      enddo  
    case(THREE_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      ny = dataset%dims(TWO_INTEGER)
      if (associated(dataset%buffer)) then
        nz = dataset%buffer%num_times_in_buffer
      else
        nz = dataset%dims(THREE_INTEGER)
      endif
      nxXny = nx*ny
      count = 0
      allocate(temp_real(nx*ny*nz))
      do i = 1, nx
        do j = 0, ny-1
          do k = 0, nz-1
            index = k*nxXny+j*nx+i
            count = count+1
            temp_real(index) = rarray(count)
          enddo
        enddo
      enddo  
    case(FOUR_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      ny = dataset%dims(TWO_INTEGER)
      nz = dataset%dims(THREE_INTEGER)
      if (associated(dataset%buffer)) then
        nt = dataset%buffer%num_times_in_buffer
      else
        nt = dataset%dims(FOUR_INTEGER)
      endif
      nxXny = nx*ny
      nxXnyXnz = nxXny*nz
      count = 0
      allocate(temp_real(nx*ny*nz*nt))
      do i = 1, nx
        do j = 0, ny-1
          do k = 0, nz-1
            do t = 0, nt-1
              index = t*nxXnyXnz+k*nxXny+j*nx+i
              count = count+1
              temp_real(index) = rarray(count)
            enddo
          enddo
        enddo
      enddo  
    case default
      write(option%io_buffer,*) dataset%ndims
      option%io_buffer = 'Dataset reordering not yet supported for rank ' // &
                         trim(adjustl(option%io_buffer)) // ' array.'
      call printErrMsgByRank(option)
  end select

  rarray = temp_real
  deallocate(temp_real)
  
end subroutine DatasetReorder

! ************************************************************************** !
!
! DatasetPrint: Prints dataset info
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetPrint(dataset,option)

  use Option_module

  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  option%io_buffer = 'TODO(geh): add DatasetPrint()'
  call printMsg(option)
            
end subroutine DatasetPrint

! ************************************************************************** !
!
! DatasetGetTimes: Fills an array of times based on a dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetGetTimes(option, dataset, max_sim_time, times)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(dataset_type) :: dataset
  PetscReal :: max_sim_time
  PetscReal, pointer :: times(:)
  
  PetscInt :: itime

  if (associated(dataset%buffer)) then
    do itime = 1, size(dataset%buffer%time_array)
      if (dataset%buffer%time_array(itime) > max_sim_time) exit
    enddo
    if (itime >= size(dataset%buffer%time_array)) then
      allocate(times(size(dataset%buffer%time_array)))
      times = dataset%buffer%time_array
    else
      allocate(times(itime))
      times(1:itime) = dataset%buffer%time_array(1:itime)
    endif
  endif
 
end subroutine DatasetGetTimes

! ************************************************************************** !
!
! DatasetBufferDestroy: Destroys a dataset buffer
! author: Glenn Hammond
! date: 10/28/11
!
! ************************************************************************** !
recursive subroutine DatasetBufferDestroy(dataset_buffer)

  implicit none
  
  type(dataset_buffer_type), pointer :: dataset_buffer
  
  if (.not.associated(dataset_buffer)) return

  if (associated(dataset_buffer%time_array)) &
    deallocate(dataset_buffer%time_array)
  nullify(dataset_buffer%time_array)
  if (associated(dataset_buffer%iarray)) deallocate(dataset_buffer%iarray)
  nullify(dataset_buffer%iarray)
  if (associated(dataset_buffer%rarray)) deallocate(dataset_buffer%rarray)
  nullify(dataset_buffer%rarray)
  
  deallocate(dataset_buffer)
  nullify(dataset_buffer)
  
end subroutine DatasetBufferDestroy

! ************************************************************************** !
!
! DatasetDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
recursive subroutine DatasetDestroy(dataset)

  implicit none
  
  type(dataset_type), pointer :: dataset
  
  if (.not.associated(dataset)) return
  
  if (associated(dataset%next)) then
    call DatasetDestroy(dataset%next)
  endif

  if (associated(dataset%iarray)) deallocate(dataset%iarray)
  nullify(dataset%iarray)
  if (associated(dataset%rarray)) deallocate(dataset%rarray)
  nullify(dataset%rarray)
  if (associated(dataset%dims)) deallocate(dataset%dims)
  nullify(dataset%dims)
  if (associated(dataset%discretization)) deallocate(dataset%discretization)
  nullify(dataset%discretization)
  if (associated(dataset%origin)) deallocate(dataset%origin)
  nullify(dataset%origin)
  if (associated(dataset%buffer)) then 
    call DatasetBufferDestroy(dataset%buffer)
  endif
  if (associated(dataset%dataset_map)) then
    if(associated(dataset%dataset_map%datatocell_ids)) &
      deallocate(dataset%dataset_map%datatocell_ids)
    if(associated(dataset%dataset_map%cell_ids_local)) &
      deallocate(dataset%dataset_map%cell_ids_local)
    if(associated(dataset%dataset_map%map)) &
      deallocate(dataset%dataset_map%map)
    nullify(dataset%dataset_map%map)
    nullify(dataset%dataset_map%datatocell_ids)
    nullify(dataset%dataset_map%cell_ids_local)
    deallocate(dataset%dataset_map)
  endif
  nullify(dataset%dataset_map)
  
  deallocate(dataset)
  nullify(dataset)
  
end subroutine DatasetDestroy

end module Dataset_Aux_module
