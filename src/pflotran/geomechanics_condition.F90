module Geomechanics_Condition_module
 
!  use Global_Aux_module
  use Dataset_Base_class
  use Time_Series_module
  
  implicit none

  private
  
#include "definitions.h"

  PetscInt, parameter                              :: NULL = 0
  PetscInt, parameter                              :: STEP = 1
  PetscInt, parameter                              :: LINEAR = 2

  type, public :: geomech_condition_dataset_type
    type(time_series_type), pointer                :: time_series
    class(dataset_base_type), pointer              ::  dataset
  end type geomech_condition_dataset_type
  
  type, public :: geomech_condition_type
    PetscInt                                       :: id    ! id from which condition can be referenced
    PetscBool                                      :: sync_time_with_update
    character(len=MAXWORDLENGTH)                   :: name    ! name of condition (e.g. boundary)
    PetscInt                                       :: num_sub_conditions
    PetscInt, pointer                              :: itype(:)
    character(len=MAXWORDLENGTH)                   :: time_units
    character(len=MAXWORDLENGTH)                   :: length_units
    type(geomech_sub_condition_type), pointer      :: displacement_x
    type(geomech_sub_condition_type), pointer      :: displacement_y
    type(geomech_sub_condition_type), pointer      :: displacement_z
    type(geomech_sub_condition_ptr_type), pointer  :: sub_condition_ptr(:)
    type(geomech_condition_type), pointer          :: next ! pointer to next condition_type for linked-lists
  end type geomech_condition_type
    
  type, public :: geomech_sub_condition_type
    PetscInt                                       :: itype ! integer describing type of condition
    PetscInt                                       :: isubtype
    character(len=MAXWORDLENGTH)                   :: ctype ! character string describing type of condition
    character(len=MAXWORDLENGTH)                   :: units      ! units
    character(len=MAXWORDLENGTH)                   :: name
    type(geomech_condition_dataset_type)           :: geomech_dataset
  end type geomech_sub_condition_type
  
  type, public :: geomech_sub_condition_ptr_type
    type(geomech_sub_condition_type), pointer      :: ptr
  end type geomech_sub_condition_ptr_type
    
  type, public :: geomech_condition_ptr_type
    type(geomech_condition_type), pointer          :: ptr
  end type geomech_condition_ptr_type
  
  type, public :: geomech_condition_list_type
    PetscInt :: num_conditions
    type(geomech_condition_type), pointer          :: first
    type(geomech_condition_type), pointer          :: last
    type(geomech_condition_type), pointer          :: array(:)    
  end type geomech_condition_list_type

  public :: GeomechConditionCreate, &
            GeomechConditionDestroy, &
            GeomechConditionRead, &
            GeomechConditionAddToList, &
            GeomechConditionInitList, &
            GeomechConditionDestroyList, &
            GeomechConditionGetPtrFromList, &
            GeomechConditionUpdate, &
            GeomechConditionPrint, &
            GeomechConditionIsTransient, &
            GeomechConditionDatasetGetTimes
            
  
contains

! ************************************************************************** !
!
! GeomechConditionCreate: Creates a condition
! author: Satish Karra, LANL
! date: 06/07/13
!
! ************************************************************************** !
function GeomechConditionCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type)                                :: option
  type(geomech_condition_type), pointer            :: GeomechConditionCreate
  
  type(geomech_condition_type), pointer            :: condition
  
  allocate(condition)
  nullify(condition%displacement_x)
  nullify(condition%displacement_y)
  nullify(condition%displacement_z)
  nullify(condition%sub_condition_ptr)
  nullify(condition%itype)
  nullify(condition%next)
  condition%sync_time_with_update = PETSC_FALSE
  condition%time_units = ''
  condition%length_units = ''
  condition%id = 0
  condition%num_sub_conditions = 0
  condition%name = ''
  
  GeomechConditionCreate => condition

end function GeomechConditionCreate

! ************************************************************************** !
!
! GeomechSubConditionCreate: Creates a sub_condition
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
function GeomechSubConditionCreate(ndof)

  use Option_module
  
  implicit none
  
  type(geomech_sub_condition_type), pointer       :: GeomechSubConditionCreate
  
  PetscInt                                        :: ndof
  
  type(geomech_sub_condition_type), pointer       :: sub_condition
  
  allocate(sub_condition)
  sub_condition%units = ''
  sub_condition%itype = 0
  sub_condition%isubtype = 0
  sub_condition%ctype = ''
  sub_condition%name = ''

  ! by default, create time series
  call GeomechConditionDatasetInit(sub_condition%geomech_dataset)
  sub_condition%geomech_dataset%time_series => TimeSeriesCreate()
  sub_condition%geomech_dataset%time_series%rank = ndof

  GeomechSubConditionCreate => sub_condition

end function GeomechSubConditionCreate

! ************************************************************************** !
!
! GeomechConditionDatasetInit: Initializes a dataset
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionDatasetInit(geomech_condition_dataset)

  implicit none
  
  type(geomech_condition_dataset_type)           :: geomech_condition_dataset

  nullify(geomech_condition_dataset%time_series)
  nullify(geomech_condition_dataset%dataset)
    
end subroutine GeomechConditionDatasetInit

! ************************************************************************** !
!
! GeomechSubConditionVerify: Verifies the data in a subcondition
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechSubConditionVerify(option, condition, sub_condition_name, &
                                  sub_condition, &
                                  default_time, &
                                  default_ctype, default_itype, &
                                  default_geomech_dataset, &
                                  destroy_if_null)

  use Option_module

  implicit none
  
  type(option_type)                                :: option
  type(geomech_condition_type)                     :: condition
  character(len=MAXWORDLENGTH)                     :: sub_condition_name
  type(geomech_sub_condition_type), pointer        :: sub_condition
  character(len=MAXWORDLENGTH)                     :: default_ctype
  PetscInt                                         :: default_itype
  PetscBool                                        :: default_cyclic
  PetscInt                                         :: default_interpolation
  PetscReal                                        :: default_time
  type(geomech_condition_dataset_type)             :: default_geomech_dataset
  PetscBool                                        :: destroy_if_null

  PetscInt :: array_size

  if (.not.associated(sub_condition)) return
  
  if (.not. (associated(sub_condition%geomech_dataset%time_series%values) .or. &
             associated(sub_condition%geomech_dataset%dataset))) then
    if (destroy_if_null) call GeomechSubConditionDestroy(sub_condition)
    return
  endif
  
  if (len_trim(sub_condition%ctype) == NULL_CONDITION) then
    option%io_buffer = 'TYPE of condition ' // trim(condition%name) // &
      ' ' // trim(sub_condition_name) // ' dataset not defined.'
    call printErrMsg(option)
  endif
  
  call GeomechConditionDatasetVerify(option,condition%name,sub_condition_name, &
                                  default_time,sub_condition%geomech_dataset, &
                                  default_geomech_dataset)

end subroutine GeomechSubConditionVerify


! ************************************************************************** !
!
! GeomechConditionDatasetVerify: Verifies the data in a dataset
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionDatasetVerify(option, condition_name, &
                                      sub_condition_name, &
                                      default_time, &
                                      dataset, default_dataset)
  use Option_module

  implicit none
  
  type(option_type)                                :: option
  character(len=MAXWORDLENGTH)                     :: condition_name
  character(len=MAXWORDLENGTH)                     :: sub_condition_name
  character(len=MAXWORDLENGTH)                     :: size1, size2
  PetscReal                                        :: default_time
  type(geomech_condition_dataset_type)             :: dataset
  type(geomech_condition_dataset_type)             :: default_dataset
  
  if (associated(dataset%dataset) .or. &
      associated(default_dataset%dataset)) then
    call TimeSeriesDestroy(dataset%time_series)
    !GB: Do not destroy default_dataset. For those modes, that have more
    !    more than 1 DOF, geomech_condition for different DOFs could be 
    !    specified via dataset or values
    !call TimeSeriesDestroy(default_dataset%time_series)
    if (associated(default_dataset%dataset) .and. &
        .not.associated(dataset%dataset)) then
      dataset%dataset => default_dataset%dataset
    endif
!geh: We cannot overwrite the default dataset pointer
!    if (associated(dataset%dataset) .and. &
!        .not.associated(default_dataset%dataset)) then
!      default_dataset%dataset => dataset%dataset
!    endif
  endif
  
  if (associated(dataset%time_series).and.associated(default_dataset%time_series)) then
    call TimeSeriesVerify(option, default_time, dataset%time_series, &
                          default_dataset%time_series)
  endif
  
end subroutine GeomechConditionDatasetVerify

! ************************************************************************** !
!
! GeomechConditionRead: Reads a condition from the input file
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionRead(condition,input,option)

  use Option_module
  use Input_module
  use String_module
  use Geomechanics_Logging_module  
  
  implicit none
  
  type(geomech_condition_type)                     :: condition
  type(input_type)                                 :: input
  type(option_type)                                :: option
  
  character(len=MAXSTRINGLENGTH)                   :: string
  character(len=MAXWORDLENGTH)                     :: word
  type(geomech_sub_condition_type), pointer :: pressure, flux, temperature, &
                                       concentration, enthalpy, rate, well,&
                                       sub_condition_ptr, saturation, &
                                       displacement_x, displacement_y, &
                                       displacement_z
  PetscReal                                        :: default_time
  PetscInt                                         :: default_iphase
  type(geomech_condition_dataset_type)             :: default_geomech_dataset
  type(geomech_condition_dataset_type)             :: default_well
  character(len=MAXWORDLENGTH)                     :: default_ctype
  PetscInt                                         :: default_itype
  PetscInt                                         :: array_size, idof
  PetscBool                                        :: found
  PetscBool                                        :: destroy_if_null
  PetscErrorCode                                   :: ierr

  call PetscLogEventBegin(geomech_logging%event_geomech_condition_read,ierr)

  default_time = 0.d0
  default_iphase = 0
  call GeomechConditionDatasetInit(default_geomech_dataset)
  default_geomech_dataset%time_series => TimeSeriesCreate()
  default_geomech_dataset%time_series%rank = 1
  default_geomech_dataset%time_series%interpolation_method = STEP
  default_geomech_dataset%time_series%is_cyclic = PETSC_FALSE

  displacement_x => GeomechSubConditionCreate(ONE_INTEGER)
  displacement_y => GeomechSubConditionCreate(ONE_INTEGER)
  displacement_z => GeomechSubConditionCreate(ONE_INTEGER)
  displacement_x%name = 'displacement_x'
  displacement_y%name = 'displacement_y'
  displacement_z%name = 'displacement_z'

  displacement_x%units = 'm'
  displacement_y%units = 'm'
  displacement_z%units = 'm'

  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC

  ! read the condition
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
      case('UNITS') ! read default units for condition arguments
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%time_units = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%length_units = trim(word)
          end select
        enddo
      case('CYCLIC')
        default_geomech_dataset%time_series%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')   
        call StringToLower(word)
        select case(word)
          case('step')
            default_geomech_dataset%time_series%interpolation_method = STEP
          case('linear') 
            default_geomech_dataset%time_series%interpolation_method = LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          call StringToUpper(word)
          select case(trim(word))
            case('PRESSURE')
            case('DISPLACEMENT_X')
              sub_condition_ptr => displacement_x
            case('DISPLACEMENT_Y')
              sub_condition_ptr => displacement_y
            case('DISPLACEMENT_Z')
              sub_condition_ptr => displacement_z
            case default
              option%io_buffer = 'keyword (' // trim(word) // &
                                 ') not recognized in condition,type'
              call printErrMsg(option)
          end select
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'TYPE','CONDITION')   
          call StringToLower(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('dirichlet')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('neumann')
              sub_condition_ptr%itype = NEUMANN_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case default
              option%io_buffer = 'bc type "' // trim(word) // &
                                 '" not recognized in condition,type'
              call printErrMsg(option)
          end select
        enddo
      case('TIME','TIMES')
        call InputReadDouble(input,option,default_time)
        call InputErrorMsg(input,option,'TIME','CONDITION')   
      case('DISPLACEMENT_X')
        call GeomechConditionReadValues(input,option,word,string, &
                                     displacement_x%geomech_dataset, &
                                     displacement_x%units)
      case('DISPLACEMENT_Y')
        call GeomechConditionReadValues(input,option,word,string, &
                                     displacement_y%geomech_dataset, &
                                     displacement_y%units) 
      case('DISPLACEMENT_Z')
        call GeomechConditionReadValues(input,option,word,string, &
                                     displacement_z%geomech_dataset, &
                                     displacement_z%units)
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                           ' not recognized in geomech condition'
        call printErrMsg(option)                                 
    end select 
  
  enddo  
  
  word = 'displacement_x'
  call GeomechSubConditionVerify(option,condition,word,displacement_x,default_time, &
                              default_ctype, default_itype, &
                              default_geomech_dataset,PETSC_TRUE)
  word = 'displacement_y'
  call GeomechSubConditionVerify(option,condition,word,displacement_y,default_time, &
                              default_ctype, default_itype, &
                              default_geomech_dataset,PETSC_TRUE)
  word = 'displacement_z'
  call GeomechSubConditionVerify(option,condition,word,displacement_z,default_time, &
                              default_ctype, default_itype, &
                              default_geomech_dataset,PETSC_TRUE)

  if (.not.associated(displacement_x)) then
    option%io_buffer = 'displacement_x condition null in condition: ' // &
                        trim(condition%name)      
    call printErrMsg(option)
  endif                         
  condition%displacement_x => displacement_x

  if (.not.associated(displacement_y)) then
    option%io_buffer = 'displacement_y condition null in condition: ' // &
                        trim(condition%name)      
    call printErrMsg(option)
  endif                         
  condition%displacement_y => displacement_y

  if (.not.associated(displacement_z)) then
    option%io_buffer = 'displacement_z condition null in condition: ' // &
                        trim(condition%name)      
    call printErrMsg(option)
  endif                         
  condition%displacement_z => displacement_z

  condition%num_sub_conditions = 3
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, 3
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo

  ! must be in this order, which matches the dofs i problem
  condition%sub_condition_ptr(ONE_INTEGER)%ptr => displacement_x
  condition%sub_condition_ptr(TWO_INTEGER)%ptr => displacement_y
  condition%sub_condition_ptr(THREE_INTEGER)%ptr => displacement_z

  allocate(condition%itype(THREE_INTEGER))
  condition%itype = 0
  condition%itype(ONE_INTEGER) = displacement_x%itype
  condition%itype(TWO_INTEGER) = displacement_y%itype
  condition%itype(THREE_INTEGER) = displacement_z%itype
      
  call GeomechConditionDatasetDestroy(default_geomech_dataset)
    
  call PetscLogEventEnd(geomech_logging%event_geomech_condition_read,ierr)

end subroutine GeomechConditionRead

! ************************************************************************** !
!
! GeomechConditionReadValues: Read the value(s) of a condition variable
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionReadValues(input,option,keyword,string, &
                                      geomech_dataset,units)

  use Input_module
  use String_module
  use Option_module
  use Geomechanics_Logging_module
  use HDF5_Aux_module
  use Units_module
  use Dataset_Base_class
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

  implicit none
  
  type(input_type)                                 :: input
  type(option_type)                                :: option
  character(len=MAXSTRINGLENGTH)                   :: string
  character(len=MAXWORDLENGTH)                     :: keyword
  type(geomech_condition_dataset_type)             :: geomech_dataset
  character(len=MAXWORDLENGTH)                     :: units
  
  character(len=MAXSTRINGLENGTH)                   :: string2
  character(len=MAXSTRINGLENGTH)                   :: filename, hdf5_path
  character(len=MAXWORDLENGTH)                     :: word, realization_word
  character(len=MAXSTRINGLENGTH)                   :: error_string
  PetscInt                                         :: length, i, icount
  PetscInt                                         :: irank
  PetscInt                                         :: ndims
  PetscInt, pointer                                :: dims(:)
  PetscReal, pointer                               :: real_buffer(:)
  type(input_type), pointer                        :: input2
  PetscErrorCode                                   :: ierr

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T)                                   :: file_id
  integer(HID_T)                                   :: prop_id
  PetscMPIInt                                      :: hdf5_err
#endif

  call PetscLogEventBegin(geomech_logging%event_geomech_condition_read_values,ierr)    

  nullify(input2)
  filename = ''
  realization_word = ''
  hdf5_path = ''
  
  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','CONDITION')
  call StringToLower(word)
  length = len_trim(word)
  if (length == FOUR_INTEGER .and. StringCompare(word,'file',FOUR_INTEGER)) then 
    input%err_buf2 = trim(keyword) // ', FILE'
    input%err_buf = 'keyword'
    call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
    if (input%ierr == 0) then
      filename = string2
    else
      option%io_buffer = 'The ability to read realization dependent datasets outside the DATASET block is no longer supported'
      call printErrMsg(option)
    endif
    
    if (len_trim(filename) < 2) then
      option%io_buffer = 'No filename listed under Geomech_Condition: ' // &
                         trim(keyword)
      call printErrMsg(option)
    endif

    if (index(filename,'.h5') > 0) then
#if !defined(PETSC_HAVE_HDF5)
      write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                               &"read HDF5 formatted geomech conditions.")')
      call printErrMsg(option)
#else   
      if (len_trim(hdf5_path) < 1) then
        option%io_buffer = 'No hdf5 path listed under Geomech_Condition: ' // &
                           trim(keyword)
        call printErrMsg(option)
      endif

      call h5open_f(hdf5_err)
      option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
      call printMsg(option)
      call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
      call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
      call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
      call h5pclose_f(prop_id,hdf5_err)

      hdf5_path = trim(hdf5_path) // trim(realization_word)
      call HDF5ReadNDimRealArray(option,file_id,hdf5_path,ndims,dims, &
                                 real_buffer)
      option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
      call printMsg(option)  
      call h5fclose_f(file_id,hdf5_err)
      call h5close_f(hdf5_err)
      
      ! dims(1) = size of array
      ! dims(2) = number of data point in time
      if (dims(1)-1 == geomech_dataset%time_series%rank) then
        ! alright, the 2d data is layed out in C-style.  now place it in
        ! the appropriate arrays
        allocate(geomech_dataset%time_series%times(dims(2)))
        geomech_dataset%time_series%times = -999.d0
        allocate(geomech_dataset%time_series%values(geomech_dataset%time_series%rank,dims(2))) 
        geomech_dataset%time_series%values = -999.d0
        icount = 1
        do i = 1, dims(2)
          geomech_dataset%time_series%times(i) = real_buffer(icount)
          icount = icount + 1
          do irank = 1, geomech_dataset%time_series%rank
            geomech_dataset%time_series%values(irank,i) = real_buffer(icount)
            icount = icount + 1
          enddo
        enddo  
      else
        option%io_buffer = 'HDF condition data set rank does not match' // &
          'rank of internal data set.  Email Glenn for additions'
        call printErrMsg(option)
      endif
      if (associated(dims)) deallocate(dims)
      nullify(dims)
      if (associated(real_buffer)) deallocate(real_buffer)
      nullify(real_buffer)
#endif    
    else
      i = index(filename,'.',PETSC_TRUE)
      if (i > 2) then
        filename = filename(1:i-1) // trim(realization_word) // filename(i:)
      else
        filename = trim(filename) // trim(realization_word)
      endif
      input2 => InputCreate(IUNIT_TEMP,filename,option)
      if (geomech_dataset%time_series%rank<=3) then
        call GeomechConditionReadValuesFromFile(input2,geomech_dataset,option)
      else
        call GeomechConditionReadValuesFromFile2(input2,geomech_dataset,option)
      endif  
      call InputDestroy(input2)
    endif
  else if (StringCompare(word,'dataset')) then
    call InputReadWord(input,option,word,PETSC_TRUE)    
    input%err_buf2 = trim(keyword) // ', DATASET'
    input%err_buf = 'dataset name'
    call InputErrorMsg(input,option)
    geomech_dataset%dataset => DatasetBaseCreate()
    geomech_dataset%dataset%name = word
  else if (length==FOUR_INTEGER .and. StringCompare(word,'list',length)) then  !sp 
    if (geomech_dataset%time_series%rank <= 3) then
      call GeomechConditionReadValuesFromFile(input,geomech_dataset,option)
    else
      call GeomechConditionReadValuesFromFile2(input,geomech_dataset,option)
    endif  
  else
    input%buf = trim(string2)
    allocate(geomech_dataset%time_series%values(geomech_dataset%time_series%rank,1))
    do irank=1,geomech_dataset%time_series%rank
      call InputReadDouble(input,option,geomech_dataset%time_series%values(irank,1))
      write(input%err_buf,'(a,i2)') trim(keyword) // ' dataset_values, irank = ', irank
      input%err_buf2 = 'CONDITION'
      call InputErrorMsg(input,option) 
    enddo
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (InputError(input)) then
      word = trim(keyword) // ' UNITS'
      call InputDefaultMsg(input,option,word)
    else
      units = trim(word)
      geomech_dataset%time_series%values(1:geomech_dataset%time_series%rank,1) = &
        UnitsConvertToInternal(units,option) * &
        geomech_dataset%time_series%values(1:geomech_dataset%time_series%rank,1)
    endif
  endif
  call PetscLogEventEnd(geomech_logging%event_geomech_condition_read_values,ierr)    

end subroutine GeomechConditionReadValues

! ************************************************************************** !
!
! GeomechConditionReadValuesFromFile: Read values from a external file
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionReadValuesFromFile(input,geomech_dataset,option)

  use Input_module
  use String_module
  use Utility_module
  use Option_module
  use Units_module

  implicit none
  
  type(input_type)                                 :: input
  type(geomech_condition_dataset_type)             :: geomech_dataset
  type(option_type)                                :: option

  character(len=MAXWORDLENGTH)                     :: time_units, data_units
  character(len=MAXSTRINGLENGTH)                   :: string
  character(len=MAXWORDLENGTH)                     :: word
  PetscReal, pointer                               :: temp_times(:), &
                                                      temp_array1(:), &
                                                      temp_array2(:), &
                                                      temp_array3(:)
  PetscReal                                        :: temp_time
  PetscReal                                        :: conversion
  PetscInt                                         :: max_size
  PetscInt                                         :: temp_max_size
  PetscInt                                         :: count, i, status
  PetscErrorCode                                   :: ierr
  
  time_units = ''
  data_units = ''
  max_size = 1000
  allocate(temp_times(max_size))
  allocate(temp_array1(max_size))
  temp_times = 0.d0
  temp_array1 = 0.d0

  if (geomech_dataset%time_series%rank > 1) then
    allocate(temp_array2(max_size))
    temp_array2 = 0.d0
  endif
  
  if (geomech_dataset%time_series%rank > 2) then
    allocate(temp_array3(max_size))
    temp_array3 = 0.d0
  endif

  count = 0
  ierr = 0
  do
    call InputReadFlotranString(input,option)
    ! reach the end of file or close out block
    if (InputError(input)) exit  ! check for end of file
    if (InputCheckExit(input,option)) exit  ! check for end of list
    ! check for units on first line
    if (count == 0) then
      string = input%buf
      call InputReadWord(string,word,PETSC_TRUE,ierr)
      call StringToUpper(word)
      select case(word)
        case('TIME_UNITS')
          call InputReadWord(string,time_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'TIME_UNITS','CONDITION FILE')
          call StringToLower(time_units) 
          cycle
        case('DATA_UNITS')
          call InputReadWord(string,data_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'DATA_UNITS','CONDITION FILE')
          call StringToLower(data_units) 
          cycle
      end select
    endif
    count = count + 1
    call InputReadDouble(input,option,temp_times(count))
    call InputErrorMsg(input,option,'time','CONDITION FILE')   
    call InputReadDouble(input,option,temp_array1(count))
    call InputErrorMsg(input,option,'array1','CONDITION FILE')
    if (geomech_dataset%time_series%rank > 1) then
      call InputReadDouble(input,option,temp_array2(count))
      call InputErrorMsg(input,option,'array2','CONDITION FILE') 
    endif
    if (geomech_dataset%time_series%rank > 2) then
      call InputReadDouble(input,option,temp_array3(count))
      call InputErrorMsg(input,option,'array3','CONDITION FILE') 
    endif
    if (count+1 > max_size) then
      temp_max_size = max_size
      call reallocateRealArray(temp_times,max_size) 
      ! careful.  reallocateRealArray double max_size every time.
      i = temp_max_size
      call reallocateRealArray(temp_array1,i) 
      if (geomech_dataset%time_series%rank > 1) then
        i = temp_max_size
        call reallocateRealArray(temp_array2,i)
      endif
      if (geomech_dataset%time_series%rank > 2) then
        i = temp_max_size
        call reallocateRealArray(temp_array3,i)
      endif
    endif  
  enddo
  
  if (associated(geomech_dataset%time_series%times)) then
    if (count /= size(geomech_dataset%time_series%times,1) .and. &
        OptionPrintToScreen(option)) then
      print *, 'Number of times (', count, ') in ', trim(input%filename), &
               ' does not match previous allocation: ', &
               size(geomech_dataset%time_series%times,1)
      stop
    endif
    do i=1,count
      if (dabs(geomech_dataset%time_series%times(i)-temp_times(i)) > 1.d-8 .and. &
          OptionPrintToScreen(option)) then
        print *, 'Time (', temp_times(i), ') in ', trim(input%filename), &
                 ' does not match previous allocation time: ', &
                 geomech_dataset%time_series%times(i), i
        stop
      endif
    enddo
  else
    allocate(geomech_dataset%time_series%times(count))
  endif

  if (associated(geomech_dataset%time_series%values)) &
    deallocate(geomech_dataset%time_series%values)
  allocate(geomech_dataset%time_series%values(geomech_dataset%time_series%rank,count))

  geomech_dataset%time_series%times(1:count) = temp_times(1:count)
  geomech_dataset%time_series%values(1,1:count) = temp_array1(1:count)
  if (geomech_dataset%time_series%rank > 1) &
    geomech_dataset%time_series%values(2,1:count) = temp_array2(1:count)
  if (geomech_dataset%time_series%rank > 2) &
    geomech_dataset%time_series%values(3,1:count) = temp_array3(1:count)
  
  deallocate(temp_times)
  deallocate(temp_array1)
  if (geomech_dataset%time_series%rank > 1) deallocate(temp_array2)
  if (geomech_dataset%time_series%rank > 2) deallocate(temp_array3)

  if (len_trim(time_units) > 0) then
    ! Times
    conversion = UnitsConvertToInternal(time_units,option)
    geomech_dataset%time_series%times(1:count) = conversion * &
                                   geomech_dataset%time_series%times(1:count)
  endif
  if (len_trim(data_units) > 0) then
    ! Data
    conversion = UnitsConvertToInternal(data_units,option)
    geomech_dataset%time_series%values(1:geomech_dataset%time_series%rank,1:count) = &
      conversion * &
      geomech_dataset%time_series%values(1:geomech_dataset%time_series%rank,1:count)
  endif
  
end subroutine GeomechConditionReadValuesFromFile
! ************************************************************************** !
! 
! GeomechConditionReadValuesFromFile: Read values from a external file with 4 more columns
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionReadValuesFromFile2(input,geomech_dataset,option)

  use Input_module
  use String_module
  use Utility_module
  use Option_module
  use Units_module

  implicit none
  
  type(input_type)                                 :: input
  type(geomech_condition_dataset_type)             :: geomech_dataset
  type(option_type)                                :: option

  character(len=MAXWORDLENGTH)                     :: time_units, data_units
  character(len=MAXSTRINGLENGTH)                   :: string
  character(len=MAXWORDLENGTH)                     :: word
  PetscReal, pointer                               :: temp_times(:), &
                                                      temp_array(:,:)
  PetscReal                                        :: temp_time
  PetscReal                                        :: conversion
  PetscInt                                         :: max_size
  PetscInt                                         :: temp_max_size
  PetscInt                                         :: count, i, status
  PetscErrorCode                                   :: ierr
  
  time_units = ''
  data_units = ''
  max_size = 1000
  allocate(temp_times(max_size))
  allocate(temp_array(geomech_dataset%time_series%rank, max_size)) 

  count = 0
  ierr = 0
  do
    call InputReadFlotranString(input,option)
    ! reach the end of file or close out block
    if (InputError(input)) exit  ! check for end of file
    if (InputCheckExit(input,option)) exit  ! check for end of list
    ! check for units on first line
    if (count == 0) then
      string = input%buf
      call InputReadWord(string,word,PETSC_TRUE,ierr)
      call StringToUpper(word)
      select case(word)
        case('TIME_UNITS')
          call InputReadWord(string,time_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'TIME_UNITS','CONDITION FILE')
          call StringToLower(time_units) 
          cycle
        case('DATA_UNITS')
          call InputReadWord(string,data_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'DATA_UNITS','CONDITION FILE')
          call StringToLower(data_units) 
          cycle
      end select
    endif
    count = count + 1
    call InputReadDouble(input,option,temp_times(count))
    call InputErrorMsg(input,option,'time','CONDITION FILE')   


   do i =1, geomech_dataset%time_series%rank
      call InputReadDouble(input,option,temp_array(i,count))
      call InputErrorMsg(input,option,'array:',' CONDITION FILE') 
   enddo

 enddo
 if (count+1 > max_size) then
   print *, 'Number of times (', count, ') in ', trim(input%filename), &
          ' exceed 1000 '
   stop
 endif 
  
 if (associated(geomech_dataset%time_series%times)) then
   if (count /= size(geomech_dataset%time_series%times,1) .and. &
        OptionPrintToScreen(option)) then
      print *, 'Number of times (', count, ') in ', trim(input%filename), &
               ' does not match previous allocation: ', &
               size(geomech_dataset%time_series%times,1)
      stop
    endif
    do i=1,count
      if (dabs(geomech_dataset%time_series%times(i)-temp_times(i)) > 1.d-8 .and. &
          OptionPrintToScreen(option)) then
        print *, 'Time (', temp_times(i), ') in ', trim(input%filename), &
                 ' does not match previous allocation time: ', &
                 geomech_dataset%time_series%times(i), i
        stop
      endif
    enddo
  else
    allocate(geomech_dataset%time_series%times(count))
  endif

  if (associated(geomech_dataset%time_series%values)) &
    deallocate(geomech_dataset%time_series%values)
  allocate(geomech_dataset%time_series%values(geomech_dataset%time_series%rank,count))

  geomech_dataset%time_series%times(1:count) = temp_times(1:count)
  geomech_dataset%time_series%values(:,1:count) = temp_array(:,1:count)
  
    
  deallocate(temp_times)
  deallocate(temp_array)

  if (len_trim(time_units) > 0) then
    ! Times
    conversion = UnitsConvertToInternal(time_units,option)
    geomech_dataset%time_series%times(1:count) = conversion * &
      geomech_dataset%time_series%times(1:count)
  endif
  if (len_trim(data_units) > 0) then
    ! Data
    conversion = UnitsConvertToInternal(data_units,option)
    geomech_dataset%time_series%values(1:geomech_dataset%time_series%rank,1:count) = &
      conversion * &
      geomech_dataset%time_series%values(1:geomech_dataset%time_series%rank,1:count)
  endif
  
end subroutine GeomechConditionReadValuesFromFile2


! ************************************************************************** !
!
! GeomechConditionDatasetGetTimes: Fills an array of times based on dataset
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionDatasetGetTimes(option, sub_condition, &
                                        max_sim_time, times)
  use Option_module
  use Time_Storage_module

  implicit none
  
  type(option_type)                                :: option
  type(geomech_sub_condition_type), pointer        :: sub_condition
  PetscReal                                        :: max_sim_time
  PetscReal, pointer                               :: times(:)
  
  type(geomech_condition_dataset_type), pointer    :: geomech_dataset
  
  geomech_dataset => sub_condition%geomech_dataset

  if (associated(geomech_dataset%time_series) .and. &
      associated(geomech_dataset%dataset)) then
    option%io_buffer = 'GeomechConditionDatasetGetTimes() currently does not ' // &
                       'support both time_series and datasets.'
    call printErrMsg(option)
  endif
  if (associated(geomech_dataset%time_series)) then
    call TimeSeriesGetTimes(option, geomech_dataset%time_series, max_sim_time, &
                            times)
  endif
  if (associated(geomech_dataset%dataset)) then
    call TimeStorageGetTimes(geomech_dataset%dataset%time_storage, option, &
                             max_sim_time, times)
  endif
 
end subroutine GeomechConditionDatasetGetTimes

! ************************************************************************** !
!
! GeomechConditionPrint: Prints geomech condition info
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionPrint(condition,option)

  use Option_module

  implicit none
  
  type(geomech_condition_type)                     :: condition
  type(option_type)                                :: option
  
  character(len=MAXSTRINGLENGTH)                   :: string
  PetscInt :: i

99 format(/,80('-'))

  write(option%fid_out,'(/,2x,''Geomech Condition: '',a)') trim(condition%name)

  if (condition%sync_time_with_update) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(4x,''Synchronize time with update: '', a)') trim(string)
  write(option%fid_out,'(4x,''Time units: '', a)') trim(condition%time_units)
  write(option%fid_out,'(4x,''Length units: '', a)') trim(condition%length_units)
  
  do i=1, condition%num_sub_conditions
    call GeomechConditionPrintSubCondition(condition%sub_condition_ptr(i)%ptr, &
                                        option)
  enddo
  write(option%fid_out,99)
  
end subroutine GeomechConditionPrint

! ************************************************************************** !
!
! GeomechConditionPrintSubCondition: Prints geomech subcondition info
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionPrintSubCondition(subcondition,option)

  use Option_module

  implicit none
  
  type(geomech_sub_condition_type)                 :: subcondition
  type(option_type)                                :: option
  
  character(len=MAXSTRINGLENGTH)                   :: string
  
  write(option%fid_out,'(/,4x,''Sub Condition: '',a)') trim(subcondition%name)
  select case(subcondition%itype)
    case(DIRICHLET_BC)
      string = 'dirichlet'
    case(NEUMANN_BC)
      string = 'neumann'
    case(ZERO_GRADIENT_BC)
      string = 'zero gradient'
  end select
  100 format(6x,'Type: ',a)  
  write(option%fid_out,100) trim(string)
  
  110 format(6x,a)  

  write(option%fid_out,110) 'Dataset:'
  if (associated(subcondition%geomech_dataset%time_series)) then
    call TimeSeriesPrint(subcondition%geomech_dataset%time_series,option)
  endif
  if (associated(subcondition%geomech_dataset%dataset)) then
!geh    call DatasetPrint(subcondition%geomech_dataset%dataset,option)
    option%io_buffer = 'TODO(geh): add DatasetPrint()'
    call printMsg(option)
  endif
            
end subroutine GeomechConditionPrintSubCondition
 
! ************************************************************************** !
!
! GeomechConditionDatasetPrint: Prints geomech condition dataset info
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionDatasetPrint(geomech_dataset,option)

  use Option_module

  implicit none
  
  type(geomech_condition_dataset_type)             :: geomech_dataset
  type(option_type)                                :: option
  
  if(associated(geomech_dataset%time_series)) then
    call TimeSeriesPrint(geomech_dataset%time_series,option)
  endif
  if(associated(geomech_dataset%dataset)) then
  endif
  
end subroutine GeomechConditionDatasetPrint

! ************************************************************************** !
!
! GeomechConditionUpdate: Updates a transient condition
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionUpdate(condition_list,option,time)

  use Option_module
  
  implicit none
  
  type(geomech_condition_list_type)                :: condition_list
  type(option_type)                                :: option
  PetscReal                                        :: time
  
  type(geomech_condition_type), pointer            :: condition
  type(geomech_sub_condition_type), pointer        :: sub_condition
  PetscInt                                         :: isub_condition   
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
    do isub_condition = 1, condition%num_sub_conditions

      sub_condition => condition%sub_condition_ptr(isub_condition)%ptr
      
      if (associated(sub_condition)) then
        call GeomechSubConditionUpdateDataset(option,time, &
                                           sub_condition%geomech_dataset)
      endif
      
    enddo
      
    condition => condition%next
    
  enddo
  
end subroutine GeomechConditionUpdate

! ************************************************************************** !
!
! GeomechSubConditionUpdateDataset: Updates a transient condition dataset
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechSubConditionUpdateDataset(option,time,geomech_condition_dataset)

  use Option_module
  use Dataset_XYZ_class
  use Dataset_Map_class
  use Dataset_Common_HDF5_class
  
  implicit none
  
  type(option_type)                              :: option
  PetscReal                                      :: time
  type(geomech_condition_dataset_type)           :: geomech_condition_dataset

  class(dataset_xyz_type), pointer               :: dataset_xyz
  class(dataset_map_type), pointer               :: dataset_map

  if (associated(geomech_condition_dataset%time_series)) then
    if (time < 1.d-40 .or. &
        geomech_condition_dataset%time_series%is_transient) then
      call TimeSeriesUpdate(option,time,geomech_condition_dataset%time_series)
    endif
  endif
  
  if (associated(geomech_condition_dataset%dataset)) then
    select type(dataset=>geomech_condition_dataset%dataset)
      class is(dataset_common_hdf5_type)
        if ((time < 1.d-40 .or. &
             dataset%is_transient) .and. &
            .not.dataset%is_cell_indexed) then
          select type(common_hdf5_dataset=>dataset)
            class is(dataset_xyz_type)
              dataset_xyz => common_hdf5_dataset
              call DatasetXYZLoad(dataset_xyz,option)
            class is(dataset_map_type)
              dataset_map => common_hdf5_dataset
              call DatasetMapLoad(dataset_map,option)
          end select
        endif
    end select
  endif
  
end subroutine GeomechSubConditionUpdateDataset

! ************************************************************************** !
!
! GeomechConditionInitList: Initializes a condition list
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionInitList(list)

  implicit none

  type(geomech_condition_list_type)                :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine GeomechConditionInitList

! ************************************************************************** !
!
! GeomechConditionAddToList: Adds a new condition to a condition list
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechConditionAddToList(new_condition,list)

  implicit none
  
  type(geomech_condition_type), pointer             :: new_condition
  type(geomech_condition_list_type)                 :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine GeomechConditionAddToList

! ************************************************************************** !
!
! GeomechConditionGetPtrFromList: Returns a pointer to the condition matching &
!                          condition_name
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
function GeomechConditionGetPtrFromList(condition_name,condition_list)

  use String_module
  
  implicit none
  
  type(geomech_condition_type), pointer     :: GeomechConditionGetPtrFromList
  character(len=MAXWORDLENGTH)              :: condition_name
  type(geomech_condition_list_type)         :: condition_list
 
  PetscInt                                  :: length
  type(geomech_condition_type), pointer     :: condition
    
  nullify(GeomechConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                      length)) then
      GeomechConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function GeomechConditionGetPtrFromList

! ************************************************************************** !
!
! GeomechConditionIsTransient: Returns PETSC_TRUE
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
function GeomechConditionIsTransient(condition)

  implicit none
  
  type(geomech_condition_type)                 :: condition
 
  PetscBool                                    :: GeomechConditionIsTransient
  
  GeomechConditionIsTransient = PETSC_FALSE

  if (GeomechSubConditionIsTransient(condition%displacement_x) .or. &
      GeomechSubConditionIsTransient(condition%displacement_y) .or. &
      GeomechSubConditionIsTransient(condition%displacement_z)) then
    GeomechConditionIsTransient = PETSC_TRUE
  endif
  
end function GeomechConditionIsTransient

! ************************************************************************** !
!
! GeomechSubConditionIsTransient: Returns PETSC_TRUE
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
function GeomechSubConditionIsTransient(sub_condition)

  implicit none
  
  type(geomech_sub_condition_type), pointer :: sub_condition
  
  PetscBool                                 :: GeomechSubConditionIsTransient
  
  GeomechSubConditionIsTransient = PETSC_FALSE

  if (associated(sub_condition)) then
    if (GeomechDatasetIsTransient(sub_condition%geomech_dataset)) then
      GeomechSubConditionIsTransient = PETSC_TRUE
    endif
  endif  
  
end function GeomechSubConditionIsTransient

! ************************************************************************** !
!
! GeomechDatasetIsTransient: Returns PETSC_TRUE
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
function GeomechDatasetIsTransient(geomech_dataset)

  use Dataset_Common_HDF5_class

  implicit none
  
  type(geomech_condition_dataset_type)           :: geomech_dataset
  
  PetscBool                                      :: GeomechDatasetIsTransient
  
  GeomechDatasetIsTransient = PETSC_FALSE

  if (associated(geomech_dataset%time_series)) then
    if (geomech_dataset%time_series%is_transient) then
      GeomechDatasetIsTransient = PETSC_TRUE
    endif
  endif
  if (associated(geomech_dataset%dataset)) then
    select type(dataset=>geomech_dataset%dataset)
      class is(dataset_common_hdf5_type)
        if (dataset%is_transient) then
          GeomechDatasetIsTransient = PETSC_TRUE
        endif
    end select
  endif
  
end function GeomechDatasetIsTransient

! ************************************************************************** !
!
! GeomechConditionDestroyList: Deallocates a list of conditions
! author: Satish Karra, LANL
! date: 06/06/13
!
! ************************************************************************** !
subroutine GeomechConditionDestroyList(condition_list)

  implicit none
  
  type(geomech_condition_list_type), pointer        :: condition_list
  
  type(geomech_condition_type), pointer             :: condition, &
                                                       prev_condition
  
  if (.not.associated(condition_list)) return
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call GeomechConditionDestroy(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine GeomechConditionDestroyList

! ************************************************************************** !
!
! GeomechConditionDestroy: Deallocates a condition
! author: Satish Karra, LANL
! date: 10/23/07
!
! ************************************************************************** !
subroutine GeomechConditionDestroy(condition)

  implicit none
  
  type(geomech_condition_type), pointer            :: condition
  
  PetscInt                                         :: i
  
  if (.not.associated(condition)) return
  
  if (associated(condition%sub_condition_ptr)) then
    do i=1,condition%num_sub_conditions
      call GeomechSubConditionDestroy(condition%sub_condition_ptr(i)%ptr)
    enddo
    deallocate(condition%sub_condition_ptr)
    nullify(condition%sub_condition_ptr)
  endif

  if (associated(condition%itype)) deallocate(condition%itype)
  nullify(condition%itype)
  
  nullify(condition%displacement_x)
  nullify(condition%displacement_y)
  nullify(condition%displacement_z)
  
  nullify(condition%next)  
  
  deallocate(condition)
  nullify(condition)

end subroutine GeomechConditionDestroy


! ************************************************************************** !
!
! GeomechSubConditionDestroy: Destroys a sub_condition
! author: Satish Karra, LANL
! date: 02/04/08
!
! ************************************************************************** !
subroutine GeomechSubConditionDestroy(sub_condition)

  implicit none
  
  type(geomech_sub_condition_type), pointer        :: sub_condition
  
  if (.not.associated(sub_condition)) return
  
  call GeomechConditionDatasetDestroy(sub_condition%geomech_dataset)
 ! call GeomechConditionDatasetDestroy(sub_condition%gradient)

  deallocate(sub_condition)
  nullify(sub_condition)

end subroutine GeomechSubConditionDestroy

! ************************************************************************** !
!
! GeomechConditionDatasetDestroy: Destroys a dataset associated with a sub_condition
! author: Satish Karra, LANL
! date: 02/04/08
!
! ************************************************************************** !
subroutine GeomechConditionDatasetDestroy(geomech_condition_dataset)

  implicit none
  
  type(geomech_condition_dataset_type)           :: geomech_condition_dataset
  
  if (associated(geomech_condition_dataset%time_series)) then
    call TimeSeriesDestroy(geomech_condition_dataset%time_series)
  endif
  nullify(geomech_condition_dataset%time_series)
  ! since datasets reside in the dataset list, we simple unlink the pointer
  ! without destroying the object.
  if (associated(geomech_condition_dataset%dataset)) &
    nullify(geomech_condition_dataset%dataset)

end subroutine GeomechConditionDatasetDestroy

end module Geomechanics_Condition_module