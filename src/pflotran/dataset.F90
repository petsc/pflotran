module Dataset_module
 
  use Dataset_Aux_module
  
  implicit none

  private
  
#include "definitions.h"

  public :: DatasetLoad, &
            DatasetProcessDatasets, &
            DatasetIsCellIndexed

contains

! ************************************************************************** !
!
! DatasetLoad: Loads a dataset from file
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetLoad(dataset,option)

  use Option_module
  use HDF5_Aux_module
  use Utility_module, only : Equal

  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option

  PetscBool :: read_dataset
  PetscBool :: interpolate_dataset
  PetscInt :: itime
  
    ! update the offset of the data array
  read_dataset = PETSC_FALSE
  interpolate_dataset = PETSC_FALSE
  if (associated(dataset%buffer)) then
    ! if we have reached the last time, no updates needed
    if (dataset%buffer%cur_time_index == -1) return
    if (.not.Equal(option%time,dataset%buffer%cur_time)) then
      !TODO(geh): modify so that interpolate_dataset is only set to
      !           true if linear interpolation (currently hardwired to step)
      !           or the current time index changes
      ! currently, I am updated every time.  Expensive!
      interpolate_dataset = PETSC_TRUE
      ! increment time index until within buffer
      itime = dataset%buffer%cur_time_index
      do
        if (itime >= dataset%buffer%num_times_total) then
          dataset%buffer%cur_time_index = dataset%buffer%num_times_total
          exit
        endif
        if (option%time >= dataset%buffer%time_array(itime+1)) then
          itime = itime + 1
        else
          dataset%buffer%cur_time_index = itime
          exit
        endif
      enddo
      ! is time outside range of buffer
      if (dataset%buffer%cur_time_index >= dataset%buffer%time_offset + &
                                           dataset%buffer%num_times_in_buffer &
          .and. dataset%buffer%cur_time_index < &
                dataset%buffer%num_times_total) then
        dataset%buffer%time_offset = dataset%buffer%cur_time_index - 1
        read_dataset = PETSC_TRUE
      endif
    endif
  else if (dataset%data_dim == DIM_NULL) then ! has not been read
    read_dataset = PETSC_TRUE
  endif
  
  if (read_dataset) then
#if !defined(PETSC_HAVE_HDF5)
    option%io_buffer = 'DataSetLoad() requires HDF5 support'
    call printErrMsg(option)
  endif
#else

#ifdef PARALLELIO_LIB
      option%io_buffer='In DataLoad: HDF5ReadDataset() not supported with ' // &
        ' PARALLELIO_LIB'
#else
    if(.not.associated(dataset%dataset_map)) then
      call HDF5ReadDataset(dataset,option)
    else
      call HDF5ReadDataset(dataset,option)
      if(dataset%dataset_map%first_time) then
        call HDF5ReadDatasetMap(dataset,option)
      endif
    endif
#endif
    call DatasetReorder(dataset,option)
    if (associated(dataset%buffer)) then
      interpolate_dataset = PETSC_TRUE ! just to be sure
    endif
  endif
  if (interpolate_dataset) then
    call DatasetInterpolateBetweenTimes(dataset,option)
    if (associated(dataset%buffer)) then
      dataset%buffer%cur_time = option%time
    endif
    ! calculate min/max values
    if (associated(dataset%rarray)) then
      dataset%rmax = maxval(dataset%rarray)
      dataset%rmin = minval(dataset%rarray)
    endif
  endif
#endif

end subroutine DatasetLoad

! *********************a***************************************************** !
!
! DatasetProcessDatasets: Determine whether a dataset is indexed by cell ids
! author: Glenn Hammond
! date: 03/26/12
!
! ************************************************************************** !
subroutine DatasetProcessDatasets(datasets,option)

  use Option_module
  
  implicit none
  
  type(dataset_type), pointer :: datasets
  type(option_type) :: option
  
  type(dataset_type), pointer :: cur_dataset
  
  cur_dataset => datasets
  do
    if (.not.associated(cur_dataset)) exit
    cur_dataset%is_cell_indexed = DatasetIsCellIndexed(cur_dataset,option)
    cur_dataset => cur_dataset%next
  enddo
  
end subroutine DatasetProcessDatasets

! ************************************************************************** !
!
! DatasetIsCellIndexed: Determine whether a dataset is indexed by cell ids
! author: Glenn Hammond
! date: 03/26/12
!
! ************************************************************************** !
function DatasetIsCellIndexed(dataset,option)

  use Option_module
  use HDF5_Aux_module

  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  PetscBool :: DatasetIsCellIndexed
  
#if defined(PETSC_HAVE_HDF5)

  DatasetIsCellIndexed = &
    .not.HDF5GroupExists(dataset%filename,dataset%h5_dataset_name,option)

#endif
 
end function DatasetIsCellIndexed

end module Dataset_module
