module Dataset_module
 
  use Dataset_Aux_module
  
  implicit none

  private
  
#include "definitions.h"

  public :: DatasetLoad

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
  use HDF5_aux_module
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
  
  
 ! Glenn:  Need to setup waypoints to match times in dataset, like done in timeseries
  
  
  
    if (.not.Equal(option%time,dataset%buffer%cur_time)) then
      interpolate_dataset = PETSC_TRUE
      ! increment time index until within buffer
      do itime = dataset%buffer%cur_time_index, dataset%buffer%num_times_total
        if (option%time < dataset%buffer%time_array(itime)) then
        else
          dataset%buffer%cur_time_index = itime
          exit
        endif
      enddo
      ! is time outside range of buffer
      if (dataset%buffer%cur_time_index > dataset%buffer%time_offset + &
                                       dataset%buffer%num_times_in_buffer) then
        read_dataset = PETSC_TRUE
      endif
    endif
  else if (dataset%data_dim == DIM_NULL) then ! has not been read
    read_dataset = PETSC_TRUE
    interpolate_dataset = PETSC_TRUE
  endif
  
  if (read_dataset) then
    call HDF5ReadDataset(dataset,option)
    call DatasetReorder(dataset,option)
    interpolate_dataset = PETSC_TRUE ! just to be sure
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

end subroutine DatasetLoad

end module Dataset_module
