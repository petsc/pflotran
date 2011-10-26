module Dataset_module
 
  use Dataset_Aux_module
  
  implicit none

  private

  public :: DatasetLoad

contains

subroutine DatasetLoad(dataset,option)

  use Option_module
  use HDF5_aux_module

  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  call HDF5ReadDataset(dataset,option)

end subroutine DatasetLoad

end module Dataset_module
