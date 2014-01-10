module Mass_Transfer_module
 
  use Dataset_Global_HDF5_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
 
  type, public :: mass_transfer_type
    PetscInt :: idof
    character(len=MAXWORDLENGTH) :: name
    class(dataset_global_hdf5_type), pointer :: dataset
    Vec :: vec
    type(mass_transfer_type), pointer :: next
  end type mass_transfer_type
  
  public :: MassTransferCreate, MassTransferDestroy, &
            MassTransferRead, MassTransferAddToList, &
            MassTransferUpdate, MassTransferInit

contains

! ************************************************************************** !

function MassTransferCreate()
  ! 
  ! Creates a mass transfer object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 
  
  implicit none

  type(mass_transfer_type), pointer :: MassTransferCreate
  
  type(mass_transfer_type), pointer :: mass_transfer
  
  allocate(mass_transfer)
  mass_transfer%idof = 0
  mass_transfer%name = ''
  nullify(mass_transfer%dataset)
  nullify(mass_transfer%next)
  mass_transfer%vec = 0
  MassTransferCreate => mass_transfer

end function MassTransferCreate

! ************************************************************************** !

subroutine MassTransferRead(mass_transfer,input,option)
  ! 
  ! Reads in contents of a mass transfer card
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  type(mass_transfer_type) :: mass_transfer
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','MASS_TRANSFER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('IDOF') 
        call InputReadInt(input,option,mass_transfer%idof)
        call InputErrorMsg(input,option,'idof','MASS_TRANSFER')
      case('DATASET')
        mass_transfer%dataset => DatasetGlobalHDF5Create()
        call InputReadNChars(input,option, &
                             mass_transfer%dataset%name,&
                             MAXWORDLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'DATASET,NAME','MASS_TRANSFER')
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in mass transfer'    
        call printErrMsg(option)
    end select
    
  enddo  

end subroutine MassTransferRead

! ************************************************************************** !

subroutine MassTransferAddToList(mass_transfer,list)
  ! 
  ! Adds a mass transfer object to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  implicit none
  
  type(mass_transfer_type), pointer :: mass_transfer
  type(mass_transfer_type), pointer :: list

  type(mass_transfer_type), pointer :: cur_mass_transfer
  
  if (associated(list)) then
    cur_mass_transfer => list
    ! loop to end of list
    do
      if (.not.associated(cur_mass_transfer%next)) exit
      cur_mass_transfer => cur_mass_transfer%next
    enddo
    cur_mass_transfer%next => mass_transfer
  else
    list => mass_transfer
  endif
  
end subroutine MassTransferAddToList

! ************************************************************************** !

recursive subroutine MassTransferInit(mass_transfer, discretization, &
                                      available_datasets, option)
  ! 
  ! Initializes mass transfer object opening dataset to
  ! set up times, vectors, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/09/13
  ! 

  use Discretization_module
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"  
  
  type(mass_transfer_type), pointer :: mass_transfer
  type(discretization_type) :: discretization
  class(dataset_base_type), pointer :: available_datasets
  type(option_type) :: option
  
  class(dataset_base_type), pointer :: dataset_base_ptr
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  if (.not.associated(mass_transfer)) return
  
  if (.not.associated(mass_transfer%dataset)) then
    option%io_buffer = 'A "global" DATASET does not exist for ' // &
      'MASS_TRANSFER object "' // trim(mass_transfer%name) // '".'
    call printErrMsg(option)
  endif

  string = 'Mass Transfer ' // trim(mass_transfer%name)
  dataset_base_ptr => &
    DatasetBaseGetPointer(available_datasets,mass_transfer%dataset%name, &
                          string,option)
  call DatasetGlobalHDF5Destroy(mass_transfer%dataset)
  select type(dataset => dataset_base_ptr)
    class is(dataset_global_hdf5_type)
      mass_transfer%dataset => dataset
    class default
      option%io_buffer = 'DATASET ' // trim(dataset%name) // 'is not of ' // &
        'GLOBAL type, which is necessary for all MASS_TRANSFER objects.'
      call printErrMsg(option)
  end select
  ! dm_wrapper is solely a pointer; it should not be allocated
  mass_transfer%dataset%dm_wrapper => discretization%dm_1dof
  mass_transfer%dataset%local_size = discretization%grid%nlmax
  mass_transfer%dataset%global_size = discretization%grid%nmax
  call DiscretizationCreateVector(discretization,ONEDOF,mass_transfer%vec, &
                                    GLOBAL,option)
  call VecZeroEntries(mass_transfer%vec,ierr)    
  
  if (.not.associated(mass_transfer%dataset%time_storage)) then
#if defined(PETSC_HAVE_HDF5)    
    call DatasetCommonHDF5ReadTimes(mass_transfer%dataset%filename, &
                                    mass_transfer%dataset%hdf5_dataset_name, &
                                    mass_transfer%dataset%time_storage,option)
#endif
    ! if time interpolation methods not set in hdf5 file, set to default of STEP
    if (mass_transfer%dataset%time_storage%time_interpolation_method == &
        INTERPOLATION_NULL) then
      mass_transfer%dataset%time_storage%time_interpolation_method = &
        INTERPOLATION_STEP
    endif
  endif 
  
  ! update the next one recursively
  if (associated(mass_transfer%next)) then
    call MassTransferInit(mass_transfer%next,discretization, &
                          available_datasets,option)
  endif  
  
end subroutine MassTransferInit

! ************************************************************************** !

recursive subroutine MassTransferUpdate(mass_transfer, grid, option)
  ! 
  ! Updates a mass transfer object transfering data from
  ! the buffer into the PETSc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  use Discretization_module
  use Grid_module
  use Option_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"  
  
  type(mass_transfer_type), pointer :: mass_transfer
  type(grid_type) :: grid
  type(option_type) :: option  
  PetscReal :: time
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  if (.not.associated(mass_transfer)) return

  call DatasetGlobalHDF5Load(mass_transfer%dataset,option)

  call VecGetArrayF90(mass_transfer%vec,vec_ptr,ierr)
  ! multiply by -1.d0 for positive contribution to residual
  vec_ptr(:) = -1.d0*mass_transfer%dataset%rarray(:)
  call VecRestoreArrayF90(mass_transfer%vec,vec_ptr,ierr)
  
  ! update the next one
  if (associated(mass_transfer%next)) then
    call MassTransferUpdate(mass_transfer%next,grid,option)
  endif
  
end subroutine MassTransferUpdate

! ************************************************************************** !

recursive subroutine MassTransferDestroy(mass_transfer)
  ! 
  ! Destroys a mass transfer object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  implicit none
  
  type(mass_transfer_type), pointer :: mass_transfer
  
  PetscErrorCode :: ierr
  
  if (.not.associated(mass_transfer)) return
  
  ! Simply nullify the pointer as the dataset resides in a list to be
  ! destroyed separately.
  nullify(mass_transfer%dataset)
  if (mass_transfer%vec /= 0) &
    call VecDestroy(mass_transfer%vec ,ierr)
  call MassTransferDestroy(mass_transfer%next)

  deallocate(mass_transfer)
  nullify(mass_transfer)
  
end subroutine MassTransferDestroy

end module Mass_Transfer_module
