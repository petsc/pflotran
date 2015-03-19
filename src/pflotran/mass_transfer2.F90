module Mass_Transfer2_module
 
  use Dataset_Global_HDF5_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
 
  type, public :: mass_transfer2_type
    VecScatter :: scatter_ctx ! scatter context from vec to residual_vec
    Vec :: vec
    type(mass_transfer2_type), pointer :: next
  end type mass_transfer2_type
  
  public :: MassTransfer2Create, MassTransfer2Destroy, &
            MassTransfer2AddToList

contains

! ************************************************************************** !

function MassTransfer2Create()
  ! 
  ! Creates a mass transfer object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 
  
  implicit none

  type(mass_transfer2_type), pointer :: MassTransfer2Create
  
  type(mass_transfer2_type), pointer :: mass_transfer
  
  allocate(mass_transfer)
  mass_transfer%vec = 0
  mass_transfer%scatter_ctx = 0
  nullify(mass_transfer%next)
  MassTransfer2Create => mass_transfer

end function MassTransfer2Create

! ************************************************************************** !

recursive subroutine MassTransfer2AddToList(mass_transfer,list)
  ! 
  ! Adds a mass transfer object to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  implicit none
  
  type(mass_transfer2_type), pointer :: mass_transfer
  type(mass_transfer2_type), pointer :: list

  if (associated(list)) then
    call MassTransfer2AddToList(mass_transfer,list%next)
  else
    list => mass_transfer
  endif
  
end subroutine MassTransfer2AddToList

! ************************************************************************** !

recursive subroutine MassTransfer2Destroy(mass_transfer)
  ! 
  ! Destroys a mass transfer object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  implicit none
  
  type(mass_transfer2_type), pointer :: mass_transfer
  
  PetscErrorCode :: ierr
  
  if (.not.associated(mass_transfer)) return
  
  ! Simply nullify the pointer as the dataset resides in a list to be
  ! destroyed separately.
  !nullify(mass_transfer%dataset)
  if (mass_transfer%scatter_ctx /= 0) then
    call VecScatterDestroy(mass_transfer%scatter_ctx,ierr);CHKERRQ(ierr)
  endif
  if (mass_transfer%vec /= 0) then
    call VecDestroy(mass_transfer%vec ,ierr);CHKERRQ(ierr)
  endif
  call MassTransfer2Destroy(mass_transfer%next)

  deallocate(mass_transfer)
  nullify(mass_transfer)
  
end subroutine MassTransfer2Destroy

end module Mass_Transfer2_module
