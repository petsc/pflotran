module Hydrogeophysics_Wrapper_module
 
  implicit none

  private

#include "definitions.h"

  public :: HydrogeophysicsWrapperInit, &
            HydrogeophysicsWrapperStart, &
            HydrogeophysicsWrapperStep, &
            HydrogeophysicsWrapperStop, &
            HydrogeophysicsWrapperDestroy

contains

! ************************************************************************** !
!
! HydrogeophysicsWrapperInit: Initializes the hydrogeophysics module
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperInit(option)
  
  use Option_module
  
  implicit none

  type(option_type) :: option
  
  call printMsg(option,'HydrogeophysicsWrapperInit()')
  
end subroutine HydrogeophysicsWrapperInit

! ************************************************************************** !
!
! HydrogeophysicsWrapperStart: Starts the hydrogeophysics forward simulation 
!                              loop
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStart(option)
  
  use Option_module

  implicit none

  type(option_type) :: option
  
  call printMsg(option,'HydrogeophysicsWrapperStart()')
  
end subroutine HydrogeophysicsWrapperStart

! ************************************************************************** !
!
! HydrogeophysicsWrapperStep: Performs a forward simulation
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStep(sigma,option)
  
  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec :: sigma
  type(option_type) :: option
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call printMsg(option,'HydrogeophysicsWrapperStep()')
  
  call VecGetArrayF90(sigma,vec_ptr,ierr)
  write(option%io_buffer,'(es13.6)') vec_ptr(365)
  call printMsg(option)
  call VecRestoreArrayF90(sigma,vec_ptr,ierr)
  
end subroutine HydrogeophysicsWrapperStep

! ************************************************************************** !
!
! HydrogeophysicsWrapperStop: Stops the hydrogeophysics forward simulation 
!                             loop
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStop(option)
  
  use Option_module

  implicit none

  type(option_type) :: option
  
  call printMsg(option,'HydrogeophysicsWrapperStop()')
  
end subroutine HydrogeophysicsWrapperStop

! ************************************************************************** !
!
! HydrogeophysicsWrapperDestroy: Destroys the contents of the hydrogeophysics 
!                                module
! author: Glenn Hammond
! date: 07/02/13
!
!
! ************************************************************************** !
recursive subroutine HydrogeophysicsWrapperDestroy(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  
  call printMsg(option,'HydrogeophysicsWrapperDestroy()')

end subroutine HydrogeophysicsWrapperDestroy

end module Hydrogeophysics_Wrapper_module
