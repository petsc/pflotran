module Process_Model_Coupler_module

  use Process_Model_Base_class

  implicit none

#include "definitions.h"
  
  private

  type, public :: process_model_coupler_type
  
    PetscInt :: variable
    PetscReal, pointer :: array(:)
    
    type(process_model_coupler_type), pointer :: next

  end type process_model_coupler_type
  
  ! add interface definitions here
! interface ProcessModelCouplerCreate
!   module procedure ProcessModelCouplerCreate1
!   module procedure ProcessModelCouplerCreate2
! end interface
  
  public :: ProcessModelCouplerCreate, &
            ProcessModelCouplerDestroy
  
contains

! ************************************************************************** !
!
! ProcessModelCouplerCreate: Allocates and initializes a new XXX process_model_coupler
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function ProcessModelCouplerCreate()

  implicit none
  
  type(process_model_coupler_type), pointer :: ProcessModelCouplerCreate
  
  type(process_model_coupler_type), pointer :: process_model_coupler
  
  allocate(process_model_coupler)
  nullify(process_model_coupler%next)
  
  ProcessModelCouplerCreate => process_model_coupler  
  
end function ProcessModelCouplerCreate

! ************************************************************************** !
!
! ProcessModelCouplerDestroy: Deallocates an XXX process_model_coupler
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine ProcessModelCouplerDestroy(process_model_coupler)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  type(process_model_coupler_type), pointer :: process_model_coupler
  
  if (.not.associated(process_model_coupler)) return

  ! destroy member process_model_couplers
! call MemberProcessModelCouplerDestroy(process_model_coupler%member_process_model_coupler)

  ! deallocate arrays
  call DeallocateArray(process_model_coupler%array)
  
  ! deallocate and nullify the process_model_coupler itself
  deallocate(process_model_coupler)
  nullify(process_model_coupler)
  
end subroutine ProcessModelCouplerDestroy
  
end module Process_Model_Coupler_module
