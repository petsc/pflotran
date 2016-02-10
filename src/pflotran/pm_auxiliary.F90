module PM_Auxiliary_class

  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(pm_base_type) :: pm_auxiliary_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    procedure(PMAuxliaryEvaluate), pointer :: Evaluate => null()
  contains
  end type pm_auxiliary_type
  
  ! interface blocks
  interface
    subroutine PMAuxliaryEvaluate(this,time,ierr)
      import :: pm_auxiliary_type
      implicit none
      class(pm_auxiliary_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr      
    end subroutine PMAuxliaryEvaluate
  end interface
  
  public :: PMAuxiliaryInit
  
contains

! ************************************************************************** !

subroutine PMAuxiliaryInit(this)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this  

  nullify(this%realization)
  nullify(this%comm1)
  
  call PMBaseInit(this)
  
end subroutine PMAuxiliaryInit

! ************************************************************************** !

subroutine PMAuxiliarySetFunctionPointer(this,string)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Option_module
  
  implicit none
  
  class(pm_auxiliary_type) :: this
  character(len=MAXSTRINGLENGTH) :: string

  select case(string)
    case('A')
      this%Evaluate => PMAuxiliaryEvaluateA
    case('B')
      this%Evaluate => PMAuxiliaryEvaluateb
    case default
      this%option%io_buffer = 'Function pointer "' // trim(string) // '" not &
        &found among available functions in PMAuxiliarySetFunctionPointer.'
      call printErrMsg(this%option)
  end select
  
end subroutine PMAuxiliarySetFunctionPointer

! ************************************************************************** !

subroutine PMAuxiliaryEvaluateA(this,time,ierr)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

 ! call InitSubsurfAssignMatIDsToRegns(this%realization)
!  call InitSubsurfAssignMatProperties(this%realization)
  
end subroutine PMAuxiliaryEvaluateA

! ************************************************************************** !

subroutine PMAuxiliaryEvaluateB(this,time,ierr)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  
  
end subroutine PMAuxiliaryEvaluateB

! ************************************************************************** !

subroutine PMAuxiliaryDestroy(this)
  ! 
  ! Destroys auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this
  
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  
end subroutine PMAuxiliaryDestroy

end module PM_Auxiliary_class
