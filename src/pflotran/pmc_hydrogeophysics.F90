module PMC_Hydrogeophysics_class

  use PMC_Base_class
  use Option_module

  implicit none

#include "definitions.h"
  
  private

  type, public, extends(pmc_base_type) :: pmc_hydrogeophysics_type
  contains
    procedure, public :: Init => PMCHydrogeophysicsInit
  end type pmc_hydrogeophysics_type
  
  public :: PMCHydrogeophysicsCreate
  
contains

! ************************************************************************** !
!
! PMCHydrogeophysicsCreate: Allocates and initializes a new  
!                           process_model_coupler object.
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function PMCHydrogeophysicsCreate()

  implicit none
  
  class(pmc_hydrogeophysics_type), pointer :: PMCHydrogeophysicsCreate
  
  class(pmc_hydrogeophysics_type), pointer :: pmc

  print *, 'PMCHydrogeophysics%Create()'
  
  allocate(pmc)
  call pmc%Init()
  
  PMCHydrogeophysicsCreate => pmc  
  
end function PMCHydrogeophysicsCreate

! ************************************************************************** !
!
! PMCHydrogeophysicsInit: Initializes a new process model coupler object.
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
subroutine PMCHydrogeophysicsInit(this)

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  print *, 'PMCHydrogeophysics%Init()'
  
  call PMCInit(this)

end subroutine PMCHydrogeophysicsInit

! ************************************************************************** !
!
! PMCHydrogeophysicsRunToTime: Runs the actual simulation.
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCHydrogeophysicsRunToTime(this,sync_time,stop_flag)

  use Output_module, only : Output
  
  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  
end subroutine PMCHydrogeophysicsRunToTime

! ************************************************************************** !
!
! PMCHydrogeophysicsFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCHydrogeophysicsFinalizeRun(this)

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%FinalizeRun()')
  
end subroutine PMCHydrogeophysicsFinalizeRun


! ************************************************************************** !
!
! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
recursive subroutine Destroy(this)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%Destroy()')
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 
  
end subroutine Destroy
  
end module PMC_Hydrogeophysics_class
