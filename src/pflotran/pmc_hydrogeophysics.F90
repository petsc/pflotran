module PMC_Hydrogeophysics_class

  use PMC_Base_class
  use Realization_class
  use Option_module

  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type, public, extends(pmc_base_type) :: pmc_hydrogeophysics_type
    class(realization_type), pointer :: realization
    Vec :: sigma
  contains
    procedure, public :: Init => PMCHydrogeophysicsInit
    procedure, public :: InitializeRun => PMCHydrogeophysicsInitializeRun
    procedure, public :: RunToTime => PMCHydrogeophysicsRunToTime
    procedure, public :: FinalizeRun => PMCHydrogeophysicsFinalizeRun
    procedure, public :: Destroy => PMCHydrogeophysicsDestroy
  end type pmc_hydrogeophysics_type
  
  public :: PMCHydrogeophysicsCreate
  
contains

! ************************************************************************** !
!
! PMCHydrogeophysicsCreate: Allocates and initializes a new  
!                           process_model_coupler object.
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
function PMCHydrogeophysicsCreate()

  implicit none
  
  class(pmc_hydrogeophysics_type), pointer :: PMCHydrogeophysicsCreate
  
  class(pmc_hydrogeophysics_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  
  PMCHydrogeophysicsCreate => pmc  
  
end function PMCHydrogeophysicsCreate

! ************************************************************************** !
!
! PMCHydrogeophysicsInit: Initializes a new process model coupler object.
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine PMCHydrogeophysicsInit(this)

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call PMCBaseInit(this)
  this%name = 'PMCHydrogeophysics'
  nullify(this%realization) 
  this%sigma = 0
  this%Synchronize1 => PMCHydrogeophysicsSynchronize

end subroutine PMCHydrogeophysicsInit


! ************************************************************************** !
!
! PMCHydrogeophysicsInitializeRun: Initializes the time stepping
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
recursive subroutine PMCHydrogeophysicsInitializeRun(this)

  use Hydrogeophysics_Wrapper_module

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%InitializeRun()')
  
  call HydrogeophysicsWrapperStart(this%option)
  
  if (associated(this%below)) then
    call this%below%InitializeRun()
  endif
  
  if (associated(this%next)) then
    call this%next%InitializeRun()
  endif

end subroutine PMCHydrogeophysicsInitializeRun

! ************************************************************************** !
!
! PMCHydrogeophysicsRunToTime: Runs the actual simulation.
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
recursive subroutine PMCHydrogeophysicsRunToTime(this,sync_time,stop_flag)

  use Hydrogeophysics_Wrapper_module

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  PetscInt :: local_stop_flag
  
  this%option%io_buffer = trim(this%name)
  call printVerboseMsg(this%option)
  
  if (associated(this%Synchronize1)) then
    call this%Synchronize1()
  endif
  
  local_stop_flag = 0
  
  call HydrogeophysicsWrapperStep(this%sigma,this%option)

  ! Run neighboring process model couplers
  if (associated(this%below)) then
    call this%below%RunToTime(sync_time,local_stop_flag)
  endif

  ! Run neighboring process model couplers
  if (associated(this%next)) then
    call this%next%RunToTime(sync_time,local_stop_flag)
  endif

  stop_flag = max(stop_flag,local_stop_flag)  
  
end subroutine PMCHydrogeophysicsRunToTime

! ************************************************************************** !
!
! PMCHydrogeophysicsSynchronize: Runs the actual simulation.
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine PMCHydrogeophysicsSynchronize(this)

  use Realization_Base_class, only : RealizationGetDataset
  use Variables_module, only : PRIMARY_MOLALITY

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_base_type), pointer :: this

  PetscErrorCode :: ierr

  select type(pmc => this)
    class is(pmc_hydrogeophysics_type)
      call RealizationGetDataset(pmc%realization,pmc%realization%field%work, &
                                 PRIMARY_MOLALITY,ONE_INTEGER,0)
      call VecCopy(pmc%realization%field%work,pmc%sigma,ierr)
  end select
  
end subroutine PMCHydrogeophysicsSynchronize

! ************************************************************************** !
!
! PMCHydrogeophysicsFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
recursive subroutine PMCHydrogeophysicsFinalizeRun(this)

  use Hydrogeophysics_Wrapper_module
  
  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%FinalizeRun()')
  
  call HydrogeophysicsWrapperStop(this%option)
  
end subroutine PMCHydrogeophysicsFinalizeRun


! ************************************************************************** !
!
! PMCHydrogeophysicsDestroy: Deallocates a process_model_coupler object
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
recursive subroutine PMCHydrogeophysicsDestroy(this)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%Destroy()')
  
  nullify(this%realization)
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 
  
end subroutine PMCHydrogeophysicsDestroy
  
end module PMC_Hydrogeophysics_class
