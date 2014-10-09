module PMC_Material_class

  use PMC_Base_class
  use Realization_class

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
  
  private

  type, public, extends(pmc_base_type) :: pmc_material_type
    class(realization_type), pointer :: realization
  contains
    procedure, public :: Init => PMCMaterialInit
    procedure, public :: RunToTime => PMCMaterialRunToTime
    procedure, public :: Destroy => PMCMaterialDestroy
  end type pmc_material_type
  
  public :: PMCMaterialCreate
  
contains

! ************************************************************************** !

function PMCMaterialCreate()
  ! 
  ! Allocates and initializes a new process_model_coupler
  ! object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pmc_material_type), pointer :: PMCMaterialCreate
  
  class(pmc_material_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCMaterial%Create()'
#endif
  
  allocate(pmc)
  call pmc%Init()
  
  PMCMaterialCreate => pmc  
  
end function PMCMaterialCreate

! ************************************************************************** !

subroutine PMCMaterialInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_material_type) :: this
  
#ifdef DEBUG
  print *, 'PMCMaterial%Init()'
#endif
  
  call PMCBaseInit(this)
  this%name = 'PMCMaterial'
  nullify(this%realization)

end subroutine PMCMaterialInit

! ************************************************************************** !

recursive subroutine PMCMaterialRunToTime(this,sync_time,stop_flag)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_Base_class
  use Subsurface_module
  use Option_module
  
  implicit none
  
#include "finclude/petscviewer.h"  

  class(pmc_material_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  PetscInt :: local_stop_flag
  PetscErrorCode :: ierr
  
  if (this%stage /= 0) then
    call PetscLogStagePush(this%stage,ierr);CHKERRQ(ierr)
  endif
  this%option%io_buffer = trim(this%name)
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  call this%GetAuxData()
  
  local_stop_flag = TS_CONTINUE
  ! if at end of simulation, skip update of material properties
  if (stop_flag /= TS_STOP_END_SIMULATION) then
    call SubsurfAssignMatIDsToRegions(this%realization)
    call SubsurfAssignMaterialProperties(this%realization)
  endif
  
  ! Run underlying process model couplers
  if (associated(this%below)) then
    ! Set data needed by process-models
    call this%SetAuxData()
    call this%below%RunToTime(this%timestepper%target_time,local_stop_flag)
    ! Get data from other process-models
    call this%GetAuxData()
  endif

  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%next)) then
    call this%next%RunToTime(sync_time,local_stop_flag)
  endif
  
  stop_flag = max(stop_flag,local_stop_flag)
  
  if (this%stage /= 0) then
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMCMaterialRunToTime

! ************************************************************************** !
!
! PMCMaterialFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCMaterialFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_material_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCMaterial%FinalizeRun()')
#endif
  
  nullify(this%realization)
  
end subroutine PMCMaterialFinalizeRun

! ************************************************************************** !

subroutine PMCMaterialStrip(this)
  !
  ! Deallocates members of PMC Material.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_material_type) :: this

  call PMCBaseStrip(this)
  nullify(this%realization)

end subroutine PMCMaterialStrip

! ************************************************************************** !

recursive subroutine PMCMaterialDestroy(this)
  ! 
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pmc_material_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCMaterial%Destroy()')
#endif

  call PMCMaterialStrip(this)
  
  if (associated(this%below)) then
    call this%below%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%below)
    nullify(this%below)
  endif 
  
  if (associated(this%next)) then
    call this%next%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%next)
    nullify(this%next)
  endif
  
end subroutine PMCMaterialDestroy
  
end module PMC_Material_class
