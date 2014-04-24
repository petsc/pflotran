module PM_Surface_class

  use PM_Base_class
  use Surface_Realization_class
  use Communicator_Base_module
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscts.h"

  type, public, extends(pm_base_type) :: pm_surface_type
    class(surface_realization_type), pointer :: surf_realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMSurfaceInit
    procedure, public :: PMSurfaceSetRealization
    procedure, public :: InitializeRun => PMSurfaceInitializeRun
    procedure, public :: PreSolve => PMSurfacePreSolve
    procedure, public :: PostSolve => PMSurfacePostSolve
    procedure, public :: Checkpoint => PMSurfaceCheckpoint
    procedure, public :: Restart => PMSurfaceRestart
    procedure, public :: UpdateAuxvars => PMSurfaceUpdateAuxvars
  end type pm_surface_type

  public :: PMSurfaceCreate, &
            PMSurfaceInit, &
            PMSurfaceUpdateSolution, &
            PMSurfaceDestroy
  
contains

! ************************************************************************** !

subroutine PMSurfaceCreate(this)
  ! 
  ! Intializes shared members of surface process models
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14

  implicit none
  
  class(pm_surface_type) :: this
  
  nullify(this%surf_realization)
  nullify(this%comm1)
  
  call PMBaseCreate(this)

end subroutine PMSurfaceCreate

! ************************************************************************** !

subroutine PMSurfaceInit(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Discretization_module
  use Unstructured_Communicator_class
  use Grid_module

  implicit none

  class(pm_surface_type) :: this

  ! set up communicator
  select case(this%surf_realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%option%io_buffer='Surface flow not supported on structured grids'
      call printErrMsg(this%option)
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select

  ! set the communicator
  call this%comm1%SetDM(this%surf_realization%discretization%dm_1dof)

end subroutine PMSurfaceInit

! ************************************************************************** !

subroutine PMSurfaceSetRealization(this, surf_realization)
  ! 
  ! Initializes relization and PETSc vectors for solution and residual.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Surface_Realization_class
  use Grid_module

  implicit none

  class(pm_surface_type) :: this
  class(surface_realization_type), pointer :: surf_realization

  this%surf_realization => surf_realization
  this%realization_base => surf_realization

  this%solution_vec = surf_realization%surf_field%flow_xx
  this%residual_vec = surf_realization%surf_field%flow_r

end subroutine PMSurfaceSetRealization

! ************************************************************************** !

recursive subroutine PMSurfaceInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  !

  implicit none

  class(pm_surface_type) :: this

end subroutine PMSurfaceInitializeRun

! ************************************************************************** !
subroutine PMSurfacePreSolve(this)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14

  use Global_module

  implicit none
  
  class(pm_surface_type) :: this
  
  this%option%io_buffer = 'PMSurfacePreSolve() must be extended.'
  call printErrMsg(this%option)  

end subroutine PMSurfacePreSolve

! ************************************************************************** !

subroutine PMSurfacePostSolve(this)
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Global_module

  implicit none
  
  class(pm_surface_type) :: this
  
  this%option%io_buffer = 'PMSurfacePostSolve() must be extended.'
  call printErrMsg(this%option)  
  
end subroutine PMSurfacePostSolve

! ************************************************************************** !

subroutine PMSurfaceUpdateSolution(this)
  !
  ! As a first step in updating the solution, update all flow-conditions.
  ! The solution will be updated by each child class of pm_surface_type.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Condition_module

  implicit none

  class(pm_surface_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE


  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%surf_realization%surf_flow_conditions, &
                           this%surf_realization%option, &
                           this%surf_realization%option%time)

  call SurfRealizAllCouplerAuxVars(this%surf_realization,force_update_flag)

end subroutine PMSurfaceUpdateSolution

! ************************************************************************** !

subroutine PMSurfaceUpdateAuxvars(this)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14

  implicit none
  
  class(pm_surface_type) :: this

  this%option%io_buffer = 'PMSurfaceUpdateAuxvars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMSurfaceUpdateAuxvars

! ************************************************************************** !

subroutine PMSurfaceCheckpoint(this,viewer)
  ! 
  ! This routine checkpoints data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Surface_Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"

  class(pm_surface_type) :: this
  PetscViewer :: viewer

  call SurfaceCheckpointProcessModel(viewer,this%surf_realization)

end subroutine PMSurfaceCheckpoint

! ************************************************************************** !

subroutine PMSurfaceRestart(this,viewer)
  ! 
  ! This routine reads checkpoint data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Surface_Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"

  class(pm_surface_type) :: this
  PetscViewer :: viewer

  call SurfaceRestartProcessModel(viewer,this%surf_realization)
  call this%UpdateAuxVars()
  call this%UpdateSolution()

end subroutine PMSurfaceRestart

! ************************************************************************** !

recursive subroutine PMSurfaceFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  implicit none

  class(pm_surface_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMSurfaceFinalizeRun

! ************************************************************************** !

subroutine PMSurfaceDestroy(this)
  ! 
  ! Destroys Surface process model
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  implicit none

  class(pm_surface_type) :: this

  call this%comm1%Destroy()

end subroutine PMSurfaceDestroy

end module PM_Surface_class
