#ifdef SURFACE_FLOW

module PM_Surface_Flow_class

  use PM_Base_class
!geh: using Richards_module here fails with gfortran (internal compiler error)
!  use Richards_module
  use Surface_Realization_class
  use Realization_class
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

  type, public, extends(pm_base_type) :: pm_surface_flow_type
    class(surface_realization_type), pointer :: surf_realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMSurfaceFlowInit
    procedure, public :: PMSurfaceFlowSetRealization
    procedure, public :: InitializeRun => PMSurfaceFlowInitializeRun
    procedure, public :: FinalizeRun => PMSurfaceFlowFinalizeRun
!    procedure, public :: InitializeTimestep => PMSurfaceFlowInitializeTimestep ! not needed for TS
!    procedure, public :: FinalizeTimestep => PMSurfaceFlowFinalizeTimeStep
!    procedure, public :: Residual => PMSurfaceFlowResidual
!    procedure, public :: Jacobian => PMSurfaceFlowJacobian
    procedure, public :: UpdateTimestep => PMSurfaceFlowUpdateTimestep
    procedure, public :: PreSolve => PMSurfaceFlowPreSolve
    procedure, public :: PostSolve => PMSurfaceFlowPostSolve
!    procedure, public :: AcceptSolution => PMSurfaceFlowAcceptSolution
!    procedure, public :: CheckUpdatePre => PMSurfaceFlowCheckUpdatePre
!    procedure, public :: CheckUpdatePost => PMSurfaceFlowCheckUpdatePost
!    procedure, public :: TimeCut => PMSurfaceFlowTimeCut
    procedure, public :: UpdateSolution => PMSurfaceFlowUpdateSolution
!    procedure, public :: MaxChange => PMSurfaceFlowMaxChange
!    procedure, public :: ComputeMassBalance => PMSurfaceFlowComputeMassBalance
    procedure, public :: Destroy => PMSurfaceFlowDestroy
    procedure, public :: RHSFunction => PMSurfaceFlowRHSFunction
    procedure, public :: Checkpoint => PMSurfaceCheckpoint
    procedure, public :: Restart => PMSurfaceRestart
  end type pm_surface_flow_type

  public :: PMSurfaceFlowCreate, &
            PMSurfaceFlowDTExplicit

contains

! ************************************************************************** !

function PMSurfaceFlowCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  implicit none

  class(pm_surface_flow_type), pointer :: PMSurfaceFlowCreate

  class(pm_surface_flow_type), pointer :: surface_flow_pm

  allocate(surface_flow_pm)
  nullify(surface_flow_pm%option)
  nullify(surface_flow_pm%output_option)
  nullify(surface_flow_pm%surf_realization)
  nullify(surface_flow_pm%comm1)

  call PMBaseCreate(surface_flow_pm)

  PMSurfaceFlowCreate => surface_flow_pm

end function PMSurfaceFlowCreate

! ************************************************************************** !

subroutine PMSurfaceFlowInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Discretization_module
  use Unstructured_Communicator_class
  use Grid_module

  implicit none

  class(pm_surface_flow_type) :: this

  ! set up communicator
  select case(this%surf_realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%option%io_buffer='Surface flow not supported on structured grids'
      call printErrMsg(this%option)
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select

  call this%comm1%SetDM(this%surf_realization%discretization%dm_1dof)

end subroutine PMSurfaceFlowInit

! ************************************************************************** !

subroutine PMSurfaceFlowPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  implicit none

  class(pm_surface_flow_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," SURFACE FLOW ",62("="))')
  endif

end subroutine PMSurfaceFlowPreSolve

! ************************************************************************** !

subroutine PMSurfaceFlowSetRealization(this, surf_realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Realization_class
  use Grid_module

  implicit none

  class(pm_surface_flow_type) :: this
  class(surface_realization_type), pointer :: surf_realization

  this%surf_realization => surf_realization
  this%realization_base => surf_realization

  this%solution_vec = surf_realization%surf_field%flow_xx
  this%residual_vec = surf_realization%surf_field%flow_r

end subroutine PMSurfaceFlowSetRealization

! ************************************************************************** !

recursive subroutine PMSurfaceFlowInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowUpdateSolution

  implicit none

  class(pm_surface_flow_type) :: this

end subroutine PMSurfaceFlowInitializeRun

! ************************************************************************** !

recursive subroutine PMSurfaceFlowFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  implicit none

  class(pm_surface_flow_type) :: this

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceFlow%FinalizeRun()')
#endif

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMSurfaceFlowFinalizeRun

! ************************************************************************** !

subroutine PMSurfaceFlowUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowUpdateSolution
  use Condition_module

  implicit none

  class(pm_surface_flow_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceFlow%UpdateSolution()')
#endif

  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%surf_realization%surf_flow_conditions, &
                           this%surf_realization%option, &
                           this%surf_realization%option%time)

  call SurfRealizAllCouplerAuxVars(this%surf_realization,force_update_flag)

  call SurfaceFlowUpdateSolution(this%surf_realization)

end subroutine PMSurfaceFlowUpdateSolution

! ************************************************************************** !

subroutine PMSurfaceFlowRHSFunction(this,ts,time,xx,ff,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowRHSFunction

  implicit none

  class(pm_surface_flow_type) :: this
  TS                                     :: ts
  PetscReal                              :: time
  Vec                                    :: xx
  Vec                                    :: ff
  type(surface_realization_type)         :: surf_realization
  PetscErrorCode                         :: ierr

  call SurfaceFlowRHSFunction(ts,time,xx,ff,this%surf_realization,ierr)

end subroutine PMSurfaceFlowRHSFunction

! ************************************************************************** !

subroutine PMSurfaceFlowUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowComputeMaxDt

  implicit none

  class(pm_surface_flow_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceFlowComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt = min(0.9d0*dt_max_glb,this%surf_realization%dt_max)

end subroutine PMSurfaceFlowUpdateTimestep

! ************************************************************************** !

subroutine PMSurfaceFlowDTExplicit(this,dt_max)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowComputeMaxDt

  implicit none

  class(pm_surface_flow_type) :: this
  PetscReal :: dt_max

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceFlowComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt_max = min(0.9d0*dt_max_glb,this%surf_realization%dt_max)

end subroutine PMSurfaceFlowDTExplicit

! ************************************************************************** !

subroutine PMSurfaceFlowPostSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Discretization_module
  use Surface_Field_module
  use Surface_Flow_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pm_surface_flow_type) :: this

  PetscReal, pointer :: xx_p(:)
  PetscInt :: local_id
  type(surface_field_type), pointer   :: surf_field 
  PetscErrorCode :: ierr

  surf_field => this%surf_realization%surf_field

  ! Ensure evolved solution is +ve
  call VecGetArrayF90(surf_field%flow_xx,xx_p,ierr)
  do local_id = 1,this%surf_realization%discretization%grid%nlmax
    if(xx_p(local_id)<1.d-15) xx_p(local_id) = 0.d0
  enddo
  call VecRestoreArrayF90(surf_field%flow_xx,xx_p,ierr)

  ! First, update the solution vector
  call DiscretizationGlobalToLocal(this%surf_realization%discretization, &
          surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! Update aux vars
  call SurfaceFlowUpdateAuxVars(this%surf_realization)

end subroutine PMSurfaceFlowPostSolve

! ************************************************************************** !

subroutine PMSurfaceCheckpoint(this,viewer)
  ! 
  ! This routine checkpoints data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

  use Surface_Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"

  class(pm_surface_flow_type) :: this
  PetscViewer :: viewer

  call SurfaceCheckpointProcessModel(viewer,this%surf_realization)

end subroutine PMSurfaceCheckpoint

! ************************************************************************** !

subroutine PMSurfaceRestart(this,viewer)
  ! 
  ! This routine reads checkpoint data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

  use Surface_Checkpoint_module
  use Surface_Flow_module, only : SurfaceFlowUpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"

  class(pm_surface_flow_type) :: this
  PetscViewer :: viewer

  call SurfaceRestartProcessModel(viewer,this%surf_realization)
  call SurfaceFlowUpdateAuxVars(this%surf_realization)
  call this%UpdateSolution()

end subroutine PMSurfaceRestart

! ************************************************************************** !

subroutine PMSurfaceFlowDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

!  use Surface_Flow_module, only : SurfaceFlowDestroy

  implicit none

  class(pm_surface_flow_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceFlowDestroy()')
#endif

#ifndef SIMPLIFY
!  call SurfaceFlowDestroy(this%surf_realization)
#endif
  call this%comm1%Destroy()

end subroutine PMSurfaceFlowDestroy

end module PM_Surface_Flow_class

#endif