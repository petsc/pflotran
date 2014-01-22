#ifdef SURFACE_FLOW

module PM_Surface_TH_class

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

  type, public, extends(pm_base_type) :: pm_surface_th_type
    class(surface_realization_type), pointer :: surf_realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMSurfaceTHInit
    procedure, public :: PMSurfaceTHSetRealization
    procedure, public :: InitializeRun => PMSurfaceTHInitializeRun
    procedure, public :: FinalizeRun => PMSurfaceTHFinalizeRun
!    procedure, public :: InitializeTimestep => PMSurfaceTHInitializeTimestep ! not needed for TS
!    procedure, public :: FinalizeTimestep => PMSurfaceTHFinalizeTimeStep
!    procedure, public :: Residual => PMSurfaceTHResidual
!    procedure, public :: Jacobian => PMSurfaceTHJacobian
    procedure, public :: UpdateTimestep => PMSurfaceTHUpdateTimestep
    procedure, public :: PreSolve => PMSurfaceTHPreSolve
    procedure, public :: PostSolve => PMSurfaceTHPostSolve
!    procedure, public :: AcceptSolution => PMSurfaceTHAcceptSolution
!    procedure, public :: CheckUpdatePre => PMSurfaceTHCheckUpdatePre
!    procedure, public :: CheckUpdatePost => PMSurfaceTHCheckUpdatePost
!    procedure, public :: TimeCut => PMSurfaceTHTimeCut
    procedure, public :: UpdateSolution => PMSurfaceTHUpdateSolution
!    procedure, public :: MaxChange => PMSurfaceTHMaxChange
!    procedure, public :: ComputeMassBalance => PMSurfaceTHComputeMassBalance
    procedure, public :: Destroy => PMSurfaceTHDestroy
    procedure, public :: RHSFunction => PMSurfaceTHRHSFunction
    procedure, public :: Checkpoint => PMSurfaceTHCheckpoint
    procedure, public :: Restart => PMSurfaceTHRestart
  end type pm_surface_th_type

  public :: PMSurfaceTHCreate, &
            PMSurfaceTHDTExplicit

contains

! ************************************************************************** !

function PMSurfaceTHCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  implicit none

  class(pm_surface_th_type), pointer :: PMSurfaceTHCreate

  class(pm_surface_th_type), pointer :: surface_th_pm

  allocate(surface_th_pm)
  nullify(surface_th_pm%option)
  nullify(surface_th_pm%output_option)
  nullify(surface_th_pm%surf_realization)
  nullify(surface_th_pm%comm1)

  call PMBaseCreate(surface_th_pm)

  PMSurfaceTHCreate => surface_th_pm

end function PMSurfaceTHCreate

! ************************************************************************** !

subroutine PMSurfaceTHInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Discretization_module
  use Unstructured_Communicator_class
  use Grid_module

  implicit none

  class(pm_surface_th_type) :: this

  ! set up communicator
  select case(this%surf_realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%option%io_buffer='Surface flow not supported on structured grids'
      call printErrMsg(this%option)
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select

  call this%comm1%SetDM(this%surf_realization%discretization%dm_1dof)

end subroutine PMSurfaceTHInit

! ************************************************************************** !

subroutine PMSurfaceTHPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  implicit none

  PetscErrorCode :: ierr
  class(pm_surface_th_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," SURFACE TH FLOW ",62("="))')
  endif

end subroutine PMSurfaceTHPreSolve

! ************************************************************************** !

subroutine PMSurfaceTHSetRealization(this, surf_realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_Realization_class
  use Grid_module

  implicit none

  class(pm_surface_th_type) :: this
  class(surface_realization_type), pointer :: surf_realization

  this%surf_realization => surf_realization
  this%realization_base => surf_realization

  this%solution_vec = surf_realization%surf_field%flow_xx
  this%residual_vec = surf_realization%surf_field%flow_r

end subroutine PMSurfaceTHSetRealization

! ************************************************************************** !

recursive subroutine PMSurfaceTHInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHUpdateSolution

  implicit none

  class(pm_surface_th_type) :: this

end subroutine PMSurfaceTHInitializeRun

! ************************************************************************** !

recursive subroutine PMSurfaceTHFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  implicit none

  class(pm_surface_th_type) :: this

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceTH%FinalizeRun()')
#endif

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMSurfaceTHFinalizeRun

! ************************************************************************** !

subroutine PMSurfaceTHUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHUpdateSolution
  use Condition_module

  implicit none

  class(pm_surface_th_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceTH%UpdateSolution()')
#endif

  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%surf_realization%surf_flow_conditions, &
                           this%surf_realization%option, &
                           this%surf_realization%option%time)

  call SurfRealizAllCouplerAuxVars(this%surf_realization,force_update_flag)

  call SurfaceTHUpdateSolution(this%surf_realization)

end subroutine PMSurfaceTHUpdateSolution

! ************************************************************************** !

subroutine PMSurfaceTHRHSFunction(this,ts,time,xx,ff,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHRHSFunction

  implicit none

  class(pm_surface_th_type) :: this
  TS                                     :: ts
  PetscReal                              :: time
  Vec                                    :: xx
  Vec                                    :: ff
  type(surface_realization_type)         :: surf_realization
  PetscErrorCode                         :: ierr

  call SurfaceTHRHSFunction(ts,time,xx,ff,this%surf_realization,ierr)

end subroutine PMSurfaceTHRHSFunction

! ************************************************************************** !

subroutine PMSurfaceTHUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHComputeMaxDt

  implicit none

  class(pm_surface_th_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceTHComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt = min(0.9d0*dt_max_glb,this%surf_realization%dt_max)

end subroutine PMSurfaceTHUpdateTimestep

! ************************************************************************** !

subroutine PMSurfaceTHDTExplicit(this,dt_max)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHComputeMaxDt

  implicit none

  class(pm_surface_th_type) :: this
  PetscReal :: dt_max

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceTHComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt_max = min(0.9d0*dt_max_glb,this%surf_realization%dt_max)

end subroutine PMSurfaceTHDTExplicit

! ************************************************************************** !

subroutine PMSurfaceTHPostSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Grid_module
  use Discretization_module
  use Surface_Field_module
  use Surface_TH_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pm_surface_th_type) :: this

  PetscReal, pointer :: xx_p(:)
  PetscInt :: local_id
  PetscInt :: istart, iend
  type(surface_field_type), pointer   :: surf_field 
  type(grid_type),pointer             :: surf_grid
  PetscErrorCode :: ierr

  surf_grid => this%surf_realization%discretization%grid
  surf_field => this%surf_realization%surf_field

  ! Ensure evolved solution is +ve
  call VecGetArrayF90(surf_field%flow_xx,xx_p,ierr)
  do local_id = 1,this%surf_realization%discretization%grid%nlmax
    iend = local_id*this%option%nflowdof
    istart = iend - this%option%nflowdof + 1
    if(xx_p(istart) < 1.d-15) then
      xx_p(istart) = 0.d0
      xx_p(iend) = this%option%reference_temperature
    endif
  enddo
  call VecRestoreArrayF90(surf_field%flow_xx,xx_p,ierr)

  ! First, update the solution vector
  call DiscretizationGlobalToLocal(this%surf_realization%discretization, &
          surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! Update aux vars
  call SurfaceTHUpdateTemperature(this%surf_realization)
  call SurfaceTHUpdateAuxVars(this%surf_realization)

end subroutine PMSurfaceTHPostSolve

! ************************************************************************** !

subroutine PMSurfaceTHCheckpoint(this,viewer)
  ! 
  ! This routine checkpoints data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

  use Surface_Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"

  class(pm_surface_th_type) :: this
  PetscViewer :: viewer

  call SurfaceCheckpointProcessModel(viewer,this%surf_realization)

end subroutine PMSurfaceTHCheckpoint

! ************************************************************************** !

subroutine PMSurfaceTHRestart(this,viewer)
  ! 
  ! This routine reads checkpoint data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

  use Surface_Checkpoint_module
  use Surface_TH_module, only : SurfaceTHUpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"

  class(pm_surface_th_type) :: this
  PetscViewer :: viewer

  call SurfaceRestartProcessModel(viewer,this%surf_realization)
  call SurfaceTHUpdateAuxVars(this%surf_realization)
  call this%UpdateSolution()

end subroutine PMSurfaceTHRestart

! ************************************************************************** !

subroutine PMSurfaceTHDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

!  use Surface_TH_module, only : SurfaceTHDestroy

  implicit none

  class(pm_surface_th_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceTHDestroy()')
#endif

#ifndef SIMPLIFY
!  call SurfaceTHDestroy(this%surf_realization)
#endif
  call this%comm1%Destroy()

end subroutine PMSurfaceTHDestroy

end module PM_Surface_TH_class

#endif