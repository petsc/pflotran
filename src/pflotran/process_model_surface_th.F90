#ifdef SURFACE_FLOW

module Process_Model_Surface_TH_class

  use Process_Model_Base_class
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
  end type pm_surface_th_type

  public :: PMSurfaceTHCreate, &
            PMSurfaceTHDTExplicit

contains

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
function PMSurfaceTHCreate()

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHInit(this)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHPreSolve(this)

  implicit none

  PetscErrorCode :: ierr
  class(pm_surface_th_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," SURFACE FLOW ",62("="))')
  endif

  call VecView(this%surf_realization%surf_field%flow_xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call printErrMsg(this%option,'PreSolve')
  

end subroutine PMSurfaceTHPreSolve

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHSetRealization(this, surf_realization)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
recursive subroutine PMSurfaceTHInitializeRun(this)

  use Surface_TH_module, only : SurfaceTHUpdateSolution

  implicit none

  class(pm_surface_th_type) :: this

end subroutine PMSurfaceTHInitializeRun

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
recursive subroutine PMSurfaceTHFinalizeRun(this)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHUpdateSolution(this)

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

  call SurfaceTHUpdateSolution(this%surf_realization)

end subroutine PMSurfaceTHUpdateSolution

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHRHSFunction(this,ts,time,xx,ff,ierr)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHDTExplicit(this,dt_max)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHPostSolve(this)

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
      xx_p(iend) = 0.d0
    endif
  enddo
  call VecRestoreArrayF90(surf_field%flow_xx,xx_p,ierr)


  ! First, update the solution vector
  call DiscretizationGlobalToLocal(this%surf_realization%discretization, &
          surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)
  call VecView(surf_field%flow_xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

  ! Update aux vars
  call SurfaceTHUpdateTemperature(this%surf_realization)
  call SurfaceTHUpdateAuxVars(this%surf_realization)

end subroutine PMSurfaceTHPostSolve


! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/23/13
! ************************************************************************** !
subroutine PMSurfaceTHDestroy(this)

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

end module Process_Model_Surface_TH_class

#endif