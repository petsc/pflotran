#ifdef GEOMECH

module PM_Geomechanics_Force_class

  use PM_Base_class
#ifdef PROCESS_MODEL
  use Geomechanics_Realization_class
#else
  use Geomechanics_Realization_module
#endif
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

  type, public, extends(pm_base_type) :: pm_geomech_force_type
    class(geomech_realization_type), pointer :: geomech_realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMGeomechForceInit
    procedure, public :: PMGeomechForceSetRealization
    procedure, public :: InitializeRun => PMGeomechForceInitializeRun
    procedure, public :: FinalizeRun => PMGeomechForceFinalizeRun
    procedure, public :: InitializeTimestep => PMGeomechForceInitializeTimestep
    procedure, public :: Residual => PMGeomechForceResidual
    procedure, public :: Jacobian => PMGeomechForceJacobian
    procedure, public :: PreSolve => PMGeomechForcePreSolve
    procedure, public :: UpdateSolution => PMGeomechForceUpdateSolution
    procedure, public :: Checkpoint => PMGeomechForceCheckpoint
    procedure, public :: Restart => PMGeomechForceRestart
    procedure, public :: Destroy => PMGeomechForceDestroy
    procedure, public :: FinalizeTimestep => PMGeomechForceFinalizeTimestep
  end type pm_geomech_force_type

  public :: PMGeomechForceCreate

contains

! ************************************************************************** !
!> This routine creates
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
function PMGeomechForceCreate()

  implicit none

  class(pm_geomech_force_type), pointer :: PMGeomechForceCreate

  class(pm_geomech_force_type), pointer :: geomech_force_pm

  allocate(geomech_force_pm)
  nullify(geomech_force_pm%option)
  nullify(geomech_force_pm%output_option)
  nullify(geomech_force_pm%geomech_realization)
  nullify(geomech_force_pm%comm1)

  call PMBaseCreate(geomech_force_pm)

  PMGeomechForceCreate => geomech_force_pm

end function PMGeomechForceCreate

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceInit(this)

  use Geomechanics_Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module

  implicit none

  class(pm_geomech_force_type) :: this

  ! set up communicator
  select case(this%geomech_realization%geomech_discretization%itype)
    case(STRUCTURED_GRID)
      this%comm1 => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select

  !call this%comm1%SetDM(this%geomech_realization%geomech_discretization%dm_1dof)

end subroutine PMGeomechForceInit

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
recursive subroutine PMGeomechForceInitializeRun(this)

  use Geomechanics_Force_module, only : GeomechUpdateSolution

  implicit none

  class(pm_geomech_force_type) :: this

end subroutine PMGeomechForceInitializeRun

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
recursive subroutine PMGeomechForceFinalizeRun(this)

  implicit none

  class(pm_geomech_force_type) :: this

#ifdef PM_GEOMECH_FORCE_DEBUG
  call printMsg(this%option,'PMGeomechForce%FinalizeRun()')
#endif

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMGeomechForceFinalizeRun

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceSetRealization(this, geomech_realization)

  use Grid_module

  implicit none

  class(pm_geomech_force_type) :: this
  class(geomech_realization_type), pointer :: geomech_realization

  this%geomech_realization => geomech_realization
  this%realization_base => geomech_realization

  this%solution_vec = geomech_realization%geomech_field%disp_xx
  this%residual_vec = geomech_realization%geomech_field%disp_r

end subroutine PMGeomechForceSetRealization

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceInitializeTimestep(this)

  use Geomechanics_Force_module, only : GeomechanicsForceInitialGuess
  use Global_module
  
  implicit none
  
  class(pm_geomech_force_type) :: this

#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%InitializeTimestep()')
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," GEOMECHANICS ",62("="))')
  endif
  
  call GeomechanicsForceInitialGuess(this%geomech_realization)
  
end subroutine PMGeomechForceInitializeTimestep

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceResidual(this,snes,xx,r,ierr)

  use Geomechanics_Force_module, only : GeomechForceResidual

  implicit none
  
  class(pm_geomech_force_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%Residual()')
#endif
  
  call GeomechForceResidual(snes,xx,r,this%geomech_realization,ierr)

end subroutine PMGeomechForceResidual

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceJacobian(this,snes,xx,A,B,flag,ierr)

  use Geomechanics_Force_module, only : GeomechForceJacobian

  implicit none
  
  class(pm_geomech_force_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%Jacobian()')
#endif
  
  call GeomechForceJacobian(snes,xx,A,B,flag,this%geomech_realization,ierr)

end subroutine PMGeomechForceJacobian

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForcePreSolve(this)

  implicit none

  class(pm_geomech_force_type) :: this

end subroutine PMGeomechForcePreSolve

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceUpdateSolution(this)

  use Geomechanics_Force_module, only : GeomechUpdateSolution, &
                                        GeomechStoreInitialDisp, &
                                        GeomechForceUpdateAuxVars
  use Condition_module

  implicit none

  class(pm_geomech_force_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_GEOMECH_FORCE_DEBUG
  call printMsg(this%option,'PMGeomechForce%UpdateSolution()')
#endif

  ! begin from RealizationUpdate()
  call GeomechUpdateSolution(this%geomech_realization)
  call GeomechStoreInitialDisp(this%geomech_realization)
  call GeomechForceUpdateAuxVars(this%geomech_realization)

end subroutine PMGeomechForceUpdateSolution

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceFinalizeTimestep(this)

  use Global_module

  implicit none
  
  class(pm_geomech_force_type) :: this
  
#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%FinalizeTimestep()')
#endif

end subroutine PMGeomechForceFinalizeTimestep

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceCheckpoint(this,viewer)

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_geomech_force_type) :: this
  PetscViewer :: viewer
  
  call printErrMsg(this%option,'add code for checkpointing Geomech in PM approach')
  
end subroutine PMGeomechForceCheckpoint

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceRestart(this,viewer)

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_geomech_force_type) :: this
  PetscViewer :: viewer
  
  call printErrMsg(this%option,'add code for restarting Geomech in PM approach')
  
end subroutine PMGeomechForceRestart

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/31/13
! ************************************************************************** !
subroutine PMGeomechForceDestroy(this)

  use Geomechanics_Realization_class, only : GeomechRealizDestroy

  implicit none
  
  class(pm_geomech_force_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_GEOMECH_FORCE_DEBUG
  call printMsg(this%option,'PMGeomechForce%Destroy()')
#endif

  !call GeomechRealizDestroy(this%geomech_realization)

  call this%comm1%Destroy()
  
end subroutine PMGeomechForceDestroy

end module PM_Geomechanics_Force_class

#endif
