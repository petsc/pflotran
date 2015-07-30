module PM_Base_class

  use Option_module
  use Output_Aux_module
  use Realization_Base_class

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

  type, public :: pm_base_type
    character(len=MAXWORDLENGTH) :: name
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    Vec :: solution_vec
    Vec :: residual_vec
    class(realization_base_type), pointer :: realization_base
    class(pm_base_type), pointer :: next
  contains
    procedure, public :: Init => PMBaseInit
    procedure, public :: Read => PMBaseRead
    procedure, public :: SetupSolvers => PMBaseSetupSolvers
    procedure, public :: InitializeRun => PMBaseThisOnly
    procedure, public :: FinalizeRun => PMBaseThisOnly
    procedure, public :: Residual => PMBaseResidual
    procedure, public :: Jacobian => PMBaseJacobian
    procedure, public :: UpdateTimestep => PMBaseUpdateTimestep
    procedure, public :: InitializeTimestep => PMBaseThisOnly
    procedure, public :: PreSolve => PMBaseThisOnly
    procedure, public :: Solve => PMBaseThisTimeError
    procedure, public :: PostSolve => PMBaseThisOnly
    procedure, public :: FinalizeTimestep => PMBaseThisOnly
    procedure, public :: AcceptSolution => PMBaseFunctionThisOnly
    procedure, public :: CheckUpdatePre => PMBaseCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMBaseCheckUpdatePost
    procedure, public :: TimeCut => PMBaseThisOnly
    procedure, public :: UpdateSolution => PMBaseThisOnly
    procedure, public :: MaxChange => PMBaseThisOnly
    procedure, public :: ComputeMassBalance => PMBaseComputeMassBalance
    procedure, public :: Destroy => PMBaseThisOnly
    procedure, public :: RHSFunction => PMBaseRHSFunction
    procedure, public :: CheckpointBinary => PMBaseCheckpointBinary
    procedure, public :: CheckpointHDF5 => PMBaseCheckpointHDF5
    procedure, public :: RestartBinary => PMBaseCheckpointBinary
  end type pm_base_type
  
  type, public :: pm_base_header_type
    PetscInt :: ndof
  end type pm_base_header_type
    
  public :: PMBaseCreate
  
  public :: PMBaseResidual
  public :: PMBaseJacobian
  public :: PMBaseRHSFunction
  
contains

! ************************************************************************** !

subroutine PMBaseCreate(this)

  implicit none
  
  class(pm_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  this%name = ''
  nullify(this%option)
  nullify(this%output_option)
  nullify(this%realization_base)
  this%solution_vec = 0
  this%residual_vec = 0
  nullify(this%next)
  
end subroutine PMBaseCreate

! ************************************************************************** !

subroutine PMBaseRead(this,input)
  use Input_Aux_module
  implicit none
  class(pm_base_type) :: this
  type(input_type) :: input
  print *, 'Must extend PMBaseRead.'
  stop
end subroutine PMBaseRead

! ************************************************************************** !

subroutine PMBaseInit(this)
  implicit none
  class(pm_base_type) :: this
  print *, 'Must extend PMBaseInit.'
  stop
end subroutine PMBaseInit

! ************************************************************************** !

subroutine PMBaseSetupSolvers(this,solver)
  use Solver_module
  implicit none
  class(pm_base_type) :: this
  type(solver_type) :: solver
  print *, 'Must extend PMBaseSetupSolvers.'
  stop
end subroutine PMBaseSetupSolvers

! ************************************************************************** !

subroutine PMBaseResidual(this,snes,xx,r,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseResidual.'
  stop
end subroutine PMBaseResidual

! ************************************************************************** !

subroutine PMBaseJacobian(this,snes,xx,A,B,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseJacobian.'
  stop
end subroutine PMBaseJacobian

! ************************************************************************** !

subroutine PMBaseUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                num_newton_iterations,tfac)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  print *, 'Must extend PMBaseUpdateTimestep.'
  stop
end subroutine PMBaseUpdateTimestep

! ************************************************************************** !

subroutine PMBaseCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  implicit none
  class(pm_base_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseCheckUpdatePre.'
  stop
end subroutine PMBaseCheckUpdatePre

! ************************************************************************** !

subroutine PMBaseCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  implicit none
  class(pm_base_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseCheckUpdatePost.'
  stop
end subroutine PMBaseCheckUpdatePost

! ************************************************************************** !

subroutine PMBaseThisOnly(this)
  implicit none
  class(pm_base_type) :: this
  print *, 'Must extend PMBaseThisOnly.'
  stop
end subroutine PMBaseThisOnly

! ************************************************************************** !

subroutine PMBaseThisTime(this,time)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: time
  print *, 'Must extend PMBaseThisTime.'
  stop
end subroutine PMBaseThisTime

! ************************************************************************** !

subroutine PMBaseThisTimeError(this,time,ierr)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: time
  PetscInt :: ierr
  print *, 'Must extend PMBaseThisTimeError.'
  stop
end subroutine PMBaseThisTimeError

! ************************************************************************** !

function PMBaseFunctionThisOnly(this)
  implicit none
  class(pm_base_type) :: this
  PetscBool ::  PMBaseFunctionThisOnly
  PMBaseFunctionThisOnly = PETSC_TRUE
  print *, 'Must extend PMBaseFunctionThisOnly.'
  stop
end function PMBaseFunctionThisOnly

! ************************************************************************** !

subroutine PMBaseComputeMassBalance(this,mass_balance_array)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: mass_balance_array(:)
  print *, 'Must extend PMBaseComputeMassBalance.'
  stop
end subroutine PMBaseComputeMassBalance

! ************************************************************************** !

subroutine PMBaseRHSFunction(this,ts,time,xx,ff,ierr)
  implicit none
  class(pm_base_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseRHSFunction.'
  stop
end subroutine PMBaseRHSFunction

! ************************************************************************** !

subroutine PMBaseCheckpointBinary(this,viewer)
  implicit none
#include "finclude/petscviewer.h"      
  class(pm_base_type) :: this
  PetscViewer :: viewer
  print *, 'Must extend PMBaseCheckpointBinary/RestartBinary.'
  stop
end subroutine PMBaseCheckpointBinary

! ************************************************************************** !

subroutine PMBaseCheckpointHDF5(this, pm_grp_id)

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  class(pm_base_type) :: this
  integer :: pm_grp_id
  print *, 'PFLOTRAN must be compiled with HDF5 to ' // &
        'write HDF5 formatted checkpoint file. Darn.'
  stop
#else

  use hdf5
  implicit none

  class(pm_base_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif
  print *, 'Must extend PMBaseCheckpointHDF5/RestartHDF5.'
  stop
#endif

end subroutine PMBaseCheckpointHDF5

end module PM_Base_class
