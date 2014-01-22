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

#ifdef ABSTRACT
  type, abstract, public :: pm_base_type
#else
  type, public :: pm_base_type
#endif
    character(len=MAXWORDLENGTH) :: name
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    Vec :: solution_vec
    Vec :: residual_vec
    class(realization_base_type), pointer :: realization_base
    class(pm_base_type), pointer :: next
  contains
#ifdef ABSTRACT  
    procedure(PMBaseInit), public, deferred :: Init
    procedure(PMBaseThisOnly), public, deferred :: InitializeRun
    procedure(PMBaseThisOnly), public, deferred :: FinalizeRun
    procedure(PMBaseResidual), public, deferred :: Residual
    procedure(PMBaseJacobian), public, deferred :: Jacobian
    procedure(PMBaseUpdateTimestep), public, deferred :: UpdateTimestep
    procedure(PMBaseThisOnly), public, deferred :: InitializeTimestep
    procedure(PMBaseThisOnly), public, deferred :: PreSolve
    procedure(PMBaseThisOnly), public, deferred :: PostSolve
    procedure(PMBaseThisOnly), public, deferred :: FinalizeTimestep
    procedure(PMBaseFunctionThisOnly), public, deferred :: AcceptSolution
    procedure(PMBaseCheckUpdatePre), public, deferred :: CheckUpdatePre
    procedure(PMBaseCheckUpdatePost), public, deferred :: CheckUpdatePost
    procedure(PMBaseThisOnly), public, deferred :: TimeCut
    procedure(PMBaseThisOnly), public, deferred :: UpdateSolution
    procedure(PMBaseThisOnly), public, deferred :: MaxChange
    procedure(PMBaseComputeMassBalance), public, deferred :: ComputeMassBalance
    procedure(PMBaseThisOnly), public, deferred :: Destroy
    procedure(PMBaseRHSFunction), public, deferred :: RHSFunction
    procedure(PMBaseCheckpoint), public, deferred :: Checkpoint
    procedure(PMBaseCheckpoint), public, deferred :: Restart
#else
    procedure, public :: Init => PMBaseInit
    procedure, public :: InitializeRun => PMBaseThisOnly
    procedure, public :: FinalizeRun => PMBaseThisOnly
    procedure, public :: Residual => PMBaseResidual
    procedure, public :: Jacobian => PMBaseJacobian
    procedure, public :: UpdateTimestep => PMBaseUpdateTimestep
    procedure, public :: InitializeTimestep => PMBaseThisOnly
    procedure, public :: PreSolve => PMBaseThisOnly
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
    procedure, public :: Checkpoint => PMBaseCheckpoint
    procedure, public :: Restart => PMBaseCheckpoint
#endif
  end type pm_base_type
  
  type, public :: pm_base_header_type
    integer*8 :: ndof
  end type pm_base_header_type
    
#ifdef ABSTRACT  
  abstract interface
    subroutine PMBaseInit(this)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
    end subroutine PMBaseInit

    subroutine PMBaseResidual(this,snes,xx,r,ierr)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      SNES :: snes
      Vec :: xx
      Vec :: r
      PetscErrorCode :: ierr
    end subroutine PMBaseResidual

    subroutine PMBaseJacobian(this,snes,xx,A,B,flag,ierr)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      SNES :: snes
      Vec :: xx
      Mat :: A, B
      MatStructure flag
      PetscErrorCode :: ierr
    end subroutine PMBaseJacobian

    subroutine PMBaseUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      PetscReal :: dt
      PetscReal :: dt_max
      PetscInt :: iacceleration
      PetscInt :: num_newton_iterations
      PetscReal :: tfac(:)
    end subroutine PMBaseUpdateTimestep
    
    subroutine PMBaseCheckUpdatePre(this,line_search,P,dP,changed,ierr)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      SNESLineSearch :: line_search
      Vec :: P
      Vec :: dP
      PetscBool :: changed
      PetscErrorCode :: ierr
    end subroutine PMBaseCheckUpdatePre
    
    subroutine PMBaseCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                      P1_changed,ierr)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      SNESLineSearch :: line_search
      Vec :: P0
      Vec :: dP
      Vec :: P1
      PetscBool :: dP_changed
      PetscBool :: P1_changed
      PetscErrorCode :: ierr
    end subroutine PMBaseCheckUpdatePost
  
    subroutine PMBasePostSolve(this)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      PetscBool :: solution_accepted
    end subroutine PMBasePostSolve
    
    subroutine PMBaseThisOnly(this)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
    end subroutine PMBaseThisOnly
    
    subroutine PMBaseThisTime(this,time)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      PetscReal :: time
    end subroutine PMBaseThisTime
    
    function PMBaseFunctionThisOnly(this)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      PetscBool ::  PMBaseFunctionThisOnly
    end function PMBaseFunctionThisOnly
    
    subroutine PMBaseComputeMassBalance(this,mass_balance_array)
      import pm_base_type
      implicit none
      class(pm_base_type) :: this
      PetscReal :: mass_balance_array(:)
    end subroutine PMBaseComputeMassBalance

    subroutine PMBaseCheckpoint(this,viewer)
      import pm_base_type
      implicit none
#include "finclude/petscviewer.h"      
      class(pm_base_type) :: this
      PetscViewer ::  viewer
    end subroutine PMBaseFunctionThisOnly
    
  end interface
#endif  
  
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

#if 0

! ************************************************************************** !

recursive subroutine RunToTime(this,time)
  ! 
  ! PMBaseRunTo: Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pm_base_type) :: this  
  PetscReal :: time
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%RunToTime(time)
  endif
  
end subroutine PMBaseRunTo
#endif

#ifndef ABSTRACT

! ************************************************************************** !

subroutine PMBaseInit(this)
  implicit none
  class(pm_base_type) :: this
  print *, 'Must extend PMBaseInit.'
  stop
end subroutine PMBaseInit

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

subroutine PMBaseJacobian(this,snes,xx,A,B,flag,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseJacobian.'
  stop
end subroutine PMBaseJacobian

! ************************************************************************** !

subroutine PMBaseUpdateTimestep(this,dt,dt_max,iacceleration, &
                                num_newton_iterations,tfac)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
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

subroutine PMBasePostSolve(this)
  implicit none
  class(pm_base_type) :: this
  PetscBool :: solution_accepted
  print *, 'Must extend PMBasePostSolve.'
  stop
end subroutine PMBasePostSolve

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

subroutine PMBaseCheckpoint(this,viewer)
  implicit none
#include "finclude/petscviewer.h"      
  class(pm_base_type) :: this
  PetscViewer :: viewer
  print *, 'Must extend PMBaseCheckpoint/Restart.'
  stop
end subroutine PMBaseCheckpoint

#endif

end module PM_Base_class
