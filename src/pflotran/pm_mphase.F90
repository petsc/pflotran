module PM_Mphase_class

  use PM_Base_class
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

  type, public, extends(pm_base_type) :: pm_mphase_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMMphaseInit
    procedure, public :: PMMphaseSetRealization
    procedure, public :: InitializeRun => PMMphaseInitializeRun
    procedure, public :: FinalizeRun => PMMphaseFinalizeRun
    procedure, public :: InitializeTimestep => PMMphaseInitializeTimestep
    procedure, public :: FinalizeTimestep => PMMphaseFinalizeTimeStep
    procedure, public :: Residual => PMMphaseResidual
    procedure, public :: Jacobian => PMMphaseJacobian
    procedure, public :: UpdateTimestep => PMMphaseUpdateTimestep
    procedure, public :: PreSolve => PMMphasePreSolve
    procedure, public :: PostSolve => PMMphasePostSolve
    procedure, public :: AcceptSolution => PMMphaseAcceptSolution
#if 0
    procedure, public :: CheckUpdatePre => PMMphaseCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMMphaseCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMMphaseTimeCut
    procedure, public :: UpdateSolution => PMMphaseUpdateSolution
    procedure, public :: MaxChange => PMMphaseMaxChange
    procedure, public :: ComputeMassBalance => PMMphaseComputeMassBalance
    procedure, public :: Checkpoint => PMMphaseCheckpoint    
    procedure, public :: Restart => PMMphaseRestart  
    procedure, public :: Destroy => PMMphaseDestroy
  end type pm_mphase_type
  
  public :: PMMphaseCreate
  
contains

! ************************************************************************** !

function PMMphaseCreate()
  ! 
  ! Creates Mphase process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_mphase_type), pointer :: PMMphaseCreate

  class(pm_mphase_type), pointer :: mphase_pm
  
#ifdef PM_MPHASE_DEBUG  
  print *, 'PMMphaseCreate()'
#endif  

  allocate(mphase_pm)
  nullify(mphase_pm%option)
  nullify(mphase_pm%output_option)
  nullify(mphase_pm%realization)
  nullify(mphase_pm%comm1)

  call PMBaseCreate(mphase_pm)
  mphase_pm%name = 'PMMphase'

  PMMphaseCreate => mphase_pm
  
end function PMMphaseCreate

! ************************************************************************** !

subroutine PMMphaseInit(this)
  ! 
  ! Initializes variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

#ifndef SIMPLIFY  
  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
#endif

  implicit none
  
  class(pm_mphase_type) :: this

#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%Init()')
#endif
  
#ifndef SIMPLIFY  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      this%comm1 => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select
  call this%comm1%SetDM(this%realization%discretization%dm_1dof)
#endif

end subroutine PMMphaseInit

! ************************************************************************** !

subroutine PMMphaseSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_mphase_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%SetRealization()')
#endif
  
  this%realization => realization
  this%realization_base => realization

  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then 
    this%solution_vec = realization%field%flow_xx_faces
    this%residual_vec = realization%field%flow_r_faces
  else
    this%solution_vec = realization%field%flow_xx
    this%residual_vec = realization%field%flow_r
  endif
  
end subroutine PMMphaseSetRealization

! ************************************************************************** !

subroutine PMMphaseInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseInitializeTimestep
  use Global_module
  
  implicit none
  
  class(pm_mphase_type) :: this

#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%InitializeTimestep()')
#endif

  this%option%flow_dt = this%option%dt

#ifndef SIMPLIFY  
  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%porosity_loc, &
                              this%realization%field%porosity_loc)
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," MPHASE FLOW ",62("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call MphaseInitializeTimestep(this%realization)
  
end subroutine PMMphaseInitializeTimestep

! ************************************************************************** !

subroutine PMMphasePreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_mphase_type) :: this
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%PreSolve()')
#endif

end subroutine PMMphasePreSolve

! ************************************************************************** !

subroutine PMMphasePostSolve(this)
  ! 
  ! PMMphaseUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_mphase_type) :: this
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%PostSolve()')
#endif
  
end subroutine PMMphasePostSolve

! ************************************************************************** !

subroutine PMMphaseFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseMaxChange
  use Global_module

  implicit none
  
  class(pm_mphase_type) :: this
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call MphaseMaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          this%option%dpmax,this%option%dtmpmax,this%option%dcmax, &
          this%option%dsmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
      this%option%dpmax,this%option%dtmpmax,this%option%dcmax, &
      this%option%dsmax
  endif  
  
end subroutine PMMphaseFinalizeTimestep

! ************************************************************************** !

function PMMphaseAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_mphase_type) :: this
  
  PetscBool :: PMMphaseAcceptSolution
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%AcceptSolution()')
#endif
  ! do nothing
  PMMphaseAcceptSolution = PETSC_TRUE
  
end function PMMphaseAcceptSolution

! ************************************************************************** !

subroutine PMMphaseUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_mphase_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
  PetscReal :: uc
  PetscReal :: uus
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%UpdateTimestep()')
#endif
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%option%dpmxe/(this%option%dpmax+0.1)
      utmp = this%option%dtmpmxe/(this%option%dtmpmax+1.d-5)
      uc = this%option%dcmxe/(this%option%dcmax+1.d-6)
      uus= this%option%dsmxe/(this%option%dsmax+1.d-6)
      ut = min(up,utmp,uc,uus)
    endif
    dtt = fac * dt * (1.d0 + ut)
  else
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    dt_tfac = tfac(ifac) * dt

    fac = 0.5d0
    up = this%option%dpmxe/(this%option%dpmax+0.1)
    dt_p = fac * dt * (1.d0 + up)

    dtt = min(dt_tfac,dt_p)
  endif
  
  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
      
  dt = dtt
  
end subroutine PMMphaseUpdateTimestep

! ************************************************************************** !

recursive subroutine PMMphaseInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Mphase_module, only : MphaseUpdateSolution

  implicit none
  
  class(pm_mphase_type) :: this
  
 
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call MphaseUpdateSolution(this%realization)
#endif  
    
end subroutine PMMphaseInitializeRun

! ************************************************************************** !

recursive subroutine PMMphaseFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pm_mphase_type) :: this
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMMphaseFinalizeRun

! ************************************************************************** !

subroutine PMMphaseResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseResidual

  implicit none
  
  class(pm_mphase_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%Residual()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call MphaseResidualMFDLP(snes,xx,r,this%realization,ierr)
    case default
      call MphaseResidual(snes,xx,r,this%realization,ierr)
  end select

end subroutine PMMphaseResidual

! ************************************************************************** !

subroutine PMMphaseJacobian(this,snes,xx,A,B,flag,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseJacobian

  implicit none
  
  class(pm_mphase_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%Jacobian()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call MphaseJacobianMFDLP(snes,xx,A,B,flag,this%realization,ierr)
    case default
      call MphaseJacobian(snes,xx,A,B,flag,this%realization,ierr)
  end select

end subroutine PMMphaseJacobian
    
#if 0

! ************************************************************************** !

subroutine PMMphaseCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseCheckUpdatePre

  implicit none
  
  class(pm_mphase_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY  
  call MphaseCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMMphaseCheckUpdatePre

! ************************************************************************** !

subroutine PMMphaseCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseCheckUpdatePost

  implicit none
  
  class(pm_mphase_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call MphaseCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMMphaseCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMMphaseTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseTimeCut

  implicit none
  
  class(pm_mphase_type) :: this
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call MphaseTimeCut(this%realization)

end subroutine PMMphaseTimeCut

! ************************************************************************** !

subroutine PMMphaseUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseUpdateSolution
  use Condition_module

  implicit none
  
  class(pm_mphase_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%UpdateSolution()')
#endif

  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%realization%flow_conditions, &
                           this%realization%option, &
                           this%realization%option%time)
  ! right now, RealizUpdateAllCouplerAuxVars only updates flow
  call RealizUpdateAllCouplerAuxVars(this%realization,force_update_flag)
  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif  
  ! end from RealizationUpdate()
  call MphaseUpdateSolution(this%realization)

end subroutine PMMphaseUpdateSolution     

! ************************************************************************** !

subroutine PMMphaseMaxChange(this)
  ! 
  ! Not needed given MphaseMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseMaxChange

  implicit none
  
  class(pm_mphase_type) :: this
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%MaxChange()')
#endif

  call MphaseMaxChange(this%realization)

end subroutine PMMphaseMaxChange

! ************************************************************************** !

subroutine PMMphaseComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseComputeMassBalance

  implicit none
  
  class(pm_mphase_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphase%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY
  !geh: currently does not include "trapped" mass
  !call MphaseComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMMphaseComputeMassBalance

! ************************************************************************** !

subroutine PMMphaseCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Mphase PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_mphase_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMMphaseCheckpoint

! ************************************************************************** !

subroutine PMMphaseRestart(this,viewer)
  ! 
  ! Restarts data associated with Mphase PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/30/13
  ! 

  use Checkpoint_module
  use Mphase_module, only : MphaseUpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_mphase_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModel(viewer,this%realization)
  call MphaseUpdateAuxVars(this%realization)
  call this%UpdateSolution()
  
end subroutine PMMphaseRestart

! ************************************************************************** !

subroutine PMMphaseDestroy(this)
  ! 
  ! Destroys Mphase process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseDestroy

  implicit none
  
  class(pm_mphase_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_MPHASE_DEBUG  
  call printMsg(this%option,'PMMphaseDestroy()')
#endif

#ifndef SIMPLIFY 
  call MphaseDestroy(this%realization)
#endif
  call this%comm1%Destroy()
  
end subroutine PMMphaseDestroy
  
end module PM_Mphase_class
