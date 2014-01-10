module PM_Richards_class

  use PM_Base_class
!geh: using Richards_module here fails with gfortran (internal compiler error)
!  use Richards_module
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

  type, public, extends(pm_base_type) :: pm_richards_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMRichardsInit
    procedure, public :: PMRichardsSetRealization
    procedure, public :: InitializeRun => PMRichardsInitializeRun
    procedure, public :: FinalizeRun => PMRichardsFinalizeRun
    procedure, public :: InitializeTimestep => PMRichardsInitializeTimestep
    procedure, public :: FinalizeTimestep => PMRichardsFinalizeTimeStep
    procedure, public :: Residual => PMRichardsResidual
    procedure, public :: Jacobian => PMRichardsJacobian
    procedure, public :: UpdateTimestep => PMRichardsUpdateTimestep
    procedure, public :: PreSolve => PMRichardsPreSolve
    procedure, public :: PostSolve => PMRichardsPostSolve
    procedure, public :: AcceptSolution => PMRichardsAcceptSolution
    procedure, public :: CheckUpdatePre => PMRichardsCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRichardsCheckUpdatePost
    procedure, public :: TimeCut => PMRichardsTimeCut
    procedure, public :: UpdateSolution => PMRichardsUpdateSolution
    procedure, public :: MaxChange => PMRichardsMaxChange
    procedure, public :: ComputeMassBalance => PMRichardsComputeMassBalance
    procedure, public :: Checkpoint => PMRichardsCheckpoint    
    procedure, public :: Restart => PMRichardsRestart  
    procedure, public :: Destroy => PMRichardsDestroy
  end type pm_richards_type
  
  public :: PMRichardsCreate
  
contains

! ************************************************************************** !

function PMRichardsCreate()
  ! 
  ! Creates Richards process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_richards_type), pointer :: PMRichardsCreate

  class(pm_richards_type), pointer :: richards_pm
  
#ifdef PM_RICHARDS_DEBUG  
  print *, 'PMRichardsCreate()'
#endif  

  allocate(richards_pm)
  nullify(richards_pm%option)
  nullify(richards_pm%output_option)
  nullify(richards_pm%realization)
  nullify(richards_pm%comm1)

  call PMBaseCreate(richards_pm)
  richards_pm%name = 'PMRichards'

  PMRichardsCreate => richards_pm
  
end function PMRichardsCreate

! ************************************************************************** !

subroutine PMRichardsInit(this)
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
  
  class(pm_richards_type) :: this

#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%Init()')
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

end subroutine PMRichardsInit

! ************************************************************************** !

subroutine PMRichardsSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_richards_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%SetRealization()')
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
  
end subroutine PMRichardsSetRealization

! ************************************************************************** !

subroutine PMRichardsInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsInitializeTimestep
  use Global_module
  
  implicit none
  
  class(pm_richards_type) :: this

#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%InitializeTimestep()')
#endif

  this%option%flow_dt = this%option%dt

#ifndef SIMPLIFY  
  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%porosity_loc, &
                              this%realization%field%porosity_loc)
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," RICHARDS FLOW ",62("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call RichardsInitializeTimestep(this%realization)
  
end subroutine PMRichardsInitializeTimestep

! ************************************************************************** !

subroutine PMRichardsPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_richards_type) :: this
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%PreSolve()')
#endif

end subroutine PMRichardsPreSolve

! ************************************************************************** !

subroutine PMRichardsPostSolve(this)
  ! 
  ! PMRichardsUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_richards_type) :: this
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%PostSolve()')
#endif
  
end subroutine PMRichardsPostSolve

! ************************************************************************** !

subroutine PMRichardsFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsMaxChange
  use Global_module

  implicit none
  
  class(pm_richards_type) :: this
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call RichardsMaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4)') this%option%dpmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') &
      this%option%dpmax
  endif  
  
end subroutine PMRichardsFinalizeTimestep

! ************************************************************************** !

function PMRichardsAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_richards_type) :: this
  
  PetscBool :: PMRichardsAcceptSolution
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%AcceptSolution()')
#endif
  ! do nothing
  PMRichardsAcceptSolution = PETSC_TRUE
  
end function PMRichardsAcceptSolution

! ************************************************************************** !

subroutine PMRichardsUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_richards_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%UpdateTimestep()')
#endif
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%option%dpmxe/(this%option%dpmax+0.1)
      ut = up
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
  
end subroutine PMRichardsUpdateTimestep

! ************************************************************************** !

recursive subroutine PMRichardsInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Richards_module, only : RichardsUpdateSolution

  implicit none
  
  class(pm_richards_type) :: this
  
 
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call RichardsUpdateSolution(this%realization)
#endif  
    
end subroutine PMRichardsInitializeRun

! ************************************************************************** !

recursive subroutine PMRichardsFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pm_richards_type) :: this
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMRichardsFinalizeRun

! ************************************************************************** !

subroutine PMRichardsResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsResidual

  implicit none
  
  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%Residual()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call RichardsResidualMFDLP(snes,xx,r,this%realization,ierr)
    case default
      call RichardsResidual(snes,xx,r,this%realization,ierr)
  end select

end subroutine PMRichardsResidual

! ************************************************************************** !

subroutine PMRichardsJacobian(this,snes,xx,A,B,flag,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsJacobian

  implicit none
  
  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%Jacobian()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call RichardsJacobianMFDLP(snes,xx,A,B,flag,this%realization,ierr)
    case default
      call RichardsJacobian(snes,xx,A,B,flag,this%realization,ierr)
  end select

end subroutine PMRichardsJacobian

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsCheckUpdatePre

  implicit none
  
  class(pm_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY  
  call RichardsCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMRichardsCheckUpdatePre

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsCheckUpdatePost

  implicit none
  
  class(pm_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call RichardsCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMRichardsCheckUpdatePost

! ************************************************************************** !

subroutine PMRichardsTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsTimeCut

  implicit none
  
  class(pm_richards_type) :: this
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call RichardsTimeCut(this%realization)

end subroutine PMRichardsTimeCut

! ************************************************************************** !

subroutine PMRichardsUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsUpdateSolution, RichardsUpdateSurfacePress
  use Condition_module

  implicit none
  
  class(pm_richards_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%UpdateSolution()')
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
  call RichardsUpdateSolution(this%realization)
#ifdef SURFACE_FLOW
  if(this%option%nsurfflowdof>0) &
    call RichardsUpdateSurfacePress(this%realization)
#endif

end subroutine PMRichardsUpdateSolution     

! ************************************************************************** !

subroutine PMRichardsMaxChange(this)
  ! 
  ! Not needed given RichardsMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsMaxChange

  implicit none
  
  class(pm_richards_type) :: this
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%MaxChange()')
#endif

  call RichardsMaxChange(this%realization)

end subroutine PMRichardsMaxChange

! ************************************************************************** !

subroutine PMRichardsComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsComputeMassBalance

  implicit none
  
  class(pm_richards_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY  
  call RichardsComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMRichardsComputeMassBalance

! ************************************************************************** !

subroutine PMRichardsCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Richards PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_richards_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMRichardsCheckpoint

! ************************************************************************** !

subroutine PMRichardsRestart(this,viewer)
  ! 
  ! Restarts data associated with Richards PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/30/13
  ! 

  use Checkpoint_module
  use Richards_module, only : RichardsUpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_richards_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModel(viewer,this%realization)
  call RichardsUpdateAuxVars(this%realization)
  call this%UpdateSolution()
  
end subroutine PMRichardsRestart

! ************************************************************************** !

subroutine PMRichardsDestroy(this)
  ! 
  ! Destroys Richards process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsDestroy

  implicit none
  
  class(pm_richards_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichardsDestroy()')
#endif

#ifndef SIMPLIFY 
  call RichardsDestroy(this%realization)
#endif
  call this%comm1%Destroy()
  
end subroutine PMRichardsDestroy
  
end module PM_Richards_class
