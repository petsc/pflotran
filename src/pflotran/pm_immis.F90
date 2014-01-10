module PM_Immis_class

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

  type, public, extends(pm_base_type) :: pm_immis_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMImmisInit
    procedure, public :: PMImmisSetRealization
    procedure, public :: InitializeRun => PMImmisInitializeRun
    procedure, public :: FinalizeRun => PMImmisFinalizeRun
    procedure, public :: InitializeTimestep => PMImmisInitializeTimestep
    procedure, public :: FinalizeTimestep => PMImmisFinalizeTimeStep
    procedure, public :: Residual => PMImmisResidual
    procedure, public :: Jacobian => PMImmisJacobian
    procedure, public :: UpdateTimestep => PMImmisUpdateTimestep
    procedure, public :: PreSolve => PMImmisPreSolve
    procedure, public :: PostSolve => PMImmisPostSolve
    procedure, public :: AcceptSolution => PMImmisAcceptSolution
#if 0
    procedure, public :: CheckUpdatePre => PMImmisCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMImmisCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMImmisTimeCut
    procedure, public :: UpdateSolution => PMImmisUpdateSolution
    procedure, public :: MaxChange => PMImmisMaxChange
    procedure, public :: ComputeMassBalance => PMImmisComputeMassBalance
    procedure, public :: Checkpoint => PMImmisCheckpoint    
    procedure, public :: Restart => PMImmisRestart  
    procedure, public :: Destroy => PMImmisDestroy
  end type pm_immis_type
  
  public :: PMImmisCreate
  
contains

! ************************************************************************** !

function PMImmisCreate()
  ! 
  ! Creates Immiscible process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type), pointer :: PMImmisCreate

  class(pm_immis_type), pointer :: immis_pm
  
#ifdef PM__DEBUG  
  print *, 'PMImmisCreate()'
#endif  

  allocate(immis_pm)
  nullify(immis_pm%option)
  nullify(immis_pm%output_option)
  nullify(immis_pm%realization)
  nullify(immis_pm%comm1)

  call PMBaseCreate(immis_pm)
  immis_pm%name = 'PMImmis'

  PMImmisCreate => immis_pm
  
end function PMImmisCreate

! ************************************************************************** !

subroutine PMImmisInit(this)
  ! 
  ! Initializes variables associated with Richard
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

#ifndef SIMPLIFY  
  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
#endif

  implicit none
  
  class(pm_immis_type) :: this

#ifdef PM_IMMIS_DEBUG
  call printMsg(this%option,'PMImmis%Init()')
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

end subroutine PMImmisInit

! ************************************************************************** !

subroutine PMImmisSetRealization(this,realization)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_immis_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%SetRealization()')
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
  
end subroutine PMImmisSetRealization

! ************************************************************************** !

subroutine PMImmisInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisInitializeTimestep
  use Global_module
  
  implicit none
  
  class(pm_immis_type) :: this

#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%InitializeTimestep()')
#endif

  this%option%flow_dt = this%option%dt

#ifndef SIMPLIFY  
  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%porosity_loc, &
                              this%realization%field%porosity_loc)
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," IMMISCIBLE FLOW ",62("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call ImmisInitializeTimestep(this%realization)
  
end subroutine PMImmisInitializeTimestep

! ************************************************************************** !

subroutine PMImmisPreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Global_module

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%PreSolve()')
#endif

end subroutine PMImmisPreSolve

! ************************************************************************** !

subroutine PMImmisPostSolve(this)
  ! 
  ! PMImmisUpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Global_module

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%PostSolve()')
#endif
  
end subroutine PMImmisPostSolve

! ************************************************************************** !

subroutine PMImmisFinalizeTimestep(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisMaxChange
  use Global_module

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call ImmisMaxChange(this%realization)
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
  
end subroutine PMImmisFinalizeTimestep

! ************************************************************************** !

function PMImmisAcceptSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
  PetscBool :: PMImmisAcceptSolution
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%AcceptSolution()')
#endif
  ! do nothing
  PMImmisAcceptSolution = PETSC_TRUE
  
end function PMImmisAcceptSolution

! ************************************************************************** !

subroutine PMImmisUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
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
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%UpdateTimestep()')
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
  
end subroutine PMImmisUpdateTimestep

! ************************************************************************** !

recursive subroutine PMImmisInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisUpdateSolution

  implicit none
  
  class(pm_immis_type) :: this
  
 
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call ImmisUpdateSolution(this%realization)
#endif  
    
end subroutine PMImmisInitializeRun

! ************************************************************************** !

recursive subroutine PMImmisFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMImmisFinalizeRun

! ************************************************************************** !

subroutine PMImmisResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisResidual

  implicit none
  
  class(pm_immis_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%Residual()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call ImmisResidualMFDLP(snes,xx,r,this%realization,ierr)
    case default
      call ImmisResidual(snes,xx,r,this%realization,ierr)
  end select

end subroutine PMImmisResidual

! ************************************************************************** !

subroutine PMImmisJacobian(this,snes,xx,A,B,flag,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisJacobian

  implicit none
  
  class(pm_immis_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%Jacobian()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call ImmisJacobianMFDLP(snes,xx,A,B,flag,this%realization,ierr)
    case default
      call ImmisJacobian(snes,xx,A,B,flag,this%realization,ierr)
  end select

end subroutine PMImmisJacobian
    
#if 0

! ************************************************************************** !

subroutine PMImmisCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisCheckUpdatePre

  implicit none
  
  class(pm_immis_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY  
  call ImmisCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMImmisCheckUpdatePre

! ************************************************************************** !

subroutine PMImmisCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisCheckUpdatePost

  implicit none
  
  class(pm_immis_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call ImmisCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMImmisCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMImmisTimeCut(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisTimeCut

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call ImmisTimeCut(this%realization)

end subroutine PMImmisTimeCut

! ************************************************************************** !

subroutine PMImmisUpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisUpdateSolution
  use Condition_module

  implicit none
  
  class(pm_immis_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%UpdateSolution()')
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
  call ImmisUpdateSolution(this%realization)

end subroutine PMImmisUpdateSolution     

! ************************************************************************** !

subroutine PMImmisMaxChange(this)
  ! 
  ! Not needed given PMImmisMaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisMaxChange

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%MaxChange()')
#endif

  call ImmisMaxChange(this%realization)

end subroutine PMImmisMaxChange

! ************************************************************************** !

subroutine PMImmisComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisComputeMassBalance

  implicit none
  
  class(pm_immis_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY
  !geh: currently does not include "trapped" mass
  !call ImmisComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMImmisComputeMassBalance

! ************************************************************************** !

subroutine PMImmisCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Immiscible PM
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_immis_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMImmisCheckpoint

! ************************************************************************** !

subroutine PMImmisRestart(this,viewer)
  ! 
  ! Restarts data associated with Immiscible PM
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Checkpoint_module
  use Immis_module, only : ImmisUpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_immis_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModel(viewer,this%realization)
  call ImmisUpdateAuxVars(this%realization)
  call this%UpdateSolution()
  
end subroutine PMImmisRestart

! ************************************************************************** !

subroutine PMImmisDestroy(this)
  ! 
  ! Destroys Immiscible process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisDestroy

  implicit none
  
  class(pm_immis_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmisDestroy()')
#endif

#ifndef SIMPLIFY 
  call ImmisDestroy(this%realization)
#endif
  call this%comm1%Destroy()
  
end subroutine PMImmisDestroy
  
end module PM_Immis_class
