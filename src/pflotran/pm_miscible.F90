module PM_Miscible_class

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

  type, public, extends(pm_base_type) :: pm_miscible_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMMiscibleInit
    procedure, public :: PMMiscibleSetRealization
    procedure, public :: InitializeRun => PMMiscibleInitializeRun
    procedure, public :: FinalizeRun => PMMiscibleFinalizeRun
    procedure, public :: InitializeTimestep => PMMiscibleInitializeTimestep
    procedure, public :: FinalizeTimestep => PMMiscibleFinalizeTimeStep
    procedure, public :: Residual => PMMiscibleResidual
    procedure, public :: Jacobian => PMMiscibleJacobian
    procedure, public :: UpdateTimestep => PMMiscibleUpdateTimestep
    procedure, public :: PreSolve => PMMisciblePreSolve
    procedure, public :: PostSolve => PMMisciblePostSolve
    procedure, public :: AcceptSolution => PMMiscibleAcceptSolution
#if 0
    procedure, public :: CheckUpdatePre => PMMiscibleCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMMiscibleCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMMiscibleTimeCut
    procedure, public :: UpdateSolution => PMMiscibleUpdateSolution
    procedure, public :: MaxChange => PMMiscibleMaxChange
    procedure, public :: ComputeMassBalance => PMMiscibleComputeMassBalance
    procedure, public :: Checkpoint => PMMiscibleCheckpoint    
    procedure, public :: Restart => PMMiscibleRestart  
    procedure, public :: Destroy => PMMiscibleDestroy
  end type pm_miscible_type
  
  public :: PMMiscibleCreate
  
contains

! ************************************************************************** !

function PMMiscibleCreate()
  ! 
  ! Creates Miscible process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type), pointer :: PMMiscibleCreate

  class(pm_miscible_type), pointer :: miscible_pm
  
#ifdef PM__DEBUG  
  print *, 'PMMiscibleCreate()'
#endif  

  allocate(miscible_pm)
  nullify(miscible_pm%option)
  nullify(miscible_pm%output_option)
  nullify(miscible_pm%realization)
  nullify(miscible_pm%comm1)

  call PMBaseCreate(miscible_pm)
  miscible_pm%name = 'PMMiscible'

  PMMiscibleCreate => miscible_pm
  
end function PMMiscibleCreate

! ************************************************************************** !

subroutine PMMiscibleInit(this)
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
  
  class(pm_miscible_type) :: this

#ifdef PM_MISCIBLE_DEBUG
  call printMsg(this%option,'PMMiscible%Init()')
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

end subroutine PMMiscibleInit

! ************************************************************************** !

subroutine PMMiscibleSetRealization(this,realization)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_miscible_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_miscible_DEBUG  
  call printMsg(this%option,'PMMiscible%SetRealization()')
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
  
end subroutine PMMiscibleSetRealization

! ************************************************************************** !

subroutine PMMiscibleInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleInitializeTimestep
  use Global_module
  
  implicit none
  
  class(pm_miscible_type) :: this

#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%InitializeTimestep()')
#endif

  this%option%flow_dt = this%option%dt

#ifndef SIMPLIFY  
  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%porosity_loc, &
                              this%realization%field%porosity_loc)
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," MISCIBLE FLOW ",62("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call MiscibleInitializeTimestep(this%realization)
  
end subroutine PMMiscibleInitializeTimestep

! ************************************************************************** !

subroutine PMMisciblePreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Global_module

  implicit none
  
  class(pm_miscible_type) :: this
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%PreSolve()')
#endif

end subroutine PMMisciblePreSolve

! ************************************************************************** !

subroutine PMMisciblePostSolve(this)
  ! 
  ! PMMiscibleUpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Global_module

  implicit none
  
  class(pm_miscible_type) :: this
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%PostSolve()')
#endif
  
end subroutine PMMisciblePostSolve

! ************************************************************************** !

subroutine PMMiscibleFinalizeTimestep(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleMaxChange
  use Global_module

  implicit none
  
  class(pm_miscible_type) :: this
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call MiscibleMaxChange(this%realization)
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
  
end subroutine PMMiscibleFinalizeTimestep

! ************************************************************************** !

function PMMiscibleAcceptSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
  
  PetscBool :: PMMiscibleAcceptSolution
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%AcceptSolution()')
#endif
  ! do nothing
  PMMiscibleAcceptSolution = PETSC_TRUE
  
end function PMMiscibleAcceptSolution

! ************************************************************************** !

subroutine PMMiscibleUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
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
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%UpdateTimestep()')
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
  
end subroutine PMMiscibleUpdateTimestep

! ************************************************************************** !

recursive subroutine PMMiscibleInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleUpdateSolution

  implicit none
  
  class(pm_miscible_type) :: this
  
 
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call MiscibleUpdateSolution(this%realization)
#endif  
    
end subroutine PMMiscibleInitializeRun

! ************************************************************************** !

recursive subroutine PMMiscibleFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMMiscibleFinalizeRun

! ************************************************************************** !

subroutine PMMiscibleResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleResidual

  implicit none
  
  class(pm_miscible_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%Residual()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call MiscibleResidualMFDLP(snes,xx,r,this%realization,ierr)
    case default
      call MiscibleResidual(snes,xx,r,this%realization,ierr)
  end select

end subroutine PMMiscibleResidual

! ************************************************************************** !

subroutine PMMiscibleJacobian(this,snes,xx,A,B,flag,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleJacobian

  implicit none
  
  class(pm_miscible_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%Jacobian()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call MiscibleJacobianMFDLP(snes,xx,A,B,flag,this%realization,ierr)
    case default
      call MiscibleJacobian(snes,xx,A,B,flag,this%realization,ierr)
  end select

end subroutine PMMiscibleJacobian
    
#if 0

! ************************************************************************** !

subroutine PMMiscibleCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleCheckUpdatePre

  implicit none
  
  class(pm_miscible_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY  
  call MiscibleCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMMiscibleCheckUpdatePre

! ************************************************************************** !

subroutine PMMiscibleCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleCheckUpdatePost

  implicit none
  
  class(pm_miscible_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call MiscibleCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMMiscibleCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMMiscibleTimeCut(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleTimeCut

  implicit none
  
  class(pm_miscible_type) :: this
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call MiscibleTimeCut(this%realization)

end subroutine PMMiscibleTimeCut

! ************************************************************************** !

subroutine PMMiscibleUpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleUpdateSolution
  use Condition_module

  implicit none
  
  class(pm_miscible_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%UpdateSolution()')
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
  call MiscibleUpdateSolution(this%realization)

end subroutine PMMiscibleUpdateSolution     

! ************************************************************************** !

subroutine PMMiscibleMaxChange(this)
  ! 
  ! Not needed given PMMiscibleMaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleMaxChange

  implicit none
  
  class(pm_miscible_type) :: this
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%MaxChange()')
#endif

  call MiscibleMaxChange(this%realization)

end subroutine PMMiscibleMaxChange

! ************************************************************************** !

subroutine PMMiscibleComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleComputeMassBalance

  implicit none
  
  class(pm_miscible_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY
  !geh: currently does not include "trapped" mass
  !call MiscibleComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMMiscibleComputeMassBalance

! ************************************************************************** !

subroutine PMMiscibleCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Miscible PM
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_miscible_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMMiscibleCheckpoint

! ************************************************************************** !

subroutine PMMiscibleRestart(this,viewer)
  ! 
  ! Restarts data associated with Miscible PM
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Checkpoint_module
  use Miscible_module, only : MiscibleUpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_miscible_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModel(viewer,this%realization)
  call MiscibleUpdateAuxVars(this%realization)
  call this%UpdateSolution()
  
end subroutine PMMiscibleRestart

! ************************************************************************** !

subroutine PMMiscibleDestroy(this)
  ! 
  ! Destroys Miscible process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleDestroy

  implicit none
  
  class(pm_miscible_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscibleDestroy()')
#endif

#ifndef SIMPLIFY 
  call MiscibleDestroy(this%realization)
#endif
  call this%comm1%Destroy()
  
end subroutine PMMiscibleDestroy
  
end module PM_Miscible_class
