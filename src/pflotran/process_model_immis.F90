module Process_Model_Immis_class

  use Process_Model_Base_class
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
!
! PMImmisCreate: Creates Immiscible process models shell
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
function PMImmisCreate()

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
!
! PMImmisInit: Initializes variables associated with Richard
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisInit(this)

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
!
! PMImmisSetRealization: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisSetRealization(this,realization)

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
! Should not need this as it is called in PreSolve.
! PMImmisInitializeTimestep: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisInitializeTimestep(this)

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
!
! PMImmisPreSolve: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisPreSolve(this)

  use Global_module

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%PreSolve()')
#endif

end subroutine PMImmisPreSolve

! ************************************************************************** !
!
! PMImmisUpdatePostSolve: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisPostSolve(this)

  use Global_module

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%PostSolve()')
#endif
  
end subroutine PMImmisPostSolve

! ************************************************************************** !
!
! PMImmisFinalizeTimestep: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisFinalizeTimestep(this)

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
!
! PMImmisAcceptSolution: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
function PMImmisAcceptSolution(this)

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
!
! PMImmisUpdateTimestep: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)

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
!
! PMImmisInitializeRun: Initializes the time stepping
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
recursive subroutine PMImmisInitializeRun(this)

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
!
! PMImmisFinalizeRun: Finalizes the time stepping
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
recursive subroutine PMImmisFinalizeRun(this)

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
!
! PMImmisResidual: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisResidual(this,snes,xx,r,ierr)

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
!
! PMImmisJacobian: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisJacobian(this,snes,xx,A,B,flag,ierr)

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
!
! PMImmisCheckUpdatePre: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisCheckUpdatePre(this,line_search,P,dP,changed,ierr)

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
!
! PMImmisCheckUpdatePost: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)

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
!
! PMImmisTimeCut: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisTimeCut(this)

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
!
! PMImmisUpdateSolution: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisUpdateSolution(this)

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
! Not needed given PMImmisMaxChange is called in PostSolve
! PMImmisMaxChange: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisMaxChange(this)

  use Immis_module, only : ImmisMaxChange

  implicit none
  
  class(pm_immis_type) :: this
  
#ifdef PM_IMMIS_DEBUG  
  call printMsg(this%option,'PMImmis%MaxChange()')
#endif

  call ImmisMaxChange(this%realization)

end subroutine PMImmisMaxChange
    
! ************************************************************************** !
!
! PMImmisComputeMassBalance: 
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisComputeMassBalance(this,mass_balance_array)

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
!
! PMImmisCheckpoint: Checkpoints data associated with Immiscible PM
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisCheckpoint(this,viewer)

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_immis_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMImmisCheckpoint


! ************************************************************************** !
!
! PMImmisRestart: Restarts data associated with Immiscible PM
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisRestart(this,viewer)

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
!
! PMImmisDestroy: Destroys Immiscible process model
! author: Gautam Bisht
! date: 11/27/13
!
! ************************************************************************** !
subroutine PMImmisDestroy(this)

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
  
end module Process_Model_Immis_class
