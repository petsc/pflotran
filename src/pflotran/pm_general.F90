module PM_General_class

  use PM_Base_class
!geh: using General_module here fails with gfortran (internal compiler error)
!  use General_module
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

  type, public, extends(pm_base_type) :: pm_general_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    PetscReal :: dPmax
    PetscReal :: dTmax
    PetscReal :: dXmax
    PetscReal :: dSmax
    PetscReal :: dPmax_allowable
    PetscReal :: dTmax_allowable
    PetscReal :: dXmax_allowable
    PetscReal :: dSmax_allowable
  contains
    procedure, public :: Init => PMGeneralInit
    procedure, public :: PMGeneralSetRealization
    procedure, public :: InitializeRun => PMGeneralInitializeRun
    procedure, public :: FinalizeRun => PMGeneralFinalizeRun
    procedure, public :: InitializeTimestep => PMGeneralInitializeTimestep
    procedure, public :: FinalizeTimestep => PMGeneralFinalizeTimeStep
    procedure, public :: Residual => PMGeneralResidual
    procedure, public :: Jacobian => PMGeneralJacobian
    procedure, public :: UpdateTimestep => PMGeneralUpdateTimestep
    procedure, public :: PreSolve => PMGeneralPreSolve
    procedure, public :: PostSolve => PMGeneralPostSolve
    procedure, public :: AcceptSolution => PMGeneralAcceptSolution
    procedure, public :: CheckUpdatePre => PMGeneralCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMGeneralCheckUpdatePost
    procedure, public :: TimeCut => PMGeneralTimeCut
    procedure, public :: UpdateSolution => PMGeneralUpdateSolution
    procedure, public :: MaxChange => PMGeneralMaxChange
    procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    procedure, public :: Checkpoint => PMGeneralCheckpoint    
    procedure, public :: Restart => PMGeneralRestart  
    procedure, public :: Destroy => PMGeneralDestroy
  end type pm_general_type
  
  public :: PMGeneralCreate
  
contains

! ************************************************************************** !

function PMGeneralCreate()
  ! 
  ! Creates General process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_general_type), pointer :: PMGeneralCreate

  class(pm_general_type), pointer :: general_pm
  
#ifdef PM_GENERAL_DEBUG  
  print *, 'PMGeneralCreate()'
#endif  

  allocate(general_pm)
  nullify(general_pm%option)
  nullify(general_pm%output_option)
  nullify(general_pm%realization)
  nullify(general_pm%comm1)
  general_pm%dPmax = 0.d0
  general_pm%dTmax = 0.d0
  general_pm%dXmax = 0.d0
  general_pm%dSmax = 0.d0
  general_pm%dPmax_allowable = 5.d5
  general_pm%dTmax_allowable = 5.d0
  general_pm%dXmax_allowable = 0.5d0
  general_pm%dSmax_allowable = 1.d0
  
  call PMBaseCreate(general_pm)
  general_pm%name = 'PMGeneral'

  PMGeneralCreate => general_pm
  
end function PMGeneralCreate

! ************************************************************************** !

subroutine PMGeneralInit(this)
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
  
  class(pm_general_type) :: this

#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%Init()')
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

end subroutine PMGeneralInit

! ************************************************************************** !

subroutine PMGeneralSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_general_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%SetRealization()')
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
  
end subroutine PMGeneralSetRealization

! ************************************************************************** !

subroutine PMGeneralInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralInitializeTimestep
  use Global_module
  use Variables_module, only : POROSITY, TORTUOSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
  
  implicit none
  
  class(pm_general_type) :: this

#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%InitializeTimestep()')
#endif

  this%option%flow_dt = this%option%dt

  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,POROSITY,0)
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY,0)
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 PERMEABILITY_X,0)
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 PERMEABILITY_Y,0)
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 PERMEABILITY_Z,0)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," GENERAL FLOW ",62("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call GeneralInitializeTimestep(this%realization)
  
end subroutine PMGeneralInitializeTimestep

! ************************************************************************** !

subroutine PMGeneralPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_general_type) :: this
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%PreSolve()')
#endif

end subroutine PMGeneralPreSolve

! ************************************************************************** !

subroutine PMGeneralPostSolve(this)
  ! 
  ! PMGeneralUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_general_type) :: this
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%PostSolve()')
#endif
  
end subroutine PMGeneralPostSolve

! ************************************************************************** !

subroutine PMGeneralFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralMaxChange
  use Global_module

  implicit none
  
  class(pm_general_type) :: this
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call GeneralMaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4)') this%option%dpmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') &
      this%option%dpmax
  endif  
  
end subroutine PMGeneralFinalizeTimestep

! ************************************************************************** !

function PMGeneralAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_general_type) :: this
  
  PetscBool :: PMGeneralAcceptSolution
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%AcceptSolution()')
#endif
  ! do nothing
  PMGeneralAcceptSolution = PETSC_TRUE
  
end function PMGeneralAcceptSolution

! ************************************************************************** !

subroutine PMGeneralUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: up, ut, ux, us, umin
  PetscReal :: dtt
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%UpdateTimestep()')
#endif
  
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    up = this%dPmax_allowable/(this%dPmax+0.1)
    ut = this%dTmax_allowable/(this%dTmax+1.d-5)
    ux = this%dXmax_allowable/(this%dXmax+1.d-5)
    us = this%dSmax_allowable/(this%dSmax+1.d-5)
    umin = min(up,ut,ux,us)
  endif
  dtt = fac * dt * (1.d0 + umin)
  dt = min(dtt,2.d0*dt,dt_max)
  
end subroutine PMGeneralUpdateTimestep

! ************************************************************************** !

recursive subroutine PMGeneralInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use General_module, only : GeneralUpdateSolution

  implicit none
  
  class(pm_general_type) :: this
  
 
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call GeneralUpdateSolution(this%realization)
#endif  
    
end subroutine PMGeneralInitializeRun

! ************************************************************************** !

recursive subroutine PMGeneralFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pm_general_type) :: this
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMGeneralFinalizeRun

! ************************************************************************** !

subroutine PMGeneralResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralResidual

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%Residual()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call GeneralResidualMFDLP(snes,xx,r,this%realization,ierr)
    case default
      call GeneralResidual(snes,xx,r,this%realization,ierr)
  end select

end subroutine PMGeneralResidual

! ************************************************************************** !

subroutine PMGeneralJacobian(this,snes,xx,A,B,flag,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralJacobian

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%Jacobian()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call GeneralJacobianMFDLP(snes,xx,A,B,flag,this%realization,ierr)
    case default
      call GeneralJacobian(snes,xx,A,B,flag,this%realization,ierr)
  end select

end subroutine PMGeneralJacobian

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralCheckUpdatePre

  implicit none
  
  class(pm_general_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY  
  call GeneralCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMGeneralCheckUpdatePre

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralCheckUpdatePost

  implicit none
  
  class(pm_general_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call GeneralCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMGeneralCheckUpdatePost

! ************************************************************************** !

subroutine PMGeneralTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralTimeCut

  implicit none
  
  class(pm_general_type) :: this
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call GeneralTimeCut(this%realization)

end subroutine PMGeneralTimeCut

! ************************************************************************** !

subroutine PMGeneralUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralUpdateSolution
  use Condition_module

  implicit none
  
  class(pm_general_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%UpdateSolution()')
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
  call GeneralUpdateSolution(this%realization)

end subroutine PMGeneralUpdateSolution     

! ************************************************************************** !

subroutine PMGeneralMaxChange(this)
  ! 
  ! Not needed given GeneralMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  use General_Aux_module

  implicit none
  
  class(pm_general_type) :: this
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: dP_liquid_max
  PetscReal :: dP_gas_max
  PetscReal :: dP_air_max
  PetscReal :: dX_air_max
  PetscReal :: dT_max
  PetscReal :: dS_gas_max
  PetscReal :: send_buffer(6), recv_buffer(6)
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscErrorCode :: ierr
  
  option => this%realization%option
  field => this%realization%field
  grid => this%realization%patch%grid

#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%MaxChange()')
#endif

  gen_auxvars => this%realization%patch%aux%General%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars
  
  dP_liquid_max = 0.d0
  dP_gas_max = 0.d0
  dP_air_max = 0.d0
  dX_air_max = 0.d0
  dT_max = 0.d0
  dS_gas_max = 0.d0
  
  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecAbs(field%flow_dxx,ierr)
  ! cannot use VecStrideMax on field%flow_dxx given the entries have no 
  ! meaning if state changes
  call VecGetArrayReadF90(field%flow_dxx,vec_ptr,ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ! flow_dxx only represents change if same primary dependent variable
    if (gen_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS) == &
        global_auxvars(ghosted_id)%istate) then
      offset = (local_id-1)*option%nflowdof
      select case(global_auxvars(ghosted_id)%istate)
        case(LIQUID_STATE)
          dP_liquid_max = max(dP_liquid_max, &
                              vec_ptr(offset+GENERAL_LIQUID_PRESSURE_DOF))
          dX_air_max = max(dX_air_max, &
                           vec_ptr(offset+ &
                                   GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF))
          dT_max = max(dT_max, &
                   vec_ptr(offset+GENERAL_LIQUID_STATE_TEMPERATURE_DOF))
        case(GAS_STATE)
          dP_gas_max = max(dP_gas_max, &
                       vec_ptr(offset+GENERAL_GAS_PRESSURE_DOF))
          dP_air_max = max(dP_air_max, &
                           vec_ptr(offset+GENERAL_AIR_PRESSURE_DOF))
          dT_max = max(dT_max, &
                       vec_ptr(offset+GENERAL_GAS_STATE_TEMPERATURE_DOF))
        case(TWO_PHASE_STATE)
          dP_gas_max = max(dP_gas_max, &
                           vec_ptr(offset+GENERAL_GAS_PRESSURE_DOF))
          dP_air_max = max(dP_air_max, &
                           vec_ptr(offset+GENERAL_AIR_PRESSURE_DOF))
          dS_gas_max = max(dS_gas_max, &
                           vec_ptr(offset+GENERAL_GAS_SATURATION_DOF))
      end select
    endif
  enddo
  call VecRestoreArrayReadF90(field%flow_dxx,vec_ptr,ierr)
  ! load the maximums into buffer for global reduction
  send_buffer(1) = dP_liquid_max
  send_buffer(2) = dP_gas_max
  send_buffer(3) = dP_air_max
  send_buffer(4) = dX_air_max
  send_buffer(5) = dT_max
  send_buffer(6) = dS_gas_max
  call MPI_Allreduce(send_buffer,recv_buffer,SIX_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  dP_liquid_max = send_buffer(1)
  dP_gas_max = send_buffer(2)
  dP_air_max = send_buffer(3)
  dX_air_max = send_buffer(4)
  dT_max = send_buffer(5)
  dS_gas_max = send_buffer(6)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4," dxa= ",1pe12.4," dt= ",1pe12.4," dsg= ",1pe12.4)') &
      dP_liquid_max,dP_gas_max,dP_air_max,dX_air_max,dT_max,dS_gas_max
  endif
  if (OptionPrintToScreen(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4," dxa= ",1pe12.4," dt= ",1pe12.4," dsg= ",1pe12.4)') &
      dP_liquid_max,dP_gas_max,dP_air_max,dX_air_max,dT_max,dS_gas_max
  endif
  this%dPmax = max(dP_liquid_max,dP_gas_max,dP_air_max)
  this%dTmax = dT_max
  this%dXmax = dX_air_max
  this%dSmax = dS_gas_max
  
end subroutine PMGeneralMaxChange

! ************************************************************************** !

subroutine PMGeneralComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralComputeMassBalance

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY  
  call GeneralComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMGeneralComputeMassBalance

! ************************************************************************** !

subroutine PMGeneralCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMGeneralCheckpoint

! ************************************************************************** !

subroutine PMGeneralRestart(this,viewer)
  ! 
  ! Restarts data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/30/13
  ! 

  use Checkpoint_module
  use General_module, only : GeneralUpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModel(viewer,this%realization)
  call GeneralUpdateAuxVars(this%realization,PETSC_FALSE)
  call this%UpdateSolution()
  
end subroutine PMGeneralRestart

! ************************************************************************** !

subroutine PMGeneralDestroy(this)
  ! 
  ! Destroys General process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralDestroy

  implicit none
  
  class(pm_general_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneralDestroy()')
#endif

#ifndef SIMPLIFY 
  call GeneralDestroy(this%realization)
#endif
  call this%comm1%Destroy()
  
end subroutine PMGeneralDestroy
  
end module PM_General_class
