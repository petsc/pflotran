module Process_Model_RT_class

  use Process_Model_Base_class
!geh: using Reactive_Transport_module here fails with gfortran (internal 
!     compiler error)
!  use Reactive_Transport_module
  use Realization_class
  use Communicator_Base_module  
  use Option_module
  
  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(pm_base_type) :: pm_rt_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    class(communicator_type), pointer :: commN
    ! local variables
    PetscBool :: steady_flow
    PetscReal :: tran_weight_t0
    PetscReal :: tran_weight_t1
  contains
    procedure, public :: Init => PMRTInit
    procedure, public :: PMRTSetRealization
    procedure, public :: InitializeRun => PMRTInitializeRun
    procedure, public :: FinalizeRun => PMRTFinalizeRun
    procedure, public :: InitializeTimestep => PMRTInitializeTimestep
    procedure, public :: FinalizeTimestep => PMRTFinalizeTimestep
    procedure, public :: Residual => PMRTResidual
    procedure, public :: Jacobian => PMRTJacobian
    procedure, public :: UpdateTimestep => PMRTUpdateTimestep
    procedure, public :: PreSolve => PMRTPreSolve
    procedure, public :: PostSolve => PMRTPostSolve
    procedure, public :: AcceptSolution => PMRTAcceptSolution
    procedure, public :: CheckUpdatePre => PMRTCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRTCheckUpdatePost
    procedure, public :: TimeCut => PMRTTimeCut
    procedure, public :: UpdateSolution => PMRTUpdateSolution
    procedure, public :: MaxChange => PMRTMaxChange
    procedure, public :: ComputeMassBalance => PMRTComputeMassBalance
    procedure, public :: SetTranWeights => SetTranWeights
    procedure, public :: Destroy => PMRTDestroy
  end type pm_rt_type
  
  public :: PMRTCreate

contains

! ************************************************************************** !
!
! PMRTCreate: Creates reactive transport process models shell
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function PMRTCreate()

  implicit none
  
  class(pm_rt_type), pointer :: PMRTCreate

  class(pm_rt_type), pointer :: rt_pm
  
#ifdef PM_RT_DEBUG  
  print *, 'PMRTCreate()'
#endif
  
  allocate(rt_pm)
  nullify(rt_pm%option)
  nullify(rt_pm%output_option)
  nullify(rt_pm%realization)
  nullify(rt_pm%comm1)
  nullify(rt_pm%commN)
  
  ! local variables
  rt_pm%steady_flow = PETSC_FALSE
  rt_pm%tran_weight_t0 = 0.d0
  rt_pm%tran_weight_t1 = 0.d0

  call PMBaseCreate(rt_pm)
  
  PMRTCreate => rt_pm
  
end function PMRTCreate

! ************************************************************************** !
!
! PMRTInit: Initializes variables associated with reactive transport
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTInit(this)

#ifndef SIMPLIFY
  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
#endif  
  
  implicit none
  
  class(pm_rt_type) :: this

#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%Init()')
#endif
  
#ifndef SIMPLIFY  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      this%comm1 => StructuredCommunicatorCreate()
      this%commN => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
      this%commN => UnstructuredCommunicatorCreate()
  end select
  call this%comm1%SetDM(this%realization%discretization%dm_1dof)
  call this%commN%SetDM(this%realization%discretization%dm_ntrandof)
#endif

end subroutine PMRTInit

! ************************************************************************** !
!
! PMRTSetRealization: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTSetRealization(this,realization)

  use Realization_class  

  implicit none
  
  class(pm_rt_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%SetRealization()')
#endif
  
  this%realization => realization
  this%realization_base => realization
  
  if (realization%reaction%use_log_formulation) then
    this%solution_vec = realization%field%tran_log_xx
  else
    this%solution_vec = realization%field%tran_xx
  endif
  this%residual_vec = realization%field%tran_r
  
end subroutine PMRTSetRealization

! ************************************************************************** !
!
! PMRTInitializeTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTInitializeTimestep(this)

  use Reactive_Transport_module, only : RTInitializeTimestep, &
                                        RTUpdateTransportCoefs
  use Global_module

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: time
 
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%InitializeTimestep()')
#endif
  
  this%option%tran_dt = this%option%dt

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," REACTIVE TRANSPORT ",57("="))')
  endif
  
#if 0  
  call DiscretizationLocalToLocal(discretization, &
                                  this%realization%field%porosity_loc, &
                                  this%realization%field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization, &
                                  this%realization%field%tortuosity_loc, &
                                  this%realization%field%tortuosity_loc,ONEDOF)
#endif  
  
  ! interpolate flow parameters/data
  ! this must remain here as these weighted values are used by both
  ! RTInitializeTimestep and RTTimeCut (which calls RTInitializeTimestep)
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    call this%SetTranWeights()
    ! set densities and saturations to t
    call GlobalUpdateDenAndSat(this%realization,this%tran_weight_t0)
  endif

  call RTInitializeTimestep(this%realization)

  !geh: this is a bug and should be moved to PreSolve()
#if 1
  ! set densities and saturations to t+dt
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    call GlobalUpdateDenAndSat(this%realization,this%tran_weight_t1)
  endif

  call RTUpdateTransportCoefs(this%realization)
#endif  

end subroutine PMRTInitializeTimestep

! ************************************************************************** !
!
! PMRTPreSolve: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTPreSolve(this)

  use Reactive_Transport_module, only : RTUpdateTransportCoefs, &
                                        RTUpdateAuxVars
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Global_module  

  implicit none
  
  class(pm_rt_type) :: this
  
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%UpdatePreSolve()')
#endif
  
#if 0
  ! set densities and saturations to t+dt
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    call GlobalUpdateDenAndSat(this%realization,this%tran_weight_t1)
  endif

  call RTUpdateTransportCoefs(this%realization)
#endif  
  
  if (this%realization%reaction%act_coef_update_frequency /= &
      ACT_COEF_FREQUENCY_OFF) then
      call RTUpdateAuxVars(this%realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
!       The below is set within RTUpdateAuxVarsPatch() when 
!         PETSC_TRUE,PETSC_TRUE,* are passed
!       patch%aux%RT%aux_vars_up_to_date = PETSC_TRUE 
  endif
  if (this%realization%reaction%use_log_formulation) then
    call VecCopy(this%realization%field%tran_xx, &
                 this%realization%field%tran_log_xx,ierr)
    call VecLog(this%realization%field%tran_log_xx,ierr)
  endif
  
end subroutine PMRTPreSolve

! ************************************************************************** !
!
! PMRTPostSolve: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTPostSolve(this)

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%PostSolve()')
#endif
  
end subroutine PMRTPostSolve

! ************************************************************************** !
!
! PMRTFinalizeTimestep: 
! author: Glenn Hammond
! date: 04/03/13
!
! ************************************************************************** !
subroutine PMRTFinalizeTimestep(this)

  use Reactive_Transport_module, only : RTMaxChange
  use Global_module

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: time  
  
#ifdef PM_RICHARDS_DEBUG  
  call printMsg(this%option,'PMRichards%FinalizeTimestep()')
#endif
  
  call RTMaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dcmx= ",1pe12.4," dc/dt= ",1pe12.4, &
            &" [mol/s]")') &
      this%option%dcmax,this%option%dcmax/this%option%tran_dt
  endif
  if (this%option%print_file_flag) then  
    write(this%option%fid_out,'("  --> max chng: dcmx= ",1pe12.4, &
                              &" dc/dt= ",1pe12.4," [mol/s]")') &
      this%option%dcmax,this%option%dcmax/this%option%tran_dt
  endif
  
end subroutine PMRTFinalizeTimestep

! ************************************************************************** !
!
! PMRichardsAcceptSolution: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function PMRTAcceptSolution(this)

  implicit none
  
  class(pm_rt_type) :: this
  
  PetscBool :: PMRTAcceptSolution
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%AcceptSolution()')
#endif
  ! do nothing
  PMRTAcceptSolution = PETSC_TRUE
  
end function PMRTAcceptSolution

! ************************************************************************** !
!
! PMRTUpdateTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTUpdateTimestep(this,dt,dt_max,iacceleration, &
                              num_newton_iterations,tfac)

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: dtt
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%UpdateTimestep()')  
#endif
  
  dtt = dt
  if (num_newton_iterations <= iacceleration) then
    if (num_newton_iterations <= size(tfac)) then
      dtt = tfac(num_newton_iterations) * dt
    else
      dtt = 0.5d0 * dt
    endif
  else
!       dtt = 2.d0 * dt
    dtt = 0.5d0 * dt
  endif

  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > dt_max) dtt = dt_max
  ! geh: see comment above under flow stepper
  dt = dtt

end subroutine PMRTUpdateTimestep

! ************************************************************************** !
!
! PMRTInitializeRun: Initializes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMRTInitializeRun(this)

  use Reactive_Transport_module, only : RTUpdateEquilibriumState, &
                                        RTJumpStartKineticSorption

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%InitializeRun()')
#endif
  
  ! restart
#if 0  
  if (transport_read .and. option%overwrite_restart_transport) then
    call CondControlAssignTranInitCond(realization)  
  endif
#endif  
  
  call RTUpdateEquilibriumState(this%realization)
  
#if 0
  if (this%option%jumpstart_kinetic_sorption .and. &
      this%option%time < 1.d-40) then
    ! only user jumpstart for a restarted simulation
    if (.not. this%option%restart_flag) then
      this%option%io_buffer = 'Only use JUMPSTART_KINETIC_SORPTION on a ' // &
        'restarted simulation.  ReactionEquilibrateConstraint() will ' // &
        'appropriately set sorbed initial concentrations for a normal ' // &
        '(non-restarted) simulation.'
      call printErrMsg(this%option)
    endif
    call RTJumpStartKineticSorption(this%realization)
  endif
  ! check on MAX_STEPS < 0 to quit after initialization.
#endif  
    
end subroutine PMRTInitializeRun

! ************************************************************************** !
!
! PMRTFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMRTFinalizeRun(this)

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%PMRTFinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMRTFinalizeRun

! ************************************************************************** !
!
! PMRTResidual: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTResidual(this,snes,xx,r,ierr)

  use Reactive_Transport_module, only : RTResidual

  implicit none
  
  class(pm_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%Residual()')  
#endif
  
  call RTResidual(snes,xx,r,this%realization,ierr)

end subroutine PMRTResidual

! ************************************************************************** !
!
! PMRTJacobian: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTJacobian(this,snes,xx,A,B,flag,ierr)

  use Reactive_Transport_module, only : RTJacobian

  implicit none
  
  class(pm_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%Jacobian()')  
#endif

  call RTJacobian(snes,xx,A,B,flag,this%realization,ierr)

end subroutine PMRTJacobian
    
! ************************************************************************** !
!
! PMRTCheckUpdatePre: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTCheckUpdatePre(this,line_search,P,dP,changed,ierr)

  use Reactive_Transport_module, only : RTCheckUpdate

  implicit none
  
  class(pm_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY 
  call RTCheckUpdate(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMRTCheckUpdatePre
    
! ************************************************************************** !
!
! PMRTCheckUpdatePost: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)

!  use Reactive_Transport_module, only : RTCheckUpdatePost

  implicit none
  
  class(pm_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%CheckUpdatePost()')
#endif
  
!  call RTCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
!                               P1_changed,this%realization,ierr)

end subroutine PMRTCheckUpdatePost
  
! ************************************************************************** !
!
! PMRTTimeCut: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTTimeCut(this)

  use Reactive_Transport_module, only : RTTimeCut

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%TimeCut()')
#endif
  
  this%option%tran_dt = this%option%dt
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    call this%SetTranWeights()
  endif
  call RTTimeCut(this%realization)

end subroutine PMRTTimeCut
    
! ************************************************************************** !
!
! PMRTUpdateSolution: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTUpdateSolution(this)

  use Reactive_Transport_module
  use Condition_module
  use Mass_Transfer_module

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%UpdateSolution()')
#endif
  
  ! begin from RealizationUpdate()
  call TranConditionUpdate(this%realization%transport_conditions, &
                           this%realization%option, &
                           this%realization%option%time)
  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif  
  ! end from RealizationUpdate()
  ! The update of status must be in this order!
  call RTUpdateEquilibriumState(this%realization)
  call RTUpdateKineticState(this%realization)
  if (this%realization%reaction%update_porosity .or. &
      this%realization%reaction%update_tortuosity .or. &
      this%realization%reaction%update_permeability .or. &
      this%realization%reaction%update_mineral_surface_area) then
    call RealizationUpdateProperties(this%realization)
  endif
  
  call MassTransferUpdate(this%realization%mass_transfer_list, &
                          this%realization%discretization, &
                          this%realization%patch%grid, &
                          this%realization%option)
  
  if (this%realization%option%compute_mass_balance_new) then
    call RTUpdateMassBalance(this%realization)
  endif  

end subroutine PMRTUpdateSolution     

! ************************************************************************** !
!
! PMRTMaxChange: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTMaxChange(this)

  use Reactive_Transport_module, only : RTMaxChange

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%MaxChange()')
#endif
  
  call RTMaxChange(this%realization)

end subroutine PMRTMaxChange
    
! ************************************************************************** !
!
! PMRTComputeMassBalance: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTComputeMassBalance(this,mass_balance_array)

  use Reactive_Transport_module, only : RTComputeMassBalance

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: mass_balance_array(:)

#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%MassBalance()')
#endif

#ifndef SIMPLIFY 
  call RTComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMRTComputeMassBalance

! ************************************************************************** !
!
! SetTranWeights: Sets the weights at t0 or t1 for transport
! author: Glenn Hammond
! date: 01/17/11; 04/03/13
!
! ************************************************************************** !
subroutine SetTranWeights(this)

  use Option_module

  implicit none
  
  class(pm_rt_type) :: this

  PetscReal :: flow_dt
  PetscReal :: flow_t0
  PetscReal :: flow_t1

  ! option%tran_time is the time at beginning of transport step
  flow_t0 = this%realization%patch%aux%Global%time_t
  flow_t1 = this%realization%patch%aux%Global%time_tpdt
  flow_dt = flow_t1-flow_t0
  this%tran_weight_t0 = max(0.d0,(this%option%time-flow_t0)/flow_dt)
  this%tran_weight_t1 = min(1.d0, &
                            (this%option%time+this%option%tran_dt-flow_t0)/ &
                            flow_dt)

end subroutine SetTranWeights

! ************************************************************************** !
!
! PMRTDestroy: Destroys RT process model
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTDestroy(this)

  use Reactive_Transport_module, only : RTDestroy

  implicit none
  
  class(pm_rt_type) :: this

#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRTDestroy()')
#endif
  
#ifndef SIMPLIFY 
  call RTDestroy(this%realization)
#endif

end subroutine PMRTDestroy
  
end module Process_Model_RT_class
