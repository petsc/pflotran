module PM_THC_class
#include "finclude/petscmat.h"
  use petscmat
  use PM_Base_class
!geh: using TH_module here fails with gfortran (internal compiler error)
!  use TH_module
  use Realization_class
  use Communicator_Base_module
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_thc_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    class(communicator_type), pointer :: commN
  contains
    procedure, public :: Init => PMTHCInit
    procedure, public :: PMTHCSetRealization
    procedure, public :: InitializeRun => PMTHCInitializeRun
    procedure, public :: FinalizeRun => PMTHCFinalizeRun
    procedure, public :: InitializeTimestep => PMTHCInitializeTimestep
    procedure, public :: FinalizeTimestep => PMTHCFinalizeTimeStep
    procedure, public :: Residual => PMTHCResidual
    procedure, public :: Jacobian => PMTHCJacobian
    procedure, public :: UpdateTimestep => PMTHCUpdateTimestep
    procedure, public :: PreSolve => PMTHCPreSolve
    procedure, public :: PostSolve => PMTHCPostSolve
    procedure, public :: AcceptSolution => PMTHCAcceptSolution
    procedure, public :: CheckUpdatePre => PMTHCCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMTHCCheckUpdatePost
    procedure, public :: TimeCut => PMTHCTimeCut
    procedure, public :: UpdateSolution => PMTHCUpdateSolution
    procedure, public :: MaxChange => PMTHCMaxChange
    procedure, public :: ComputeMassBalance => PMTHCComputeMassBalance
    procedure, public :: Destroy => PMTHCDestroy
  end type pm_thc_type
  
  public :: PMTHCCreate
  
contains

! ************************************************************************** !

function PMTHCCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_thc_type), pointer :: PMTHCCreate

  class(pm_thc_type), pointer :: thc_pm
  
#ifdef PM_THC_DEBUG
  print *, 'PMTHCCreate()'
#endif  

  allocate(thc_pm)
  nullify(thc_pm%option)
  nullify(thc_pm%output_option)
  nullify(thc_pm%realization)
  nullify(thc_pm%comm1)
  nullify(thc_pm%commN)

  call PMBaseCreate(thc_pm)
  thc_pm%name = 'PMTHC'

  PMTHCCreate => thc_pm
  
end function PMTHCCreate

! ************************************************************************** !

subroutine PMTHCInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

#ifndef SIMPLIFY  
  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
#endif

  implicit none
  
  class(pm_thc_type) :: this

#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHC%Init()')
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
  call this%commN%SetDM(this%realization%discretization%dm_nflowdof)
#endif

end subroutine PMTHCInit

! ************************************************************************** !

subroutine PMTHCSetRealization(this,realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_thc_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHCSetRealization%SetRealization()')
#endif
  
  this%realization => realization
  this%realization_base => realization

  this%solution_vec = realization%field%flow_xx
  this%residual_vec = realization%field%flow_r
  
end subroutine PMTHCSetRealization

! ************************************************************************** !

subroutine PMTHCInitializeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCInitializeTimestep
  use Global_module
  
  implicit none
  
  class(pm_thc_type) :: this

#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHCInitializeTimestep%InitializeTimestep()')
#endif

  this%option%flow_dt = this%option%dt

#ifndef SIMPLIFY  
  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%porosity_loc, &
                              this%realization%field%porosity_loc)
  call this%comm1%LocalToLocal(this%realization%field%tortuosity_loc, &
                              this%realization%field%tortuosity_loc)
  call this%comm1%LocalToLocal(this%realization%field%icap_loc, &
                              this%realization%field%icap_loc)
  call this%comm1%LocalToLocal(this%realization%field%ithrm_loc, &
                              this%realization%field%ithrm_loc)
  call this%comm1%LocalToLocal(this%realization%field%iphas_loc, &
                              this%realization%field%iphas_loc)
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," THC FLOW ",68("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call THCInitializeTimestep(this%realization)
  
end subroutine PMTHCInitializeTimestep

! ************************************************************************** !

subroutine PMTHCPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Global_module

  implicit none
  
  class(pm_thc_type) :: this

#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHCPreSolve%PreSolve()')
#endif

end subroutine PMTHCPreSolve

! ************************************************************************** !

subroutine PMTHCPostSolve(this)
  ! 
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_thc_type) :: this
  
#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHC%PostSolve()')
#endif
  
end subroutine PMTHCPostSolve

! ************************************************************************** !

subroutine PMTHCFinalizeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCMaxChange
  use Global_module

  implicit none
  
  class(pm_thc_type) :: this
  
#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHC%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call THCMaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4, " dtmpmx= ",1pe12.4)') this%option%dpmax, this%option%dtmpmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
      & "dtmpmx= ",1pe12.4)') &
      this%option%dpmax, this%option%dtmpmax
  endif  
  
end subroutine PMTHCFinalizeTimestep

! ************************************************************************** !

function PMTHCAcceptSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_thc_type) :: this
  
  PetscBool :: PMTHCAcceptSolution
  
#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHCs%AcceptSolution()')
#endif
  ! do nothing
  PMTHCAcceptSolution = PETSC_TRUE
  
end function PMTHCAcceptSolution

! ************************************************************************** !

subroutine PMTHCUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                              num_newton_iterations,tfac)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_thc_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
  PetscReal :: uus
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHC%UpdateTimestep()')
#endif
  
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    ut = 0.d0
  else
    up = this%option%dpmxe/(this%option%dpmax+0.1)
    utmp = this%option%dtmpmxe/(this%option%dtmpmax+1.d-5)
    uus= this%option%dsmxe/(this%option%dsmax+1.d-6)
    ut = min(up,utmp,uus)
  endif
  dtt = fac * dt * (1.d0 + ut)

  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
  dtt = max(dtt,dt_min)
  dt = dtt
  
end subroutine PMTHCUpdateTimestep

! ************************************************************************** !

recursive subroutine PMTHCInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCUpdateSolution

  implicit none
  
  class(pm_thc_type) :: this
  
 
#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHC%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call THCUpdateSolution(this%realization)
#endif  
    
end subroutine PMTHCInitializeRun

! ************************************************************************** !

recursive subroutine PMTHCFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_thc_type) :: this
  
#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHC%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMTHCFinalizeRun

! ************************************************************************** !

subroutine PMTHCResidual(this,snes,xx,r,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCResidual

  implicit none
  
  class(pm_thc_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHC%Residual()')
#endif
  
  call THCResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTHCResidual

! ************************************************************************** !

subroutine PMTHCJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCJacobian

  implicit none
  
  class(pm_thc_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
#ifdef PM_THC_DEBUG
  call printMsg(this%option,'PMTHC%Jacobian()')
#endif
  
  call THCJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMTHCJacobian

! ************************************************************************** !

subroutine PMTHCCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCCheckUpdatePre

  implicit none
  
  class(pm_thc_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHC%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY
  call THCCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMTHCCheckUpdatePre

! ************************************************************************** !

subroutine PMTHCCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCCheckUpdatePost

  implicit none
  
  class(pm_thc_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHC%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call THCCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMTHCCheckUpdatePost

! ************************************************************************** !

subroutine PMTHCTimeCut(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCTimeCut

  implicit none
  
  class(pm_thc_type) :: this
  
#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHC%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call THCTimeCut(this%realization)

end subroutine PMTHCTimeCut

! ************************************************************************** !

subroutine PMTHCUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCUpdateSolution
  use Condition_module

  implicit none
  
  class(pm_thc_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHC%UpdateSolution()')
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
  call THCUpdateSolution(this%realization)

end subroutine PMTHCUpdateSolution     

! ************************************************************************** !

subroutine PMTHCMaxChange(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCMaxChange

  implicit none
  
  class(pm_thc_type) :: this
  
#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHC%MaxChange()')
#endif

  call THCMaxChange(this%realization)

end subroutine PMTHCMaxChange

! ************************************************************************** !

subroutine PMTHCComputeMassBalance(this,mass_balance_array)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCComputeMassBalance

  implicit none
  
  class(pm_thc_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHC%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY  
  call THCComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMTHCComputeMassBalance

! ************************************************************************** !

subroutine PMTHCDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use THC_module, only : THCDestroy

  implicit none
  
  class(pm_thc_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_THC_DEBUG  
  call printMsg(this%option,'PMTHCDestroy()')
#endif

#ifndef SIMPLIFY 
  call THCDestroy(this%realization%patch)
#endif
  call this%comm1%Destroy()
  call this%commN%Destroy()

end subroutine PMTHCDestroy

end module PM_THC_class
