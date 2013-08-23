module Process_Model_TH_class

  use Process_Model_Base_class
!geh: using TH_module here fails with gfortran (internal compiler error)
!  use TH_module
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

  type, public, extends(pm_base_type) :: pm_th_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    class(communicator_type), pointer :: commN
  contains
    procedure, public :: Init => PMTHInit
    procedure, public :: PMTHSetRealization
    procedure, public :: InitializeRun => PMTHInitializeRun
    procedure, public :: FinalizeRun => PMTHFinalizeRun
    procedure, public :: InitializeTimestep => PMTHInitializeTimestep
    procedure, public :: FinalizeTimestep => PMTHFinalizeTimeStep
    procedure, public :: Residual => PMTHResidual
    procedure, public :: Jacobian => PMTHJacobian
    procedure, public :: UpdateTimestep => PMTHUpdateTimestep
    procedure, public :: PreSolve => PMTHPreSolve
    procedure, public :: PostSolve => PMTHPostSolve
    procedure, public :: AcceptSolution => PMTHAcceptSolution
    procedure, public :: CheckUpdatePre => PMTHCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMTHCheckUpdatePost
    procedure, public :: TimeCut => PMTHTimeCut
    procedure, public :: UpdateSolution => PMTHUpdateSolution
    procedure, public :: MaxChange => PMTHMaxChange
    procedure, public :: ComputeMassBalance => PMTHComputeMassBalance
    procedure, public :: Destroy => PMTHDestroy
  end type pm_th_type
  
  public :: PMTHCreate
  
contains


! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
function PMTHCreate()

  implicit none
  
  class(pm_th_type), pointer :: PMTHCreate

  class(pm_th_type), pointer :: th_pm
  
#ifdef PM_TH_DEBUG
  print *, 'PMTHCreate()'
#endif  

  allocate(th_pm)
  nullify(th_pm%option)
  nullify(th_pm%output_option)
  nullify(th_pm%realization)
  nullify(th_pm%comm1)
  nullify(th_pm%commN)

  call PMBaseCreate(th_pm)
  th_pm%name = 'PMTH'

  PMTHCreate => th_pm
  
end function PMTHCreate

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHInit(this)

#ifndef SIMPLIFY  
  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
#endif

  implicit none
  
  class(pm_th_type) :: this

#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%Init()')
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

end subroutine PMTHInit

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHSetRealization(this,realization)

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_th_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTHSetRealization%SetRealization()')
#endif
  
  this%realization => realization
  this%realization_base => realization

  this%solution_vec = realization%field%flow_xx
  this%residual_vec = realization%field%flow_r
  
end subroutine PMTHSetRealization

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHInitializeTimestep(this)

  use TH_module, only : THInitializeTimestep
  use Global_module
  
  implicit none
  
  class(pm_th_type) :: this

#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTHInitializeTimestep%InitializeTimestep()')
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
    write(*,'(/,2("=")," TH FLOW ",62("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call THInitializeTimestep(this%realization)
  
end subroutine PMTHInitializeTimestep

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHPreSolve(this)

  use Global_module

  implicit none
  
  class(pm_th_type) :: this

#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTHPreSolve%PreSolve()')
#endif

end subroutine PMTHPreSolve

! ************************************************************************** !
!
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMTHPostSolve(this)

  use Global_module

  implicit none
  
  class(pm_th_type) :: this
  
#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTH%PostSolve()')
#endif
  
end subroutine PMTHPostSolve

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHFinalizeTimestep(this)

  use TH_module, only : THMaxChange
  use Global_module

  implicit none
  
  class(pm_th_type) :: this
  
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call THMaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4)') this%option%dpmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') &
      this%option%dpmax
  endif  
  
end subroutine PMTHFinalizeTimestep

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
function PMTHAcceptSolution(this)

  implicit none
  
  class(pm_th_type) :: this
  
  PetscBool :: PMTHAcceptSolution
  
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTHs%AcceptSolution()')
#endif
  ! do nothing
  PMTHAcceptSolution = PETSC_TRUE
  
end function PMTHAcceptSolution

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHUpdateTimestep(this,dt,dt_max,iacceleration, &
                              num_newton_iterations,tfac)

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
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
  
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%UpdateTimestep()')
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
      
  dt = dtt
  
end subroutine PMTHUpdateTimestep

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
recursive subroutine PMTHInitializeRun(this)

  use TH_module, only : THUpdateSolution

  implicit none
  
  class(pm_th_type) :: this
  
 
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call THUpdateSolution(this%realization)
#endif  
    
end subroutine PMTHInitializeRun

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
recursive subroutine PMTHFinalizeRun(this)

  implicit none
  
  class(pm_th_type) :: this
  
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMTHFinalizeRun

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHResidual(this,snes,xx,r,ierr)

  use TH_module, only : THResidual

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%Residual()')
#endif
  
  call THResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTHResidual

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHJacobian(this,snes,xx,A,B,flag,ierr)

  use TH_module, only : THJacobian

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%Jacobian()')
#endif
  
  call THJacobian(snes,xx,A,B,flag,this%realization,ierr)

end subroutine PMTHJacobian

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHCheckUpdatePre(this,line_search,P,dP,changed,ierr)

  use TH_module, only : THCheckUpdatePre

  implicit none
  
  class(pm_th_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTH%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY  
  call THCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMTHCheckUpdatePre
    
! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)

  use TH_module, only : THCheckUpdatePost

  implicit none
  
  class(pm_th_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTH%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call THCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMTHCheckUpdatePost
  
! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHTimeCut(this)

  use TH_module, only : THTimeCut

  implicit none
  
  class(pm_th_type) :: this
  
#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTH%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call THTimeCut(this%realization)

end subroutine PMTHTimeCut
    
! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHUpdateSolution(this)

  use TH_module, only : THUpdateSolution
  use Condition_module

  implicit none
  
  class(pm_th_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTH%UpdateSolution()')
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
  call THUpdateSolution(this%realization)

end subroutine PMTHUpdateSolution     

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHMaxChange(this)

  use TH_module, only : THMaxChange

  implicit none
  
  class(pm_th_type) :: this
  
#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTH%MaxChange()')
#endif

  call THMaxChange(this%realization)

end subroutine PMTHMaxChange
    
! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHComputeMassBalance(this,mass_balance_array)

  use TH_module, only : THComputeMassBalance

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTH%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY  
  call THComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMTHComputeMassBalance

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/90/13
! ************************************************************************** !
subroutine PMTHDestroy(this)

  use TH_module, only : THDestroy

  implicit none
  
  class(pm_th_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_TH_DEBUG  
  call printMsg(this%option,'PMTHDestroy()')
#endif

#ifndef SIMPLIFY 
  call THDestroy(this%realization%patch)
#endif
  call this%comm1%Destroy()
  call this%commN%Destroy()

end subroutine PMTHDestroy

end module Process_Model_TH_class
