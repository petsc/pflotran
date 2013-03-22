module Process_Model_RT_class

  use Process_Model_Base_class
  use Reactive_Transport_module
  
  use Realization_class
  use Communicator_Base_module  
  
  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(process_model_base_type) :: process_model_rt_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    class(communicator_type), pointer :: commN
  contains
    procedure, public :: Init => PMRTInit
    procedure, public :: PMRTSetRealization
    procedure, public :: InitializeRun => PMRTInitializeRun
    procedure, public :: FinalizeRun => PMRTFinalizeRun
    procedure, public :: InitializeTimeStep => PMRTInitializeTimestep
    procedure, public :: Residual => PMRTResidual
    procedure, public :: Jacobian => PMRTJacobian
    procedure, public :: UpdateTimestep => PMRTUpdateTimestep
    procedure, public :: UpdatePreSolve => PMRTUpdatePreSolve
    procedure, public :: CheckUpdatePre => PMRTCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRTCheckUpdatePost
    procedure, public :: TimeCut => PMRTTimeCut
    procedure, public :: UpdateSolution => PMRTUpdateSolution
    procedure, public :: MaxChange => PMRTMaxChange
    procedure, public :: ComputeMassBalance => PMRTComputeMassBalance
    procedure, public :: Destroy => PMRTDestroy
  end type process_model_rt_type

  type :: realization_type
    PetscInt :: i
  end type realization_type
  
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
  
  class(process_model_rt_type), pointer :: PMRTCreate

  class(process_model_rt_type), pointer :: rt_pm
  
  allocate(rt_pm)
  nullify(rt_pm%realization)
  nullify(rt_pm%comm1)
  nullify(rt_pm%commN)

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

  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
  
  implicit none
  
  class(process_model_rt_type) :: this

#ifdef SIMPLIFY  
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

  use Realization_Base_class  

  implicit none
  
  class(process_model_rt_type) :: this
  class(realization_type), pointer :: realization

  this%realization => realization
  
end subroutine PMRTSetRealization

! ************************************************************************** !
!
! PMRTInitializeTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTInitializeTimestep(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
#ifdef SIMPLIFY 
  call RTSetup(this%realization)
#endif

end subroutine PMRTInitializeTimestep

! ************************************************************************** !
!
! PMRTUpdatePreSolve: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTUpdatePreSolve(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
#ifdef SIMPLIFY 
  ! update tortuosity
  call this%comm1%LocalToLocal(this%realization%field%tortuosity_loc, &
                               this%realization%field%tortuosity_loc)
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," REACTIVE TRANSPORT ",58("="))')
  endif
  
end subroutine PMRTUpdatePreSolve

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
  
  class(process_model_rt_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: dtt
  
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

  implicit none
  
  class(process_model_rt_type) :: this
  
  ! restart
#if 0  
  if (transport_read .and. option%overwrite_restart_transport) then
    call CondControlAssignTranInitCond(realization)  
  endif
#endif  
  
#ifdef SIMPLIFY 
  call RTUpdateSolution(this%realization)
  
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
#endif
  ! check on MAX_STEPS < 0 to quit after initialization.
    
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
  
  class(process_model_rt_type) :: this
  
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

  implicit none
  
  class(process_model_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef SIMPLIFY 
  call RTResidual(snes,xx,r,this%realization,ierr)
#endif

end subroutine PMRTResidual

! ************************************************************************** !
!
! PMRTJacobian: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTJacobian(this,snes,xx,A,B,flag,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef SIMPLIFY 
  call RTJacobian(snes,xx,A,B,flag,this%realization,ierr)
#endif

end subroutine PMRTJacobian
    
! ************************************************************************** !
!
! PMRTCheckUpdatePre: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTCheckUpdatePre(this,line_search,P,dP,changed,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef SIMPLIFY 
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
  implicit none
  
  class(process_model_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
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

  implicit none
  
  class(process_model_rt_type) :: this
  
#ifdef SIMPLIFY 
  call RTTimeCut(this%realization)
#endif

end subroutine PMRTTimeCut
    
! ************************************************************************** !
!
! PMRTUpdateSolution: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTUpdateSolution(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
#ifdef SIMPLIFY 
  call RTUpdateSolution(this%realization)
#endif

end subroutine PMRTUpdateSolution     

! ************************************************************************** !
!
! PMRTMaxChange: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTMaxChange(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
#ifdef SIMPLIFY 
  call RTMaxChange(this%realization)
#endif

end subroutine PMRTMaxChange
    
! ************************************************************************** !
!
! PMRTComputeMassBalance: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTComputeMassBalance(this,mass_balance_array)

  implicit none
  
  class(process_model_rt_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef SIMPLIFY 
  call RTComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMRTComputeMassBalance

! ************************************************************************** !
!
! PMRTDestroy: Destroys RT process model
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTDestroy(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
#ifdef SIMPLIFY 
  call RTDestroy(this%realization)
#endif

end subroutine PMRTDestroy
  
end module Process_Model_RT_class
