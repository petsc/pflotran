module Process_Model_Richards_class

  use Process_Model_Base_class
  use Richards_module
#ifdef SIMPLIFY  
  use Realization_class
#endif  
  use Communicator_Base_module
  
  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(process_model_base_type) :: process_model_richards_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm
  contains
    procedure, public :: Init => PMRichardsInit
    procedure, public :: PMRichardsSetRealization
    procedure, public :: InitializeRun => PMRichardsInitializeRun
    procedure, public :: FinalizeRun => PMRichardsFinalizeRun
    procedure, public :: InitializeTimeStep => PMRichardsInitializeTimestep
    procedure, public :: Residual => PMRichardsResidual
    procedure, public :: Jacobian => PMRichardsJacobian
    procedure, public :: UpdateTimestep => PMRichardsUpdateTimestep
    procedure, public :: UpdatePreSolve => PMRichardsUpdatePreSolve
    procedure, public :: CheckUpdatePre => PMRichardsCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRichardsCheckUpdatePost
    procedure, public :: TimeCut => PMRichardsTimeCut
    procedure, public :: UpdateSolution => PMRichardsUpdateSolution
    procedure, public :: MaxChange => PMRichardsMaxChange
    procedure, public :: ComputeMassBalance => PMRichardsComputeMassBalance
    procedure, public :: Destroy => PMRichardsDestroy
  end type process_model_richards_type
  
  type :: realization_type
    PetscInt :: i
  end type realization_type
  
  public :: PMRichardsCreate
  
contains

! ************************************************************************** !
!
! PMRichardsCreate: Creates Richards process models shell
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function PMRichardsCreate()

  implicit none
  
  class(process_model_richards_type), pointer :: PMRichardsCreate

  class(process_model_richards_type), pointer :: richards_pm
  
  allocate(richards_pm)
  nullify(richards_pm%realization)
  nullify(richards_pm%comm)

  PMRichardsCreate => richards_pm
  
end function PMRichardsCreate

! ************************************************************************** !
!
! PMRichardsInit: Initializes variables associated with Richard
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsInit(this)

  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
  
  implicit none
  
  class(process_model_richards_type) :: this

#ifdef SIMPLIFY  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      this%comm => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%comm => UnstructuredCommunicatorCreate()
  end select
  call this%comm%SetDM(this%realization%discretization%dm_1dof)
#endif

end subroutine PMRichardsInit

! ************************************************************************** !
!
! PMRichardsSetRealization: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsSetRealization(this,realization)

  use Realization_Base_class  

  implicit none
  
  class(process_model_richards_type) :: this
  class(realization_type), pointer :: realization

  this%realization => realization
  
end subroutine PMRichardsSetRealization

! ************************************************************************** !
!
! PMRichardsInitializeTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsInitializeTimestep(this)

  implicit none
  
  class(process_model_richards_type) :: this

#ifdef SIMPLIFY  
  call RichardsSetup(this%realization)
#endif
  
end subroutine PMRichardsInitializeTimestep

! ************************************************************************** !
!
! PMRichardsUpdatePreSolve: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsUpdatePreSolve(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
#ifdef SIMPLIFY  
  ! update porosity
  call this%comm%LocalToLocal(this%realization%field%porosity_loc, &
                              this%realization%field%porosity_loc)
#endif
  
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," RICHARDS FLOW ",63("="))')
  endif
  
end subroutine PMRichardsUpdatePreSolve

! ************************************************************************** !
!
! PMRichardsUpdateTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)

  implicit none
  
  class(process_model_richards_type) :: this
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
!
! PMRichardsInitializeRun: Initializes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMRichardsInitializeRun(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
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
!
! PMRichardsFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMRichardsFinalizeRun(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMRichardsFinalizeRun

! ************************************************************************** !
!
! PMRichardsResidual: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsResidual(this,snes,xx,r,ierr)

  implicit none
  
  class(process_model_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef SIMPLIFY  
  call RichardsResidual(snes,xx,r,this%realization,ierr)
#endif

end subroutine PMRichardsResidual

! ************************************************************************** !
!
! PMRichardsJacobian: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsJacobian(this,snes,xx,A,B,flag,ierr)

  implicit none
  
  class(process_model_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef SIMPLIFY  
  call RichardsJacobian(snes,xx,A,B,flag,this%realization,ierr)
#endif

end subroutine PMRichardsJacobian
    
! ************************************************************************** !
!
! PMRichardsCheckUpdatePre: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsCheckUpdatePre(this,line_search,P,dP,changed,ierr)

  implicit none
  
  class(process_model_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef SIMPLIFY  
  call RichardsCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMRichardsCheckUpdatePre
    
! ************************************************************************** !
!
! PMRichardsCheckUpdatePost: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  implicit none
  
  class(process_model_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef SIMPLIFY  
  call RichardsCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMRichardsCheckUpdatePost
  
! ************************************************************************** !
!
! PMRichardsTimeCut: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsTimeCut(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
#ifdef SIMPLIFY  
  call RichardsTimeCut(this%realization)
#endif

end subroutine PMRichardsTimeCut
    
! ************************************************************************** !
!
! PMRichardsUpdateSolution: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsUpdateSolution(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
#ifdef SIMPLIFY  
  call RichardsUpdateSolution(this%realization)
#endif

end subroutine PMRichardsUpdateSolution     

! ************************************************************************** !
!
! PMRichardsMaxChange: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsMaxChange(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
#ifdef SIMPLIFY  
  call RichardsMaxChange(this%realization)
#endif

end subroutine PMRichardsMaxChange
    
! ************************************************************************** !
!
! PMRichardsComputeMassBalance: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsComputeMassBalance(this,mass_balance_array)

  implicit none
  
  class(process_model_richards_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef SIMPLIFY  
  call RichardsComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMRichardsComputeMassBalance

! ************************************************************************** !
!
! PMRichardsDestroy: Destroys Richards process model
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsDestroy(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef SIMPLIFY 
  call RichardsDestroy(this%realization)
#endif
  call this%comm%Destroy()
  
end subroutine PMRichardsDestroy
  
end module Process_Model_Richards_class
