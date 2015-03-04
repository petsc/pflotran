module PM_TH_class

  use PM_Base_class
  use PM_Subsurface_class
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

  type, public, extends(pm_subsurface_type) :: pm_th_type
    class(communicator_type), pointer :: commN
  contains
    procedure, public :: Init => PMTHInit
    procedure, public :: Read => PMTHRead
    procedure, public :: InitializeTimestep => PMTHInitializeTimestep
    procedure, public :: Residual => PMTHResidual
    procedure, public :: Jacobian => PMTHJacobian
    procedure, public :: UpdateTimestep => PMTHUpdateTimestep
    procedure, public :: PreSolve => PMTHPreSolve
    procedure, public :: PostSolve => PMTHPostSolve
    procedure, public :: CheckUpdatePre => PMTHCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMTHCheckUpdatePost
    procedure, public :: TimeCut => PMTHTimeCut
    procedure, public :: UpdateSolution => PMTHUpdateSolution
    procedure, public :: UpdateAuxvars => PMTHUpdateAuxvars
    procedure, public :: MaxChange => PMTHMaxChange
    procedure, public :: ComputeMassBalance => PMTHComputeMassBalance
    procedure, public :: Destroy => PMTHDestroy
  end type pm_th_type
  
  public :: PMTHCreate
  
contains

! ************************************************************************** !

function PMTHCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_th_type), pointer :: PMTHCreate

  class(pm_th_type), pointer :: th_pm
  
#ifdef PM_TH_DEBUG
  print *, 'PMTHCreate()'
#endif  

  allocate(th_pm)

  nullify(th_pm%commN)

  call PMSubsurfaceCreate(th_pm)
  th_pm%name = 'PMTH'

  PMTHCreate => th_pm
  
end function PMTHCreate

! ************************************************************************** !

subroutine PMTHRead(this,input)
  ! 
  ! Reads input file parameters associated with the TH process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
 
  implicit none
  
  class(pm_th_type) :: this
  type(input_type) :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  option => this%option
  
  error_string = 'TH Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    select case(trim(word))
      case('FREEZING')
        option%use_th_freezing = PETSC_TRUE
        option%io_buffer = ' TH: using FREEZING submode!'
        call printMsg(option)
        ! Override the default setting for TH-mode with freezing
        call EOSWaterSetDensityPainter()
        call EOSWaterSetEnthalpyPainter()
      case('ICE_MODEL')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('PAINTER_EXPLICIT')
            option%ice_model = PAINTER_EXPLICIT
          case ('PAINTER_KARRA_IMPLICIT')
            option%ice_model = PAINTER_KARRA_IMPLICIT
          case ('PAINTER_KARRA_EXPLICIT')
            option%ice_model = PAINTER_KARRA_EXPLICIT
          case ('PAINTER_KARRA_EXPLICIT_NOCRYO')
            option%ice_model = PAINTER_KARRA_EXPLICIT_NOCRYO
          case ('DALL_AMICO')
            option%ice_model = DALL_AMICO
          case default
            option%io_buffer = 'Cannot identify the specificed ice model.' // &
             'Specify PAINTER_EXPLICIT or PAINTER_KARRA_IMPLICIT' // &
             ' or PAINTER_KARRA_EXPLICIT or PAINTER_KARRA_EXPLICIT_NOCRYO ' // &
             ' or DALL_AMICO.'
            call printErrMsg(option)
          end select
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMTHRead

! ************************************************************************** !

subroutine PMTHInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module 

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceInit(this)
  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%commN => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%commN => UnstructuredCommunicatorCreate()
  end select
  call this%commN%SetDM(this%realization%discretization%dm_nflowdof)

end subroutine PMTHInit

! ************************************************************************** !

subroutine PMTHSetupSolvers(this,solver)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/14

  use TH_module, only : THCheckUpdatePre, THCheckUpdatePost
  use Solver_module
  
  implicit none
  
  class(pm_subsurface_type) :: this
  type(solver_type) :: solver
  
  SNESLineSearch :: linesearch
  PetscErrorCode :: ierr
  
  call PMSubsurfaceSetupSolvers(this,solver)

  call SNESGetLineSearch(solver%snes, linesearch, ierr);CHKERRQ(ierr)
  if (dabs(this%option%pressure_dampening_factor) > 0.d0 .or. &
      dabs(this%option%pressure_change_limit) > 0.d0 .or. &
      dabs(this%option%temperature_change_limit) > 0.d0) then
    call SNESLineSearchSetPreCheck(linesearch, &
                                   THCheckUpdatePre, &
                                   this%realization,ierr);CHKERRQ(ierr)
  endif
  if (solver%check_post_convergence) then
    call SNESLineSearchSetPostCheck(linesearch, &
                                    THCheckUpdatePost, &
                                    this%realization,ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMTHSetupSolvers

! ************************************************************************** !

subroutine PMTHInitializeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THInitializeTimestep

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceInitializeTimestepA(this)

  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%icap_loc, &
                               this%realization%field%icap_loc)
  call this%comm1%LocalToLocal(this%realization%field%ithrm_loc, &
                               this%realization%field%ithrm_loc)
  call this%comm1%LocalToLocal(this%realization%field%iphas_loc, &
                               this%realization%field%iphas_loc)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," TH FLOW ",69("="))')
  endif
  
  call THInitializeTimestep(this%realization)
  call PMSubsurfaceInitializeTimestepB(this)
  
end subroutine PMTHInitializeTimestep

! ************************************************************************** !

subroutine PMTHPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Global_module

  implicit none
  
  class(pm_th_type) :: this

end subroutine PMTHPreSolve

! ************************************************************************** !

subroutine PMTHPostSolve(this)
  ! 
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_th_type) :: this
  
end subroutine PMTHPostSolve

! ************************************************************************** !

subroutine PMTHUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                              num_newton_iterations,tfac)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
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
    ut = min(up,utmp)
  endif
  dtt = fac * dt * (1.d0 + ut)

  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
  dtt = max(dtt,dt_min)
  dt = dtt
  
end subroutine PMTHUpdateTimestep

! ************************************************************************** !

subroutine PMTHResidual(this,snes,xx,r,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THResidual

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call THResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTHResidual

! ************************************************************************** !

subroutine PMTHJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THJacobian

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call THJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMTHJacobian

! ************************************************************************** !

subroutine PMTHCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THCheckUpdatePre

  implicit none
  
  class(pm_th_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call THCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

end subroutine PMTHCheckUpdatePre

! ************************************************************************** !

subroutine PMTHCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

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
  
  call THCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)

end subroutine PMTHCheckUpdatePost

! ************************************************************************** !

subroutine PMTHTimeCut(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THTimeCut

  implicit none
  
  class(pm_th_type) :: this
  
  call PMSubsurfaceTimeCut(this)
  call THTimeCut(this%realization)

end subroutine PMTHTimeCut

! ************************************************************************** !

subroutine PMTHUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THUpdateSolution, THUpdateSurfaceBC

  implicit none
  
  class(pm_th_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call THUpdateSolution(this%realization)
  if (this%option%surf_flow_on) &
    call THUpdateSurfaceBC(this%realization)

end subroutine PMTHUpdateSolution     

! ************************************************************************** !

subroutine PMTHUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use TH_module, only : THUpdateAuxVars
  
  implicit none
  
  class(pm_th_type) :: this

  call THUpdateAuxVars(this%realization)

end subroutine PMTHUpdateAuxvars   

! ************************************************************************** !

subroutine PMTHMaxChange(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THMaxChange

  implicit none
  
  class(pm_th_type) :: this
  
  call THMaxChange(this%realization)
    if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4," dtmpmx= ",1pe12.4)') &
      this%option%dpmax,this%option%dtmpmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4)') &
      this%option%dpmax,this%option%dtmpmax
  endif 

end subroutine PMTHMaxChange

! ************************************************************************** !

subroutine PMTHComputeMassBalance(this,mass_balance_array)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THComputeMassBalance

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call THComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMTHComputeMassBalance

! ************************************************************************** !

subroutine PMTHDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THDestroy

  implicit none
  
  class(pm_th_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call this%commN%Destroy()

  ! preserve this ordering
  call THDestroy(this%realization%patch)
  call PMSubsurfaceDestroy(this)

end subroutine PMTHDestroy

end module PM_TH_class
