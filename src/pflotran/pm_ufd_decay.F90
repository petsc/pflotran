module PM_UFD_Decay_class

  use PM_Base_class
  use Realization_class
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  
  type, public, extends(pm_base_type) :: pm_ufd_decay_type
    class(realization_type), pointer :: realization
  contains
!geh: commented out subroutines can only be called externally
!    procedure, public :: Init => PMUFDDecayInit
    procedure, public :: Read => PMUFDDecayRead
!    procedure, public :: SetupSolvers => PMUFDDecaySetupSolvers
    procedure, public :: PMUFDDecaySetRealization
    procedure, public :: InitializeRun => PMUFDDecayInitializeRun
!!    procedure, public :: FinalizeRun => PMUFDDecayFinalizeRun
    procedure, public :: InitializeTimestep => PMUFDDecayInitializeTimestep
    procedure, public :: FinalizeTimestep => PMUFDDecayFinalizeTimestep
!    procedure, public :: PreSolve => PMUFDDecayPreSolve
    procedure, public :: Solve => PMUFDDecaySolve
!    procedure, public :: PostSolve => PMUFDDecayPostSolve
!    procedure, public :: AcceptSolution => PMUFDDecayAcceptSolution
!    procedure, public :: TimeCut => PMUFDDecayTimeCut
!    procedure, public :: UpdateSolution => PMUFDDecayUpdateSolution
!    procedure, public :: UpdateAuxvars => PMUFDDecayUpdateAuxvars
!    procedure, public :: Checkpoint => PMUFDDecayCheckpoint    
!    procedure, public :: Restart => PMUFDDecayRestart  
    procedure, public :: Destroy => PMUFDDecayDestroy
  end type pm_ufd_decay_type
  
  public :: PMUFDDecayCreate, &
            PMUFDDecayInit !, &
!            PMUFDDecaySetupSolvers, &
!            PMUFDDecayInitializeTimestepA, &
!            PMUFDDecayInitializeTimestepB, &
!            PMUFDDecayInitializeRun, &
!            PMUFDDecayUpdateSolution, &
!            PMUFDDecayUpdatePropertiesNI, &
!            PMUFDDecayTimeCut, &
!            PMUFDDecayDestroy
  
contains

! ************************************************************************** !

function PMUFDDecayCreate()
  ! 
  ! Creates the UFD decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type), pointer :: PMUFDDecayCreate
  
  allocate(PMUFDDecayCreate)
  nullify(PMUFDDecayCreate%realization)

  call PMBaseInit(PMUFDDecayCreate)

end function PMUFDDecayCreate

! ************************************************************************** !

subroutine PMUFDDecayRead(this,input)
  ! 
  ! Reads input file parameters associated with the ufd decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  type(input_type) :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  option => this%option
  
  error_string = 'UFD Decay'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    select case(trim(word))
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMUFDDecayRead

! ************************************************************************** !

subroutine PMUFDDecayInit(this)
  ! 
  ! Initializes variables associated with the UFD decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Reaction_Aux_module
  use Option_module

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  
  grid => this%realization%patch%grid
  option => this%realization%option
  reaction => this%realization%reaction
  
end subroutine PMUFDDecayInit

! ************************************************************************** !

subroutine PMUFDDecaySetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Realization_class

  implicit none
  
  class(pm_ufd_decay_type) :: this
  class(realization_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMUFDDecaySetRealization

! ************************************************************************** !

recursive subroutine PMUFDDecayInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none

  class(pm_ufd_decay_type) :: this
  

end subroutine PMUFDDecayInitializeRun

! ************************************************************************** !

subroutine PMUFDDecayInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Global_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," USED FUEL DISPOSITION DECAY MODEL ",43("="))')
  endif

end subroutine PMUFDDecayInitializeTimestep

! ************************************************************************** !

subroutine PMUFDDecayPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Grid_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  
end subroutine PMUFDDecayPreSolve

! ************************************************************************** !

subroutine PMUFDDecaySolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  !
  use Option_module
  use Reaction_Aux_module
  use Patch_module
  use Grid_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none

  class(pm_ufd_decay_type) :: this
  
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  ierr = 0
  
  option => this%realization%option
  reaction => this%realization%reaction
  patch => this%realization%patch
  grid => patch%grid
  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    
    ! perform the magic here
    
  enddo
  
end subroutine PMUFDDecaySolve

! ************************************************************************** !

subroutine PMUFDDecayPostSolve(this)
  ! 
  ! PMUFDDecayUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  ! 
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayPostSolve

! ************************************************************************** !

function PMUFDDecayAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscBool :: PMUFDDecayAcceptSolution
  
  ! do nothing
  PMUFDDecayAcceptSolution = PETSC_TRUE
  
end function PMUFDDecayAcceptSolution

! ************************************************************************** !

subroutine PMUFDDecayUpdatePropertiesTS(this)
  ! 
  ! Updates parameters/properties at each Newton iteration
  !
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayUpdatePropertiesTS

! ************************************************************************** !

subroutine PMUFDDecayTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine PMUFDDecayTimeCut

! ************************************************************************** !

subroutine PMUFDDecayFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayFinalizeTimestep

! ************************************************************************** !

subroutine PMUFDDecayUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscErrorCode :: ierr

end subroutine PMUFDDecayUpdateSolution  

! ************************************************************************** !

subroutine PMUFDDecayUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this

  this%option%io_buffer = 'PMUFDDecayUpdateAuxvars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMUFDDecayUpdateAuxvars   

! ************************************************************************** !

subroutine PMUFDDecayCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_ufd_decay_type) :: this
  PetscViewer :: viewer
  
end subroutine PMUFDDecayCheckpoint

! ************************************************************************** !

subroutine PMUFDDecayRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_ufd_decay_type) :: this
  PetscViewer :: viewer
  
end subroutine PMUFDDecayRestart

! ************************************************************************** !

recursive subroutine PMUFDDecayFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMUFDDecayFinalizeRun

! ************************************************************************** !

subroutine PMUFDDecayDestroy(this)
  ! 
  ! Destroys Subsurface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayDestroy
  
end module PM_UFD_Decay_class
