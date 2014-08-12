module PM_Subsurface_class

  use PM_Base_class
!geh: using Subsurface_module here fails with gfortran (internal compiler error)
!  use Subsurface_module
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

  type, public, extends(pm_base_type) :: pm_subsurface_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    PetscBool :: transient_permeability
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Init => PMSubsurfaceInit
    procedure, public :: PMSubsurfaceSetRealization
    procedure, public :: InitializeRun => PMSubsurfaceInitializeRun
!    procedure, public :: FinalizeRun => PMSubsurfaceFinalizeRun
!    procedure, public :: InitializeTimestep => PMSubsurfaceInitializeTimestep
    procedure, public :: FinalizeTimestep => PMSubsurfaceFinalizeTimestep
    procedure, public :: PreSolve => PMSubsurfacePreSolve
    procedure, public :: PostSolve => PMSubsurfacePostSolve
    procedure, public :: AcceptSolution => PMSubsurfaceAcceptSolution
!    procedure, public :: TimeCut => PMSubsurfaceTimeCut
!    procedure, public :: UpdateSolution => PMSubsurfaceUpdateSolution
    procedure, public :: UpdateAuxvars => PMSubsurfaceUpdateAuxvars
    procedure, public :: Checkpoint => PMSubsurfaceCheckpoint    
    procedure, public :: Restart => PMSubsurfaceRestart  
!    procedure, public :: Destroy => PMSubsurfaceDestroy
  end type pm_subsurface_type
  
  public :: PMSubsurfaceCreate, &
            PMSubsurfaceInit, &
            PMSubsurfaceInitializeTimestep, &
            PMSubsurfaceInitializeRun, &
            PMSubsurfaceUpdateSolution, &
            PMSubsurfaceUpdatePropertiesNI, &
            PMSubsurfaceTimeCut, &
            PMSubsurfaceDestroy
  
contains

! ************************************************************************** !

subroutine PMSubsurfaceCreate(this)
  ! 
  ! Intializes shared members of subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_type) :: this
  
  nullify(this%realization)
  nullify(this%comm1)
  this%transient_permeability = PETSC_FALSE
  
  call PMBaseCreate(this)

end subroutine PMSubsurfaceCreate

! ************************************************************************** !

subroutine PMSubsurfaceInit(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 

  implicit none
  
  class(pm_subsurface_type) :: this
  
  PetscErrorCode :: ierr

  ! set the communicator
  this%comm1 => this%realization%comm1
  
end subroutine PMSubsurfaceInit

! ************************************************************************** !

subroutine PMSubsurfaceSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_subsurface_type) :: this
  class(realization_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then 
    this%solution_vec = realization%field%flow_xx_faces
    this%residual_vec = realization%field%flow_r_faces
  else
    this%solution_vec = realization%field%flow_xx
    this%residual_vec = realization%field%flow_r
  endif
  
end subroutine PMSubsurfaceSetRealization

! ************************************************************************** !

recursive subroutine PMSubsurfaceInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 
  use Condition_Control_module

  implicit none
  
  class(pm_subsurface_type) :: this

  ! overridden in pm_general only
  
  ! restart
  if (this%option%restart_flag .and. this%option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(this%realization)
!geh: for testing only.  In general, we only revert parameter, not flow.
!    call CondControlAssignFlowInitCond(this%realization)
!    call this%UpdateAuxVars()
  endif

  call this%UpdateSolution()  
  
    
end subroutine PMSubsurfaceInitializeRun

! ************************************************************************** !

subroutine PMSubsurfaceInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_subsurface_type) :: this

  this%option%flow_dt = this%option%dt
!geh:remove
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,POROSITY,0)
                                 
  if (this%transient_permeability) then
  !geh:remove
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
  endif

  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
end subroutine PMSubsurfaceInitializeTimestep

! ************************************************************************** !

subroutine PMSubsurfacePreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module

  implicit none
  
  class(pm_subsurface_type) :: this
  
  this%option%io_buffer = 'PMSubsurfacePreSolve() must be extended.'
  call printErrMsg(this%option)  

end subroutine PMSubsurfacePreSolve

! ************************************************************************** !

subroutine PMSubsurfacePostSolve(this)
  ! 
  ! PMSubsurfaceUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_subsurface_type) :: this
  
  this%option%io_buffer = 'PMSubsurfacePostSolve() must be extended.'
  call printErrMsg(this%option)  
  
end subroutine PMSubsurfacePostSolve

! ************************************************************************** !

function PMSubsurfaceAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_type) :: this
  
  PetscBool :: PMSubsurfaceAcceptSolution
  
  ! do nothing
  PMSubsurfaceAcceptSolution = PETSC_TRUE
  
end function PMSubsurfaceAcceptSolution

! ************************************************************************** !

subroutine PMSubsurfaceUpdatePropertiesNI(this)
  ! 
  ! Updates parameters/properties at each Newton iteration
  !
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_type) :: this
  
  call RealizationUpdatePropertiesNI(this%realization)

end subroutine PMSubsurfaceUpdatePropertiesNI

! ************************************************************************** !

subroutine PMSubsurfaceTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 

  implicit none
  
  class(pm_subsurface_type) :: this
  
  this%option%flow_dt = this%option%dt

end subroutine PMSubsurfaceTimeCut

! ************************************************************************** !

subroutine PMSubsurfaceFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module

  implicit none
  
  class(pm_subsurface_type) :: this

  if (this%option%ntrandof > 0) then 
    ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call this%MaxChange()
  
end subroutine PMSubsurfaceFinalizeTimestep

! ************************************************************************** !

subroutine PMSubsurfaceUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Condition_module

  implicit none
  
  class(pm_subsurface_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

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

end subroutine PMSubsurfaceUpdateSolution  

! ************************************************************************** !

subroutine PMSubsurfaceUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_type) :: this

  this%option%io_buffer = 'PMSubsurfaceUpdateAuxvars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMSubsurfaceUpdateAuxvars   

! ************************************************************************** !

subroutine PMSubsurfaceCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_subsurface_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMSubsurfaceCheckpoint

! ************************************************************************** !

subroutine PMSubsurfaceRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_subsurface_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModel(viewer,this%realization)
  call this%UpdateAuxVars()
  call this%UpdateSolution()
  
end subroutine PMSubsurfaceRestart

! ************************************************************************** !

recursive subroutine PMSubsurfaceFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMSubsurfaceFinalizeRun

! ************************************************************************** !

subroutine PMSubsurfaceDestroy(this)
  ! 
  ! Destroys Subsurface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_type) :: this
  
  ! destroyed in realization
  nullify(this%comm1)
  
end subroutine PMSubsurfaceDestroy
  
end module PM_Subsurface_class
