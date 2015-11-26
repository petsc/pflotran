module PM_Subsurface_class

  use PM_Base_class
!geh: using Init_Subsurface_module here fails with gfortran (internal compiler error)
!  use Init_Subsurface_module
  use Realization_class
  use Communicator_Base_module
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"

  type, public, extends(pm_base_type) :: pm_subsurface_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    PetscBool :: transient_permeability
    PetscBool :: store_porosity_for_ts_cut
    PetscBool :: store_porosity_for_transport
    PetscBool :: check_post_convergence
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Setup => PMSubsurfaceSetup
    procedure, public :: SetupSolvers => PMSubsurfaceSetupSolvers
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
    procedure, public :: CheckpointBinary => PMSubsurfaceCheckpointBinary
    procedure, public :: CheckpointHDF5 => PMSubsurfaceCheckpointHDF5
    procedure, public :: RestartBinary => PMSubsurfaceRestartBinary
    procedure, public :: RestartHDF5 => PMSubsurfaceRestartHDF5
!    procedure, public :: Destroy => PMSubsurfaceDestroy
  end type pm_subsurface_type
  
  public :: PMSubsurfaceCreate, &
            PMSubsurfaceSetup, &
            PMSubsurfaceSetupSolvers, &
            PMSubsurfaceInitializeTimestepA, &
            PMSubsurfaceInitializeTimestepB, &
            PMSubsurfaceInitializeRun, &
            PMSubsurfaceUpdateSolution, &
            PMSubsurfaceUpdatePropertiesNI, &
            PMSubsurfaceTimeCut, &
            PMSubsurfaceCheckpointBinary, &
            PMSubsurfaceCheckpointHDF5, &
            PMSubsurfaceRestartBinary, &
            PMSubsurfaceRestartHDF5, &
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
  this%store_porosity_for_ts_cut = PETSC_FALSE
  this%store_porosity_for_transport = PETSC_FALSE
  this%check_post_convergence = PETSC_FALSE
  
  call PMBaseInit(this)

end subroutine PMSubsurfaceCreate

! ************************************************************************** !

subroutine PMSubsurfaceSetup(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module 

  implicit none
  
  class(pm_subsurface_type) :: this
  
  PetscErrorCode :: ierr

  ! set the communicator
  this%comm1 => this%realization%comm1
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity .or. &
        this%realization%reaction%update_tortuosity .or. &
        this%realization%reaction%update_permeability .or. &
        this%realization%reaction%update_mnrl_surf_with_porosity) then
      this%store_porosity_for_ts_cut = PETSC_TRUE
      this%store_porosity_for_transport = PETSC_TRUE
    endif
  endif
  if (this%option%flow%transient_porosity) then
    this%store_porosity_for_ts_cut = PETSC_TRUE
    if (this%option%ntrandof > 0) then
      this%store_porosity_for_transport = PETSC_TRUE
    endif
  endif
  
end subroutine PMSubsurfaceSetup

! ************************************************************************** !

subroutine PMSubsurfaceSetupSolvers(this,solver)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/14

  use Solver_module
  
  implicit none
  
  class(pm_subsurface_type) :: this
  type(solver_type) :: solver
  
  PetscErrorCode :: ierr
  
end subroutine PMSubsurfaceSetupSolvers

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

  this%solution_vec = realization%field%flow_xx
  this%residual_vec = realization%field%flow_r
  
end subroutine PMSubsurfaceSetRealization

! ************************************************************************** !

recursive subroutine PMSubsurfaceInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 
  use Condition_Control_module
  use Material_module
  use Variables_module, only : POROSITY
  use Material_Aux_class, only : POROSITY_MINERAL, POROSITY_CURRENT

  implicit none
  
  class(pm_subsurface_type) :: this
  PetscBool :: update_initial_porosity

  ! overridden in pm_general only
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity) then
      call RealizationCalcMineralPorosity(this%realization)
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc, &
                                   POROSITY,POROSITY_MINERAL)
      call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                    this%realization%field%porosity0)
    endif
  endif
  
  ! restart
  if (this%option%restart_flag .and. this%option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(this%realization)
!geh: for testing only.  In general, we only revert parameter, not flow.
!    call CondControlAssignFlowInitCond(this%realization)
!    call this%UpdateAuxVars()
  endif
  ! update material properties that are a function of mineral vol fracs
  update_initial_porosity = PETSC_TRUE
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity .or. &
        this%realization%reaction%update_tortuosity .or. &
        this%realization%reaction%update_permeability .or. &
        this%realization%reaction%update_mineral_surface_area) then
      call RealizationUpdatePropertiesTS(this%realization)
      update_initial_porosity = PETSC_FALSE
    endif
  endif
  if (update_initial_porosity) then
    call this%comm1%GlobalToLocal(this%realization%field%porosity0, &
                                  this%realization%field%work_loc)
    ! push values to porosity_base
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_MINERAL)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_CURRENT)
  endif  

  call this%PreSolve()
  call this%UpdateAuxVars()
  call this%UpdateSolution() 
    
end subroutine PMSubsurfaceInitializeRun

! ************************************************************************** !

subroutine PMSubsurfaceInitializeTimestepA(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
  use Material_module
  use Material_Aux_class, only : POROSITY_MINERAL
  
  implicit none
  
  class(pm_subsurface_type) :: this

  this%option%flow_dt = this%option%dt

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

  if (this%store_porosity_for_ts_cut) then
    ! store base properties for reverting at time step cut
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,POROSITY, &
                                 POROSITY_MINERAL)
    call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                  this%realization%field%porosity_base_store)
  endif

end subroutine PMSubsurfaceInitializeTimestepA

! ************************************************************************** !

subroutine PMSubsurfaceInitializeTimestepB(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
  use Material_module
  use Material_Aux_class, only : POROSITY_CURRENT
  
  implicit none
  
  class(pm_subsurface_type) :: this

  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
    if (this%store_porosity_for_transport) then
      ! store time t properties for transport
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc,POROSITY, &
                                   POROSITY_CURRENT)
      call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                    this%realization%field%porosity_t)
    endif
  endif 
  
  ! update porosity for time t+dt
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity .or. &
        this%realization%reaction%update_tortuosity .or. &
        this%realization%reaction%update_permeability .or. &
        this%realization%reaction%update_mineral_surface_area) then
      call RealizationUpdatePropertiesTS(this%realization)
    endif
  endif
  
end subroutine PMSubsurfaceInitializeTimestepB

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
  use Material_module
  use Variables_module, only : POROSITY
  use Material_Aux_class, only : POROSITY_MINERAL
  
  implicit none
  
  class(pm_subsurface_type) :: this
  
  PetscErrorCode :: ierr
  
  this%option%flow_dt = this%option%dt
  call VecCopy(this%realization%field%flow_yy, &
               this%realization%field%flow_xx,ierr);CHKERRQ(ierr)
  if (this%store_porosity_for_transport) then
    ! store base properties for reverting at time step cut
    call this%comm1%GlobalToLocal(this%realization%field%porosity_base_store, &
                                  this%realization%field%work_loc)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,POROSITY, &
                                 POROSITY_MINERAL)
  endif             

end subroutine PMSubsurfaceTimeCut

! ************************************************************************** !

subroutine PMSubsurfaceFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14
  use Material_module
  use Global_module
  use Variables_module, only : POROSITY
  use Material_Aux_class, only : POROSITY_CURRENT

  implicit none
  
  class(pm_subsurface_type) :: this

  if (this%option%ntrandof > 0) then 
    ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
    if (this%store_porosity_for_transport) then
      ! store time t properties for transport
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc,POROSITY, &
                                   POROSITY_CURRENT)
      call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                    this%realization%field%porosity_tpdt)
    endif    
  endif
  
  call this%MaxChange()
  
end subroutine PMSubsurfaceFinalizeTimestep

! ************************************************************************** !

subroutine PMSubsurfaceUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Condition_module
  use Integral_Flux_module
  use SrcSink_Sandbox_module

  implicit none
  
  class(pm_subsurface_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE
  PetscErrorCode :: ierr

  call VecCopy(this%realization%field%flow_xx, &
               this%realization%field%flow_yy,ierr);CHKERRQ(ierr)
  
  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%realization%flow_conditions, &
                           this%realization%option, &
                           this%realization%option%time)
  call SSSandboxUpdate(ss_sandbox_list,this%realization%option%time, &
                       this%realization%option)
  ! right now, RealizUpdateAllCouplerAuxVars only updates flow
  call RealizUpdateAllCouplerAuxVars(this%realization,force_update_flag)
  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif
  call IntegralFluxUpdate(this%realization%patch%integral_flux_list, &
                          this%realization%patch%internal_flow_fluxes, &
                          this%realization%patch%boundary_flow_fluxes, &
                          INTEGRATE_FLOW,this%option)
  ! end from RealizationUpdate()

end subroutine PMSubsurfaceUpdateSolution  

! ************************************************************************** !

subroutine PMSubsurfaceUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_type) :: this

  this%option%io_buffer = 'PMSubsurfaceUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMSubsurfaceUpdateAuxVars   

! ************************************************************************** !

subroutine PMSubsurfaceCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_subsurface_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModelBinary(viewer,this%realization) 
  
end subroutine PMSubsurfaceCheckpointBinary

! ************************************************************************** !

subroutine PMSubsurfaceRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_subsurface_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModelBinary(viewer,this%realization)
  call this%UpdateAuxVars()
  call this%UpdateSolution()
  
end subroutine PMSubsurfaceRestartBinary

! ************************************************************************** !

subroutine PMSubsurfaceCheckpointHDF5(this, pm_grp_id)
  !
  ! Checkpoints data associated with Subsurface PM
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  class(pm_subsurface_type) :: this
  integer :: pm_grp_id
  print *, 'PFLOTRAN must be compiled with HDF5 to ' // &
        'write HDF5 formatted checkpoint file. Darn.'
  stop
#else

  use Checkpoint_module
  use hdf5

  implicit none

  class(pm_subsurface_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif

  call CheckpointFlowProcessModelHDF5(pm_grp_id, this%realization)

#endif

end subroutine PMSubsurfaceCheckpointHDF5

! ************************************************************************** !

subroutine PMSubsurfaceRestartHDF5(this, pm_grp_id)
  !
  ! Checkpoints data associated with Subsurface PM
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  class(pm_subsurface_type) :: this
  integer :: pm_grp_id
  print *, 'PFLOTRAN must be compiled with HDF5 to ' // &
        'write HDF5 formatted checkpoint file. Darn.'
  stop
#else

  use Checkpoint_module
  use hdf5

  implicit none

  class(pm_subsurface_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif

  call RestartFlowProcessModelHDF5(pm_grp_id, this%realization)
  call this%UpdateAuxVars()
  call this%UpdateSolution()

#endif

end subroutine PMSubsurfaceRestartHDF5

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
