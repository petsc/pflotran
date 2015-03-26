module Factory_Geomechanics_module

  use Simulation_Geomechanics_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: GeomechanicsInitialize

contains

! ************************************************************************** !

subroutine GeomechanicsInitialize(simulation_base,pm_list,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Option_module 
  use Simulation_Base_class
  use PM_Base_class

  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  class(pm_base_type), pointer :: pm_list
  type(option_type), pointer :: option

  class(geomechanics_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => GeomechanicsSimulationCreate(option)
  simulation%process_model_list => pm_list
  call GeomechanicsInitializePostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine GeomechanicsInitialize

! ************************************************************************** !

subroutine GeomechanicsInitializePostPETSc(simulation, option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 02/10/15
  ! 

  use Simulation_Geomechanics_class
  use Simulation_Subsurface_class
  use Factory_Subsurface_module
  use Init_Common_module
  use Geomechanics_Init_module
  use Option_module
  use PM_Base_class
  use PM_Base_Pointer_module
  use PM_Geomechanics_Force_class
  use PMC_Base_class
  use PMC_Geomechanics_class
  use PFLOTRAN_Constants_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module
  use Geomechanics_Realization_class
  use Simulation_Aux_module
  use Realization_class
  use Timestepper_Geomechanics_class
  use Input_Aux_module
  use Logging_module

  implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(geomechanics_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  class(realization_type), pointer :: subsurf_realization
  class(geomech_realization_type), pointer :: geomech_realization
  class(pmc_base_type), pointer :: cur_process_model_coupler
  type(gmdm_ptr_type), pointer  :: dm_ptr
  class(pm_base_type), pointer :: cur_pm, prev_pm
  class(pm_geomech_force_type), pointer :: pm_geomech
  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(timestepper_geomechanics_type), pointer :: timestepper
  character(len=MAXSTRINGLENGTH) :: string
  type(waypoint_type), pointer :: waypoint
  type(input_type), pointer :: input
  PetscErrorCode :: ierr

  nullify(prev_pm)
  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_geomech_force_type)
        pm_geomech => cur_pm
        if (associated(prev_pm)) then
          prev_pm%next => cur_pm%next
        else
          simulation%process_model_list => cur_pm%next
        endif
        exit
      class default
    end select
    prev_pm => cur_pm
    cur_pm => cur_pm%next
  enddo
  
  call SubsurfaceInitializePostPetsc(simulation,option)  
  simulation%process_model_coupler_list%is_master = PETSC_TRUE
    
  if (option%geomech_on) then
    simulation%geomech_realization => GeomechRealizCreate(option)
    geomech_realization => simulation%geomech_realization
    subsurf_realization => simulation%realization
    geomech_realization%input => InputCreate(IN_UNIT,option%input_filename,option)
    call GeomechicsInitReadRequiredCards(geomech_realization)
    geomech_realization%waypoint_list => WaypointListCreate()
    pm_geomech%output_option => simulation%geomech_realization%output_option
    pmc_geomech => PMCGeomechanicsCreate()
    pmc_geomech%name = 'PMCGeomech'
    simulation%geomech_process_model_coupler => pmc_geomech
    pmc_geomech%option => option
    pmc_geomech%pms => pm_geomech
    pmc_geomech%pm_ptr%ptr => pm_geomech
    pmc_geomech%geomech_realization => simulation%geomech_realization
    pmc_geomech%subsurf_realization => simulation%realization
    timestepper => TimestepperGeomechanicsCreate()
    pmc_geomech%timestepper => timestepper
    ! set up logging stage
    string = trim(pmc_geomech%name) // 'Geomechanics'
    call LoggingCreateStage(string,pmc_geomech%stage)

    input => InputCreate(IN_UNIT,option%input_filename,option)    
    string = 'GEOMECHANICS'
    call InputFindStringInFile(input,option,string)
    call InputFindStringErrorMsg(input,option,string)  
    call GeomechanicsInitReadInput(geomech_realization, &
                              timestepper%solver,input,option)

    ! Add first waypoint
    waypoint => WaypointCreate()
    waypoint%time = 0.d0
    call WaypointInsertInList(waypoint,geomech_realization%waypoint_list)
 
    ! Add final_time waypoint to geomech_realization
    waypoint => WaypointCreate()
    waypoint%final = PETSC_TRUE
    waypoint%time = simulation%realization%waypoint_list%last%time
    waypoint%print_output = PETSC_TRUE
    call WaypointInsertInList(waypoint,geomech_realization%waypoint_list)   

    if (associated(simulation%geomech_process_model_coupler)) then
      if (associated(simulation%geomech_process_model_coupler% &
                     timestepper)) then
        simulation%geomech_process_model_coupler%timestepper%cur_waypoint => &
          subsurf_realization%waypoint_list%first
      endif
    endif
 
    call InitGeomechSetupRealization(geomech_realization,subsurf_realization)
    call InitGeomechSetupSolvers(geomech_realization,subsurf_realization, &
                                  timestepper%solver)
                                  

    call pm_geomech%PMGeomechForceSetRealization(geomech_realization)
    call pm_geomech%Init()
    call SNESSetFunction(timestepper%solver%snes, &
                         pm_geomech%residual_vec, &
                         PMResidual, &
                         pmc_geomech%pm_ptr, &
                         ierr);CHKERRQ(ierr)
    call SNESSetJacobian(timestepper%solver%snes, &
                         timestepper%solver%J, &
                         timestepper%solver%Jpre, &
                         PMJacobian, &
                         pmc_geomech%pm_ptr, &
                         ierr);CHKERRQ(ierr)
                                  
     nullify(simulation%process_model_coupler_list)                             
   endif
  ! sim_aux: Create PETSc Vectors and VectorScatters
  if (option%ngeomechdof > 0) then

    call GeomechCreateGeomechSubsurfVec(subsurf_realization, &
                                        geomech_realization)
    call SimAuxCopySubsurfVec(simulation%sim_aux,subsurf_realization%field%work)

    call GeomechCreateSubsurfStressStrainVec(subsurf_realization, &
                                             geomech_realization)
    call SimAuxCopySubsurfGeomechVec(simulation%sim_aux, &
          geomech_realization%geomech_field%strain_subsurf)

    call GeomechRealizMapSubsurfGeomechGrid(subsurf_realization, &
                                            geomech_realization, &
                                            option)

    dm_ptr => GeomechDiscretizationGetDMPtrFromIndex( &
                geomech_realization%geomech_discretization, ONEDOF)

    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                              SUBSURF_TO_GEOMECHANICS)
    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                              GEOMECHANICS_TO_SUBSURF)
  endif

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if (associated(simulation%rt_process_model_coupler)) &
    simulation%rt_process_model_coupler%sim_aux => simulation%sim_aux
  if (option%ngeomechdof>0 .and. &
     associated(simulation%geomech_process_model_coupler)) &
    simulation%geomech_process_model_coupler%sim_aux => simulation%sim_aux
 
  ! set geomech as not master
  simulation%geomech_process_model_coupler%is_master = PETSC_FALSE
  ! link geomech and master
  simulation%process_model_coupler_list => &
    simulation%geomech_process_model_coupler
  ! link subsurface flow as peer
  simulation%process_model_coupler_list%peer => &
    simulation%flow_process_model_coupler  
    

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  if (associated(cur_process_model_coupler%peer)) then
    cur_process_model_coupler => cur_process_model_coupler%peer
    call cur_process_model_coupler%GetAuxData()
    call cur_process_model_coupler%SetAuxData()
    select type(pmc => cur_process_model_coupler)
      class is(pmc_geomechanics_type)
        call GeomechStoreInitialPressTemp(pmc%geomech_realization)
    end select
  endif    

  call GeomechanicsJumpStart(simulation)
  
end subroutine GeomechanicsInitializePostPETSc

! ************************************************************************** !

subroutine GeomechanicsJumpStart(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Geomechanics_Realization_class
  use Option_module
  use Timestepper_Geomechanics_class
  use Output_Aux_module
  use Output_module, only : Output, OutputPrintCouplers
  use Output_Geomechanics_module
  use Logging_module
  use Condition_Control_module

  implicit none

  type(geomechanics_simulation_type) :: simulation

  class(geomech_realization_type), pointer :: geomch_realization
  class(timestepper_geomechanics_type), pointer :: master_timestepper
  class(timestepper_geomechanics_type), pointer :: geomech_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, transient_plot_flag
  PetscBool :: geomech_read
  PetscBool :: failure
  PetscErrorCode :: ierr

  geomch_realization => simulation%geomech_realization

  select type(ts => simulation%geomech_process_model_coupler%timestepper)
    class is(timestepper_geomechanics_type)
      geomech_timestepper => ts
  end select
  nullify(master_timestepper)

  option => geomch_realization%option
  output_option => geomch_realization%output_option

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", &
                           failure, ierr);CHKERRQ(ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported for surface-flow.'
    call printErrMsg(option)
    return
  endif
  
  master_timestepper => geomech_timestepper

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  geomech_read = PETSC_FALSE
  failure = PETSC_FALSE
  
  call OutputGeomechInit(master_timestepper%steps)


end subroutine GeomechanicsJumpStart

! ************************************************************************** !

end module Factory_Geomechanics_module
