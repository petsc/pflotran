module Factory_Subsurface_module

  use Simulation_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: SubsurfaceInitialize, &
            SubsurfaceInitializePostPETSc, &
            HijackSimulation, &
            SubsurfaceJumpStart, &
            ! move to init_subsurface
            SubsurfaceReadFlowPM, &
            SubsurfaceReadRTPM

contains

! ************************************************************************** !

subroutine SubsurfaceInitialize(simulation_base,pm_list,option)
  ! 
  ! Sets up PFLOTRAN subsurface simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  use Option_module
  use Input_Aux_module
  use Simulation_Base_class
  use PM_Base_class
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  class(pm_base_type), pointer :: pm_list
  type(option_type), pointer :: option

  class(subsurface_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => SubsurfaceSimulationCreate(option)
  simulation%process_model_list => pm_list
  call SubsurfaceInitializePostPetsc(simulation,option)
  
  ! set first process model coupler as the master
  simulation%process_model_coupler_list%is_master = PETSC_TRUE
  
  simulation_base => simulation

end subroutine SubsurfaceInitialize

! ************************************************************************** !

subroutine SubsurfaceInitializePostPetsc(simulation, option)
  ! 
  ! Sets up PFLOTRAN subsurface simulation
  ! framework after to PETSc initialization
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  ! 

  use Option_module
  use PM_Subsurface_class
  use PM_Base_class
  use PM_RT_class
  use Timestepper_BE_class
  use Realization_class
  use Logging_module
#ifndef INIT_REFACTOR
  use Simulation_module
  use Init_Common_module
#else  
  use Simulation_Subsurface_class
  use PMC_Subsurface_class
  use Solver_module
  use Waypoint_module
  use Init_Subsurface_module
  use Input_Aux_module
#endif  
  
  implicit none
  
  class(subsurface_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  class(pm_subsurface_type), pointer :: pm_flow
  class(pm_rt_type), pointer :: pm_rt
  class(pm_base_type), pointer :: cur_pm
  class(realization_type), pointer :: realization
  class(timestepper_BE_type), pointer :: timestepper
  character(len=MAXSTRINGLENGTH) :: string
  
#ifndef INIT_REFACTOR
  type(simulation_type), pointer :: simulation_old
  
  ! process command line arguments specific to subsurface
  call SubsurfInitCommandLineSettings(option)

  simulation_old => SimulationCreate(option)
  call Init(simulation_old)
  call HijackSimulation(simulation_old,simulation)
  
  ! no longer need simulation
  ! nullify realization and regression so that it is not destroyed
  nullify(simulation_old%realization)
  nullify(simulation_old%regression)
  call SimulationDestroy(simulation_old)
#else
  nullify(pm_flow)
  nullify(pm_rt)
  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_subsurface_type)
        pm_flow => cur_pm
      class is(pm_rt_type)
        pm_rt => cur_pm
      class default
        option%io_buffer = &
         'PM Class unrecogmized in SubsurfaceInitializePostPetsc.'
        call printErrMsg(option)
    end select
    cur_pm => cur_pm%next
  enddo
  call SubsurfaceSetFlowMode(pm_flow,option)
  realization => RealizationCreate(option)
  simulation%realization => realization
  realization%waypoint_list => WaypointListCreate()
  if (associated(pm_flow)) then
    simulation%flow_process_model_coupler => PMCSubsurfaceCreate()
    simulation%process_model_coupler_list => simulation%flow_process_model_coupler
    simulation%flow_process_model_coupler%option => option
    simulation%flow_process_model_coupler%pms => pm_flow
    simulation%flow_process_model_coupler%pm_ptr%ptr => pm_flow
    simulation%flow_process_model_coupler%realization => realization
    ! set up logging stage
    string = trim(pm_flow%name) // 'Flow'
    call LoggingCreateStage(string,simulation%flow_process_model_coupler%stage)
!    timestepper => TimestepperBECreate()
!    timestepper%solver => SolverCreate()
!    simulation%flow_process_model_coupler%timestepper => timestepper
  endif
  if (associated(pm_rt)) then
    simulation%rt_process_model_coupler => PMCSubsurfaceCreate()
    if (.not.associated(simulation%process_model_coupler_list)) then
      simulation%process_model_coupler_list => simulation%rt_process_model_coupler
    endif
    simulation%rt_process_model_coupler%option => option
    simulation%rt_process_model_coupler%pms => pm_rt
    simulation%rt_process_model_coupler%pm_ptr%ptr => pm_rt
    simulation%rt_process_model_coupler%realization => realization
    ! set up logging stage
    string = 'Reactive Transport'
    call LoggingCreateStage(string,simulation%rt_process_model_coupler%stage)
!    timestepper => TimestepperBECreate()
!    timestepper%solver => SolverCreate()
!    simulation%rt_process_model_coupler%timestepper => timestepper
  endif
  
  realization%input => InputCreate(IN_UNIT,option%input_filename,option)
  call InitSubsurfaceReadRequiredCards(realization)
  call InitSubsurfaceReadInput(simulation)
  call InputDestroy(realization%input)
  call InitSubsurfaceSimulation(simulation)
  
#endif
  call SubsurfaceJumpStart(simulation)

end subroutine SubsurfaceInitializePostPetsc

! ************************************************************************** !

subroutine SubsurfInitCommandLineSettings(option)
  ! 
  ! Initializes PFLTORAN subsurface output
  ! filenames, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  ! 

  use Option_module
  use Input_Aux_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = MULTISIMULATION_SIM_TYPE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = STOCHASTIC_SIM_TYPE
  endif
  
end subroutine SubsurfInitCommandLineSettings

! ************************************************************************** !

subroutine SubsurfaceSetFlowMode(pm_flow,option)
  ! 
  ! Sets the flow mode (richards, vadose, mph, etc.)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  use Option_module
  use PM_Subsurface_class
  use PM_Base_class
  use PM_Flash2_class
  use PM_General_class
  use PM_Immis_class
  use PM_Miscible_class
  use PM_Mphase_class
  use PM_Richards_class
  use PM_TH_class

  implicit none 

  type(option_type) :: option
  class(pm_subsurface_type), pointer :: pm_flow
  
  if (.not.associated(pm_flow)) then
    option%nphase = 1
    option%liquid_phase = 1
    ! assume default isothermal when only transport
    option%use_isothermal = PETSC_TRUE  
    return
  endif
  
  select type(pm_flow)
    class is (pm_flash2_type)
      option%iflowmode = FLASH2_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%use_isothermal = PETSC_FALSE
    class is (pm_general_type)
      option%iflowmode = G_MODE
      option%nphase = 2
      option%liquid_phase = 1  ! liquid_pressure
      option%gas_phase = 2     ! gas_pressure

      option%air_pressure_id = 3
      option%capillary_pressure_id = 4
      option%vapor_pressure_id = 5
      option%saturation_pressure_id = 6

      option%water_id = 1
      option%air_id = 2
      option%energy_id = 3

      option%nflowdof = 3
      option%nflowspec = 2
      option%use_isothermal = PETSC_FALSE
    class is (pm_immis_type)
      option%iflowmode = IMS_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%io_buffer = 'Material Auxvars must be refactored for IMMIS.'
      call printErrMsg(option)
    class is (pm_miscible_type)
      option%iflowmode = MIS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 2
      option%nflowspec = 2
      option%io_buffer = 'Material Auxvars must be refactored for MISCIBLE.'
      call printErrMsg(option)
    class is (pm_mphase_type)
      option%iflowmode = MPH_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2 ! read CO2DATA0.dat
!     option%itable = 1 ! create CO2 database: co2data.dat
      option%use_isothermal = PETSC_FALSE
    class is (pm_richards_type)
      option%iflowmode = RICHARDS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    class is (pm_th_type)
      option%iflowmode = TH_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2
      option%nflowdof = 2
      option%nflowspec = 1
      option%use_isothermal = PETSC_FALSE
      option%flow%store_fluxes = PETSC_TRUE
    class default
  end select
  
end subroutine SubsurfaceSetFlowMode

! ************************************************************************** !

subroutine SubsurfaceReadFlowPM(input, option, pm)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module
  
  use PMC_Base_class
  
  use PM_Base_class
  use PM_Flash2_class
  use PM_General_class
  use PM_Immis_class
  use PM_Miscible_class
  use PM_Mphase_class
  use PM_Richards_class
  use PM_TH_class
  
  use Init_Common_module

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  class(pm_base_type), pointer :: pm
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SIMULATION,PROCESS_MODEL'

  word = ''
  do   
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('MODE')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'flow mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('GENERAL')
            pm => PMGeneralCreate()
          case('MPHASE')
            pm => PMMphaseCreate()
          case('FLASH2')
            pm => PMFlash2Create()
          case('IMS','IMMIS','THS')
            pm => PMImmisCreate()
          case('MIS','MISCIBLE')
            pm => PMMiscibleCreate()
          case('RICHARDS')
            pm => PMRichardsCreate()
          case('TH')
            pm => PMTHCreate()
          case default
            option%io_buffer = 'FLOW PM "' // trim(word) // '" not recognized.'
            call printErrMsg(option)
        end select
      case default
    end select
  enddo
  
end subroutine SubsurfaceReadFlowPM

! ************************************************************************** !

subroutine SubsurfaceReadRTPM(input, option, pm)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module
  
  use PMC_Base_class
  use PM_Base_class
  use PM_RT_class
  
  use Init_Common_module

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  class(pm_base_type), pointer :: pm
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SIMULATION,PROCESS_MODEL'

  pm => PMRTCreate()
  
  word = ''
  do   
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('OPERATOR_SPLIT','OPERATOR_SPLITTING')
      case default ! includes 'GLOBAL_IMPLICIT'
    end select
  enddo
  
end subroutine SubsurfaceReadRTPM


! ************************************************************************** !

subroutine InitSubsurfaceSimulation(simulation)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Realization_class
  use Option_module
  use Output_module, only : Output
  use Output_Aux_module
  
  
  use Global_module
  use Init_Common_module
  use Init_Subsurface_module
  use Init_Subsurface_Flow_module
  use Init_Subsurface_Tran_module
  use Waypoint_module
  use Strata_module
  use Regression_module
  
  use PMC_Subsurface_class
  use PMC_Material_class
  use PMC_Base_class
  use PM_Base_class
  use PM_Subsurface_class
  use PM_RT_class
  use Timestepper_BE_class
  
  implicit none
  
#include "finclude/petscsnes.h"  
  
  class(subsurface_simulation_type) :: simulation
  
  class(pmc_subsurface_type), pointer :: flow_process_model_coupler
  class(pmc_subsurface_type), pointer :: tran_process_model_coupler
  class(pmc_material_type), pointer :: material_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pm_base_type), pointer :: cur_process_model
  
  class(realization_type), pointer :: realization
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  SNESLineSearch :: linesearch
  PetscErrorCode :: ierr
  
  realization => simulation%realization
  option => realization%option

! begin from old Init()  
  call InitSubsurfSetupRealization(realization)
  
  !TODO(geh): refactor
  if (associated(simulation%flow_process_model_coupler%timestepper)) then
    simulation%flow_process_model_coupler%timestepper%cur_waypoint => &
      realization%waypoint_list%first
  endif
  if (associated(simulation%rt_process_model_coupler%timestepper)) then
    simulation%rt_process_model_coupler%timestepper%cur_waypoint => &
      realization%waypoint_list%first
  endif
  
  !TODO(geh): refactor
  ! initialize global auxiliary variable object
  call GlobalSetup(realization)
  
  ! always call the flow side since a velocity field still has to be
  ! set if no flow exists
  call InitSubsurfFlowSetupRealization(realization)
  if (option%ntrandof > 0) call InitSubsurfTranSetupRealization(realization)
  call OutputVariableAppendDefaults(realization%output_option% &
                                      output_variable_list,option)
    ! check for non-initialized data sets, e.g. porosity, permeability
  call RealizationNonInitializedData(realization)

  if (realization%debug%print_couplers) then
    call InitCommonVerifyAllCouplers(realization)
  endif
  if (realization%debug%print_waypoints) then
    call WaypointListPrint(realization%waypoint_list,option, &
                           realization%output_option)
  endif  

  if (option%nflowdof > 0) then
    select type(ts => simulation%flow_process_model_coupler%timestepper)
      class is (timestepper_BE_type)
        call InitSubsurfFlowSetupSolvers(realization,ts%solver)
    end select
  endif
  if (option%ntrandof > 0) then
    select type(ts => simulation%rt_process_model_coupler%timestepper)
      class is (timestepper_BE_type)
        call InitSubsurfFlowSetupSolvers(realization,ts%solver)
    end select
  endif
  call RegressionCreateMapping(simulation%regression,realization)
! end from old Init()
  
  simulation%waypoint_list => RealizCreateSyncWaypointList(realization)

  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!
  simulation%output_option => realization%output_option

  
  if (StrataEvolves(realization%patch%strata_list)) then
    material_process_model_coupler => PMCMaterialCreate()
    material_process_model_coupler%option => option
    material_process_model_coupler%realization => realization
    ! place the material process model as %peer for the top pmc
    simulation%process_model_coupler_list%peer => &
      material_process_model_coupler
  endif

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! - SNES functions (Residual/Jacobain), or TS function (RHSFunction)
  cur_process_model_coupler_top => simulation%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler_top%waypoint_list => realization%waypoint_list
    cur_process_model_coupler => cur_process_model_coupler_top
    do
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%pms
      do
        if (.not.associated(cur_process_model)) exit
        ! set realization
        select type(cur_process_model)
          class is (pm_subsurface_type)
            call cur_process_model%PMSubsurfaceSetRealization(realization)
          class is (pm_rt_type)
            call cur_process_model%PMRTSetRealization(realization)
        end select
        ! set time stepper
        select type(cur_process_model)
          class is (pm_rt_type)
            cur_process_model_coupler%timestepper%dt = option%tran_dt
          class default ! otherwise flow
            cur_process_model_coupler%timestepper%dt = option%flow_dt
        end select
        cur_process_model%output_option => realization%output_option
        call cur_process_model%Init()
        ! Until classes are resolved as user-defined contexts in PETSc, 
        ! we cannot use SetupSolvers.  Therefore, everything has to be
        ! explicitly defined here.  This may be easier in the long
        ! run as it creates an intermediate refactor in pulling 
        ! functionality in Init() into the factories - geh
#if 0        
        select type(ts => cur_process_model_coupler%timestepper)
          class is(timestepper_BE_type)
            call cur_process_model%SetupSolvers(ts%solver)
        end select
#endif
#if 0
        select type(ts => cur_process_model_coupler%timestepper)
          class is(timestepper_BE_type)
            call SNESGetLineSearch(ts%solver%snes, linesearch, ierr);CHKERRQ(ierr)
            select type(cur_process_model)
              class is(pm_richards_type)
                if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
                    dabs(option%saturation_change_limit) > 0.d0) then
                  call SNESLineSearchSetPreCheck(linesearch, &
                                                 RichardsCheckUpdatePre, &
                                                 realization,ierr);CHKERRQ(ierr)
                endif              
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch, &
                                                  RichardsCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)        
                endif
              class is(pm_general_type)
                call SNESLineSearchSetPreCheck(linesearch, &
                                               GeneralCheckUpdatePre, &
                                               realization,ierr);CHKERRQ(ierr)              
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch, &
                                                  GeneralCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)        
                endif
              class is(pm_th_type)
                if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
                    dabs(option%pressure_change_limit) > 0.d0 .or. &
                    dabs(option%temperature_change_limit) > 0.d0) then
                  call SNESLineSearchSetPreCheck(linesearch, &
                                                 THCheckUpdatePre, &
                                                 realization,ierr);CHKERRQ(ierr)
                endif 
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch, &
                                                  THCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)        
                endif
              class is(pm_rt_type)
                if (realization%reaction%check_update) then
                  call SNESLineSearchSetPreCheck(linesearch,RTCheckUpdatePre, &
                                                 realization,ierr);CHKERRQ(ierr)
                endif
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch,RTCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)
                endif        
              class default
            end select
        end select
#endif
#if 0        
        select type(cur_process_model)
          class default
            select type(ts => cur_process_model_coupler%timestepper)
              class is(timestepper_BE_type)
                call SNESSetFunction(ts%solver%snes, &
                                     cur_process_model%residual_vec, &
                                     PMResidual, &
                                     cur_process_model_coupler%pm_ptr, &
                                     ierr);CHKERRQ(ierr)
                call SNESSetJacobian(ts%solver%snes, &
                                     ts%solver%J, &
                                     ts%solver%Jpre, &
                                     PMJacobian, &
                                     cur_process_model_coupler%pm_ptr, &
                                     ierr);CHKERRQ(ierr)
            end select
        end select
#endif            
        cur_process_model => cur_process_model%next
      enddo
      ! has to be called after realizations are set above
      call cur_process_model_coupler%SetupSolvers()
      cur_process_model_coupler => cur_process_model_coupler%child
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%peer
  enddo
  
  ! point the top process model coupler to Output
  simulation%process_model_coupler_list%Output => Output

end subroutine InitSubsurfaceSimulation

! ************************************************************************** !

subroutine HijackSimulation(simulation_old,simulation)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Simulation_module
  use Realization_class
  use Option_module
  use Output_module, only : Output
  use Output_Aux_module
  
  use PMC_Base_class
  use PMC_Subsurface_class  
  use PMC_Material_class
  use Simulation_Base_class
  use PM_Base_class
  use PM_General_class
  use PM_Flash2_class
  use PM_Immis_class
  use PM_Mphase_class
  use PM_Miscible_class
  use PM_Richards_class
  use PM_RT_class
  use PM_Subsurface_class
  use PM_TH_class
  use PM_Base_Pointer_module
  use Timestepper_BE_class
  use Logging_module
  use Strata_module
  
  use General_module
  use TH_module
  use Richards_module
  use Reactive_Transport_module
  
  use Global_module
  use Init_Common_module
  use Init_Subsurface_module
  use Init_Subsurface_Flow_module
  use Init_Subsurface_Tran_module
  use Waypoint_module
  use Regression_module
  
  implicit none
  
#include "finclude/petscsnes.h"  
  
  type(simulation_type) :: simulation_old
  class(subsurface_simulation_type) :: simulation
  
  class(pmc_subsurface_type), pointer :: flow_process_model_coupler
  class(pmc_subsurface_type), pointer :: tran_process_model_coupler
  class(pmc_material_type), pointer :: material_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pm_base_type), pointer :: cur_process_model
  
  class(realization_type), pointer :: realization
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  SNESLineSearch :: linesearch
  PetscErrorCode :: ierr
  
  realization => simulation_old%realization
  option => realization%option

! begin from old Init()  
  call InitSubsurfSetupRealization(realization)
  
  !TODO(geh): refactor
  if (associated(simulation_old%flow_timestepper)) then
    simulation_old%flow_timestepper%cur_waypoint => &
      realization%waypoint_list%first
  endif
  if (associated(simulation_old%tran_timestepper)) then
    simulation_old%tran_timestepper%cur_waypoint => &
      realization%waypoint_list%first
  endif
  
  !TODO(geh): refactor
  ! initialize global auxiliary variable object
  call GlobalSetup(realization)
  
  ! always call the flow side since a velocity field still has to be
  ! set if no flow exists
  call InitSubsurfFlowSetupRealization(realization)
  if (option%ntrandof > 0) call InitSubsurfTranSetupRealization(realization)
  call OutputVariableAppendDefaults(realization%output_option% &
                                      output_variable_list,option)
    ! check for non-initialized data sets, e.g. porosity, permeability
  call RealizationNonInitializedData(realization)

  if (realization%debug%print_couplers) then
    call InitCommonVerifyAllCouplers(realization)
  endif
  if (realization%debug%print_waypoints) then
    call WaypointListPrint(realization%waypoint_list,option, &
                           realization%output_option)
  endif  

  if (option%nflowdof > 0) call InitSubsurfFlowSetupSolvers(realization, &
                                         simulation_old%flow_timestepper%solver)
  if (option%ntrandof > 0) call InitSubsurfTranSetupSolvers(realization, &
                                         simulation_old%tran_timestepper%solver)
  call RegressionCreateMapping(simulation_old%regression,realization)
! end from old Init()
  
  simulation%waypoint_list => RealizCreateSyncWaypointList(realization)

  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!
  simulation%output_option => realization%output_option
  simulation%option => realization%option
  nullify(cur_process_model)

  nullify(flow_process_model_coupler)
  nullify(tran_process_model_coupler)
  
  ! Create Subsurface-flow ProcessModel & ProcessModelCoupler
  if (option%nflowdof > 0) then
    select case(option%iflowmode)
      case(G_MODE)
        cur_process_model => PMGeneralCreate()
      case(FLASH2_MODE)
        cur_process_model => PMFlash2Create()
      case(IMS_MODE)
        cur_process_model => PMImmisCreate()
      case(MPH_MODE)
        cur_process_model => PMMphaseCreate()
      case(MIS_MODE)
        cur_process_model => PMMiscibleCreate()
      case(RICHARDS_MODE)
        cur_process_model => PMRichardsCreate()
      case(TH_MODE)
        cur_process_model => PMTHCreate()
    end select
    cur_process_model%option => realization%option
    cur_process_model%output_option => realization%output_option

    flow_process_model_coupler => PMCSubsurfaceCreate()
    flow_process_model_coupler%option => option
    flow_process_model_coupler%pms => cur_process_model
    flow_process_model_coupler%pm_ptr%ptr => cur_process_model
!    flow_process_model_coupler%timestepper => simulation_old%flow_timestepper
    flow_process_model_coupler%realization => realization
    call HijackTimestepper(simulation_old%flow_timestepper, &
                           flow_process_model_coupler%timestepper)
    flow_process_model_coupler%timestepper%name = 'FLOW'
    ! set up logging stage
    string = trim(cur_process_model%name) // 'Flow'
    call LoggingCreateStage(string,flow_process_model_coupler%stage)
    nullify(cur_process_model)
  endif

  ! Create Subsurface transport ProcessModel & ProcessModelCoupler
  if (option%ntrandof > 0) then
    cur_process_model => PMRTCreate()
    cur_process_model%output_option => realization%output_option
    cur_process_model%option => realization%option
   
    tran_process_model_coupler => PMCSubsurfaceCreate()
    tran_process_model_coupler%option => option
    tran_process_model_coupler%pms => cur_process_model
    tran_process_model_coupler%pm_ptr%ptr => cur_process_model
!    tran_process_model_coupler%timestepper => simulation_old%tran_timestepper
    tran_process_model_coupler%realization => realization
    call HijackTimestepper(simulation_old%tran_timestepper, &
                           tran_process_model_coupler%timestepper)
    tran_process_model_coupler%timestepper%name = 'TRAN'
    ! set up logging stage
    string = 'Reactive Transport'
    call LoggingCreateStage(string,tran_process_model_coupler%stage)
    nullify(cur_process_model)
  endif

  ! Add the ProcessModelCouplers in a list
  if (associated(flow_process_model_coupler)) then
    simulation%process_model_coupler_list => &
      flow_process_model_coupler%CastToBase()
    if (associated(tran_process_model_coupler)) then
      flow_process_model_coupler%child => &
        tran_process_model_coupler%CastToBase()
    endif
  else
    simulation%process_model_coupler_list => &
      tran_process_model_coupler%CastToBase()
  endif
  
  if (StrataEvolves(realization%patch%strata_list)) then
    material_process_model_coupler => PMCMaterialCreate()
    material_process_model_coupler%option => option
    material_process_model_coupler%realization => realization
    ! place the material process model as %peer for the top pmc
    simulation%process_model_coupler_list%peer => &
      material_process_model_coupler
  endif

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! - SNES functions (Residual/Jacobain), or TS function (RHSFunction)
  cur_process_model_coupler_top => simulation%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler_top%waypoint_list => realization%waypoint_list
    cur_process_model_coupler => cur_process_model_coupler_top
    do
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%pms
      do
        if (.not.associated(cur_process_model)) exit
        ! set realization
        select type(cur_process_model)
          class is (pm_subsurface_type)
            call cur_process_model%PMSubsurfaceSetRealization(realization)
          class is (pm_rt_type)
            call cur_process_model%PMRTSetRealization(realization)
        end select
        ! set time stepper
        select type(cur_process_model)
          class is (pm_rt_type)
            call cur_process_model_coupler%SetTimestepper( &
                   tran_process_model_coupler%timestepper)
            tran_process_model_coupler%timestepper%dt = option%tran_dt
          class default ! otherwise flow
            call cur_process_model_coupler%SetTimestepper( &
                   flow_process_model_coupler%timestepper)
            flow_process_model_coupler%timestepper%dt = option%flow_dt
        end select

        call cur_process_model%Init()
        ! Until classes are resolved as user-defined contexts in PETSc, 
        ! we cannot use SetupSolvers.  Therefore, everything has to be
        ! explicitly defined here.  This may be easier in the long
        ! run as it creates an intermediate refactor in pulling 
        ! functionality in Init() into the factories - geh
#if 0        
        select type(ts => cur_process_model_coupler%timestepper)
          class is(timestepper_BE_type)
            call cur_process_model%SetupSolvers(ts%solver)
        end select
#endif
#if 0
        select type(ts => cur_process_model_coupler%timestepper)
          class is(timestepper_BE_type)
            call SNESGetLineSearch(ts%solver%snes, linesearch, ierr);CHKERRQ(ierr)
            select type(cur_process_model)
              class is(pm_richards_type)
                if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
                    dabs(option%saturation_change_limit) > 0.d0) then
                  call SNESLineSearchSetPreCheck(linesearch, &
                                                 RichardsCheckUpdatePre, &
                                                 realization,ierr);CHKERRQ(ierr)
                endif              
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch, &
                                                  RichardsCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)        
                endif
              class is(pm_general_type)
                call SNESLineSearchSetPreCheck(linesearch, &
                                               GeneralCheckUpdatePre, &
                                               realization,ierr);CHKERRQ(ierr)              
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch, &
                                                  GeneralCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)        
                endif
              class is(pm_th_type)
                if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
                    dabs(option%pressure_change_limit) > 0.d0 .or. &
                    dabs(option%temperature_change_limit) > 0.d0) then
                  call SNESLineSearchSetPreCheck(linesearch, &
                                                 THCheckUpdatePre, &
                                                 realization,ierr);CHKERRQ(ierr)
                endif 
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch, &
                                                  THCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)        
                endif
              class is(pm_rt_type)
                if (realization%reaction%check_update) then
                  call SNESLineSearchSetPreCheck(linesearch,RTCheckUpdatePre, &
                                                 realization,ierr);CHKERRQ(ierr)
                endif
                if (ts%solver%check_post_convergence) then
                  call SNESLineSearchSetPostCheck(linesearch,RTCheckUpdatePost, &
                                                  realization,ierr);CHKERRQ(ierr)
                endif        
              class default
            end select
        end select
#endif
#if 0        
        select type(cur_process_model)
          class default
            select type(ts => cur_process_model_coupler%timestepper)
              class is(timestepper_BE_type)
                call SNESSetFunction(ts%solver%snes, &
                                     cur_process_model%residual_vec, &
                                     PMResidual, &
                                     cur_process_model_coupler%pm_ptr, &
                                     ierr);CHKERRQ(ierr)
                call SNESSetJacobian(ts%solver%snes, &
                                     ts%solver%J, &
                                     ts%solver%Jpre, &
                                     PMJacobian, &
                                     cur_process_model_coupler%pm_ptr, &
                                     ierr);CHKERRQ(ierr)
            end select
        end select
#endif            
        cur_process_model => cur_process_model%next
      enddo
      ! has to be called after realizations are set above
      call cur_process_model_coupler%SetupSolvers()
      cur_process_model_coupler => cur_process_model_coupler%child
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%peer
  enddo
  
  simulation%realization => realization
  simulation%flow_process_model_coupler => flow_process_model_coupler
  simulation%rt_process_model_coupler => tran_process_model_coupler
  simulation%regression => simulation_old%regression
  
  ! point the top process model coupler to Output
  simulation%process_model_coupler_list%Output => Output

end subroutine HijackSimulation

! ************************************************************************** !

subroutine SubsurfaceJumpStart(simulation)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Realization_class
  use Option_module
  use Timestepper_Base_class
  use Timestepper_BE_class
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Condition_Control_module
  use Reactive_Transport_module, only : RTJumpStartKineticSorption  

  implicit none

  type(subsurface_simulation_type) :: simulation
  
  class(realization_type), pointer :: realization
  class(timestepper_base_type), pointer :: master_timestepper
  class(timestepper_BE_type), pointer :: flow_timestepper
  class(timestepper_BE_type), pointer :: tran_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, transient_plot_flag
  PetscBool :: activity_coefs_read
  PetscBool :: flow_read
  PetscBool :: transport_read
  PetscBool :: failure
  PetscErrorCode :: ierr

  realization => simulation%realization
  
  if (associated(simulation%flow_process_model_coupler)) then
    select type(ts => simulation%flow_process_model_coupler%timestepper)
      class is(timestepper_BE_type)
        flow_timestepper => ts
    end select
  else
    nullify(flow_timestepper)
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    select type(ts => simulation%rt_process_model_coupler%timestepper)
      class is(timestepper_BE_type)
        tran_timestepper => ts
    end select
  else
    nullify(tran_timestepper)
  endif
  nullify(master_timestepper)
  
  option => realization%option
  output_option => realization%output_option

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr);CHKERRQ(ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported in refactored code.'
    call printErrMsg(option)
#if 0    
    call StepperRunSteadyState(realization,flow_timestepper,tran_timestepper)
#endif    
    ! do not want to run through time stepper
    option%status = DONE
    return 
  endif
  
  if (associated(flow_timestepper)) then
    master_timestepper => flow_timestepper
  else
    master_timestepper => tran_timestepper
  endif

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  activity_coefs_read = PETSC_FALSE
  flow_read = PETSC_FALSE
  transport_read = PETSC_FALSE
  failure = PETSC_FALSE

#if 0
!geh: moved to within PMInitialize routines
!geh: removed 8/11
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
    call CondControlAssignFlowInitCond(realization)
  endif

  if (transport_read .and. option%overwrite_restart_transport) then
    call CondControlAssignTranInitCond(realization)  
  endif

  ! turn on flag to tell RTUpdateSolution that the code is not timestepping
  if (associated(simulation%flow_process_model_coupler)) then
    call simulation%flow_process_model_coupler%UpdateSolution()
  endif
#endif
 
!geh: now performed in PMRTInitializeRun()
!  if (associated(simulation%rt_process_model_coupler)) then
!    call simulation%rt_process_model_coupler%UpdateSolution()
!  endif

  if (option%transport%jumpstart_kinetic_sorption .and. &
      option%time < 1.d-40) then
    ! only user jumpstart for a restarted simulation
    if (.not. option%restart_flag) then
      option%io_buffer = 'Only use JUMPSTART_KINETIC_SORPTION on a ' // &
        'restarted simulation.  ReactionEquilibrateConstraint() will ' // &
        'appropriately set sorbed initial concentrations for a normal ' // &
        '(non-restarted) simulation.'
      call printErrMsg(option)
    endif
    call RTJumpStartKineticSorption(realization)
  endif
#if 0  
!geh: removed 8/11
  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_timestepper%max_time_step < 0) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
    option%status = DONE
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputInit(master_timestepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_timestepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
    call Output(realization,plot_flag,transient_plot_flag)
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_timestepper%max_time_step < 1) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'') 
    option%status = DONE
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_timestepper)) then
    if (.not.associated(flow_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      flow_timestepper%dt_max = flow_timestepper%cur_waypoint%dt_max
    endif
  endif
  if (associated(tran_timestepper)) then
    if (.not.associated(tran_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start ' // &
        'time or simulation time needs to be extended on a restart.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      tran_timestepper%dt_max = tran_timestepper%cur_waypoint%dt_max
    endif
  endif
           
  if (associated(flow_timestepper)) &
    flow_timestepper%start_time_step = flow_timestepper%steps + 1
  if (associated(tran_timestepper)) &
    tran_timestepper%start_time_step = tran_timestepper%steps + 1
  
  if (realization%debug%print_couplers) then
    call OutputPrintCouplers(realization,ZERO_INTEGER)
  endif
#endif
end subroutine SubsurfaceJumpStart

! ************************************************************************** !

subroutine HijackTimestepper(timestepper_old,timestepper_base)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Timestepper_BE_class
  use Timestepper_Base_class
  use Timestepper_module
  use Solver_module

  implicit none
  
  type(timestepper_type), pointer :: timestepper_old
  class(timestepper_base_type), pointer :: timestepper_base
  
  class(timestepper_BE_type), pointer :: timestepper

  timestepper => TimestepperBECreate()
  
  timestepper%steps = timestepper_old%steps
  timestepper%num_newton_iterations = timestepper_old%num_newton_iterations
  timestepper%num_linear_iterations = timestepper_old%num_linear_iterations
  timestepper%num_constant_time_steps = timestepper_old%num_constant_time_steps

  timestepper%max_time_step = timestepper_old%max_time_step
  timestepper%max_time_step_cuts = timestepper_old%max_time_step_cuts
  timestepper%constant_time_step_threshold = timestepper_old%constant_time_step_threshold
  timestepper%iaccel = timestepper_old%iaccel

  timestepper%cumulative_newton_iterations = timestepper_old%cumulative_newton_iterations
  timestepper%cumulative_linear_iterations = timestepper_old%cumulative_linear_iterations
  timestepper%cumulative_time_step_cuts = timestepper_old%cumulative_time_step_cuts
  timestepper%cumulative_solver_time = timestepper_old%cumulative_solver_time

  timestepper%start_time = timestepper_old%start_time
  timestepper%start_time_step = timestepper_old%start_time_step
  timestepper%time_step_tolerance = timestepper_old%time_step_tolerance
  timestepper%target_time = timestepper_old%target_time
  
  timestepper%prev_dt = timestepper_old%prev_dt
!  stepper%dt = timestepper_old%dt
  timestepper%dt_init = timestepper_old%dt_init
  timestepper%dt_min = timestepper_old%dt_min
  timestepper%dt_max = timestepper_old%dt_max
  timestepper%cfl_limiter = timestepper_old%cfl_limiter
  timestepper%cfl_limiter_ts = timestepper_old%cfl_limiter_ts
  
  timestepper%time_step_cut_flag = timestepper_old%time_step_cut_flag

  timestepper%ntfac = timestepper_old%ntfac
  
  ! we destroy the new tfac here as the timestepper is pointed to the legacy
  ! one below.
  deallocate(timestepper%tfac)
  timestepper%tfac => timestepper_old%tfac
  nullify(timestepper_old%tfac)
  
  timestepper%init_to_steady_state = timestepper_old%init_to_steady_state
  timestepper%steady_state_rel_tol = timestepper_old%steady_state_rel_tol
  timestepper%run_as_steady_state = timestepper_old%run_as_steady_state

  ! we destroy the new solver here as the timestepper is pointed to the legacy
  ! one below.  Remove use Solver_module statement above when removed.
  call SolverDestroy(timestepper%solver)
  timestepper%solver => timestepper_old%solver
  nullify(timestepper_old%solver)

  timestepper%convergence_context => timestepper_old%convergence_context
  nullify(timestepper_old%convergence_context)
  timestepper%cur_waypoint => timestepper_old%cur_waypoint
  nullify(timestepper_old%cur_waypoint)
!  stepper%prev_waypoint => timestepper_old%prev_waypoint
!  nullify(timestepper_old%prev_waypoint)
  
!  stepper%revert_dt = timestepper_old%revert_dt
!  stepper%num_contig_revert_due_to_sync = &
!  timestepper_old%num_contig_revert_due_to_sync

  timestepper_base => timestepper
  
end subroutine HijackTimestepper

end module Factory_Subsurface_module
