module Factory_Subsurface_module

  use Simulation_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: SubsurfaceInitialize, &
            SubsurfaceInitializePostPETSc, &
            SubsurfaceJumpStart, &
            ! move to init_subsurface
            SubsurfaceReadFlowPM, &
            SubsurfaceReadRTPM, &
            SubsurfaceReadWasteFormPM

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
  use Simulation_Base_class
  use PM_Base_class
  use Creep_Closure_module
  use Klinkenberg_module
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  class(pm_base_type), pointer :: pm_list
  type(option_type), pointer :: option

  class(subsurface_simulation_type), pointer :: simulation

  ! Modules that must be initialized
  call CreepClosureInit()
  call KlinkenbergInit()
  
  ! NOTE: PETSc must already have been initialized here!
  simulation => SubsurfaceSimulationCreate(option)
  simulation%process_model_list => pm_list
  call SubsurfaceInitializePostPetsc(simulation,option)
  
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
  use PM_Waste_Form_class
  use PMC_Subsurface_class
  use PMC_Third_Party_class
  use Timestepper_BE_class
  use Realization_class
  use Logging_module
  use Simulation_Subsurface_class
  use Solver_module
  use Waypoint_module
  use Init_Common_module
  use Init_Subsurface_module
  use Input_Aux_module
  
  implicit none
  
  class(subsurface_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  class(pmc_subsurface_type), pointer :: pmc_subsurface
  class(pmc_third_party_type), pointer :: pmc_third_party
  class(pm_subsurface_type), pointer :: pm_flow
  class(pm_rt_type), pointer :: pm_rt
  class(pm_fmdm_type), pointer :: pm_waste_form
  class(pm_base_type), pointer :: cur_pm, prev_pm
  class(realization_type), pointer :: realization
  class(timestepper_BE_type), pointer :: timestepper
  character(len=MAXSTRINGLENGTH) :: string
  
  ! process command line arguments specific to subsurface
  call SubsurfInitCommandLineSettings(option)
  nullify(pm_flow)
  nullify(pm_rt)
  nullify(pm_waste_form)
  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_subsurface_type)
        pm_flow => cur_pm
      class is(pm_rt_type)
        pm_rt => cur_pm
      class is (pm_fmdm_type)
        pm_waste_form => cur_pm
      class default
        option%io_buffer = &
         'PM Class unrecognized in SubsurfaceInitializePostPetsc.'
        call printErrMsg(option)
    end select
    prev_pm => cur_pm
    cur_pm => cur_pm%next
    ! we must destroy the linkage between pms so that they are in independent
    ! lists among pmcs
    nullify(prev_pm%next)
  enddo
  call SubsurfaceSetFlowMode(pm_flow,option)
  realization => RealizationCreate(option)
  simulation%realization => realization
  realization%waypoint_list => WaypointListCreate()
  if (associated(pm_flow)) then
    pmc_subsurface => PMCSubsurfaceCreate()
    pmc_subsurface%option => option
    pmc_subsurface%pms => pm_flow
    pmc_subsurface%pm_ptr%ptr => pm_flow
    pmc_subsurface%realization => realization
    ! set up logging stage
    string = trim(pm_flow%name) // 'Flow'
    call LoggingCreateStage(string,pmc_subsurface%stage)
!    timestepper => TimestepperBECreate()
!    timestepper%solver => SolverCreate()
!    simulation%flow_process_model_coupler%timestepper => timestepper
    simulation%flow_process_model_coupler => pmc_subsurface
    simulation%process_model_coupler_list => simulation%flow_process_model_coupler
    nullify(pmc_subsurface)
  endif
  if (associated(pm_rt)) then
    pmc_subsurface => PMCSubsurfaceCreate()
    pmc_subsurface%option => option
    pmc_subsurface%pms => pm_rt
    pmc_subsurface%pm_ptr%ptr => pm_rt
    pmc_subsurface%realization => realization
    ! set up logging stage
    string = 'Reactive Transport'
    call LoggingCreateStage(string,pmc_subsurface%stage)
!    timestepper => TimestepperBECreate()
!    timestepper%solver => SolverCreate()
!    simulation%rt_process_model_coupler%timestepper => timestepper
    simulation%rt_process_model_coupler => pmc_subsurface
    if (.not.associated(simulation%process_model_coupler_list)) then
      simulation%process_model_coupler_list => pmc_subsurface
    else
      simulation%flow_process_model_coupler%child => pmc_subsurface
    endif
    nullify(pmc_subsurface)
  endif
  if (associated(pm_waste_form)) then
    if (.not.associated(simulation%rt_process_model_coupler)) then
      option%io_buffer = 'The waste form process model requires reactive ' // &
        'transport.'
      call printErrMsg(option)
    endif
    simulation%misc_process_model_coupler => PMCThirdPartyCreate()
    simulation%rt_process_model_coupler%child => &
      simulation%misc_process_model_coupler%CastToBase()
    simulation%misc_process_model_coupler%option => option
    simulation%misc_process_model_coupler%pms => pm_waste_form
    simulation%misc_process_model_coupler%pm_ptr%ptr => pm_waste_form
    simulation%misc_process_model_coupler%realization => realization
    ! set up logging stage
    string = 'Waste Form'
    call LoggingCreateStage(string,simulation%misc_process_model_coupler%stage)
  endif

  realization%input => InputCreate(IN_UNIT,option%input_filename,option)
  call InitSubsurfaceReadRequiredCards(simulation)
  call InitSubsurfaceReadInput(simulation)
  if (associated(pm_waste_form)) then
    string = 'FMDM'
    call InputFindStringInFile(realization%input,option,string)
    call InputFindStringErrorMsg(realization%input,option,string)
    call pm_waste_form%Read(realization%input)
  endif
  call InputDestroy(realization%input)
  call InitSubsurfaceSimulation(simulation)
  
  if (associated(pm_waste_form)) then
    pmc_third_party => PMCThirdPartyCreate()
    pmc_third_party%option => option
    pmc_third_party%pms => pm_waste_form
    pmc_third_party%pm_ptr%ptr => pm_waste_form
    pmc_third_party%realization => realization
    ! set up logging stage
    string = 'FMDM'
    call LoggingCreateStage(string,pmc_third_party%stage)
    simulation%rt_process_model_coupler%child => pmc_third_party
    nullify(pmc_third_party)
  endif

  ! clean up waypoints
  if (.not.option%steady_state) then
    ! fill in holes in waypoint data
    call WaypointListFillIn(option,realization%waypoint_list)
    call WaypointListRemoveExtraWaypnts(option,realization%waypoint_list)
  endif


  ! debugging output
  if (realization%debug%print_couplers) then
    call InitCommonVerifyAllCouplers(realization)
  endif
  if (realization%debug%print_waypoints) then
    call WaypointListPrint(realization%waypoint_list,option, &
                           realization%output_option)
  endif  

  call SubsurfaceJumpStart(simulation)
  ! set first process model coupler as the master
  simulation%process_model_coupler_list%is_master = PETSC_TRUE

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
      option%water_id = 1
      option%air_id = 2
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
      option%io_buffer = 'Material AuxVars must be refactored for IMMIS.'
      call printErrMsg(option)
    class is (pm_miscible_type)
      option%iflowmode = MIS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 2
      option%nflowspec = 2
      option%io_buffer = 'Material AuxVars must be refactored for MISCIBLE.'
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
      option%water_id = 1
      option%air_id = 2
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

  use General_module

  implicit none
  
  type(input_type) :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SIMULATION,PROCESS_MODEL'

  nullify(pm)
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
            call InputKeywordUnrecognized(word, &
                     'SIMULATION,PROCESS_MODELS,SUBSURFACE_FLOW,MODE',option)
        end select
        pm%option => option
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'MODE keyword must be read first under ' // &
            'SUBSURFACE_FLOW PROCESS_MODEL.'
          call printErrMsg(option)
        endif
        select type(pm)
          class is(pm_general_type)
            ! inorder to not immediately return out of GeneralRead
            !TODO(geh): remove dummy word
            input%buf = 'dummy_word'
            call pm%Read(input)
          class is(pm_th_type)
            call pm%Read(input)
          class default
            option%io_buffer = 'OPTIONS not set up for PM.'
            call printErrMsg(option)
        end select
      case default
        error_string = trim(error_string) // ',SUBSURFACE_FLOW'
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
  if (.not.associated(pm)) then
    option%io_buffer = 'A flow MODE (card) must be included in the ' // &
      'SUBSURFACE_FLOW block in ' // trim(error_string) // '.'
    call printErrMsg(option)
  endif
  
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
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SIMULATION,PROCESS_MODEL'

  pm => PMRTCreate()
  pm%option => option
  
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

subroutine SubsurfaceReadWasteFormPM(input, option, pm)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module
  
  use PM_Base_class
  use PM_Waste_Form_class

  implicit none
  
  type(input_type) :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SIMULATION,PROCESS_MODEL'

  pm => PMWasteFormCreate()
  pm%option => option
  
  word = ''
  do   
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case default
    end select
  enddo
  
end subroutine SubsurfaceReadWasteFormPM

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
  use PM_Base_Pointer_module
  use PM_Subsurface_class
  
  use PM_General_class
  use PM_Richards_class
  use PM_TH_class
  use PM_RT_class
  use PM_Waste_Form_class

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
  if (associated(simulation%flow_process_model_coupler)) then
    if (associated(simulation%flow_process_model_coupler%timestepper)) then
      simulation%flow_process_model_coupler%timestepper%cur_waypoint => &
        realization%waypoint_list%first
    endif
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    if (associated(simulation%rt_process_model_coupler%timestepper)) then
      simulation%rt_process_model_coupler%timestepper%cur_waypoint => &
        realization%waypoint_list%first
    endif
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

  if (option%nflowdof > 0) then
    select type(ts => simulation%flow_process_model_coupler%timestepper)
      class is (timestepper_BE_type)
        call InitSubsurfFlowSetupSolvers(realization,ts%convergence_context, &
                                         ts%solver)
    end select
  endif
  if (option%ntrandof > 0) then
    select type(ts => simulation%rt_process_model_coupler%timestepper)
      class is (timestepper_BE_type)
        call InitSubsurfTranSetupSolvers(realization,ts%convergence_context, &
                                         ts%solver)
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
          class is (pm_fmdm_type)
            call cur_process_model%PMwasteFormSetRealization(realization)
        end select
        ! set time stepper
        select type(cur_process_model)
          class is (pm_subsurface_type)
            cur_process_model_coupler%timestepper%dt = option%flow_dt
          class is (pm_rt_type)
            cur_process_model_coupler%timestepper%dt = option%tran_dt
        end select
        cur_process_model%output_option => realization%output_option
        call cur_process_model%Init()
        if (associated(cur_process_model_coupler%timestepper)) then
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
          select type(ts => cur_process_model_coupler%timestepper)
            class is(timestepper_BE_type)
              call SNESGetLineSearch(ts%solver%snes,linesearch,ierr);CHKERRQ(ierr)
              ! Post
              select type(cur_process_model)
                ! flow solutions
                class is(pm_subsurface_type)
                  if (ts%solver%check_post_convergence) then
                    call SNESLineSearchSetPostCheck(linesearch, &
                                                    PMCheckUpdatePostPtr, &
                                               cur_process_model_coupler%pm_ptr, &
                                                    ierr);CHKERRQ(ierr)
                  endif
                class is(pm_rt_type)
                  if (ts%solver%check_post_convergence .or. option%use_mc) then
                    call SNESLineSearchSetPostCheck(linesearch, &
                                                    PMCheckUpdatePostPtr, &
                                               cur_process_model_coupler%pm_ptr, &
                                                    ierr);CHKERRQ(ierr)
                  endif
              end select
              ! Pre
              select type(cur_process_model)
                class is(pm_richards_type)
                  if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
                      dabs(option%saturation_change_limit) > 0.d0) then
                    call SNESLineSearchSetPreCheck(linesearch, &
                                                   PMCheckUpdatePrePtr, &
                                               cur_process_model_coupler%pm_ptr, &
                                                   ierr);CHKERRQ(ierr)
                  endif              
                class is(pm_general_type)
                  call SNESLineSearchSetPreCheck(linesearch, &
                                                 PMCheckUpdatePrePtr, &
                                               cur_process_model_coupler%pm_ptr, &
                                                 ierr);CHKERRQ(ierr)
                class is(pm_th_type)
                  if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
                      dabs(option%pressure_change_limit) > 0.d0 .or. &
                      dabs(option%temperature_change_limit) > 0.d0) then
                    call SNESLineSearchSetPreCheck(linesearch, &
                                                   PMCheckUpdatePrePtr, &
                                               cur_process_model_coupler%pm_ptr, &
                                                   ierr);CHKERRQ(ierr)
                  endif 
                class is(pm_rt_type)
                  if (realization%reaction%check_update) then
                    call SNESLineSearchSetPreCheck(linesearch, &
                                                   PMCheckUpdatePrePtr, &
                                               cur_process_model_coupler%pm_ptr, &
                                                   ierr);CHKERRQ(ierr)
                  endif
                class default
              end select
          end select
          select type(cur_process_model)
            class default
              select type(ts => cur_process_model_coupler%timestepper)
                class is(timestepper_BE_type)
                  call SNESSetFunction(ts%solver%snes, &
                                       cur_process_model%residual_vec, &
                                       PMResidual, &
                                       cur_process_model, &
                                       ierr);CHKERRQ(ierr)
                  call SNESSetJacobian(ts%solver%snes, &
                                       ts%solver%J, &
                                       ts%solver%Jpre, &
                                       PMJacobian, &
                                       cur_process_model, &
                                       ierr);CHKERRQ(ierr)
              end select
          end select
        endif ! if associated(cur_process_model_coupler%timestepper)
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

#if 0
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
#endif  

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

end subroutine SubsurfaceJumpStart

end module Factory_Subsurface_module
