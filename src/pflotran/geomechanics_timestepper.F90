#ifdef GEOMECH

module Geomechanics_Timestepper_module
 
  use Solver_module
  use Waypoint_module 
  use Convergence_module 
  use Timestepper_module
  use PFLOTRAN_Constants_module
 
  implicit none

  private
  
#include "finclude/petscsys.h"

  public :: GeomechTimestepperInitializeRun, &
            GeomechTimestepperExecuteRun, &
            GeomechTimestepperFinalizeRun
            
contains

! ************************************************************************** !
!
! GeomechTimestepperInitializeRun: Initializes timestepping run the time step loop
! author: Glenn Hammond
! date: 03/11/13
!
! ************************************************************************** !
subroutine GeomechTimestepperInitializeRun(realization,geomech_realization, &
                                    master_stepper, &
                                    flow_stepper,tran_stepper, &
                                    geomech_stepper, &
                                    init_status)

  use Realization_class

  use Option_module
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Condition_Control_module
#ifdef GEOMECH  
  use Geomechanics_Realization_module
  use Geomechanics_Logging_module  
  use Output_Geomechanics_module
  use Geomechanics_Force_module
#endif

  implicit none

  type(realization_type) :: realization
  type(geomech_realization_type) :: geomech_realization
  type(stepper_type), pointer :: master_stepper
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(stepper_type), pointer :: geomech_stepper
  PetscInt :: init_status

  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, transient_plot_flag
  PetscBool :: activity_coefs_read
  PetscBool :: flow_read
  PetscBool :: transport_read
  PetscBool :: failure
  PetscErrorCode :: ierr
#ifdef GEOMECH  
  PetscBool :: geomech_plot_flag, geomech_transient_plot_flag
#endif

  option => realization%option
  output_option => realization%output_option

  nullify(master_stepper)
  init_status = TIMESTEPPER_INIT_PROCEED

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr)
                             
  if (option%steady_state) then
    call StepperRunSteadyState(realization,flow_stepper,tran_stepper)
    ! do not want to run through time stepper
    init_status = TIMESTEPPER_INIT_DONE
    return 
  endif
  
  if (associated(flow_stepper)) then
    master_stepper => flow_stepper
  else
    master_stepper => tran_stepper
  endif

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  activity_coefs_read = PETSC_FALSE
  flow_read = PETSC_FALSE
  transport_read = PETSC_FALSE
  failure = PETSC_FALSE
  
  if (option%restart_flag) then
    call TimestepperRestart(realization,flow_stepper,tran_stepper, &
                            flow_read,transport_read,activity_coefs_read)


  else if (master_stepper%init_to_steady_state) then
    option%print_screen_flag = OptionPrintToScreen(option)
    option%print_file_flag = OptionPrintToFile(option)
    if (associated(flow_stepper)) then
      if (flow_stepper%init_to_steady_state) then
        option%flow_dt = master_stepper%dt_min
        call FlowStepperStepToSteadyState(realization,flow_stepper,failure)
        if (failure) then ! if flow solve fails, exit
          if (OptionPrintToScreen(option)) then
            write(*,*) ' ERROR: steady state solve failed!!!'
          endif
          if (OptionPrintToFile(option)) then
            write(option%fid_out,*) ' ERROR: steady state solve failed!!!'
          endif
          init_status = TIMESTEPPER_INIT_FAIL
          return 
        endif
        option%flow_dt = master_stepper%dt_min
      endif
      if (flow_stepper%run_as_steady_state) then
        master_stepper => tran_stepper
      endif
    endif

    if (associated(tran_stepper)) then
      if (tran_stepper%init_to_steady_state) then
    ! not yet functional
    !    step_to_steady_state = PETSC_TRUE
    !    option%tran_dt = master_stepper%dt_min
    !    call StepperStepTransportDT(realization,tran_stepper, &
    !                                tran_timestep_cut_flag, &
    !                                idum,step_to_steady_state,failure)
    !    if (failure) return ! if flow solve fails, exit
    !    option%tran_dt = master_stepper%dt_min
      endif
    endif
  endif

  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif

  if (transport_read .and. option%overwrite_restart_transport) then
    call CondControlAssignTranInitCond(realization)  
  endif

  ! turn on flag to tell RTUpdateSolution that the code is not timestepping
  call StepperUpdateSolution(realization,PETSC_FALSE)

  if (option%jumpstart_kinetic_sorption .and. option%time < 1.d-40) then
    ! only user jumpstart for a restarted simulation
    if (.not. option%restart_flag) then
      option%io_buffer = 'Only use JUMPSTART_KINETIC_SORPTION on a ' // &
        'restarted simulation.  ReactionEquilibrateConstraint() will ' // &
        'appropriately set sorbed initial concentrations for a normal ' // &
        '(non-restarted) simulation.'
      call printErrMsg(option)
    endif
    call StepperJumpStart(realization)
  endif
  
  ! pushed in Init()
  call PetscLogStagePop(ierr)
  option%init_stage = PETSC_FALSE

  ! popped in TimeStepperFinalizeRun()
  call PetscLogStagePush(geomech_logging%stage(GEOMECH_TS_STAGE),ierr)

  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_stepper%max_time_step < 0) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
    init_status = TIMESTEPPER_INIT_DONE
    return
  endif
  
#ifdef GEOMECH       
    if (option%ngeomechdof > 0) then
      call GeomechUpdateFromSubsurf(realization,geomech_realization)
      call GeomechStoreInitialPressTemp(geomech_realization)
      call StepperSolveGeomechSteadyState(geomech_realization,geomech_stepper, &
                                          failure)
      call GeomechUpdateSolution(geomech_realization)
      call GeomechStoreInitialDisp(geomech_realization)
      call GeomechForceUpdateAuxVars(geomech_realization)
      if (option%geomech_subsurf_coupling == GEOMECH_TWO_WAY_COUPLED) then 
        call GeomechStoreInitialPorosity(realization,geomech_realization)
        call GeomechUpdateSubsurfFromGeomech(realization,geomech_realization)
        call GeomechUpdateSubsurfPorosity(realization,geomech_realization)
      endif
    endif
#endif  

  ! print initial condition output if not a restarted sim
  call OutputInit(master_stepper%steps)
#ifdef GEOMECH
  if (option%ngeomechdof > 0) then
    call OutputGeomechInit(geomech_realization,master_stepper%steps)
  endif
#endif

  if (output_option%plot_number == 0 .and. &
      master_stepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
    call Output(realization,plot_flag,transient_plot_flag)
#ifdef GEOMECH
    if (option%ngeomechdof > 0) then
      geomech_plot_flag = PETSC_TRUE
      geomech_transient_plot_flag = PETSC_TRUE
      call OutputGeomechanics(geomech_realization,geomech_plot_flag, &
                              geomech_transient_plot_flag)
    endif
#endif    
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_stepper%max_time_step < 1) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'') 
    init_status = TIMESTEPPER_INIT_DONE
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_stepper)) then
    if (.not.associated(flow_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      init_status = TIMESTEPPER_INIT_FAIL
      return
    else
      flow_stepper%dt_max = flow_stepper%cur_waypoint%dt_max
    endif
  endif
  if (associated(tran_stepper)) then
    if (.not.associated(tran_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start ' // &
        'time or simulation time needs to be extended on a restart.'
      call printMsg(option)
      init_status = TIMESTEPPER_INIT_FAIL
      return
    else
      tran_stepper%dt_max = tran_stepper%cur_waypoint%dt_max
    endif
  endif
           
  if (associated(flow_stepper)) &
    flow_stepper%start_time_step = flow_stepper%steps + 1
  if (associated(tran_stepper)) &
    tran_stepper%start_time_step = tran_stepper%steps + 1
  
  if (realization%debug%print_couplers) then
    call OutputPrintCouplers(realization,ZERO_INTEGER)
  endif

end subroutine GeomechTimestepperInitializeRun

! ************************************************************************** !
!
! GeomechTimestepperExecuteRun: Runs the time step loop
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GeomechTimestepperExecuteRun(realization,geomech_realization, &
                                        master_stepper,flow_stepper, &
                                        tran_stepper,geomech_stepper)

  use Realization_class

  use Option_module
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Logging_module  
  use Discretization_module
  use Condition_Control_module
#ifdef GEOMECH
  use Geomechanics_Realization_module
  use Output_Geomechanics_module, only : OutputGeomechanics
  use Geomechanics_Force_module
#endif  

  implicit none
  
#include "finclude/petscdef.h"
#include "finclude/petsclog.h"
#include "finclude/petscsys.h"
#include "finclude/petscviewer.h"

  type(realization_type) :: realization
  type(geomech_realization_type) :: geomech_realization
  type(stepper_type), pointer :: master_stepper
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper  
  type(stepper_type), pointer :: null_stepper
  type(stepper_type), pointer :: geomech_stepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(waypoint_type), pointer :: prev_waypoint  
     

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, stop_flag, transient_plot_flag
  PetscBool :: step_to_steady_state
  PetscBool :: run_flow_as_steady_state
  PetscBool :: failure, surf_failure
  PetscReal :: tran_dt_save, flow_t0
  PetscReal :: flow_to_tran_ts_ratio
  PetscBool :: geomech_plot_flag
  PetscBool :: checkpoint_flag

  PetscLogDouble :: stepper_start_time, current_time, average_step_time
  PetscErrorCode :: ierr

  option => realization%option
  output_option => realization%output_option

  nullify(null_stepper)

  stop_flag = PETSC_FALSE
  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  failure = PETSC_FALSE
  run_flow_as_steady_state = PETSC_FALSE
  step_to_steady_state = PETSC_FALSE
  if (associated(flow_stepper)) then
    run_flow_as_steady_state = flow_stepper%run_as_steady_state
  endif
  checkpoint_flag = PETSC_FALSE

  call PetscTime(master_stepper%start_time, ierr)

  do

    if (OptionPrintToScreen(option) .and. &
        mod(master_stepper%steps,output_option%screen_imod) == 0) then
      option%print_screen_flag = PETSC_TRUE
    else
      option%print_screen_flag = PETSC_FALSE
    endif

    if (OptionPrintToFile(option) .and. &
        mod(master_stepper%steps,output_option%output_file_imod) == 0) then
      option%print_file_flag = PETSC_TRUE
    else
      option%print_file_flag = PETSC_FALSE
    endif

    prev_waypoint => master_stepper%cur_waypoint
    plot_flag = PETSC_FALSE
    transient_plot_flag = PETSC_FALSE

    call StepperSetTargetTimes(flow_stepper,tran_stepper, &
                               option,plot_flag, &
                               transient_plot_flag, &
                               checkpoint_flag)

    ! flow solution
    if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then

      flow_t0 = option%flow_time
      call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)
      call StepperStepFlowDT(realization,flow_stepper,failure)
      call PetscLogStagePop(ierr)
      if (failure) return ! if flow solve fails, exit
      option%flow_time = flow_stepper%target_time
    endif

    ! (reactive) transport solution
    if (associated(tran_stepper)) then
      call PetscLogStagePush(logging%stage(TRAN_STAGE),ierr)
      tran_dt_save = -999.d0
      ! reset transpor time step if flow time step is cut
      if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then
        if (flow_stepper%time_step_cut_flag) then
          tran_stepper%target_time = flow_stepper%target_time
          option%tran_dt = min(option%tran_dt,option%flow_dt)
        endif
      endif
      ! CFL Limiting
      if (tran_stepper%cfl_limiter > 0.d0) then
        call TimestepperCheckCFLLimit(tran_stepper,realization)
        call TimestepperEnforceCFLLimit(tran_stepper,option,output_option)
      endif
      ! if transport time step is less than flow step, but greater than
      ! half the flow time step, might as well cut it down to half the 
      ! flow time step size
      if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then
        flow_to_tran_ts_ratio = option%flow_dt / option%tran_dt
        if (flow_to_tran_ts_ratio > 1.d0 .and. &
            flow_to_tran_ts_ratio < 2.d0) then
          option%tran_dt = option%flow_dt * 0.5d0
        endif
      endif
      do ! loop on transport until it reaches the target time
        if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
          !global implicit
          call StepperStepTransportDT_GI(realization,tran_stepper, &
                                         run_flow_as_steady_state, &
                                         flow_t0,option%flow_time,failure)
        else
          !operator splitting
          call StepperStepTransportDT_OS(realization,tran_stepper, &
                                         run_flow_as_steady_state, &
                                         flow_t0,option%flow_time,failure)
        endif
        if (failure) then ! if transport solve fails, exit
          call PetscLogStagePop(ierr)
          return 
        endif
        
        ! update transport time
        option%tran_time = option%tran_time + option%tran_dt

        ! if target time reached, we are done
        if ((tran_stepper%target_time-option%tran_time) / &
            option%tran_time < 1.d-10) exit

        call StepperUpdateTransportSolution(realization,PETSC_TRUE)
        
        ! if still stepping, update the time step size based on convergence
        ! criteria
        call StepperUpdateDT(null_stepper,tran_stepper,option)
        ! check CFL limit if turned on
        if (tran_stepper%cfl_limiter > 0.d0) then
          call TimestepperEnforceCFLLimit(tran_stepper,option,output_option)
        endif         
        ! if current step exceeds (factoring in tolerance) exceeds the target
        ! time, set the step size to the difference
        if ((option%tran_time + &
             (1.d0+tran_stepper%time_step_tolerance)*option%tran_dt) > &
            tran_stepper%target_time) then
          ! tran_dt_save enables us to regain the original tran_dt if
          ! it is cut to match the target time
          tran_dt_save = option%tran_dt
          option%tran_dt = tran_stepper%target_time - option%tran_time
        endif
      enddo
      ! if substepping occured and time step size was modified to 
      ! to match the stepper target time, need to set the transport
      ! step size back to the pre-modified value
      !geh: moved below
      !geh: if (tran_dt_save > -998.d0) option%tran_dt = tran_dt_save
      option%tran_time = tran_stepper%target_time
      call PetscLogStagePop(ierr)
    endif

    if (realization%debug%print_couplers) then
      call OutputPrintCouplers(realization,master_stepper%steps)
    endif  

    ! update solution variables
    
    option%time = master_stepper%target_time
    call StepperUpdateSolution(realization,PETSC_TRUE)

    if (associated(tran_stepper)) then
      !geh: must revert tran_dt back after update of solution.  Otherwise,
      !     explicit update of mineral vol fracs, etc. will be updated based
      !     on incorrect dt
      if (tran_dt_save > -998.d0) option%tran_dt = tran_dt_save
    endif

    ! if a time step cut has occured, need to set the below back to original values
    ! if they changed. 
    if (master_stepper%time_step_cut_flag) then
      master_stepper%cur_waypoint => prev_waypoint
      plot_flag = PETSC_FALSE
    endif
    ! however, if we are using the modulus of the output_option%imod, we may still print
    if (mod(master_stepper%steps, &
            output_option%periodic_output_ts_imod) == 0) then
      plot_flag = PETSC_TRUE
    endif
    if (plot_flag .or. mod(master_stepper%steps, &
                           output_option%periodic_tr_output_ts_imod) == 0) then
      transient_plot_flag = PETSC_TRUE
    endif

!   deprecated - geh
!    if (option%compute_mass_balance) then
!      call MassBalanceUpdate(realization,flow_stepper%solver, &
!                             tran_stepper%solver)
!    endif
#ifdef GEOMECH       
    if (option%ngeomechdof > 0) &
      geomech_plot_flag = plot_flag
#endif 
    call Output(realization,plot_flag,transient_plot_flag)
    
    call StepperUpdateDTMax(flow_stepper,tran_stepper,option)
    call StepperUpdateDT(flow_stepper,tran_stepper,option)
  
#ifdef GEOMECH       
    if (option%ngeomechdof > 0) then
      select case (option%geomech_subsurf_coupling)
        case (GEOMECH_ONE_WAY_COUPLED) ! call geomech only at plot times
          if (geomech_plot_flag) then
            call GeomechUpdateFromSubsurf(realization,geomech_realization)
            call StepperSolveGeomechSteadyState(geomech_realization,geomech_stepper, &
                                                failure)
            call GeomechUpdateSolution(geomech_realization)
          endif
        case (GEOMECH_TWO_WAY_COUPLED)
          call GeomechUpdateFromSubsurf(realization,geomech_realization)
          call StepperSolveGeomechSteadyState(geomech_realization,geomech_stepper, &
                                              failure)
          call GeomechUpdateSolution(geomech_realization)  
          call GeomechUpdateSubsurfFromGeomech(realization,geomech_realization)
          call GeomechUpdateSubsurfPorosity(realization,geomech_realization)
        case default
          call StepperSolveGeomechSteadyState(geomech_realization,geomech_stepper, &
                                              failure)
          call GeomechUpdateSolution(geomech_realization)
      end select
      call OutputGeomechanics(geomech_realization,geomech_plot_flag, &
                              transient_plot_flag)
    endif
#endif    

    ! if a simulation wallclock duration time is set, check to see that the
    ! next time step will not exceed that value.  If it does, print the
    ! checkpoint and exit.
    if (option%wallclock_stop_flag) then
      call PetscTime(current_time, ierr)
      average_step_time = (current_time-master_stepper%start_time)/ &
                          real(master_stepper%steps-&
                               master_stepper%start_time_step+1) &
                          *2.d0  ! just to be safe, double it
      if (average_step_time + current_time > option%wallclock_stop_time) then
        call printMsg(option,"Wallclock stop time exceeded.  Exiting!!!")
        call printMsg(option,"")
        stop_flag = PETSC_TRUE
      endif
    endif

    if (option%checkpoint_flag .and. &
        mod(master_stepper%steps,option%checkpoint_frequency) == 0) then
      call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                             master_stepper%steps)  
    endif
    
    ! if at end of waypoint list (i.e. cur_waypoint = null), we are done!
    if (.not.associated(master_stepper%cur_waypoint) .or. stop_flag) exit

  enddo

end subroutine GeomechTimestepperExecuteRun

! ************************************************************************** !
!
! GeomechTimestepperFinalizeRun: Finalizes timestepping runs the time step loop
! author: Glenn Hammond
! date: 03/11/13
!
! ************************************************************************** !
subroutine GeomechTimestepperFinalizeRun(realization,geomech_realization, &
                                  master_stepper,flow_stepper, &
                                  tran_stepper,geomech_stepper)


  use Realization_class

  use Option_module
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Logging_module  
  use Geomechanics_Realization_module

  implicit none
  
#include "finclude/petscdef.h"
#include "finclude/petsclog.h"
#include "finclude/petscsys.h"
#include "finclude/petscviewer.h"

  type(realization_type) :: realization
  type(geomech_realization_type) :: geomech_realization
  type(stepper_type), pointer :: master_stepper
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(stepper_type), pointer :: geomech_stepper
  
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  PetscErrorCode :: ierr

  option => realization%option
  output_option => realization%output_option
  
  if (master_stepper%steps >= master_stepper%max_time_step) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
  endif

  if (option%checkpoint_flag) then
    call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                           NEG_ONE_INTEGER)  

  endif

  if (OptionPrintToScreen(option)) then
    if (option%nflowdof > 0) then
      write(*,'(/," FLOW steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            flow_stepper%steps,flow_stepper%cumulative_newton_iterations, &
            flow_stepper%cumulative_linear_iterations,flow_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') flow_stepper%cumulative_solver_time
      write(*,*) 'FLOW SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif
    if (option%ntrandof > 0) then
      write(*,'(/," TRAN steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            tran_stepper%steps,tran_stepper%cumulative_newton_iterations, &
            tran_stepper%cumulative_linear_iterations,tran_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') tran_stepper%cumulative_solver_time
      write(*,*) 'TRAN SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif            
  endif
  
  if (OptionPrintToFile(option)) then
    if (option%nflowdof > 0) then
      write(option%fid_out,'(/," FLOW steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            flow_stepper%steps,flow_stepper%cumulative_newton_iterations, &
            flow_stepper%cumulative_linear_iterations,flow_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') flow_stepper%cumulative_solver_time
      write(option%fid_out,*) 'FLOW SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif
    if (option%ntrandof > 0) then
      write(option%fid_out,'(/," TRAN steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            tran_stepper%steps,tran_stepper%cumulative_newton_iterations, &
            tran_stepper%cumulative_linear_iterations,tran_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') tran_stepper%cumulative_solver_time
      write(option%fid_out,*) 'TRAN SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif            
  endif

  ! pushed in TimeStepperInitializeRun
  call PetscLogStagePop(ierr)

end subroutine GeomechTimestepperFinalizeRun

! ************************************************************************** !
!
! StepperSolveGeomechSteadyState: Solves the steady-state force balance equation
! author: Satish Karra
! date: 06/19/13
!
! ************************************************************************** !
subroutine StepperSolveGeomechSteadyState(geomech_realization,stepper,failure)
  
  use Geomechanics_Realization_module
  use Geomechanics_Discretization_module
  use Option_module
  use Solver_module
  use Geomechanics_Field_module
  use Geomechanics_Force_module, only : GeomechanicsForceInitialGuess

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(geomech_realization_type) :: geomech_realization
  type(stepper_type) :: stepper
  PetscBool :: failure

  PetscErrorCode :: ierr
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscInt :: snes_reason
  PetscReal :: fnorm
  PetscReal :: inorm
  PetscReal :: scaled_fnorm

  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: field 
  type(solver_type), pointer :: solver
  type(geomech_discretization_type), pointer :: geomech_discretization

  option => geomech_realization%option
  geomech_discretization => geomech_realization%geomech_discretization
  field => geomech_realization%geomech_field
  solver => stepper%solver

  sum_newton_iterations = 0
  sum_linear_iterations = 0

    
  if (option%print_screen_flag) write(*,'(/,2("=")," GEOMECHANICS ",65("="))')

  call GeomechanicsForceInitialGuess(geomech_realization)


  call SNESSolve(solver%snes,PETSC_NULL_OBJECT,field%disp_xx,ierr)
     
  call SNESGetIterationNumber(solver%snes,num_newton_iterations,ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,ierr)
  call SNESGetConvergedReason(solver%snes,snes_reason,ierr)

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in GEOMECHANICS, reason: ', &
                snes_reason
    endif
    failure = PETSC_TRUE
    return
  endif
  
!  stepper%cumulative_newton_iterations = num_newton_iterations
!  stepper%cumulative_linear_iterations = num_linear_iterations

  ! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(field%disp_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    if (associated(geomech_discretization%grid)) then
       scaled_fnorm = fnorm/geomech_discretization%grid%nmax_node
    else
       scaled_fnorm = fnorm
    endif
    write(*,*) ''
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  

  sum_newton_iterations = sum_newton_iterations + num_newton_iterations
  sum_linear_iterations = sum_linear_iterations + num_linear_iterations
  
  stepper%steps = stepper%steps + 1      
  stepper%cumulative_newton_iterations = &
    stepper%cumulative_newton_iterations + sum_newton_iterations
  stepper%cumulative_linear_iterations = &
    stepper%cumulative_linear_iterations + sum_linear_iterations  
  
  if (option%print_screen_flag) print *, ""
  
  if (option%print_file_flag) then
    write(option%fid_out, '(" GEOMECHANICS ",i6," snes_conv_reason: ",i4,/, &
      &"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]")') &
      stepper%steps, &
      snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations
  endif  

end subroutine StepperSolveGeomechSteadyState

end module Geomechanics_Timestepper_module
#endif
