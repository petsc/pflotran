module Timestepper_module
 
  use Solver_module
  use Waypoint_module 
  use Convergence_module 
 
  implicit none

  private
  
#include "definitions.h"
 
  type, public :: stepper_type
  
    PetscInt :: steps         ! The number of time-steps taken by the code.
    PetscInt :: nstepmax      ! Maximum number of timesteps taken by the code.
    PetscInt :: icut_max      ! Maximum number of timestep cuts within one time step.
    PetscInt :: ndtcmx        ! Steps needed after cutting to increase time step
    PetscInt :: newton_cum       ! Total number of Newton iterations
    PetscInt :: linear_cum     ! Total number of linear iterations
    PetscInt :: icutcum       ! Total number of cuts in the timestep taken.    
    PetscInt :: iaccel        ! Accelerator index
    
    
    PetscReal :: dt_min
    PetscReal :: dt_max
    
    PetscReal :: cumulative_solver_time
    
    PetscTruth :: init_to_steady_state
    PetscTruth :: run_as_steady_state
    PetscReal :: steady_state_rel_tol

    ! An array of multiplicative factors that specify how to increase time step.
    PetscReal, pointer :: tfac(:)
        
    type(solver_type), pointer :: solver
    
    type(waypoint_type), pointer :: cur_waypoint

    type(convergence_context_type), pointer :: convergence_context
    
  end type stepper_type
  
  public :: TimestepperCreate, TimestepperDestroy, StepperRun, &
            TimestepperRead, TimestepperPrintInfo

contains

! ************************************************************************** !
!
! TimestepperCreate: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function TimestepperCreate()

  implicit none
  
  type(stepper_type), pointer :: TimestepperCreate
  
  type(stepper_type), pointer :: stepper
  
  allocate(stepper)
  stepper%steps = 0
  stepper%nstepmax = 999999
  
  stepper%icut_max = 16
  stepper%ndtcmx = 5
  stepper%newton_cum = 0
  stepper%linear_cum = 0
  stepper%icutcum = 0    
  stepper%iaccel = 5
  
  stepper%dt_min = 1.d0
  stepper%dt_max = 3.1536d6 ! One-tenth of a year.  
  
  stepper%cumulative_solver_time = 0.d0

  allocate(stepper%tfac(13))
  stepper%tfac(1)  = 2.0d0; stepper%tfac(2)  = 2.0d0
  stepper%tfac(3)  = 2.0d0; stepper%tfac(4)  = 2.0d0
  stepper%tfac(5)  = 2.0d0; stepper%tfac(6)  = 1.8d0
  stepper%tfac(7)  = 1.6d0; stepper%tfac(8)  = 1.4d0
  stepper%tfac(9)  = 1.2d0; stepper%tfac(10) = 1.0d0
  stepper%tfac(11) = 1.0d0; stepper%tfac(12) = 1.0d0
  stepper%tfac(13) = 1.0d0
  
  stepper%init_to_steady_state = PETSC_FALSE
  stepper%steady_state_rel_tol = 1.d-30
  stepper%run_as_steady_state = PETSC_FALSE
  
  nullify(stepper%solver)
  nullify(stepper%convergence_context)
  nullify(stepper%cur_waypoint)
  
  stepper%solver => SolverCreate()
  
  TimeStepperCreate => stepper
  
end function TimestepperCreate 

! ************************************************************************** !
!
! TimestepperRead: Reads parameters associated with time stepper
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine TimestepperRead(stepper,input,option)

  use Option_module
  use String_module
  use Input_module
  
  implicit none

  type(stepper_type) :: stepper
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','TIMESTEPPER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))

      case('NUM_STEPS_AFTER_TS_CUT')
        call InputReadInt(input,option,stepper%ndtcmx)
        call InputDefaultMsg(input,option,'num_steps_after_ts_cut')

      case('MAX_STEPS')
        call InputReadInt(input,option,stepper%nstepmax)
        call InputDefaultMsg(input,option,'nstepmax')
  
      case('TS_ACCELERATION')
        call InputReadInt(input,option,stepper%iaccel)
        call InputDefaultMsg(input,option,'iaccel')

      case('MAX_TS_CUTS')
        call InputReadInt(input,option,stepper%icut_max)
        call InputDefaultMsg(input,option,'icut_max')

      case('INITIALIZE_TO_STEADY_STATE')
        stepper%init_to_steady_state = PETSC_TRUE
        call InputReadDouble(input,option,stepper%steady_state_rel_tol)
        call InputDefaultMsg(input,option,'steady state convergence relative tolerance')

      case('RUN_AS_STEADY_STATE')
        stepper%run_as_steady_state = PETSC_TRUE

      case('MAX_PRESSURE_CHANGE')
        call InputReadDouble(input,option,option%dpmxe)
        call InputDefaultMsg(input,option,'dpmxe')

      case('MAX_TEMPERATURE_CHANGE')
        call InputReadDouble(input,option,option%dtmpmxe)
        call InputDefaultMsg(input,option,'dtmpmxe')
  
      case('MAX_CONCENTRATION_CHANGE')
        call InputReadDouble(input,option,option%dcmxe)
        call InputDefaultMsg(input,option,'dcmxe')

      case('MAX_SATURATION_CHANGE')
        call InputReadDouble(input,option,option%dsmxe)
        call InputDefaultMsg(input,option,'dsmxe')

      case default
        option%io_buffer = 'Timestepper option: '//trim(keyword)// &
                           ' not recognized.'
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine TimestepperRead

! ************************************************************************** !
!
! StepperRun: Runs the time step loop
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine StepperRun(realization,flow_stepper,tran_stepper)

  use Realization_module

  use Option_module
  use Output_module, only : Output, OutputInit, OutputVectorTecplot
  use Logging_module  
  use Mass_Balance_module
  use Discretization_module
  
  implicit none
  
#include "finclude/petscdef.h"
#include "finclude/petsclog.h"
#include "finclude/petscsys.h"
#include "finclude/petscviewer.h"

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(stepper_type), pointer :: master_stepper
  type(stepper_type), pointer :: flow_stepper_save
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(waypoint_type), pointer :: prev_waypoint  

  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: plot_flag, stop_flag, transient_plot_flag
  PetscTruth :: master_timestep_cut_flag
  PetscTruth :: flow_timestep_cut_flag, tran_timestep_cut_flag
  PetscInt :: istep, start_step
  PetscInt :: num_const_timesteps
  PetscInt :: num_newton_iterations, idum, idum2
  PetscTruth :: activity_coefs_read
  PetscTruth :: flow_read
  PetscTruth :: transport_read
  PetscTruth :: step_to_steady_state
  PetscTruth :: run_flow_as_steady_state
  PetscTruth :: failure
  PetscLogDouble :: start_time, end_time
  
  PetscLogDouble :: stepper_start_time, current_time, average_step_time
  PetscErrorCode :: ierr

  option => realization%option
  output_option => realization%output_option
  
  if (option%steady_state) then
    call StepperRunSteadyState(realization,flow_stepper,tran_stepper)
    return 
  endif
  
  if (associated(flow_stepper)) then
    master_stepper => flow_stepper
  else
    master_stepper => tran_stepper
  endif

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  master_timestep_cut_flag = PETSC_FALSE
  flow_timestep_cut_flag = PETSC_FALSE
  tran_timestep_cut_flag = PETSC_FALSE
  stop_flag = PETSC_FALSE
  activity_coefs_read = PETSC_FALSE
  flow_read = PETSC_FALSE
  transport_read = PETSC_FALSE
  step_to_steady_state = PETSC_FALSE
  run_flow_as_steady_state = PETSC_FALSE
  failure = PETSC_FALSE
  num_const_timesteps = 0  

  if (option%restart_flag) then
    call StepperRestart(realization,flow_stepper,tran_stepper, &
                        num_const_timesteps,num_newton_iterations, &
                        flow_read,transport_read,activity_coefs_read)
    if (associated(flow_stepper)) flow_stepper%cur_waypoint => &
      WaypointSkipToTime(realization%waypoints,option%time)
    if (associated(tran_stepper)) tran_stepper%cur_waypoint => &
      WaypointSkipToTime(realization%waypoints,option%time)

    if (flow_read) then
      call StepperUpdateFlowAuxVars(realization)
    endif

  else if (master_stepper%init_to_steady_state) then
    option%print_screen_flag = OptionPrintToScreen(option)
    option%print_file_flag = OptionPrintToFile(option)
    if (associated(flow_stepper)) then
      if (flow_stepper%init_to_steady_state) then
        step_to_steady_state = PETSC_TRUE
        option%flow_dt = master_stepper%dt_min
        call StepperStepFlowDT(realization,flow_stepper, &
                               flow_timestep_cut_flag, &
                               num_newton_iterations, &
                               step_to_steady_state,failure)
        if (failure) then ! if flow solve fails, exit
          if (OptionPrintToScreen(option)) then
            write(*,*) ' ERROR: steady state solve failed!!!'
          endif
          if (OptionPrintToFile(option)) then
            write(option%fid_out,*) ' ERROR: steady state solve failed!!!'
          endif
          return 
        endif
        option%flow_dt = master_stepper%dt_min
        run_flow_as_steady_state = flow_stepper%run_as_steady_state
      endif
      if (run_flow_as_steady_state) then
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
    call RealizAssignTransportInitCond(realization)  
  endif

  ! turn on flag to tell RTUpdateSolution that the code is not timestepping
  call StepperUpdateSolution(realization)

  if (option%jumpstart_kinetic_sorption .and. option%time < 1.d-40) then
    call StepperJumpStart(realization)
  endif
  
  call PetscLogStagePop(ierr)
  option%init_stage = PETSC_FALSE
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr)

  ! print initial condition output if not a restarted sim
  call OutputInit(realization)
  if (output_option%plot_number == 0 .and. &
      master_stepper%nstepmax >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
    call Output(realization,plot_flag,transient_plot_flag)
    if (output_option%print_permeability) then
      if (len_trim(option%group_prefix) > 1) then
        string = 'permeability-' // trim(option%group_prefix) // '.tec'
      else
        string = 'permeability.tec'
      endif
      call DiscretizationLocalToGlobal(realization%discretization, &
                                       realization%field%perm_xx_loc, &
                                       realization%field%work,ONEDOF)
      call OutputVectorTecplot(string,string,realization,realization%field%work)
    endif
#if 0
! this now occurs in the standard output file    
    if (output_option%print_porosity) then
      if (len_trim(option%group_prefix) > 1) then
        string = 'porosity-' // trim(option%group_prefix) // '.tec'
      else
        string = 'porosity.tec'
      endif
      call DiscretizationLocalToGlobal(realization%discretization, &
                                       realization%field%porosity_loc, &
                                       realization%field%work,ONEDOF)
      call OutputVectorTecplot(string,string,realization,realization%field%work)
    endif
#endif    
  endif
  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_stepper)) then
    if (.not.associated(flow_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      return
    else
      flow_stepper%dt_max = flow_stepper%cur_waypoint%dt_max
    endif
  endif
  if (associated(tran_stepper)) then
    if (.not.associated(tran_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start time.'
      call printMsg(option)
      return
    else
      tran_stepper%dt_max = tran_stepper%cur_waypoint%dt_max
    endif
  endif
           
  ! ensure that steady_state flag is off
  step_to_steady_state = PETSC_FALSE
  call PetscGetTime(stepper_start_time, ierr)
  start_step = master_stepper%steps+1
  do istep = start_step, master_stepper%nstepmax

    if (OptionPrintToScreen(option) .and. &
        mod(master_stepper%steps,output_option%screen_imod) == 0) then
      option%print_screen_flag = PETSC_TRUE
    else
      option%print_screen_flag = PETSC_FALSE
    endif

    if (OptionPrintToFile(option)) then
      option%print_file_flag = PETSC_TRUE
    else
      option%print_file_flag = PETSC_FALSE
    endif

    prev_waypoint => master_stepper%cur_waypoint
    flow_timestep_cut_flag = PETSC_FALSE
    tran_timestep_cut_flag = PETSC_FALSE
    plot_flag = PETSC_FALSE
    transient_plot_flag = PETSC_FALSE
    
    call StepperSetTargetTimes(flow_stepper,tran_stepper,option,plot_flag, &
                               transient_plot_flag)

    ! flow solution
    if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then
      call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)
      call StepperStepFlowDT(realization,flow_stepper, &
                             flow_timestep_cut_flag, &
                             num_newton_iterations, &
                             step_to_steady_state,failure)
      call PetscLogStagePop(ierr)
      if (failure) return ! if flow solve fails, exit
    endif
    ! (reactive) transport solution
    if (associated(tran_stepper)) then
      call PetscLogStagePush(logging%stage(TRAN_STAGE),ierr)
      call StepperStepTransportDT(realization,tran_stepper, &
                                  flow_timestep_cut_flag, &
                                  tran_timestep_cut_flag, &
                                  idum,failure)
      call PetscLogStagePop(ierr)
      if (failure) return ! if flow solve fails, exit
      if (.not.associated(flow_stepper) .or. run_flow_as_steady_state) &
        num_newton_iterations = idum
    endif

    if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then
      master_timestep_cut_flag = flow_timestep_cut_flag
    else
      master_timestep_cut_flag = tran_timestep_cut_flag
    endif

    ! update solution variables
    call StepperUpdateSolution(realization)
    
    ! if a time step cut has occured, need to set the below back to original values
    ! if they changed. 
    if (master_timestep_cut_flag) then
      master_stepper%cur_waypoint => prev_waypoint
      plot_flag = PETSC_FALSE
    endif
    ! however, if we are using the modulus of the output_option%imod, we may still print
    if (mod(istep,output_option%periodic_output_ts_imod) == 0) then
      plot_flag = PETSC_TRUE
    endif
    if (plot_flag .or. mod(istep,output_option%periodic_tr_output_ts_imod) == 0) then
      transient_plot_flag = PETSC_TRUE
    endif

!   deprecated - geh
!    if (option%compute_mass_balance) then
!      call MassBalanceUpdate(realization,flow_stepper%solver, &
!                             tran_stepper%solver)
!    endif
    call Output(realization,plot_flag,transient_plot_flag)
  
    call StepperUpdateDT(flow_stepper,tran_stepper,option, &
                         master_timestep_cut_flag, &
                         num_const_timesteps,num_newton_iterations)

    ! if a simulation wallclock duration time is set, check to see that the
    ! next time step will not exceed that value.  If it does, print the
    ! checkpoint and exit.
    if (option%wallclock_stop_flag) then
      call PetscGetTime(current_time, ierr)
      average_step_time = (current_time-stepper_start_time)/ &
                          real(istep-start_step+1) &
                          *2.d0  ! just to be safe, double it
      if (average_step_time + current_time > option%wallclock_stop_time) then
        call printMsg(option,"Wallclock stop time exceeded.  Exiting!!!")
        call printMsg(option,"")
        stop_flag = PETSC_TRUE
      endif
    endif

    if (option%checkpoint_flag .and. &
        mod(istep,option%checkpoint_frequency) == 0) then
      call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                             num_const_timesteps,num_newton_iterations, &
                             istep)  
    endif
    
   
    ! if at end of waypoint list (i.e. cur_waypoint = null), we are done!
    if (.not.associated(master_stepper%cur_waypoint) .or. stop_flag) exit

  enddo

  if (option%checkpoint_flag) then
    call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                           num_const_timesteps,num_newton_iterations, &
                           NEG_ONE_INTEGER)  
  endif

  if (OptionPrintToScreen(option)) then
    if (option%nflowdof > 0) then
      write(*,'(/," FLOW steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            flow_stepper%steps,flow_stepper%newton_cum, &
            flow_stepper%linear_cum,flow_stepper%icutcum
      write(string,'(f12.1)') flow_stepper%cumulative_solver_time
      write(*,*) 'FLOW SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif
    if (option%ntrandof > 0) then
      write(*,'(/," TRAN steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            tran_stepper%steps,tran_stepper%newton_cum, &
            tran_stepper%linear_cum,tran_stepper%icutcum
      write(string,'(f12.1)') tran_stepper%cumulative_solver_time
      write(*,*) 'TRAN SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif            
  endif
  
  if (OptionPrintToFile(option)) then
    if (option%nflowdof > 0) then
      write(option%fid_out,'(/," FLOW steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            flow_stepper%steps,flow_stepper%newton_cum, &
            flow_stepper%linear_cum,flow_stepper%icutcum
      write(string,'(f12.1)') flow_stepper%cumulative_solver_time
      write(option%fid_out,*) 'FLOW SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif
    if (option%ntrandof > 0) then
      write(option%fid_out,'(/," TRAN steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            tran_stepper%steps,tran_stepper%newton_cum, &
            tran_stepper%linear_cum,tran_stepper%icutcum
      write(string,'(f12.1)') tran_stepper%cumulative_solver_time
      write(option%fid_out,*) 'TRAN SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif            
  endif

  call PetscLogStagePop(ierr)

end subroutine StepperRun

! ************************************************************************** !
!
! StepperUpdateDT: Updates time step
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperUpdateDT(flow_stepper,tran_stepper,option,timestep_cut_flag, &
                           num_const_timesteps,num_newton_iterations)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(option_type) :: option
  PetscTruth :: timestep_cut_flag
  PetscInt :: num_const_timesteps
  PetscInt :: num_newton_iterations
  
  type(stepper_type), pointer :: stepper
  PetscReal :: time, dt
  PetscReal :: fac,dtt,up,utmp,uc,ut,uus
  PetscTruth :: update_dt_with_flow_stepper

  if (num_const_timesteps > 0) num_const_timesteps = num_const_timesteps + 1
  if (timestep_cut_flag) num_const_timesteps = 1

  update_dt_with_flow_stepper = PETSC_TRUE
  if (associated(flow_stepper)) then
    if (flow_stepper%run_as_steady_state) then
      update_dt_with_flow_stepper = PETSC_FALSE
    endif
  else
    update_dt_with_flow_stepper = PETSC_FALSE
  endif

  if (update_dt_with_flow_stepper) then
    if (num_const_timesteps > flow_stepper%ndtcmx) then
      num_const_timesteps = 0
    else if (num_const_timesteps > 0) then
      return
    endif
    stepper => flow_stepper
    time = option%flow_time
    dt = option%flow_dt
  else
    if (num_const_timesteps > tran_stepper%ndtcmx) then
      num_const_timesteps = 0
    else if (num_const_timesteps > 0) then
      return
    endif
    stepper => tran_stepper
    time = option%tran_time
    dt = option%tran_dt
  endif

  if (stepper%iaccel == 0) return

  select case(option%iflowmode)
    case(IMS_MODE)   
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = option%dpmxe/(option%dpmax+0.1)
        utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
        uus= option%dsmxe/(option%dsmax+1.d-6)
        ut = min(up,utmp,uus)
      endif
      dtt = fac * dt * (1.d0 + ut)
    case(MPH_MODE)   
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = option%dpmxe/(option%dpmax+0.1)
        utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
        uc = option%dcmxe/(option%dcmax+1.d-6)
        uus= option%dsmxe/(option%dsmax+1.d-6)
        ut = min(up,utmp,uc,uus)
      endif
      dtt = fac * dt * (1.d0 + ut)
    case(THC_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = option%dpmxe/(option%dpmax+0.1)
        utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
        uus= option%dsmxe/(option%dsmax+1.d-6)
        ut = min(up,utmp,uus)
      endif
      dtt = fac * dt * (1.d0 + ut)
    case(RICHARDS_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = option%dpmxe/(option%dpmax+0.1)
        ut = up
      endif
      dtt = fac * dt * (1.d0 + ut)
    case default
      dtt = dt
      if (num_newton_iterations <= stepper%iaccel .and. &
          num_newton_iterations <= size(stepper%tfac)) then
        if (num_newton_iterations == 0) then
          dtt = stepper%tfac(1) * dt
        else
          dtt = stepper%tfac(num_newton_iterations) * dt
        endif
      endif
      
  end select
  
  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > stepper%dt_max) dtt = stepper%dt_max
  if (dtt>.25d0*time .and. time>5.d2) dtt=.25d0*time
  dt = dtt

  if (update_dt_with_flow_stepper) then
    option%flow_dt = dt
  else
    option%tran_dt = dt
  endif

end subroutine StepperUpdateDT

! ************************************************************************** !
!
! StepperSetTargetTimes: Sets target time for flow and transport solvers
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperSetTargetTimes(flow_stepper,tran_stepper,option,plot_flag, &
                                 transient_plot_flag)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: flow_stepper, tran_stepper
  type(option_type) :: option
  PetscTruth :: plot_flag
  PetscTruth :: transient_plot_flag
  
  PetscReal :: time
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: steps
  PetscInt :: nstepmax
  PetscTruth :: set_target_with_flow_stepper
  type(waypoint_type), pointer :: cur_waypoint

  ! target time will always be dictated by the flow solver, if present
  set_target_with_flow_stepper = PETSC_TRUE
  if (associated(flow_stepper)) then
    if (flow_stepper%run_as_steady_state) then
      set_target_with_flow_stepper = PETSC_FALSE
    endif
  else
    set_target_with_flow_stepper = PETSC_FALSE
  endif

  if (set_target_with_flow_stepper) then
    time = option%flow_time + option%flow_dt
    dt = option%flow_dt
    dt_max = flow_stepper%dt_max
    cur_waypoint => flow_stepper%cur_waypoint
    steps = flow_stepper%steps
    nstepmax = flow_stepper%nstepmax
  else
    time = option%tran_time + option%tran_dt
    dt = option%tran_dt
    dt_max = tran_stepper%dt_max
    cur_waypoint => tran_stepper%cur_waypoint
    steps = tran_stepper%steps
    nstepmax = tran_stepper%nstepmax
  endif
  
  ! FOr the case where the second waypoint is a printout after the first time step
  ! we must increment the waypoint beyond the first (time=0.) waypoint.  Otherwise
  ! the second time step will be zero. - geh
  if (cur_waypoint%time < 1.d-40) then
    cur_waypoint => cur_waypoint%next
  endif
  
  if (option%match_waypoint) then
    time = time - dt
    dt = option%prev_dt
    time = time + dt
    option%match_waypoint = PETSC_FALSE
  else
    option%prev_dt = dt
  endif

! If a waypoint calls for a plot or change in src/sinks, adjust time step to match waypoint
  if (time + 0.2*dt >= cur_waypoint%time .and. &
      (cur_waypoint%update_srcs .or. &
       cur_waypoint%print_output .or. &
       cur_waypoint%print_tr_output)) then
    time = time - dt
    dt = cur_waypoint%time - time
    if (dt > dt_max .and. dabs(dt-dt_max) > 1.d0) then ! 1 sec tolerance to avoid cancellation
      dt = dt_max                    ! error from waypoint%time - time
      time = time + dt
    else
      time = cur_waypoint%time
      if (cur_waypoint%print_output) plot_flag = PETSC_TRUE
      if (cur_waypoint%print_tr_output) transient_plot_flag = PETSC_TRUE
      option%match_waypoint = PETSC_TRUE
      cur_waypoint => cur_waypoint%next
      if (associated(cur_waypoint)) &
        dt_max = cur_waypoint%dt_max
    endif
  else if (time > cur_waypoint%time) then
    cur_waypoint => cur_waypoint%next
    if (associated(cur_waypoint)) &
      dt_max = cur_waypoint%dt_max
  else if (steps >= nstepmax) then
    plot_flag = PETSC_TRUE
    nullify(cur_waypoint)
  endif
    
  ! target time will always be dictated by the flow solver, if present
  if (associated(flow_stepper)) then
    option%flow_time = time
    option%flow_dt = dt
    flow_stepper%dt_max = dt_max
    flow_stepper%cur_waypoint => cur_waypoint
  endif
  if (associated(tran_stepper)) then
    option%tran_time = time
    option%tran_dt = dt
    tran_stepper%dt_max = dt_max
    tran_stepper%cur_waypoint => cur_waypoint
  endif
  
  option%time = time
  option%dt = dt
  
end subroutine StepperSetTargetTimes

! ************************************************************************** !
!
! StepperStepFlowDT: Steps forward one step in time
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperStepFlowDT(realization,stepper,timestep_cut_flag, &
                             num_newton_iterations,step_to_steady_state, &
                             failure)
  use MPHASE_module, only : MphaseMaxChange, MphaseInitializeTimestep, &
                           MphaseTimeCut, MPhaseUpdateReason
  use Immis_module, only : ImmisMaxChange, ImmisInitializeTimestep, &
                           ImmisTimeCut, ImmisUpdateReason
  use Richards_module, only : RichardsMaxChange, RichardsInitializeTimestep, &
                             RichardsTimeCut
  use THC_module, only : THCMaxChange, THCInitializeTimestep, THCTimeCut
  use Global_module

  use Output_module, only : Output
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper

  PetscTruth :: timestep_cut_flag
  PetscInt :: num_newton_iterations
  PetscTruth :: step_to_steady_state
  PetscTruth :: failure
  
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason
  PetscInt :: sum_newton_iterations, sum_linear_iterations, num_linear_iterations
  PetscReal :: fnorm, scaled_fnorm, inorm, prev_norm, dif_norm, rel_norm
  Vec :: update_vec
  PetscTruth :: plot_flag
  PetscTruth :: transient_plot_flag
  PetscLogDouble :: log_start_time, log_end_time

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(discretization_type), pointer :: discretization 
  type(solver_type), pointer :: solver

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  num_newton_iterations = 0
  icut = 0
  sum_newton_iterations = 0
  sum_linear_iterations = 0
  prev_norm = 1.d20

  do ! this loop is for steady state initial condition
  
    ! Perform some global-to-local scatters to update the ghosted vectors.
    ! We have to do this so that the routines for calculating the residual
    ! and the Jacobian will have the ghost points they need.
    ! Note that we don't do the global-to-local scatter for the pressure 
    ! vector, as that needs to be done within the residual calculation routine
    ! because that routine may get called several times during one Newton step
    ! if a method such as line search is being used.
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%tor_loc, &
                                    field%tor_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                    field%icap_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                    field%ithrm_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                    field%iphas_loc,ONEDOF)
    
    if (option%print_screen_flag) write(*,'(/,2("=")," FLOW ",52("="))')

    if (option%ntrandof > 0) then ! store initial saturations for transport
      call GlobalUpdateAuxVars(realization,TIME_T)
    endif
    
    select case(option%iflowmode)
      case(THC_MODE)
        call THCInitializeTimestep(realization)
      case(RICHARDS_MODE)
        call RichardsInitializeTimestep(realization)
      case(MPH_MODE)
        call MphaseInitializeTimestep(realization)
      case(IMS_MODE)
        call ImmisInitializeTimestep(realization)
    end select
    
    do
      
      call PetscGetTime(log_start_time, ierr)
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
          call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%flow_xx, ierr)
      end select
      call PetscGetTime(log_end_time, ierr)
      stepper%cumulative_solver_time = stepper%cumulative_solver_time + (log_end_time - log_start_time)

  ! do we really need all this? - geh 
      call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
      call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
      call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

      sum_newton_iterations = sum_newton_iterations + num_newton_iterations
      sum_linear_iterations = sum_linear_iterations + num_linear_iterations
      update_reason = 1
      
      if (snes_reason >= 0) then
        select case(option%iflowmode)
          case(IMS_MODE)
            call ImmisUpdateReason(update_reason,realization)
          case(MPH_MODE)
            call MPhaseUpdateReason(update_reason,realization)
          case(THC_MODE)
            update_reason=1
          case(RICHARDS_MODE)
            update_reason=1
        end select   
        if (option%print_screen_flag) print *,'update_reason: ',update_reason
      endif
   
  !******************************************************************
      
      if (snes_reason <= 0 .or. update_reason <= 0) then
        ! The Newton solver diverged, so try reducing the time step.
        icut = icut + 1
        timestep_cut_flag = PETSC_TRUE

        if (icut > stepper%icut_max .or. option%flow_dt<1.d-20) then
          if (option%print_screen_flag) then
            print *,"--> icut_max exceeded: icut/icutmax= ",icut, &
                    stepper%icut_max, "t= ", &
                    option%flow_time/realization%output_option%tconv, " dt= ", &
                    option%flow_dt/realization%output_option%tconv
            print *,"Stopping execution!"
          endif
          realization%output_option%plot_name = 'flow_cut_to_failure'
          plot_flag = PETSC_TRUE
          transient_plot_flag = PETSC_FALSE
          call Output(realization,plot_flag,transient_plot_flag)
          failure = PETSC_TRUE
          return
        endif

        option%flow_time = option%flow_time - option%flow_dt
        option%flow_dt = 0.5d0 * option%flow_dt
      
        if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
          &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.4, '' dt= '', &
          &   1pe12.4,l2)')  snes_reason,icut,stepper%icutcum, &
              option%flow_time/realization%output_option%tconv, &
              option%flow_dt/realization%output_option%tconv,timestep_cut_flag

        option%flow_time = option%flow_time + option%flow_dt

        select case(option%iflowmode)
          case(THC_MODE)
            call THCTimeCut(realization)
          case(RICHARDS_MODE)
            call RichardsTimeCut(realization)
          case(MPH_MODE)
            call MphaseTimeCut(realization)
          case(IMS_MODE)
            call ImmisTimeCut(realization)
        end select
        call VecCopy(field%iphas_old_loc, field%iphas_loc, ierr)

      else
        ! The Newton solver converged, so we can exit.
        exit
      endif
    enddo
    
    if (.not.step_to_steady_state) then
      exit
    else
    
      call SNESGetSolutionUpdate(solver%snes,update_vec,ierr)
      call VecStrideNorm(update_vec,ZERO_INTEGER,NORM_INFINITY,inorm,ierr)
      dif_norm = inorm-prev_norm
      rel_norm = dif_norm/prev_norm
      if ((sum_newton_iterations > 100 .and. &
           dabs(dif_norm) < stepper%steady_state_rel_tol) .or. &
          dabs(rel_norm) < stepper%steady_state_rel_tol) then
        if (option%print_file_flag) then
          write(option%fid_out,*) 'Steady state solve converged after ', &
            sum_newton_iterations, ' iterations'
        endif
        exit
      endif
      
      prev_norm = inorm
      
      if (sum_newton_iterations > 1000) then
        option%io_buffer = 'Steady-state solve exceeded 1000 Newton iterations.'
        call printMsg(option)
        if (option%print_file_flag) then
          write(option%fid_out,*) trim(option%io_buffer)
          write(option%fid_out,*) 'inorm: ', inorm
          write(option%fid_out,*) 'prev_norm: ', prev_norm
          write(option%fid_out,*) 'rel_norm: ', rel_norm
          write(option%fid_out,*) 'dif_norm: ', dif_norm
          write(option%fid_out,*) 'relative tolerance: ', &
                                  stepper%steady_state_rel_tol
        endif
        failure = PETSC_TRUE
        exit
      endif

      ! zero out dt to prevent potential error in mass balance calc while 
      ! updating solution
      fnorm = option%flow_dt
      option%flow_dt = 0.d0  
      ! take next step at larger dt      
      call StepperUpdateFlowSolution(realization)
      option%flow_dt = fnorm
      option%flow_dt = option%flow_dt*2.d0
      if (option%print_file_flag) then
        write(option%fid_out,*) 'Dt: ', option%flow_dt
      endif
      if (option%print_screen_flag) then
        write(*,*) 'Dt: ', option%flow_dt
      endif
    endif

  enddo ! end of steady-state loop
  
  if (.not.step_to_steady_state) then
    stepper%steps = stepper%steps + 1      
  endif
  stepper%newton_cum = stepper%newton_cum + sum_newton_iterations
  stepper%linear_cum = stepper%linear_cum + sum_linear_iterations
  stepper%icutcum = stepper%icutcum + icut

  if (option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(realization,TIME_TpDT)
  endif
    
! print screen output
  call VecNorm(field%flow_r,NORM_2,fnorm,ierr) 
  call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    write(*, '(/," FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps,option%flow_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%newton_cum,sum_linear_iterations,stepper%linear_cum,icut, &
      stepper%icutcum

    ! the grid pointer is null if we are working with SAMRAI
    if(associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  if (option%print_file_flag) then
    write(option%fid_out, '(" FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') stepper%steps, &
      option%flow_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%newton_cum,sum_linear_iterations,stepper%linear_cum,icut, &
      stepper%icutcum
  endif
  
  if (option%iflowmode == THC_MODE) then
     call THCMaxChange(realization)
    if (option%print_screen_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax, option%dcmax
    endif
    if (option%print_file_flag) then 
      write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dcmax
    endif
  else if (option%iflowmode == RICHARDS_MODE) then
    call RichardsMaxChange(realization)
    if (option%print_screen_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
    endif
    if (option%print_file_flag) then
      write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
    endif
  else if (option%iflowmode == MPH_MODE) then
    call MphaseMaxChange(realization)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (option%print_screen_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
    endif
    if (option%print_file_flag) then  
      write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
    endif
  else if (option%iflowmode == IMS_MODE) then
    call ImmisMaxChange(realization)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (option%print_screen_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dsmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dsmax
    endif
    if (option%print_file_flag) then  
      write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dsmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dsmax
    endif
  endif

  if (option%print_screen_flag) print *, ""

end subroutine StepperStepFlowDT

! ************************************************************************** !
!
! StepperStepTransportDT: Steps forward one step in time
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperStepTransportDT(realization,stepper,flow_timestep_cut_flag, &
                                  tran_timestep_cut_flag, &
                                  num_newton_iterations,failure)
  
  use Reactive_Transport_module
  use Output_module, only : Output
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  use Grid_module
  use Level_module
  use Patch_module
  use Global_module  
  
  implicit none

#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper

  PetscTruth :: flow_timestep_cut_flag
  PetscTruth :: tran_timestep_cut_flag
  PetscInt :: num_newton_iterations
  PetscTruth :: failure
  
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason
  PetscInt :: sum_newton_iterations, sum_linear_iterations, num_linear_iterations
  PetscInt :: n, nmax_inf
  PetscReal :: fnorm, scaled_fnorm, inorm
  PetscReal :: start_time, end_time, dt_orig
  PetscTruth :: plot_flag  
  PetscTruth :: transient_plot_flag  
  PetscReal, pointer :: r_p(:), xx_p(:), log_xx_p(:)
  PetscReal, parameter :: time_tol = 1.d-10
  PetscLogDouble :: log_start_time, log_end_time

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver
  type(patch_type), pointer :: cur_patch
  type(level_type), pointer :: cur_level

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  call DiscretizationLocalToLocal(discretization,field%porosity_loc,field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tor_loc,field%tor_loc,ONEDOF)

  if (flow_timestep_cut_flag) then
    option%tran_time = option%flow_time
    option%tran_dt = option%flow_dt
    option%time = option%flow_time
  endif
  
  end_time = option%tran_time
  start_time = end_time-option%tran_dt
  dt_orig = option%tran_dt
  
  ! test
!  option%tran_time = option%tran_time - option%tran_dt
!  option%tran_dt = option%tran_dt * 0.5d0
!  option%tran_time = option%tran_time + option%tran_dt
  
  do
  
    if (option%nflowdof > 0) then
      option%tran_weight_t0 = (option%tran_time-option%tran_dt-start_time)/ &
                              (end_time-start_time)
      option%tran_weight_t1 = (option%tran_time-start_time)/ &
                              (end_time-start_time)
      ! set densities and saturations to t
      call GlobalUpdateDenAndSat(realization,option%tran_weight_t0)
    endif

    call RTInitializeTimestep(realization)

    ! set densities and saturations to t+dt
    if (option%nflowdof > 0) then
      call GlobalUpdateDenAndSat(realization,option%tran_weight_t1)
    endif

    sum_newton_iterations = 0
    sum_linear_iterations = 0
    icut = 0

    if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT ",47("="))')

    do
     
      if (realization%reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
        call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE)
! The below is set within RTUpdateAuxVarsPatch() when PETSC_TRUE,PETSC_TRUE are passed
!        patch%aux%RT%aux_vars_up_to_date = PETSC_TRUE 
      endif
      if (realization%reaction%use_log_formulation) then
        if (associated(realization%patch%grid%structured_grid) .and. &
            (.not.(realization%patch%grid%structured_grid%p_samr_patch.eq.0))) then
          cur_level => realization%level_list%first
          do 
            if (.not.associated(cur_level)) exit
            cur_patch => cur_level%patch_list%first
            do
              if (.not.associated(cur_patch)) exit
              call GridVecGetArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
              call GridVecGetArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
              log_xx_p(:) = log(xx_p(:))
              call GridVecRestoreArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
              call GridVecRestoreArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
              cur_patch => cur_patch%next
            enddo
            cur_level => cur_level%next
          enddo
        else
          call VecCopy(field%tran_xx,field%tran_log_xx,ierr)
          call VecLog(field%tran_log_xx,ierr)
        endif

        call PetscGetTime(log_start_time, ierr)
        call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_log_xx, ierr)
        call PetscGetTime(log_end_time, ierr)
        stepper%cumulative_solver_time = stepper%cumulative_solver_time + (log_end_time - log_start_time)          
          
        if (associated(realization%patch%grid%structured_grid) .and. &
            (.not.(realization%patch%grid%structured_grid%p_samr_patch.eq.0))) then
          cur_level => realization%level_list%first
          do 
            if (.not.associated(cur_level)) exit
            cur_patch => cur_level%patch_list%first
            do
              if (.not.associated(cur_patch)) exit
              call GridVecGetArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
              call GridVecGetArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
              xx_p(:) = exp(log_xx_p(:))
              call GridVecRestoreArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
              call GridVecRestoreArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
              cur_patch => cur_patch%next
            enddo
            cur_level => cur_level%next
          enddo
        else
          call VecCopy(field%tran_log_xx,field%tran_xx,ierr)
          call VecExp(field%tran_xx,ierr)
        endif
      else
        call PetscGetTime(log_start_time, ierr)
        call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_xx, ierr)
        call PetscGetTime(log_end_time, ierr)
        stepper%cumulative_solver_time = stepper%cumulative_solver_time + (log_end_time - log_start_time)          
      endif

  ! do we really need all this? - geh 
      call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
      call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
      call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

      sum_newton_iterations = sum_newton_iterations + num_newton_iterations
      sum_linear_iterations = sum_linear_iterations + num_linear_iterations
      
      if (snes_reason <= 0) then
        ! The Newton solver diverged, so try reducing the time step.
        icut = icut + 1
        tran_timestep_cut_flag = PETSC_TRUE

        if (icut > stepper%icut_max .or. option%tran_dt<1.d-20) then
          if (option%print_screen_flag) then
            print *,"--> icut_max exceeded: icut/icutmax= ",icut,stepper%icut_max, &
                    "t= ",option%tran_time/realization%output_option%tconv, " dt= ", &
                    option%tran_dt/realization%output_option%tconv
            print *,"Stopping execution!"
          endif
          realization%output_option%plot_name = 'tran_cut_to_failure'
          plot_flag = PETSC_TRUE
          transient_plot_flag = PETSC_FALSE
          call Output(realization,plot_flag,transient_plot_flag)
          failure = PETSC_TRUE
          return
        endif

        option%tran_time = option%tran_time - option%tran_dt
        option%tran_dt = 0.5d0 * option%tran_dt
      
        if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
          &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.4, '' dt= '', &
          &   1pe12.4,l2)')  snes_reason,icut,stepper%icutcum, &
              option%tran_time/realization%output_option%tconv, &
              option%tran_dt/realization%output_option%tconv, &
              tran_timestep_cut_flag

        option%tran_time = option%tran_time + option%tran_dt

        ! recompute weights
        if (option%nflowdof > 0) then
          option%tran_weight_t0 = (option%tran_time-option%tran_dt-start_time)/ &
                                  (end_time-start_time)
          option%tran_weight_t1 = (option%tran_time-start_time)/ &
                                  (end_time-start_time)
        endif
        call RTTimeCut(realization)

      else
        ! increment time steps number
        stepper%steps = stepper%steps + 1      
        ! The Newton solver converged, so we can exit.
        exit
      endif
    enddo

    stepper%newton_cum = stepper%newton_cum + sum_newton_iterations
    stepper%linear_cum = stepper%linear_cum + sum_linear_iterations
    stepper%icutcum = stepper%icutcum + icut

  ! print screen output
    call VecNorm(field%tran_r,NORM_2,fnorm,ierr) 
    call VecNorm(field%tran_r,NORM_INFINITY,inorm,ierr)
    if (option%print_screen_flag) then

      write(*, '(/," TRAN ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
        & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
        & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
        stepper%steps,option%tran_time/realization%output_option%tconv, &
        option%tran_dt/realization%output_option%tconv, &
        realization%output_option%tunit,snes_reason,sum_newton_iterations, &
        stepper%newton_cum,sum_linear_iterations,stepper%linear_cum,icut, &
        stepper%icutcum

      ! the grid pointer is null if we are working with SAMRAI
      if(associated(discretization%grid)) then
         scaled_fnorm = fnorm/discretization%grid%nmax   
      else
         scaled_fnorm = fnorm
      endif
      print *,' --> SNES Linear/Non-Linear Iterations = ', &
               num_linear_iterations,' / ',num_newton_iterations
      write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
    endif

    if (option%print_file_flag) then
      write(option%fid_out, '(" TRAN ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1, &
        & "]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
        & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') stepper%steps, &
        option%tran_time/realization%output_option%tconv, &
        option%tran_dt/realization%output_option%tconv, &
        realization%output_option%tunit,snes_reason,sum_newton_iterations, &
        stepper%newton_cum,sum_linear_iterations,stepper%linear_cum,icut, &
        stepper%icutcum
    endif
    
    call RTMaxChange(realization)
    if (option%print_screen_flag) then
      write(*,'("  --> max chng: dcmx= ",1pe12.4," dc/dt= ",1pe12.4," [mol/s]")') &
        option%dcmax,option%dcmax/option%tran_dt
    endif
    if (option%print_file_flag) then  
      write(option%fid_out,'("  --> max chng: dcmx= ",1pe12.4," dc/dt= ",1pe12.4," [mol/s]")') &
        option%dcmax,option%dcmax/option%tran_dt
    endif

    if (option%print_screen_flag) print *, ""

    ! get out if not simulating flow or time is met
    if (option%nflowdof == 0 .or. &
        option%flow_time - option%tran_time <= time_tol*option%flow_time) exit

    ! taking intermediate time step, need to update solution
    call StepperUpdateTransportSolution(realization)

    ! if dt is smaller than dt_orig/4, try growing it by 0.25d0
    if (icut == 0 .and. option%tran_dt < 0.25d0*dt_orig) then
      option%tran_dt = 1.25d0*option%tran_dt
    endif

    ! compute next time step
    ! can't exceed the flow step
    if (option%tran_time + 1.0001d0*option%tran_dt >= option%flow_time) then
      option%tran_dt = option%flow_time - option%tran_time
      option%tran_time = option%flow_time
    else
      option%tran_time = option%tran_time + option%tran_dt
    endif

  enddo

end subroutine StepperStepTransportDT

! ************************************************************************** !
!
! StepperRunSteadyState: Solves steady state solution for flow and transport
! author: Glenn Hammond
! date: 03/10/09
!
! ************************************************************************** !
subroutine StepperRunSteadyState(realization,flow_stepper,tran_stepper)

  use Realization_module

  use Option_module
  use Output_module, only : Output, OutputInit, OutputVectorTecplot
  use Logging_module
  use Discretization_module

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper

  PetscTruth :: transient_plot_flag
  PetscTruth :: plot_flag
  PetscTruth :: failure
  PetscLogDouble :: start_time, end_time
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  PetscErrorCode :: ierr

  option => realization%option
  output_option => realization%output_option

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  failure = PETSC_FALSE

  call PetscLogStagePush(logging%stage(TS_STAGE),ierr)

  ! print initial condition output if not a restarted sim
  call OutputInit(realization)
  transient_plot_flag = PETSC_FALSE
  plot_flag = PETSC_TRUE
  if (output_option%print_initial) then
    transient_plot_flag = PETSC_FALSE
    call Output(realization,plot_flag,transient_plot_flag)
    if (output_option%print_permeability) then
      if (len_trim(option%group_prefix) > 1) then
        string = 'permeability-' // trim(option%group_prefix) // '.tec'
      else
        string = 'permeability.tec'
      endif
      call DiscretizationLocalToGlobal(realization%discretization, &
                                       realization%field%perm_xx_loc, &
                                       realization%field%work,ONEDOF)
      call OutputVectorTecplot(string,string,realization,realization%field%work)
    endif
    if (output_option%print_porosity) then
      if (len_trim(option%group_prefix) > 1) then
        string = 'porosity-' // trim(option%group_prefix) // '.tec'
      else
        string = 'porosity.tec'
      endif
      call DiscretizationLocalToGlobal(realization%discretization, &
                                       realization%field%porosity_loc, &
                                       realization%field%work,ONEDOF)
      call OutputVectorTecplot(string,string,realization,realization%field%work)
    endif
  endif

  if (OptionPrintToScreen(option)) then
    option%print_screen_flag = PETSC_TRUE
  else
    option%print_screen_flag = PETSC_FALSE
  endif

  if (OptionPrintToFile(option)) then
    option%print_file_flag = PETSC_TRUE
  else
    option%print_file_flag = PETSC_FALSE
  endif

  plot_flag = PETSC_TRUE
    
  if (associated(flow_stepper)) then
    call PetscGetTime(start_time, ierr)
    call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)
    call StepperSolveFlowSteadyState(realization,flow_stepper,failure)
    call PetscLogStagePop(ierr)
    call PetscGetTime(end_time, ierr)
    if (OptionPrintToScreen(option)) then
      write(*, &
         &  '(/,1pe12.4," seconds to solve steady state flow problem",/)') &
        end_time-start_time
    endif
    if (OptionPrintToFile(option)) then
      write(option%fid_out, &
         &  '(/,1pe12.4," seconds to solve steady state flow problem",/)') &
        end_time-start_time
    endif
    if (failure) return ! if flow solve fails, exit
  endif

  if (associated(tran_stepper)) then
    call PetscGetTime(start_time, ierr)
    call PetscLogStagePush(logging%stage(TRAN_STAGE),ierr)
    call StepperSolveTranSteadyState(realization,tran_stepper,failure)
    call PetscLogStagePop(ierr)
    call PetscGetTime(end_time, ierr)
    if (OptionPrintToScreen(option)) then
      write(*, &
         &'(/,1pe12.4," seconds to solve steady state transport problem",/)') &
        end_time-start_time
    endif
    if (OptionPrintToFile(option)) then
      write(option%fid_out, &
         &'(/,1pe12.4," seconds to solve steady state transport problem",/)') &
        end_time-start_time
    endif
    if (failure) return ! if transport solve fails, exit
  endif

  if (output_option%print_initial) then
    output_option%plot_number = 1
    call Output(realization,plot_flag,transient_plot_flag)
  endif
   
  if (option%checkpoint_flag) then
    call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                           ZERO_INTEGER,ZERO_INTEGER,NEG_ONE_INTEGER)  
  endif

  if (OptionPrintToScreen(option)) then
    if (option%nflowdof > 0) then
      write(*,'(/," FLOW newton = ",i8," linear = ",i10)') &
            flow_stepper%newton_cum,flow_stepper%linear_cum
    endif
    if (option%ntrandof > 0) then
      write(*,'(/," TRAN newton = ",i8," linear = ",i10)') &
            tran_stepper%newton_cum,tran_stepper%linear_cum
    endif            
  endif

  if (OptionPrintToFile(option)) then
    if (option%nflowdof > 0) then
      write(option%fid_out,'(/," FLOW newton = ",i8," linear = ",i10)') &
            flow_stepper%steps,flow_stepper%linear_cum
    endif
    if (option%ntrandof > 0) then
      write(option%fid_out,'(/," TRAN newton = ",i8," linear = ",i10)') &
            tran_stepper%newton_cum,tran_stepper%linear_cum
    endif            
  endif

  call PetscLogStagePop(ierr)


end subroutine StepperRunSteadyState
 
! ************************************************************************** !
!
! StepperSolveFlowSteadyState: Solves the steady-state flow equation
! author: Glenn Hammond
! date: 03/10/09
!
! ************************************************************************** !
subroutine StepperSolveFlowSteadyState(realization,stepper,failure)

  use Global_module, only : GlobalUpdateAuxVars
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscTruth :: failure

  PetscErrorCode :: ierr
  PetscInt :: num_newton_iterations
  PetscInt :: num_linear_iterations
  PetscInt :: snes_reason
  PetscReal :: fnorm
  PetscReal :: inorm
  PetscReal :: scaled_fnorm

  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(solver_type), pointer :: solver
  type(discretization_type), pointer :: discretization 

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

  num_newton_iterations = 0
  num_linear_iterations = 0

  call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                  field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tor_loc, &
                                  field%tor_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                  field%icap_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                  field%ithrm_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)
    
  if (option%print_screen_flag) write(*,'(/,2("=")," FLOW (STEADY STATE) ",37("="))')

  call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%flow_xx, ierr)

  call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in FLOW, reason: ', snes_reason
    endif
    failure = PETSC_TRUE
    return
  endif
  
  stepper%newton_cum = num_newton_iterations
  stepper%linear_cum = num_linear_iterations

  if (option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(realization,TIME_T)
    call GlobalUpdateAuxVars(realization,TIME_TpDT)
  endif
    
! print screen output
  call VecNorm(field%flow_r,NORM_2,fnorm,ierr) 
  call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    ! the grid pointer is null if we are working with SAMRAI
    if(associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  
  if (option%print_screen_flag) print *, ""

end subroutine StepperSolveFlowSteadyState

! ************************************************************************** !
!
! StepperSolveTranSteadyState: Steps forward one step in time
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperSolveTranSteadyState(realization,stepper,failure)
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  use Patch_module
  use Level_module
  use Grid_module
    
  use Global_module, only : GlobalUpdateDenAndSat
  use Reactive_Transport_module, only : RTUpdateAuxVars  

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscTruth :: failure
  
  PetscErrorCode :: ierr
  SNESConvergedReason :: snes_reason 
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscReal :: fnorm, scaled_fnorm, inorm
  PetscReal, pointer :: r_p(:), xx_p(:), log_xx_p(:)

  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver
  type(patch_type), pointer :: cur_patch
  type(level_type), pointer :: cur_level

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  call DiscretizationLocalToLocal(discretization,field%porosity_loc,field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tor_loc,field%tor_loc,ONEDOF)

  call GlobalUpdateDenAndSat(realization,1.d0)
  num_newton_iterations = 0
  num_linear_iterations = 0

  if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT (STEADY STATE) ",32("="))')

  if (realization%reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE)
  endif

  if (realization%reaction%use_log_formulation) then
    if (associated(realization%patch%grid%structured_grid) .and. &
        (.not.(realization%patch%grid%structured_grid%p_samr_patch.eq.0))) then
      cur_level => realization%level_list%first
      do 
        if (.not.associated(cur_level)) exit
        cur_patch => cur_level%patch_list%first
        do
          if (.not.associated(cur_patch)) exit
          call GridVecGetArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
          call GridVecGetArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
          log_xx_p(:) = log(xx_p(:))
          call GridVecRestoreArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
          call GridVecRestoreArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
          cur_patch => cur_patch%next
        enddo
        cur_level => cur_level%next
      enddo
    else
      call VecCopy(field%tran_xx,field%tran_log_xx,ierr)
      call VecLog(field%tran_log_xx,ierr)
    endif
      
    call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_log_xx, ierr)
      
    if (associated(realization%patch%grid%structured_grid) .and. &
        (.not.(realization%patch%grid%structured_grid%p_samr_patch.eq.0))) then
      cur_level => realization%level_list%first
      do 
        if (.not.associated(cur_level)) exit
        cur_patch => cur_level%patch_list%first
        do
          if (.not.associated(cur_patch)) exit
          call GridVecGetArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
          call GridVecGetArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
          xx_p(:) = exp(log_xx_p(:))
          call GridVecRestoreArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
          call GridVecRestoreArrayF90(cur_patch%grid,field%tran_log_xx,log_xx_p,ierr)
          cur_patch => cur_patch%next
        enddo
        cur_level => cur_level%next
      enddo
    else
      call VecCopy(field%tran_log_xx,field%tran_xx,ierr)
      call VecExp(field%tran_xx,ierr)
    endif
  else
    call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_xx, ierr)
  endif

  ! do we really need all this? - geh 
  call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in TRAN, reason: ', snes_reason
    endif
    failure = PETSC_TRUE
    return
  endif
      
  stepper%newton_cum = num_newton_iterations
  stepper%linear_cum = num_linear_iterations

  ! print screen output
  call VecNorm(field%tran_r,NORM_2,fnorm,ierr) 
  call VecNorm(field%tran_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    ! the grid pointer is null if we are working with SAMRAI
    if(associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax   
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif

  if (option%print_screen_flag) print *, ""

end subroutine StepperSolveTranSteadyState

! ************************************************************************** !
!
! StepperUpdateSolution: Updates the solution variables
! author: Glenn Hammond
! date: 02/19/08 
!
! ************************************************************************** !
subroutine StepperUpdateSolution(realization)

  use Realization_module
  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  ! update solution variables
  call RealizationUpdate(realization)
  if (realization%option%nflowdof > 0) &
    call StepperUpdateFlowSolution(realization)
  if (realization%option%ntrandof > 0) &
    call StepperUpdateTransportSolution(realization)
    
end subroutine StepperUpdateSolution

! ************************************************************************** !
!
! StepperUpdateFlowSolution: Updates the flow solution variables
! author: Glenn Hammond
! date: 02/19/08 
!
! ************************************************************************** !
subroutine StepperUpdateFlowSolution(realization)
  
  use MPHASE_module, only: MphaseUpdateSolution
  use Immis_module, only: ImmisUpdateSolution
  use Richards_module, only : RichardsUpdateSolution
  use THC_module, only : THCUpdateSolution

  use Realization_module
  use Option_module

  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr
  
  option => realization%option
  
  select case(option%iflowmode)
    case(MPH_MODE)
      call MphaseUpdateSolution(realization)
    case(IMS_MODE)
      call ImmisUpdateSolution(realization)
    case(THC_MODE)
      call THCUpdateSolution(realization)
    case(RICHARDS_MODE)
      call RichardsUpdateSolution(realization)
  end select    

end subroutine StepperUpdateFlowSolution

! ************************************************************************** !
!
! StepperUpdateTransportSolution: Updates the transport solution variables
! author: Glenn Hammond
! date: 02/19/08 
!
! ************************************************************************** !
subroutine StepperUpdateTransportSolution(realization)

  use Realization_module
  use Reactive_Transport_module, only : RTUpdateSolution, RTCalculatePorosity

  implicit none

  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call RTUpdateSolution(realization)
  if (realization%option%calculate_porosity) then
    call RTCalculatePorosity(realization)
  endif

end subroutine StepperUpdateTransportSolution

! ************************************************************************** !
!
! StepperJumpStart: Sets kinetic sorbed concentrations
! author: Glenn Hammond
! date: 08/05/09 
!
! ************************************************************************** !
subroutine StepperJumpStart(realization)

  use Realization_module
  use Reactive_Transport_module, only : RTJumpStartKineticSorption

  implicit none

  type(realization_type) :: realization

  call RTJumpStartKineticSorption(realization)

end subroutine StepperJumpStart

! ************************************************************************** !
!
! StepperUpdateFlowAuxVars: Updates the flow auxilliary variables
! author: Glenn Hammond
! date: 10/11/08 
!
! ************************************************************************** !
subroutine StepperUpdateFlowAuxVars(realization)
  
  use MPHASE_module, only: MphaseUpdateAuxVars
  use Immis_module, only: ImmisUpdateAuxVars
  use Richards_module, only : RichardsUpdateAuxVars
  use THC_module, only : THCUpdateAuxVars

  use Realization_module
  use Option_module

  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr
  
  option => realization%option
  
  select case(option%iflowmode)
    case(IMS_MODE)
      call ImmisUpdateAuxVars(realization)
    case(MPH_MODE)
      call MphaseUpdateAuxVars(realization)
    case(THC_MODE)
      call THCUpdateAuxVars(realization)
    case(RICHARDS_MODE)
      call RichardsUpdateAuxVars(realization)
  end select    

end subroutine StepperUpdateFlowAuxVars

! ************************************************************************** !
!
! StepperCheckpoint: Calls appropriate routines to write a checkpoint file
! author: Glenn Hammond
! date: 03/07/08 
!
! ************************************************************************** !
subroutine StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                             num_const_timesteps,num_newton_iterations,id)

  use Realization_module
  use Checkpoint_module
  use Option_module

  implicit none

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  PetscInt :: num_const_timesteps, num_newton_iterations  
  PetscInt :: id

  type(option_type), pointer :: option
  PetscInt :: flow_steps, flow_newton_cum, flow_icutcum, flow_linear_cum, &
              flow_num_const_timesteps, flow_num_newton_iterations
  PetscInt :: tran_steps, tran_newton_cum, tran_icutcum, tran_linear_cum, &
              tran_num_const_timesteps, tran_num_newton_iterations
  PetscReal :: flow_cumulative_solver_time
  PetscReal :: tran_cumulative_solver_time
  
  option => realization%option

  if (associated(flow_stepper)) then
    flow_steps = flow_stepper%steps
    flow_newton_cum = flow_stepper%newton_cum
    flow_icutcum = flow_stepper%icutcum
    flow_linear_cum = flow_stepper%linear_cum
    flow_num_const_timesteps = num_const_timesteps
    flow_num_newton_iterations = num_newton_iterations
    flow_cumulative_solver_time = flow_stepper%cumulative_solver_time
  endif
  if (associated(tran_stepper)) then
    tran_steps = tran_stepper%steps
    tran_newton_cum = tran_stepper%newton_cum
    tran_icutcum = tran_stepper%icutcum
    tran_linear_cum = tran_stepper%linear_cum
    tran_num_const_timesteps = num_const_timesteps
    tran_num_newton_iterations = num_newton_iterations
    tran_cumulative_solver_time = tran_stepper%cumulative_solver_time
  endif
  
  call Checkpoint(realization, &
                  flow_steps,flow_newton_cum,flow_icutcum,flow_linear_cum, &
                  flow_num_const_timesteps,flow_num_newton_iterations, &
                  flow_cumulative_solver_time, &
                  tran_steps,tran_newton_cum,tran_icutcum,tran_linear_cum, &
                  tran_num_const_timesteps,tran_num_newton_iterations, &
                  tran_cumulative_solver_time, &
                  id)
                      
end subroutine StepperCheckpoint

! ************************************************************************** !
!
! StepperRestart: Calls appropriate routines to read checkpoint file and
!                 restart
! author: Glenn Hammond
! date: 03/07/08 
!
! ************************************************************************** !
subroutine StepperRestart(realization,flow_stepper,tran_stepper, &
                          num_const_timesteps,num_newton_iterations, &
                          flow_read,transport_read,activity_coefs_read)

  use Realization_module
  use Checkpoint_module
  use Option_module

  implicit none

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  PetscInt :: num_const_timesteps, num_newton_iterations
  PetscTruth :: activity_coefs_read
  PetscTruth :: flow_read
  PetscTruth :: transport_read

  type(option_type), pointer :: option
  PetscInt :: flow_steps, flow_newton_cum, flow_icutcum, flow_linear_cum ,&
              flow_num_const_timesteps, flow_num_newton_iterations
  PetscReal :: flow_cum_solver_time
  PetscInt :: tran_steps, tran_newton_cum, tran_icutcum, tran_linear_cum, &
              tran_num_const_timesteps, tran_num_newton_iterations
  PetscReal :: tran_cum_solver_time
  
  option => realization%option

  call Restart(realization, &
               flow_steps,flow_newton_cum,flow_icutcum,flow_linear_cum, &
               flow_num_const_timesteps,flow_num_newton_iterations, &
               flow_cum_solver_time, &
               tran_steps,tran_newton_cum,tran_icutcum,tran_linear_cum, &
               tran_num_const_timesteps,tran_num_newton_iterations, &
               tran_cum_solver_time, &
               flow_read,transport_read,activity_coefs_read)
  if (option%restart_time < -998.d0) then
    option%time = max(option%flow_time,option%tran_time)
    if (associated(flow_stepper) .and. flow_read) then
      flow_stepper%steps = flow_steps
      flow_stepper%newton_cum = flow_newton_cum
      flow_stepper%icutcum = flow_icutcum
      flow_stepper%linear_cum = flow_linear_cum
      flow_stepper%cumulative_solver_time = flow_cum_solver_time
    endif
    if (associated(tran_stepper) .and. transport_read) then
      tran_stepper%steps = tran_steps
      tran_stepper%newton_cum = tran_newton_cum
      tran_stepper%icutcum = tran_icutcum
      tran_stepper%linear_cum = tran_linear_cum
      tran_stepper%cumulative_solver_time = tran_cum_solver_time
    endif
    if (flow_read .and. .not.flow_stepper%run_as_steady_state) then
      num_const_timesteps = flow_num_const_timesteps
      num_newton_iterations = flow_num_newton_iterations
    else if (transport_read) then
      num_const_timesteps = tran_num_const_timesteps
      num_newton_iterations = tran_num_newton_iterations
    else
      num_const_timesteps = 0
      num_newton_iterations = 0
    endif
  else
    option%time = option%restart_time
    option%flow_time = option%restart_time
    option%tran_time = option%restart_time
    if (associated(flow_stepper)) then
      option%flow_dt = flow_stepper%dt_min
      flow_stepper%steps = 0
      flow_stepper%newton_cum = 0
      flow_stepper%icutcum = 0
      flow_stepper%linear_cum = 0
      flow_stepper%cumulative_solver_time = 0.d0
    endif
    if (associated(tran_stepper)) then
      option%tran_dt = tran_stepper%dt_min
      tran_stepper%steps = 0
      tran_stepper%newton_cum = 0
      tran_stepper%icutcum = 0
      tran_stepper%linear_cum = 0
      tran_stepper%cumulative_solver_time = 0.d0
    endif
    num_const_timesteps = 0
    num_newton_iterations = 0
    realization%output_option%plot_number = 0
  endif
    
end subroutine StepperRestart

! ************************************************************************** !
!
! TimestepperPrintInfo: Prints information about time stepper
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine TimestepperPrintInfo(stepper,fid,header,option)

  use Option_module
  
  implicit none
  
  type(stepper_type) :: stepper
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header  
  type(option_type) :: option
  
  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(header)
    write(*,'("max steps:",i8)') stepper%nstepmax
    write(*,'("max const steps:",i4)') stepper%ndtcmx
    write(*,'("max cuts:",i4)') stepper%icut_max
  endif
  if (OptionPrintToFile(option)) then
    write(fid,*) 
    write(fid,'(a)') trim(header)
    write(fid,'("max steps:",i8)') stepper%nstepmax
    write(fid,'("max const steps:",i4)') stepper%ndtcmx
    write(fid,'("max cuts:",i4)') stepper%icut_max
  endif    

end subroutine TimestepperPrintInfo

! ************************************************************************** !
!
! TimestepperDestroy: Deallocates a time stepper
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine TimestepperDestroy(stepper)

  implicit none
  
  type(stepper_type), pointer :: stepper
  
  if (.not.associated(stepper)) return
    
  call SolverDestroy(stepper%solver)
  call ConvergenceContextDestroy(stepper%convergence_context)
  
  if (associated(stepper%tfac)) deallocate(stepper%tfac)
  nullify(stepper%tfac)
  deallocate(stepper)
  nullify(stepper)
  
end subroutine TimestepperDestroy

end module Timestepper_module
