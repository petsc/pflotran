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
subroutine TimestepperRead(stepper,fid,option)

  use Option_module
  use Fileio_module
  
  implicit none

  type(stepper_type) :: stepper
  PetscInt :: fid
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string, error_string
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  PetscErrorCode :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)

    if (fiCheckExit(string)) exit  

    call fiReadWord(string,keyword,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'keyword','TIMESTEPPER', ierr)
    call fiWordToUpper(keyword)   
      
    select case(trim(keyword))

      case('NUM_STEPS_AFTER_TS_CUT')
        call fiReadInt(string,stepper%ndtcmx,ierr)
        call fiDefaultMsg(option%myrank,'num_steps_after_ts_cut',ierr)

      case('MAX_STEPS')
        call fiReadInt(string,stepper%nstepmax,ierr)
        call fiDefaultMsg(option%myrank,'nstepmax',ierr)
  
      case('TS_ACCELERATION')
        call fiReadInt(string,stepper%iaccel,ierr)
        call fiDefaultMsg(option%myrank,'iaccel',ierr)

      case('MAX_TS_CUTS')
        call fiReadInt(string,stepper%icut_max,ierr)
        call fiDefaultMsg(option%myrank,'icut_max',ierr)

      case('MAX_PRESSURE_CHANGE')
        call fiReadDouble(string,option%dpmxe,ierr)
        call fiDefaultMsg(option%myrank,'dpmxe',ierr)

      case('MAX_TEMPERATURE_CHANGE')
        call fiReadDouble(string,option%dtmpmxe,ierr)
        call fiDefaultMsg(option%myrank,'dtmpmxe',ierr)
  
      case('MAX_CONCENTRATION_CHANGE')
        call fiReadDouble(string,option%dcmxe,ierr)
        call fiDefaultMsg(option%myrank,'dcmxe',ierr)

      case('MAX_SATURATION_CHANGE')
        call fiReadDouble(string,option%dsmxe,ierr)
        call fiDefaultMsg(option%myrank,'dsmxe',ierr)

      case default
        call printErrMsg(option,'Timestepper option: '//trim(word)// &
                         ' not recognized.')
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
  use Output_module
  use Logging_module  
  use Mass_Balance_module
  
  implicit none
  
#include "finclude/petscdef.h"
#include "finclude/petsclog.h"
#include "finclude/petscsys.h"
#include "finclude/petscviewer.h"

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(stepper_type), pointer :: master_stepper
  
  type(option_type), pointer :: option
  type(waypoint_type), pointer :: prev_waypoint  

  PetscTruth :: plot_flag, timestep_cut_flag, stop_flag
  PetscInt :: istep, start_step
  PetscInt :: num_const_timesteps
  PetscInt :: num_newton_iterations, idum
  
  PetscLogDouble :: stepper_start_time, current_time, average_step_time
  PetscErrorCode :: ierr

  option => realization%option
  
  if (associated(flow_stepper)) then
    master_stepper => flow_stepper
  else
    master_stepper => tran_stepper
  endif

  plot_flag = PETSC_FALSE
  timestep_cut_flag = PETSC_FALSE
  stop_flag = PETSC_FALSE
  num_const_timesteps = 0  

  if (option%restart_flag) then
    call StepperRestart(realization,flow_stepper,tran_stepper, &
                        num_const_timesteps,num_newton_iterations)
    if (associated(flow_stepper)) flow_stepper%cur_waypoint => &
      WaypointSkipToTime(realization%waypoints,option%time)
    if (associated(tran_stepper)) tran_stepper%cur_waypoint => &
      WaypointSkipToTime(realization%waypoints,option%time)

    if (option%nflowdof > 0) then
      call StepperUpdateFlowAuxVars(realization)
    endif

  endif

  if (option%overwrite_restart_transport .and. option%ntrandof > 0) then
    call RealizAssignTransportInitCond(realization)  
  endif
  
  if (option%overwrite_restart_flow_params) then
    call RealizationRevertFlowParameters(realization)
  endif

  call PetscLogStagePop(ierr)
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr)
  
  call StepperUpdateSolution(realization)

  ! print initial condition output if not a restarted sim
  if (realization%output_option%plot_number == 0 .and. master_stepper%nstepmax >= 0) then
    plot_flag = PETSC_TRUE
    call Output(realization,plot_flag)
  endif
           
  call PetscGetTime(stepper_start_time, ierr)
  start_step = master_stepper%steps+1
  do istep = start_step, master_stepper%nstepmax

    if (option%myrank == 0 .and. &
        mod(master_stepper%steps,option%imod) == 0) then
      option%print_flag = PETSC_TRUE
    else
      option%print_flag = PETSC_FALSE
    endif

    prev_waypoint => master_stepper%cur_waypoint
    timestep_cut_flag = PETSC_FALSE
    plot_flag = PETSC_FALSE
    call StepperSetTargetTimes(flow_stepper,tran_stepper,option,plot_flag)

    if (associated(flow_stepper)) then
      call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)
      call StepperStepFlowDT(realization,flow_stepper,timestep_cut_flag, &
                             num_newton_iterations)
      call PetscLogStagePop(ierr)
    endif
    if (associated(tran_stepper)) then
      call PetscLogStagePush(logging%stage(TRAN_STAGE),ierr)
      call StepperStepTransportDT(realization,tran_stepper,timestep_cut_flag, &
                                  idum)
      call PetscLogStagePop(ierr)
      if (.not.associated(flow_stepper)) num_newton_iterations = idum
    endif

    ! update solution variables
    call StepperUpdateSolution(realization)
    
    ! if a time step cut has occured, need to set the below back to original values
    ! if they changed. 
    if (timestep_cut_flag) then
      master_stepper%cur_waypoint => prev_waypoint
      plot_flag = PETSC_FALSE
    endif

    if (option%compute_mass_balance) then
      call MassBalanceUpdate(realization,flow_stepper%solver, &
                             tran_stepper%solver)
    endif
    call Output(realization,plot_flag)
  
    call StepperUpdateDT(flow_stepper,tran_stepper,option,timestep_cut_flag, &
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

  if (option%myrank == 0) then
    if (option%nflowdof > 0) then
      write(*,'(/," FLOW steps = ",i6," newton = ",i6," linear = ",i6, &
            & " cuts = ",i6)') &
!            istep,flow_stepper%newton_cum,flow_stepper%linear_cum, &
            flow_stepper%steps,flow_stepper%newton_cum, &
            flow_stepper%linear_cum,flow_stepper%icutcum

      write(option%fid_out,'(/," FLOW steps = ",i6," newton = ",i6," linear = ",i6, &
            & " cuts = ",i6)') &
!            istep,flow_stepper%newton_cum,flow_stepper%linear_cum, &
            flow_stepper%steps,flow_stepper%newton_cum, &
            flow_stepper%linear_cum,flow_stepper%icutcum
    endif
    if (option%ntrandof > 0) then
      write(*,'(/," TRAN steps = ",i6," newton = ",i6," linear = ",i6, &
            & " cuts = ",i6)') &
!            istep,tran_stepper%newton_cum,tran_stepper%linear_cum, &
            tran_stepper%steps,tran_stepper%newton_cum, &
            tran_stepper%linear_cum,tran_stepper%icutcum

      write(option%fid_out,'(/," TRAN steps = ",i6," newton = ",i6," linear = ",i6, &
            & " cuts = ",i6)') &
!            istep,tran_stepper%newton_cum,tran_stepper%linear_cum, &
            tran_stepper%steps,tran_stepper%newton_cum, &
            tran_stepper%linear_cum,tran_stepper%icutcum
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

  if (num_const_timesteps > 0) num_const_timesteps = num_const_timesteps + 1
  if (timestep_cut_flag) num_const_timesteps = 1

  if (associated(flow_stepper)) then
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
          num_newton_iterations <= size(option%tfac)) then
        if (num_newton_iterations == 0) then
          dtt = option%tfac(1) * dt
        else
          dtt = option%tfac(num_newton_iterations) * dt
        endif
      endif
  end select
  
  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > stepper%dt_max) dtt = stepper%dt_max
  if (dtt>.25d0*time .and. time>5.d2) dtt=.25d0*time
  dt = dtt

  if (associated(flow_stepper)) then
    option%flow_time = time
    option%flow_dt = dt
  else
    option%tran_time = time
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
subroutine StepperSetTargetTimes(flow_stepper,tran_stepper,option,plot_flag)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: flow_stepper, tran_stepper
  type(option_type) :: option
  PetscTruth :: plot_flag
  
  PetscReal :: time
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: steps
  PetscInt :: nstepmax
  type(waypoint_type), pointer :: cur_waypoint

  ! target time will always be dictated by the flow solver, if present
  if (associated(flow_stepper)) then
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
       cur_waypoint%print_output)) then
    time = time - dt
    dt = cur_waypoint%time - time
    if (dt > dt_max .and. dabs(dt-dt_max) > 1.d0) then ! 1 sec tolerance to avoid cancellation
      dt = dt_max                    ! error from waypoint%time - time
      time = time + dt
    else
      time = cur_waypoint%time
      if (cur_waypoint%print_output) plot_flag = PETSC_TRUE
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
                             num_newton_iterations)
  
  use MPHASE_module
  use Richards_module
  use THC_module
  use Output_module
  use Global_module
  
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

  PetscTruth :: timestep_cut_flag
  PetscInt :: num_newton_iterations
  
  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason
  PetscInt :: sum_newton_iterations, sum_linear_iterations, num_linear_iterations
  PetscReal :: fnorm, scaled_fnorm, inorm
  Vec :: global_vec
  PetscTruth :: plot_flag
  
  PetscInt, save :: linear_solver_divergence_count = 0

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
  
  if (option%print_flag) write(*,'(/,2("=")," FLOW ",52("="))')

  if (option%ntrandof > 0) then ! store initial saturations for transport
    call RealizationGetDataset(realization,field%work,LIQUID_DENSITY, &
                               ZERO_INTEGER)
    call DiscretizationGlobalToLocal(realization%discretization, &
                                     field%work,field%work_loc,ONEDOF)
    call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_DENSITY,TIME_T)                                     

    call RealizationGetDataset(realization,field%work,LIQUID_SATURATION, &
                               ZERO_INTEGER)
    call DiscretizationGlobalToLocal(realization%discretization, &
                                     field%work,field%work_loc,ONEDOF)
    call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_SATURATION,TIME_T)                                     
  endif
  
  select case(option%iflowmode)
    case(THC_MODE)
      call THCInitializeTimestep(realization)
    case(RICHARDS_MODE)
      call RichardsInitializeTimestep(realization)
    case(MPH_MODE)
      call MphaseInitializeTimestep(realization)
  end select
  
  sum_newton_iterations = 0
  sum_linear_iterations = 0
    
  do
    
    select case(option%iflowmode)
      case(MPH_MODE,THC_MODE,RICHARDS_MODE)
        call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%flow_xx, ierr)
    end select

! do we really need all this? - geh 
    call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
    call VecNorm(field%flow_r,NORM_2,fnorm,ierr) 
    call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr)
    ! the grid pointer is null if we are working with SAMRAI
    if(associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif

    call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
    update_reason = 1
    
    if (snes_reason >= 0) then
      select case(option%iflowmode)
        case(MPH_MODE)
          call MPhaseUpdateReason(update_reason,realization)
        case(THC_MODE)
          update_reason=1
        case(RICHARDS_MODE)
          update_reason=1
      end select   
      if (option%print_flag) print *,'update_reason: ',update_reason
    endif

!******************************************************************
    
    if (snes_reason <= 0 .or. update_reason <= 0) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      timestep_cut_flag = PETSC_TRUE

      if (icut > stepper%icut_max .or. option%flow_dt<1.d-20) then
        if (option%print_flag) then
          print *,"--> icut_max exceeded: icut/icutmax= ",icut, &
                  stepper%icut_max, "t= ", &
                  option%flow_time/realization%output_option%tconv, " dt= ", &
                  option%flow_dt/realization%output_option%tconv
          print *,"Stopping execution!"
        endif
        realization%output_option%plot_name = 'cut_to_failure'
        plot_flag = PETSC_TRUE
        call Output(realization,plot_flag)
        call PetscFinalize(ierr)
        stop
      endif

      option%flow_time = option%flow_time - option%flow_dt
      option%flow_dt = 0.5d0 * option%flow_dt
    
      if (option%print_flag) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.4, '' dt= '', &
        &   1pe12.4,i3)')  snes_reason,icut,stepper%icutcum, &
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
      end select
      call VecCopy(field%iphas_old_loc, field%iphas_loc, ierr)

    else
      ! increment time steps number
      stepper%steps = stepper%steps + 1      
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  ! for debugging
!  call THCInitializeTimestep(realization)    
!  call SNESComputeFunction(stepper%solver%snes,field%flow_xx,field%flow_r,ierr)
!  call SNESComputeJacobian(stepper%solver%snes,field%flow_xx,stepper%solver%J, &
!                           stepper%solver%Jpre,PETSC_NULL_INTEGER,ierr)

  stepper%newton_cum = stepper%newton_cum + sum_newton_iterations
  stepper%linear_cum = stepper%linear_cum + sum_linear_iterations
  stepper%icutcum = stepper%icutcum + icut

  if (option%ntrandof > 0) then ! store initial saturations for transport
    call RealizationGetDataset(realization,field%work,LIQUID_DENSITY, &
                               ZERO_INTEGER)
    call DiscretizationGlobalToLocal(realization%discretization, &
                                     field%work,field%work_loc,ONEDOF)
    call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_DENSITY,TIME_TpDT)                                     

    call RealizationGetDataset(realization,field%work,LIQUID_SATURATION, &
                               ZERO_INTEGER)
    call DiscretizationGlobalToLocal(realization%discretization, &
                                     field%work,field%work_loc,ONEDOF)
    call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_SATURATION,TIME_TpDT)                                     
  endif
  
! print screen output
  if (option%print_flag) then
    write(*, '(/," FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i2," [",i6,"]", &
      & " linear = ",i5," [",i8,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps,option%flow_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,num_newton_iterations, &
      stepper%newton_cum,num_linear_iterations,stepper%linear_cum,icut, &
      stepper%icutcum

    print *,' --> SNES Linear/Non-Linear Interations = ', &
             num_linear_iterations,num_newton_iterations
    print *,' --> SNES Residual: ', fnorm, scaled_fnorm, inorm 
     
    write(option%fid_out, '(" FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newton = ",i2," [",i6,"]", &
      & " linear = ",i5," [",i8,"]"," cuts = ",i2," [",i4,"]")') stepper%steps, &
      option%flow_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,num_newton_iterations, &
      stepper%newton_cum,num_linear_iterations,stepper%linear_cum,icut, &
      stepper%icutcum
  endif
  
  if (option%iflowmode == THC_MODE) then
     call THCMaxChange(realization)
    if (option%print_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax, option%dcmax
      
      write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dcmax
    endif
  else if (option%iflowmode == RICHARDS_MODE) then
    call RichardsMaxChange(realization)
    if (option%print_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
        
      write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
    endif
  else if (option%iflowmode == MPH_MODE) then
    call MphaseMaxChange(realization)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (option%print_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
        
      write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
        & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
        option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
    endif
  endif

  if (option%print_flag) print *, ""
  
end subroutine StepperStepFlowDT

! ************************************************************************** !
!
! StepperStepTransportDT: Steps forward one step in time
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperStepTransportDT(realization,stepper,timestep_cut_flag, &
                                  num_newton_iterations)
  
  use Reactive_Transport_module
  use Output_module
  
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

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper

  character(len=MAXSTRINGLENGTH) :: string, string2, string3

  PetscTruth :: timestep_cut_flag
  PetscInt :: num_newton_iterations
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason
  PetscInt :: sum_newton_iterations, sum_linear_iterations, num_linear_iterations
  PetscInt :: n, nmax_inf
  PetscReal :: fnorm, scaled_fnorm, inorm
  PetscReal :: start_time, end_time, dt_orig
  PetscTruth :: plot_flag  
  PetscReal, pointer :: r_p(:), xx_p(:), log_xx_p(:)
  PetscReal, parameter :: time_tol = 1.d-10

  PetscInt, save :: linear_solver_divergence_count = 0

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

  if (timestep_cut_flag) then
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

    if (option%print_flag) write(*,'(/,2("=")" TRANSPORT ",47("="))')

    do
     
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
      call VecNorm(field%tran_r,NORM_2,fnorm,ierr) 
      call VecNorm(field%tran_r,NORM_INFINITY,inorm,ierr)
      ! the grid pointer is null if we are working with SAMRAI
      if(associated(discretization%grid)) then
         scaled_fnorm = fnorm/discretization%grid%nmax   
      else
         scaled_fnorm = fnorm
      endif

      call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

      sum_newton_iterations = sum_newton_iterations + num_newton_iterations
      sum_linear_iterations = sum_linear_iterations + num_linear_iterations
      
      if (snes_reason <= 0) then
        ! The Newton solver diverged, so try reducing the time step.
        icut = icut + 1
        timestep_cut_flag = PETSC_TRUE

        if (icut > stepper%icut_max .or. option%tran_dt<1.d-20) then
          if (option%print_flag) then
            print *,"--> icut_max exceeded: icut/icutmax= ",icut,stepper%icut_max, &
                    "t= ",option%tran_time/realization%output_option%tconv, " dt= ", &
                    option%tran_dt/realization%output_option%tconv
            print *,"Stopping execution!"
          endif
          realization%output_option%plot_name = 'cut_to_failure'
          plot_flag = PETSC_TRUE
          call Output(realization,plot_flag)
          call PetscFinalize(ierr)
          stop
        endif

        option%tran_time = option%tran_time - option%tran_dt
        option%tran_dt = 0.5d0 * option%tran_dt
      
        if (option%print_flag) write(*,'('' -> Cut time step: snes='',i3, &
          &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.4, '' dt= '', &
          &   1pe12.4,i3)')  snes_reason,icut,stepper%icutcum, &
              option%tran_time/realization%output_option%tconv, &
              option%tran_dt/realization%output_option%tconv,timestep_cut_flag

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
    if (option%print_flag) then

      write(*, '(/," TRAN ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
        & " snes_conv_reason: ",i4,/,"  newton = ",i2," [",i6,"]", &
        & " linear = ",i5," [",i8,"]"," cuts = ",i2," [",i4,"]")') &
        stepper%steps,option%tran_time/realization%output_option%tconv, &
        option%tran_dt/realization%output_option%tconv, &
        realization%output_option%tunit,snes_reason,num_newton_iterations, &
        stepper%newton_cum,num_linear_iterations,stepper%linear_cum,icut, &
        stepper%icutcum

      print *,' --> SNES Linear/Non-Linear Interations = ', &
               num_linear_iterations,num_newton_iterations
      print *,' --> SNES Residual: ', fnorm, scaled_fnorm, inorm 
       
      write(option%fid_out, '(" TRAN ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1, &
        & "]"," snes_conv_reason: ",i4,/,"  newton = ",i2," [",i6,"]", &
        & " linear = ",i5," [",i8,"]"," cuts = ",i2," [",i4,"]")') stepper%steps, &
        option%tran_time/realization%output_option%tconv, &
        option%tran_dt/realization%output_option%tconv, &
        realization%output_option%tunit,snes_reason,num_newton_iterations, &
        stepper%newton_cum,num_linear_iterations,stepper%linear_cum,icut, &
        stepper%icutcum
    endif
    
    call RTMaxChange(realization)
    if (option%print_flag) then
      write(*,'("  --> max chng: dcmx= ",1pe12.4)') option%dcmax
        
      write(option%fid_out,'("  --> max chng: dcmx= ",1pe12.4)') option%dcmax
    endif

    if (option%print_flag) print *, ""

    ! get out if not simulating flow or time is met
    if (option%nflowdof == 0 .or. &
        option%flow_time - option%tran_time <= time_tol*option%flow_time) exit

    ! if dt is smaller than dt_orig/4, try growing it by 0.25d0
    if (option%tran_dt < 0.25d0*dt_orig) then
      option%tran_dt = 1.25d0*option%tran_dt
    endif

    ! compute next time step
    if (option%tran_time + 1.2d0*option%tran_dt >= option%flow_dt) then
      option%tran_dt = option%flow_time - option%tran_time
      option%tran_time = option%flow_time
    else
      option%tran_time = option%tran_time + option%tran_dt
    endif
    
  enddo

end subroutine StepperStepTransportDT

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
  use Reactive_Transport_module, only : RTUpdateSolution

  implicit none

  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call RTUpdateSolution(realization)

end subroutine StepperUpdateTransportSolution

! ************************************************************************** !
!
! StepperUpdateFlowAuxVars: Updates the flow auxilliary variables
! author: Glenn Hammond
! date: 10/11/08 
!
! ************************************************************************** !
subroutine StepperUpdateFlowAuxVars(realization)
  
  use MPHASE_module, only: MphaseUpdateAuxVars
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
  
  option => realization%option

  if (associated(flow_stepper)) then
    flow_steps = flow_stepper%steps
    flow_newton_cum = flow_stepper%newton_cum
    flow_icutcum = flow_stepper%icutcum
    flow_linear_cum = flow_stepper%linear_cum
    flow_num_const_timesteps = num_const_timesteps
    flow_num_newton_iterations = num_newton_iterations
  endif
  if (associated(tran_stepper)) then
    tran_steps = tran_stepper%steps
    tran_newton_cum = tran_stepper%newton_cum
    tran_icutcum = tran_stepper%icutcum
    tran_linear_cum = tran_stepper%linear_cum
    tran_num_const_timesteps = num_const_timesteps
    tran_num_newton_iterations = num_newton_iterations
  endif
  
  call Checkpoint(realization, &
                  flow_steps,flow_newton_cum,flow_icutcum,flow_linear_cum, &
                  flow_num_const_timesteps,flow_num_newton_iterations, &
                  tran_steps,tran_newton_cum,tran_icutcum,tran_linear_cum, &
                  tran_num_const_timesteps,tran_num_newton_iterations, &
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
                          num_const_timesteps,num_newton_iterations)

  use Realization_module
  use Checkpoint_module
  use Option_module

  implicit none

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  PetscInt :: num_const_timesteps, num_newton_iterations

  type(option_type), pointer :: option
  PetscInt :: flow_steps, flow_newton_cum, flow_icutcum, flow_linear_cum ,&
              flow_num_const_timesteps, flow_num_newton_iterations
  PetscInt :: tran_steps, tran_newton_cum, tran_icutcum, tran_linear_cum, &
              tran_num_const_timesteps, tran_num_newton_iterations
  
  option => realization%option

  call Restart(realization, &
               flow_steps,flow_newton_cum,flow_icutcum,flow_linear_cum, &
               flow_num_const_timesteps,flow_num_newton_iterations, &
               tran_steps,tran_newton_cum,tran_icutcum,tran_linear_cum, &
               tran_num_const_timesteps,tran_num_newton_iterations)
  if (option%restart_time < -998.d0) then
    option%time = max(option%flow_time,option%tran_time)
    if (associated(flow_stepper)) then
      flow_stepper%steps = flow_steps
      flow_stepper%newton_cum = flow_newton_cum
      flow_stepper%icutcum = flow_icutcum
      flow_stepper%linear_cum = flow_linear_cum
    endif
    if (associated(tran_stepper)) then
      tran_stepper%steps = tran_steps
      tran_stepper%newton_cum = tran_newton_cum
      tran_stepper%icutcum = tran_icutcum
      tran_stepper%linear_cum = tran_linear_cum
    endif
    num_const_timesteps = flow_num_const_timesteps
    num_newton_iterations = flow_num_newton_iterations
!    num_const_timesteps = tran_num_const_timesteps
!    num_newton_iterations = tran_num_newton_iterations
  else
    option%time = option%restart_time
    option%flow_time = option%restart_time
    option%flow_dt = flow_stepper%dt_min
    option%tran_time = option%restart_time
    option%tran_dt = tran_stepper%dt_min
    if (associated(flow_stepper)) then
      flow_stepper%steps = 0
      flow_stepper%newton_cum = 0
      flow_stepper%icutcum = 0
      flow_stepper%linear_cum = 0
    endif
    if (associated(tran_stepper)) then
      tran_stepper%steps = 0
      tran_stepper%newton_cum = 0
      tran_stepper%icutcum = 0
      tran_stepper%linear_cum = 0
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
  
  character(len=MAXSTRINGLENGTH) :: string

  if (option%myrank == 0) then
    write(*,*) 
    write(fid,*) 
    write(*,'(a)') trim(header)
    write(fid,'(a)') trim(header)
    write(*,'("max steps:",i6)') stepper%nstepmax
    write(fid,'("max steps:",i6)') stepper%nstepmax
    write(*,'("max const steps:",i4)') stepper%ndtcmx
    write(fid,'("max const steps:",i4)') stepper%ndtcmx
    write(*,'("max cuts:",i4)') stepper%icut_max
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

  deallocate(stepper)
  nullify(stepper)
  
end subroutine TimestepperDestroy
  
end module Timestepper_module
