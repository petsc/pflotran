module Timestepper_module
 
  use Solver_module
  use Option_module
  use Waypoint_module 
  use Convergence_module 
 
  implicit none

  private
  
#include "definitions.h"
 
  type, public :: stepper_type
  
    PetscInt :: flowsteps     ! The number of time-steps taken by the flow code.
    PetscInt :: stepmax       ! The maximum number of time-steps taken by the flow code.
    PetscInt :: nstpmax       ! The maximum number of time-step increments.
    PetscInt :: newton_max    ! Max number of Newton steps for one time step.
    PetscInt :: icut_max      ! Max number of dt cuts for one time step.
    PetscInt :: ndtcmx        ! Steps needed after cutting to increase time step
    PetscInt :: newtcum       ! Total number of Newton steps taken.
    PetscInt :: icutcum       ! Total number of cuts in the timestep taken.    
    PetscInt :: iaccel    
    
    PetscReal :: dt_min
    PetscReal :: dt_max
        
    type(solver_type), pointer :: solver
    type(waypoint_list_type), pointer :: waypoints
    type(waypoint_type), pointer :: cur_waypoint
    PetscReal, pointer :: steady_eps(:)  ! tolerance for stead state convergence
    type(convergence_context_type), pointer :: convergence_context
    
  end type stepper_type
  
  public :: TimestepperCreate, TimestepperDestroy, StepperRun
  
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
  stepper%flowsteps = 0
  stepper%stepmax = 0
  stepper%nstpmax = 0
  
  stepper%newton_max = 0
  stepper%icut_max = 0
  stepper%ndtcmx = 5
  stepper%newtcum = 0
  stepper%icutcum = 0    
  stepper%iaccel = 1
  
  stepper%dt_min = 1.d0
  stepper%dt_max = 3.1536d6 ! One-tenth of a year.  
      
  nullify(stepper%solver)
  nullify(stepper%convergence_context)
  nullify(stepper%cur_waypoint)
  
  stepper%solver => SolverCreate()
  stepper%waypoints => WaypointListCreate()
  
  TimeStepperCreate => stepper
  
end function TimestepperCreate 

! ************************************************************************** !
!
! StepperRun: Runs the time step loop
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine StepperRun(realization,stepper)

  use Realization_module
  use Option_module
  use Output_module
  use pflow_checkpoint
  
  implicit none
  
#include "include/finclude/petscdef.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  
  type(option_type), pointer :: option

  logical :: plot_flag, timestep_cut_flag, stop_flag
  PetscInt :: istep, num_timestep_cuts, start_step
  PetscInt :: num_newton_iterations
  
  PetscLogDouble :: stepper_start_time, current_time, average_step_time
  PetscErrorCode :: ierr
  
  option => realization%option

  plot_flag = .false.
  num_timestep_cuts = 0
  timestep_cut_flag = .false.
  stop_flag = .false.

  ! add waypoints associated with boundary conditions, source/sinks etc. to list
  call RealizationAddWaypointsToList(realization,stepper%waypoints)
  ! fill in holes in waypoint data
  call WaypointListFillIn(option,stepper%waypoints)
  call WaypointListRemoveExtraWaypnts(option,stepper%waypoints)
  ! convert times from in put time to seconds
  call WaypointConvertTimes(stepper%waypoints,realization%output_option%tconv)
  stepper%cur_waypoint => stepper%waypoints%first

  if(option%restart_flag == PETSC_TRUE) then
    call pflowGridRestart(realization,stepper%flowsteps,stepper%newtcum, &
                          stepper%icutcum, &
                          timestep_cut_flag,num_timestep_cuts, &
                          num_newton_iterations)
    stepper%cur_waypoint => WaypointSkipToTime(stepper%waypoints,option%time)
    call StepperUpdateSolution(realization)
  endif

  ! print initial condition output if not a restarted sim
  if (option%restart_flag == PETSC_FALSE) then
    call Output(realization)
    call OutputBreakthrough(realization)
  endif
           
  call PetscGetTime(stepper_start_time, ierr)
  start_step = stepper%flowsteps+1
  do istep = start_step, stepper%stepmax
  
    call StepperStepDT(realization,stepper,plot_flag,timestep_cut_flag, &
                       num_timestep_cuts,num_newton_iterations)
    call StepperUpdateSolution(realization)

    if (plot_flag) then
      call Output(realization)
      plot_flag = .false.
    endif
    call OutputBreakthrough(realization)
  
    if (.not.timestep_cut_flag) &
      call StepperUpdateDT(stepper,option,num_newton_iterations)

    ! if a simulation wallclock duration time is set, check to see that the
    ! next time step will not exceed that value.  If it does, print the
    ! checkpoint and exit.
    if (option%wallclock_stop_flag == PETSC_TRUE) then
      call PetscGetTime(current_time, ierr)
      average_step_time = (current_time-stepper_start_time)/ &
                          real(istep-start_step+1) &
                          *2.d0  ! just to be safe, double it
      if (average_step_time + current_time > option%wallclock_stop_time) then
        call printMsg(option,"Wallclock stop time exceeded.  Exiting!!!")
        call printMsg(option,"")
        stop_flag = .true.
      endif
    endif

    if (option%checkpoint_flag == PETSC_TRUE .and. &
        mod(istep,option%checkpoint_frequency) == 0) then
      call pflowGridCheckpoint(realization,stepper%flowsteps,stepper%newtcum, &
                               stepper%icutcum,timestep_cut_flag, &
                               num_timestep_cuts,num_newton_iterations,istep)
    endif
    
   
    ! if at end of waypoint list (i.e. cur_waypoint = null), we are done!
    if (.not.associated(stepper%cur_waypoint) .or. stop_flag) exit

  enddo

  if (option%checkpoint_flag == PETSC_TRUE) then
    call pflowGridCheckpoint(realization,stepper%flowsteps,stepper%newtcum, &
                             stepper%icutcum,timestep_cut_flag, &
                             num_timestep_cuts,num_newton_iterations,NEG_ONE_INTEGER)
  endif

  if (option%myrank == 0) then
    write(*,'(/," PFLOW steps = ",i6," newton = ",i6," cuts = ",i6)') &
          istep-1,stepper%newtcum,stepper%icutcum

    write(IUNIT2,'(/," PFLOW steps = ",i6," newton = ",i6," cuts = ",i6)') &
          istep-1,stepper%newtcum,stepper%icutcum
  endif

end subroutine StepperRun

! ************************************************************************** !
!
! StepperUpdateDT: Updates time step
! author: 
! date: 
!
! ************************************************************************** !
subroutine StepperUpdateDT(stepper,option,num_newton_iterations)

  use Option_module
  
  implicit none

  type(stepper_type) :: stepper
  type(option_type) :: option
  PetscInt, intent(in) :: num_newton_iterations
  
  PetscReal :: fac,dtt,up,utmp,uc,ut,uus
  
  if (stepper%iaccel == 0) return

  select case(option%iflowmode)
    case(THC_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) fac = 0.33d0
      up = option%dpmxe/(option%dpmax+0.1)
      utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
      uc = option%dcmxe/(option%dcmax+1.d-6)
      ut = min(up,utmp,uc)
      dtt = fac * option%dt * (1.d0 + ut)
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
      dtt = fac * option%dt * (1.d0 + ut)
    case(RICHARDS_MODE)
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
      dtt = fac * option%dt * (1.d0 + ut)
    case(RICHARDS_LITE_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = option%dpmxe/(option%dpmax+0.1)
        ut = up
      endif
      dtt = fac * option%dt * (1.d0 + ut)
    case default
      if (num_newton_iterations <= stepper%iaccel .and. &
          num_newton_iterations <= size(option%tfac)) then
        if (num_newton_iterations == 0) then
          dtt = option%tfac(1) * option%dt
        else
          dtt = option%tfac(num_newton_iterations) * option%dt
        endif
      endif
  end select
  
  if (dtt > 2.d0 * option%dt) dtt = 2.d0 * option%dt 
  if (dtt > stepper%dt_max) dtt = stepper%dt_max
  if (dtt>.25d0*option%time .and. option%time>1.d-2) dtt=.25d0*option%time
  option%dt = dtt

end subroutine StepperUpdateDT

! ************************************************************************** !
!
! StepperStepDT: Steps forward one step in time
! author: 
! date: 
!
! ************************************************************************** !
subroutine StepperStepDT(realization,stepper,plot_flag,timestep_cut_flag, &
                         num_timestep_cuts,num_newton_iterations)
  
  use MPHASE_module
  use Richards_Lite_module
  use Richards_Analytical_module
  use THC_module
  use Output_module
  
  use Realization_module
  use Grid_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper

  character(len=MAXSTRINGLENGTH) :: string, string2, string3

  logical :: plot_flag, timestep_cut_flag
  PetscInt :: num_timestep_cuts,num_newton_iterations
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason, it_linear=0, it_snes
  PetscInt :: n, nmax_inf
  PetscReal m_r2norm, s_r2norm, norm_inf, s_r2norm0, norm_inf0, r2norm
  PetscReal :: tsrc
  PetscReal, pointer :: r_p(:)  

  PetscInt, save :: linear_solver_divergence_count = 0

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver

  option => realization%option
  grid => realization%grid
  field => realization%field
  solver => stepper%solver

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  num_newton_iterations = 0
  icut = 0

  ! Perform some global-to-local scatters to update the ghosted vectors.
  ! We have to do this so that the routines for calculating the residual
  ! and the Jacobian will have the ghost points they need.
  ! Note that we don't do the global-to-local scatter for the ppressure 
  ! vector, as that needs to be done within the residual calculation routine
  ! because that routine may get called several times during one Newton step
  ! if a method such as line search is being used.
  call GridLocalToLocal(grid,field%porosity_loc,field%porosity_loc,ONEDOF)
  call GridLocalToLocal(grid,field%tor_loc,field%tor_loc,ONEDOF)
  call GridLocalToLocal(grid,field%icap_loc,field%icap_loc,ONEDOF)
  call GridLocalToLocal(grid,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_loc,ONEDOF)

  option%time = option%time + option%dt
  stepper%flowsteps = stepper%flowsteps + 1

!print *, 'pflow_step:1:',  ntstep, option%dt
 

! If a waypoint calls for a plot or change in src/sinks, adjust time step to match waypoint
  if (option%time + 0.2*option%dt >= stepper%cur_waypoint%time .and. &
      (stepper%cur_waypoint%update_srcs .or. &
       stepper%cur_waypoint%print_output)) then
    option%time = option%time - option%dt
    option%dt = stepper%cur_waypoint%time - option%time
    if (option%dt > stepper%dt_max .and. &
        dabs(option%dt-stepper%dt_max) > 1.d0) then ! 1 sec tolerance to avoid cancellation
      option%dt = stepper%dt_max                    ! error from waypoint%time - option%time
      option%time = option%time + option%dt
    else
      option%time = stepper%cur_waypoint%time
      if (stepper%cur_waypoint%print_output) plot_flag = .true.
      stepper%cur_waypoint => stepper%cur_waypoint%next
      if (associated(stepper%cur_waypoint)) &
        stepper%dt_max = stepper%cur_waypoint%dt_max
    endif
  else if (option%time > stepper%cur_waypoint%time) then
    stepper%cur_waypoint => stepper%cur_waypoint%next
    if (associated(stepper%cur_waypoint)) &
      stepper%dt_max = stepper%cur_waypoint%dt_max
  else if (stepper%flowsteps == stepper%stepmax) then
    plot_flag = .true.
    nullify(stepper%cur_waypoint)
  endif
  
   ! print *, 'pflow_step:3:',  ntstep, option%dt
  if (timestep_cut_flag) then
    num_timestep_cuts = num_timestep_cuts + 1
    if (num_timestep_cuts > stepper%ndtcmx) then
      timestep_cut_flag = .false.
      num_timestep_cuts = 0
    endif
  endif
  
  if (option%myrank == 0) then
    write(*,'(/,60("="))')
  endif
  
  do
   
    
    select case(option%iflowmode)
      case(THC_MODE,MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
        call SNESSolve(solver%snes, PETSC_NULL, field%flow_xx, ierr)
    end select

! do we really need all this? - geh 
    call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
    it_snes = num_newton_iterations
    call VecNorm(field%flow_r, NORM_2, r2norm, ierr) 
    call VecGetArrayF90(field%flow_r, r_p, ierr)
    
    s_r2norm = 0.D0 ; norm_inf = -1.D0 ; nmax_inf =-1
    do n=1, grid%nlmax
       s_r2norm = s_r2norm + r_p(n)* r_p(n)
       if(dabs(r_p(n))> norm_inf)then
          norm_inf = dabs(r_p(n))
          nmax_inf = grid%nL2A(n)
       endif     
    enddo 
   call VecRestoreArrayF90(field%flow_r, r_p, ierr)
   
    if(option%commsize >1)then 
    call MPI_REDUCE(s_r2norm, s_r2norm0,ONE_INTEGER, MPI_DOUBLE_PRECISION ,MPI_SUM,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(norm_inf, norm_inf0,ONE_INTEGER, MPI_DOUBLE_PRECISION, MPI_MAX,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
    if(option%myrank==0) then
      s_r2norm =s_r2norm0
      norm_inf =norm_inf0
    endif
   endif
    s_r2norm = dsqrt(s_r2norm)
    m_r2norm = s_r2norm/grid%nmax   
#if (PETSC_VERSION_RELEASE == 0 || PETSC_VERSION_SUBMINOR == 3)      
    call SNESGetLinearSolveIterations(solver%snes, it_linear, ierr)
#endif      
    call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

    update_reason = 1
    
    if (snes_reason >= 0) then
      select case(option%iflowmode)
        case(MPH_MODE)
          call MPhase_Update_Reason(update_reason,realization)
        case(RICHARDS_MODE)
          update_reason=1
        case(RICHARDS_LITE_MODE)
          update_reason=1
      end select   
      if (option%myrank==0) print *,'update_reason: ',update_reason
    endif

!******************************************************************
    
    if (snes_reason <= 0 .or. update_reason <= 0) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      timestep_cut_flag = 1

      if (icut > stepper%icut_max .or. option%dt<1.d-20) then
        if (option%myrank == 0) then
          print *,"--> icut_max exceeded: icut/icutmax= ",icut,stepper%icut_max, &
                  "t= ",option%time/realization%output_option%tconv, " dt= ", &
                  option%dt/realization%output_option%tconv
          print *,"Stopping execution!"
        endif
        realization%output_option%plot_name = 'cut_to_failure'
        call Output(realization)
        realization%output_option%plot_name = ''
 !       call pflowgrid_destroy(grid)
        call PetscFinalize(ierr)
        stop
      endif

      option%time = option%time - option%dt
      option%dt = 0.5d0 * option%dt
      option%time = option%time + option%dt
    
      if (option%myrank == 0) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.4, '' dt= '', &
        &   1pe12.4,i2)')  snes_reason,icut,stepper%icutcum, &
            option%time/realization%output_option%tconv, &
            option%dt/realization%output_option%tconv,timestep_cut_flag

      select case(option%iflowmode)
        case(THC_MODE)
          call THCTimeCut(realization)
        case(RICHARDS_MODE)
          call RichardsTimeCut(realization)
        case(RICHARDS_LITE_MODE)
          call RichardsLiteTimeCut(realization)
        case(MPH_MODE)
          call pflow_mphase_timecut(realization)
      end select
      call VecCopy(field%iphas_old_loc, field%iphas_loc, ierr)

    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  stepper%newtcum = stepper%newtcum + num_newton_iterations
  stepper%icutcum = stepper%icutcum + icut

! print screen output
  if (option%myrank == 0) then
    if (mod(stepper%flowsteps,option%imod) == 0 .or. stepper%flowsteps == 1) then
      write(*, '(/," FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newt= ",i2," [",i6,"]"," cut= ",i2," [",i4,"]")') &
      stepper%flowsteps,option%time/realization%output_option%tconv, &
      option%dt/realization%output_option%tconv,realization%output_option%tunit, &
      snes_reason,num_newton_iterations,stepper%newtcum,icut,stepper%icutcum

      print *,' --> SNES Linear/Non-Linear Interations = ',it_linear,it_snes
      print *,' --> SNES Residual: ', r2norm, s_r2norm, m_r2norm, norm_inf 
       
      write(IUNIT2, '(" FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newt= ",i2," [",i6,"]"," cut= ",i2," [",i4, &
      & "]")') stepper%flowsteps,option%time/realization%output_option%tconv, &
      option%dt/realization%output_option%tconv, &
      realization%output_option%tunit, snes_reason,num_newton_iterations,stepper%newtcum,icut,stepper%icutcum
    endif
  endif
  
  if (option%iflowmode == RICHARDS_MODE) then
     call RichardsMaxChange(realization)
    if (option%myrank==0) then
      if (mod(stepper%flowsteps,option%imod) == 0 .or. stepper%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax, option%dcmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax
      endif
    endif
  else if (option%iflowmode == RICHARDS_LITE_MODE) then
    call RichardsLiteMaxChange(realization)
    if (option%myrank==0) then
      if (mod(stepper%flowsteps,option%imod) == 0 .or. stepper%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
      endif
    endif
  else if (option%iflowmode == MPH_MODE) then
     call MphaseMaxChange(realization)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (option%myrank==0) then
      if (mod(stepper%flowsteps,option%imod) == 0 .or. stepper%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
      endif
    endif
  endif

  if (option%myrank == 0 .and. mod(stepper%flowsteps,option%imod) == 0) then
    print *, ""
  endif

end subroutine StepperStepDT

! ************************************************************************** !
!
! StepperUpdateSolution: Updates the realization and realization-dependent variables
! author: 
! date: 
!
! ************************************************************************** !
subroutine StepperUpdateSolution(realization)
  
  use MPHASE_module, only: pflow_update_mphase
  use Richards_Lite_module, only : RichardsLiteUpdateSolution
  use Richards_Analytical_module, only : RichardsUpdateSolution

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr
  PetscInt :: m, n
  PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:), phis_p(:)
  
  option => realization%option
  grid => realization%grid
  field => realization%field
  
  select case(option%iflowmode)
    case(MPH_MODE)
      call pflow_update_mphase(realization)
    case(RICHARDS_MODE)
      call RichardsUpdateSolution(realization)
    case(RICHARDS_LITE_MODE)
      call RichardsLiteUpdateSolution(realization)
  end select    
  
  ! update solution variables
  call RealizationUpdate(realization)

end subroutine StepperUpdateSolution

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
