module Timestepper_module
 
  use Solver_module
  use Option_module
  use Waypoint_module  
 
  implicit none

  private
 
  type, public :: stepper_type
  
    integer :: flowsteps     ! The number of time-steps taken by the flow code.
    integer :: stepmax       ! The maximum number of time-steps taken by the flow code.
    integer :: nstpmax       ! The maximum number of time-step increments.
    integer :: newton_max    ! Max number of Newton steps for one time step.
    integer :: icut_max      ! Max number of dt cuts for one time step.
    integer :: ndtcmx        ! Steps needed after cutting to increase time step
    integer :: newtcum       ! Total number of Newton steps taken.
    integer :: icutcum       ! Total number of cuts in the timestep taken.    
    integer :: iaccel    
    
    real*8 :: dt_min
    real*8 :: dt_max
!    real*8, pointer :: tstep(:)
!    real*8, pointer :: dtstep(:)    
        
    type(solver_type), pointer :: solver
    type(waypoint_list_type), pointer :: waypoints
    type(waypoint_type), pointer :: cur_waypoint
    real*8, pointer :: steady_eps(:)  ! tolerance for stead state convergence
    
  end type stepper_type
  
  public :: TimestepperCreate, StepperUpdateDT, StepperStepDT, StepperUpdateSolution, &
            TimestepperDestroy, StepperRun, &
            WaypointCreate, WaypointInsertInList
  
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
subroutine StepperRun(solution,stepper,stage)

  use Solution_module
  use Option_module
  use Output_module
  
  implicit none
  
#include "definitions.h"
#include "include/finclude/petsc.h"
#include "include/finclude/petscdef.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

  type(solution_type) :: solution
  type(stepper_type) :: stepper
  PetscInt :: stage(*)
  
  type(option_type), pointer :: option

  real*8 dt_cur
  real*8, pointer :: dxdt(:)
  integer :: idx, ista=0  
  
  logical :: plot_flag, timestep_cut_flag
  integer :: istep, num_timestep_cuts
  integer :: num_newton_iterations
  
  PetscErrorCode :: ierr
  
  option => solution%option

  plot_flag = .false.
  num_timestep_cuts = 0
  timestep_cut_flag = .false.

  call SolutionAddWaypointsToList(solution,stepper%waypoints)
  call WaypointListFillIn(option,stepper%waypoints)
  call WaypointListRemoveExtraWaypnts(option,stepper%waypoints)
  call WaypointConvertTimes(stepper%waypoints,solution%output_option%tconv)
  stepper%cur_waypoint => stepper%waypoints%first

  allocate(dxdt(1:option%ndof))  

  do istep = stepper%flowsteps+1, stepper%stepmax
  
    call StepperStepDT(solution,stepper,plot_flag,timestep_cut_flag, &
                       num_timestep_cuts,num_newton_iterations)
    call StepperUpdateSolution(solution)

#if 0
    ! needs to be modularized
    dt_cur = option%dt 
       
    if(option%imode == THC_MODE)then
      dxdt(1)=option%dpmax/dt_cur
      dxdt(2)=option%dtmpmax/dt_cur
      dxdt(3)=option%dcmax/dt_cur
    endif  
#endif    

    call PetscLogStagePush(stage(2), ierr)
    if (plot_flag) then
      if(option%imode /= OWG_MODE) then
        call Output(solution,istep)
      else
 !       call pflow_var_output(grid,kplt,iplot)
      endif
      plot_flag = .false.
    endif
    call PetscLogStagePop(ierr)
  
    if (.not.timestep_cut_flag) &
      call StepperUpdateDT(stepper,option,num_newton_iterations)

#if 0
    ! still needs implementation
    call PetscLogStagePush(stage(2), ierr)
    if(chkptflag == PETSC_TRUE .and. mod(steps, chkptfreq) == 0) then
      call pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                               num_newton_iterations, steps)
    endif
    call PetscLogStagePop(ierr)
#endif
    
#if 0    
    ! needs to be modularized
    ista=0
    if(option%imode == THC_MODE)then
      do idx = 1, option%ndof
        if(dxdt(idx) < stepper%steady_eps(idx)) ista=ista+1
      enddo 
      
      if(ista >= option%ndof)then
        solution%output_option%plot_number=option%kplot; iplot=1     
      endif
    endif         
#endif
    
    ! if at end of waypoint list (i.e. cur_waypoint = null), we are done!
    if (.not.associated(stepper%cur_waypoint)) exit

  enddo

#if 0
  ! still needs implementation
  if(chkptflag == PETSC_TRUE .and. mod(steps, chkptfreq) /= 0) then
    call pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                             num_newton_iterations, steps)
  endif
#endif  

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

#include "include/finclude/petsc.h"
#include "definitions.h"

  type(stepper_type) :: stepper
  type(option_type) :: option
  integer, intent(in) :: num_newton_iterations
  
  real*8 :: fac,dtt,up,utmp,uc,ut,uus
  
  if (stepper%iaccel == 0) return

  select case(option%imode)
    case(THC_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) fac = 0.33d0
      up = option%dpmxe/(option%dpmax+0.1)
      utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
      uc = option%dcmxe/(option%dcmax+1.d-6)
      ut = min(up,utmp,uc)
      dtt = fac * option%dt * (1.d0 + ut)
    case(TWOPH_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = option%dpmxe/(option%dpmax+0.1)
        utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
        uc = option%dcmxe/(option%dcmax+1.d-6)
        uus=(0.01D0/(option%dsmax+1.d-6))**2
        ut = min(up,utmp,uc)
      endif
      dtt = fac * option%dt * (1.d0 + ut)
    case(MPH_MODE,FLASH_MODE,OWG_MODE,VADOSE_MODE)   
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
    case(TH_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) fac = 0.33d0
      up = option%dpmxe/(option%dpmax+0.1)
      utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
      ut = min(up,utmp)
      dtt = fac * option%dt * (1.d0 + ut)
    case(COND_MODE)
      fac = 0.5d0
      if (num_newton_iterations >= stepper%iaccel) fac = 0.33d0
      ut = option%dtmpmxe/(option%dpmax+1.e-5)
      dtt = fac * option%dt * (1.d0 + ut)
      if (dtt > 2.d0 * option%dt) dtt = 2.d0 * option%dt 
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
subroutine StepperStepDT(solution,stepper,plot_flag,timestep_cut_flag, &
                         num_timestep_cuts,num_newton_iterations)
  
  use translator_mph_module, only : translator_mph_step_maxchange
  use translator_owg_module, only : translator_owg_step_maxchange
  use translator_vad_module, only : translator_vad_step_maxchange
  use translator_flash_module, only : translator_flash_step_maxchange
  use translator_Richards_module, only : translator_ric_step_maxchange
  use TTPHASE_module
  use MPHASE_module
  use Flash_module
  use OWG_module
  use vadose_module
  use Richards_module
  use pflow_solv_module
  
  use Solution_module
  use Grid_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

#include "definitions.h"
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscsnes.h"

  type(solution_type) :: solution
  type(stepper_type) :: stepper

  logical :: plot_flag, timestep_cut_flag
  integer :: num_timestep_cuts,num_newton_iterations
  integer :: ierr
  integer :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  integer :: update_reason, it_linear=0, it_snes
  integer n, nmax_inf
  real*8 m_r2norm, s_r2norm, norm_inf, s_r2norm0, norm_inf0, r2norm
  real*8 :: tsrc
  real*8, pointer :: r_p(:)  

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver

  option => solution%option
  grid => solution%grid
  field => solution%field
  solver => stepper%solver

! real*8, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

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
    if (option%dt > stepper%dt_max) then
      option%dt = stepper%dt_max
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
    plot_flag = 1
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
   
    
    option%iphch=0
    select case(option%imode)
      case(COND_MODE)
        call SNESSolve(solver%snes, PETSC_NULL, field%ttemp, ierr)
      case(TH_MODE)
        call SNESSolve(solver%snes, PETSC_NULL, field%xx, ierr)
      case(THC_MODE)
        call SNESSolve(solver%snes, PETSC_NULL, field%xx, ierr)
      case(TWOPH_MODE)
     ! call  TTPhase_Update(field%xx,grid)
        call SNESSolve(solver%snes, PETSC_NULL, field%xx, ierr)
      case(MPH_MODE)
        if (option%use_ksp == PETSC_TRUE) then
          call pflow_solve(solution,num_newton_iterations, &
                           stepper%newton_max,snes_reason,ierr)
        else 
          call SNESSolve(solver%snes, PETSC_NULL, field%xx, ierr)
        endif
      case(RICHARDS_MODE)
        if (option%use_ksp == PETSC_TRUE) then
          call pflow_solve(solution,num_newton_iterations, &
                           stepper%newton_max,snes_reason,ierr)
        else 
          call SNESSolve(solver%snes, PETSC_NULL, field%xx, ierr)
        endif
      case(FLASH_MODE)
        if (option%use_ksp == PETSC_TRUE) then
          call pflow_solve(solution,num_newton_iterations, &
                           stepper%newton_max,snes_reason,ierr)
        else 
          call SNESSolve(solver%snes, PETSC_NULL, field%xx, ierr)
        endif
      case(VADOSE_MODE)
        if (option%use_ksp == PETSC_TRUE) then
          call pflow_solve(solution,num_newton_iterations, &
                           stepper%newton_max,snes_reason,ierr)
        else 
          call SNESSolve(solver%snes, PETSC_NULL, field%xx, ierr)
        endif
      case(OWG_MODE)
        if (option%use_ksp == PETSC_TRUE) then
          call pflow_solve(solution,num_newton_iterations, &
                           stepper%newton_max,snes_reason,ierr)
        else 
          call SNESSolve(solver%snes,PETSC_NULL, field%xx, ierr)
        endif
      case default
        call SNESSolve(solver%snes,PETSC_NULL, field%ppressure, ierr)
    end select

  ! print *,'pflow_step, finish SNESSolve'
    call MPI_Barrier(PETSC_COMM_WORLD,ierr)
    if (option%use_ksp /= PETSC_TRUE) then
      call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
      it_snes = num_newton_iterations
!     call SNESGetFunctionNorm(solver%snes,r2norm, ierr)
      call VecNorm(field%r, NORM_2, r2norm, ierr) 
      call VecGetArrayF90(field%r, r_p, ierr)
      
      s_r2norm = 0.D0 ; norm_inf = -1.D0 ; nmax_inf =-1
      do n=1, grid%nlmax
         s_r2norm = s_r2norm + r_p(n)* r_p(n)
         if(dabs(r_p(n))> norm_inf)then
            norm_inf = dabs(r_p(n))
            nmax_inf = grid%nL2A(n)
         endif     
      enddo 
     call VecRestoreArrayF90(field%r, r_p, ierr)
     
      if(option%commsize >1)then 
      call MPI_REDUCE(s_r2norm, s_r2norm0,1, MPI_DOUBLE_PRECISION ,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
      call MPI_REDUCE(norm_inf, norm_inf0,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PETSC_COMM_WORLD,ierr)
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
!     if (option%myrank == 0) print *,'SNES R2Norm = ',r2norm, &
!   ' linear Interations = ',it_linear
      call SNESGetConvergedReason(solver%snes, snes_reason, ierr)
    endif
!   parameter (SNES_CONVERGED_ITERATING         =  0)
!   parameter (SNES_CONVERGED_FNORM_ABS         =  2)
!   parameter (SNES_CONVERGED_FNORM_RELATIVE    =  3)
!   parameter (SNES_CONVERGED_PNORM_RELATIVE    =  4)
!   parameter (SNES_CONVERGED_GNORM_ABS         =  5)
!   parameter (SNES_CONVERGED_TR_REDUCTION      =  6)
!   parameter (SNES_CONVERGED_TR_DELTA          =  7)

!   parameter (SNES_DIVERGED_FUNCTION_COUNT     = -2)
!   parameter (SNES_DIVERGED_FNORM_NAN          = -4)
!   parameter (SNES_DIVERGED_MAX_IT             = -5)
!   parameter (SNES_DIVERGED_LS_FAILURE         = -6)
!   parameter (SNES_DIVERGED_TR_REDUCTION       = -7)
!   parameter (SNES_DIVERGED_LOCAL_MIN          = -8)


!******************************************************************
! since calculation on saturation has special requirments, need special treatment
! update_reason
    !call PETScBarrier(PETSC_NULL_OBJECT, ierr)
    
    update_reason = 1
    
    if (snes_reason >= 0) then
    
      select case(option%imode)
#if 0      
        case(TWOPH_MODE)
          call TTPhase_Update_Reason(update_reason,solution)
        case(MPH_MODE)
          call MPhase_Update_Reason(update_reason,solution)
        case(FLASH_MODE)
          call flash_Update_Reason(update_reason,solution)
        case(VADOSE_MODE)
          call Vadose_Update_Reason(update_reason,solution)
#endif          
        case(RICHARDS_MODE)
          update_reason=1
         !call Richards_Update_Reason(update_reason,solution)
#if 0         
        case(OWG_MODE)
          call OWG_Update_Reason(update_reason,solution)
#endif          
      end select   
      if (option%myrank==0) print *,'update_reason: ',update_reason
    endif

!******************************************************************
    
    if (snes_reason <= 0 .or. update_reason <= 0) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      timestep_cut_flag = 1

      if (icut > stepper%icut_max .or. option%dt<1.d-20) then
!       call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)
        if (option%myrank == 0) then
!         t = pflowgrid_get_t(grid)
          print *,"--> icut_max exceeded: icut/icutmax= ",icut,stepper%icut_max, &
                  "t= ",option%time/solution%output_option%tconv, " dt= ", &
                  option%dt/solution%output_option%tconv
          print *,"Stopping execution!"
        endif
        plot_flag = .true.
 !       call pflow_output(grid%ppressure,field%ttemp,grid%conc,grid%phis,grid%porosity, &
 !       grid%perm_xx, grid%perm_yy, grid%perm_zz, &
 !       grid%porosity, grid%sat, grid%vl, &
 !       grid%c_nat,grid%vl_nat,grid%p_nat,option%t_nat,grid%s_nat,grid%phis_nat,grid%por_nat, &
 !       grid%ibrkface, grid%jh2o, grid%nphase, grid%nmax, &
 !       solver%snes, &
 !       option%t, option%dt, option%tconv, stepper%flowsteps, option%rk, &
 !       grid%k1brk,grid%k2brk,grid%j1brk,grid%j2brk,grid%i1brk,grid%i2brk, &
 !       grid%nx,grid%ny,grid%nz,grid%nxy,grid%dx0,grid%dy0,grid%dz0,grid%x,grid%y,grid%z, &
 !       grid%da_nphase_dof,grid%da_1_dof,grid%da_3np_dof,grid%da_ndof,option%ndof, &
 !       kplt,iplot,grid%iprint,grid%ibrkcrv, &
 !       grid%itecplot,grid%write_init,option%myrank)
    
  !      call pflow_output(grid,kplt,iplot)
        ! The above line won't work when scope is restricted properly!
        ! Replace this with a different function call!
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
            option%time/solution%output_option%tconv, &
            option%dt/solution%output_option%tconv,timestep_cut_flag

      if (option%ndof == 1) then
        ! VecCopy(x,y): y=x
        call VecCopy(field%pressure, field%ppressure, ierr)
        call VecCopy(field%temp, field%ttemp, ierr)
      else
        select case(option%imode)
#if 0        
          case(OWG_MODE)
            call pflow_owg_timecut(grid)
          case(MPH_MODE)
            call pflow_mphase_timecut(grid)
          case(FLASH_MODE)
            call pflow_flash_timecut(grid)
#endif            
          case(RICHARDS_MODE)
            call pflow_richards_timecut(solution)
#if 0            
          case(VADOSE_MODE)
            call pflow_vadose_timecut(grid)
          case default
            call VecCopy(grid%h, grid%hh, ierr)
            call VecCopy(grid%yy, field%xx, ierr)
            call VecCopy(grid%density, grid%ddensity, ierr)
#endif            
        end select
        call VecCopy(field%iphas_old_loc, field%iphas_loc, ierr)
      endif

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
      stepper%flowsteps,option%time/solution%output_option%tconv, &
      option%dt/solution%output_option%tconv,solution%output_option%tunit, &
      snes_reason,num_newton_iterations,stepper%newtcum,icut,stepper%icutcum

      if (option%use_ksp /= PETSC_TRUE) then
        print *,' --> SNES Linear/Non-Linear Interations = ',it_linear,it_snes
        print *,' --> SNES Residual: ', r2norm, s_r2norm, m_r2norm, norm_inf 
      endif
       
       
      write(IUNIT2, '(" FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newt= ",i2," [",i6,"]"," cut= ",i2," [",i4, &
      & "]")') stepper%flowsteps,option%time/solution%output_option%tconv, &
      option%dt/solution%output_option%tconv, &
      solution%output_option%tunit, snes_reason,num_newton_iterations,stepper%newtcum,icut,stepper%icutcum
    endif
  endif
  
  ! calculate maxium changes in fields over a time step
#if 0
  if (option%ndof == 1 .and. option%use_cond == PETSC_TRUE) then
  
    call VecWAXPY(field%dp,-1.d0,field%ttemp,field%temp,ierr)
    call VecStrideNorm(field%dp,0,NORM_INFINITY,option%dpmax,ierr)
    if (option%myrank==0) then
      if (mod(stepper%flowsteps,option%imod) == 0) then
        write(*,'("  --> max chng: dTmx= ",1pe12.4)') option%dpmax
        write(IUNIT2,'("  --> max chng: dTmx= ",1pe12.4)') option%dpmax
      endif
    endif
    
  else if (option%ndof == 2 .and. option%use_th == PETSC_TRUE) then
  
    call VecWAXPY(field%dxx,-1.d0,field%xx,field%yy,ierr)
!   call VecAbs(field%dxx,ierr)
!   call VecMax(field%dxx,option%idxxmax,field%dxxmax,ierr)
!   call VecStrideMax(field%dxx,1,PETSC_NULL,field%dxxmax,ierr)
    call VecStrideNorm(field%dxx,0,NORM_INFINITY,option%dpmax,ierr)
    call VecStrideNorm(field%dxx,1,NORM_INFINITY,option%dtmpmax,ierr)
!   call VecStrideNorm(field%dxx,2,NORM_INFINITY,option%dcmax,ierr)
!   n = option%idxxmax/option%ndof+1
!   kz = n/grid%nxy+1
!   jy = n/grid%nx+1 - (kz-1)*grid%ny
!   ix = n - (jy-1)*grid%nx - (kz-1)*grid%nxy
!   j = option%idxxmax - (n-1)*option%ndof
!   if (option%myrank==0) print *,'  --> max chng: ',option%idxxmax,field%dxxmax, &
!     ' @ node ',ix,jy,kz,j,option%ndof,grid%nx,grid%nxy
    if (option%myrank==0) then
      if (mod(stepper%flowsteps,option%imod) == 0 .or. stepper%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4)') option%dpmax,option%dtmpmax
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4)') option%dpmax,option%dtmpmax
      endif
    endif
    
  else if (option%ndof == 3 .and. option%use_thc == PETSC_TRUE) then
  
    call VecWAXPY(field%dxx,-1.d0,field%xx,field%yy,ierr)
!   call VecAbs(field%dxx,ierr)
!   call VecMax(field%dxx,option%idxxmax,field%dxxmax,ierr)
!   call VecStrideMax(field%dxx,1,PETSC_NULL,field%dxxmax,ierr)
    call VecStrideNorm(field%dxx,0,NORM_INFINITY,option%dpmax,ierr)
    call VecStrideNorm(field%dxx,1,NORM_INFINITY,option%dtmpmax,ierr)
    call VecStrideNorm(field%dxx,2,NORM_INFINITY,option%dcmax,ierr)
!   n = option%idxxmax/option%ndof+1
!   kz = n/grid%nxy+1
!   jy = n/grid%nx+1 - (kz-1)*grid%ny
!   ix = n - (jy-1)*grid%nx - (kz-1)*grid%nxy
!   j = option%idxxmax - (n-1)*option%ndof
!   if (option%myrank==0) print *,'  --> max chng: ',option%idxxmax,field%dxxmax, &
!     ' @ node ',ix,jy,kz,j,option%ndof,grid%nx,grid%nxy
    if (option%myrank==0) then
      if (mod(stepper%flowsteps,option%imod) == 0 .or. stepper%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax
      endif
    endif
    
  else if (option%ndof == 4 .and. option%use_2ph == PETSC_TRUE) then
  
    call VecWAXPY(field%dxx,-1.d0,field%xx,field%yy,ierr)
    call VecStrideNorm(field%dxx,0,NORM_INFINITY,option%dpmax,ierr)
    call VecStrideNorm(field%dxx,1,NORM_INFINITY,option%dtmpmax,ierr)
    call VecStrideNorm(field%dxx,2,NORM_INFINITY,option%dcmax,ierr)
    call VecStrideNorm(field%dxx,3,NORM_INFINITY,option%dsmax,ierr)
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

  else if (option%use_mph == PETSC_TRUE) then
     call translator_mph_step_maxchange(grid)
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

  else if (option%use_richards == PETSC_TRUE) then
#endif
  if (option%imode == RICHARDS_MODE) then
     call translator_ric_step_maxchange(solution)
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
#if 0

  else if (option%use_flash == PETSC_TRUE) then
     call translator_flash_step_maxchange(grid)
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


  else if (option%use_vadose == PETSC_TRUE) then
     call translator_vad_step_maxchange(grid)
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

  else if (option%use_owg == PETSC_TRUE) then
    call translator_owg_step_maxchange(grid)
   
      
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
    
  else ! use_liquid
  
    call VecWAXPY(field%dp,-1.d0,field%ppressure,field%pressure,ierr)
!   call VecAbs(field%dp,ierr)
!   call VecMax(field%dp,idpmax,dpmax,ierr)
    call VecStrideNorm(field%dp,0,NORM_INFINITY,option%dpmax,ierr)
    if (option%myrank==0) then
      if (mod(stepper%flowsteps,option%imod) == 0 .or. stepper%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
      endif
    endif
#endif    
  endif

  if (option%myrank == 0 .and. mod(stepper%flowsteps,option%imod) == 0) then
    print *, ""
  endif

end subroutine StepperStepDT

! ************************************************************************** !
!
! StepperUpdateSolution: Updates the solution and solution-dependent variables
! author: 
! date: 
!
! ************************************************************************** !
subroutine StepperUpdateSolution(solution)
  
  use pflow_vector_ops_module
  use TTPHASE_module
  use MPHASE_module
  use Flash_module
  use OWG_module
  use Vadose_module
  use Richards_module, only: pflow_update_richards
  use hydrostat_module, only: recondition_bc

  use Solution_module
  use Option_module
  use Grid_module

  implicit none

#include "include/finclude/petsc.h"  

  type(solution_type) :: solution

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  integer :: ierr
  integer*4 m, n
  real*8, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:), phis_p(:)
  
  option => solution%option
  grid => solution%grid
  field => solution%field
  
! update solution vector and physical properties (VecCopy(x,y): y=x)
  if (option%ndof == 1) then
    call VecCopy(field%ppressure, field%pressure, ierr)
    call VecCopy(field%ttemp, field%temp, ierr)
    call VecCopy(field%ddensity, field%density, ierr)
  else
    if (option%ndof <= 3 .and. &
        option%imode /= MPH_MODE .and. &
        option%imode /= OWG_MODE .and. &
        option%imode /= VADOSE_MODE .and. &
        option%imode /= FLASH_MODE .and. &
        option%imode /= RICHARDS_MODE) then
      call VecCopy(field%xx, field%yy, ierr)
      call VecCopy(field%hh, field%h, ierr)
      call VecCopy(field%ddensity, field%density, ierr)
     ! if (option%use_thc == PETSC_TRUE) call recondition_bc(grid)
    endif    
  endif
 
  select case(option%imode)
#if 0  
    case(OWG_MODE)
      call pflow_update_owg(grid)
    case(MPH_MODE)
      call pflow_update_mphase(grid)
#endif      
    case(RICHARDS_MODE)
      call pflow_update_richards(solution)
#if 0      
    case(FLASH_MODE)
      call pflow_update_flash(grid)
    case(VADOSE_MODE)
      call pflow_update_vadose(grid)
    case(TWOPH_MODE)
      call pflow_update_2phase(grid)  
    case default
      if (option%ndof > 1) then
        call VecGetArrayF90(field%xx, xx_p, ierr)
        call VecGetArrayF90(grid%pressure, press_p, ierr)
        call VecGetArrayF90(field%temp, temp_p, ierr)
        if (option%ndof == 3) call VecGetArrayF90(grid%conc, conc_p, ierr)
        do m = 1, grid%nlmax
          press_p(m) = xx_p(1+(m-1)*option%ndof)
          temp_p(m) = xx_p(2+(m-1)*option%ndof)
          if (option%ndof == 3) conc_p(m) = xx_p(3+(m-1)*option%ndof)
        enddo
        call VecRestoreArrayF90(field%xx, xx_p, ierr)
        call VecRestoreArrayF90(grid%pressure, press_p, ierr)
        call VecRestoreArrayF90(field%temp, temp_p, ierr)
        if (option%ndof == 3) call VecRestoreArrayF90(grid%conc, conc_p, ierr)
#endif        
  end select    

  if (option%run_coupled == PETSC_TRUE) then
    field%xphi_co2 = field%xxphi_co2
    field%den_co2 = field%dden_co2
  endif
  
  !integrate solid volume fraction using explicit finite difference
  if (option%rk > 0.d0) then
    call VecGetArrayF90(field%phis,phis_p,ierr)
    do n = 1, grid%nlmax
      phis_p(n) = phis_p(n) + option%dt * option%vbars * option%rate(n)
      if (phis_p(n) < 0.d0) phis_p(n) = 0.d0
      option%area_var(n) = (phis_p(n)/option%phis0)**option%pwrsrf
      
!     print *,'update: ',n,phis_p(n),option%rate(n),grid%area_var(n)
    enddo
    call VecRestoreArrayF90(field%phis,phis_p,ierr)
  endif
  
  ! update solutoin variables
  call SolutionUpdate(solution)

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

  deallocate(stepper)
  nullify(stepper)
  
end subroutine TimestepperDestroy
  
end module Timestepper_module
