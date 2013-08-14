module Timestepper_BE_class
 
  use Solver_module
  use Convergence_module
  use Timestepper_Base_class
 
  implicit none

  private
  
#include "definitions.h"
 
  type, public, extends(stepper_base_type) :: stepper_BE_type
  
    PetscInt :: num_newton_iterations ! number of Newton iterations in a time step
    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step
    PetscInt :: cumulative_newton_iterations       ! Total number of Newton iterations
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations

    PetscInt :: iaccel        ! Accelerator index
    ! An array of multiplicative factors that specify how to increase time step.
    PetscReal, pointer :: tfac(:)
    PetscInt :: ntfac             ! size of tfac
            
    type(solver_type), pointer :: solver
    type(convergence_context_type), pointer :: convergence_context
  
  contains
    
    procedure, public :: ReadInput => TimestepperBERead
    procedure, public :: Init => TimestepperBEInit
!    procedure, public :: SetTargetTime => TimeStepperBaseSetTargetTime
    procedure, public :: StepDT => TimeStepperBEStepDT
    procedure, public :: UpdateDT => TimeStepperBEUpdateDT
    procedure, public :: FinalizeRun => TimestepperBEFinalizeRun
    procedure, public :: Destroy => TimeStepperBEDestroy
    
  end type stepper_BE_type
  
  public :: TimeStepperBECreate, TimeStepperBEPrintInfo, &
            TimestepperBEInit

contains

! ************************************************************************** !
!
! TimeStepperBECreate: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
function TimeStepperBECreate()

  implicit none
  
  class(stepper_BE_type), pointer :: TimeStepperBECreate
  
  class(stepper_BE_type), pointer :: stepper
  
  allocate(stepper)
  call stepper%Init()
  
  stepper%solver => SolverCreate()
  
  TimeStepperBECreate => stepper
  
end function TimeStepperBECreate

! ************************************************************************** !
!
! TimeStepperBEInit: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimeStepperBEInit(stepper)

  implicit none
  
  class(stepper_BE_type) :: stepper
  
  call TimeStepperBaseInit(stepper)
  
  stepper%num_newton_iterations = 0
  stepper%num_linear_iterations = 0

  stepper%cumulative_newton_iterations = 0
  stepper%cumulative_linear_iterations = 0

  stepper%iaccel = 5
  stepper%ntfac = 13
  allocate(stepper%tfac(13))
  stepper%tfac(1)  = 2.0d0; stepper%tfac(2)  = 2.0d0
  stepper%tfac(3)  = 2.0d0; stepper%tfac(4)  = 2.0d0
  stepper%tfac(5)  = 2.0d0; stepper%tfac(6)  = 1.8d0
  stepper%tfac(7)  = 1.6d0; stepper%tfac(8)  = 1.4d0
  stepper%tfac(9)  = 1.2d0; stepper%tfac(10) = 1.0d0
  stepper%tfac(11) = 1.0d0; stepper%tfac(12) = 1.0d0
  stepper%tfac(13) = 1.0d0
  
  nullify(stepper%solver)
  nullify(stepper%convergence_context)
  
end subroutine TimeStepperBEInit

! ************************************************************************** !
!
! TimeStepperBERead: Reads parameters associated with time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimeStepperBERead(stepper,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none

  class(stepper_BE_type) :: stepper
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','TIMESTEPPER_BE')
    call StringToUpper(keyword)   

    select case(trim(keyword))
  
      case('TS_ACCELERATION')
        call InputReadInt(input,option,stepper%iaccel)
        call InputDefaultMsg(input,option,'iaccel')

      case('DT_FACTOR')
        string='time_step_factor'
        call UtilityReadArray(stepper%tfac,NEG_ONE_INTEGER,string,input, &
            option)
        stepper%ntfac = size(stepper%tfac)

      case default
        call TimeStepperBaseProcessKeyword(stepper,input,option,keyword)
    end select 
  
  enddo  

end subroutine TimeStepperBERead


! ************************************************************************** !
!
! TimeStepperBEUpdateDT: Updates time step
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimeStepperBEUpdateDT(timestepper,process_model)

  use Process_Model_Base_class
  
  implicit none

  class(stepper_BE_type) :: timestepper
  class(pm_base_type) :: process_model
  
  PetscBool :: update_time_step
  
  update_time_step = PETSC_TRUE

  if (timestepper%time_step_cut_flag) then
    timestepper%num_constant_time_steps = 1
  else if (timestepper%num_constant_time_steps > 0) then
    ! otherwise, only increment if the constant time step counter was
    ! initialized to 1
    timestepper%num_constant_time_steps = &
      timestepper%num_constant_time_steps + 1
  endif

  ! num_constant_time_steps = 0: normal time stepping with growing steps
  ! num_constant_time_steps > 0: restriction of constant time steps until
  !                              constant_time_step_threshold is met
  if (timestepper%num_constant_time_steps > &
      timestepper%constant_time_step_threshold) then
    timestepper%num_constant_time_steps = 0
  else if (timestepper%num_constant_time_steps > 0) then
    ! do not increase time step size
    update_time_step = PETSC_FALSE
  endif
    
  if (update_time_step .and. timestepper%iaccel /= 0) then
      
    call process_model%UpdateTimestep(timestepper%dt, &
                                      timestepper%dt_max, &
                                      timestepper%iaccel, &
                                      timestepper%num_newton_iterations, &
                                      timestepper%tfac)
    
  endif

end subroutine TimeStepperBEUpdateDT

! ************************************************************************** !
!
! TimeStepperBEStepDT: Steps forward one step in time
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimeStepperBEStepDT(timestepper,process_model,stop_flag)

  use Process_Model_Base_class
  use Option_module
  use Output_module, only : Output
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"

  class(stepper_BE_type) :: timestepper
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag
  
  SNESConvergedReason :: snes_reason
  PetscInt :: icut
  
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscInt :: num_newton_iterations
  PetscInt :: num_linear_iterations
  PetscInt :: sum_newton_iterations
  PetscInt :: sum_linear_iterations
  character(len=2) :: tunit
  PetscReal :: tconv
  PetscReal :: fnorm, inorm, scaled_fnorm
  PetscBool :: plot_flag, transient_plot_flag
  PetscErrorCode :: ierr
  
  solver => timestepper%solver
  option => process_model%option
  
  write(process_model%option%io_buffer,'(es12.5)') timestepper%dt
  process_model%option%io_buffer = 'StepperStepDT(' // &
    trim(adjustl(process_model%option%io_buffer)) // ')'
  call printMsg(process_model%option)  

  tconv = process_model%output_option%tconv
  tunit = process_model%output_option%tunit
  sum_linear_iterations = 0
  sum_newton_iterations = 0
  icut = 0
  
  option%dt = timestepper%dt
  option%time = timestepper%target_time-timestepper%dt

  call process_model%InitializeTimestep()
    
  do
      
    call process_model%PreSolve()
    
    call PetscTime(log_start_time, ierr)

    call SNESSolve(solver%snes,PETSC_NULL_OBJECT, &
                   process_model%solution_vec,ierr)

    call PetscTime(log_end_time, ierr)

    timestepper%cumulative_solver_time = &
      timestepper%cumulative_solver_time + &
      (log_end_time - log_start_time)

    call SNESGetIterationNumber(solver%snes,num_newton_iterations,ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,ierr)
    call SNESGetConvergedReason(solver%snes,snes_reason,ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
  
    if (snes_reason <= 0 .or. .not. process_model%AcceptSolution()) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      timestepper%time_step_cut_flag = PETSC_TRUE

      if (icut > timestepper%max_time_step_cuts .or. &
          timestepper%dt < 1.d-20) then
        if (option%print_screen_flag) then
          print *,"--> max_time_step_cuts exceeded: icut/icutmax= ",icut, &
                  timestepper%max_time_step_cuts, "t= ", &
                  timestepper%target_time/tconv, &
                  " dt= ", &
                  timestepper%dt/tconv
          print *,"Stopping execution!"
        endif
        process_model%output_option%plot_name = 'flow_cut_to_failure'
        plot_flag = PETSC_TRUE
        transient_plot_flag = PETSC_FALSE
        call Output(process_model%realization_base,plot_flag,transient_plot_flag)
        stop_flag = 2
        return
      endif
 
      timestepper%target_time = timestepper%target_time - timestepper%dt

      timestepper%dt = 0.5d0 * timestepper%dt  
      
      if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
        &   1pe12.5)')  snes_reason,icut,timestepper%cumulative_time_step_cuts, &
            option%time/tconv, &
            timestepper%dt/tconv

      timestepper%target_time = timestepper%target_time + timestepper%dt
      option%dt = timestepper%dt
      call process_model%TimeCut()
  
    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  timestepper%steps = timestepper%steps + 1      
  timestepper%cumulative_newton_iterations = &
    timestepper%cumulative_newton_iterations + sum_newton_iterations
  timestepper%cumulative_linear_iterations = &
    timestepper%cumulative_linear_iterations + sum_linear_iterations
  timestepper%cumulative_time_step_cuts = &
    timestepper%cumulative_time_step_cuts + icut

  timestepper%num_newton_iterations = num_newton_iterations
  timestepper%num_linear_iterations = num_linear_iterations  
  
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(process_model%residual_vec,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    write(*, '(/," Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      timestepper%steps, &
      timestepper%target_time/tconv, &
      timestepper%dt/tconv, &
      tunit,snes_reason,sum_newton_iterations, &
      timestepper%cumulative_newton_iterations,sum_linear_iterations, &
      timestepper%cumulative_linear_iterations,icut, &
      timestepper%cumulative_time_step_cuts

#if 0    
    if (associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif
#endif
    scaled_fnorm = fnorm

    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  if (option%print_file_flag) then
    write(option%fid_out, '(" Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      & " [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      timestepper%steps, &
      timestepper%target_time/tconv, &
      timestepper%dt/tconv, &
      tunit,snes_reason,sum_newton_iterations, &
      timestepper%cumulative_newton_iterations,sum_linear_iterations, &
      timestepper%cumulative_linear_iterations,icut, &
      timestepper%cumulative_time_step_cuts
  endif  
  
  option%time = timestepper%target_time
  call process_model%FinalizeTimestep()
  
  if (option%print_screen_flag) print *, ""  
  
end subroutine TimeStepperBEStepDT

! ************************************************************************** !
!
! TimeStepperBEPrintInfo: Prints information about time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimeStepperBEPrintInfo(stepper,fid,header,option)

  use Option_module
  
  implicit none
  
  class(stepper_BE_type) :: stepper
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  
  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(header)
    write(string,*) stepper%max_time_step
    write(*,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) stepper%constant_time_step_threshold
    write(*,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) stepper%max_time_step_cuts
    write(*,'("max cuts:",x,a)') trim(adjustl(string))
  endif
  if (OptionPrintToFile(option)) then
    write(fid,*) 
    write(fid,'(a)') trim(header)
    write(string,*) stepper%max_time_step
    write(fid,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) stepper%constant_time_step_threshold
    write(fid,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) stepper%max_time_step_cuts
    write(fid,'("max cuts:",x,a)') trim(adjustl(string))
  endif    

end subroutine TimeStepperBEPrintInfo

! ************************************************************************** !
!
! TimestepperCheckpoint: Checkpoints parameters/variables associated with a 
!                        time stepper.
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperCheckpoint(stepper,viewer)

  implicit none

#include "finclude/petscviewer.h"

  class(stepper_BE_type) :: stepper
  PetscViewer :: viewer
    
end subroutine TimestepperCheckpoint

! ************************************************************************** !
!
! TimestepperBEFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
recursive subroutine TimestepperBEFinalizeRun(this,option)

  use Option_module
  
  implicit none
  
  class(stepper_BE_type) :: this
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  call printMsg(option,'TSBE%FinalizeRun()')
  
  if (OptionPrintToScreen(option)) then
    write(*,'(/," TS BE steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            this%steps, &
            this%cumulative_newton_iterations, &
            this%cumulative_linear_iterations, &
            this%cumulative_time_step_cuts
    write(string,'(f12.1)') this%cumulative_solver_time
    write(*,*) 'TS BE SNES time = ' // trim(adjustl(string)) // ' seconds'
  endif
  
end subroutine TimestepperBEFinalizeRun

! ************************************************************************** !
!
! TimeStepperBEStrip: Deallocates members of a time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimeStepperBEStrip(stepper)

  implicit none
  
  class(stepper_BE_type) :: stepper
  
  call SolverDestroy(stepper%solver)
  call ConvergenceContextDestroy(stepper%convergence_context)

  if (associated(stepper%tfac)) deallocate(stepper%tfac)
  nullify(stepper%tfac)
  
end subroutine TimeStepperBEStrip

! ************************************************************************** !
!
! TimeStepperBEDestroy: Deallocates a time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimeStepperBEDestroy(stepper)

  implicit none
  
  class(stepper_BE_type) :: stepper
  
  call TimeStepperBaseStrip(stepper)
  
!  deallocate(stepper)
!  nullify(stepper)
  
end subroutine TimeStepperBEDestroy

end module Timestepper_BE_class
