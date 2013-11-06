module Timestepper_BE_class
 
  use Solver_module
  use Convergence_module
  use Timestepper_Base_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "finclude/petscsys.h"
 
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
!    procedure, public :: SetTargetTime => TimestepperBaseSetTargetTime
    procedure, public :: StepDT => TimestepperBEStepDT
    procedure, public :: UpdateDT => TimestepperBEUpdateDT
    procedure, public :: Checkpoint => TimestepperBECheckpoint
    procedure, public :: Restart => TimestepperBERestart
    procedure, public :: FinalizeRun => TimestepperBEFinalizeRun
    procedure, public :: Destroy => TimestepperBEDestroy
    
  end type stepper_BE_type
  
  ! For checkpointing
  type, public, extends(stepper_base_header_type) :: stepper_BE_header_type
    integer*8 :: cumulative_newton_iterations
    integer*8 :: cumulative_linear_iterations
    integer*8 :: num_newton_iterations
  end type stepper_BE_header_type
  PetscSizeT, parameter, private :: bagsize = 88 ! 64 (base) + 24 (BE)

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: stepper_BE_header_type
      implicit none
#include "finclude/petscbag.h"      
      PetscBag :: bag
      class(stepper_BE_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData  

  public :: TimestepperBECreate, TimestepperBEPrintInfo, &
            TimestepperBEInit

contains

! ************************************************************************** !
!
! TimestepperBECreate: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
function TimestepperBECreate()

  implicit none
  
  class(stepper_BE_type), pointer :: TimestepperBECreate
  
  class(stepper_BE_type), pointer :: stepper
  
  allocate(stepper)
  call stepper%Init()
  
  stepper%solver => SolverCreate()
  
  TimestepperBECreate => stepper
  
end function TimestepperBECreate

! ************************************************************************** !
!
! TimestepperBEInit: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBEInit(this)

  implicit none
  
  class(stepper_BE_type) :: this
  
  call TimestepperBaseInit(this)
  
  this%num_newton_iterations = 0
  this%num_linear_iterations = 0

  this%cumulative_newton_iterations = 0
  this%cumulative_linear_iterations = 0

  this%iaccel = 5
  this%ntfac = 13
  allocate(this%tfac(13))
  this%tfac(1)  = 2.0d0; this%tfac(2)  = 2.0d0
  this%tfac(3)  = 2.0d0; this%tfac(4)  = 2.0d0
  this%tfac(5)  = 2.0d0; this%tfac(6)  = 1.8d0
  this%tfac(7)  = 1.6d0; this%tfac(8)  = 1.4d0
  this%tfac(9)  = 1.2d0; this%tfac(10) = 1.0d0
  this%tfac(11) = 1.0d0; this%tfac(12) = 1.0d0
  this%tfac(13) = 1.0d0
  
  nullify(this%solver)
  nullify(this%convergence_context)
  
end subroutine TimestepperBEInit

! ************************************************************************** !
!
! TimestepperBERead: Reads parameters associated with time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBERead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none

  class(stepper_BE_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','TIMESTEPPER_BE')
    call StringToUpper(keyword)   

    select case(trim(keyword))
  
      case('TS_ACCELERATION')
        call InputReadInt(input,option,this%iaccel)
        call InputDefaultMsg(input,option,'iaccel')

      case('DT_FACTOR')
        string='time_step_factor'
        call UtilityReadArray(this%tfac,NEG_ONE_INTEGER,string,input, &
            option)
        this%ntfac = size(this%tfac)

      case default
        call TimestepperBaseProcessKeyword(this,input,option,keyword)
    end select 
  
  enddo  

end subroutine TimestepperBERead


! ************************************************************************** !
!
! TimestepperBEUpdateDT: Updates time step
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBEUpdateDT(this,process_model)

  use Process_Model_Base_class
  
  implicit none

  class(stepper_BE_type) :: this
  class(pm_base_type) :: process_model
  
  PetscBool :: update_time_step
  
  update_time_step = PETSC_TRUE

  if (this%time_step_cut_flag) then
    this%num_constant_time_steps = 1
  else if (this%num_constant_time_steps > 0) then
    ! otherwise, only increment if the constant time step counter was
    ! initialized to 1
    this%num_constant_time_steps = &
      this%num_constant_time_steps + 1
  endif

  ! num_constant_time_steps = 0: normal time stepping with growing steps
  ! num_constant_time_steps > 0: restriction of constant time steps until
  !                              constant_time_step_threshold is met
  if (this%num_constant_time_steps > &
      this%constant_time_step_threshold) then
    this%num_constant_time_steps = 0
  else if (this%num_constant_time_steps > 0) then
    ! do not increase time step size
    update_time_step = PETSC_FALSE
  endif
    
  if (update_time_step .and. this%iaccel /= 0) then
      
    call process_model%UpdateTimestep(this%dt, &
                                      this%dt_max, &
                                      this%iaccel, &
                                      this%num_newton_iterations, &
                                      this%tfac)
    
  endif

end subroutine TimestepperBEUpdateDT

! ************************************************************************** !
!
! TimestepperBEStepDT: Steps forward one step in time
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBEStepDT(this,process_model,stop_flag)

  use Process_Model_Base_class
  use Option_module
  use Output_module, only : Output
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"

  class(stepper_BE_type) :: this
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
  
  solver => this%solver
  option => process_model%option
  
!geh: for debugging
!  write(process_model%option%io_buffer,'(es12.5)') this%dt
!  process_model%option%io_buffer = 'StepperStepDT(' // &
!    trim(adjustl(process_model%option%io_buffer)) // ')'
!  call printMsg(process_model%option)  

  tconv = process_model%output_option%tconv
  tunit = process_model%output_option%tunit
  sum_linear_iterations = 0
  sum_newton_iterations = 0
  icut = 0
  
  option%dt = this%dt
  option%time = this%target_time-this%dt

  call process_model%InitializeTimestep()
    
  do
      
    call process_model%PreSolve()
    
    call PetscTime(log_start_time, ierr)

    call SNESSolve(solver%snes,PETSC_NULL_OBJECT, &
                   process_model%solution_vec,ierr)

    call PetscTime(log_end_time, ierr)

    this%cumulative_solver_time = &
      this%cumulative_solver_time + &
      (log_end_time - log_start_time)

    call SNESGetIterationNumber(solver%snes,num_newton_iterations,ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,ierr)
    call SNESGetConvergedReason(solver%snes,snes_reason,ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
  
    if (snes_reason <= 0 .or. .not. process_model%AcceptSolution()) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      this%time_step_cut_flag = PETSC_TRUE

      if (icut > this%max_time_step_cuts .or. &
          this%dt < 1.d-20) then
        if (option%print_screen_flag) then
          print *,"--> max_time_step_cuts exceeded: icut/icutmax= ",icut, &
                  this%max_time_step_cuts, "t= ", &
                  this%target_time/tconv, &
                  " dt= ", &
                  this%dt/tconv
          print *,"Stopping execution!"
        endif
        process_model%output_option%plot_name = 'flow_cut_to_failure'
        plot_flag = PETSC_TRUE
        transient_plot_flag = PETSC_FALSE
        call Output(process_model%realization_base,plot_flag,transient_plot_flag)
        stop_flag = 2
        return
      endif
 
      this%target_time = this%target_time - this%dt

      this%dt = 0.5d0 * this%dt  
      
      if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
        &   1pe12.5)')  snes_reason,icut,this%cumulative_time_step_cuts, &
            option%time/tconv, &
            this%dt/tconv

      this%target_time = this%target_time + this%dt
      option%dt = this%dt
      call process_model%TimeCut()
  
    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  this%steps = this%steps + 1      
  this%cumulative_newton_iterations = &
    this%cumulative_newton_iterations + sum_newton_iterations
  this%cumulative_linear_iterations = &
    this%cumulative_linear_iterations + sum_linear_iterations
  this%cumulative_time_step_cuts = &
    this%cumulative_time_step_cuts + icut

  this%num_newton_iterations = num_newton_iterations
  this%num_linear_iterations = num_linear_iterations  
  
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(process_model%residual_vec,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    write(*, '(/," Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      this%steps, &
      this%target_time/tconv, &
      this%dt/tconv, &
      tunit,snes_reason,sum_newton_iterations, &
      this%cumulative_newton_iterations,sum_linear_iterations, &
      this%cumulative_linear_iterations,icut, &
      this%cumulative_time_step_cuts

    if (associated(process_model%realization_base%discretization%grid)) then
       scaled_fnorm = fnorm/process_model%realization_base% &
                        discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif

    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  if (option%print_file_flag) then
    write(option%fid_out, '(" Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      & " [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      this%steps, &
      this%target_time/tconv, &
      this%dt/tconv, &
      tunit,snes_reason,sum_newton_iterations, &
      this%cumulative_newton_iterations,sum_linear_iterations, &
      this%cumulative_linear_iterations,icut, &
      this%cumulative_time_step_cuts
  endif  
  
  option%time = this%target_time
  call process_model%FinalizeTimestep()
  
  if (option%print_screen_flag) print *, ""  
  
end subroutine TimestepperBEStepDT

! ************************************************************************** !
!
! TimestepperBEPrintInfo: Prints information about time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBEPrintInfo(this,fid,header,option)

  use Option_module
  
  implicit none
  
  class(stepper_BE_type) :: this
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  
  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(header)
    write(string,*) this%max_time_step
    write(*,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) this%constant_time_step_threshold
    write(*,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) this%max_time_step_cuts
    write(*,'("max cuts:",x,a)') trim(adjustl(string))
  endif
  if (OptionPrintToFile(option)) then
    write(fid,*) 
    write(fid,'(a)') trim(header)
    write(string,*) this%max_time_step
    write(fid,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) this%constant_time_step_threshold
    write(fid,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) this%max_time_step_cuts
    write(fid,'("max cuts:",x,a)') trim(adjustl(string))
  endif    

end subroutine TimestepperBEPrintInfo

! ************************************************************************** !
!
! TimestepperBECheckpoint: Checkpoints parameters/variables associated with 
!                          a time stepper.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBECheckpoint(this,viewer,option)

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"
  
  class(stepper_BE_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  class(stepper_BE_header_type), pointer :: header
  PetscBag :: bag
  PetscErrorCode :: ierr

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr)
  call PetscBagGetData(bag,header,ierr)
  call TimestepperBERegisterHeader(this,bag,header)
  call TimestepperBESetHeader(this,bag,header)
  call PetscBagView(bag,viewer,ierr)
  call PetscBagDestroy(bag,ierr)  

end subroutine TimestepperBECheckpoint

! ************************************************************************** !
!
! TimestepperBERegisterHeader: Register header entries.
! author: Glenn Hammond
! date: 07/30/13
!
! ************************************************************************** !
subroutine TimestepperBERegisterHeader(this,bag,header)

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

  class(stepper_BE_type) :: this
  class(stepper_BE_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  ! bagsize = 3 * 8 bytes = 24 bytes
  call PetscBagRegisterInt(bag,header%cumulative_newton_iterations,0, &
                           "cumulative_newton_iterations","",ierr)
  call PetscBagRegisterInt(bag,header%cumulative_linear_iterations,0, &
                           "cumulative_linear_iterations","",ierr)
  call PetscBagRegisterInt(bag,header%num_newton_iterations,0, &
                           "num_newton_iterations","",ierr)

  call TimestepperBaseRegisterHeader(this,bag,header)
  
end subroutine TimestepperBERegisterHeader

! ************************************************************************** !
!
! TimestepperBESetHeader: Sets values in checkpoint header.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBESetHeader(this,bag,header)

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

  class(stepper_BE_type) :: this
  class(stepper_BE_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  header%cumulative_newton_iterations = this%cumulative_newton_iterations
  header%cumulative_linear_iterations = this%cumulative_linear_iterations
  header%num_newton_iterations = this%num_newton_iterations

  call TimestepperBaseSetHeader(this,bag,header)
  
end subroutine TimestepperBESetHeader

! ************************************************************************** !
!
! TimestepperBERestart: Checkpoints parameters/variables associated with 
!                          a time stepper.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBERestart(this,viewer,option)

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

  class(stepper_BE_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  class(stepper_BE_header_type), pointer :: header
  PetscBag :: bag
  PetscErrorCode :: ierr
  
  call PetscBagCreate(option%mycomm,bagsize,bag,ierr)
  call PetscBagGetData(bag,header,ierr)
  call TimestepperBERegisterHeader(this,bag,header)
  call PetscBagLoad(viewer,bag,ierr)
  call TimestepperBEGetHeader(this,header)
  call PetscBagDestroy(bag,ierr)  

end subroutine TimestepperBERestart

! ************************************************************************** !
!
! TimestepperBEGetHeader: Gets values in checkpoint header.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBEGetHeader(this,header)

  use Option_module

  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

  class(stepper_BE_type) :: this
  class(stepper_BE_header_type) :: header
  
  this%cumulative_newton_iterations = header%cumulative_newton_iterations
  this%cumulative_linear_iterations = header%cumulative_linear_iterations
  this%num_newton_iterations = header%num_newton_iterations

  call TimestepperBaseGetHeader(this,header)
  
end subroutine TimestepperBEGetHeader

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
! TimestepperBEStrip: Deallocates members of a time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBEStrip(this)

  implicit none
  
  class(stepper_BE_type) :: this
  
  call SolverDestroy(this%solver)
  call ConvergenceContextDestroy(this%convergence_context)

  if (associated(this%tfac)) deallocate(this%tfac)
  nullify(this%tfac)
  
end subroutine TimestepperBEStrip

! ************************************************************************** !
!
! TimestepperBEDestroy: Deallocates a time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBEDestroy(this)

  implicit none
  
  class(stepper_BE_type) :: this
  
  call TimestepperBaseStrip(this)
  
!  deallocate(this)
!  nullify(this)
  
end subroutine TimestepperBEDestroy

end module Timestepper_BE_class
