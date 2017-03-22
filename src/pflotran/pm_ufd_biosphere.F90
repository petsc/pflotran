module PM_UFD_Biosphere_class

  use PM_Base_class
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: ERB_base_type
    class(ERB_base_type), pointer :: next
    character(len=MAXWORDLENGTH) :: name
  contains
  end type ERB_base_type
  
  type, public, extends(ERB_base_type) :: ERB_1A_type
  contains
  end type ERB_1A_type
  
  type, public, extends(ERB_base_type) :: ERB_1B_type
  contains
  end type ERB_1B_type

  type, public, extends(pm_base_type) :: pm_ufd_biosphere_type
    class(realization_subsurface_type), pointer :: realization
    class(ERB_base_type), pointer :: ERB_list
  contains
    procedure, public :: PMUFDBSetRealization
    procedure, public :: Setup => PMUFDBSetup
    procedure, public :: Read => PMUFDBRead
    procedure, public :: InitializeRun => PMUFDBInitializeRun
    procedure, public :: InitializeTimestep => PMUFDBInitializeTimestep
    procedure, public :: FinalizeTimestep => PMUFDBFinalizeTimestep
    procedure, public :: Solve => PMUFDBSolve
    procedure, public :: Output => PMUFDBOutput
    procedure, public :: InputRecord => PMUFDBInputRecord
    procedure, public :: Destroy => PMUFDBDestroy
  end type pm_ufd_biosphere_type

  public :: PMUFDBCreate, &
            ERB_1A_Create, &
            ERB_1B_Create
  
contains

! *************************************************************************** !

function PMUFDBCreate()
  !
  ! Creates and initializes the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none
  
  class(pm_ufd_biosphere_type), pointer :: PMUFDBCreate
  
  allocate(PMUFDBCreate)
  nullify(PMUFDBCreate%realization)
  PMUFDBCreate%name = 'ufd biosphere'

  call PMBaseInit(PMUFDBCreate)
  
end function PMUFDBCreate

! *************************************************************************** !

subroutine ERBInit(ERB_model)
  !
  ! Initializes an ERB type model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  class(ERB_base_type) :: ERB_model
  
  ERB_model%name = ''
  nullify(ERB_model%next)

end subroutine ERBInit

! *************************************************************************** !

function ERB_1A_Create()
  !
  ! Creates and initializes an ERB_1A model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  type(ERB_1A_type), pointer :: ERB_1A_Create
  
  allocate(ERB_1A_Create)

  call ERBInit(ERB_1A_Create)
  
end function ERB_1A_Create

! *************************************************************************** !

function ERB_1B_Create()
  !
  ! Creates and initializes an ERB_1B model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  type(ERB_1B_type), pointer :: ERB_1B_Create
  
  allocate(ERB_1B_Create)

  call ERBInit(ERB_1B_Create)
  
end function ERB_1B_Create

! *************************************************************************** !

subroutine PMUFDBSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  use Realization_Subsurface_class

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMUFDBSetRealization

! *************************************************************************** !

subroutine PMUFDBRead(this,input)
  !
  ! Reads input file parameters for the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  character(len=MAXSTRINGLENGTH) :: error_string
  type(ERB_1B_type), pointer :: new_ERB1B
  type(ERB_1A_type), pointer :: new_ERB1A
  class(ERB_base_type), pointer :: cur_ERB
  PetscBool :: added
  
  option => this%option
  input%ierr = 0
  option%io_buffer = 'pflotran card:: UFD_BIOSPHERE'
  call printMsg(option)
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    error_string = 'UFD_BIOSPHERE'
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------------
    !-----------------------------------------
      case('ERB_1A')
    !-----------------------------------------
    !-----------------------------------------
      case('ERB_1B')
        error_string = trim(error_string) // ',ERB_1B'
        allocate(new_ERB1B)
        new_ERB1B => ERB_1B_Create()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_ERB1B%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_ERB1B%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('a')

          !-----------------------------------
            case('b')

          !-----------------------------------    
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages ---------------------
        added = PETSC_FALSE
        if (.not.associated(this%ERB_list)) then
          this%ERB_list => new_ERB1B
        else
          cur_ERB => this%ERB_list
          do
            if (.not.associated(cur_ERB)) exit
            if (.not.associated(cur_ERB%next)) then
              cur_ERB%next => new_ERB1B
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_ERB => cur_ERB%next
          enddo
        endif
        nullify(new_ERB1B)
    !-----------------------------------------
    !-----------------------------------------
      case default
        call InputKeywordUnrecognized(word,'WIPP_SOURCE_SINK',option)
    !-----------------------------------------
    end select  
  enddo
  
end subroutine PMUFDBRead

! *************************************************************************** !

subroutine PMUFDBSetup(this)
  !
  ! 
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBSetup

! ************************************************************************** !

subroutine PMUFDBInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBInitializeRun

! *************************************************************************** !

subroutine PMUFDBInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
  PetscReal :: dt

  dt = this%option%flow_dt
  
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," UFD BIOSPHERE MODEL ",45("="))')
  endif

  
end subroutine PMUFDBInitializeTimestep

! *************************************************************************** !

 subroutine PMUFDBSolve(this,time,ierr)
  ! 
  ! 
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

end subroutine PMUFDBSolve

! ************************************************************************** !

subroutine PMUFDBFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBFinalizeTimestep

! *************************************************************************** !

 subroutine PMUFDBOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBOutput

! *************************************************************************** !

subroutine PMUFDBInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  ! 
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMUFDBInputRecord

! *************************************************************************** !

subroutine PMUFDBDestroy(this)
  ! 
  ! Strips and destroys the UFD Biosphere process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
  call PMUFDBStrip(this)
  
end subroutine PMUFDBDestroy

! ************************************************************************** !

subroutine PMUFDBStrip(this)
  ! 
  ! Strips the UFD Biosphere process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this


end subroutine PMUFDBStrip

! ************************************************************************** !

end module PM_UFD_Biosphere_class