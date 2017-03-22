module PM_UFD_Biosphere_class

  use PM_Base_class
  use Region_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: unsupported_rad_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: supported_parent
    character(len=MAXWORDLENGTH) :: element
    PetscReal :: dcf
    PetscReal :: emanation_factor
    PetscReal :: kd
    PetscReal :: kd_supp_parent
    type(unsupported_rad_type), pointer :: next
  end type unsupported_rad_type

  type, public :: ERB_base_type
    class(ERB_base_type), pointer :: next
    character(len=MAXWORDLENGTH) :: name
    type(region_type), pointer :: region
    character(len=MAXWORDLENGTH) :: region_name
    PetscReal :: indv_consumption_rate
    PetscBool :: incl_unsupported_rads
  contains
  end type ERB_base_type
  
  type, public, extends(ERB_base_type) :: ERB_1A_type
  contains
  end type ERB_1A_type
  
  type, public, extends(ERB_base_type) :: ERB_1B_type
    PetscReal :: dilution_factor
  contains
  end type ERB_1B_type

  type, public, extends(pm_base_type) :: pm_ufd_biosphere_type
    class(realization_subsurface_type), pointer :: realization
    class(ERB_base_type), pointer :: ERB_list
    type(unsupported_rad_type), pointer :: unsupported_rad_list
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
            ERB_1B_Create, &
            UnsuppRadCreate
  
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
  ERB_model%region_name = ''
  ERB_model%indv_consumption_rate = UNINITIALIZED_DOUBLE
  ERB_model%incl_unsupported_rads = PETSC_FALSE
  nullify(ERB_model%region)
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
  
  ERB_1B_Create%dilution_factor = UNINITIALIZED_DOUBLE
  
  call ERBInit(ERB_1B_Create)
  
end function ERB_1B_Create

! *************************************************************************** !

function UnsuppRadCreate()
  !
  ! Creates and initializes an unsupported radionuclide type.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  type(unsupported_rad_type), pointer :: UnsuppRadCreate
  
  allocate(UnsuppRadCreate)
  
  UnsuppRadCreate%name = ''
  UnsuppRadCreate%supported_parent = ''
  UnsuppRadCreate%element = ''
  UnsuppRadCreate%dcf = UNINITIALIZED_DOUBLE  ! Sv/Bq
  UnsuppRadCreate%emanation_factor = 1.d0     ! default value
  UnsuppRadCreate%kd = UNINITIALIZED_DOUBLE
  UnsuppRadCreate%kd_supp_parent = UNINITIALIZED_DOUBLE
  nullify(UnsuppRadCreate%next)
  
end function UnsuppRadCreate

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
        error_string = trim(error_string) // ',ERB_1A'
        allocate(new_ERB1A)
        new_ERB1A => ERB_1A_Create()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_ERB1A%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_ERB1A%name)
        call ReadERBmodel(this,input,option,new_ERB1A,error_string)
        ! add new ERB_1A model to ERB_list
        added = PETSC_FALSE
        if (.not.associated(this%ERB_list)) then
          this%ERB_list => new_ERB1A
        else
          cur_ERB => this%ERB_list
          do
            if (.not.associated(cur_ERB)) exit
            if (.not.associated(cur_ERB%next)) then
              cur_ERB%next => new_ERB1A
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_ERB => cur_ERB%next
          enddo
        endif
        nullify(new_ERB1A)
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
        call ReadERBmodel(this,input,option,new_ERB1B,error_string)
        ! add new ERB_1B model to ERB_list
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
      case('UNSUPPORTED_RADIONUCLIDES')
        error_string = trim(error_string) // ',UNSUPPORTED_RADIONUCLIDES'
        call ReadUnsuppRad(this,input,option,error_string)      
    !-----------------------------------------
    !-----------------------------------------
      case default
        call InputKeywordUnrecognized(word,'UFD_BIOSPHERE',option)
    !-----------------------------------------
    end select  
  enddo
  
end subroutine PMUFDBRead

! *************************************************************************** !

subroutine ReadERBmodel(this,input,option,ERB_model,error_string)
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
  class(ERB_base_type) :: ERB_model
  character(len=MAXSTRINGLENGTH) :: error_string
  
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------
      case('REGION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'region assignment',error_string)
        ERB_model%region_name = trim(word)
    !-----------------------------------
      case('DILUTION_FACTOR')
        select type(ERB_model)
          type is(ERB_1A_type)
            option%io_buffer = 'DILUTION_FACTOR cannot be specified for ' &
                               // trim(error_string)
            call printErrMsg(option)
          type is(ERB_1B_type)
            call InputReadDouble(input,option,ERB_model%dilution_factor)
            call InputErrorMsg(input,option,'DILUTION_FACTOR',error_string)
        end select
    !-----------------------------------
      case('INDIVIDUAL_CONSUMPTION_RATE')
          call InputReadDouble(input,option,double)
          call InputErrorMsg(input,option,'INDIVIDUAL_CONSUMPTION_RATE', &
                             error_string)
          call InputReadAndConvertUnits(input,double,'L/day', &
               trim(error_string) // ',INDIVIDUAL_CONSUMPTION_RATE units', &
               option)
          ERB_model%indv_consumption_rate = double
    !-----------------------------------
      case('INCLUDE_UNSUPPORTED_RADS')
        ERB_model%incl_unsupported_rads = PETSC_TRUE
    !-----------------------------------    
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !----------------------------------- 
    end select
  enddo

  ! error messages
  select type(ERB_model)
    type is(ERB_1B_type)
      if (uninitialized(ERB_model%dilution_factor)) then
        option%io_buffer = 'DILUTION_FACTOR must be specified in the ' // &
                           trim(error_string) // ' block.'
        call printErrMsg(option)
      endif
  end select
  
  if (ERB_model%region_name == '') then
    option%io_buffer = 'REGION must be specified in the ' // &
                       trim(error_string) // ' block.'
    call printErrMsg(option)
  endif
  if (uninitialized(ERB_model%indv_consumption_rate)) then
    option%io_buffer = 'INDIVIDUAL_CONSUMPTION_RATE must be specified &
                       &in the ' // trim(error_string) // ' block.'
    call printErrMsg(option)
  endif
  
end subroutine ReadERBmodel

! *************************************************************************** !

subroutine ReadUnsuppRad(this,input,option,error_string)
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
  character(len=MAXSTRINGLENGTH) :: error_string
  
  type(unsupported_rad_type), pointer :: new_unsupp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: added
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------
      case('RADIONUCLIDE')
        error_string = trim(error_string) // ',RADIONUCLIDE'
        allocate(new_unsupp_rad)
        new_unsupp_rad => UnsuppRadCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'radionuclide name',error_string)
        new_unsupp_rad%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_unsupp_rad%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('ELEMENT_KD')
              call InputReadDouble(input,option,new_unsupp_rad%kd)
              call InputErrorMsg(input,option,'ELEMENT_KD',error_string)
          !-----------------------------------
            case('SUPPORTED_PARENT')
              call InputReadWord(input,option,new_unsupp_rad%supported_parent, &
                                 PETSC_TRUE)
              call InputErrorMsg(input,option,'SUPPORTED_PARENT',error_string)
          !-----------------------------------
            case('INGESTION_DOSE_COEF')
              call InputReadDouble(input,option,new_unsupp_rad%dcf)
              call InputErrorMsg(input,option,'INGESTION_DOSE_COEF', &
                                 error_string)
          !-----------------------------------
            case('EMANATION_FACTOR')
              call InputReadDouble(input,option,new_unsupp_rad%emanation_factor)
              call InputErrorMsg(input,option,'EMANATION_FACTOR',error_string)
          !-----------------------------------  
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages
        if (uninitialized(new_unsupp_rad%kd)) then
          option%io_buffer = 'ELEMENT_KD must be specified in the ' // &
                             trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        if (uninitialized(new_unsupp_rad%dcf)) then
          option%io_buffer = 'INGESTION_DOSE_COEF must be specified in the ' &
                             // trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        if (new_unsupp_rad%supported_parent == '') then
          option%io_buffer = 'SUPPORTED_PARENT must be specified in the ' // &
                             trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        ! add new unsupported radionuclide to list
        added = PETSC_FALSE
        if (.not.associated(this%unsupported_rad_list)) then
          this%unsupported_rad_list => new_unsupp_rad
        else
          cur_unsupp_rad => this%unsupported_rad_list
          do
            if (.not.associated(cur_unsupp_rad)) exit
            if (.not.associated(cur_unsupp_rad%next)) then
              cur_unsupp_rad%next => new_unsupp_rad
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_unsupp_rad => cur_unsupp_rad%next
          enddo
        endif
        nullify(new_unsupp_rad)
    !-----------------------------------
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !-----------------------------------
    end select
  enddo
  
end subroutine ReadUnsuppRad

! *************************************************************************** !

subroutine AssociateRegion(this,region_list)
  ! 
  ! Associates the ERB model to its assigned region via the REGION keyword.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !

  use Region_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: cur_region
  class(ERB_base_type), pointer :: cur_ERB
  type(option_type), pointer :: option
  PetscBool :: matched
  
  option => this%option
  
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
      cur_region => region_list%first     
      do
        if (.not.associated(cur_region)) exit
        matched = PETSC_FALSE
        if (StringCompare(trim(cur_region%name), &
                          trim(cur_ERB%region_name))) then
          cur_ERB%region => cur_region
          matched = PETSC_TRUE
        endif
        if (matched) exit
        cur_region => cur_region%next
      enddo      
      if (.not.associated(cur_ERB%region)) then
        option%io_buffer = 'ERB model (' // trim(cur_ERB%name) // ') REGION ' &
                           // trim(cur_ERB%region_name) // ' not found among &
                           &defined regions.'
        call printErrMsg(option)
      endif
    cur_ERB => cur_ERB%next
  enddo
  
end subroutine AssociateRegion

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
  
  ! point the ERB model's regions to the desired regions 
  call AssociateRegion(this,this%realization%patch%region_list)
  
  ! check if all unsupported rads have KDs
  
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
    write(*,'(/,2("=")," UFD BIOSPHERE MODEL ",57("="))')
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