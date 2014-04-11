module SrcSink_Sandbox_WIPP_Gas_class

! Sandbox srcsink for WIPP gas generation source terms
  
  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_wipp_gas_type
  contains
    procedure, public :: ReadInput => WIPPGasGenerationRead
    procedure, public :: Setup => WIPPGasGenerationSetup
    procedure, public :: Evaluate => WIPPGasGenerationSrcSink
    procedure, public :: Destroy => WIPPGasGenerationDestroy
  end type srcsink_sandbox_wipp_gas_type

  public :: WIPPGasGenerationCreate

contains

! ************************************************************************** !

function WIPPGasGenerationCreate()
  ! 
  ! Allocates WIPP gas generation src/sink object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none
  
  class(srcsink_sandbox_wipp_gas_type), pointer :: WIPPGasGenerationCreate

  allocate(WIPPGasGenerationCreate)
  nullify(WIPPGasGenerationCreate%next)  
      
end function WIPPGasGenerationCreate

! ************************************************************************** !

subroutine WIPPGasGenerationRead(this,input,option)
  ! 
  ! Reads input deck for WIPP gas generation src/sink parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'SRCSINK_SANDBOX,WIPP')
    call StringToUpper(word)   

    select case(trim(word))

!      case('AQUEOUS_SPECIES_NAME')
!        call InputReadWord(input,option,this%aqueous_species_name,PETSC_TRUE)  
!        call InputReadDouble(input,option,this%rate_constant)
!        call InputErrorMsg(input,option,'aqueous species_name', &
!                           'CHEMISTRY,REACTION_SANDBOX,UFD_WP')    
!        ! Read the units
!        call InputReadWord(input,option,word,PETSC_TRUE)
!        if (InputError(input)) then
!          ! If units do not exist, assume default units of 1/s which are the
!          ! standard internal PFLOTRAN units for this rate constant.
!          input%err_buf = 'REACTION_SANDBOX,UFD-WP,RATE CONSTANT UNITS'
!          call InputDefaultMsg(input,option)
!        else              
!          ! If units exist, convert to internal units of 1/s
!          this%rate_constant = this%rate_constant * &
!            UnitsConvertToInternal(word,option)
!        endif
      case default
        option%io_buffer = 'SRCSINK_SANDBOX,WIPP keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine WIPPGasGenerationRead

! ************************************************************************** !

subroutine WIPPGasGenerationSetup(this,option)
  ! 
  ! Sets up the WIPP gas generation src/sink
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module

  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this
  type(option_type) :: option

end subroutine WIPPGasGenerationSetup 

! ************************************************************************** !

subroutine WIPPGasGenerationSrcSink(this,Residual,Jacobian, &
                                    compute_derivative, &
                                    material_auxvar,option)
  ! 
  ! Evaluates src/sink storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this  
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  class(material_auxvar_type) :: material_auxvar
  
  ! residual something
  
  if (compute_derivative) then
    
    ! jacobian something

  endif
  
end subroutine WIPPGasGenerationSrcSink

! ************************************************************************** !

subroutine WIPPGasGenerationDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this  

end subroutine WIPPGasGenerationDestroy

end module SrcSink_Sandbox_WIPP_Gas_class
