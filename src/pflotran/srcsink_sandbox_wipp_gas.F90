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
    procedure, public :: Evaluate2 => WIPPGasGenerationSrcSink
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
  call SSSandboxBaseInit(WIPPGasGenerationCreate)  
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
  PetscBool :: found
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'SRCSINK_SANDBOX,WIPP')
    call StringToUpper(word)   

    ! reads the REGION
    call SSSandboxBaseRead(this,input,option,word,found)
    if (found) cycle
    
    select case(trim(word))
      case default
        option%io_buffer = 'SRCSINK_SANDBOX,WIPP keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo

end subroutine WIPPGasGenerationRead

! ************************************************************************** !

subroutine WIPPGasGenerationSetup(this,region_list,option)
  ! 
  ! Sets up the WIPP gas generation src/sink
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14

  use Option_module
  use Region_module

  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this
  type(region_list_type) :: region_list  
  type(option_type) :: option
  
  call SSSandboxBaseSetup(this,region_list,option)

end subroutine WIPPGasGenerationSetup 

! ************************************************************************** !

subroutine Base_SrcSink1(this,Residual,Jacobian,compute_derivative, &
                         material_auxvar,option)
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  class(material_auxvar_type) :: material_auxvar
      
end subroutine Base_SrcSink1

! ************************************************************************** !

subroutine WIPPGasGenerationSrcSink(this,Residual,Jacobian, &
                                    compute_derivative, &
                                    material_auxvar,aux_real,option)
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
  PetscReal :: aux_real
  
  PetscReal :: water_saturation
  PetscReal :: r_ci, r_mi
  PetscReal :: f_c, f_m
  PetscReal :: s_H2_Fe, s_H2_CH2O
  PetscReal :: r_ch, r_mh
  PetscReal :: q_rc, q_rm
  PetscReal :: q_h2
  
  water_saturation = aux_real
  r_ci = 3.d-8  ! mol Fe/m^3/s
  r_mi = 1.5d-7 ! mol CH2O/m^3/s
  f_c = 1.d-3
  f_m = 0.2d0
  r_ch = f_c * r_ci
  r_mh = f_m * r_mi
  q_rc = r_ci * water_saturation + r_ch * (1.d0-water_saturation)
  q_rm = r_mi * water_saturation + r_mh * (1.d0-water_saturation)
  s_H2_Fe = 1.3081d0
  s_H2_CH2O = 1.1100d0
  q_h2 = s_H2_Fe * q_rc + s_H2_CH2O * q_rm
  
  ! gas production is a negative source (same as injection)
  Residual(TWO_INTEGER) = q_h2
  
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
  
  call SSSandboxBaseDestroy(this)

end subroutine WIPPGasGenerationDestroy

end module SrcSink_Sandbox_WIPP_Gas_class
