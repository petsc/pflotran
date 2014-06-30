module SrcSink_Sandbox_Base_class
  
  use PFLOTRAN_Constants_module
  use Region_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, abstract, public :: srcsink_sandbox_base_type
    character(len=MAXWORDLENGTH) :: region_name
    type(region_type), pointer :: region
    class(srcsink_sandbox_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => Base_Read
    procedure, public :: Setup => SSSandboxBaseSetup
    procedure, public :: Evaluate => Base_SrcSink
    procedure, public :: Destroy => SSSandboxBaseDestroy    
  end type srcsink_sandbox_base_type
  
  public :: SSSandboxBaseInit, &
            SSSandboxBaseSetup, &
            SSSandboxBaseRead, &
            SSSandboxBaseDestroy
  
contains

! ************************************************************************** !

subroutine SSSandboxBaseInit(this)
    
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
    
  this%region_name = ''
  nullify(this%region)
  nullify(this%next)
  
end subroutine SSSandboxBaseInit

! ************************************************************************** !

subroutine SSSandboxBaseSetup(this,region_list,option)
    
  use Option_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(region_list_type) :: region_list
  type(option_type) :: option
  
  this%region => &
    RegionGetPtrFromList(this%region_name,region_list)
  if (.not.associated(this%region)) then
    option%io_buffer = 'Source/Sink Sandbox region "' // &
                       trim(this%region_name) // &
                         '" not found in list of regions.'
    call printErrMsg(option)
  endif  
  
end subroutine SSSandboxBaseSetup 

! ************************************************************************** !

subroutine SSSandboxBaseRead(this,input,option,keyword,found)
    
  use Option_module
  use Input_Aux_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  
  found = PETSC_TRUE
  select case(trim(keyword))
    case('REGION')
      call InputReadWord(input,option,this%region_name,PETSC_TRUE)
      call InputErrorMsg(input,option,'REGION','SOURCE_SINK_SANDBOX')
    case default
      found = PETSC_FALSE
  end select   
  
end subroutine SSSandboxBaseRead

! ************************************************************************** !

subroutine Base_Read(this,input,option)
    
  use Option_module
  use Input_Aux_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
end subroutine Base_Read   

! ************************************************************************** !

subroutine Base_SrcSink(this,Residual,Jacobian,compute_derivative, &
                        material_auxvar,aux_real,option)
  
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)
      
end subroutine Base_SrcSink

! ************************************************************************** !

subroutine SSSandboxBaseDestroy(this)

  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  
  nullify(this%region)

end subroutine SSSandboxBaseDestroy  

end module SrcSink_Sandbox_Base_class
