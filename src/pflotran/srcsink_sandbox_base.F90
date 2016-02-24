module SrcSink_Sandbox_Base_class
  
  use PFLOTRAN_Constants_module
  use Region_module
  use Geometry_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, abstract, public :: srcsink_sandbox_base_type
    character(len=MAXWORDLENGTH) :: region_name
    type(region_type), pointer :: region
    PetscInt :: local_cell_id
    type(point3d_type) :: coordinate    
    PetscBool :: mass_balance
    PetscReal, pointer :: instantaneous_mass_rate(:)
    PetscReal, pointer :: cumulative_mass(:)
    class(srcsink_sandbox_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => SSSandboxBaseRead
    procedure, public :: Setup => SSSandboxBaseSetup
    procedure, public :: Update => SSSandboxBaseUpdate
    procedure, public :: Evaluate => SSSandboxBaseEvaluate
    procedure, public :: Destroy => SSSandboxBaseDestroy    
  end type srcsink_sandbox_base_type
  
  public :: SSSandboxBaseInit, &
            SSSandboxBaseSetup, &
            SSSandboxBaseRead, &
            SSSandboxBaseSelectCase, &
            SSSandboxBaseDestroy
  
contains

! ************************************************************************** !

subroutine SSSandboxBaseInit(this)
    
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
    
  this%region_name = ''
  this%coordinate%x = UNINITIALIZED_DOUBLE
  this%coordinate%y = UNINITIALIZED_DOUBLE
  this%coordinate%z = UNINITIALIZED_DOUBLE
  this%local_cell_id = UNINITIALIZED_INTEGER
  this%mass_balance = PETSC_FALSE
  nullify(this%region)
  nullify(this%instantaneous_mass_rate)
  nullify(this%cumulative_mass)
  nullify(this%next)
  
end subroutine SSSandboxBaseInit

! ************************************************************************** !

subroutine SSSandboxBaseSetup(this,region_list,grid,option)
    
  use Option_module
  use Grid_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(region_list_type) :: region_list
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscInt :: local_id

  if (len_trim(this%region_name) > 0) then
    this%region => &
      RegionGetPtrFromList(this%region_name,region_list)
    if (.not.associated(this%region)) then
      option%io_buffer = 'Source/Sink Sandbox region "' // &
                         trim(this%region_name) // &
                           '" not found in list of regions.'
      call printErrMsg(option)
    endif
  else if (Initialized(this%coordinate%x)) then
    call GridGetLocalIDFromCoordinate(grid,this%coordinate,option,local_id)
    if (local_id > 0) then
      this%local_cell_id = local_id
    endif
  else
    option%io_buffer = 'Source/sink in SSSandbox not associate with the &
      &domain thorugh either a REGION or COORDINATE.'
    call printErrMsg(option)
  endif  

  ! perhaps we need a global reduction here to catch the lack of a connection.
  
end subroutine SSSandboxBaseSetup 

! ************************************************************************** !

subroutine SSSandboxBaseRead(this,input,option)
    
  use Option_module
  use Input_Aux_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
end subroutine SSSandboxBaseRead  

! ************************************************************************** !

subroutine SSSandboxBaseSelectCase(this,input,option,keyword,found)
    
  use Option_module
  use Input_Aux_module
  use Geometry_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SOURCE_SINK_SANDBOX'
  
  found = PETSC_TRUE
  select case(trim(keyword))
    case('REGION')
      call InputReadWord(input,option,this%region_name,PETSC_TRUE)
      call InputErrorMsg(input,option,'REGION',error_string)
    case('COORDINATE')
      call GeometryReadCoordinate(input,option,this%coordinate,error_string)
    case('MASS_BALANCE')
      this%mass_balance = PETSC_TRUE
    case default
      found = PETSC_FALSE
  end select   
  
end subroutine SSSandboxBaseSelectCase

! ************************************************************************** !

subroutine SSSandboxBaseUpdate(this,time,option)
    
  use Option_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  PetscReal :: time
  type(option_type) :: option
  
  if (associated(this%cumulative_mass)) then
    this%cumulative_mass(:) = this%cumulative_mass(:) + &
      option%flow_dt*this%instantaneous_mass_rate(:)
  endif
  
end subroutine SSSandboxBaseUpdate   

! ************************************************************************** !

subroutine SSSandboxBaseEvaluate(this,Residual,Jacobian,compute_derivative, &
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
      
end subroutine SSSandboxBaseEvaluate

! ************************************************************************** !

subroutine SSSandboxBaseDestroy(this)

  use Utility_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  
  nullify(this%region)
  call DeallocateArray(this%instantaneous_mass_rate)
  call DeallocateArray(this%cumulative_mass)

end subroutine SSSandboxBaseDestroy  

end module SrcSink_Sandbox_Base_class
