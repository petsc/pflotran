module SrcSink_Sandbox_Base_class
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, abstract, public :: srcsink_sandbox_base_type
    class(srcsink_sandbox_base_type), pointer :: next
  contains
#if 1  
    procedure(Base_Read), public, deferred :: ReadInput
    procedure(Base_Setup), public, deferred :: Setup 
    procedure(Base_SrcSink), public, deferred :: Evaluate
    procedure(Base_Destroy), public, deferred :: Destroy
#else
    procedure, public :: ReadInput => Base_Read
    procedure, public :: Setup => Base_Setup
    procedure, public :: Evaluate => Base_SrcSink
    procedure, public :: Destroy => Base_Destroy    
#endif
  end type srcsink_sandbox_base_type
  
! for some reason cannot use the interfaces when passing in "this"
! with Intel
#if 1 
  abstract interface
  
    subroutine Base_Setup(this,option)
    
      use Option_module
  
      import srcsink_sandbox_base_type
    
      implicit none
  
      class(srcsink_sandbox_base_type) :: this
      type(option_type) :: option
  
    end subroutine Base_Setup 

    subroutine Base_Read(this,input,option)
    
      use Option_module
      use Input_Aux_module
  
      import srcsink_sandbox_base_type
    
      implicit none
  
      class(srcsink_sandbox_base_type) :: this
      type(input_type) :: input
      type(option_type) :: option
  
    end subroutine Base_Read 
    
    subroutine Base_SkipBlock(this,input,option)
    
      use Option_module
      use Input_Aux_module
  
      import srcsink_sandbox_base_type
    
      implicit none
  
      class(srcsink_sandbox_base_type) :: this
      type(input_type) :: input
      type(option_type) :: option
  
    end subroutine Base_SkipBlock 
    
    subroutine Base_SrcSink(this,Residual,Jacobian,compute_derivative, &
                            material_auxvar,option)

      use Option_module
      use Material_Aux_class
  
      import srcsink_sandbox_base_type
    
      implicit none
  
      class(srcsink_sandbox_base_type) :: this
      type(option_type) :: option
      PetscBool :: compute_derivative
      PetscReal :: Residual(option%nflowdof)
      PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
      class(material_auxvar_type) :: material_auxvar
      
    end subroutine
    
    subroutine Base_Destroy(this)

      import srcsink_sandbox_base_type
    
      implicit none
  
      class(srcsink_sandbox_base_type) :: this

    end subroutine Base_Destroy   
    
  end interface

#else

contains

! ************************************************************************** !

  subroutine Base_Setup(this,option)
    
    use Option_module
  
    implicit none
  
    class(srcsink_sandbox_base_type) :: this
    type(option_type) :: option
  
  end subroutine Base_Setup 

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

  subroutine Base_SkipBlock(this,input,option)
    
    use Option_module
    use Input_Aux_module
  
    implicit none
  
    class(srcsink_sandbox_base_type) :: this
    type(input_type) :: input
    type(option_type) :: option
  
  end subroutine Base_SkipBlock   

! ************************************************************************** !

  subroutine Base_SrcSink(this,Residual,Jacobian,compute_derivative, &
                          material_auxvar,option)
    use Option_module
    use Material_Aux_class
  
    implicit none
  
    class(srcsink_sandbox_base_type) :: this
    type(option_type) :: option
    PetscBool :: compute_derivative
    PetscReal :: Residual(option%nflowdof)
    PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
    class(material_auxvar_type) :: material_auxvar
      
  end subroutine

! ************************************************************************** !

  subroutine Base_Destroy(this)

    implicit none
  
    class(srcsink_sandbox_base_type) :: this

  end subroutine Base_Destroy  
#endif

end module SrcSink_Sandbox_Base_class
