module test_module

  implicit none
  
  private
  
  type, public :: auxvar_type
  contains
    procedure, public :: TestProcedure => Procedure1
  end type auxvar_type
  
  type, public :: container_type
    class(auxvar_type), pointer :: aux_vars(:)
  end type container_type  
  
contains

subroutine Procedure1(aux_var)

  implicit none
  
  class(auxvar_type) :: aux_var

end subroutine Procedure1

subroutine test(container)

  implicit none
  
  type(container_type) :: container
  
  class(auxvar_type), pointer :: aux_vars(:)
  
  aux_vars => container%aux_vars

end subroutine test  

end module test_module

program test

  use test_module
  
  implicit none
  
  type(container_type) :: container
  
  class(auxvar_type), pointer :: aux_vars(:)
  
  aux_vars => container%aux_vars
  
  print *, 'Hello World!'
  
end program test

