module cf90interface_module

public
  
  type, public :: f90ptrwrap
     real*8, pointer, dimension(:) :: f90ptr
  end type f90ptrwrap
  
end module cf90interface_module
