module shape_class
  
  implicit none
  
  private

  type, public :: shape
  
  ! empty instance variables
  
  contains
  
    procedure, public :: area => calc_area_fn
    procedure, public :: perimeter => calc_perimeter_fn
    procedure, public :: to_string => to_string_fn

  end type shape
  
contains

real function calc_area_fn(this)
  
  implicit none
  
  class(shape) :: this
  
  calc_area_fn = 0.
  
end function calc_area_fn

real function calc_perimeter_fn(this)
  
  implicit none
  
  class(shape) :: this
  
  calc_perimeter_fn = 0.
  
end function calc_perimeter_fn

character(len=50) function to_string_fn(this)
  
  implicit none
  
  class(shape) :: this
  
  to_string_fn = ''
  
end function to_string_fn

end module shape_class

!**********************************************************************!
module circle_class

  use shape_class

  implicit none
  
  private

  type, public, extends(shape) :: circle
  
    real :: r = 0.
  
  contains
  
    procedure, public :: initialize => initialize_sub
    procedure, public :: area => get_area_fn
    procedure, public :: perimeter => get_perimeter_fn
    procedure, public :: to_string => to_string_fn
  
  end type circle
  
  real, parameter :: PI = 3.141593
  
contains

subroutine initialize_sub(this,r)

  implicit none
  
  class(circle) :: this
  real :: r
  
  this%r = r

end subroutine initialize_sub

real function get_area_fn(this)

  implicit none

  class(circle) :: this
  
  get_area_fn = PI * this%r**2.

end function get_area_fn

real function get_perimeter_fn(this)

  implicit none

  class(circle) :: this
  
  get_perimeter_fn = 2. * PI * this%r

end function get_perimeter_fn

character(len=50) function to_string_fn(this)
  
  implicit none
  
  class(circle) :: this
  
  write(to_string_fn,'(a,f6.2)') 'Circle of radius ', this%r
  
end function to_string_fn

end module circle_class

!**********************************************************************!
module triangle_class

  use shape_class

  implicit none
  
  private

  type, public, extends(shape) :: triangle
  
    real :: s = 0.
  
  contains
  
    procedure, public :: initialize => initialize_sub
    procedure, public :: area => get_area_fn
    procedure, public :: perimeter => get_perimeter_fn
    procedure, public :: to_string => to_string_fn
  
  end type triangle
  
contains

subroutine initialize_sub(this,s)

  implicit none
  
  class(triangle) :: this
  real :: s
  
  this%s = s

end subroutine initialize_sub

real function get_area_fn(this)

  implicit none

  class(triangle) :: this
  
  get_area_fn = sqrt(3.) / 4. * this%s**2.

end function get_area_fn

real function get_perimeter_fn(this)

  implicit none

  class(triangle) :: this
  
  get_perimeter_fn = 3. * this%s

end function get_perimeter_fn

character(len=50) function to_string_fn(this)
  
  implicit none
  
  class(triangle) :: this
  
  write(to_string_fn,'(a,f6.2)') 'Equilateral triangle of side ', this%s
  
end function to_string_fn

end module triangle_class

!**********************************************************************!
module rectangle_class

  use shape_class

  implicit none
  
  private

  type, public, extends(shape) :: rectangle
  
    real :: l = 0.
    real :: w = 0.
  
  contains
  
    procedure, public :: initialize => initialize_sub
    procedure, public :: area => get_area_fn
    procedure, public :: perimeter => get_perimeter_fn
    procedure, public :: to_string => to_string_fn
  
  end type rectangle
  
contains

subroutine initialize_sub(this,l,w)

  implicit none
  
  class(rectangle) :: this
  real :: l
  real :: w
  
  this%l = l
  this%w = w

end subroutine initialize_sub

real function get_area_fn(this)

  implicit none

  class(rectangle) :: this
  
  get_area_fn = this%l * this%w

end function get_area_fn

real function get_perimeter_fn(this)

  implicit none

  class(rectangle) :: this
  
  get_perimeter_fn = 2. * this%l + 2. * this%w

end function get_perimeter_fn

character(len=50) function to_string_fn(this)
  
  implicit none
  
  class(rectangle) :: this
  
  write(to_string_fn,'(a,f6.2,a,f6.2)') 'Rectangle of length ',  &
    this%l, ' and width ', this%w
  
end function to_string_fn

end module rectangle_class

!**********************************************************************!
module square_class

  use rectangle_class

  implicit none
  
  private

  type, public, extends(rectangle) :: square
  
  ! no new variables
  
  contains
  
    procedure, public :: to_string => to_string_fn
  
  end type square
  
contains

character(len=50) function to_string_fn(this)
  
  implicit none
  
  class(square) :: this
  
  write(to_string_fn,'(a,f6.2)') 'Square of length ', this%l
  
end function to_string_fn

end module square_class

!**********************************************************************!
module pentagon_class

  use shape_class

  implicit none
  
  private

  type, public, extends(shape) :: pentagon
  
    real :: s = 0.
  
  contains
  
    procedure, public :: initialize => initialize_sub
    procedure, public :: area => get_area_fn
    procedure, public :: perimeter => get_perimeter_fn
    procedure, public :: to_string => to_string_fn
  
  end type pentagon
  
contains

subroutine initialize_sub(this,s)

  implicit none
  
  class(pentagon) :: this
  real :: s

  this%s = s

end subroutine initialize_sub

real function get_area_fn(this)

  implicit none

  class(pentagon) :: this
  
  get_area_fn = 1.25 * this%s**2. / 0.72654253

end function get_area_fn

real function get_perimeter_fn(this)

  implicit none

  class(pentagon) :: this
  
  get_perimeter_fn = 5. * this%s

end function get_perimeter_fn

character(len=50) function to_string_fn(this)
  
  implicit none
  
  class(pentagon) :: this
  
  write(to_string_fn,'(a,f6.2,a,f6.2)') 'Pentagon of side ', this%s
  
end function to_string_fn

end module pentagon_class

!**********************************************************************!
!**********************************************************************!

program test_shape

  use circle_class
  use square_class
  use rectangle_class
  use pentagon_class
  use triangle_class
  use shape_class
  
  implicit none
  
  type(circle), pointer :: cir
  type(square), pointer :: squ
  type(rectangle), pointer :: rec
  type(pentagon), pointer :: pen
  type(triangle), pointer :: tri
  
  type :: shape_ptr
    class(shape), pointer :: p
  end type shape_ptr
  
  integer :: i
  character(len=50) :: id_string
  type(shape_ptr), dimension(5) :: shapes
  
  write(*,'(a)') 'Beginning of Fortran 2003 test...'
  
  allocate(cir)
  call cir%initialize(2.)
  
  allocate(squ)
  call squ%initialize(2.,2.)
  
  allocate(rec)
  call rec%initialize(2.,1.)
  
  allocate(tri)
  call tri%initialize(2.)
  
  allocate(pen)
  call pen%initialize(2.)
  
  shapes(1)%p => cir
  shapes(2)%p => squ
  shapes(3)%p => rec
  shapes(4)%p => tri
  shapes(5)%p => pen

  do i = 1, 5
    id_string = shapes(i)%p%to_string()
    write(*,'(/a)') id_string
    write(*,'(a,f8.4)') 'Area      = ', shapes(i)%p%area()
    write(*,'(a,f8.4)') 'Perimeter = ', shapes(i)%p%perimeter()
  enddo
  
  write(*,'(/a)') 'Successful completion of Fortran 2003 test!'

end program test_shape
