module Integral_Flux_module

  use Geometry_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public :: integral_flux_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    type(point3d_type), pointer :: coordinates(:)
    PetscInt, pointer :: connections(:)
    type(integral_flux_type), pointer :: next
  end type integral_flux_type
  
  type, public :: integral_flux_list_type
    PetscInt :: num_integral_fluxes
    type(integral_flux_type), pointer :: first
    type(integral_flux_type), pointer :: last
    type(integral_flux_type), pointer :: array(:)
  end type integral_flux_list_type

  public :: IntegralFluxCreate, &
            IntegralFluxDestroy, &
            IntegralFluxRead, &
            IntegralFluxAddToList, &
            IntegralFluxInitList, &
            IntegralFluxDestroyList, &
            IntegralFluxGetPtrFromList, &
            IntegralFluxRemoveFromList

  interface IntegralFluxCreate
    module procedure IntegralFluxCreate1
    module procedure IntegralFluxCreateFromIntegralFlux
  end interface
    
contains

! ************************************************************************** !

function IntegralFluxCreate1()
  ! 
  ! Create object that stores integral flux information
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: IntegralFluxCreate1
  
  type(integral_flux_type), pointer :: integral_flux
  
  allocate(integral_flux)
  
  integral_flux%name = ""
  integral_flux%id = 0
  nullify(integral_flux%connections)
  nullify(integral_flux%coordinates)
  nullify(integral_flux%next)
  
  IntegralFluxCreate1 => integral_flux

end function IntegralFluxCreate1

! ************************************************************************** !

function IntegralFluxCreateFromIntegralFlux(integral_flux)
  ! 
  ! Create object that stores integral flux information
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: IntegralFluxCreateFromIntegralFlux
  type(integral_flux_type), pointer :: integral_flux

  type(integral_flux_type), pointer :: new_integral_flux
  
  new_integral_flux => IntegralFluxCreate1()
  
  new_integral_flux%name = integral_flux%name
  new_integral_flux%id = integral_flux%id
  
  nullify(new_integral_flux%next)
  
  IntegralFluxCreateFromIntegralFlux => new_integral_flux

end function IntegralFluxCreateFromIntegralFlux

! ************************************************************************** !

subroutine IntegralFluxRead(integral_flux,input,option)
  ! 
  ! Reads integral flux data from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(integral_flux_type) :: integral_flux
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','integral_flux')   
      
    select case(trim(keyword))
      case('COORDINATES')
        call GeometryReadCoordinates(input,option,integral_flux%name, &
                                     integral_flux%coordinates)
      case default
        option%io_buffer = 'Keyword (' // trim(keyword) // &
                           ') not recognized under integral_flux.'
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine IntegralFluxRead

! ************************************************************************** !

subroutine IntegralFluxInitList(list)
  ! 
  ! Initializes a integral flux list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none

  type(integral_flux_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_integral_fluxes = 0

end subroutine IntegralFluxInitList

! ************************************************************************** !

subroutine IntegralFluxAddToList(new_integral_flux,list)
  ! 
  ! Adds a new integral_flux to a integral flux list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: new_integral_flux
  type(integral_flux_list_type) :: list
  
  list%num_integral_fluxes = list%num_integral_fluxes + 1
  new_integral_flux%id = list%num_integral_fluxes
  if (.not.associated(list%first)) list%first => new_integral_flux
  if (associated(list%last)) list%last%next => new_integral_flux
  list%last => new_integral_flux
  
end subroutine IntegralFluxAddToList

! ************************************************************************** !

subroutine IntegralFluxRemoveFromList(integral_flux,list)
  ! 
  ! Removes a integral_flux from a integral flux list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: integral_flux
  type(integral_flux_list_type) :: list
  
  type(integral_flux_type), pointer :: cur_integral_flux, prev_integral_flux
  
  cur_integral_flux => list%first
  nullify(prev_integral_flux)
  
  do
    if (.not.associated(cur_integral_flux)) exit
    if (associated(cur_integral_flux,integral_flux)) then
      if (associated(prev_integral_flux)) then
        prev_integral_flux%next => cur_integral_flux%next
      else
        list%first => cur_integral_flux%next
      endif
      if (.not.associated(cur_integral_flux%next)) then
        list%last => prev_integral_flux
      endif
      list%num_integral_fluxes = list%num_integral_fluxes-1
      call IntegralFluxDestroy(cur_integral_flux)
      return
    endif
    prev_integral_flux => cur_integral_flux
    cur_integral_flux => cur_integral_flux%next
  enddo
  
end subroutine IntegralFluxRemoveFromList

! ************************************************************************** !

function IntegralFluxGetPtrFromList(integral_flux_name,integral_flux_list)
  ! 
  ! Returns a pointer to the integral flux matching integral_flux_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use String_module

  implicit none
  
  type(integral_flux_type), pointer :: IntegralFluxGetPtrFromList
  character(len=MAXWORDLENGTH) :: integral_flux_name
  type(integral_flux_list_type) :: integral_flux_list
 
  PetscInt :: length
  type(integral_flux_type), pointer :: integral_flux
    
  nullify(IntegralFluxGetPtrFromList)
  integral_flux => integral_flux_list%first
  
  do 
    if (.not.associated(integral_flux)) exit
    length = len_trim(integral_flux_name)
    if (length == len_trim(integral_flux%name) .and. &
        StringCompare(integral_flux%name,integral_flux_name, &
                        length)) then
      IntegralFluxGetPtrFromList => integral_flux
      return
    endif
    integral_flux => integral_flux%next
  enddo
  
end function IntegralFluxGetPtrFromList

! ************************************************************************** !

subroutine IntegralFluxDestroyList(integral_flux_list)
  ! 
  ! Deallocates a list of integral fluxes
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_list_type), pointer :: integral_flux_list
  
  type(integral_flux_type), pointer :: integral_flux, prev_integral_flux
  
  if (.not.associated(integral_flux_list)) return
  
  integral_flux => integral_flux_list%first
  do 
    if (.not.associated(integral_flux)) exit
    prev_integral_flux => integral_flux
    integral_flux => integral_flux%next
    call IntegralFluxDestroy(prev_integral_flux)
  enddo
  
  integral_flux_list%num_integral_fluxes = 0
  nullify(integral_flux_list%first)
  nullify(integral_flux_list%last)
  if (associated(integral_flux_list%array)) deallocate(integral_flux_list%array)
  nullify(integral_flux_list%array)
  
  deallocate(integral_flux_list)
  nullify(integral_flux_list)

end subroutine IntegralFluxDestroyList

! ************************************************************************** !

subroutine IntegralFluxDestroy(integral_flux)
  ! 
  ! Deallocates a integral flux
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: integral_flux
  
  PetscInt :: i
  
  if (.not.associated(integral_flux)) return
  
  if (associated(integral_flux%coordinates)) &
    deallocate(integral_flux%coordinates)
  nullify(integral_flux%coordinates)
  if (associated(integral_flux%connections)) &
    deallocate(integral_flux%connections)
  nullify(integral_flux%connections)
  deallocate(integral_flux)
  nullify(integral_flux)

end subroutine IntegralFluxDestroy

end module Integral_Flux_module
