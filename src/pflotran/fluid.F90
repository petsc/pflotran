module Fluid_module
 
  implicit none

  private

#include "definitions.h"
 
  type, public :: fluid_property_type
    PetscReal :: tort_bin_diff
    PetscReal :: vap_air_diff_coef
    PetscReal :: exp_binary_diff
    PetscReal :: enh_binary_diff_coef
    PetscReal :: diff_base
    PetscReal :: diff_exp
    character(len=MAXWORDLENGTH) :: phase_name
    PetscInt :: phase_id
    PetscReal :: diff_coef
    type(fluid_property_type), pointer :: next
  end type fluid_property_type
  
  public :: FluidPropertyCreate, FluidPropertyDestroy, &
            FluidPropertyRead, FluidPropertyAddToList

contains

! ************************************************************************** !
!
! FluidPropertyCreate: Creates a fluid property object
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
function FluidPropertyCreate()
  
  implicit none

  type(fluid_property_type), pointer :: FluidPropertyCreate
  
  type(fluid_property_type), pointer :: fluid_property
  
  allocate(fluid_property)
  fluid_property%tort_bin_diff = 0.d0
  fluid_property%vap_air_diff_coef = 0.d0
  fluid_property%exp_binary_diff = 0.d0
  fluid_property%enh_binary_diff_coef = 0.d0
  fluid_property%diff_base = 0.d0
  fluid_property%diff_exp = 0.d0
  fluid_property%phase_name = ''
  fluid_property%phase_id = 0
  fluid_property%diff_coef = 0.d0
  nullify(fluid_property%next)
  FluidPropertyCreate => fluid_property

end function FluidPropertyCreate

! ************************************************************************** !
!
! FluidPropertyRead: Reads in contents of a fluid property card
! author: Glenn Hammond
! date: 01/21/09
! 
! ************************************************************************** !
subroutine FluidPropertyRead(fluid_property,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  type(fluid_property_type) :: fluid_property
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','FLUID_PROPERTY')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('PHASE') 
        call InputReadWord(input,option,fluid_property%phase_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'phase','FLUID_PROPERTY')
      case('DIFFUSION_COEFFICIENT') 
        call InputReadDouble(input,option,fluid_property%diff_coef)
        call InputErrorMsg(input,option,'diffusion coefficient','FLUID_PROPERTY')

      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in fluid property'    
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine FluidPropertyRead

! ************************************************************************** !
!
! FluidPropertyAddToList: Adds a thermal property to linked list
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
subroutine FluidPropertyAddToList(fluid_property,list)

  implicit none
  
  type(fluid_property_type), pointer :: fluid_property
  type(fluid_property_type), pointer :: list

  type(fluid_property_type), pointer :: cur_fluid_property
  
  if (associated(list)) then
    cur_fluid_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_fluid_property%next)) exit
      cur_fluid_property => cur_fluid_property%next
    enddo
    cur_fluid_property%next => fluid_property
  else
    list => fluid_property
  endif
  
end subroutine FluidPropertyAddToList

! ************************************************************************** !
!
! FluidPropertyDestroy: Destroys a fluid property
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
recursive subroutine FluidPropertyDestroy(fluid_property)

  implicit none
  
  type(fluid_property_type), pointer :: fluid_property
  
  if (.not.associated(fluid_property)) return
  
  call FluidPropertyDestroy(fluid_property%next)

  deallocate(fluid_property)
  nullify(fluid_property)
  
end subroutine FluidPropertyDestroy

end module Fluid_module
