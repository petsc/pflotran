module Subcontinuum_module

  implicit none

  private

#include "definitions.h"

  type, public :: subcontinuum_property_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: num_subgrids
    PetscReal :: radius
    type(subcontinuum_property_type), pointer :: next
  end type subcontinuum_property_type

  type, public :: subcontinuum_property_ptr_type
    type(subcontinuum_property_type), pointer :: ptr
  end type subcontinuum_property_ptr_type
  
  type, public :: subcontinuum_field_typec
    ! NOTE (Jitu): Check the required fields again and add/delete as needed
    PetscInt :: num_subgrids
    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: tran_xx, tran_dxx, tran_yy, tran_accum
    ! Vectors for operator splitting
    Vec :: tran_rhs, tran_rhs_coef

    Vec :: tran_log_xx
    Vec :: tran_ts_mass_balance, tran_total_mass_balance

  end type subcontinuum_field_typec

  type, public :: subcontinuum_field_typen
    PetscInt :: num_continuum
    type(subcontinuum_field_typec), pointer :: sub_field_continuum
  end type subcontinuum_field_typen

  type, public :: subcontinuum_field_typep
    PetscInt :: num_cells
    type(subcontinuum_field_typen), pointer :: sub_field_node
  end type subcontinuum_field_typep  

  public :: SubcontinuumPropertyCreate, &
            SubcontinuumPropertyDestroy, &
            SubcontinuumPropertyAddToList, &
            SubcontinuumPropGetPtrFromList, &
            SubcontinuumPropGetPtrFromArray, &
            SubcontinuumPropConvertListToArray, &
            SubcontinuumPropertyRead

contains

! ************************************************************************ !
!
! SubcontinuumPropertyCreate: Creates a subcontinuum property
! author: Jitendra Kumar
! date: 10/04/2010
!
! ************************************************************************ !
function SubcontinuumPropertyCreate()

  implicit none

  type(subcontinuum_property_type), pointer :: SubcontinuumPropertyCreate

  type(subcontinuum_property_type), pointer :: subcontinuum_property
  
  allocate(subcontinuum_property)
  subcontinuum_property%id = 0
  subcontinuum_property%num_subgrids = 0
  subcontinuum_property%name = ''
  subcontinuum_property%radius = 0.d0
  nullify(subcontinuum_property%next)
  SubcontinuumPropertyCreate => subcontinuum_property

end function SubcontinuumPropertyCreate


! ************************************************************************ !
!
! SubcontinuumPropertyRead: Reads in contents of a subcontinuum property card
! author: Jitendra Kumar
! date: 10/04/2010
!
! ************************************************************************ !
subroutine SubcontinuumPropertyRead(subcontinuum_property,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none

  type(subcontinuum_property_type) :: subcontinuum_property
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXWORDLENGTH) :: string

  PetscInt :: length
  
  input%err = 0

  do

    call InputReadFlotranString(input,option)
    
    if(InputCheckExit(input,option)) exit

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SUBCONTINUUM_PROPERTY')
    call StringToUpper(keyword)

    select case(trim(keyword))

      case('NAME')
        call InputReadWord(input,option,subcontinuum_property%name,PETSC_TRUE)
        call InputErrormsg(input,option,'name','SUBCONTINUUM_PROPERTY')
      case('ID')
        call InputReadInt(input,option,subcontinuum_property%id)
        call InputErrorMsg(input,option,'id','SUBCONTINUUM_PROPERTY')
      case('RADIUS')
        call InputReadDouble(input,option,subcontinuum_property%radius)
        call InuputErrorMsg(input,option,'radius','SUBCONTINUUM_PROPERTY')
      case('SUBGRID')
        call InputReadInt(input,option,subcontinuum_property%radius)
        call InuputErrorMsg(input,option,'num_subgrids','SUBCONTINUUM_PROPERTY')
      case default
        option%io_buffer = 'Keyword ('// trim(keyword) // ') not recognized in
        subcontinuum_property'
        call printErrMsg(option)
    end select

  enddo   


end subroutine SubcontinuumPropertyRead   



! ************************************************************************ !
!
! SubcontinuumPropertyAddToList: Adds a subcontinuum property to linked list
! author: Jitendra Kumar
! date: 10/04/2010
!
! ************************************************************************ !
subroutine SubcontinuumPropertyAddToList(subcontinuum_property,list)

  implicit none

  type(subcontinuum_property_type), pointer :: subcontinuum_property
  type(subcontinuum_property_type), pointer :: list

  type(subcontinuum_property_type), pointer :: cur_subcontinuum_property

  if(associated(list)) then
    cur_material_property => list
    ! loop to the end of list
    do
      if (.not.associated(cur_subcontinuum_property%next)) exit
      cur_subcontinuum_property => cur_subcontinuum_property%next
    enddo
    cur_subcontinuum_property%next => subcontinuum_property
  else
    list => subcontinuum_property
  endif

end subroutine SubcontinuumPropertyAddToList

! ************************************************************************ !
!
! SubcontinuumPropConvertListToArray: Creates an array of pointers to the 
!                                subcontinuum_properties in the list
! author: Jitendra Kumar 
! date: 10/04/2010  
!
! ************************************************************************ !
subroutine SubcontinuumPropConvertListToArray(list,array)

  implicit none
  
  type(subcontinuum_property_type), pointer :: list
  type(subcontinuum_property_ptr_type), pointer :: array(:)
    
  type(subcontinuum_property_type), pointer :: cur_subcontinuum_property
  type(subcontinuum_property_type), pointer :: prev_subcontinuum_property
  type(subcontinuum_property_type), pointer :: next_subcontinuum_property
  PetscInt :: i, max_id

  max_id = 0
  cur_subcontinuum_property => list
  do 
    if (.not.associated(cur_subcontinuum_property)) exit
    max_id = max(max_id,cur_subcontinuum_property%id)
    cur_subcontinuum_property => cur_subcontinuum_property%next
  enddo
  
  allocate(array(max_id))
  do i = 1, max_id
    nullify(array(i)%ptr)
  enddo
  
  cur_subcontinuum_property => list
  do 
    if (.not.associated(cur_subcontinuum_property)) exit
    array(cur_subcontinuum_property%id)%ptr => cur_subcontinuum_property
    cur_subcontinuum_property => cur_subcontinuum_property%next
  enddo

end subroutine SubcontinuumPropConvertListToArray

! ************************************************************************ !
!
! SubcontinuumPropGetPtrFromList: Returns a pointer to the subcontinuum 
!                                 property matching subcontinuum_name
!                             
! author: Jitendra Kumar 
! date: 10/04/2010
!
! ************************************************************************ !
function SubcontinuumPropGetPtrFromList(subcontinuum_property_name, &
                                    subcontinuum_property_list)

  use String_module
  
  implicit none
  
  type(material_property_type), pointer :: SubcontinuumPropGetPtrFromList
  character(len=MAXWORDLENGTH) :: subcontinuum_property_name
  type(subcontinuum_property_type), pointer :: subcontinuum_property_list
  PetscInt :: length
  type(subcontinuum_property_type), pointer :: subcontinuum_property
    
  nullify(SubcontinuumPropGetPtrFromList)
  subcontinuum_property => subcontinuum_property_list
  
  do 
    if (.not.associated(subcontinuum_property)) exit
    length = len_trim(subcontinuum_property_name)
    if (length == len_trim(subcontinuum_property%name) .and. &
        StringCompare(subcontinuum_property%name,subcontinuum_property_name,length)) then
      SubcontinuumPropGetPtrFromList => subcontinuum_property
      return
    endif
    subcontinuum_property => subcontinuum_property%next
  enddo
  
end function SubcontinuumPropGetPtrFromList

! ************************************************************************ !
!
! SubcontinuumPropGetPtrFromArray: Returns a pointer to the 
!                     subcontinuum property matching subcontinuum_name
!
! author: Jitendra Kumar 
! date: 10/04/2010
!
! ************************************************************************ !
function SubcontinuumPropGetPtrFromArray(subcontinuum_property_name, &
                                     subcontinuum_property_array)

  use String_module

  implicit none
  
  type(subcontinuum_property_type), pointer :: SubcontinuumPropGetPtrFromArray
  character(len=MAXWORDLENGTH) :: subcontinuum_property_name
  type(subcontinuum_property_ptr_type), pointer :: subcontinuum_property_array(:)
  PetscInt :: length
  PetscInt :: isubcontinuum_property
    
  nullify(SubcontinuumPropGetPtrFromArray)
  
  do isubcontinuum_property = 1, size(subcontinuum_property_array)
    length = len_trim(subcontinuum_property_name)
    if (length == len_trim(subcontinuum_property_array(isubcontinuum_property)%ptr%name) .and. &
        StringCompare(subcontinuum_property_array(isubcontinuum_property)%ptr%name, &
                        subcontinuum_property_name,length)) then
      SubcontinuumPropGetPtrFromArray => subcontinuum_property_array(isubcontinuum_property)%ptr
      return
    endif
  enddo
  
end function SubcontinuumPropGetPtrFromArray

! ************************************************************************ !
!
! SubcontinuumPropertyDestroy: Destroys a subcontinuum_property
! author: Jitendra Kumar 
! date: 10/04/2010
!
! ************************************************************************ !
recursive subroutine SubcontinuumPropertyDestroy(subcontinuum_property)

  implicit none
  
  type(subcontinuum_property_type), pointer :: subcontinuum_property
  
  if (.not.associated(subcontinuum_property)) return
  
  call SubcontinuumPropertyDestroy(subcontinuum_property%next)
  
  deallocate(subcontinuum_property)
  nullify(subcontinuum_property)
  
end subroutine SubcontinuumPropertyDestroy

end module Subcontinuum_module































