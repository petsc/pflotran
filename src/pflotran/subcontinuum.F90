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
  
  type, public :: subcontinuum_field_types
    ! NOTE (Jitu): Check the required fields again and add/delete as needed
    PetscInt :: num_subgrids
    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: tran_xx, tran_dxx, tran_yy, tran_accum
    ! Vectors for operator splitting
    Vec :: tran_rhs, tran_rhs_coef

    Vec :: tran_log_xx
    Vec :: tran_ts_mass_balance, tran_total_mass_balance

  end type subcontinuum_field_typec

  type, public :: subcontinuum_field_typec
    PetscInt :: num_continuum
    type(subcontinuum_field_types), pointer :: subcontinuum_field_continuum
  end type subcontinuum_field_typec

  type, public :: subcontinuum_field_typep
    PetscInt :: num_cells
    type(subcontinuum_field_typec), pointer :: subcontinuum_field_cell
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


! ************************************************************************ !
!
! SubcontinuumFieldCreatePatch: Create subcontinuum field at a patch 
! author: Jitendra Kumar 
! date: 11/29/2010
!
! ************************************************************************ !
subroutine SubcontinuumFieldCreatePatch(patch)
  
  implicit none

  type(patch_type), pointer :: patch

  type(subcontinuum_field_typep), pointer :: subcontinuum_field_patch

  PetscInt :: icell, isub, num_cells, num_continuum, num_subgrids, offset

  subcontinuum_field_patch => patch%subcontinuum_field_patch 
  
  ! Save no. of local grid cells in subcontinuum_field_patch object  
  num_cells = patch%nlmax
  subcontinuum_field_patch%num_cells = num_cells
  ! Allocate storage for subcontinuum_field_node pointers for every cell
  allocate(subcontinuum_field_patch%subcontinuum_field_cell(num_cells))
  
  ! Loop over each cell and allocate subcontinuum_field_cell
  
  do icell=1,num_cells
    num_continuum = patch%num_subcontinuum_type(icell,1)
    subcontinuum_field_patch%subcontinuum_field_cell%num_continuum = &
                                                         num_continuum
    allocate(subcontinuum_field_patch%subcontinuum_field_cell(num_continuum)
    
    ! Loop over each subcontinuum and allocate the Petsc vectors for the
    ! based on the subgrid mesh information
    offset = patch%subcontinuum_grid_offset(icell,2)
    do isub=1,num_continuum
      ! Get the num_subgrid associated with the current subcontinuum
      num_subgrids = patch%subcontinuum_grid(offset+isub-1)
      ! Allocate Petsc vectors for this subgrid
      VecCreateSeq(PETSC_COMM_SELF, num_subgrids,  &
          subcontinuum_field_patch%subcontinuum_field_cell(icell)%subcontinuum_field_continuum(isub)%tran_xx)
      
    enddo
  enddo

end subroutine SubcontinuumFieldCreatePatch


! ************************************************************************ !
!
! SubcontinuumFieldCreateRealization: Create subcontinuum field for a
! realization 
! author: Jitendra Kumar 
! date: 11/29/2010
!
! ************************************************************************ !
subroutine SubcontinuumFieldCreateRealization(realization)
  
  use Realization_module
  use Patch_module
  use Subcontinuum_module
  use Level_module

  implicit none

  type(realization_type), pointer :: realization
  type(patch_type), pointer :: patch
  type(level_type), pointer :: level
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
    
  ! Jump to last/finest level
  cur_level => realization%level_list%last
  if(.not.associated(cur_level)) exit

  ! Loop over all patches at this level and create the subcontinuum field
  cur_patch => cur_level%patch_list%first
  do 
    call SubcontinuumFieldCreatePatch(cur_patch)
    cur_patch => cur_patch%next
  enddo

end subroutine SubcontinuumFieldCreateRealization

! ********************************************************************** !
!
! GridPopulateSubcontinuum: Populate subcontinuum discretization info in
!                           the grid object
! author: Jitendra Kumar
! date: 11/22/2010
!
! *********************************************************************** !
subroutine GridPopulateSubcontinuum(realization)

  use Realization_module
  use Discretization_module
  use Material_module
  use Grid_module
  use Patch_module
  use Level_module

  implicit none

  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: patch 
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(subcontinuum_type), pointer :: subcontinuum_type

  PetscInt :: icell, ssub

  ! do this only at the last and finest level
  cur_level => realization%level_list%last
  if(.not.associated(cur_level)) exit
  cur_patch => cur_level%patch_list%first
  do 
    grid => cur_patch%grid
    ! Allocate storage for subcontinuum_grid_offset
    allocate(cur_patch%grid%subcontinuum_grid_offset(cur_patch%grid%nlmax,2)
    
    ! Loop through all cells in the patch, copy the no. of subcontinuum 
    ! and offset from patch%num_subcontinuum_type into
    ! grid%subcontinuum_grid_offset.
    ssub = 0
    do icell=1, cur_patch%grid%nlmax
      cur_patch%grid%subcontinuum_grid_offset(icell,1) =   &
                          cur_patch%num_subcontinuum_type(icell,1)
      cur_patch%grid%subcontinuum_grid_offset(icell,2) =   &
                          cur_patch%num_subcontinuum_type(icell,2)
      ssub = ssub + cur_patch%grid%subcontinuum_grid_offset(icell,1)
    enddo

    ! Allocate storage for grid%subcontinuum_grid and store the subgrid
    ! information
    allocate(cur_patch%grid%subcontinuum_grid(ssub))
    
    ! Loop through all the cells-> all subcontinuum and set the subgrid
    ! information
    do icell=1, cur_patch%grid%nlmax
      if(associated(region)) the
        local_id = region%cell_ids(icell)
      else
        local_id = icell
      endif

      ! Loop over all subcontinua
      offset = cur_patch%grid%subcontinuum_grid_offset(icell,2)
      do isub = 1, cur_patch%grid%subcontinuum_grid_offset(icell,1)
        cur_patch%grid%subcontinuum_grid(offset) =  &
        realization%subcontinuum_properties(cur_patch%subcontinuum_type_ids(offset))%num_subgrids
        offset = offset + 1
      enddo
    enddo
    cur_patch => cur_patch%next
  enddo           
end subroutine GridPopulateSubcontinuum


end module Subcontinuum_module































