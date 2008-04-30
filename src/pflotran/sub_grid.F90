module Sub_Grid_module

  use Auxilliary_module
  
  implicit none

#include "definitions.h"

  private

  type, public :: sub_grid_type
    PetscInt :: id
    PetscInt :: num_cells
    PetscReal :: total_volume
    PetscReal :: total_distance
    PetscReal, pointer :: dist_up(:)
    PetscReal, pointer :: dist_dn(:)
    PetscReal, pointer :: area(:)
    PetscReal, pointer :: vol(:)
    type(auxilliary_type), pointer :: aux
    type(sub_grid_type), pointer :: next
  end type sub_grid_type


  ! pointer data structure required for making an array of region pointers in F90
  type, public :: sub_grid_ptr_type
    type(sub_grid_type), pointer :: ptr           ! pointer to the sub_grid_type
  end type sub_grid_ptr_type 
  
  type, public :: sub_grid_list_type
    PetscInt :: num_sub_grid_objects
    type(sub_grid_type), pointer :: first
    type(sub_grid_type), pointer :: last
    type(sub_grid_ptr_type), pointer :: array(:)
  end type sub_grid_list_type
  
  public :: SubGridCreate, SubGridAddToList, &
            SubGridAllocateLists, &
            SubGridGetNumberInList, &
            SubGridInitList, SubGridDestroyList, SubGridDestroy
  
contains

! ************************************************************************** !
!
! SubGridCreate: Allocates and initializes a new sub_grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function SubGridCreate(ncells,ndof,volume)

  implicit none
  
  PetscInt :: ncells
  PetscInt :: ndof
  PetscReal :: volume
  
  type(sub_grid_type), pointer :: SubGridCreate

  type(sub_grid_type), pointer :: sub_grid

  allocate(sub_grid)
  sub_grid%id = 0
  sub_grid%num_cells = 0
  nullify(sub_grid%dist_up)
  nullify(sub_grid%dist_dn)
  nullify(sub_grid%area)
  nullify(sub_grid%vol)
  nullify(sub_grid%aux)
  
  nullify(sub_grid%next)
  
  SubGridCreate => sub_grid

end function SubGridCreate

! ************************************************************************** !
!
! SubGridRead: Reads parameters associated with linear solver
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
subroutine SubGridRead(sub_grid,fid,myrank)

  use Fileio_module
  
  implicit none

  type(sub_grid_type) :: sub_grid
  PetscInt :: fid
  PetscMPIInt :: myrank
  
  character(len=MAXSTRINGLENGTH) :: string, error_string
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  PetscErrorCode :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadWord(string,keyword,.true.,ierr)
    call fiErrorMsg(myrank,'keyword','LINEAR SOLVER', ierr)
    call fiWordToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('NUMCELLS')
        call fiReadInteger(string,sub_grid%num_cells,ierr)
        call fiDefaultMsg(myrank,'num_cells',ierr)
      case('VOLUME')
        call fiReadDouble(string,sub_grid%volume,ierr)
        call fiDefaultMsg(myrank,'volume',ierr)

        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(myrank,'ksp_type','SUBGRID', ierr)   
        call fiWordToUpper(word)
        select case(trim(word))
          case('NONE','PREONLY')
            solver%ksp_type = KSPPREONLY
          case('GMRES')
            solver%ksp_type = KSPGMRES
          case('BCGS','BICGSTAB','BI-CGSTAB')
            solver%ksp_type = KSPBCGS
          case default
            string  = 'ERROR: Krylov solver type: ' // trim(word) // ' unknown.'
            if (myrank == 0) print *, string
            stop
        end select
      case('PRECONDITIONER_TYPE','PRECONDITIONER','PC','PC_TYPE')
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(myrank,'pc_type','SOLVER', ierr)   
        call fiWordToUpper(word)
        select case(trim(word))
          case('ILU','PCILU')
            solver%pc_type = PCILU
          case('LU','PCLU')
            solver%pc_type = PCLU
          case('BJACOBI','BLOCK_JACOBI')
            solver%pc_type = PCBJACOBI
          case('ASM','ADDITIVE_SCHWARTZ')
            solver%pc_type = PCASM
          case default
            string  = 'ERROR: Preconditioner type: ' // trim(word) // ' unknown.'
            if (myrank == 0) print *, string
            stop
        end select

      case('ATOL')
        call fiReadDouble(string,solver%linear_atol,ierr)
        call fiDefaultMsg(myrank,'linear_atol',ierr)

      case('RTOL')
        call fiReadDouble(string,solver%linear_rtol,ierr)
        call fiDefaultMsg(myrank,'linear_rtol',ierr)

      case('DTOL')
        call fiReadDouble(string,solver%linear_dtol,ierr)
        call fiDefaultMsg(myrank,'linear_dtol',ierr)
   
      case('MAXIT')
        call fiReadInt(string,solver%linear_maxit,ierr)
        call fiDefaultMsg(myrank,'linear_maxit',ierr)

      case default
        if (myrank == 0) print *, 'Keyword: '//keyword// &
                                  &' not recognized in linear solver'    

    end select 
  
  enddo  

end subroutine SubGridRead

! ************************************************************************** !
!
! SubGridGetNumberInList: Returns the number of sub_grids in a list
! author: Glenn Hammond
! date: 11/19/07
!
! ************************************************************************** !
function SubGridGetNumberInList(list)

  implicit none
  
  type(sub_grid_list_type) :: list

  PetscInt :: SubGridGetNumberInList
  type(sub_grid_type), pointer :: cur_sub_grid
  
  SubGridGetNumberInList = 0
  cur_sub_grid => list%first
  do
    if (.not.associated(cur_sub_grid)) exit
    SubGridGetNumberInList = SubGridGetNumberInList + 1
    cur_sub_grid => cur_sub_grid%next
  enddo

end function SubGridGetNumberInList

! ************************************************************************** !
!
! SubGridAllocateLists: Allocates sub_grids lists
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine SubGridAllocateLists(list)

  implicit none
  
  type(sub_grid_list_type), pointer :: list

  allocate(list)
  call SubGridInitList(list)
  
end subroutine

! ************************************************************************** !
!
! InitSubGridModule: Initializes module variables, lists, arrays.
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine SubGridInitList(list)

  implicit none

  type(sub_grid_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_sub_grid_objects = 0

end subroutine SubGridInitList

! ************************************************************************** !
!
! SubGridAddToList: Adds a new sub_grid of the module global list of 
!                      sub_grids
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine SubGridAddToList(new_sub_grid,list)

  implicit none
  
  type(sub_grid_type), pointer :: new_sub_grid
  type(sub_grid_list_type) :: list
  
  list%num_sub_grid_objects = list%num_sub_grid_objects + 1
  new_sub_grid%id = list%num_sub_grid_objects
  if (.not.associated(list%first)) list%first => new_sub_grid
  if (associated(list%last)) list%last%next => new_sub_grid
  list%last => new_sub_grid
  
end subroutine SubGridAddToList

! ************************************************************************** !
!
! SubGridConvertListToArray: Creates an array of pointers to the 
!                               sub_grids in the sub_grid list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine SubGridConvertListToArray(list)

  implicit none
  
  type(sub_grid_list_type) :: list
    
  PetscInt :: count
  type(sub_grid_type), pointer :: cur_sub_grid
  
  
  allocate(list%array(list%num_sub_grid_objects))
  
  cur_sub_grid => list%first
  do 
    if (.not.associated(cur_sub_grid)) exit
    list%array(cur_sub_grid%id)%ptr => cur_sub_grid
    cur_sub_grid => cur_sub_grid%next
  enddo

end subroutine SubGridConvertListToArray

! ************************************************************************** !
!
! SubGridDestroy: Deallocates a sub_grid
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine SubGridDestroy(sub_grid)

  implicit none
  
  type(sub_grid_type), pointer :: sub_grid
  
  if (.not.associated(sub_grid)) return
  
  if (associated(sub_grid%dist_up)) deallocate(sub_grid%dist_up)
  nullify(sub_grid%dist_up)
  if (associated(sub_grid%dist_dn)) deallocate(sub_grid%dist_dn)
  nullify(sub_grid%dist_dn)
  if (associated(sub_grid%area)) deallocate(sub_grid%area)
  nullify(sub_grid%area)
  if (associated(sub_grid%vol)) deallocate(sub_grid%vol)
  nullify(sub_grid%vol)
  
  call AuxDestroy(sub_grid%aux)
  
  nullify(sub_grid%next)
  
  deallocate(sub_grid)
  nullify(sub_grid)

end subroutine SubGridDestroy

! ************************************************************************** !
!
! SubGridDestroyList: Deallocates the module global list and array of regions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine SubGridDestroyList(list)

  implicit none
  
  type(sub_grid_list_type), pointer :: list
    
  type(sub_grid_type), pointer :: cur_sub_grid, prev_sub_grid
  
  if (.not.associated(list)) return
  
  if (associated(list%array)) deallocate(list%array)
  nullify(list%array)
  
  cur_sub_grid => list%first
  do 
    if (.not.associated(cur_sub_grid)) exit
    prev_sub_grid => cur_sub_grid
    cur_sub_grid => cur_sub_grid%next
    call SubGridDestroy(prev_sub_grid)
  enddo
  
  nullify(list%first)
  nullify(list%last)
  list%num_sub_grid_objects = 0
  
  deallocate(list)
  nullify(list)

end subroutine SubGridDestroyList

end module Sub_Grid_module
