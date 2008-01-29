module Connection_module

  implicit none

#include "definitions.h"

  private

  type, public :: connection_type
    PetscInt :: id
    PetscInt :: itype                  ! connection type (boundary, internal, source sink
    PetscInt :: num_connections
    PetscInt, pointer :: id_up(:)      ! list of ids of upwind cells
    PetscInt, pointer :: id_dn(:)      ! list of ids of downwind cells
    PetscReal, pointer :: dist(:,:)      ! list of distance vectors, size(-1:3,num_connections) where
                                      !   -1 = fraction upwind
                                      !   0 = magnitude of distance 
                                      !   1-3 = components of unit vector
    PetscReal, pointer :: area(:)        ! list of areas of faces normal to distance vectors
!    PetscReal, pointer :: velocity(:,:)  ! velocity scalars for each phase
    type(connection_type), pointer :: next
  end type connection_type


  ! pointer data structure required for making an array of region pointers in F90
  type, public :: connection_ptr_type
    type(connection_type), pointer :: ptr           ! pointer to the connection_type
  end type connection_ptr_type 
  
  type, public :: connection_list_type
    PetscInt :: num_connection_objects
    type(connection_type), pointer :: first
    type(connection_type), pointer :: last
    type(connection_ptr_type), pointer :: array(:)
  end type connection_list_type
  
  type(connection_list_type), pointer, private :: internal_connection_list, &
                                                  boundary_connection_list


  public :: ConnectionCreate, ConnectionAddToList, &
            ConnectionAllocateLists, &
            ConnectionGetNumberInList, &
            ConnectionInitList, ConnectionDestroyList, ConnectionDestroy
  
contains

! ************************************************************************** !
!
! ConnectionCreate: Allocates and initializes a new connection
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function ConnectionCreate(num_connections,num_dof,connection_itype)

  implicit none
  
  PetscInt :: num_connections
  PetscInt :: num_dof
  PetscInt :: connection_itype
  
  type(connection_type), pointer :: ConnectionCreate

  type(connection_type), pointer :: connection

  allocate(connection)
  connection%id = 0
  connection%itype = connection_itype
  connection%num_connections = num_connections
  nullify(connection%id_up)
  nullify(connection%id_dn)
  nullify(connection%dist)
  nullify(connection%area)
!  nullify(connection%velocity)
  select case(connection_itype)
    case(INTERNAL_CONNECTION_TYPE)
      allocate(connection%id_up(num_connections))
      allocate(connection%id_dn(num_connections))
      allocate(connection%dist(-1:3,num_connections))
      allocate(connection%area(num_connections))
!      allocate(connection%velocity(num_dof,num_connections))
      connection%id_up = 0
      connection%id_dn = 0
      connection%dist = 0.d0
      connection%area = 0.d0
!      connection%velocity = 0.d0
    case(BOUNDARY_CONNECTION_TYPE)
      allocate(connection%id_dn(num_connections))
      allocate(connection%dist(-1:3,num_connections))
      allocate(connection%area(num_connections))
!      allocate(connection%velocity(num_dof,num_connections))
      connection%id_dn = 0
      connection%dist = 0.d0
      connection%area = 0.d0
!      connection%velocity = 0.d0
    case(SRC_SINK_CONNECTION_TYPE,INITIAL_CONNECTION_TYPE)
      allocate(connection%id_dn(num_connections))
      connection%id_dn = 0
  end select
  nullify(connection%next)
  
  ConnectionCreate => connection

end function ConnectionCreate

! ************************************************************************** !
!
! ConnectionGetNumberInList: Returns the number of connections in a list
! author: Glenn Hammond
! date: 11/19/07
!
! ************************************************************************** !
function ConnectionGetNumberInList(list)

  implicit none
  
  type(connection_list_type) :: list

  PetscInt :: ConnectionGetNumberInList
  type(connection_type), pointer :: cur_connection_set
  
  ConnectionGetNumberInList = 0
  cur_connection_set => list%first
  do
    if (.not.associated(cur_connection_set)) exit
    ConnectionGetNumberInList = ConnectionGetNumberInList + &
                                cur_connection_set%num_connections
    cur_connection_set => cur_connection_set%next
  enddo

end function ConnectionGetNumberInList

! ************************************************************************** !
!
! ConnectionAllocateLists: Allocates connections lists
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine ConnectionAllocateLists()

  implicit none
  
  allocate(internal_connection_list)
  call ConnectionInitList(internal_connection_list)
  allocate(boundary_connection_list)
  call ConnectionInitList(boundary_connection_list)
  
end subroutine

! ************************************************************************** !
!
! InitConnectionModule: Initializes module variables, lists, arrays.
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine ConnectionInitList(list)

  implicit none

  type(connection_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_connection_objects = 0

end subroutine ConnectionInitList

! ************************************************************************** !
!
! ConnectionAddToList: Adds a new connection of the module global list of 
!                      connections
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine ConnectionAddToList(new_connection,list)

  implicit none
  
  type(connection_type), pointer :: new_connection
  type(connection_list_type) :: list
  
  list%num_connection_objects = list%num_connection_objects + 1
  new_connection%id = list%num_connection_objects
  if (.not.associated(list%first)) list%first => new_connection
  if (associated(list%last)) list%last%next => new_connection
  list%last => new_connection
  
end subroutine ConnectionAddToList

! ************************************************************************** !
!
! ConnectionConvertListToArray: Creates an array of pointers to the 
!                               connections in the connection list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine ConnectionConvertListToArray(list)

  implicit none
  
  type(connection_list_type) :: list
    
  PetscInt :: count
  type(connection_type), pointer :: cur_connection
  
  
  allocate(list%array(list%num_connection_objects))
  
  cur_connection => list%first
  do 
    if (.not.associated(cur_connection)) exit
    list%array(cur_connection%id)%ptr => cur_connection
    cur_connection => cur_connection%next
  enddo

end subroutine ConnectionConvertListToArray

! ************************************************************************** !
!
! ConnectionDestroy: Deallocates a connection
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine ConnectionDestroy(connection)

  implicit none
  
  type(connection_type), pointer :: connection
  
  if (.not.associated(connection)) return
  
  if (associated(connection%id_up)) deallocate(connection%id_up)
  nullify(connection%id_up)
  if (associated(connection%id_dn)) deallocate(connection%id_dn)
  nullify(connection%id_dn)
  if (associated(connection%dist)) deallocate(connection%dist)
  nullify(connection%dist)
  if (associated(connection%area)) deallocate(connection%area)
  nullify(connection%area)
!  if (associated(connection%velocity)) deallocate(connection%velocity)
!  nullify(connection%velocity)
  nullify(connection%next)
  
  deallocate(connection)
  nullify(connection)

end subroutine ConnectionDestroy

! ************************************************************************** !
!
! ConnectionDestroyList: Deallocates the module global list and array of regions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine ConnectionDestroyList(list)

  implicit none
  
  type(connection_list_type), pointer :: list
    
  type(connection_type), pointer :: cur_connection, prev_connection
  
  if (.not.associated(list)) return
  
  if (associated(list%array)) deallocate(list%array)
  nullify(list%array)
  
  cur_connection => list%first
  do 
    if (.not.associated(cur_connection)) exit
    prev_connection => cur_connection
    cur_connection => cur_connection%next
    call ConnectionDestroy(prev_connection)
  enddo
  
  nullify(list%first)
  nullify(list%last)
  list%num_connection_objects = 0
  
  deallocate(list)
  nullify(list)

end subroutine ConnectionDestroyList

end module Connection_module
