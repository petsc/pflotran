module Connection_module

  implicit none

#include "definitions.h"

  private

  type, public :: connection_type
    integer :: id
    integer :: num_connections
    integer, pointer :: id_up(:)      ! list of ids of upwind cells
    integer, pointer :: id_dn(:)      ! list of ids of downwind cells
    real*8, pointer :: dist(:,:)      ! list of distance vectors, size(-1:3,num_connections) where
                                      !   -1 = fraction upwind
                                      !   0 = magnitude of distance 
                                      !   1-3 = components of unit vector
    real*8, pointer :: area(:)        ! list of areas of faces normal to distance vectors
    real*8, pointer :: velocity(:,:)  ! velocity scalars for each phase
    type(connection_type), pointer :: next
  end type connection_type


  ! pointer data structure required for making an array of region pointers in F90
  type, public :: connection_ptr_type
    type(connection_type), pointer :: ptr           ! pointer to the connection_type
  end type connection_ptr_type 
  
  type, public :: connection_list_type
    integer :: num_connection_objects
    type(connection_type), pointer :: first
    type(connection_type), pointer :: last
    type(connection_ptr_type), pointer :: array(:)
  end type connection_list_type
  
  type(connection_list_type), pointer, private :: internal_connection_list, &
                                                  boundary_connection_list

  public :: allocateConnectionLists, createConnection, addConnectionToList, &
            getInternalConnectionList, getBoundaryConnectionList, &
            getNumberOfInternalConnections, getNumberOfBoundaryConnections
  
contains

! ************************************************************************** !
!
! getInternalConnectionList: Returns pointer to internal_connection_list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function getInternalConnectionList()

  implicit none
  
  type(connection_list_type), pointer :: getInternalConnectionList
  getInternalConnectionList => internal_connection_list

end function getInternalConnectionList

! ************************************************************************** !
!
! getBoundaryConnectionList: Returns pointer to boundary_connection_list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function getBoundaryConnectionList()

  implicit none
  
  type(connection_list_type), pointer :: getBoundaryConnectionList
  getBoundaryConnectionList => boundary_connection_list

end function getBoundaryConnectionList

! ************************************************************************** !
!
! getNumberOfBoundaryConnections: Returns pointer to boundary_connection_list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function getNumberOfBoundaryConnections()

  implicit none
  
  integer :: getNumberOfBoundaryConnections
  getNumberOfBoundaryConnections = boundary_connection_list%first%num_connections

end function getNumberOfBoundaryConnections

! ************************************************************************** !
!
! getNumberOfInternalConnections: Returns pointer to internal_connection_list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function getNumberOfInternalConnections()

  implicit none
  
  integer :: getNumberOfInternalConnections
  getNumberOfInternalConnections = internal_connection_list%first%num_connections

end function getNumberOfInternalConnections

! ************************************************************************** !
!
! allocateConnectionLists: Allocates connections lists
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine allocateConnectionLists()

  implicit none
  
  allocate(internal_connection_list)
  call initConnectionList(internal_connection_list)
  allocate(boundary_connection_list)
  call initConnectionList(boundary_connection_list)
  
end subroutine

! ************************************************************************** !
!
! InitConnectionModule: Initializes module variables, lists, arrays.
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine initConnectionList(list)

  implicit none

  type(connection_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_connection_objects = 0

end subroutine InitConnectionList

! ************************************************************************** !
!
! createConnection: Allocates and initializes a new connection
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function createConnection(num_connections,num_dof)

  implicit none
  
  integer :: num_connections
  integer :: num_dof
  
  type(connection_type), pointer :: createConnection
  allocate(createConnection)
  createConnection%id = 0
  createConnection%num_connections = num_connections
  allocate(createConnection%id_up(num_connections))
  allocate(createConnection%id_dn(num_connections))
  allocate(createConnection%dist(-1:3,num_connections))
  allocate(createConnection%area(num_connections))
  allocate(createConnection%velocity(num_dof,num_connections))
  createConnection%id_up = 0
  createConnection%id_dn = 0
  createConnection%dist = 0.d0
  createConnection%area = 0.d0
  createConnection%velocity = 0.d0
  nullify(createConnection%next)

end function createConnection

! ************************************************************************** !
!
! addConnectionToList: Adds a new connection of the module global list of 
!                      connections
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine addConnectionToList(new_connection,list)

  implicit none
  
  type(connection_type), pointer :: new_connection
  type(connection_list_type) :: list
  
  list%num_connection_objects = list%num_connection_objects + 1
  new_connection%id = list%num_connection_objects
  if (.not.associated(list%first)) list%first => new_connection
  if (associated(list%last)) list%last%next => new_connection
  list%last => new_connection
  
end subroutine addConnectionToList

! ************************************************************************** !
!
! convertConnectionListToArray: Creates an array of pointers to the 
!                               connections in the connection list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine convertConnectionListToArray(list)

  implicit none
  
  type(connection_list_type) :: list
    
  integer :: count
  type(connection_type), pointer :: cur_connection
  
  
  allocate(list%array(list%num_connection_objects))
  
  cur_connection => list%first
  do 
    if (.not.associated(cur_connection)) exit
    list%array(cur_connection%id)%ptr => cur_connection
    cur_connection => cur_connection%next
  enddo

end subroutine convertConnectionListToArray

! ************************************************************************** !
!
! destroyConnectionList: Deallocates the module global list and array of regions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine destroyConnectionList(list)

  implicit none
  
  type(connection_list_type) :: list
    
  type(connection_type), pointer :: cur_connection, prev_connection
  
  deallocate(list%array)
  nullify(list%array)
  
  cur_connection => list%first
  do 
    if (.not.associated(cur_connection)) exit
    if (associated(cur_connection%id_up)) deallocate(cur_connection%id_up)
    if (associated(cur_connection%id_dn)) deallocate(cur_connection%id_dn)
    if (associated(cur_connection%dist)) deallocate(cur_connection%dist)
    if (associated(cur_connection%area)) deallocate(cur_connection%area)
    if (associated(cur_connection%velocity)) &
      deallocate(cur_connection%velocity)
    prev_connection => cur_connection
    cur_connection => cur_connection%next
    deallocate(prev_connection)
  enddo
  
  nullify(list%first)
  nullify(list%last)
  list%num_connection_objects = 0

end subroutine destroyConnectionList

end module Connection_module