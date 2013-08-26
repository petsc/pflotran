module Connection_module

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  private

  type, public :: connection_set_type
    PetscInt :: id
    PetscInt :: itype                  ! connection type (boundary, internal, source sink
    PetscInt :: num_connections
    PetscInt :: offset
    PetscInt, pointer :: local(:)      ! 1 if connection is local, 0 if connection is ghosted
    PetscInt, pointer :: id_up(:)      ! list of ids of upwind cells
    PetscInt, pointer :: id_dn(:)      ! list of ids of downwind cells
    PetscInt, pointer :: id_up2(:)     ! list of ids of 2nd upwind cells
    PetscInt, pointer :: id_dn2(:)     ! list of ids of 2nd downwind cells
    PetscReal, pointer :: dist(:,:)    ! list of distance vectors, size(-1:3,num_connections) where
                                       !   -1 = fraction upwind
                                       !   0 = magnitude of distance 
                                       !   1-3 = components of unit vector
    PetscReal, pointer :: intercp(:,:) ! x,y,z location of intercept between the line connecting
                                       ! upwind and downwind cells with the face shared by the cells
    PetscReal, pointer :: area(:)      ! list of areas of faces normal to distance vectors
    PetscReal, pointer :: cntr(:,:)    ! coordinates (1:3, num_connections) of the mass center of the face
    PetscInt, pointer :: face_id(:)    ! list of ids of faces (in local order)
    type(connection_set_type), pointer :: next
  end type connection_set_type


  ! pointer data structure required for making an array of region pointers in F90
  type, public :: connection_set_ptr_type
    type(connection_set_type), pointer :: ptr           ! pointer to the connection_set_type
  end type connection_set_ptr_type 
  
  type, public :: connection_set_list_type
    PetscInt :: num_connection_objects
    type(connection_set_type), pointer :: first
    type(connection_set_type), pointer :: last
    type(connection_set_ptr_type), pointer :: array(:)
  end type connection_set_list_type
  
  public :: ConnectionCreate, ConnectionAddToList, &
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
function ConnectionCreate(num_connections,connection_itype)

  implicit none
  
  PetscInt :: num_connections
  PetscInt :: num_dof
  PetscInt :: connection_itype
  
  type(connection_set_type), pointer :: ConnectionCreate

  type(connection_set_type), pointer :: connection

  allocate(connection)
  connection%id = 0
  connection%itype = connection_itype
  connection%offset = 0
  connection%num_connections = num_connections
  nullify(connection%local)
  nullify(connection%id_up)
  nullify(connection%id_dn)
  nullify(connection%id_up2)
  nullify(connection%id_dn2)
  nullify(connection%face_id)
  nullify(connection%dist)
  nullify(connection%intercp)
  nullify(connection%area)
  nullify(connection%cntr)
  select case(connection_itype)
    case(INTERNAL_CONNECTION_TYPE)
      allocate(connection%id_up(num_connections))
      allocate(connection%id_dn(num_connections))
      allocate(connection%dist(-1:3,num_connections))
      allocate(connection%intercp(1:3,num_connections))
      allocate(connection%area(num_connections))
      allocate(connection%face_id(num_connections))
#ifdef DASVYAT
      allocate(connection%cntr(1:3, num_connections))
      allocate(connection%local(num_connections))
      connection%local = 0
#endif
      connection%id_up = 0
      connection%id_dn = 0
      connection%face_id = 0
      connection%dist = 0.d0
      connection%intercp = 0.d0
      connection%area = 0.d0
    case(BOUNDARY_CONNECTION_TYPE)
      allocate(connection%id_dn(num_connections))
      allocate(connection%dist(-1:3,num_connections))
      allocate(connection%intercp(1:3,num_connections))
      allocate(connection%area(num_connections))
      allocate(connection%face_id(num_connections))
#ifdef DASVYAT
      allocate(connection%cntr(1:3, num_connections))
      connection%cntr = 0.d0
#endif
      connection%id_dn = 0
      connection%dist = 0.d0
      connection%intercp = 0.d0
      connection%area = 0.d0
    case(SRC_SINK_CONNECTION_TYPE,INITIAL_CONNECTION_TYPE)
      allocate(connection%id_dn(num_connections))
#ifdef DASVYAT
      allocate(connection%cntr(1:3, num_connections))
#endif
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
  
  type(connection_set_list_type) :: list

  PetscInt :: ConnectionGetNumberInList
  type(connection_set_type), pointer :: cur_connection_set
  
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
! InitConnectionModule: Initializes module variables, lists, arrays.
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine ConnectionInitList(list)

  implicit none

  type(connection_set_list_type) :: list
  
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
subroutine ConnectionAddToList(new_connection_set,list)

  implicit none
  
  type(connection_set_type), pointer :: new_connection_set
  type(connection_set_list_type) :: list
  
  list%num_connection_objects = list%num_connection_objects + 1
  new_connection_set%id = list%num_connection_objects
  if (.not.associated(list%first)) list%first => new_connection_set
  if (associated(list%last)) list%last%next => new_connection_set
  list%last => new_connection_set
  
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
  
  type(connection_set_list_type) :: list
    
  PetscInt :: count
  type(connection_set_type), pointer :: cur_connection_set
  
  
  allocate(list%array(list%num_connection_objects))
  
  cur_connection_set => list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    list%array(cur_connection_set%id)%ptr => cur_connection_set
    cur_connection_set => cur_connection_set%next
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
  
  type(connection_set_type), pointer :: connection
  
  if (.not.associated(connection)) return
  
  if (associated(connection%local)) deallocate(connection%local)
  nullify(connection%local)
  if (associated(connection%id_up)) deallocate(connection%id_up)
  nullify(connection%id_up)
  if (associated(connection%id_dn)) deallocate(connection%id_dn)
  nullify(connection%id_dn)
  if (associated(connection%id_up2)) deallocate(connection%id_up2)
  nullify(connection%id_up2)
  if (associated(connection%id_dn2)) deallocate(connection%id_dn2)
  nullify(connection%id_dn2)
  if (associated(connection%face_id)) deallocate(connection%face_id)
  nullify(connection%face_id)
  if (associated(connection%dist)) deallocate(connection%dist)
  nullify(connection%dist)
  if (associated(connection%intercp)) deallocate(connection%intercp)
  nullify(connection%intercp)
  if (associated(connection%area)) deallocate(connection%area)
  nullify(connection%area)
  if (associated(connection%cntr)) deallocate(connection%cntr)
  nullify(connection%cntr)
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
  
  type(connection_set_list_type), pointer :: list
    
  type(connection_set_type), pointer :: cur_connection_set, prev_connection_set
  
  if (.not.associated(list)) return
  
  if (associated(list%array)) deallocate(list%array)
  nullify(list%array)
  
  cur_connection_set => list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    prev_connection_set => cur_connection_set
    cur_connection_set => cur_connection_set%next
    call ConnectionDestroy(prev_connection_set)
  enddo
  
  nullify(list%first)
  nullify(list%last)
  list%num_connection_objects = 0
  
  deallocate(list)
  nullify(list)

end subroutine ConnectionDestroyList

end module Connection_module
