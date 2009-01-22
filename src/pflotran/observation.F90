module Observation_module

  use Region_module
  use Connection_module
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public :: observation_type
    ! all added variables must be included in ObservationCreateFromObservation
    PetscInt :: id
    PetscInt :: itype
    PetscTruth :: print_velocities
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: linkage_name
    type(connection_set_type), pointer :: connection_set
    type(region_type), pointer :: region
    type(observation_type), pointer :: next
  end type observation_type
  
  type, public :: observation_list_type
    PetscInt :: num_observations
    type(observation_type), pointer :: first
    type(observation_type), pointer :: last
    type(observation_type), pointer :: array(:)
  end type observation_list_type

  public :: ObservationCreate, ObservationDestroy, ObservationRead, &
            ObservationAddToList, ObservationInitList, ObservationDestroyList, &
            ObservationGetPtrFromList, ObservationRemoveFromList

  interface ObservationCreate
    module procedure ObservationCreate1
    module procedure ObservationCreateFromObservation
  end interface
    
contains

! ************************************************************************** !
!
! ObservationCreate1: Create object that stores observation regions
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function ObservationCreate1()

  implicit none
  
  type(observation_type), pointer :: ObservationCreate1
  
  type(observation_type), pointer :: observation
  
  allocate(observation)
  
  observation%name = ""
  observation%linkage_name = ""
  observation%id = 0
  observation%itype = OBSERVATION_SCALAR
  observation%print_velocities = PETSC_FALSE
  nullify(observation%region)
  nullify(observation%next)
  
  ObservationCreate1 => observation

end function ObservationCreate1

! ************************************************************************** !
!
! ObservationCreate: Create object that stores observation regions
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function ObservationCreateFromObservation(observation)

  implicit none
  
  type(observation_type), pointer :: ObservationCreateFromObservation
  type(observation_type), pointer :: observation

  type(observation_type), pointer :: new_observation
  
  new_observation => ObservationCreate1()
  
  new_observation%name = observation%name
  new_observation%linkage_name = observation%linkage_name
  new_observation%id = observation%id
  new_observation%itype = observation%itype
  new_observation%print_velocities = observation%print_velocities
  ! keep these null for now to catch bugs
  nullify(new_observation%region)
  nullify(new_observation%next)
  
  ObservationCreateFromObservation => new_observation

end function ObservationCreateFromObservation

! ************************************************************************** !
!
! ObservationRead: Reads observation data from the input file
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine ObservationRead(observation,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none
  
  type(observation_type) :: observation
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','OBSERVATION')   
      
    select case(trim(keyword))
    
      case('BOUNDARY_CONDITION')
        call InputReadWord(input,option,observation%linkage_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'boundary condition name','OBSERVATION')
        option%store_solute_fluxes = PETSC_TRUE
        observation%itype = OBSERVATION_FLUX
      case('REGION')
        call InputReadWord(input,option,observation%linkage_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'region name','OBSERVATION')
        observation%itype = OBSERVATION_SCALAR
      case('VELOCITY')
        observation%print_velocities = PETSC_TRUE
      case default
        option%io_buffer = 'Keyword (' // trim(keyword) // &
                           ') not recognized under' // &
                           ' OBSERVATION.'
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine ObservationRead

! ************************************************************************** !
!
! ObservationInitList: Initializes a observation list
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine ObservationInitList(list)

  implicit none

  type(observation_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_observations = 0

end subroutine ObservationInitList

! ************************************************************************** !
!
! ObservationAddToList: Adds a new observation to a observation list
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine ObservationAddToList(new_observation,list)

  implicit none
  
  type(observation_type), pointer :: new_observation
  type(observation_list_type) :: list
  
  list%num_observations = list%num_observations + 1
  new_observation%id = list%num_observations
  if (.not.associated(list%first)) list%first => new_observation
  if (associated(list%last)) list%last%next => new_observation
  list%last => new_observation
  
end subroutine ObservationAddToList

! ************************************************************************** !
!
! ObservationRemoveFromList: Removes a observation from a observation list
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine ObservationRemoveFromList(observation,list)

  implicit none
  
  type(observation_type), pointer :: observation
  type(observation_list_type) :: list
  
  type(observation_type), pointer :: cur_observation, prev_observation
  
  cur_observation => list%first
  nullify(prev_observation)
  
  do
    if (.not.associated(cur_observation)) exit
    if (associated(cur_observation,observation)) then
      if (associated(prev_observation)) then
        prev_observation%next => cur_observation%next
      else
        list%first => cur_observation%next
      endif
      if (.not.associated(cur_observation%next)) then
        list%last => prev_observation
      endif
      list%num_observations = list%num_observations-1
      call ObservationDestroy(cur_observation)
      return
    endif
    prev_observation => cur_observation
    cur_observation => cur_observation%next
  enddo
  
end subroutine ObservationRemoveFromList

! ************************************************************************** !
!
! ObservationGetPtrFromList: Returns a pointer to the observation matching &
!                             observation_name
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function ObservationGetPtrFromList(observation_name,observation_list)

  use String_module

  implicit none
  
  type(observation_type), pointer :: ObservationGetPtrFromList
  character(len=MAXWORDLENGTH) :: observation_name
  type(observation_list_type) :: observation_list
 
  PetscInt :: length
  type(observation_type), pointer :: observation
    
  nullify(ObservationGetPtrFromList)
  observation => observation_list%first
  
  do 
    if (.not.associated(observation)) exit
    length = len_trim(observation_name)
    if (length == len_trim(observation%name) .and. &
        StringCompare(observation%name,observation_name, &
                        length)) then
      ObservationGetPtrFromList => observation
      return
    endif
    observation => observation%next
  enddo
  
end function ObservationGetPtrFromList

! ************************************************************************** !
!
! ObservationDestroyList: Deallocates a list of observations
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine ObservationDestroyList(observation_list)

  implicit none
  
  type(observation_list_type), pointer :: observation_list
  
  type(observation_type), pointer :: observation, prev_observation
  
  if (.not.associated(observation_list)) return
  
  observation => observation_list%first
  do 
    if (.not.associated(observation)) exit
    prev_observation => observation
    observation => observation%next
    call ObservationDestroy(prev_observation)
  enddo
  
  observation_list%num_observations = 0
  nullify(observation_list%first)
  nullify(observation_list%last)
  if (associated(observation_list%array)) deallocate(observation_list%array)
  nullify(observation_list%array)
  
  deallocate(observation_list)
  nullify(observation_list)

end subroutine ObservationDestroyList

! ************************************************************************** !
!
! ObservationDestroy: Deallocates a observation
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine ObservationDestroy(observation)

  implicit none
  
  type(observation_type), pointer :: observation
  
  PetscInt :: i
  
  if (.not.associated(observation)) return
  
  nullify(observation%region)
  nullify(observation%connection_set)
  deallocate(observation)
  nullify(observation)

end subroutine ObservationDestroy

end module Observation_module
