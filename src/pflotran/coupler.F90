module Coupler_module
 
  use Condition_module
  use Connection_module
  use Region_module
 
  implicit none

  private
 
#include "definitions.h"
 
  type, public :: coupler_type
    integer :: id                                       ! id of coupler
    integer :: itype                                    ! integer defining type
    character(len=MAXWORDLENGTH) :: ctype               ! character string definign type
    character(len=MAXWORDLENGTH) :: condition_name      ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: region_name         ! character string defining name of region to be applied
    integer :: icondition                               ! id of condition in condition array/list
    integer :: iregion                                  ! id of region in region array/list
    integer :: iface                                    ! for structured grids only
    integer, pointer :: aux_int_var(:,:)                ! auxilliary array for integer value
    real*8, pointer :: aux_real_var(:,:)                ! auxilliary array for real values
    type(condition_type), pointer :: condition          ! pointer to condition in condition array/list
    type(region_type), pointer :: region                ! pointer to region in region array/list
    type(connection_type), pointer :: connection        ! pointer to an array/list of connections
    type(coupler_type), pointer :: next                 ! pointer to next coupler
  end type coupler_type
  
  type, public :: coupler_ptr_type
    type(coupler_type), pointer :: ptr
  end type coupler_ptr_type
    
  type, public :: coupler_list_type
    integer :: num_couplers
    type(coupler_type), pointer :: first
    type(coupler_type), pointer :: last
    type(coupler_ptr_type), pointer :: array(:)    
  end type coupler_list_type
  
  integer, save :: num_couplers = 0
  
  public :: CouplerCreate, CouplerDestroy, CouplerInitList, CouplerAddToList, &
            CouplerRead, CouplerDestroyList, CouplerGetNumConnectionsInList, &
            CouplerUpdateAuxVars

  
  interface CouplerCreate
    module procedure CouplerCreate1
    module procedure CouplerCreate2
  end interface
    
contains

! ************************************************************************** !
!
! CouplerCreate: Creates a coupler
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function CouplerCreate1()

  implicit none

  type(coupler_type), pointer :: CouplerCreate1
  
  type(coupler_type), pointer :: coupler
  
  allocate(coupler)
  coupler%id = 0
  coupler%itype = BOUNDARY_COUPLER_TYPE
  coupler%ctype = "boundary"
  coupler%condition_name = ""
  coupler%region_name = ""
  coupler%icondition = 0
  coupler%iregion = 0
  coupler%iface = 0
  nullify(coupler%aux_int_var)
  nullify(coupler%aux_real_var)
  nullify(coupler%condition)
  nullify(coupler%region)
  nullify(coupler%connection)
  nullify(coupler%next)
  
  num_couplers = num_couplers + 1
  coupler%id = num_couplers
  
  CouplerCreate1 => coupler

end function CouplerCreate1

! ************************************************************************** !
!
! CouplerCreate2: Creates a coupler
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function CouplerCreate2(itype)

  implicit none

  integer :: itype
  
  type(coupler_type), pointer :: CouplerCreate2
  
  type(coupler_type), pointer :: coupler
  
  coupler => CouplerCreate1()
  coupler%itype = itype
  select case(itype)
    case(INITIAL_COUPLER_TYPE)
      coupler%ctype = 'initial'
    case(BOUNDARY_COUPLER_TYPE)
      coupler%ctype = 'boundary'
    case(SRC_SINK_COUPLER_TYPE)
      coupler%ctype = 'source_sink'
  end select

  CouplerCreate2 => coupler

end function CouplerCreate2

! ************************************************************************** !
!
! CouplerInitList: Initializes a coupler list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine CouplerInitList(list)

  implicit none

  type(coupler_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_couplers = 0

end subroutine CouplerInitList

! ************************************************************************** !
!
! CouplerRead: Reads a coupler from the input file
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine CouplerRead(coupler,fid)

  use Fileio_module
  
  implicit none
  
  type(coupler_type) :: coupler
  integer :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  integer :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg('keyword','COUPLER', ierr)   
      
    select case(trim(word))
    
      case('REGION')
        call fiReadWord(string,coupler%region_name,.true.,ierr)
      case('CONDITION')
        call fiReadWord(string,coupler%condition_name,.true.,ierr)
      case('TYPE')
        call fiReadWord(string,coupler%ctype,.true.,ierr)
        call fiCharsToLower(coupler%ctype,len_trim(coupler%ctype))
        select case(trim(coupler%ctype))
          case('initial')
            coupler%itype = INITIAL_COUPLER_TYPE
          case('boundary')
            coupler%itype = BOUNDARY_COUPLER_TYPE
          case('source_sink')
            coupler%itype = SRC_SINK_COUPLER_TYPE
          case default
            print *, 'ERROR: TYPE option (', trim(coupler%ctype), ') not recognized.'
            stop
        end select    
      case('FACE')
        call fiReadWord(string,word,.true.,ierr)
        call fiCharsToUpper(word,len_trim(word))
        select case(word)
          case('WEST')
            coupler%iface = 1
          case('EAST')
            coupler%iface = 2
          case('SOUTH')
            coupler%iface = 3
          case('NORTH')
            coupler%iface = 4
          case('BOTTOM')
            coupler%iface = 5
          case('TOP')
            coupler%iface = 6
          case default
            print *, 'ERROR: FACE option (', trim(word), ') not recognized.'
            stop
        end select
      case('END')
        exit
    end select 
  
  enddo  

end subroutine CouplerRead

! ************************************************************************** !
!
! CouplerUpdateAuxVars: Updates auxilliary variables associated with a coupler
! author: Glenn Hammond
! date: 11/19/07
!
! ************************************************************************** !
subroutine CouplerUpdateAuxVars(coupler,option)

  use Option_module

  implicit none
  
  type(coupler_type) :: coupler
  type(option_type) :: option
  
  integer :: idof, num_connections
  
  num_connections = coupler%connection%num_connections
  select case(option%imode)
    case(RICHARDS_MODE,MPH_MODE)
      coupler%aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
        coupler%condition%iphase
  end select
  
  do idof = 1, option%ndof
  
    select case(coupler%condition%itype(idof))
      case(DIRICHLET_BC,NEUMANN_BC,MASS_RATE)
        coupler%aux_real_var(idof,1:num_connections) = &
          coupler%condition%cur_value(idof)
      case(HYDROSTATIC_BC)
      ! call hydrostatic subroutines
    end select 
  enddo

end subroutine CouplerUpdateAuxVars

! ************************************************************************** !
!
! CouplerAddToList: Adds a new coupler to a coupler list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine CouplerAddToList(new_coupler,list)

  implicit none
  
  type(coupler_type), pointer :: new_coupler
  type(coupler_list_type) :: list
  
  list%num_couplers = list%num_couplers + 1
  new_coupler%id = list%num_couplers
  if (.not.associated(list%first)) list%first => new_coupler
  if (associated(list%last)) list%last%next => new_coupler
  list%last => new_coupler
  
end subroutine CouplerAddToList

! ************************************************************************** !
!
! CouplerGetNumConnectionsInList: Returns the number of connections associated
!                                 with all couplers in the list
! author: Glenn Hammond
! date: 11/19/07
!
! ************************************************************************** !
function CouplerGetNumConnectionsInList(list)

  implicit none
  
  type(coupler_list_type) :: list
  
  integer :: CouplerGetNumConnectionsInList
  type(coupler_type), pointer :: coupler
  
  CouplerGetNumConnectionsInList = 0
  coupler => list%first
  
  do
    if (.not.associated(coupler)) exit
    CouplerGetNumConnectionsInList = CouplerGetNumConnectionsInList + &
                                     coupler%connection%num_connections
    coupler => coupler%next
  enddo

end function CouplerGetNumConnectionsInList

! ************************************************************************** !
!
! CouplerDestroyList: Deallocates a list of couplers
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine CouplerDestroyList(coupler_list)

  implicit none
  
  type(coupler_list_type), pointer :: coupler_list
  
  type(coupler_type), pointer :: coupler, prev_coupler
  
  
  coupler => coupler_list%first
  do 
    if (.not.associated(coupler)) exit
    prev_coupler => coupler
    coupler => coupler%next
    call CouplerDestroy(prev_coupler)
  enddo
  
  coupler_list%num_couplers = 0
  nullify(coupler_list%first)
  nullify(coupler_list%last)
  if (associated(coupler_list%array)) deallocate(coupler_list%array)
  nullify(coupler_list%array)
  
  deallocate(coupler_list)
  nullify(coupler_list)

end subroutine CouplerDestroyList
  
! ************************************************************************** !
!
! CouplerDestroy: Destroys a coupler
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine CouplerDestroy(coupler)

  implicit none
  
  type(coupler_type), pointer :: coupler
  
  if (.not.associated(coupler)) return
  
  ! since the below are simply pointers to objects in list that have already
  ! or will be deallocated from the list, nullify instead of destroying
  
  nullify(coupler%condition)     ! since these are simply pointers to 
  nullify(coupler%region)        ! conditoins in list, nullify

  if (associated(coupler%aux_int_var)) deallocate(coupler%aux_int_var)
  if (associated(coupler%aux_real_var)) deallocate(coupler%aux_real_var)

  call ConnectionDestroy(coupler%connection)
  nullify(coupler%connection)
  
  deallocate(coupler)
  nullify(coupler)

end subroutine CouplerDestroy

end module Coupler_module
