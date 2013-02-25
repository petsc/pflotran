module Coupler_module
 
  use Condition_module
  use Connection_module
  use Region_module
 
  implicit none

  private
 
#include "definitions.h"

  ! coupler types
  PetscInt, parameter, public :: INITIAL_COUPLER_TYPE = 1
  PetscInt, parameter, public :: BOUNDARY_COUPLER_TYPE = 2
  PetscInt, parameter, public :: SRC_SINK_COUPLER_TYPE = 3
  PetscInt, parameter, public :: COUPLER_IPHASE_INDEX = 1

  type, public :: coupler_type
    PetscInt :: id                                      ! id of coupler
    character(len=MAXWORDLENGTH) :: name                ! name of coupler
    PetscInt :: itype                                   ! integer defining type
    character(len=MAXWORDLENGTH) :: ctype               ! character string defining type
    character(len=MAXWORDLENGTH) :: flow_condition_name ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: tran_condition_name ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: region_name         ! character string defining name of region to be applied
    PetscInt :: iflow_condition                         ! id of condition in condition array/list
    PetscInt :: itran_condition                         ! id of condition in condition array/list
    PetscInt :: iregion                                 ! id of region in region array/list
    PetscInt :: iface                                   ! for structured grids only
    PetscInt, pointer :: flow_aux_int_var(:,:)          ! auxiliary array for integer value
    PetscReal, pointer :: flow_aux_real_var(:,:)        ! auxiliary array for real values
    type(flow_condition_type), pointer :: flow_condition     ! pointer to condition in condition array/list
    type(tran_condition_type), pointer :: tran_condition     ! pointer to condition in condition array/list
    type(region_type), pointer :: region                ! pointer to region in region array/list
    type(connection_set_type), pointer :: connection_set ! pointer to an array/list of connections
    PetscInt, pointer :: faces_set(:)                   ! ids of the elements of array grid%faces. Include local and ghosted faces. Doesn't require additional allocation of connections. Implemented in MIMETIC mode.
    PetscInt :: numfaces_set
    type(coupler_type), pointer :: next                 ! pointer to next coupler
  end type coupler_type
  
  type, public :: coupler_ptr_type
    type(coupler_type), pointer :: ptr
  end type coupler_ptr_type
    
  type, public :: coupler_list_type
    PetscInt :: num_couplers
    type(coupler_type), pointer :: first
    type(coupler_type), pointer :: last
    type(coupler_ptr_type), pointer :: array(:)    
  end type coupler_list_type
  
  public :: CouplerCreate, CouplerDestroy, CouplerInitList, CouplerAddToList, &
            CouplerRead, CouplerDestroyList, CouplerGetNumConnectionsInList, &
            CouplerListComputeConnections, CouplerGetPtrFromList,&
            CouplerAssignBCtoCells, CouplerGetNumBoundConnectionsInListMFD

  
  interface CouplerCreate
    module procedure CouplerCreate1
    module procedure CouplerCreate2
    module procedure CouplerCreateFromCoupler
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
  coupler%name = ''
  coupler%itype = BOUNDARY_COUPLER_TYPE
  coupler%ctype = "boundary"
  coupler%flow_condition_name = ""
  coupler%tran_condition_name = ""
  coupler%region_name = ""
  coupler%iflow_condition = 0
  coupler%itran_condition = 0
  coupler%iregion = 0
  coupler%iface = 0
  nullify(coupler%flow_aux_int_var)
  nullify(coupler%flow_aux_real_var)
  nullify(coupler%flow_condition)
  nullify(coupler%tran_condition)
  nullify(coupler%region)
  nullify(coupler%connection_set)
  nullify(coupler%faces_set)
  nullify(coupler%next)
  
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

  PetscInt :: itype
  
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
! CouplerCreateFromCoupler: Creates a coupler
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function CouplerCreateFromCoupler(coupler)

  implicit none
  
  type(coupler_type), pointer :: coupler
  
  type(coupler_type), pointer :: CouplerCreateFromCoupler
  type(coupler_type), pointer :: new_coupler

  new_coupler => CouplerCreate1()

  new_coupler%id = coupler%id
  new_coupler%name = coupler%name
  new_coupler%itype = coupler%itype
  new_coupler%ctype = coupler%ctype
  new_coupler%flow_condition_name = coupler%flow_condition_name
  new_coupler%tran_condition_name = coupler%tran_condition_name
  new_coupler%region_name = coupler%region_name
  new_coupler%iflow_condition = coupler%iflow_condition
  new_coupler%itran_condition = coupler%itran_condition
  new_coupler%iregion = coupler%iregion
  new_coupler%iface = coupler%iface

  ! these must remain null  
  nullify(coupler%flow_condition)
  nullify(coupler%tran_condition)
  nullify(coupler%region)
  nullify(coupler%flow_aux_int_var)
  nullify(coupler%flow_aux_real_var)
  nullify(coupler%connection_set)
  nullify(coupler%next)

  CouplerCreateFromCoupler => new_coupler

end function CouplerCreateFromCoupler

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
subroutine CouplerRead(coupler,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(coupler_type) :: coupler
  type(input_type) :: input
  
  character(len=MAXWORDLENGTH) :: word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','COUPLER')   
    call StringToUpper(word)      
    
    select case(trim(word))
    
      case('REGION','SURF_REGION')
        call InputReadWord(input,option,coupler%region_name,PETSC_TRUE)
      case('FLOW_CONDITION','SURF_FLOW_CONDITION')
        call InputReadWord(input,option,coupler%flow_condition_name,PETSC_TRUE)
      case('TRANSPORT_CONDITION')
        call InputReadWord(input,option,coupler%tran_condition_name,PETSC_TRUE)
      case default
        option%io_buffer = 'Coupler card (' // trim(word) // ') not recognized.'
        call printErrMsg(option)        
    end select 
  
  enddo  

end subroutine CouplerRead

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
! CouplerListComputeConnections: computes connectivity for a list of couplers
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine CouplerListComputeConnections(grid,option,coupler_list)

  use Option_module
  use Grid_module
  
  implicit none
 
  type(grid_type) :: grid
  type(option_type) :: option
  type(coupler_list_type), pointer :: coupler_list
  
  type(coupler_type), pointer :: coupler
  PetscInt :: offset
  
  if (.not.associated(coupler_list)) return
  
  offset = 0
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit 
    if (grid%itype == STRUCTURED_GRID_MIMETIC .and. &
        (coupler%itype == INITIAL_COUPLER_TYPE .or. &
         coupler%itype == BOUNDARY_COUPLER_TYPE)) then  
       call CouplerComputeConnections(grid,option,coupler)
       call CouplerComputeConnectionsFaces(grid,option,coupler)      
       call CouplerAssignBCtoCells(grid,option,coupler)
 
    else
       call CouplerComputeConnections(grid,option,coupler)
    end if
    if (associated(coupler%connection_set)) then
      coupler%connection_set%offset = offset
      offset = offset + coupler%connection_set%num_connections
    endif
    coupler => coupler%next
  enddo

end subroutine CouplerListComputeConnections

! ************************************************************************** !
!
! CouplerComputeConnections: computes connectivity coupler to a grid
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine CouplerComputeConnections(grid,option,coupler)

  use Connection_module
  use Option_module
  use Region_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Explicit_module, only : ExplicitUGridSetBoundaryConnect, &
                                           ExplicitUGridSetConnections
  
  implicit none
 
  type(grid_type) :: grid
  type(option_type) :: option
  type(coupler_type), pointer :: coupler_list
  
  PetscInt :: iconn
  PetscInt :: cell_id_local, cell_id_ghosted
  PetscInt :: connection_itype
  PetscInt :: iface
  type(connection_set_type), pointer :: connection_set
  type(region_type), pointer :: region
  type(coupler_type), pointer :: coupler
  PetscBool :: nullify_connection_set
  PetscErrorCode :: ierr

  if (.not.associated(coupler)) return
  
  nullify_connection_set = PETSC_FALSE
  select case(coupler%itype)
    case(INITIAL_COUPLER_TYPE)
      if (associated(coupler%flow_condition)) then
        if (associated(coupler%flow_condition%pressure)) then
          if (coupler%flow_condition%pressure%itype /= HYDROSTATIC_BC .and. &
              coupler%flow_condition%pressure%itype /= SEEPAGE_BC .and. &
              coupler%flow_condition%pressure%itype /= CONDUCTANCE_BC) then
            nullify_connection_set = PETSC_TRUE
          endif
        else if (associated(coupler%flow_condition%concentration)) then
          ! need to calculate connection set
        endif
      else
        nullify_connection_set = PETSC_TRUE
      endif
      connection_itype = INITIAL_CONNECTION_TYPE
    case(SRC_SINK_COUPLER_TYPE)
      connection_itype = SRC_SINK_CONNECTION_TYPE
    case(BOUNDARY_COUPLER_TYPE)
      connection_itype = BOUNDARY_CONNECTION_TYPE
  end select
  
  if (nullify_connection_set) then
    nullify(coupler%connection_set)
    return
  endif
  
  region => coupler%region

  select case(grid%itype)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      if (associated(region%explicit_faceset)) then
        connection_set => &
          ExplicitUGridSetBoundaryConnect(grid%unstructured_grid% &
                                            explicit_grid, &
                                          region%cell_ids, &
                                     region%explicit_faceset%face_centroids, &
                                     region%explicit_faceset%face_areas, &
                                     option)
      else
        connection_set => &
          ExplicitUGridSetConnections(grid%unstructured_grid% &
                                        explicit_grid, &
                                      region%cell_ids, &
                                      connection_itype,option)
      endif
    case default
      connection_set => ConnectionCreate(region%num_cells,connection_itype)
    
      ! if using higher order advection, allocate associated arrays
      if (option%itranmode == EXPLICIT_ADVECTION .and. &
          option%tvd_flux_limiter /= 1 .and. &  ! 1 = upwind
          connection_set%itype == BOUNDARY_CONNECTION_TYPE) then
        ! connections%id_up2 should remain null as it will not be used
        allocate(connection_set%id_dn2(size(connection_set%id_dn)))
        connection_set%id_dn2 = 0
      endif  

      iface = coupler%iface
      do iconn = 1,region%num_cells
    
        cell_id_local = region%cell_ids(iconn)
        if (associated(region%faces)) iface = region%faces(iconn)
    
        connection_set%id_dn(iconn) = cell_id_local

        call GridPopulateConnection(grid,connection_set,iface,iconn, &
                                    cell_id_local,option)
      enddo
  end select

  coupler%connection_set => connection_set
  nullify(connection_set)
 
end subroutine CouplerComputeConnections

! ************************************************************************** !
!
! CouplerComputeConnections: computes connectivity coupler to a grid
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine CouplerComputeConnectionsFaces(grid,option,coupler)

  use Connection_module
  use Option_module
  use Region_module
  use Grid_module
  use MFD_Aux_module
  
  implicit none
 
  type(grid_type) :: grid
  type(option_type) :: option
  type(coupler_type), pointer :: coupler
  
#ifdef DASVYAT
  PetscInt :: iconn, reg_numconn
  PetscInt :: cell_id_local, cell_id_ghosted, face_id_ghosted
  PetscInt :: face_id_local
  PetscInt :: connection_itype
  PetscInt :: iface, icell
  type(connection_set_type), pointer :: connection_set
  type(region_type), pointer :: region

  type(mfd_type), pointer :: mfd_aux
  type(mfd_auxvar_type), pointer :: aux_var

  PetscErrorCode :: ierr
  PetscInt, pointer :: local_faces(:)
  PetscInt :: conn_id, stride, iface_type
  type(connection_set_type), pointer :: conn_set_ptr
  PetscBool :: bnd_face

  if (.not.associated(coupler)) return
  

  select case(coupler%itype)
    case(INITIAL_COUPLER_TYPE)
      if (associated(coupler%flow_condition)) then
        if (coupler%flow_condition%pressure%itype /= HYDROSTATIC_BC .and. &
            coupler%flow_condition%pressure%itype /= SEEPAGE_BC .and. &
            coupler%flow_condition%pressure%itype /= CONDUCTANCE_BC) then
          nullify(coupler%connection_set)
          return
        endif
      else
        nullify(coupler%connection_set)
        return
      endif
!      connection_itype = INITIAL_CONNECTION_TYPE
      connection_itype = INTERNAL_CONNECTION_TYPE
    case(SRC_SINK_COUPLER_TYPE)
!      connection_itype = SRC_SINK_CONNECTION_TYPE
      connection_itype = INTERNAL_CONNECTION_TYPE
    case(BOUNDARY_COUPLER_TYPE)
      connection_itype = BOUNDARY_CONNECTION_TYPE
  end select
  
  region => coupler%region
  mfd_aux => grid%MFD
 
  coupler%numfaces_set = 0
 

  allocate(local_faces(grid%nlmax_faces))

  local_faces = 0


  if (coupler%itype == BOUNDARY_COUPLER_TYPE) then
    coupler%numfaces_set = region%num_cells
  else 
    do icell = 1, region%num_cells
      
      cell_id_local = region%cell_ids(icell)
  
      aux_var => mfd_aux%aux_vars(cell_id_local)
  
 
      do iface = 1,aux_var%numfaces
        face_id_ghosted = aux_var%face_id_gh(iface)
        face_id_local = grid%fG2L(face_id_ghosted)
 
        if (face_id_local > 0) then
            if (local_faces(face_id_local)==0) then
              coupler%numfaces_set = coupler%numfaces_set + 1
              local_faces(face_id_local) = 1
            end if
        end if
      end do
    end do
  end if 

 allocate(coupler%faces_set(coupler%numfaces_set))

!    connection_set => ConnectionCreate(coupler%numfaces_set, &
!                                     connection_itype)
!	stop
! else
    connection_set => ConnectionCreate(ZERO_INTEGER, &
                                      connection_itype)
! end if
 local_faces = 0
  

 iconn = 1 
 do icell = 1, region%num_cells
    
    cell_id_local = region%cell_ids(icell)

    aux_var => mfd_aux%aux_vars(cell_id_local)

    do iface = 1,aux_var%numfaces
      face_id_ghosted = aux_var%face_id_gh(iface)
      face_id_local = grid%fG2L(face_id_ghosted)

      if (coupler%itype == BOUNDARY_COUPLER_TYPE) then
        conn_set_ptr => grid%faces(face_id_ghosted)%conn_set_ptr
        conn_id = grid%faces(face_id_ghosted)%id
 
!        write(*,*) "cell ", cell_id_local, "face", face_id_ghosted       
 

        if (conn_set_ptr%itype /= BOUNDARY_CONNECTION_TYPE) cycle
        if (associated(region%faces)) iface_type = region%faces(icell)
        bnd_face = PETSC_FALSE
        select case (iface_type)
             case(WEST_FACE)
               if ((conn_set_ptr%dist(1,conn_id) == 1).and.&
                   (conn_set_ptr%dist(2,conn_id) == 0).and.&
                   (conn_set_ptr%dist(3,conn_id) == 0)) bnd_face = PETSC_TRUE
             case(EAST_FACE)
               if ((conn_set_ptr%dist(1,conn_id) == -1).and.&
                   (conn_set_ptr%dist(2,conn_id) == 0).and.&
                   (conn_set_ptr%dist(3,conn_id) == 0)) bnd_face = PETSC_TRUE
             case(SOUTH_FACE)
               if ((conn_set_ptr%dist(1,conn_id) == 0).and.&
                   (conn_set_ptr%dist(2,conn_id) == 1).and.&
                   (conn_set_ptr%dist(3,conn_id) == 0)) bnd_face = PETSC_TRUE
             case(NORTH_FACE)
               if ((conn_set_ptr%dist(1,conn_id) == 0).and.&
                   (conn_set_ptr%dist(2,conn_id) == -1).and.&
                   (conn_set_ptr%dist(3,conn_id) == 0)) bnd_face = PETSC_TRUE
             case(TOP_FACE)
               if (grid%structured_grid%invert_z_axis) then
                if ((conn_set_ptr%dist(1,conn_id) == 0).and.&
                    (conn_set_ptr%dist(2,conn_id) == 0).and.&
                    (conn_set_ptr%dist(3,conn_id) == 1)) bnd_face = PETSC_TRUE
               else
                if ((conn_set_ptr%dist(1,conn_id) == 0).and.&
                    (conn_set_ptr%dist(2,conn_id) == 0).and.&
                    (conn_set_ptr%dist(3,conn_id) == -1)) bnd_face = PETSC_TRUE
               end if
             case(BOTTOM_FACE)
               if (grid%structured_grid%invert_z_axis) then
                if ((conn_set_ptr%dist(1,conn_id) == 0).and.&
                    (conn_set_ptr%dist(2,conn_id) == 0).and.&
                    (conn_set_ptr%dist(3,conn_id) == -1)) bnd_face = PETSC_TRUE
               else
                if ((conn_set_ptr%dist(1,conn_id) == 0).and.&
                    (conn_set_ptr%dist(2,conn_id) == 0).and.&
                    (conn_set_ptr%dist(3,conn_id) == 1)) bnd_face = PETSC_TRUE
               end if
        end select


       if (bnd_face)   then
!          write(*,*) face_id_ghosted, conn_set_ptr%cntr(1,conn_id), conn_set_ptr%cntr(2,conn_id),conn_set_ptr%cntr(3,conn_id)
           if (face_id_local > 0) then
              coupler%faces_set(iconn) = face_id_ghosted
              iconn = iconn + 1
           end if
       end if

      else 

        if (face_id_local > 0) then
            if (local_faces(face_id_local)==0) then
  
              coupler%faces_set(iconn) = face_id_ghosted
              iconn = iconn + 1
              local_faces(face_id_local) = 1
  
            end if
        end if

      end if
    end do
 end do 

  select case(coupler%itype)
    case(INITIAL_COUPLER_TYPE)
      connection_set%itype = INITIAL_CONNECTION_TYPE
    case(SRC_SINK_COUPLER_TYPE)
      connection_set%itype = SRC_SINK_CONNECTION_TYPE
    case(BOUNDARY_COUPLER_TYPE)
      connection_set%itype = BOUNDARY_CONNECTION_TYPE
  end select
   
!  coupler%connection_set => connection_set

  deallocate(local_faces)

#endif
  
!  nullify(connection_set)


!  stop
 
end subroutine CouplerComputeConnectionsFaces




subroutine CouplerAssignBCtoCells(grid,option,coupler)

  use Connection_module
  use Option_module
  use Region_module
  use Grid_module
  use MFD_Aux_module

  
  implicit none
 
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  type(option_type) :: option
  type(coupler_type), pointer :: coupler

#ifdef DASVYAT
  
  PetscInt :: iconn, reg_numconn
  PetscInt :: cell_id_local, cell_id_ghosted, face_id_ghosted
  PetscInt :: face_id_local
  PetscInt :: connection_itype
  PetscInt :: iface, icell
  type(connection_set_type), pointer :: connection_set
  type(region_type), pointer :: region
  type(mfd_type), pointer :: mfd_aux
  PetscErrorCode :: ierr
  PetscInt, pointer :: local_faces(:)
  type(mfd_auxvar_type), pointer :: aux_var
  PetscInt :: conn_id, stride, iface_type, e2n_size
  PetscScalar, pointer :: e2n_local(:)
  type(connection_set_type), pointer :: conn_set_ptr
  PetscBool :: bnd_face


  stride = 6 !hex only

  if (.not.associated(coupler)) return
  
  if (coupler%itype /= BOUNDARY_COUPLER_TYPE) return 


  call VecGetSize(grid%e2n, e2n_size, ierr)

  if (e2n_size > 0) then 
      call VecGetArrayF90(grid%e2n, e2n_local, ierr)
  end if

  region => coupler%region

  if (.not.associated(region)) then
      option%io_buffer = '.not.associated(region)'
      call printErrMsg(option) 
      stop
  end if

  mfd_aux => grid%MFD

  if (.not.associated(mfd_aux)) then
      option%io_buffer = '.not.associated(mfd_aux)'
      call printErrMsg(option)  
     stop
  end if


  do icell = 1, region%num_cells

    cell_id_local = region%cell_ids(icell)

    aux_var => mfd_aux%aux_vars(cell_id_local)

    do iface = 1,aux_var%numfaces
      face_id_ghosted = aux_var%face_id_gh(iface)
      if (coupler%faces_set(icell) == face_id_ghosted) then
        e2n_local((cell_id_local-1)*stride + iface) = -coupler%flow_condition%itype(RICHARDS_PRESSURE_DOF)
      end if
    end do
  end do

!  do icell = 1, region%num_cells
!
!    cell_id_local = region%cell_ids(icell)
!
!    do iface = 1,aux_var%numfaces
!        write(*,*) e2n_local((cell_id_local-1)*stride + iface)
!    end do
!    write(*,*)
!  end do

  if (e2n_size > 0) call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)

#endif 

!  stop

end subroutine CouplerAssignBCtoCells



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
  
  PetscInt :: CouplerGetNumConnectionsInList
  type(coupler_type), pointer :: coupler
  
  CouplerGetNumConnectionsInList = 0
  coupler => list%first
  
  do
    if (.not.associated(coupler)) exit
    CouplerGetNumConnectionsInList = CouplerGetNumConnectionsInList + &
                                     coupler%connection_set%num_connections
    coupler => coupler%next
  enddo

end function CouplerGetNumConnectionsInList

! ************************************************************************** !
!
! CouplerGetNumBoundConnectionsInListMFD: Returns the number of boundary connections associated
!	                                    with all couplers in the list. Establish connections between
!   	                                 local face_id and bound_face_id.
!										 (Since boundary fluxes allocated only for active boundary faces
!											they have different indexing
! author: Daniil Svyatskiy
! date: 11/04/10
!
! ************************************************************************** !
function CouplerGetNumBoundConnectionsInListMFD(grid, list, option)

  use Grid_module
  use Option_module

  implicit none
  

  type(coupler_list_type) :: list
  type(grid_type) :: grid
  type(option_type) :: option
  PetscInt :: CouplerGetNumBoundConnectionsInListMFD
 
  type(coupler_type), pointer :: coupler

  PetscInt :: numfaces_set, iconn, local_face_id, i
  
  iconn = 0

  coupler => list%first

  allocate(grid%fL2B(grid%nlmax_faces))

  grid%fL2B = 0

  do
    if (.not.associated(coupler)) exit

    numfaces_set = coupler%numfaces_set
    do i = 1, numfaces_set
       iconn = iconn + 1
       local_face_id = grid%fG2L(coupler%faces_set(i))
       if (local_face_id<=0) then
          write(*,*) "local_face_id for boundary face is <=0"
!          call printMsg(option)
       end if
       grid%fL2B(local_face_id) = iconn
    enddo
    coupler => coupler%next
  enddo

  CouplerGetNumBoundConnectionsInListMFD = iconn

end function CouplerGetNumBoundConnectionsInListMFD 

! ************************************************************************** !
!
! CouplerGetPtrFromList: Returns a pointer to the coupler matching 
!                        coupler_name
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
function CouplerGetPtrFromList(coupler_name,coupler_list)

  use String_module

  implicit none
  
  type(coupler_type), pointer :: CouplerGetPtrFromList
  character(len=MAXWORDLENGTH) :: coupler_name
  PetscInt :: length
  type(coupler_list_type) :: coupler_list

  type(coupler_type), pointer :: coupler
    
  nullify(CouplerGetPtrFromList)

  coupler => coupler_list%first
  do 
    if (.not.associated(coupler)) exit
    length = len_trim(coupler_name)
    if (length == len_trim(coupler%name) .and. &
        StringCompare(coupler%name,coupler_name,length)) then
      CouplerGetPtrFromList => coupler
      return
    endif
    coupler => coupler%next
  enddo
  
end function CouplerGetPtrFromList

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
  
  if (.not.associated(coupler_list)) return
  
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
  
  nullify(coupler%flow_condition)     ! since these are simply pointers to 
  nullify(coupler%tran_condition)     ! since these are simply pointers to 
  nullify(coupler%region)        ! conditoins in list, nullify

  if (associated(coupler%flow_aux_int_var)) &
    deallocate(coupler%flow_aux_int_var)
  nullify(coupler%flow_aux_int_var)
  if (associated(coupler%flow_aux_real_var)) &
    deallocate(coupler%flow_aux_real_var)
  nullify(coupler%flow_aux_real_var)

  call ConnectionDestroy(coupler%connection_set)
  nullify(coupler%connection_set)

  if (associated(coupler%faces_set)) &
   deallocate(coupler%faces_set)
  nullify(coupler%faces_set)
  
  deallocate(coupler)
  nullify(coupler)

end subroutine CouplerDestroy

end module Coupler_module
