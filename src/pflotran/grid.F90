module Grid_module

  use Structured_Grid_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Polyhedra_module
  use Connection_module
  use MFD_Aux_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
 
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type, public :: grid_type 
  
    character(len=MAXWORDLENGTH) :: ctype
    PetscInt :: itype  ! type of grid (e.g. structured_grid, implicit_unstructured_grid, etc.)
    
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt :: global_offset ! Offset of first cell on process in petsc ordering
    PetscInt :: nlmax_faces  ! Total number of non-ghosted faces in local domain.
    PetscInt :: ngmax_faces  ! Number of ghosted & non-ghosted faces in local domain.
    PetscInt :: nmax_faces  ! Number of ghosted & non-ghosted faces in local domain.
    PetscInt :: global_cell_offset, global_faces_offset  ! offsets for LP formulation
   
    ! Below, we define several arrays used for mapping between different 
    ! types of array indices.  Our terminology is as follows:
    !
    ! 'Local' indices are used to access arrays containing values that are 
    ! entirely local to the MPI process -- these arrays contain no "ghost" 
    ! entries used to hold copies of values that are owned by neighboring 
    ! processes.
    !
    ! 'Ghosted local' (or simply 'ghost') indices are used to access arrays 
    ! that contain additional entries that hold copies of values that are 
    ! owned by neighboring processes.  (These entries are filled in by 
    ! DMGlobalToLocalBegin/End() in the structured grid case.)
    !
    ! Entries of a vector created with DMCreateGlobalVector() should be 
    ! indexed using 'local' indices.  The array returned from a call to 
    ! VecGetArrayF90() on such a vector consists of local entries only and 
    ! NO ghost points.
    !
    ! Entries of a vector created with DMCreateLocalVector() should be 
    ! indexed using 'ghosted local' indices.  The array returned from a call 
    ! to VecGetArrayF90() on such a vector contains the truly3 local entries 
    ! as well as ghost points.
    !
    ! The index mapping arrays are the following:
    ! nL2G :  not collective, local processor: local  =>  ghosted local  
    ! nG2L :  not collective, local processor:  ghosted local => local  
    ! nG2A :  not collective, ghosted local => natural

    PetscInt, pointer :: nL2G(:), nG2L(:)
    PetscInt, pointer :: nG2A(:), nG2P(:), nG2LP(:)

    PetscInt, pointer :: fL2G(:), fG2L(:), fG2P(:), fL2P(:)
    PetscInt, pointer :: fL2B(:)
    Vec :: e2f             ! global vector to establish connection between global face_id and cell_id
    Vec :: e2n, e2n_LP     ! global cell connectivity vector

    ! For mapping faces between unstrutured grid and mimetic discretization
    PetscInt, pointer :: fU2M(:,:),fM2U(:)
    PetscInt :: discretization_itype

    PetscReal, pointer :: x(:), y(:), z(:) ! coordinates of ghosted grid cells

    PetscReal :: x_min_global, y_min_global, z_min_global
    PetscReal :: x_max_global, y_max_global, z_max_global
    PetscReal :: x_min_local, y_min_local, z_min_local
    PetscReal :: x_max_local, y_max_local, z_max_local

    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash_bins

    type(structured_grid_type), pointer :: structured_grid
    type(unstructured_grid_type), pointer :: unstructured_grid
    
    type(connection_set_list_type), pointer :: internal_connection_set_list
    type(connection_set_list_type), pointer :: boundary_connection_set_list
    type(face_type), pointer :: faces(:)
    type(mfd_type), pointer :: MFD

    ! For "Least Square Method" to compute flux
    ! This vector has information regarding how far away a ghost cell is from
    ! a local cell.
    ! 1) Zero-value represent a local cell.
    ! 2) If SNES stencil_width = 2, the maximum value of this vector can be 2.
    PetscInt,pointer :: ghosted_level(:)

    ! Displacement: i-th row entry = (x_i-x_0, y_i-y_0, z_i-z_0)
    ! Minv = inverse(disp^T * disp).
    ! These quantities depends on grid.
    !
    ! TODO[GB]: Create one sparse matrix for the local part, rather than
    !          individual matrix for each control volume.
    Mat,pointer      :: dispT(:)
    Mat,pointer      :: Minv(:)
    PetscReal, pointer :: jacfac(:,:,:)
    
    ! PETSC_TRUE -- if a cell is boundary cell,
    ! PETSC_FALSE -- if a cell is interior cell
    PetscBool, pointer :: bnd_cell(:)

  end type grid_type
  
  type, public :: face_type
    type(connection_set_type), pointer :: conn_set_ptr
    PetscInt :: id
  end type face_type

  public :: GridCreate, &
            GridDestroy, &
            GridComputeInternalConnect, &
            GridMapIndices, &
            GridComputeSpacing, &
            GridComputeCoordinates, &
            GridComputeVolumes, &
            GridComputeAreas, &
            GridLocalizeRegions, &
            GridLocalizeRegionsFromCellIDsUGrid, &
            GridPopulateConnection, &
            GridPopulateFaces, &
            GridCopyIntegerArrayToVec, &
            GridCopyRealArrayToVec, &
            GridCopyVecToIntegerArray, &
            GridCopyVecToRealArray, &
            GridCreateNaturalToGhostedHash, &
            GridDestroyHashTable, &
            GridGetLocalGhostedIdFromHash, &
            GridIndexToCellID, &
            GridComputeCell2FaceConnectivity, &
            GridComputeGlobalCell2FaceConnectivity, &
            GridGetGhostedNeighbors, &
            GridGetGhostedNeighborsWithCorners, &
            GridComputeNeighbors, &
            GridComputeMinv, &
            GridSaveBoundaryCellInfo
  
!geh  public :: VecGetArrayReadF90, VecRestoreArrayReadF90
  
contains

#if 0
! ************************************************************************** !
!
! VecGetArrayReadF90: Redirects VecGetArrayReadF90 to VecGetArrayF90 for
!                     older petsc-dev that lacks stubs for the read version
! author: Glenn Hammond
! date: 03/04/13
!
! ************************************************************************** !
subroutine VecGetArrayReadF90(vec, f90ptr, ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(vec, f90ptr, ierr)

end subroutine VecGetArrayReadF90
      
! ************************************************************************** !
!
! VecRestoreArrayReadF90: Redirects VecRestoreArrayReadF90 to VecGetArrayF90 
!                         older petsc-dev that lacks stubs for the read 
!                         for version
! author: Glenn Hammond
! date: 03/04/13
!
! ************************************************************************** !
subroutine VecRestoreArrayReadF90(vec, f90ptr, ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  call VecRestoreArrayF90(vec, f90ptr, ierr)
  
end subroutine VecRestoreArrayReadF90
#endif

! ************************************************************************** !
!
! GridCreate: Creates a structured or unstructured grid
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function GridCreate()

  implicit none
  
  type(grid_type), pointer :: GridCreate
  
  type(grid_type), pointer :: grid
  
  allocate(grid)
  grid%ctype = ''
  grid%itype = 0

  nullify(grid%structured_grid)
  nullify(grid%unstructured_grid)

  nullify(grid%internal_connection_set_list)
#ifdef DASVYAT
  nullify(grid%boundary_connection_set_list)
#endif

  nullify(grid%nL2G)
  nullify(grid%nG2L)
  nullify(grid%nG2A)
  nullify(grid%nG2P)
  nullify(grid%ghosted_level)
  nullify(grid%Minv)
  nullify(grid%dispT)
  nullify(grid%jacfac)
  nullify(grid%bnd_cell)

#ifdef DASVYAT
  nullify(grid%fL2G)
  nullify(grid%fG2L)
  nullify(grid%fG2P)
  nullify(grid%fL2P)
  nullify(grid%fL2B)
#endif

  nullify(grid%x)
  nullify(grid%y)
  nullify(grid%z)

  grid%x_min_global = 1.d20
  grid%x_max_global = -1.d20
  grid%y_min_global = 1.d20
  grid%y_max_global = -1.d20
  grid%z_min_global = 1.d20
  grid%z_max_global = -1.d20

  grid%x_min_local = 1.d20
  grid%x_max_local = -1.d20
  grid%y_min_local = 1.d20
  grid%y_max_local = -1.d20
  grid%z_min_local = 1.d20
  grid%z_max_local = -1.d20

  grid%nmax = 0
  grid%nlmax = 0 
  grid%ngmax = 0
  grid%global_offset = 0

#ifdef DASVYAT  
  nullify(grid%faces)
  nullify(grid%MFD)
  grid%e2f = 0
  grid%e2n = 0
  grid%e2n_LP = 0
#endif

  nullify(grid%fU2M)
  nullify(grid%fM2U)
  nullify(grid%hash)
  grid%num_hash_bins = 1000
  grid%discretization_itype = 0

  GridCreate => grid

end function GridCreate

! ************************************************************************** !
!
! GridComputeInternalConnect: computes internal connectivity of a grid
! author: Glenn Hammond
! date: 10/17/07
!
! sp modified December 2010 
! ************************************************************************** !
subroutine GridComputeInternalConnect(grid,option,ugdm)

  use Connection_module
  use Option_module
  use Unstructured_Explicit_module
  use Unstructured_Polyhedra_module
    
  implicit none
  
  PetscInt ierr

  type(grid_type) :: grid
  type(option_type) :: option
  type(ugdm_type), optional :: ugdm
  
  type(connection_set_type), pointer :: connection_set, connection_bound_set
  type(connection_set_type), pointer :: connection_set_2
  nullify(connection_set); nullify(connection_bound_set)
  
  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      connection_set => &
        StructGridComputeInternConnect( grid%structured_grid, grid%x, grid%y, &
                                    grid%z, option)
    case(IMPLICIT_UNSTRUCTURED_GRID) 
      connection_set => &
        UGridComputeInternConnect(grid%unstructured_grid,grid%x,grid%y, &
                                  grid%z,option)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      connection_set => &
        ExplicitUGridSetInternConnect(grid%unstructured_grid%explicit_grid, &
                                        option)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      connection_set => &
        PolyhedraUGridComputeInternConnect(grid%unstructured_grid, &
                                           grid%x, grid%y, grid%z, &
                                           option)
  end select
  
  allocate(grid%internal_connection_set_list)
  call ConnectionInitList(grid%internal_connection_set_list)
  call ConnectionAddToList(connection_set,grid%internal_connection_set_list)


  select case(grid%itype)
    case(STRUCTURED_GRID_MIMETIC)
#ifdef DASVYAT
      connection_bound_set => &
        StructGridComputeBoundConnect(grid%structured_grid, grid%x, grid%y, grid%z, option)
        grid%nlmax_faces = grid%structured_grid%nlmax_faces;
        grid%ngmax_faces = grid%structured_grid%ngmax_faces;
#endif
    case(IMPLICIT_UNSTRUCTURED_GRID) 
!      connection_bound_set => &
!        UGridComputeBoundConnect(grid%unstructured_grid,option)
  end select

#ifdef DASVYAT  
  if (associated(connection_bound_set)) then
    allocate(grid%boundary_connection_set_list)
    call ConnectionInitList(grid%boundary_connection_set_list)
    call ConnectionAddToList(connection_bound_set,grid%boundary_connection_set_list)
  endif
  if ((grid%itype==STRUCTURED_GRID_MIMETIC)) then
    allocate(grid%fL2G(grid%nlmax_faces))
    allocate(grid%fG2L(grid%ngmax_faces))
  
    grid%fL2G = 0
    grid%fG2L = 0
    call GridPopulateFaces(grid, option)
  endif
#endif

#ifdef MFD_UGRID
  if(grid%itype==IMPLICIT_UNSTRUCTURED_GRID) then

    call GridPopulateFaces(grid, option)

    ! Note: Allocation of memory has to happen after call to GridPopulateFaces()
    allocate(grid%fL2G(grid%nlmax_faces))
    allocate(grid%fG2L(grid%ngmax_faces))
    grid%fL2G = 0
    grid%fG2L = 0
  endif
#endif

end subroutine GridComputeInternalConnect

! ************************************************************************** !
!
! GridPopulateConnection: computes connectivity coupler to a grid
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine GridPopulateConnection(grid,connection,iface,iconn,cell_id_local, &
                                  option)

  use Connection_module
  use Structured_Grid_module
  use Option_module
  
  implicit none
 
  type(grid_type) :: grid
  type(connection_set_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: cell_id_local
  type(option_type) :: option
  
  PetscInt :: cell_id_ghosted
  
  cell_id_ghosted = grid%nL2G(cell_id_local)
  ! Use ghosted index to access dx, dy, dz because we have
  ! already done a global-to-local scatter for computing the
  ! interior node connections.
  
  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call StructGridPopulateConnection(grid%x,grid%structured_grid,connection, &
                                        iface,iconn,cell_id_ghosted,option)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridPopulateConnection(grid%unstructured_grid,connection,iface,&
                                   iconn,cell_id_ghosted,option)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      call PolyhedraUGridPopulateConnection(grid%unstructured_grid,connection,iface, &
                                            iconn,cell_id_ghosted,option)
  end select

end subroutine GridPopulateConnection


! ************************************************************************** !
!
! GridPopulateFaces: allocate and populate array of faces
! author: Daniil Svyatskiy
! date: 02/04/10
!
! ************************************************************************** !
subroutine GridPopulateFaces(grid, option)

   use Connection_module
   use Option_module
    
   type(grid_type) :: grid
   type(option_type) :: option

#ifdef DASVYAT
   PetscInt :: total_faces, face_id, iconn
   type(connection_set_type), pointer :: cur_connection_set
   type(connection_set_list_type), pointer :: connection_set_list
   type(face_type), pointer :: faces(:)
   character(len=MAXWORDLENGTH) :: filename

#ifdef MFD_UGRID
  if (grid%itype==IMPLICIT_UNSTRUCTURED_GRID) then
    call GridPopulateFacesForUGrid(grid,option)
    return
  endif
#endif

   total_faces = grid%ngmax_faces
   allocate(faces(total_faces))
   face_id = 0

   connection_set_list => grid%internal_connection_set_list
   cur_connection_set => connection_set_list%first
   do 
     if (.not.associated(cur_connection_set)) exit
     do iconn = 1, cur_connection_set%num_connections 
       face_id = face_id + 1
       faces(face_id)%conn_set_ptr => cur_connection_set
       faces(face_id)%id = iconn
     enddo
     cur_connection_set => cur_connection_set%next
   enddo

   connection_set_list => grid%boundary_connection_set_list
   cur_connection_set => connection_set_list%first
   do 
     if (.not.associated(cur_connection_set)) exit
     do iconn = 1, cur_connection_set%num_connections
       face_id = face_id + 1
       faces(face_id)%conn_set_ptr => cur_connection_set
       faces(face_id)%id = iconn
     enddo
     cur_connection_set => cur_connection_set%next
   enddo

   grid%faces => faces

#endif
end subroutine GridPopulateFaces 

! ************************************************************************** !
! ************************************************************************** !
subroutine GridComputeCell2FaceConnectivity(grid, MFD_aux, option)

  use MFD_Aux_module
  use Option_module

  implicit none

  type(grid_type) :: grid
  type(mfd_type), pointer :: MFD_aux
  
  type(option_type) :: option

#ifdef DASVYAT
  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  PetscInt :: icount, icell, iface, local_id
  PetscInt :: local_id_dn, local_id_up, ghosted_id_dn, ghosted_id_up
  character(len=MAXWORDLENGTH) :: filename

  PetscInt, pointer :: numfaces(:)

#ifdef MFD_UGRID
  if (grid%itype==IMPLICIT_UNSTRUCTURED_GRID) then
    call GridComputeCell2FaceForUGrid(grid,MFD_aux,option)
    return
  endif
#endif

  MFD_aux => MFDAuxCreate()
  grid%MFD => MFD_aux
 
  call MFDAuxInit(MFD_aux, grid%nlmax, option)
  allocate(numfaces(grid%nlmax))

  numfaces = 6

  do icell = 1, grid%nlmax
    aux_var => MFD_aux%aux_vars(icell)
    call MFDAuxVarInit(aux_var, numfaces(icell), option)
  enddo

  local_id = 1
  do icount = 1, grid%ngmax_faces
    conn => grid%faces(icount)%conn_set_ptr
    iface = grid%faces(icount)%id

    if (conn%itype==BOUNDARY_CONNECTION_TYPE) then
      ghosted_id_dn = conn%id_dn(iface)
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      if (local_id_dn>0) then
        aux_var => MFD_aux%aux_vars(local_id_dn)
        call MFDAuxAddFace(aux_var,option, icount)
        grid%fG2L(icount)=local_id
        grid%fL2G(local_id) = icount
        local_id = local_id + 1
      endif
        
    else if (conn%itype==INTERNAL_CONNECTION_TYPE) then 
      ghosted_id_up = conn%id_up(iface)
      ghosted_id_dn = conn%id_dn(iface)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      if (conn%local(iface) == 1) then
        grid%fG2L(icount)=local_id
        grid%fL2G(local_id) = icount
        local_id = local_id + 1
      endif

      if (local_id_dn>0) then
        aux_var => MFD_aux%aux_vars(local_id_dn)
        call MFDAuxAddFace(aux_var,option, icount)
      endif
      if (local_id_up>0) then
        aux_var => MFD_aux%aux_vars(local_id_up)
        call MFDAuxAddFace(aux_var,option, icount)
      endif
    endif
  enddo

  if (associated(numfaces)) deallocate(numfaces)

#endif

end subroutine GridComputeCell2FaceConnectivity


! ************************************************************************** !
! ************************************************************************** !
subroutine GridComputeGlobalCell2FaceConnectivity( grid, MFD_aux, sgdm, DOF, option)

  use MFD_Aux_module
  use Option_module
 
  implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"
#include "finclude/petscsys.h"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"

  type(grid_type) :: grid
  type(mfd_type), pointer :: MFD_aux
  DM :: sgdm
  PetscInt :: DOF
  type(option_type) :: option

#ifdef DASVYAT

  Vec :: ghosted_e2f
  Vec :: ghosted_e2n

  PetscErrorCode :: ierr
  PetscInt :: ndof
  PetscInt :: global_offset, stride
  PetscInt :: local_cell_id, ghost_cell_id, global_cell_id
  PetscInt :: local_face_id, ghost_face_id, global_face_id
  PetscInt :: face_id_gh, global_neigh_id
  PetscInt :: iface, icell, icount, jcount
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscInt :: num_ghosted_upd
  PetscInt :: istart, iend, temp_int
  IS :: is_a2g_bl
  IS :: is_local_bl
  PetscViewer :: viewer

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn

  PetscInt, allocatable :: ghosted_ids(:)
  PetscInt, allocatable :: strided_indices_local(:)
  PetscInt, allocatable :: strided_indices_ghosted(:)
  PetscInt, allocatable :: int_array(:)
  PetscInt, pointer :: int_ptr(:)

  PetscScalar, pointer :: e2f_local_values(:)
  PetscScalar, pointer :: e2n_local_values(:)
!  PetscScalar, pointer :: lp_cell_ids(:), lp_cell_ids_loc(:)

  PetscScalar, pointer :: e2f_ghosted_values(:)
  PetscScalar, pointer :: e2n_ghosted_values(:)
   
  PetscScalar, pointer :: vec_ptr_e2f(:)
  PetscScalar, pointer :: vec_ptr_e2n(:)
    
  PetscScalar, pointer :: vec_ptr_e2f_gh(:)
  PetscScalar, pointer :: vec_ptr_e2n_gh(:)

  VecScatter :: VC_global2ghosted

#ifdef MFD_UGRID
  if (grid%itype==IMPLICIT_UNSTRUCTURED_GRID) then
    call GridSetGlobalCell2FaceForUGrid(grid,MFD_aux,DOF,option)
    return
  endif
#endif

  select case(DOF)
    case(ONEDOF)
      ndof = 1
    case(NFLOWDOF)
      ndof = option%nflowdof
    case(NTRANDOF)
      ndof = option%ntrandof
  end select

  MFD_aux%ndof = ndof

  allocate(grid%fL2P(grid%nlmax_faces))
  allocate(grid%fG2P(grid%ngmax_faces))

  global_offset = 0
  do iface = 1,grid%nlmax_faces
    grid%fL2P(iface)=grid%global_faces_offset + grid%global_cell_offset + iface - 1
  enddo

  global_offset = grid%global_faces_offset + grid%global_cell_offset

  stride = 6                ! Only for hexagons

  call VecCreate(option%mycomm, grid%e2f, ierr)
  call VecSetSizes(grid%e2f, grid%nlmax*stride, PETSC_DECIDE, ierr)
  call VecSetFromOptions(grid%e2f, ierr)

  call VecDuplicate(grid%e2f, grid%e2n, ierr)

  call VecGetArrayF90(grid%e2f, e2f_local_values, ierr)
  call VecGetArrayF90(grid%e2n, e2n_local_values, ierr)

  e2f_local_values = 0
  e2n_local_values = 0

  do icell = 1, grid%nlmax
    aux_var => MFD_aux%aux_vars(icell)
    do icount = 1, aux_var%numfaces
      ghost_face_id = aux_var%face_id_gh(icount)
      local_face_id = grid%fG2L(ghost_face_id)
      conn => grid%faces(ghost_face_id)%conn_set_ptr
      iface = grid%faces(ghost_face_id)%id

      if (conn%itype==INTERNAL_CONNECTION_TYPE) then

        if (local_face_id > 0) then

          e2f_local_values((icell-1)*stride + icount) = grid%fL2P(local_face_id) + 1
          ghosted_id_up = conn%id_up(iface)
          ghosted_id_dn = conn%id_dn(iface)

          local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
          local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
          if (local_id_up==icell) then
            e2n_local_values((icell-1)*stride + icount) = grid%nG2P(ghosted_id_dn) + 1
          else if (local_id_dn==icell) then
            e2n_local_values((icell-1)*stride + icount) = grid%nG2P(ghosted_id_up) + 1
          endif

        else if (local_face_id == 0) then

          ghosted_id_up = conn%id_up(iface)
          ghosted_id_dn = conn%id_dn(iface)

          local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
          local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
          if (local_id_up==icell) then
            e2n_local_values((icell-1)*stride + icount) = grid%nG2P(ghosted_id_dn) + 1
          else if (local_id_dn==icell) then
            e2n_local_values((icell-1)*stride + icount) = grid%nG2P(ghosted_id_up) + 1
          endif
        endif

      else if (conn%itype==BOUNDARY_CONNECTION_TYPE) then
        e2n_local_values((icell-1)*stride + icount) = 0
        if (local_face_id > 0) e2f_local_values((icell-1)*stride + icount) = grid%fL2P(local_face_id) + 1
      endif
    enddo
  enddo

  call VecRestoreArrayF90(grid%e2f, e2f_local_values, ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local_values, ierr)

  allocate(ghosted_ids(grid%ngmax_faces - grid%nlmax_faces))
  ghosted_ids = 0

  num_ghosted_upd = 0
  do icount = 1, grid%ngmax_faces
    conn => grid%faces(icount)%conn_set_ptr
    iface = grid%faces(icount)%id
    if (conn%itype==INTERNAL_CONNECTION_TYPE) then
      if (conn%local(iface) == 0) then
        num_ghosted_upd = num_ghosted_upd + 1
        ghosted_id_up = conn%id_up(iface)
        ghosted_id_dn = conn%id_dn(iface)
        if (grid%nG2L(ghosted_id_up)==0) then    ! up_cell is ghosted
          ghosted_ids(num_ghosted_upd) = grid%nG2P(ghosted_id_up) + 1          ! +1 since global indexes strats from 0
        else if (grid%nG2L(ghosted_id_dn)==0) then    ! down_cell is ghosted
          ghosted_ids(num_ghosted_upd) = grid%nG2P(ghosted_id_dn) + 1
        endif
      endif
    endif
  enddo

  call VecCreateSeq(PETSC_COMM_SELF, num_ghosted_upd*stride, ghosted_e2f, ierr)
  call VecCreateSeq(PETSC_COMM_SELF, num_ghosted_upd*stride, ghosted_e2n, ierr)
         
  allocate(strided_indices_local(num_ghosted_upd))
  allocate(strided_indices_ghosted(num_ghosted_upd))

  do icount = 1, num_ghosted_upd
    strided_indices_local(icount) = (icount -1)
    strided_indices_ghosted(icount) = (ghosted_ids(icount)-1)
  enddo

  call ISCreateBlock(option%mycomm, stride, num_ghosted_upd, strided_indices_local, &
                    PETSC_COPY_VALUES, is_local_bl, ierr)
  call ISCreateBlock(option%mycomm, stride, num_ghosted_upd, strided_indices_ghosted,&
                    PETSC_COPY_VALUES, is_a2g_bl, ierr)

  call VecScatterCreate(grid%e2f, is_a2g_bl, ghosted_e2f, is_local_bl, VC_global2ghosted, ierr)

  call ISDestroy(is_local_bl, ierr)
  call ISDestroy(is_a2g_bl, ierr)

  call VecScatterBegin(VC_global2ghosted, grid%e2f, ghosted_e2f, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(VC_global2ghosted, grid%e2f, ghosted_e2f, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecScatterBegin(VC_global2ghosted, grid%e2n, ghosted_e2n, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(VC_global2ghosted, grid%e2n, ghosted_e2n, &
                    INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(ghosted_e2n, vec_ptr_e2n_gh, ierr)
  call VecGetArrayF90(ghosted_e2f, vec_ptr_e2f_gh, ierr)

  call VecGetArrayF90(grid%e2n, vec_ptr_e2n, ierr)
  call VecGetArrayF90(grid%e2f, vec_ptr_e2f, ierr)
 
  do icell = 1, grid%nlmax
    aux_var => MFD_aux%aux_vars(icell)
    do icount = 1, aux_var%numfaces
      ghost_face_id = aux_var%face_id_gh(icount)
      local_face_id = grid%fG2L(ghost_face_id)

      if (local_face_id == 0) then
        conn => grid%faces(ghost_face_id)%conn_set_ptr
        iface = grid%faces(ghost_face_id)%id

        ghosted_id_up = conn%id_up(iface)
        ghosted_id_dn = conn%id_dn(iface)

        local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
        local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
        if (icell==local_id_up) then
          global_neigh_id = grid%nG2P(ghosted_id_dn) + 1
        else if (icell==local_id_dn) then
          global_neigh_id = grid%nG2P(ghosted_id_up) + 1
        endif
        do jcount = 1, num_ghosted_upd
          if (ghosted_ids(jcount)==global_neigh_id) exit
        enddo
        do iface=1, stride             ! assumption that cell has 6 faces
          if (vec_ptr_e2n_gh((jcount-1)*stride + iface) == grid%nG2P(grid%nL2G(icell)) + 1) then
            vec_ptr_e2f((icell -1)*stride +icount) = vec_ptr_e2f_gh((jcount-1)*stride + iface)
            grid%fG2P(ghost_face_id) = int(vec_ptr_e2f_gh((jcount-1)*stride + iface)) - 1
            exit
          endif
        end do
      else
        grid%fG2P(ghost_face_id) = grid%fL2P(local_face_id)
      endif
    enddo
  enddo

  call VecRestoreArrayF90(grid%e2n, vec_ptr_e2n, ierr)
  call VecRestoreArrayF90(grid%e2f, vec_ptr_e2f, ierr)
  call VecRestoreArrayF90(ghosted_e2n, vec_ptr_e2n_gh, ierr)
  call VecRestoreArrayF90(ghosted_e2f, vec_ptr_e2f_gh, ierr)

  call VecDestroy(ghosted_e2n, ierr)
  call VecDestroy(ghosted_e2f, ierr)
  call VecScatterDestroy(VC_global2ghosted, ierr)

  call CreateMFDStruct4LP(grid, MFD_aux, ndof, option)

  deallocate(ghosted_ids)
  deallocate(strided_indices_local)
  deallocate(strided_indices_ghosted)
#endif

end subroutine GridComputeGlobalCell2FaceConnectivity

! ************************************************************************** !
! ************************************************************************** !
subroutine CreateMFDStruct4Faces(grid, MFD_aux, ndof, option)

    
  use MFD_Aux_module
  use Option_module
 
  implicit none

  type(grid_type) :: grid
  type(mfd_type), pointer :: MFD_aux
  PetscInt :: ndof
  type(option_type) :: option

  PetscInt, allocatable :: int_array(:)
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: iface, icell, icount, jcount
  PetscInt :: istart, iend, temp_int
  PetscErrorCode :: ierr
  PetscInt :: global_offset, stride

  Vec :: global_vec
  Vec :: local_vec

  PetscViewer :: viewer

  ! create global vec
  call VecCreateMPI(option%mycomm, grid%nlmax_faces*ndof, &
                    PETSC_DETERMINE, global_vec,ierr)
  call VecSetBlockSize(global_vec,ndof,ierr)  
  
  ! create local vec
  call VecCreateSeq(PETSC_COMM_SELF, grid%ngmax_faces*ndof, &
                    local_vec,ierr)
  call VecSetBlockSize(local_vec,ndof,ierr)

  ! IS for global numbering of local, non-ghosted cells
  call VecGetOwnershipRange(global_vec,istart,iend,ierr)


  allocate(int_array(grid%nlmax_faces))
  do iface = 1, grid%nlmax_faces
    int_array(iface) = (iface-1) + istart/ndof
  enddo

  call ISCreateBlock(option%mycomm, ndof, grid%nlmax_faces, &
                     int_array, PETSC_COPY_VALUES, MFD_aux%is_local_petsc_faces, ierr)
  deallocate(int_array)


  allocate(int_array(grid%ngmax_faces - grid%nlmax_faces))
  do iface =  1, grid%ngmax_faces - grid%nlmax_faces
     int_array(iface) = (iface + grid%nlmax_faces  - 1)
  end do

  call ISCreateBlock(option%mycomm,ndof,grid%ngmax_faces - grid%nlmax_faces, &
                     int_array, PETSC_COPY_VALUES, MFD_aux%is_ghosts_local_faces,ierr)
  deallocate(int_array)

  allocate(int_array(grid%ngmax_faces - grid%nlmax_faces))
  do iface =  1, grid%ngmax_faces - grid%nlmax_faces
     int_array(iface) = (grid%fG2P(grid%nlmax_faces + iface))
  end do
  
  call ISCreateBlock(option%mycomm,ndof,grid%ngmax_faces - grid%nlmax_faces, &
                     int_array, PETSC_COPY_VALUES, MFD_aux%is_ghosts_petsc_faces,ierr)

  deallocate(int_array)

  allocate(int_array(grid%nlmax_faces))
  do iface = 1, grid%nlmax_faces
    int_array(iface)=(iface - 1)
  end do

  call ISCreateBlock(option%mycomm,ndof,grid%nlmax_faces, &
      int_array, PETSC_COPY_VALUES, MFD_aux%is_local_local_faces, ierr)
  deallocate(int_array)

  allocate(int_array(grid%ngmax_faces))
  do iface = 1, grid%ngmax_faces
    int_array(iface)=(iface - 1)
  end do

  call ISCreateBlock(option%mycomm,ndof,grid%ngmax_faces, &
      int_array, PETSC_COPY_VALUES, MFD_aux%is_ghosted_local_faces, ierr)
  deallocate(int_array)

   allocate(int_array(grid%ngmax_faces))
   do iface = 1, grid%ngmax_faces
      int_array(iface) = (grid%fG2P(iface))
   end do

  call ISCreateBlock(option%mycomm,ndof,grid%ngmax_faces, &
       int_array, PETSC_COPY_VALUES, MFD_aux%is_ghosted_petsc_faces, ierr)
  deallocate(int_array)

  call ISLocalToGlobalMappingCreateIS(MFD_aux%is_ghosted_petsc_faces, &
                                    MFD_aux%mapping_ltog_faces,ierr)

  call ISLocalToGlobalMappingBlock(MFD_aux%mapping_ltog_faces,ndof, &
                                  MFD_aux%mapping_ltogb_faces,ierr)

  call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_petsc.out',viewer,ierr)
  call ISView(MFD_aux%is_ghosted_petsc_faces,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call VecScatterCreate(local_vec,MFD_aux%is_local_local_faces,global_vec, &
                        MFD_aux%is_local_petsc_faces,MFD_aux%scatter_ltog_faces,ierr)
  call VecScatterCreate(global_vec,MFD_aux%is_ghosted_petsc_faces,local_vec, &
                        MFD_aux%is_ghosted_local_faces, MFD_aux%scatter_gtol_faces, ierr)

  ! Create local to local scatter.  Essentially remap the global to local as
  ! PETSc does in daltol.c
  call VecScatterCopy(MFD_aux%scatter_gtol_faces, MFD_aux%scatter_ltol_faces, ierr)
  call ISGetIndicesF90(MFD_aux%is_local_local_faces,int_ptr,ierr)
  call VecScatterRemap(MFD_aux%scatter_ltol_faces,int_ptr,PETSC_NULL_INTEGER,ierr)
  call ISRestoreIndicesF90(MFD_aux%is_local_local_faces,int_ptr,ierr)

  call VecDestroy(local_vec, ierr)
  call VecDestroy(global_vec, ierr)

end subroutine CreateMFDStruct4Faces

! ************************************************************************** !
! ************************************************************************** !
subroutine CreateMFDStruct4LP(grid, MFD_aux, ndof, option)

  use MFD_Aux_module
  use Option_module
 
  implicit none

  type(grid_type) :: grid
  type(mfd_type), pointer :: MFD_aux
  type(option_type) :: option
  PetscInt :: ndof
  PetscViewer :: viewer

  PetscInt, allocatable :: int_array(:)
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: icount, jcount, id
  PetscInt :: istart, iend, temp_int
  PetscErrorCode :: ierr
  PetscInt :: global_offset, stride
  PetscInt :: NG, NL

  Vec :: global_vec_LP
  Vec :: local_vec_LP

  NG = grid%ngmax_faces + grid%ngmax
  NL = grid%nlmax_faces + grid%nlmax

  ! create global vec_LP
  call VecCreateMPI(option%mycomm, NL*ndof, &
                    PETSC_DETERMINE, global_vec_LP,ierr)
  call VecSetBlockSize(global_vec_LP,ndof,ierr)
  ! create local vec
  call VecCreateSeq(PETSC_COMM_SELF, NG*ndof, &
                    local_vec_LP,ierr)
  call VecSetBlockSize(local_vec_LP,ndof,ierr)

  ! IS for global numbering of local, non-ghosted cells
  call VecGetOwnershipRange(global_vec_LP,istart,iend,ierr)

  allocate(int_array(NL))
  do id = 1, NL
    int_array(id) = (id-1) + istart/ndof
  enddo

  call ISCreateBlock(option%mycomm, ndof, NL, &
                     int_array, PETSC_COPY_VALUES, MFD_aux%is_local_petsc_LP, ierr)
  deallocate(int_array)


  allocate(int_array(NG - NL))
  do id =  1, grid%ngmax_faces - grid%nlmax_faces
     int_array(id) = (id + grid%nlmax_faces  - 1)
  end do

  do id =  1, grid%ngmax - grid%nlmax
     int_array(id + grid%ngmax_faces - grid%nlmax_faces) = (grid%ngmax_faces + grid%nlmax + id - 1)
  end do

  call ISCreateBlock(option%mycomm,ndof,NG - NL, &
                     int_array, PETSC_COPY_VALUES, MFD_aux%is_ghosts_local_LP,ierr)
  deallocate(int_array)

  allocate(int_array(NL))
  do id = 1, NL
    int_array(id)=(id - 1)
  end do

  call ISCreateBlock(option%mycomm,ndof, NL, &
       int_array, PETSC_COPY_VALUES, MFD_aux%is_local_local_LP, ierr)
  deallocate(int_array)

  allocate(int_array(NG))
  do id = 1, NG
    int_array(id)=(id - 1)
  end do

  call ISCreateBlock(option%mycomm,ndof, NG, &
       int_array, PETSC_COPY_VALUES, MFD_aux%is_ghosted_local_LP, ierr)
  deallocate(int_array)

  call VecScatterCreate(local_vec_LP,MFD_aux%is_local_local_LP,global_vec_LP, &
                        MFD_aux%is_local_petsc_LP,MFD_aux%scatter_ltog_LP,ierr)


  allocate(int_array(NG))
  do id = 1, grid%ngmax_faces
    int_array(id) = (grid%fG2P(id))
  end do

  do id = 1, grid%ngmax
    int_array(grid%ngmax_faces + id) =  grid%nG2LP(id)
  end do

  call ISCreateBlock(option%mycomm,ndof, NG, &
       int_array, PETSC_COPY_VALUES, MFD_aux%is_ghosted_petsc_LP, ierr)
  deallocate(int_array)

  call ISLocalToGlobalMappingCreateIS(MFD_aux%is_ghosted_petsc_LP, &
                                      MFD_aux%mapping_ltog_LP,ierr)

  call ISLocalToGlobalMappingBlock(MFD_aux%mapping_ltog_LP, ndof, &
                                   MFD_aux%mapping_ltogb_LP, ierr)

  call VecScatterCreate(global_vec_LP,MFD_aux%is_ghosted_petsc_LP,local_vec_LP, &
                        MFD_aux%is_ghosted_local_LP, MFD_aux%scatter_gtol_LP, ierr)

 ! Create local to local scatter.  Essentially remap the global to local as
 ! PETSc does in daltol.c
  call VecScatterCopy(MFD_aux%scatter_gtol_LP, MFD_aux%scatter_ltol_LP, ierr)
  call ISGetIndicesF90(MFD_aux%is_local_local_LP,int_ptr,ierr)
  call VecScatterRemap(MFD_aux%scatter_ltol_LP,int_ptr,PETSC_NULL_INTEGER,ierr)
  call ISRestoreIndicesF90(MFD_aux%is_local_local_LP,int_ptr,ierr)

  call VecDestroy(local_vec_LP, ierr)
  call VecDestroy(global_vec_LP, ierr)

end subroutine CreateMFDStruct4LP

! ************************************************************************** !
!
! GridMapIndices: maps global, local and natural indices of cells 
!                 to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridMapIndices(grid, sgdm, stencil_type, lsm_flux_method, option)

use Option_module

  implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"


  
  type(grid_type) :: grid
  DM :: sgdm
  PetscInt :: stencil_type
  PetscBool :: lsm_flux_method
  type(option_type) :: option

  PetscInt :: ierr, icount
  PetscInt, allocatable :: int_tmp(:)
! PetscInt, pointer :: int_tmp(:)
  PetscInt :: n
  PetscOffset :: i_da
  
  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call StructGridMapIndices(grid%structured_grid,stencil_type, &
                                    lsm_flux_method, &
                                    grid%nG2L,grid%nL2G,grid%nG2A, &
                                    grid%ghosted_level,option)
#ifdef DASVYAT
      if ((grid%itype==STRUCTURED_GRID_MIMETIC)) then
        allocate(grid%nG2P(grid%ngmax))
        allocate(int_tmp(grid%ngmax))
!geh     call DMDAGetGlobalIndicesF90(sgdm, n, int_tmp, ierr)
        call DMDAGetGlobalIndices(sgdm,  grid%ngmax, int_tmp, i_da, ierr)
        do icount = 1, grid%ngmax
!geh         write(*,*) icount,  int_tmp(icount + i_da)
          grid%nG2P(icount) = int_tmp(icount + i_da)
!             write(*,*) icount,  int_tmp(icount)
!geh        grid%nG2P(icount) = int_tmp(icount)
        enddo
        deallocate(int_tmp)
      endif
#endif
    case(IMPLICIT_UNSTRUCTURED_GRID)
  end select
 
 
end subroutine GridMapIndices

! ************************************************************************** !
!
! GridComputeSpacing: Computes grid spacing (only for structured grid
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine GridComputeSpacing(grid,option)

  use Option_module

  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call StructGridComputeSpacing(grid%structured_grid,option)
    case(IMPLICIT_UNSTRUCTURED_GRID)
  end select
  
end subroutine GridComputeSpacing

! ************************************************************************** !
!
! GridComputeCoordinates: Computes x,y,z coordinates of grid cells
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridComputeCoordinates(grid,origin_global,option,ugdm)

  use Option_module
  use Unstructured_Explicit_module
  use Unstructured_Polyhedra_module
  
  implicit none

  type(grid_type) :: grid
  PetscReal :: origin_global(3)
  type(option_type) :: option
  type(ugdm_type), optional :: ugdm ! sp 
  PetscInt :: icell
  
  PetscErrorCode :: ierr  

  allocate(grid%x(grid%ngmax))
  grid%x = 0.d0
  allocate(grid%y(grid%ngmax))
  grid%y = 0.d0
  allocate(grid%z(grid%ngmax))
  grid%z = 0.d0
  
  select case(grid%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      call StructGridComputeCoord(grid%structured_grid,option, &
                                      origin_global, &
                                      grid%x,grid%y,grid%z, &
                                      grid%x_min_local,grid%x_max_local, &
                                      grid%y_min_local,grid%y_max_local, &
                                      grid%z_min_local,grid%z_max_local)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridComputeCoord(grid%unstructured_grid,option, &
                             grid%x,grid%y,grid%z, &
                             grid%x_min_local,grid%x_max_local, &
                             grid%y_min_local,grid%y_max_local, &
                             grid%z_min_local,grid%z_max_local)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      call ExplicitUGridSetCellCentroids(grid%unstructured_grid% &
                                         explicit_grid, &
                                         grid%x,grid%y,grid%z, &
                             grid%x_min_local,grid%x_max_local, &
                             grid%y_min_local,grid%y_max_local, &
                             grid%z_min_local,grid%z_max_local)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      call PolyhedraUGridSetCellCentroids(grid%unstructured_grid%polyhedra_grid, &
                                          grid%x,grid%y,grid%z, &
                                          grid%x_min_local,grid%x_max_local, &
                                          grid%y_min_local,grid%y_max_local, &
                                          grid%z_min_local,grid%z_max_local,option)

  end select

  if (associated(grid%structured_grid)) then
    ! compute global max/min from the local max/in
    call MPI_Allreduce(grid%x_min_local,grid%x_min_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
    call MPI_Allreduce(grid%y_min_local,grid%y_min_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
    call MPI_Allreduce(grid%z_min_local,grid%z_min_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
    call MPI_Allreduce(grid%x_max_local,grid%x_max_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    call MPI_Allreduce(grid%y_max_local,grid%y_max_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    call MPI_Allreduce(grid%z_max_local,grid%z_max_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
 endif

  if (associated(grid%unstructured_grid)) then
     ! compute global max/min from the local max/in
     call MPI_Allreduce(grid%x_min_local,grid%x_min_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
     call MPI_Allreduce(grid%y_min_local,grid%y_min_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
     call MPI_Allreduce(grid%z_min_local,grid%z_min_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
     call MPI_Allreduce(grid%x_max_local,grid%x_max_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
     call MPI_Allreduce(grid%y_max_local,grid%y_max_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
     call MPI_Allreduce(grid%z_max_local,grid%z_max_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
   !endif
 endif

end subroutine GridComputeCoordinates

! ************************************************************************** !
!
! GridComputeVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GridComputeVolumes(grid,volume,option)

  use Option_module
  use Unstructured_Explicit_module
  use Unstructured_Polyhedra_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(grid_type) :: grid
  type(option_type) :: option
  Vec :: volume
  
  select case(grid%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      call StructGridComputeVolumes(grid%x,grid%structured_grid,option, &
                                        grid%nL2G,volume)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridComputeVolumes(grid%unstructured_grid,option,volume)
      call UGridComputeQuality(grid%unstructured_grid,option)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      call ExplicitUGridComputeVolumes(grid%unstructured_grid, &
                                       option,volume)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      call PolyhedraUGridComputeVolumes(grid%unstructured_grid,option,volume)
  end select

end subroutine GridComputeVolumes

! ************************************************************************** !
!
! GridComputeAreas: Computes the areas for 2D-mesh
! author: Gautam Bisht
! date: 03/07/2012
!
! ************************************************************************** !
subroutine GridComputeAreas(grid,area,option)

  use Option_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(grid_type) :: grid
  type(option_type) :: option
  Vec :: area
  
  select case(grid%itype)
    !case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      !call StructGridComputeVolumes(grid%x,grid%structured_grid,option, &
      !                                  grid%nL2G,volume)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridComputeAreas(grid%unstructured_grid,option,area)
      call UGridComputeQuality(grid%unstructured_grid,option)
    case default
      option%io_buffer = 'ERROR: GridComputeAreas only implemented for Unstructured grid'
      call printErrMsg(option)
  end select

end subroutine GridComputeAreas

! ************************************************************************** !
!
! GridLocalizeRegions: Resticts regions to cells local to processor
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine GridLocalizeRegions(grid,region_list,option)

  use Option_module
  use Region_module

  implicit none
  
  type(region_list_type), pointer :: region_list
  type(grid_type), pointer :: grid
  type(option_type) :: option
  
  type(region_type), pointer :: region
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscReal :: del_x, del_y, del_z
  PetscInt :: iflag
  PetscBool :: same_point
  PetscErrorCode :: ierr
  
  iflag = 0
  region => region_list%first
  do
  
    if (.not.associated(region)) exit
    
    if (.not.(associated(region%cell_ids) .or. &
              associated(region%sideset) .or. &
              associated(region%polygonal_volume) .or. &
              associated(region%vertex_ids))) then
      ! i, j, k block
      if (region%i1 > 0 .and. region%i2 > 0 .and. &
          region%j1 > 0 .and. region%j2 > 0 .and. &
          region%k1 > 0 .and. region%k2 > 0) then

        call GridLocalizeRegionFromBlock(grid,region,option)

      else if (associated(region%coordinates)) then
        call GridLocalizeRegionFromCoordinates(grid,region,option)
      endif
    else if (associated(region%sideset)) then
      call UGridMapSideSet(grid%unstructured_grid, & 
                           region%sideset%face_vertices, & 
                           region%sideset%nfaces,region%name, & 
                           option,region%cell_ids,region%faces) 
      region%num_cells = size(region%cell_ids)
    else if (associated(region%polygonal_volume)) then
      call UGridMapBoundFacesInPolVol(grid%unstructured_grid, &
                                      region%polygonal_volume, &
                                      region%name,option, &
                                      region%cell_ids,region%faces)
      region%num_cells = size(region%cell_ids)
    else if (associated(region%cell_ids)) then
      select case(grid%itype) 
        case(IMPLICIT_UNSTRUCTURED_GRID)
          call GridLocalizeRegionsFromCellIDsUGrid(grid,region,option)

        case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
!sp following was commented out 
!sp remove? 
!geh: Not for now.  The below maps a list of natural ids to local.  This is now done
!     in hdf5.F90 for lists of cells from hdf5 files.  If we elect to support ascii
!     files, we will need this functionality here.
#if 0
          allocate(temp_int_array(region%num_cells))
          do count=1,region%num_cells
            i = mod(region%cell_ids(count),grid%structured_grid%nx) - &
                  grid%structured_grid%lxs
            j = mod((region%cell_ids(count)-1)/grid%structured_grid%nx, &
                    grid%structured_grid%ny)+1 - &
                  grid%structured_grid%lys
            k = ((region%cell_ids(count)-1)/grid%structured_grid%nxy)+1 - &
                  grid%structured_grid%lzs
            if (i > 0 .and. i <= grid%structured_grid%nlx .and. &
                j > 0 .and. j <= grid%structured_grid%nly .and. &
                k > 0 .and. k <= grid%structured_grid%nlz) then
              temp_int_array(local_count) = &
                  i + (j-1)*grid%structured_grid%nlx + &
                  (k-1)*grid%structured_grid%nlxy
              local_count = local_count + 1
            endif
          enddo
#endif
        case(EXPLICIT_UNSTRUCTURED_GRID)
          call GridLocalizeExplicitFaceset(grid%unstructured_grid,region, &
                                           option)
        case default
          option%io_buffer = 'GridLocalizeRegions: define region by list ' // &
            'of cells not implemented: ' // trim(region%name)
          call printErrMsg(option)
      end select
      !sp end 
    endif
    
    if (region%num_cells == 0 .and. associated(region%cell_ids)) then
      deallocate(region%cell_ids)
      nullify(region%cell_ids)
    endif
    if (region%num_cells == 0 .and. associated(region%faces)) then
      deallocate(region%faces)
      nullify(region%faces)
    endif
    region => region%next
    
  enddo

end subroutine GridLocalizeRegions

! ************************************************************************** !
!
! GridLocalizeRegionsFromCellIDsUGrid: Resticts regions to cells local 
!    to processor for unstrucutred mesh when the region is defined by a
!    list of cell ids (in natural order)
! author: Gautam Bisht
! date: 5/30/2011
!
! ************************************************************************** !
subroutine GridLocalizeRegionsFromCellIDsUGrid(grid, region, option)

  use Option_module
  use Region_module

  implicit none
  
#include "finclude/petsclog.h"
#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscmat.h"

  type(grid_type)                 :: grid
  type(region_type)               :: region
  type(option_type)               :: option

  ! local
  type(unstructured_grid_type),pointer    :: ugrid
  Vec                             :: vec_cell_ids,vec_cell_ids_loc
  Vec                             :: vec_face_ids,vec_face_ids_loc
  IS                              :: is_from, is_to
  Mat                             :: adj, adj_t,adj_d, adj_o
  VecScatter                      :: vec_scat
  PetscErrorCode                  :: ierr
  PetscViewer                     :: viewer
  PetscInt                        :: ii,jj,kk,count
  PetscInt                        :: istart, iend
  PetscInt                        :: ghosted_id,local_id,natural_id
  PetscInt,pointer                :: tmp_int_array(:),tmp_int_array2(:)
  PetscScalar,pointer             :: v_loc_p(:),v_loc2_p(:)
  PetscScalar,pointer             :: tmp_scl_array(:)


  PetscInt, pointer               :: ia_p(:), ja_p(:)
  PetscInt                        :: n,rstart,rend,icol(1)
  PetscInt                        :: index
  PetscInt                        :: vertex_id
  PetscOffset                     :: iia,jja,aaa,iicol
  PetscBool                       :: done,found
  PetscScalar                     :: aa(1)
  ! PetscScalar, pointer            :: aa(:)
  ! Would like to use the above, but I have to fix MatGetArrayF90() first. --RTM
  
  Mat                 :: mat_vert2cell, mat_vert2cell_diag, mat_vert2cell_offdiag
  Vec                 :: vec_vert2cell, vec_cell2facevert
  Vec                 :: vec_vert2cell_reg_subset, vec_cell2facevert_reg_subset
  PetscInt            :: vert_id_loc, vert_id_nat, counter1, counter2
  PetscInt,pointer    :: cell_count(:), cell_ids(:)
  PetscInt,pointer    :: cell_ids_for_face(:), face_ids_for_face(:)
  PetscScalar,pointer :: vert2cell_array(:)

  ugrid => grid%unstructured_grid
  
  if (associated(region%cell_ids)) then
    
    call VecCreateMPI(option%mycomm, ugrid%nlmax, PETSC_DECIDE, &
                      vec_cell_ids, ierr)
    call VecCreateMPI(option%mycomm, ugrid%nlmax, PETSC_DECIDE, &
                      vec_cell_ids_loc, ierr)
    
    call VecZeroEntries(vec_cell_ids, ierr)
    
    allocate(tmp_int_array(region%num_cells))
    allocate(tmp_scl_array(region%num_cells))

    count = 0
    do ii = 1, region%num_cells
      count = count + 1
      tmp_int_array(count) = region%cell_ids(ii)
      tmp_scl_array(count) = 1.d0
    enddo

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_cell_ids_bef.out', &
                              viewer, ierr)
    call VecView(vec_cell_ids, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    call VecSetValues(vec_cell_ids, region%num_cells, tmp_int_array, &
                      tmp_scl_array, ADD_VALUES, ierr)
    
    deallocate(tmp_int_array)
    deallocate(tmp_scl_array)

    call VecAssemblyBegin(vec_cell_ids, ierr)
    call VecAssemblyEnd(vec_cell_ids, ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_cell_ids_aft.out', &
                              viewer, ierr)
    call VecView(vec_cell_ids, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    allocate(tmp_int_array(ugrid%nlmax))
    count = 0
    do ghosted_id = 1, ugrid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (local_id < 1) cycle
      count = count + 1
      natural_id = grid%nG2A(ghosted_id)
      tmp_int_array(count) = natural_id
    enddo
    
    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, ugrid%nlmax, &
                        tmp_int_array, PETSC_COPY_VALUES, is_from, ierr)
    
    call VecGetOwnershipRange(vec_cell_ids_loc,istart,iend,ierr)
    do ii=1,ugrid%nlmax
      tmp_int_array(ii) = ii + istart
    enddo

    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, ugrid%nlmax, &
                        tmp_int_array, PETSC_COPY_VALUES, is_to, ierr)
    deallocate(tmp_int_array)
    
    call VecScatterCreate(vec_cell_ids,is_from,vec_cell_ids_loc,is_to, &
                          vec_scat, ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)
    
    call VecScatterBegin(vec_scat, vec_cell_ids, vec_cell_ids_loc, &
                          INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, vec_cell_ids, vec_cell_ids_loc, &
                        INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_cell_ids_loc.out', &
                              viewer, ierr)
    call VecView(vec_cell_ids_loc, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    call VecGetArrayF90(vec_cell_ids_loc, v_loc_p, ierr)
    count = 0
    do ii=1, ugrid%nlmax
      if (v_loc_p(ii) == 1) count = count + 1
    enddo
    
    region%num_cells = count
    if (count > 0) then
      allocate(tmp_int_array(count))
      count = 0
      do ii =1, ugrid%nlmax
        if (v_loc_p(ii) == 1) then
          count = count + 1
          tmp_int_array(count) = ii
        endif
      enddo

      deallocate(region%cell_ids)
      allocate(region%cell_ids(region%num_cells))
      region%cell_ids = tmp_int_array
      deallocate(tmp_int_array)
    endif
    
    call VecRestoreArrayF90(vec_cell_ids_loc,v_loc_p,ierr)
    
    call VecDestroy(vec_cell_ids,ierr)
    call VecDestroy(vec_cell_ids_loc,ierr)

  endif

end subroutine GridLocalizeRegionsFromCellIDsUGrid

! ************************************************************************** !
!
! GridLocalizeExplicitFaceset
! author: Glenn Hammond
! date: 10/10/12
!
! ************************************************************************** !
subroutine GridLocalizeExplicitFaceset(ugrid,region,option)

  use Region_module
  use Option_module

  implicit none
  
  type(unstructured_grid_type) :: ugrid
  type(region_type) :: region
  type(option_type) :: option
  Vec :: volume

  type(unstructured_explicit_type), pointer :: explicit_grid
  type(region_explicit_face_type), pointer :: faceset
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icell, count
  PetscInt, allocatable :: int_array(:)
  PetscReal, allocatable :: real_array_2d(:,:)
  PetscErrorCode :: ierr

  explicit_grid => ugrid%explicit_grid
  faceset => region%explicit_faceset
  
  ! convert ids to petsc
  region%cell_ids = region%cell_ids - 1
  call AOApplicationToPetsc(ugrid%ao_natural_to_petsc,size(region%cell_ids), &
                            region%cell_ids,ierr)
  region%cell_ids = region%cell_ids + 1
  
  ! if petsc ids are below global_offset or above global_offset + nlmax, 
  ! they are off processor; otherwise, local

  allocate(int_array(size(region%cell_ids)))
  ! negate off processor ids
  int_array = -999
  ! only if faceset exists
  if (associated(faceset)) then
    allocate(real_array_2d(4,size(region%cell_ids)))
    real_array_2d = -999.d0
  endif
  count = 0
  do icell = 1, size(region%cell_ids)
    if (region%cell_ids(icell) > ugrid%global_offset .and. &
        region%cell_ids(icell) <= ugrid%global_offset + ugrid%nlmax) then
      count = count + 1
      ! local cell id
      int_array(count) = region%cell_ids(icell) - ugrid%global_offset
      if (associated(faceset)) then
        real_array_2d(1,count) = faceset%face_centroids(icell)%x
        real_array_2d(2,count) = faceset%face_centroids(icell)%y
        real_array_2d(3,count) = faceset%face_centroids(icell)%z
        real_array_2d(4,count) = faceset%face_areas(icell)
      endif
    endif
  enddo
  
  deallocate(region%cell_ids)
  allocate(region%cell_ids(count))
  region%cell_ids = int_array(1:count)
  deallocate(int_array)

  if (associated(faceset)) then
    deallocate(faceset%face_centroids)
    deallocate(faceset%face_areas)
    allocate(faceset%face_centroids(count))
    allocate(faceset%face_areas(count))
  
    do icell = 1, count
      faceset%face_centroids(icell)%x = real_array_2d(1,icell)
      faceset%face_centroids(icell)%y = real_array_2d(2,icell)
      faceset%face_centroids(icell)%z = real_array_2d(3,icell)
      faceset%face_areas(icell) = real_array_2d(4,icell)
    enddo
    deallocate(real_array_2d)
  endif
  

  region%num_cells = count
  
  if (region%num_cells == 0) then
    deallocate(region%cell_ids)
    nullify(region%cell_ids)
    if (associated(faceset)) then
      deallocate(faceset%face_centroids)
      nullify(faceset%face_centroids)
      deallocate(faceset%face_areas)
      nullify(faceset%face_areas)
      ! note that have to use full reference
      deallocate(region%explicit_faceset)
      nullify(region%explicit_faceset)
    endif
  endif

#if UGRID_DEBUG
  if (region%num_cells > 0) then
    write(string,*) option%myrank
    string = 'region_faceset_' // trim(region%name) // trim(adjustl(string)) // '.out'
    open(unit=86,file=trim(string))
    if (associated(faceset)) then
      do icell = 1, region%num_cells
        write(86,'(i5,4f7.3)') region%cell_ids(icell), &
                    faceset%face_centroids(icell)%x, &
                    faceset%face_centroids(icell)%y, &
                    faceset%face_centroids(icell)%z, &
                    faceset%face_areas(icell)
      enddo
    else
      do icell = 1, region%num_cells
        write(86,'(i5)') region%cell_ids(icell)
      enddo
    endif
    close(86)
  endif
#endif  

end subroutine GridLocalizeExplicitFaceset

! ************************************************************************** !
!
! GridCopyIntegerArrayToVec: Copies values from an integer array into a 
!                                 PETSc Vec
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !

subroutine GridCopyIntegerArrayToVec(grid, array,vector,num_values)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(grid_type) :: grid
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90( vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90( vector,vec_ptr,ierr)
  
end subroutine GridCopyIntegerArrayToVec

! ************************************************************************** !
!
! GridCopyRealArrayToVec: Copies values from an integer array into a 
!                              PETSc Vec
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyRealArrayToVec(grid,array,vector,num_values)

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
    
  type(grid_type) :: grid
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyRealArrayToVec

! ************************************************************************** !
!
! GridCopyVecToIntegerArray: Copies values from a PETSc Vec to an  
!                                 integer array
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyVecToIntegerArray(grid,array,vector,num_values)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(grid_type) :: grid
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  do i=1,num_values
    if (vec_ptr(i) > 0.d0) then
      array(i) = int(vec_ptr(i)+1.d-4)
    else
      array(i) = int(vec_ptr(i)-1.d-4)
    endif
  enddo
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyVecToIntegerArray

! ************************************************************************** !
!
! GridCopyVecToRealArray: Copies values from a PETSc Vec to an integer 
!                              array
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyVecToRealArray(grid,array,vector,num_values)

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
    
  type(grid_type) :: grid
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  array(1:num_values) = vec_ptr(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyVecToRealArray

! ************************************************************************** !
!
! GridCreateNaturalToGhostedHash: Creates a hash table for looking up the  
!                                 local ghosted id of a natural id, if it 
!                                 exists
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine GridCreateNaturalToGhostedHash(grid,option)

  use Option_module
  use Logging_module  
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: local_ghosted_id, natural_id
  PetscInt :: num_in_hash, num_ids_per_hash, hash_id, id, ierr, hash_id_2
  PetscInt :: max_num_ids_per_hash
  PetscInt, pointer :: hash(:,:,:), temp_hash(:,:,:)

  if (associated(grid%hash)) return

  call PetscLogEventBegin(logging%event_hash_create,ierr)
                          
  max_num_ids_per_hash = 0
  ! initial guess of 10% of ids per hash
  ! must be at least 5 so that reallocation (*1.2) works below
  num_ids_per_hash = max(grid%nlmax/(grid%num_hash_bins/10),5)

  allocate(hash(2,0:num_ids_per_hash,grid%num_hash_bins))
  hash(:,:,:) = 0

  
  do local_ghosted_id = 1, grid%ngmax
    natural_id = grid%nG2A(local_ghosted_id) !nG2A is 1-based
    hash_id = mod(natural_id,grid%num_hash_bins)+1 
    num_in_hash = hash(1,0,hash_id)
    num_in_hash = num_in_hash+1
    if (num_in_hash > max_num_ids_per_hash) max_num_ids_per_hash = num_in_hash
    ! if a hash runs out of space reallocate
    if (num_in_hash > num_ids_per_hash) then 
      allocate(temp_hash(2,0:num_ids_per_hash,grid%num_hash_bins))
      ! copy old hash
      temp_hash(1:2,0:num_ids_per_hash,1:grid%num_hash_bins) = &
                             hash(1:2,0:num_ids_per_hash,1:grid%num_hash_bins)
      deallocate(hash)
      ! recompute hash 20% larger
      num_ids_per_hash = int(dble(num_ids_per_hash)*1.2)
      allocate(hash(2,0:num_ids_per_hash,grid%num_hash_bins))
      ! copy old to new
      do hash_id_2 = 1, grid%num_hash_bins
        do id = 1, temp_hash(1,0,hash_id_2)
          hash(1:2,id,hash_id_2) = temp_hash(1:2,id,hash_id_2)
        enddo
        hash(1,0,hash_id_2) = temp_hash(1,0,hash_id_2)
      enddo
      deallocate(temp_hash)
    endif
    hash(1,0,hash_id) = num_in_hash
    hash(1,num_in_hash,hash_id) = natural_id
    hash(2,num_in_hash,hash_id) = local_ghosted_id
  enddo

  grid%hash => hash
  
!  call GridPrintHashTable(grid)
  call MPI_Allreduce(max_num_ids_per_hash,num_in_hash,ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  write(option%io_buffer,'("max_num_ids_per_hash: ",i5)') num_in_hash
  call printMsg(option)

  call PetscLogEventEnd(logging%event_hash_create,ierr)

end subroutine GridCreateNaturalToGhostedHash

! ************************************************************************** !
!
! GetLocalIdFromNaturalId: Returns the local id corresponding to a natural
!                          id or 0, if the natural id is off-processor
! WARNING: Extremely inefficient for large jobs
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalIdFromNaturalId(grid,natural_id)

  implicit none

  type(grid_type) :: grid

  PetscInt :: natural_id, local_id
  
  do local_id = 1, grid%nlmax
    if (natural_id == grid%nG2A(grid%nL2G(local_id))) then
      GridGetLocalIdFromNaturalId = local_id
      return
    endif
  enddo
  GridGetLocalIdFromNaturalId = 0

end function GridGetLocalIdFromNaturalId

! ************************************************************************** !
!
! GridGetLocalGhostedIdFromNatId: Returns the local ghosted id corresponding 
!                                 to a natural id or 0, if the natural id 
!                                 is off-processor
! WARNING: Extremely inefficient for large jobs
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalGhostedIdFromNatId(grid,natural_id)

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: local_ghosted_id
  
  do local_ghosted_id = 1, grid%ngmax
    !geh: nG2A is 1-based
    if (natural_id == grid%nG2A(local_ghosted_id)) then
      GridGetLocalGhostedIdFromNatId = local_ghosted_id
      return 
    endif
  enddo
  GridGetLocalGhostedIdFromNatId = 0

end function GridGetLocalGhostedIdFromNatId

! ************************************************************************** !
!
! GridGetLocalGhostedIdFromHash: Returns the local ghosted id of a natural 
!                                id, if it exists.  Otherwise 0 is returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalGhostedIdFromHash(grid,natural_id)

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: hash_id, id

  GridGetLocalGhostedIdFromHash = 0
  hash_id = mod(natural_id,grid%num_hash_bins)+1 
  do id = 1, grid%hash(1,0,hash_id)
    if (grid%hash(1,id,hash_id) == natural_id) then
      GridGetLocalGhostedIdFromHash = grid%hash(2,id,hash_id)
      return
    endif
  enddo

end function GridGetLocalGhostedIdFromHash

! ************************************************************************** !
!
! GridDestroyHashTable: Deallocates the hash table
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine GridDestroyHashTable(grid)

  implicit none

  type(grid_type), pointer :: grid
  
  if (associated(grid%hash)) deallocate(grid%hash)
  
#ifdef DASVYAT
  if (associated(grid%faces)) deallocate(grid%faces)
  nullify(grid%faces)
#endif

  nullify(grid%hash)
  grid%num_hash_bins = 100

end subroutine GridDestroyHashTable

! ************************************************************************** !
!
! UnstructGridPrintHashTable: Prints the hashtable for viewing
! author: Glenn Hammond
! date: 03/09/07
!
! ************************************************************************** !
subroutine GridPrintHashTable(grid)

  implicit none

  type(grid_type) :: grid
  
  PetscInt :: ihash, id, fid

  fid = 87 
  open(fid,file='hashtable.dat',action='write')
  do ihash=1,grid%num_hash_bins
    write(fid,'(a4,i3,a,i5,a2,x,200(i6,x))') 'Hash',ihash,'(', &
                         grid%hash(1,0,ihash), &
                         '):', &
                         (grid%hash(1,id,ihash),id=1,grid%hash(1,0,ihash))
  enddo
  close(fid)

end subroutine GridPrintHashTable

! ************************************************************************** !
!
! GridGetNeighbors: Returns an array of neighboring cells
! author: Glenn Hammond
! date: 01/28/11
!
! ************************************************************************** !
subroutine GridGetGhostedNeighbors(grid,ghosted_id,stencil_type, &
                                   stencil_width_i,stencil_width_j, &
                                   stencil_width_k,x_count,y_count, &
                                   z_count, &
                                   ghosted_neighbors,option)

  use Option_module

  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscInt :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: x_count
  PetscInt :: y_count
  PetscInt :: z_count
  PetscInt :: ghosted_neighbors(*)
  
  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call StructGridGetGhostedNeighbors(grid%structured_grid, &
                                         ghosted_id,stencil_type, &
                                         stencil_width_i, &
                                         stencil_width_j,stencil_width_k, &
                                         x_count,y_count,z_count, &
                                         ghosted_neighbors,option)
    case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID) 
      option%io_buffer = 'GridGetNeighbors not currently supported for ' // &
        'unstructured grids.'
      call printErrMsg(option)
  end select

end subroutine GridGetGhostedNeighbors

! ************************************************************************** !
!
! GridGetNeighborsWithCorners: Returns an array of neighboring cells along with corner
! cells
! author: Satish Karra, LANL
! date: 02/19/12
!
! ************************************************************************** !
subroutine GridGetGhostedNeighborsWithCorners(grid,ghosted_id,stencil_type, &
                                   stencil_width_i,stencil_width_j, &
                                   stencil_width_k,icount, &
                                   ghosted_neighbors,option)

  use Option_module

  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscInt :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: icount
  PetscInt :: ghosted_neighbors(*)
  
  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call StructGridGetGhostedNeighborsCorners(grid%structured_grid, &
                                         ghosted_id,stencil_type, &
                                         stencil_width_i, &
                                         stencil_width_j, &
                                         stencil_width_k, &
                                         icount, &
                                         ghosted_neighbors,option)
    case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID) 
      option%io_buffer = 'GridGetNeighbors not currently supported for ' // &
        'unstructured grids.'
      call printErrMsg(option)
  end select

end subroutine GridGetGhostedNeighborsWithCorners


! ************************************************************************** !
!
! GridDestroy: Deallocates a grid
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine GridDestroy(grid)

  implicit none
  
  type(grid_type), pointer :: grid
  PetscErrorCode :: ierr
  PetscInt :: ghosted_id
    
  if (.not.associated(grid)) return
      
  if (associated(grid%nL2G)) deallocate(grid%nL2G)
  nullify(grid%nL2G)
  if (associated(grid%nG2L)) deallocate(grid%nG2L)
  nullify(grid%nG2L)
  if (associated(grid%nG2A)) deallocate(grid%nG2A)
  nullify(grid%nG2A)
  if (associated(grid%nG2P)) deallocate(grid%nG2P)
  nullify(grid%nG2P)

  !Note: Destroying for ghosted_level<TWO_INTEGER assumes that max_stencil_width
  !      was TWO_INTEGER.
  if (associated(grid%dispT)) then
    do ghosted_id=1,grid%ngmax
      if(grid%ghosted_level(ghosted_id)<TWO_INTEGER) then
        call MatDestroy(grid%dispT(ghosted_id),ierr)
      endif
    enddo
  endif
  nullify(grid%dispT)
  if (associated(grid%Minv)) then
    do ghosted_id=1,grid%ngmax
      if(grid%ghosted_level(ghosted_id)<TWO_INTEGER) then
        call MatDestroy(grid%Minv(ghosted_id),ierr)
      endif
    enddo
  endif
  nullify(grid%Minv)
  if (associated(grid%ghosted_level)) deallocate(grid%ghosted_level)
  nullify(grid%ghosted_level)
  if (associated(grid%jacfac)) deallocate(grid%jacfac)
  nullify(grid%jacfac)
  if (associated(grid%bnd_cell)) deallocate(grid%bnd_cell)
  nullify(grid%bnd_cell)

#ifdef DASVYAT
  if (associated(grid%fL2G)) deallocate(grid%fL2G)
  nullify(grid%fL2G)
  if (associated(grid%fG2L)) deallocate(grid%fG2L)
  nullify(grid%fG2L)
  if (associated(grid%fG2P)) deallocate(grid%fG2P)
  nullify(grid%fG2P)
  if (associated(grid%fL2P)) deallocate(grid%fL2P)
  nullify(grid%fL2P)
  if (associated(grid%fL2B)) deallocate(grid%fL2B)
  nullify(grid%fL2B)

  if (grid%e2f /= 0) Call VecDestroy(grid%e2f, ierr)
  if (grid%e2n /= 0) Call VecDestroy(grid%e2n, ierr)
  if (grid%e2n_LP /= 0) Call VecDestroy(grid%e2n_LP, ierr)

  call MFDAuxDestroy(grid%MFD)
#endif

  if (associated(grid%x)) deallocate(grid%x)
  nullify(grid%x)
  if (associated(grid%y)) deallocate(grid%y)
  nullify(grid%y)
  if (associated(grid%z)) deallocate(grid%z)
  nullify(grid%z)
  
  if (associated(grid%hash)) call GridDestroyHashTable(grid)
  
  call UGridDestroy(grid%unstructured_grid)    
  call StructGridDestroy(grid%structured_grid)
                                           
  call ConnectionDestroyList(grid%internal_connection_set_list)

end subroutine GridDestroy
      
! ************************************************************************** !
!
! GridIndexToCellID: Returns the local grid cell id of a Petsc Vec index
! author: Glenn Hammond
! date: 01/07/09
!
! ************************************************************************** !
function GridIndexToCellID(vec,index,grid,vec_type)

  implicit none
  
  Vec :: vec
  PetscInt :: index
  type(grid_type) :: grid
  PetscInt :: vec_type
  
  PetscInt :: GridIndexToCellID
  
  PetscInt :: low
  PetscInt :: high
  PetscInt :: ndof
  PetscInt :: cell_id
  PetscErrorCode :: ierr

  
  cell_id = -1
  call VecGetOwnershipRange(vec,low,high,ierr)
  call VecGetBlockSize(vec,ndof,ierr)
  if (index >= low .and. index < high) then
    cell_id = (index-low)/ndof+1
    if (vec_type == GLOBAL) then
      cell_id = grid%nG2A(grid%nL2G(cell_id))
    else if (vec_type == LOCAL) then
      cell_id = grid%nG2A(cell_id) !nG2A is 1-based
    endif
  endif
  
  call MPI_Allreduce(cell_id,GridIndexToCellID,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,PETSC_COMM_WORLD,ierr)
                     
end function GridIndexToCellID

! ************************************************************************** !
!> This routine computes the indices of neighbours for all ghosted cells
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 08/24/12
! ************************************************************************** !
subroutine GridComputeNeighbors(grid,is_bnd_vec,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  Vec :: is_bnd_vec
  type(option_type) :: option

  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call StructGridComputeNeighbors(grid%structured_grid,grid%nG2L,is_bnd_vec,option)
    case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID) 
      option%io_buffer = 'GridComputeNeighbors not currently supported for ' // &
        'unstructured grids.'
      call printErrMsg(option)
  end select

end subroutine GridComputeNeighbors

! ************************************************************************** !
!> This routine resticts regions to cells local to processor when the region
!! was defined using a BLOCK from inputfile.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/04/12
! ************************************************************************** !
subroutine GridLocalizeRegionFromBlock(grid,region,option)

  use Option_module
  use Region_module

  implicit none
  
  type(region_type), pointer :: region
  type(grid_type), pointer   :: grid
  type(option_type)          :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscReal :: del_x, del_y, del_z
  PetscInt :: iflag
  PetscBool :: same_point
  PetscErrorCode :: ierr

  if(grid%itype /= STRUCTURED_GRID .and. &
     grid%itype /= STRUCTURED_GRID_MIMETIC) then
     
     option%io_buffer='Region definition using BLOCK is only supported for ' //&
       ' structured grids'
     call printErrMsg(option)
  endif
  
  ! convert indexing from global (entire domain) to local processor
  region%i1 = region%i1 - grid%structured_grid%lxs
  region%i2 = region%i2 - grid%structured_grid%lxs
  region%j1 = region%j1 - grid%structured_grid%lys
  region%j2 = region%j2 - grid%structured_grid%lys
  region%k1 = region%k1 - grid%structured_grid%lzs
  region%k2 = region%k2 - grid%structured_grid%lzs
          
  ! clip region to within local processor domain
  region%i1 = max(region%i1,1)
  region%i2 = min(region%i2,grid%structured_grid%nlx)
  region%j1 = max(region%j1,1)
  region%j2 = min(region%j2,grid%structured_grid%nly)
  region%k1 = max(region%k1,1)
  region%k2 = min(region%k2,grid%structured_grid%nlz)
   
  count = 0  
  if (region%i1 <= region%i2 .and. &
      region%j1 <= region%j2 .and. &
      region%k1 <= region%k2) then
    region%num_cells = (region%i2-region%i1+1)* &
                       (region%j2-region%j1+1)* &
                       (region%k2-region%k1+1)
    allocate(region%cell_ids(region%num_cells))
    if (region%iface /= 0) then
      allocate(region%faces(region%num_cells))
      region%faces = region%iface
    endif
    region%cell_ids = 0
    do k=region%k1,region%k2
      do j=region%j1,region%j2
        do i=region%i1,region%i2
          count = count + 1
          region%cell_ids(count) = &
                 i + (j-1)*grid%structured_grid%nlx + &
                 (k-1)*grid%structured_grid%nlxy
        enddo
      enddo
    enddo
!   if (region%num_cells > 0) then
!     region%coordinates(1)%x = grid%x(region%cell_ids(ONE_INTEGER))
!     region%coordinates(1)%y = grid%y(region%cell_ids(ONE_INTEGER))
!     region%coordinates(1)%z = grid%z(region%cell_ids(ONE_INTEGER))
!   endif
  else
    region%num_cells = 0
  endif

  if (count /= region%num_cells) then
    option%io_buffer = 'Mismatch in number of cells in block region'
    call printErrMsg(option)
  endif

end subroutine GridLocalizeRegionFromBlock

! ************************************************************************** !
!> This routine resticts regions to cells local to processor when the region
!! was defined using COORDINATES from inputfile.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/04/12
! ************************************************************************** !
subroutine GridLocalizeRegionFromCoordinates(grid,region,option)

  use Option_module
  use Region_module

  implicit none
  
  type(region_type), pointer :: region
  type(grid_type), pointer   :: grid
  type(option_type)          :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscReal :: del_x, del_y, del_z
  PetscInt :: iflag
  PetscBool :: same_point
  PetscErrorCode :: ierr

  iflag = 0
  if (size(region%coordinates) > TWO_INTEGER) then
    option%io_buffer = 'GridLocalizeRegions: more than 2 coordinates' // &
                       ' not supported in region object'
    call printErrMsg(option)
  endif

  same_point = PETSC_FALSE
  ! if two coordinates, determine whether they are the same point
  if (size(region%coordinates) == TWO_INTEGER) then
    if (dabs(region%coordinates(ONE_INTEGER)%x - &
             region%coordinates(TWO_INTEGER)%x) < tol .and. &
        dabs(region%coordinates(ONE_INTEGER)%y - &
             region%coordinates(TWO_INTEGER)%y) < tol .and. &
        dabs(region%coordinates(ONE_INTEGER)%z - &
             region%coordinates(TWO_INTEGER)%z) < tol) then
      same_point = PETSC_TRUE
    endif
  endif
   
  ! treat two identical coordinates the same as a single coordinate
  if (size(region%coordinates) == ONE_INTEGER .or. same_point) then
    if (region%coordinates(ONE_INTEGER)%x >= grid%x_min_global .and. &
        region%coordinates(ONE_INTEGER)%x <= grid%x_max_global .and. &
        region%coordinates(ONE_INTEGER)%y >= grid%y_min_global .and. &
        region%coordinates(ONE_INTEGER)%y <= grid%y_max_global .and. &
        region%coordinates(ONE_INTEGER)%z >= grid%z_min_global .and. &
        region%coordinates(ONE_INTEGER)%z <= grid%z_max_global) then
      ! If a point is on the corner of 4 or 8 patches in AMR, the region
      ! will be assigned to all 4/8 patches...a problem.  To avoid this, 
      ! we are going to perturb all point coordinates slightly upwind, as
      ! long as they are not on a global boundary (i.e. boundary condition)
      ! -- shift the coorindate slightly upwind
      x_shift = region%coordinates(ONE_INTEGER)%x - &
                pert*(grid%x_max_global-grid%x_min_global)
      y_shift = region%coordinates(ONE_INTEGER)%y - &
                pert*(grid%y_max_global-grid%y_min_global)
      z_shift = region%coordinates(ONE_INTEGER)%z - &
                pert*(grid%z_max_global-grid%z_min_global)
      ! if the coodinate is shifted out of the global domain or 
      ! onto an exterior edge, set it back to the original value
      if (x_shift - grid%x_min_global < tol) &
        x_shift = region%coordinates(ONE_INTEGER)%x
      if (y_shift - grid%y_min_global < tol) &
        y_shift = region%coordinates(ONE_INTEGER)%y
      if (z_shift - grid%z_min_global < tol) &
        z_shift = region%coordinates(ONE_INTEGER)%z
      select case(grid%itype)
        case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
          call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                              x_shift,y_shift,z_shift, &
                                              i,j,k)
          if (i > 0 .and. j > 0 .and. k > 0) then
            region%num_cells = 1
            allocate(region%cell_ids(region%num_cells))
            if (region%iface /= 0) then
              allocate(region%faces(region%num_cells))
              region%faces = region%iface
            endif
            region%cell_ids = 0
            region%cell_ids(1) = i + (j-1)*grid%structured_grid%nlx + &
                                (k-1)*grid%structured_grid%nlxy
          else
            region%num_cells = 0
          endif
          ! the next test as designed will only work on a uniform grid
          call MPI_Allreduce(region%num_cells,count, &
                              ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                              option%mycomm,ierr)   
          if (count == 0) then
            write(option%io_buffer,*) 'Region: (coord)', &
                  region%coordinates(ONE_INTEGER)%x, &
                  region%coordinates(ONE_INTEGER)%y, &
                  region%coordinates(ONE_INTEGER)%z, &
                  ' not found in global domain.', count
            call printErrMsg(option)
          else if (count > 1) then
            write(option%io_buffer,*) 'Region: (coord)', &
                  region%coordinates(ONE_INTEGER)%x, &
                  region%coordinates(ONE_INTEGER)%y, &
                  region%coordinates(ONE_INTEGER)%z, &
                  ' duplicated across ', count, &
                  ' procs in global domain.'
            call printErrMsg(option)
          endif
        case(IMPLICIT_UNSTRUCTURED_GRID)
          !geh: must check each cell individually
          call UGridGetCellFromPoint(region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     grid%unstructured_grid,option,local_id)
          if (local_id > 0) then
            region%num_cells = 1
            allocate(region%cell_ids(region%num_cells))
            region%cell_ids(1) = local_id
          else
            region%num_cells = 0
          endif
        case(POLYHEDRA_UNSTRUCTURED_GRID)
           option%io_buffer = 'add code POLYHDERA in GridLocalizeRegionFromCoordinates'
           call printErrMsg(option)
      end select
    endif
  else ! 2 coordinates
    x_min = min(region%coordinates(ONE_INTEGER)%x, &
                region%coordinates(TWO_INTEGER)%x)
    x_max = max(region%coordinates(ONE_INTEGER)%x, &
                region%coordinates(TWO_INTEGER)%x)
    y_min = min(region%coordinates(ONE_INTEGER)%y, &
                region%coordinates(TWO_INTEGER)%y)
    y_max = max(region%coordinates(ONE_INTEGER)%y, &
                region%coordinates(TWO_INTEGER)%y)
    z_min = min(region%coordinates(ONE_INTEGER)%z, &
                region%coordinates(TWO_INTEGER)%z)
    z_max = max(region%coordinates(ONE_INTEGER)%z, &
                region%coordinates(TWO_INTEGER)%z)
                
    if (grid%itype == STRUCTURED_GRID .or. &
        grid%itype == STRUCTURED_GRID_MIMETIC) then
      ! shift box slightly inward
      x_shift = 1.d-8*(grid%x_max_global-grid%x_min_global)
      x_min = x_min+x_shift            
      x_max = x_max-x_shift
      y_shift = 1.d-8*(grid%y_max_global-grid%y_min_global)
      y_min = y_min+y_shift            
      y_max = y_max-y_shift
      z_shift = 1.d-8*(grid%z_max_global-grid%z_min_global)
      z_min = z_min+z_shift            
      z_max = z_max-z_shift
         
      ! if plane or line, ensure it is within the grid cells     
      if (x_max-x_min < 1.d-10) then
        x_max = region%coordinates(ONE_INTEGER)%x
        x_shift = 1.d-8*(grid%x_max_global-grid%x_min_global)
        if (region%iface == WEST_FACE) then
          x_max = x_max + x_shift
        else if (region%iface == EAST_FACE) then
          x_max = x_max - x_shift
        ! otherwise, shift upwind, unless at upwind physical boundary
        else
          if (x_max > grid%x_min_global + x_shift) then
            x_max = x_max - x_shift
          else
            x_max = x_max + x_shift
          endif
        endif
        x_min = x_max
      endif
      if (y_max-y_min < 1.d-10) then
        y_max = region%coordinates(ONE_INTEGER)%y
        y_shift = 1.d-8*(grid%y_max_global-grid%y_min_global)
        if (region%iface == SOUTH_FACE) then
          y_max = y_max + y_shift
        else if (region%iface == NORTH_FACE) then
          y_max = y_max - y_shift
        ! otherwise, shift upwind, unless at upwind physical boundary
        else
          if (y_max > grid%y_min_global + y_shift) then
            y_max = y_max - y_shift
          else
            y_max = y_max + y_shift
          endif
        endif
        y_min = y_max
      endif
      if (z_max-z_min < 1.d-10) then
        z_max = region%coordinates(ONE_INTEGER)%z
        z_shift = 1.d-8*(grid%z_max_global-grid%z_min_global)
        if (region%iface == BOTTOM_FACE) then
          z_max = z_max + z_shift
        else if (region%iface == TOP_FACE) then
          z_max = z_max - z_shift
        ! otherwise, shift upwind, unless at upwind physical boundary
        else
          if (z_max > grid%z_min_global + z_shift) then
            z_max = z_max - z_shift
          else
            z_max = z_max + z_shift
          endif
        endif
        z_min = z_max
      endif
    endif   
             
    ! ensure overlap
    if (x_min <= grid%x_max_local .and. &
        x_max >= grid%x_min_local .and. &
        y_min <= grid%y_max_local .and. &
        y_max >= grid%y_min_local .and. &
        z_min <= grid%z_max_local .and. &
        z_max >= grid%z_min_local) then
        
      ! get I,J,K bounds
      select case(grid%itype)
        case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
          ! local, non-ghosted i,j,k's are returned
          call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                    max(x_min,grid%x_min_local+x_shift), &
                                    max(y_min,grid%y_min_local+y_shift), &
                                    max(z_min,grid%z_min_local+z_shift), &
                                              i_min,j_min,k_min)
          call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                    min(x_max,grid%x_max_local-x_shift), &
                                    min(y_max,grid%y_max_local-y_shift), &
                                    min(z_max,grid%z_max_local-z_shift), &
                                              i_max,j_max,k_max)
          if (i_min > 0 .and. j_min > 0 .and. k_min > 0 .and. &
              i_max > 0 .and. j_max > 0 .and. k_max > 0) then
            region%num_cells = (i_max-i_min+1)*(j_max-j_min+1)*(k_max-k_min+1)
            allocate(region%cell_ids(region%num_cells))
            if (region%iface /= 0) then
              allocate(region%faces(region%num_cells))
              region%faces = region%iface
            endif
            region%cell_ids = 0
            count = 0
            do k = k_min, k_max
              do j = j_min, j_max
                do i = i_min, i_max
                  count = count+1
                  region%cell_ids(count) = i + (j-1)*grid%structured_grid%nlx + &
                                      (k-1)*grid%structured_grid%nlxy
                enddo
              enddo
            enddo
          else
            iflag = 1
          endif
        case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID,POLYHEDRA_UNSTRUCTURED_GRID)
          del_x = x_max-x_min
          del_y = y_max-y_min
          del_z = z_max-z_min
          ! 3D box
          if (del_x > 1.d-10 .and. &
              del_y > 1.d-10 .and. &
              del_z > 1.d-10) then
            ! geh: if the coordinates define a 3D box, add all cell centers
            ! that reside within the box
            count = 0
            do local_id = 1, grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              if (grid%x(ghosted_id) >= x_min .and. &
                  grid%x(ghosted_id) <= x_max .and. &
                  grid%y(ghosted_id) >= y_min .and. &
                  grid%y(ghosted_id) <= y_max .and. &
                  grid%z(ghosted_id) >= z_min .and. &
                  grid%z(ghosted_id) <= z_max) then
                count = count + 1
              endif
            enddo
            allocate(region%cell_ids(count))
            region%cell_ids = 0
            count = 0
            do local_id = 1, grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              if (grid%x(ghosted_id) >= x_min .and. &
                  grid%x(ghosted_id) <= x_max .and. &
                  grid%y(ghosted_id) >= y_min .and. &
                  grid%y(ghosted_id) <= y_max .and. &
                  grid%z(ghosted_id) >= z_min .and. &
                  grid%z(ghosted_id) <= z_max) then
                count = count + 1
                region%cell_ids(count) = local_id
              endif
            enddo
            region%num_cells = count
          ! 2D plane
          elseif ((del_x < 1.d-10 .and. del_y > 1.d-10 .and. &
                   del_z > 1.d-10) .or. &
                  (del_x > 1.d-10 .and. del_y < 1.d-10 .and. &
                   del_z > 1.d-10) .or. &
                  (del_x > 1.d-10 .and. del_y > 1.d-10 .and. &
                   del_z < 1.d-10)) then
            if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID .or. &
                grid%itype == EXPLICIT_UNSTRUCTURED_GRID) then
              call UGridGetCellsInRectangle(x_min,x_max,y_min,y_max, &
                                            z_min,z_max, &
                                            grid%unstructured_grid,option, &
                                            region%num_cells,region%cell_ids, &
                                            region%faces)
            else
              call PolyhedraUGridGetCellsInRectangle(x_min,x_max,y_min,y_max, &
                                                     z_min,z_max, &
                                                     grid%unstructured_grid,option, &
                                                     region%num_cells,region%cell_ids, &
                                                     region%faces)
            endif

          endif
      end select
    endif

    call MPI_Allreduce(iflag,i,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX, &
                       option%mycomm,ierr)
    iflag = i
    if (iflag > 0) then
      option%io_buffer = 'GridLocalizeRegions, between two points'
      call printErrMsg(option)
    endif
  endif

end subroutine GridLocalizeRegionFromCoordinates

! ************************************************************************** !
!> This routine computes the following:
!! - Displacement matrix (D), and
!! - Inverse of D^T*D
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 11/20/12
! ************************************************************************** !
subroutine GridComputeMinv(grid,max_stencil_width,option)

  use Option_module
  use Utility_module

  implicit none

  type(grid_type) :: grid
  type(option_type) :: option
  PetscInt :: max_stencil_width

  PetscInt, pointer :: cell_neighbors(:,:)
  PetscInt :: ghosted_id,nid
  PetscReal :: dx,dy,dz
  PetscErrorCode :: ierr
  Mat :: A,B
  PetscScalar, pointer :: xx_v(:,:)
  PetscScalar :: b_v(3)
  PetscInt :: cols(3), ncol
  PetscInt :: INDX(3)
  PetscInt :: D,ii,jj
  PetscReal :: disp_mat(3,3)
  PetscReal :: identity(3)
  PetscInt :: max_neighbors
  Vec :: iden_vec,C
  PetscReal, pointer :: v_p(:)

  select case(grid%itype)
    case(STRUCTURED_GRID)
      cell_neighbors => grid%structured_grid%cell_neighbors
    case(UNSTRUCTURED_GRID)
      option%io_buffer='GridComputeMinv() not implemented for unstructured grid.'
      call printErrMsg(option)
  end select

  allocate(grid%dispT(grid%ngmax))
  allocate(grid%Minv(grid%ngmax))

  max_neighbors = maxval(cell_neighbors(0,:))
  allocate(grid%jacfac(grid%ngmax,0:max_neighbors,THREE_INTEGER))

  call VecCreateSeq(PETSC_COMM_SELF,THREE_INTEGER,C,ierr)

  do ghosted_id = 1,grid%ngmax

    if(grid%ghosted_level(ghosted_id)<max_stencil_width) then

      ! Create the disp matrix
      call MatCreate(MPI_COMM_SELF,grid%dispT(ghosted_id),ierr)
      call MatSetSizes(grid%dispT(ghosted_id),3,cell_neighbors(0,ghosted_id),3,cell_neighbors(0,ghosted_id),ierr)
      call MatSetType(grid%dispT(ghosted_id),MATSEQDENSE,ierr)
      call MatSetUp(grid%dispT(ghosted_id),ierr)
      do nid = 1,cell_neighbors(0,ghosted_id)
        dx = grid%x(cell_neighbors(nid,ghosted_id)) - grid%x(ghosted_id)
        dy = grid%y(cell_neighbors(nid,ghosted_id)) - grid%y(ghosted_id)
        dz = grid%z(cell_neighbors(nid,ghosted_id)) - grid%z(ghosted_id)
        call MatSetValue(grid%dispT(ghosted_id),0,nid-1,dx,INSERT_VALUES,ierr)
        call MatSetValue(grid%dispT(ghosted_id),1,nid-1,dy,INSERT_VALUES,ierr)
        call MatSetValue(grid%dispT(ghosted_id),2,nid-1,dz,INSERT_VALUES,ierr)
      enddo
      call MatAssemblyBegin(grid%dispT(ghosted_id),MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(  grid%dispT(ghosted_id),MAT_FINAL_ASSEMBLY,ierr)

      ! Compute transpose of disp matrix
      call MatTranspose(grid%dispT(ghosted_id),MAT_INITIAL_MATRIX, &
                        A,ierr)

      ! B = disp_mat^T * disp_mat
      call MatMatMult(grid%dispT(ghosted_id),A, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL_PRECISION,B,ierr)

      ! Pack the values of B in disp_mat for obtaining the inverse of matrix
      do ii=0,2
        call MatGetRow(B,ii,ncol,cols,b_v,ierr)
        disp_mat(ii+1,:) = b_v(:)
        call MatRestoreRow(B,ii,ncol,cols,b_v,ierr)
      enddo

      ! LU decomposition of disp_mat
      call ludcmp(disp_mat,THREE_INTEGER,INDX,D)

      ! Save the inverse matrix
      call MatCreate(MPI_COMM_SELF,grid%Minv(ghosted_id),ierr)
      call MatSetSizes(grid%Minv(ghosted_id),3,3,3,3,ierr)
      call MatSetType(grid%Minv(ghosted_id),MATSEQDENSE,ierr)
      call MatSetUp(grid%Minv(ghosted_id),ierr)

      ! Find inverse matrix column-by-column
      do ii=1,3
        identity = 0
        identity(ii) = 1
        call lubksb(disp_mat,THREE_INTEGER,INDX,identity)
        call MatSetValue(grid%Minv(ghosted_id),0,ii-1,identity(1),INSERT_VALUES,ierr)
        call MatSetValue(grid%Minv(ghosted_id),1,ii-1,identity(2),INSERT_VALUES,ierr)
        call MatSetValue(grid%Minv(ghosted_id),2,ii-1,identity(3),INSERT_VALUES,ierr)
      enddo

      call MatAssemblyBegin(grid%Minv(ghosted_id),MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(  grid%Minv(ghosted_id),MAT_FINAL_ASSEMBLY,ierr)

      call MatDestroy(A,ierr)
      call MatDestroy(B,ierr)

      ! Save values used in Jacobian computation

      ! Compute Minv * dispT
      call MatMatMult(grid%Minv(ghosted_id),grid%dispT(ghosted_id), &
                      MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL_PRECISION,A,ierr)

      call VecCreateSeq(PETSC_COMM_SELF,cell_neighbors(0,ghosted_id),iden_vec,ierr)

      ! Save for 'ghosted_id'
      call VecGetArrayF90(iden_vec,v_p,ierr)
      v_p = -1.d0
      call VecRestoreArrayF90(iden_vec,v_p,ierr)
      call MatMult(A,iden_vec,C,ierr)

      call VecGetArrayF90(C,v_p,ierr)
      grid%jacfac(ghosted_id,0,:) = v_p(:)
      call VecGetArrayF90(C,v_p,ierr)

      ! Save for neighbors of 'ghosted_id'
      do nid = 1,cell_neighbors(0,ghosted_id)

        call VecGetArrayF90(iden_vec,v_p,ierr)
        v_p = 0.d0
        v_p(nid) = 1.d0
        call VecRestoreArrayF90(iden_vec,v_p,ierr)

        call MatMult(A,iden_vec,C,ierr)

        call VecGetArrayF90(C,v_p,ierr)
        grid%jacfac(ghosted_id,nid,:) = v_p(:)
        call VecGetArrayF90(C,v_p,ierr)

      enddo

      call VecDestroy(iden_vec,ierr)
      call MatDestroy(A,ierr)

    endif
  enddo

  call VecDestroy(C,ierr)

end subroutine GridComputeMinv

! ************************************************************************** !
!> This routine saves information regarding a cell being boundary or 
!! interior cell.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/19/12
! ************************************************************************** !
subroutine GridSaveBoundaryCellInfo(grid,is_bnd_vec,option)

  use Option_module

  implicit none

  type(grid_type) :: grid
  Vec :: is_bnd_vec
  type(option_type) :: option
  
  PetscInt:: ghosted_id
  PetscScalar,pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  allocate(grid%bnd_cell(grid%ngmax))

  call VecGetArrayF90(is_bnd_vec,vec_ptr,ierr)

  do ghosted_id=1,grid%ngmax
    if(vec_ptr(ghosted_id)==0.d0) then
      grid%bnd_cell(ghosted_id) = PETSC_FALSE
    else
      grid%bnd_cell(ghosted_id) = PETSC_TRUE
    endif
  enddo

  call VecRestoreArrayF90(is_bnd_vec,vec_ptr,ierr)

end subroutine GridSaveBoundaryCellInfo

! ************************************************************************** !
!> This routine populates faces of a unstructured grid for MIMETIC
!! discretization. It does the following:
!! 1) Calculates nlmax_faces, ngmax_faces
!! 2) If boundary faces are present, adds boundary connection in 
!!    grid%boundary_connection_set_list
!! 3) Lastly, save information about ghosted faces: faces(ngmax_faces)
!!
!! Note: This subroutine performs functions of StructGridComputeInternConnect(),
!!       StructGridComputeBoundConnect(), and GridPopulateFaces() for a
!!       STRUCTURED_GRID_MIMETIC.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/29/13
! ************************************************************************** !
subroutine GridPopulateFacesForUGrid(grid,option)

  use Option_module
  use Unstructured_Cell_module
  use Connection_module

  implicit none

  type(grid_type) :: grid
  type(option_type) :: option

#if defined(MFD_UGRID)
  type(face_type), pointer :: faces(:)
  type(unstructured_grid_type),pointer :: ugrid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(connection_set_type), pointer :: bnd_connections
  type(connection_set_type), pointer :: conn

  PetscInt :: iconn
  PetscInt :: nconn
  PetscInt :: icell
  PetscInt :: iface
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ghosted_id_up
  PetscInt :: ghosted_id_dn
  PetscInt :: face_id
  PetscInt :: nfaces
  PetscInt :: nfaces_intrn_loc
  PetscInt :: nfaces_intrn_nonloc
  PetscInt :: nfaces_bnd
  PetscInt :: offset
  PetscInt,pointer::int_array(:)
  PetscReal, pointer :: vec_ptr(:)

  Vec :: proc_id_loc
  Vec :: proc_id_nat
  Vec :: proc_id_ghosted
  VecScatter :: vec_scat

  IS :: is_tmp1, is_tmp2

  PetscErrorCode :: ierr
  PetscViewer :: viewer

  ugrid => grid%unstructured_grid

  nfaces_bnd = 0
  nfaces_intrn_loc = 0
  nfaces_intrn_nonloc = 0

  ! Compute number of boundary faces
  do local_id = 1, ugrid%nlmax
    nfaces = UCellGetNFaces(ugrid%cell_type(local_id),option)
    do iface = 1, nfaces
      face_id = ugrid%cell_to_face_ghosted(iface,local_id)
      if (ugrid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        nfaces_bnd = nfaces_bnd + 1
      endif
    enddo
  enddo

  ! For local+ghost cells, determine processor id on which the cell resides
  call VecCreateMPI(option%mycomm,ugrid%nlmax,PETSC_DETERMINE,proc_id_loc,ierr)
  call VecCreateMPI(option%mycomm,ugrid%nlmax,PETSC_DETERMINE,proc_id_nat,ierr)
  call VecCreateMPI(option%mycomm,ugrid%ngmax,PETSC_DETERMINE,proc_id_ghosted,ierr)

  call VecGetArrayF90(proc_id_loc,vec_ptr,ierr)
  do local_id = 1,ugrid%nlmax
    vec_ptr(local_id)=option%myrank
  enddo
  call VecRestoreArrayF90(proc_id_loc,vec_ptr,ierr)

  allocate(int_array(ugrid%nlmax))
  do local_id=1,ugrid%nlmax
    int_array(local_id)=local_id-1+ugrid%global_offset
  enddo
  call ISCreateGeneral(option%mycomm,ugrid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr)

  do local_id=1,ugrid%nlmax
    int_array(local_id)=grid%nG2A(grid%nL2G(local_id))-1
  enddo
  call ISCreateGeneral(option%mycomm,ugrid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr)
  deallocate(int_array)

  call VecScatterCreate(proc_id_loc,is_tmp1,proc_id_nat,is_tmp2,vec_scat,ierr)
  call ISDestroy(is_tmp1,ierr)
  call ISDestroy(is_tmp2,ierr)

  call VecScatterBegin(vec_scat,proc_id_loc,proc_id_nat, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scat,proc_id_loc,proc_id_nat, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scat,ierr)

  offset=0
  call MPI_Exscan(ugrid%ngmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  allocate(int_array(ugrid%ngmax))
  do ghosted_id=1,ugrid%ngmax
    int_array(ghosted_id)=ghosted_id-1+offset
  enddo
  call ISCreateGeneral(option%mycomm,ugrid%ngmax, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr)

  do ghosted_id=1,ugrid%ngmax
    int_array(ghosted_id)=grid%nG2A(ghosted_id)-1
  enddo
  call ISCreateGeneral(option%mycomm,ugrid%ngmax, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr)
  deallocate(int_array)

  call VecScatterCreate(proc_id_nat,is_tmp1,proc_id_ghosted,is_tmp2,vec_scat,ierr)
  call ISDestroy(is_tmp1,ierr)
  call ISDestroy(is_tmp2,ierr)

  call VecScatterBegin(vec_scat,proc_id_nat,proc_id_ghosted, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scat,proc_id_nat,proc_id_ghosted, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scat,ierr)

  ! Compute number of the two types of internal faces:
  ! - nfaces_intrn_loc: Both cells sharing the face are local cells
  ! - nfaces_intrn_nonloc: One of the cells sharing the face is a ghost cell
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  call VecGetArrayF90(proc_id_ghosted,vec_ptr,ierr)
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      if (ghosted_id_up<=ugrid%nlmax .and. ghosted_id_dn<=ugrid%nlmax) then
        ! Both cells sharing the face are local cells
        nfaces_intrn_loc = nfaces_intrn_loc + 1
        cur_connection_set%local(iconn) = 1
      else
        if((ghosted_id_dn>ugrid%nlmax.and.vec_ptr(ghosted_id_dn)>option%myrank).or. &
           (ghosted_id_up>ugrid%nlmax.and.vec_ptr(ghosted_id_up)>option%myrank)) then
          ! If the ghost cell resides on a proc whose rank is greater than
          ! mine, the face is local
          nfaces_intrn_loc = nfaces_intrn_loc + 1
          cur_connection_set%local(iconn) = 1
        else
          ! Otherwise, the face is non-local
          nfaces_intrn_nonloc = nfaces_intrn_nonloc + 1
        endif
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  call VecRestoreArrayF90(proc_id_ghosted,vec_ptr,ierr)

  ! Save information about number of faces
  grid%nlmax_faces = nfaces_intrn_loc+nfaces_bnd
  grid%ngmax_faces = nfaces_intrn_loc+nfaces_bnd+nfaces_intrn_nonloc

  ! If there are any boundary faces, set up Boundary connections for the grid
  ! Note: Equivalent subroutine for STRUCTURED_GRID_MIMETIC is 
  ! StructGridComputeBoundConnect()
  if (nfaces_bnd>0) then
    ! allocate memory
    nullify(bnd_connections)
    bnd_connections => ConnectionCreate(nfaces_bnd,BOUNDARY_CONNECTION_TYPE)

    ! Populate 'bnd_connections'
    nconn = 0
    do local_id = 1, ugrid%nlmax
      nfaces = UCellGetNFaces(ugrid%cell_type(local_id),option)
      do iface = 1, nfaces
        face_id = ugrid%cell_to_face_ghosted(iface,local_id)
        if (ugrid%face_to_cell_ghosted(2,face_id) < 1) then
          ! boundary face, since not connected to 2 cells
          nconn=nconn+1
          bnd_connections%id_dn(nconn) = local_id
          call GridPopulateConnection(grid,bnd_connections,iface,nconn, &
                                    local_id,option)
          bnd_connections%cntr(:,nconn)=bnd_connections%intercp(:,nconn)
        endif
      enddo
    enddo

    ! Add to list
    allocate(grid%boundary_connection_set_list)
    call ConnectionInitList(grid%boundary_connection_set_list)
    call ConnectionAddToList(bnd_connections,grid%boundary_connection_set_list)
  endif

  ! Allocate memory to save 'faces'
  allocate(faces(grid%ngmax_faces))

  ! Add all internal-faces
  face_id = 0
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      face_id = face_id + 1
      faces(face_id)%conn_set_ptr => cur_connection_set
      faces(face_id)%id = iconn
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Add any boundary-faces
  connection_set_list => grid%boundary_connection_set_list
  cur_connection_set => connection_set_list%first
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      face_id = face_id + 1
      faces(face_id)%conn_set_ptr => cur_connection_set
      faces(face_id)%id = iconn
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Save it
  grid%faces => faces

  ! Free up memory
  call VecDestroy(proc_id_loc,ierr)
  call VecDestroy(proc_id_nat,ierr)
  call VecDestroy(proc_id_ghosted,ierr)
#endif

end subroutine GridPopulateFacesForUGrid

! ************************************************************************** !
!> This routine sets up cell to face connection for an unstructure grid for
!! MIMETIC discretization.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/29/13
! ************************************************************************** !
subroutine GridComputeCell2FaceForUGrid(grid,MFD,option)

  use MFD_Aux_module
  use Option_module
  use Unstructured_Cell_module

  implicit none

  type(grid_type) :: grid
  type(mfd_type), pointer :: MFD
  type(option_type) :: option

#ifdef DASVYAT

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  type(unstructured_grid_type),pointer :: ugrid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: icell
  PetscInt :: iface
  PetscInt :: local_id
  PetscInt :: local_id_up
  PetscInt :: local_id_dn
  PetscInt :: ghosted_id
  PetscInt :: ghosted_id_up
  PetscInt :: ghosted_id_dn
  PetscInt :: face_id
  PetscInt :: nfaces
  PetscInt :: nfaces_intrn_loc
  PetscInt :: nfaces_intrn_nonloc
  PetscInt :: nfaces_bnd
  PetscInt :: offset
  PetscInt :: iface2
  PetscInt :: face_id2
  PetscInt :: count

  PetscInt,pointer :: bnd_count(:)
  
  PetscBool :: found
  PetscErrorCode :: ierr
  
  character(len=MAXWORDLENGTH) :: filename

  MFD => MFDAuxCreate()
  grid%MFD => MFD
  ugrid => grid%unstructured_grid
 
  ! Test
  allocate(grid%fU2M(MAX_FACE_PER_CELL,grid%nlmax))
  allocate(grid%fM2U(grid%ngmax_faces))
  allocate(bnd_count(grid%nlmax))
  bnd_count=0
  grid%fM2U = -1
  grid%fU2M = -1

  ! Find offset for boundary faces
  offset=0
  nfaces_bnd=0
  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  do
    if (.not.associated(cur_connection_set)) exit
    offset=offset+cur_connection_set%num_connections
    cur_connection_set => cur_connection_set%next
  enddo

  call MFDAuxInit(MFD,grid%nlmax,option)

  ! For each local cell, allocate memory for mfd_auxvar_type that depends on
  ! on number of faces for a given cell.
  do local_id = 1, grid%nlmax
    aux_var => MFD%aux_vars(local_id)
    nfaces = UCellGetNFaces(ugrid%cell_type(local_id),option)
    call MFDAuxVarInit(aux_var,nfaces,option)
  enddo

  ! Compute fG2L and fL2G mapping
  local_id = 1
  do iface = 1, grid%ngmax_faces
    conn => grid%faces(iface)%conn_set_ptr
    face_id = grid%faces(iface)%id

    if (conn%itype==INTERNAL_CONNECTION_TYPE) then

      ghosted_id_up = conn%id_up(face_id)
      ghosted_id_dn = conn%id_dn(face_id)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      if (conn%local(face_id) == 1) then
        grid%fG2L(iface)=local_id
        grid%fL2G(local_id) = iface
        local_id = local_id + 1
      endif

      if (local_id_dn>0) then
        aux_var => MFD%aux_vars(local_id_dn)
        call MFDAuxAddFace(aux_var,option,iface)
      endif

      if (local_id_up>0) then
        aux_var => MFD%aux_vars(local_id_up)
        call MFDAuxAddFace(aux_var,option,iface)
      endif

      ! For 'local_id_up', find the face-id that is shared by cells
      ! local_id_up ----- local_id_dn
      found=PETSC_FALSE
      do iface2=1,MAX_FACE_PER_CELL
        face_id2=ugrid%cell_to_face_ghosted(iface2,local_id_up)
        if( (ugrid%face_to_cell_ghosted(1,face_id2)==ghosted_id_up.and. &
             ugrid%face_to_cell_ghosted(2,face_id2)==ghosted_id_dn).or. &
            (ugrid%face_to_cell_ghosted(1,face_id2)==ghosted_id_dn.and. &
             ugrid%face_to_cell_ghosted(2,face_id2)==ghosted_id_up)) then
          found=PETSC_TRUE
          grid%fU2M(iface2,local_id_up)=iface
          grid%fM2U(iface)=face_id2
          exit
        endif
      enddo
      if(.not.found) then
        option%io_buffer='1) UGRID face not found to match face in MIMETIC discretization'
        call printErrMsg(option)
      endif
      
      ! For 'local_id_dn', find the face-id that is shared by cells
      ! local_id_up ----- local_id_dn
      if(ghosted_id_dn<=grid%nlmax) then
        found=PETSC_FALSE
        do iface2=1,MAX_FACE_PER_CELL
          face_id2=ugrid%cell_to_face_ghosted(iface2,local_id_dn)
          if( (ugrid%face_to_cell_ghosted(1,face_id2)==ghosted_id_up.and. &
               ugrid%face_to_cell_ghosted(2,face_id2)==ghosted_id_dn).or. &
              (ugrid%face_to_cell_ghosted(1,face_id2)==ghosted_id_dn.and. &
               ugrid%face_to_cell_ghosted(2,face_id2)==ghosted_id_up)) then
            found=PETSC_TRUE
            grid%fU2M(iface2,local_id_dn)=iface
            grid%fM2U(iface)=face_id2
            exit
          endif
        enddo
        if(.not.found) then
          option%io_buffer='2) UGRID face not found to match face in MIMETIC discretization'
          call printErrMsg(option)
        endif
      endif

    else if (conn%itype==BOUNDARY_CONNECTION_TYPE) then

      ghosted_id_dn = conn%id_dn(face_id)
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
   
      if (local_id_dn>0) then
        aux_var => MFD%aux_vars(local_id_dn)
        call MFDAuxAddFace(aux_var,option,iface)
        grid%fG2L(iface)=local_id
        grid%fL2G(local_id) = iface
        local_id = local_id + 1
      endif

      nfaces=UCellGetNFaces(ugrid%cell_type(local_id_dn),option)
      count=0
      found=PETSC_FALSE
      do iface2=1,nfaces
        face_id2=ugrid%cell_to_face_ghosted(iface2,local_id_dn)
        if (ugrid%face_to_cell_ghosted(2,face_id2) < 1) then
          ! boundary face, since not connected to 2 cells
          count=count+1
          if(count>bnd_count(local_id_dn)) then
            bnd_count(local_id_dn)=bnd_count(local_id_dn)+1
            found=PETSC_TRUE
            grid%fM2U(iface)=face_id2
            grid%fU2M(iface2,local_id_dn)=iface
            exit
          endif
        endif
        if(found) exit
      enddo
      if(.not.found) then
        option%io_buffer='UGRID boundary face not found to match face in MIMETIC discretization'
        call printErrMsg(option)
      endif

    endif

  enddo

#endif

end subroutine GridComputeCell2FaceForUGrid

! ************************************************************************** !
!> This routine sets up cell to face connection for an unstructure grid for
!! MIMETIC discretization.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/29/13
! ************************************************************************** !
subroutine GridSetGlobalCell2FaceForUGrid(grid,MFD,DOF,option)

  use MFD_Aux_module
  use Option_module
  use Unstructured_Cell_module
 
  implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"
#include "finclude/petscsys.h"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"

  type(grid_type) :: grid
  type(mfd_type), pointer :: MFD
  PetscInt :: DOF
  type(option_type) :: option

  type(unstructured_grid_type),pointer :: ugrid
  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn

  PetscInt :: local_id
  PetscInt :: local_id_up
  PetscInt :: local_id_dn
  PetscInt :: ghosted_id
  PetscInt :: ghosted_id_up
  PetscInt :: ghosted_id_dn
  PetscInt :: iface
  PetscInt :: face_id
  PetscInt :: ghost_face_id
  PetscInt :: local_face_id
  PetscInt :: nfaces
  PetscInt :: ndof
  PetscInt :: global_offset
  PetscInt :: num_ghosted_upd
  PetscInt :: icount,icell,stride,global_neigh_id,jcount
  IS :: is_from
  IS :: is_to
  PetscViewer :: viewer
  VecScatter :: scatter

  PetscInt, pointer :: num_faces_cumm(:)

  PetscInt, allocatable :: ghosted_ids(:)
  PetscInt, allocatable :: indices_to(:)
  PetscInt, allocatable :: indices_from(:)
  PetscInt, allocatable :: int_array(:)
  PetscInt, pointer :: int_ptr(:)

  PetscScalar, pointer :: e2f_local_values(:)
  PetscScalar, pointer :: e2n_local_values(:)
  PetscScalar, pointer :: e2f_ghosted_values(:)
  PetscScalar, pointer :: e2n_ghosted_values(:)
  PetscScalar, pointer :: vec_ptr_e2f(:)
  PetscScalar, pointer :: vec_ptr_e2n(:)
  PetscScalar, pointer :: vec_ptr_e2f_gh(:)
  PetscScalar, pointer :: vec_ptr_e2n_gh(:)

  Vec :: e2f
  Vec :: e2n
  Vec :: e2f_ghosted
  Vec :: e2n_ghosted

  PetscErrorCode :: ierr

  select case(DOF)
    case(ONEDOF)
      ndof = 1
    case(NFLOWDOF)
      ndof = option%nflowdof
    case(NTRANDOF)
      ndof = option%ntrandof
  end select

  MFD%ndof = ndof
  ugrid => grid%unstructured_grid

  ! Allocate memory
  allocate(grid%fL2P(grid%nlmax_faces))
  allocate(grid%fG2P(grid%ngmax_faces))

  ! Set fL2P
  global_offset = 0
  do iface=1,grid%nlmax_faces
    grid%fL2P(iface)=grid%global_faces_offset+grid%global_cell_offset+iface-1
  enddo

  global_offset = grid%global_faces_offset + grid%global_cell_offset

  ! Create e2f and e2n
  stride=MAX_FACE_PER_CELL
  call VecCreate(option%mycomm,e2f,ierr)
  call VecSetSizes(e2f,grid%nlmax*stride,PETSC_DECIDE,ierr)
  call VecSetFromOptions(e2f,ierr)
  call VecDuplicate(e2f,e2n,ierr)

  call VecGetArrayF90(e2f,e2f_local_values,ierr)
  call VecGetArrayF90(e2n,e2n_local_values,ierr)

  e2f_local_values = 0
  e2n_local_values = 0

  do local_id = 1, grid%nlmax
    aux_var => MFD%aux_vars(local_id)
    do iface = 1, aux_var%numfaces
      ghost_face_id = aux_var%face_id_gh(iface)
      local_face_id = grid%fG2L(ghost_face_id)
      conn => grid%faces(ghost_face_id)%conn_set_ptr
      face_id = grid%faces(ghost_face_id)%id

      if (conn%itype==INTERNAL_CONNECTION_TYPE) then

        if(local_face_id>0) then
          ! Face is a 'local' internal connection
          e2f_local_values((local_id-1)*stride+iface)=grid%fL2P(local_face_id)+1

          ghosted_id_up=conn%id_up(face_id)
          ghosted_id_dn=conn%id_dn(face_id)
          local_id_up=grid%nG2L(ghosted_id_up)
          local_id_dn=grid%nG2L(ghosted_id_dn)

          if (local_id_up==local_id) then
            e2n_local_values((local_id-1)*stride+iface)=grid%nG2P(ghosted_id_dn)+1
          else if (local_id_dn==local_id) then
            e2n_local_values((local_id-1)*stride+iface)=grid%nG2P(ghosted_id_up)+1
          endif

        else if(local_face_id==0) then
          ! Face is a 'non-local' internal connection
          ghosted_id_up=conn%id_up(face_id)
          ghosted_id_dn=conn%id_dn(face_id)
          local_id_up=grid%nG2L(ghosted_id_up)
          local_id_dn=grid%nG2L(ghosted_id_dn)

          if (local_id_up==local_id) then
            e2n_local_values((local_id-1)*stride+iface)=grid%nG2P(ghosted_id_dn)+1
          else if (local_id_dn==local_id) then
            e2n_local_values((local_id-1)*stride+iface)=grid%nG2P(ghosted_id_up)+1
          endif
        endif

      else if (conn%itype == BOUNDARY_CONNECTION_TYPE) then

        if (local_face_id>0) then
          e2f_local_values((local_id-1)*stride+iface)=grid%fL2P(local_face_id)+1
        endif

        e2n_local_values((local_id-1)*stride+iface)=0
      endif
    enddo
  enddo

  call VecRestoreArrayF90(e2f,e2f_local_values,ierr)
  call VecRestoreArrayF90(e2n,e2n_local_values,ierr)

  ! Find PETSc ID of cells with which a non-local internal connections is shared
  allocate(ghosted_ids(grid%ngmax_faces-grid%nlmax_faces))
  ghosted_ids = 0

  num_ghosted_upd = 0
  do iface = 1, grid%ngmax_faces
    conn => grid%faces(iface)%conn_set_ptr
    face_id=grid%faces(iface)%id

    if (conn%itype==INTERNAL_CONNECTION_TYPE) then
      if (conn%local(face_id) == 0) then
        ! Non-local internal connection
        num_ghosted_upd = num_ghosted_upd + 1
        ghosted_id_up = conn%id_up(face_id)
        ghosted_id_dn = conn%id_dn(face_id)

        ! Check if upwind or downwind is a ghost cell
        if (grid%nG2L(ghosted_id_up)==0) then
          ghosted_ids(num_ghosted_upd)=grid%nG2P(ghosted_id_up)+1 ! +1 since global indexes strats from 0
        else if (grid%nG2L(ghosted_id_dn)==0) then
          ghosted_ids(num_ghosted_upd)=grid%nG2P(ghosted_id_dn)+1
        endif
      endif
    endif
  enddo

  ! Create sequential vectors with a stride
  call VecCreateSeq(PETSC_COMM_SELF,num_ghosted_upd*stride,e2f_ghosted,ierr)
  call VecCreateSeq(PETSC_COMM_SELF,num_ghosted_upd*stride,e2n_ghosted,ierr)
         
  allocate(indices_to(num_ghosted_upd))
  allocate(indices_from(num_ghosted_upd))

  do icount=1,num_ghosted_upd
    indices_to(icount)=icount-1
    indices_from(icount)=ghosted_ids(icount)-1
  enddo

  ! Create IS with a block size of 'stride'
  call ISCreateBlock(option%mycomm,stride,num_ghosted_upd,indices_to, &
                    PETSC_COPY_VALUES,is_to,ierr)
  call ISCreateBlock(option%mycomm,stride,num_ghosted_upd,indices_from,&
                    PETSC_COPY_VALUES,is_from,ierr)
  deallocate(indices_from)
  deallocate(indices_to)

  ! Create scatter context
  call VecScatterCreate(e2f,is_from,e2f_ghosted,is_to,scatter,ierr)
  call ISDestroy(is_to,ierr)
  call ISDestroy(is_from,ierr)

  ! Scatter forward the e2f and e2n
  call VecScatterBegin(scatter,e2f,e2f_ghosted,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter,e2f,e2f_ghosted,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterBegin(scatter,e2n,e2n_ghosted,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter,e2n,e2n_ghosted,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(scatter,ierr)

  call VecGetArrayF90(e2n_ghosted,vec_ptr_e2n_gh,ierr)
  call VecGetArrayF90(e2f_ghosted,vec_ptr_e2f_gh,ierr)
  call VecGetArrayF90(e2n,vec_ptr_e2n,ierr)
  call VecGetArrayF90(e2f,vec_ptr_e2f,ierr)
 
  jcount=0
  do local_id=1,grid%nlmax
    aux_var => MFD%aux_vars(local_id)
    do iface=1,aux_var%numfaces
      ghost_face_id=aux_var%face_id_gh(iface)
      local_face_id=grid%fG2L(ghost_face_id)

      if (local_face_id==0) then
        conn => grid%faces(ghost_face_id)%conn_set_ptr
        face_id=grid%faces(ghost_face_id)%id

        ghosted_id_up=conn%id_up(face_id)
        ghosted_id_dn=conn%id_dn(face_id)

        local_id_up=grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
        local_id_dn=grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

        if (local_id==local_id_up) then
          global_neigh_id=grid%nG2P(ghosted_id_dn)+1
        else if (local_id==local_id_dn) then
          global_neigh_id=grid%nG2P(ghosted_id_up)+1
        endif

        ! GB: Can avoid the this loop
        !do jcount=1,num_ghosted_upd
        !  if (ghosted_ids(jcount)==global_neigh_id) exit
        !enddo
        jcount=jcount+1

        do face_id=1,stride
          if (vec_ptr_e2n_gh((jcount-1)*stride+face_id)==grid%nG2P(grid%nL2G(local_id))+1) then
            vec_ptr_e2f((local_id-1)*stride+iface)=vec_ptr_e2f_gh((jcount-1)*stride+face_id)
            grid%fG2P(ghost_face_id)=int(vec_ptr_e2f_gh((jcount-1)*stride+face_id))-1
            exit
          endif
        end do
      else
        grid%fG2P(ghost_face_id)=grid%fL2P(local_face_id)
      endif
    enddo
  enddo

  call VecRestoreArrayF90(e2n_ghosted,vec_ptr_e2n_gh,ierr)
  call VecRestoreArrayF90(e2f_ghosted,vec_ptr_e2f_gh,ierr)
  call VecDestroy(e2f_ghosted,ierr)
  call VecDestroy(e2n_ghosted,ierr)

  ! Count number of faces for all local cells
  nfaces=0
  allocate(num_faces_cumm(ugrid%nlmax))
  num_faces_cumm=0
  do local_id=1,ugrid%nlmax
    nfaces=nfaces+UCellGetNFaces(ugrid%cell_type(local_id),option)
    if(local_id<ugrid%nlmax) then
      num_faces_cumm(local_id+1)=num_faces_cumm(local_id)+nfaces
    endif
  enddo

  ! Allocate memory for grid%e2f and grid%e2n
  call VecCreate(option%mycomm,grid%e2f,ierr)
  call VecSetSizes(grid%e2f,nfaces,PETSC_DECIDE,ierr)
  call VecSetFromOptions(grid%e2f,ierr)
  call VecDuplicate(grid%e2f,grid%e2n,ierr)

  call VecGetArrayF90(grid%e2f,e2f_local_values,ierr)
  call VecGetArrayF90(grid%e2n,e2n_local_values,ierr)

  nfaces=0
  do local_id=1,grid%nlmax
    do iface=1,UCellGetNFaces(ugrid%cell_type(local_id),option)
      e2f_local_values(nfaces+iface)=vec_ptr_e2f((local_id-1)*stride+iface)
      e2n_local_values(nfaces+iface)=vec_ptr_e2n((local_id-1)*stride+iface)
    enddo
    nfaces=nfaces+UCellGetNFaces(ugrid%cell_type(local_id),option)
  enddo

  call VecRestoreArrayF90(grid%e2f,e2f_local_values,ierr)
  call VecRestoreArrayF90(grid%e2n,e2n_local_values,ierr)
  call VecRestoreArrayF90(e2n,vec_ptr_e2n,ierr)
  call VecRestoreArrayF90(e2f,vec_ptr_e2f,ierr)

  call VecDestroy(e2f,ierr)
  call VecDestroy(e2n,ierr)

  call CreateMFDStruct4LP(grid,MFD,ndof,option)

  deallocate(num_faces_cumm)


end subroutine GridSetGlobalCell2FaceForUGrid

end module Grid_module
