module Grid_module

  use Structured_Grid_module
  use Unstructured_Grid_module
  use Connection_module
  use MFD_Aux_module
 
  implicit none

  private
 
#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type, public :: grid_type 
  
    character(len=MAXWORDLENGTH) :: ctype
    PetscInt :: itype  ! type of grid (e.g. structured, unstructured, etc.)
    
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
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

    !geh: do we even need nL2A???
    
    ! nL2A :   collective, local => natural index, used for initialization   
    !                               and source/sink setup  (zero-based)
    PetscInt, pointer :: nL2G(:), nG2L(:), nL2A(:)
    PetscInt, pointer :: nG2A(:), nG2P(:), nG2LP(:)

    PetscInt, pointer :: fL2G(:), fG2L(:), fG2P(:), fL2P(:)
    PetscInt, pointer :: fL2B(:)
    Vec :: e2f             ! global vector to establish connection between global face_id and cell_id
    Vec :: e2n, e2n_LP     ! global cell connectivity vector

    
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
#ifdef SURFACE_FLOW
    !type(connection_set_list_type), pointer :: surface_internal_connection_set_list
#endif
    type(face_type), pointer :: faces(:)
    type(mfd_type), pointer :: MFD

#ifdef SUBCONTINUUM_MODEL
    ! Save no. of subcontinuum subgrids for all subcontinua patched 
    ! together in one array
    PetscInt, pointer :: subcontinuum_grid(:)
    ! Offsets to access subcontinuum_grid array. No. of rows = no. of cells
    ! in the patch. First column = no. of subcontinua at the cell
    ! Second column = offset to access subcontinuum_grid for the cell
    PetscInt, pointer :: subcontinuum_grid_offset(:,:)
#endif

  end type grid_type

  type, public :: face_type
    type(connection_set_type), pointer :: conn_set_ptr
    PetscInt :: id
  end type face_type

  interface GridVecGetArrayF90
     module procedure GridVecGetArrayCellF90
     module procedure GridVecGetArraySideF90
  end interface

  public :: GridCreate, &
            GridDestroy, &
            GridComputeInternalConnect, &
            GridMapIndices, &
            GridComputeSpacing, &
            GridComputeCoordinates, &
            GridComputeVolumes, &
            GridLocalizeRegions, &
            GridLocalizeRegionsForUGrid, &
            GridPopulateConnection, &
            GridPopulateFaces, &
            GridCopyIntegerArrayToVec, &
            GridCopyRealArrayToVec, &
            GridCopyVecToIntegerArray, &
            GridCopyVecToRealArray, &
            GridCreateNaturalToGhostedHash, &
            GridDestroyHashTable, &
            GridGetLocalGhostedIdFromHash, &
            GridVecGetMaskArrayCellF90, &
            GridVecGetArrayF90, &
            GridVecRestoreArrayF90, &
            GridIndexToCellID, &
            GridComputeCell2FaceConnectivity, &
            GridComputeGlobalCell2FaceConnectivity, &
            GridGetGhostedNeighbors
contains

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

#ifdef SURFACE_FLOW
  !nullify(grid%surface_internal_connection_set_list)
#endif

  nullify(grid%nL2G)
  nullify(grid%nG2L)
  nullify(grid%nL2A)
  nullify(grid%nG2A)
  nullify(grid%nG2P)

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

#ifdef DASVYAT  
  nullify(grid%faces)
  nullify(grid%MFD)
#endif

  nullify(grid%hash)
  grid%num_hash_bins = 1000

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
    case(UNSTRUCTURED_GRID) 
      connection_set => &
        UGridComputeInternConnect(grid%unstructured_grid,grid%x,grid%y, &
                                  grid%z,ugdm%scatter_ltol,option)
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
    case(UNSTRUCTURED_GRID) 
!      connection_bound_set => &
!        UGridComputeBoundConnect(grid%unstructured_grid,option)
  end select

#ifdef DASVYAT  
  if (associated(connection_bound_set)) then
    allocate(grid%boundary_connection_set_list)
    call ConnectionInitList(grid%boundary_connection_set_list)
    call ConnectionAddToList(connection_bound_set,grid%boundary_connection_set_list)
  end if
  if ((grid%itype==STRUCTURED_GRID_MIMETIC)) then
    allocate(grid%fL2G(grid%nlmax_faces))
    allocate(grid%fG2L(grid%ngmax_faces))
  
    grid%fL2G = 0
    grid%fG2L = 0
    call GridPopulateFaces(grid, option)
  end if
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
                                        iface,iconn,cell_id_ghosted)
    case(UNSTRUCTURED_GRID)
      call UGridPopulateConnection(grid%unstructured_grid,connection,iface,&
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

subroutine GridComputeCell2FaceConnectivity(grid, MFD_aux, option)


  use MFD_Aux_module
  use Option_module

  implicit none

  type(grid_type) :: grid
  type(mfd_type), pointer :: MFD_aux
  

!  type(auxilliary_type) :: aux
  type(option_type) :: option

#ifdef DASVYAT

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  PetscInt :: icount, icell, iface, local_id
  PetscInt :: local_id_dn, local_id_up, ghosted_id_dn, ghosted_id_up
  character(len=MAXWORDLENGTH) :: filename


  PetscInt, pointer :: numfaces(:)





  MFD_aux => MFDAuxCreate()
  grid%MFD => MFD_aux
 
  call MFDAuxInit(MFD_aux, grid%nlmax, option)
  allocate(numfaces(grid%nlmax))

  numfaces = 6
  
!  do icount = 1, grid%ngmax_faces
!    conn => grid%faces(icount)%conn_set_ptr
!    iface = grid%faces(icount)%id
!    if (conn%itype==BOUNDARY_CONNECTION_TYPE) then
!        ghosted_id_dn = conn%id_dn(iface)
!        local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
!
!        if (local_id_dn>0) then
!           aux_var => MFD_aux%aux_vars(local_id_dn)
!           numfaces(local_id_dn) = numfaces(local_id_dn) + 1
!        end if
!        
!    else if (conn%itype==INTERNAL_CONNECTION_TYPE) then 
!        ghosted_id_up = conn%id_up(iface)
!        ghosted_id_dn = conn%id_dn(iface)
!
!        local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
!        local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
!
!        if (local_id_dn>0) then
!           aux_var => MFD_aux%aux_vars(local_id_dn)
!           numfaces(local_id_dn) = numfaces(local_id_dn) + 1
!        end if
!        if (local_id_up>0) then
!           aux_var => MFD_aux%aux_vars(local_id_up)
!           numfaces(local_id_up) = numfaces(local_id_up) + 1
!        end if
!    end if
!  end do

  do icell = 1, grid%nlmax
    aux_var => MFD_aux%aux_vars(icell)
    call MFDAuxVarInit(aux_var, numfaces(icell), option)
  end do


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
        end if
        
    else if (conn%itype==INTERNAL_CONNECTION_TYPE) then 
        ghosted_id_up = conn%id_up(iface)
        ghosted_id_dn = conn%id_dn(iface)

        local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
        local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

        if (conn%local(iface) == 1) then
           grid%fG2L(icount)=local_id
           grid%fL2G(local_id) = icount
           local_id = local_id + 1
        end if
            


        if (local_id_dn>0) then
           aux_var => MFD_aux%aux_vars(local_id_dn)
           call MFDAuxAddFace(aux_var,option, icount)
        end if
        if (local_id_up>0) then
           aux_var => MFD_aux%aux_vars(local_id_up)
           call MFDAuxAddFace(aux_var,option, icount)
        end if
    end if
    
  end do

 
    
!  do icell = 1, grid%nlmax
!    aux_var => MFD_aux%aux_vars(icell)
!    write(9,*) option%myrank, icell
!    write(9,*) option%myrank, "+++++++++++++++", (aux_var%face_id_gh(icount),icount=1,6)
!  end do
  
  if (associated(numfaces)) deallocate(numfaces)

#endif

end subroutine GridComputeCell2FaceConnectivity


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


#include "definitions.h"

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
!    Vec :: vec_LP_cell_id
!    Vec :: vec_LP_cell_id_loc



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
!    PetscScalar, pointer :: lp_cell_ids(:), lp_cell_ids_loc(:)

    PetscScalar, pointer :: e2f_ghosted_values(:)
    PetscScalar, pointer :: e2n_ghosted_values(:)
   
    PetscScalar, pointer :: vec_ptr_e2f(:)
    PetscScalar, pointer :: vec_ptr_e2n(:)
    
    PetscScalar, pointer :: vec_ptr_e2f_gh(:)
    PetscScalar, pointer :: vec_ptr_e2n_gh(:)

    VecScatter :: VC_global2ghosted

    

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
!    grid%global_faces_offset = 0
!    grid%global_cell_offset = 0
!
!    call MPI_Exscan(grid%nlmax_faces, grid%global_faces_offset, &
!                      ONE_INTEGER,MPI_INTEGER,MPI_SUM,option%mycomm,ierr)
!
!    call MPI_Exscan(grid%nlmax, grid%global_cell_offset, &
!                      ONE_INTEGER,MPI_INTEGER,MPI_SUM,option%mycomm,ierr)


    do iface = 1,grid%nlmax_faces
       grid%fL2P(iface)=grid%global_faces_offset + grid%global_cell_offset + iface - 1
    end do

    global_offset = grid%global_faces_offset + grid%global_cell_offset

!    write(*,*) option%myrank, global_offset, grid%global_faces_offset, grid%global_cell_offset


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
              end if

          else if (local_face_id == 0) then

              ghosted_id_up = conn%id_up(iface)
              ghosted_id_dn = conn%id_dn(iface)

              local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
              local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
              if (local_id_up==icell) then
                 e2n_local_values((icell-1)*stride + icount) = grid%nG2P(ghosted_id_dn) + 1
              else if (local_id_dn==icell) then
                 e2n_local_values((icell-1)*stride + icount) = grid%nG2P(ghosted_id_up) + 1
              end if

         end if

       else if (conn%itype==BOUNDARY_CONNECTION_TYPE) then
            e2n_local_values((icell-1)*stride + icount) = 0
            if (local_face_id > 0) e2f_local_values((icell-1)*stride + icount) = grid%fL2P(local_face_id) + 1
        end if
      end do
   end do


   call VecRestoreArrayF90(grid%e2f, e2f_local_values, ierr)
   call VecRestoreArrayF90(grid%e2n, e2n_local_values, ierr)


!   call PetscViewerASCIIOpen(option%mycomm,'vec_e2f_before.out',viewer,ierr)
!   call VecView(grid%e2f,viewer,ierr)
!   call VecView(vec_LP_cell_id,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)

!   call PetscViewerASCIIOpen(option%mycomm,'vec_LP_cell_id.out',viewer,ierr)
!   call VecView(vec_LP_cell_id,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)

!   stop


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
           end if
         end if
       end if
     end do


        call VecCreateSeq(PETSC_COMM_SELF, num_ghosted_upd*stride, ghosted_e2f, ierr)
        call VecCreateSeq(PETSC_COMM_SELF, num_ghosted_upd*stride, ghosted_e2n, ierr)
         
        allocate(strided_indices_local(num_ghosted_upd))
        allocate(strided_indices_ghosted(num_ghosted_upd))

        do icount = 1, num_ghosted_upd
          strided_indices_local(icount) = (icount -1)
          strided_indices_ghosted(icount) = (ghosted_ids(icount)-1)
        end do


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
 
!    write(*,*) option%myrank, num_ghosted_upd
!    do icount = 1, num_ghosted_upd
!      write(*,*) option%myrank, ghosted_ids(icount)," f: ", (vec_ptr_e2f_gh((icount - 1)*stride + icell), icell =1,6) 
!      write(*,*) option%myrank, ghosted_ids(icount)," n: ", (vec_ptr_e2n_gh((icount - 1)*stride + icell), icell =1,6) 
!    end do


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
            end if
            do jcount = 1, num_ghosted_upd
                if (ghosted_ids(jcount)==global_neigh_id) exit
            end do                    
            do iface=1, stride             ! assumption that cell has 6 faces
               if (vec_ptr_e2n_gh((jcount-1)*stride + iface) == grid%nG2P(grid%nL2G(icell)) + 1) then
                   vec_ptr_e2f((icell -1)*stride +icount) = vec_ptr_e2f_gh((jcount-1)*stride + iface)
                   grid%fG2P(ghost_face_id) = vec_ptr_e2f_gh((jcount-1)*stride + iface) - 1
                   exit
               end if
            end do
          else 
             grid%fG2P(ghost_face_id) = grid%fL2P(local_face_id)
          end if

       end do
     end do






!    do iface = 1, grid%ngmax_faces
!      write(*,*) option%myrank, iface, grid%fG2P(iface)
!    end do

!    do icell = 1, grid%ngmax
!      grid%nG2LP(icell) = lp_cell_ids_loc(icell) - 1
!      if (option%myrank == 0) write(*,*) option%myrank, icell, grid%nG2LP(icell)
!    end do



    call VecRestoreArrayF90(grid%e2n, vec_ptr_e2n, ierr)
    call VecRestoreArrayF90(grid%e2f, vec_ptr_e2f, ierr)

    call VecRestoreArrayF90(ghosted_e2n, vec_ptr_e2n_gh, ierr)
    call VecRestoreArrayF90(ghosted_e2f, vec_ptr_e2f_gh, ierr)

  
!    call PetscViewerASCIIOpen(option%mycomm,'vec_e2f_after.out',viewer,ierr)
!    call VecView(grid%e2f,viewer,ierr)
!    call PetscViewerDestroy(viewer,ierr)


!    call PetscViewerASCIIOpen(option%mycomm,'nG2LP.out',viewer,ierr)
!    call VecView(grid%nG2LP, viewer,ierr)
!    call PetscViewerDestroy(viewer,ierr)

    call VecDestroy(ghosted_e2n, ierr)
    call VecDestroy(ghosted_e2f, ierr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!    call CreateMFDStruct4Faces(grid, MFD_aux, ndof, option)
    call CreateMFDStruct4LP(grid, MFD_aux, ndof, option)

    deallocate(ghosted_ids)
    deallocate(strided_indices_local)
    deallocate(strided_indices_ghosted)
!    stop


#endif

end subroutine GridComputeGlobalCell2FaceConnectivity

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
!
!   write(*,*) 'istart', istart, 'ndof', ndof, 'grid%nlmax_faces', grid%nlmax_faces
!  do iface = 1, grid%nlmax_faces
!     write(*,*) iface, int_array(iface)
!  enddo



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

   

!   call PetscViewerASCIIOpen(option%mycomm,'is_local_petsc.out',viewer,ierr)
!   call ISView(MFD_aux%is_local_petsc_faces,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)


   call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_petsc.out',viewer,ierr)
   call ISView(MFD_aux%is_ghosted_petsc_faces,viewer,ierr)
   call PetscViewerDestroy(viewer,ierr)


!   call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_local.out',viewer,ierr)
!   call ISView(MFD_aux%is_ghosted_local_faces,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)


    call VecScatterCreate(local_vec,MFD_aux%is_local_local_faces,global_vec, &
                        MFD_aux%is_local_petsc_faces,MFD_aux%scatter_ltog_faces,ierr)


    call VecScatterCreate(global_vec,MFD_aux%is_ghosted_petsc_faces,local_vec, &
                        MFD_aux%is_ghosted_local_faces, MFD_aux%scatter_gtol_faces, ierr)

    
!   call PetscViewerASCIIOpen(option%mycomm,'scatter_gtol_faces.out',viewer,ierr)
!   call VecScatterView(MFD_aux%scatter_gtol_faces,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)

 ! Create local to local scatter.  Essentially remap the global to local as
 ! PETSc does in daltol.c
  call VecScatterCopy(MFD_aux%scatter_gtol_faces, MFD_aux%scatter_ltol_faces, ierr)
  call ISGetIndicesF90(MFD_aux%is_local_local_faces,int_ptr,ierr)
  call VecScatterRemap(MFD_aux%scatter_ltol_faces,int_ptr,PETSC_NULL_INTEGER,ierr)
  call ISRestoreIndicesF90(MFD_aux%is_local_local_faces,int_ptr,ierr)


  call VecDestroy(local_vec, ierr)
  call VecDestroy(global_vec, ierr)
  

end subroutine CreateMFDStruct4Faces

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

 !    Compute global number of faces
 !     call MPI_Allreduce(grid%nlmax_faces, grid%nmax_faces,ONE_INTEGER_MPI, &
 !                       MPI_INT,MPI_SUM,option%mycomm,ierr)


    NG = grid%ngmax_faces + grid%ngmax
    NL = grid%nlmax_faces + grid%nlmax


!     write(*,*) "Faces", option%myrank, "ngmax_faces", grid%ngmax_faces, "nlmax_faces", grid%nlmax_faces
!     write(*,*) "Cells", option%myrank, "ngmax", grid%ngmax, "nlmax", grid%nlmax
!     write(*,*) 'NG', NG, 'NL', NL

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

!   call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_petsc_LP.out',viewer,ierr)
!   call ISView(MFD_aux%is_ghosted_petsc_LP,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)


!   call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_petsc.out',viewer,ierr)
!   call ISView(MFD_aux%is_ghosted_petsc_faces,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)


!   call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_local.out',viewer,ierr)
!   call ISView(MFD_aux%is_ghosted_local_faces,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)


    call VecScatterCreate(local_vec_LP,MFD_aux%is_local_local_LP,global_vec_LP, &
                        MFD_aux%is_local_petsc_LP,MFD_aux%scatter_ltog_LP,ierr)


    call VecScatterCreate(global_vec_LP,MFD_aux%is_ghosted_petsc_LP,local_vec_LP, &
                        MFD_aux%is_ghosted_local_LP, MFD_aux%scatter_gtol_LP, ierr)


    
!   call PetscViewerASCIIOpen(option%mycomm,'scatter_gtol_LP.out',viewer,ierr)
!   call VecScatterView(MFD_aux%scatter_gtol_LP,viewer,ierr)
!   call PetscViewerDestroy(viewer,ierr)

 ! Create local to local scatter.  Essentially remap the global to local as
 ! PETSc does in daltol.c
  call VecScatterCopy(MFD_aux%scatter_gtol_LP, MFD_aux%scatter_ltol_LP, ierr)
  call ISGetIndicesF90(MFD_aux%is_local_local_LP,int_ptr,ierr)
  call VecScatterRemap(MFD_aux%scatter_ltol_LP,int_ptr,PETSC_NULL_INTEGER,ierr)
  call ISRestoreIndicesF90(MFD_aux%is_local_local_LP,int_ptr,ierr)


  call VecDestroy(local_vec_LP, ierr)
  call VecDestroy(global_vec_LP, ierr)

!   write(*,*) "End CreateMFDStruct4LP"
  
!   stop

end subroutine CreateMFDStruct4LP
   



! ************************************************************************** !
!
! GridMapIndices: maps global, local and natural indices of cells 
!                 to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridMapIndices(grid, sgdm)


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

  PetscInt :: ierr, icount
  PetscInt, allocatable :: int_tmp(:)
  PetscOffset :: i_da
  
  select case(grid%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call StructuredGridMapIndices(grid%structured_grid,grid%nG2L,grid%nL2G, &
                                    grid%nL2A,grid%nG2A)
#ifdef DASVYAT
      if ((grid%itype==STRUCTURED_GRID_MIMETIC)) then
         allocate(grid%nG2P(grid%ngmax))
         allocate(int_tmp(grid%ngmax))
         call DMDAGetGlobalIndices(sgdm,  grid%ngmax, int_tmp, i_da, ierr)
         do icount = 1, grid%ngmax
         !   write(*,*) icount, int_tmp(icount + i_da)
            grid%nG2P(icount) = int_tmp(icount + i_da)
         end do
        deallocate(int_tmp)
      end if
#endif
    case(UNSTRUCTURED_GRID)
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
      call StructuredGridComputeSpacing(grid%structured_grid,option)
    case(UNSTRUCTURED_GRID)
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
  
  implicit none

  type(grid_type) :: grid
  PetscReal :: origin_global(3)
  type(option_type) :: option
  type(ugdm_type), optional :: ugdm ! sp  
  
  PetscErrorCode :: ierr  
  
  select case(grid%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      allocate(grid%x(grid%ngmax))
      grid%x = 0.d0
      allocate(grid%y(grid%ngmax))
      grid%y = 0.d0
      allocate(grid%z(grid%ngmax))
      grid%z = 0.d0
      call StructuredGridComputeCoord(grid%structured_grid,option, &
                                      origin_global, &
                                      grid%x,grid%y,grid%z, &
                                      grid%x_min_local,grid%x_max_local, &
                                      grid%y_min_local,grid%y_max_local, &
                                      grid%z_min_local,grid%z_max_local)
    case(UNSTRUCTURED_GRID)
      allocate(grid%x(grid%ngmax))
      grid%x = 0.d0
      allocate(grid%y(grid%ngmax))
      grid%y = 0.d0
      allocate(grid%z(grid%ngmax))
      grid%z = 0.d0
      call UGridComputeCoord(grid%unstructured_grid,option, &
                             ugdm%scatter_ltol, & !sp 
                             grid%x,grid%y,grid%z, &
                             grid%x_min_local,grid%x_max_local, &
                             grid%y_min_local,grid%y_max_local, &
                             grid%z_min_local,grid%z_max_local)
  end select

  if (associated(grid%structured_grid)) then
    if (grid%structured_grid%p_samr_patch==0) then
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
 endif

  if (associated(grid%unstructured_grid)) then
    !if (grid%unstructured_grid%p_samr_patch==0) then
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
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(grid_type) :: grid
  type(option_type) :: option
  Vec :: volume
  
  select case(grid%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      call StructuredGridComputeVolumes(grid%x,grid%structured_grid,option, &
                                        grid%nL2G,volume)
    case(UNSTRUCTURED_GRID)
      call UGridComputeVolumes(grid%unstructured_grid,option,volume)
  end select

end subroutine GridComputeVolumes

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

        ! convert indexing from global (entire domain) to local processor
        region%i1 = region%i1 - grid%structured_grid%nxs
        region%i2 = region%i2 - grid%structured_grid%nxs
        region%j1 = region%j1 - grid%structured_grid%nys
        region%j2 = region%j2 - grid%structured_grid%nys
        region%k1 = region%k1 - grid%structured_grid%nzs
        region%k2 = region%k2 - grid%structured_grid%nzs
          
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
!          if (region%num_cells > 0) then
!            region%coordinates(1)%x = grid%x(region%cell_ids(ONE_INTEGER))
!            region%coordinates(1)%y = grid%y(region%cell_ids(ONE_INTEGER))
!            region%coordinates(1)%z = grid%z(region%cell_ids(ONE_INTEGER))
!          endif
        else
          region%num_cells = 0
        endif
     
        if (count /= region%num_cells) then
          option%io_buffer = 'Mismatch in number of cells in block region'
          call printErrMsg(option)
        endif

      else if (associated(region%coordinates)) then
      
        ! 1 or 2 coordinates
        if (size(region%coordinates) <= TWO_INTEGER) then
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
                  if (.not. (option%use_samr)) then
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
                  endif
                case(UNSTRUCTURED_GRID)
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
                case(UNSTRUCTURED_GRID)
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
                    call UGridGetCellsInRectangle(x_min,x_max,y_min,y_max, &
                                                  z_min,z_max, &
                                                  grid%unstructured_grid,option, &
                                                  region%num_cells,region%cell_ids, &
                                                  region%faces)
                  endif
              end select
            endif

            if (.not. (option%use_samr)) then
              call MPI_Allreduce(iflag,i,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX, &
                                 option%mycomm,ierr)
            else
              i = 0
            endif

            iflag = i
            if (iflag > 0) then
              option%io_buffer = 'GridLocalizeRegions, between two points'
              call printErrMsg(option)
            endif
          endif  
        else
          option%io_buffer = 'GridLocalizeRegions: more than 2 coordinates' // &
                             ' not supported in region object'
          call printErrMsg(option)
        endif  
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
        case(UNSTRUCTURED_GRID)
          allocate(temp_int_array(region%num_cells))
          temp_int_array = 0
          local_count=0
          do count=1,region%num_cells
            local_id = GridGetLocalIdFromNaturalId( grid, region%cell_ids(count))
            if (local_id <= 0) cycle
            if (local_id > grid%unstructured_grid%nlmax) cycle
            local_count = local_count + 1
            temp_int_array(local_count) = local_id
          enddo
          if (local_count /= region%num_cells) then
            deallocate(region%cell_ids)
            allocate(region%cell_ids(local_count))
            region%num_cells = local_count 
          endif
          region%cell_ids(1:local_count) = temp_int_array(1:local_count)
          deallocate(temp_int_array)
        case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
!sp following was commented out 
!sp remove? 
!geh: Not for now.  The below maps a list of natural ids to local.  This is now done
!     in hdf5.F90 for lists of cells from hdf5 files.  If we elect to support ascii
!     files, we will need this functionality here.
#if 0
          do count=1,region%num_cells
            i = mod(region%cell_ids(count),grid%structured_grid%nx) - &
                  grid%structured_grid%nxs
            j = mod((region%cell_ids(count)-1)/grid%structured_grid%nx, &
                    grid%structured_grid%ny)+1 - &
                  grid%structured_grid%nys
            k = ((region%cell_ids(count)-1)/grid%structured_grid%nxy)+1 - &
                  grid%structured_grid%nzs
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
        case default
          option%io_buffer = 'GridLocalizeRegions: define region by list ' // &
            'of cells not implemented: ' // trim(region%name)
          call printErrMsg(option)
      end select
      !sp end 
    else if (associated(region%vertex_ids)) then
      call GridLocalizeRegionsForUGrid(grid, region, option)
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
! GridLocalizeRegionsForUGrid: Resticts regions to cells local 
!    to processor for unstrucutred mesh
! author: Gautam Bisht
! date: 5/30/2011
!
! ************************************************************************** !
subroutine GridLocalizeRegionsForUGrid(grid, region, option)

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

      if (region%iface /= 0) then
        allocate(region%faces(region%num_cells))
          
        !
        ! In unstrucutred mesh, faces for hex are saved in the following order:
        !   (1) SOUTH
        !   (2) EAST
        !   (3) NORTH
        !   (4) WEST
        !   (5) BOTTOM
        !   (6) TOP
        ! PS: This will not work for WEDGE
        !
        select case (region%iface)
          case(WEST_FACE)
            region%faces = 4
          case(EAST_FACE)
            region%faces = 2
          case(SOUTH_FACE)
            region%faces = 1
          case(NORTH_FACE)
            region%faces = 3
          case(BOTTOM_FACE)
            region%faces = 5
          case(TOP_FACE)
            region%faces = 6
          case(NULL_FACE)
            region%faces = NULL_FACE
        end select
      endif
    else
      nullify(region%cell_ids)
    endif
    
    call VecRestoreArrayF90(vec_cell_ids_loc,v_loc_p,ierr)
      

  endif

  
  !
  !  Region is defined as a collection of faces. Each face is identified by 
  !  a list of vertices forming it.
  !
  if (associated(region%vertex_ids)) then
    !
    ! Create a sparse matrix: mat_vert2cell
    !
    !  - size(mat_vert2cell) = num_vertices_global x num_cells_global
    !  - i-th row of mat_vert2cell corresponds to global vertex id in natural
    !    indexing
    !  - j-th column of mat_vert2cell corresponds to global cell id in 
    !    natural indexing
    !  - mat_vert2cell(i,j) = j (i.e. i-th global vertex constitutes j-th
    !    cell
    !
    call MatCreateMPIAIJ(option%mycomm, PETSC_DECIDE, PETSC_DECIDE, &
                          ugrid%num_vertices_global, ugrid%nmax, &
                          MAX_CELLS_SHARING_A_VERTEX, PETSC_NULL_INTEGER, &
                          MAX_CELLS_SHARING_A_VERTEX, PETSC_NULL_INTEGER, &
                          mat_vert2cell, ierr)
    
    do ghosted_id = 1, ugrid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (local_id < 1) cycle
      natural_id = grid%nG2A(ghosted_id)
      do ii = 1, ugrid%cell_vertices(0, local_id)
        vertex_id = ugrid%cell_vertices(ii, local_id)
!geh: I believe that this is incorrect since MatSetValues uses petsc ordering,
!     unless teh matrix is a local MATSEQXXX matrix
        call MatSetValues(mat_vert2cell, &
                          1, &
                          ugrid%vertex_ids_natural(vertex_id)-1, &
                          1, &
                          natural_id-1, &
                          natural_id-1.0d0, &
                          INSERT_VALUES, &
                          ierr)
      enddo
    enddo
      
    call MatAssemblyBegin(mat_vert2cell, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(  mat_vert2cell, MAT_FINAL_ASSEMBLY, ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'mat_vert2cell.out', viewer, ierr)
    call MatView(mat_vert2cell, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif      

    call MatGetOwnershipRange(mat_vert2cell, rstart, rend, ierr)
    if (option%mycommsize > 1) then
      call MatMPIAIJGetSeqAIJ(mat_vert2cell, mat_vert2cell_diag, &
                              mat_vert2cell_offdiag, icol, iicol, ierr)
      call MatGetRowIJF90(mat_vert2cell_diag, 1, PETSC_FALSE, PETSC_FALSE, &
                          n, ia_p, ja_p, done, ierr)
      call MatGetArray(mat_vert2cell_diag, aa, aaa, ierr)
      ! call MatGetArrayF90(mat_vert2cell_diag, aa, ierr)
    else 
      call MatGetRowIJF90(mat_vert2cell, 1, PETSC_FALSE, PETSC_FALSE, n, &
                          ia_p, ja_p, done, ierr)
      call MatGetArray(mat_vert2cell, aa, aaa, ierr)
      !call MatGetArrayF90(mat_vert2cell, aa, ierr)
    endif
      
    !
    ! vert2cell_array 
    !  Let MCSV = MAX_CELLS_SHARING_A_VERTEX
    !
    !  - size(vert2cell_array) = No. rows of mat_vert2cell owned by proc x MCSV
    !  - vert2cell_array((i-1)*MCSV + 1: i*MCSV ) = Cell ids (in natural 
    !    index) which contain i-th vertex as one of its vertices
    !
    allocate(vert2cell_array(1: (rend - rstart)*MAX_CELLS_SHARING_A_VERTEX))
    vert2cell_array = -1
    do ii = 1, n
      count = ia_p(ii + 1) - ia_p(ii)
      do jj = ia_p(ii), ia_p(ii + 1) - 1
        found = PETSC_FALSE
        do kk = 1, MAX_CELLS_SHARING_A_VERTEX
          index = (ii - 1)*MAX_CELLS_SHARING_A_VERTEX + kk
          if (vert2cell_array(index) == -1) then
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (.not.found) then
          option%io_buffer = 'Increase the value of ' // &
            'MAX_CELLS_SHARING_A_VERTEX within the code.'
          call printErrMsg(option)
        endif
        vert2cell_array(index) = aa(aaa+ jj)
      enddo
    enddo
    if (option%mycommsize > 1) then
      call MatRestoreRowIJF90(mat_vert2cell_diag, 1, PETSC_FALSE, &
                              PETSC_FALSE, n, ia_p, ja_p, done, ierr)
      call MatRestoreArray(mat_vert2cell_diag, aa, aaa, ierr)
      !call MatRestoreArrayF90(mat_vert2cell_diag, aa, ierr)
    else
      call MatRestoreRowIJF90(mat_vert2cell, 1, PETSC_FALSE, PETSC_FALSE, n, &
                              ia_p, ja_p, done, ierr)
      call MatRestoreArray(mat_vert2cell, aa, aaa, ierr)
      ! call MatRestoreArrayF90(mat_vert2cell, aa, ierr)
    endif
      
    if (option%mycommsize > 1) then
      call MatGetRowIJF90(mat_vert2cell_offdiag, 1, PETSC_FALSE, &
                          PETSC_FALSE, n, ia_p, ja_p, done, ierr)
      call MatGetArray(mat_vert2cell_offdiag, aa, aaa, ierr)
      !call MatGetArrayF90(mat_vert2cell_offdiag,aa,ierr)
      do ii=1, n
        count = ia_p(ii+1) - ia_p(ii)
        do jj = ia_p(ii), ia_p(ii+1)-1
          found = PETSC_FALSE
          do kk = 1,MAX_CELLS_SHARING_A_VERTEX
            if (vert2cell_array( (ii-1)*MAX_CELLS_SHARING_A_VERTEX + kk) &
                                  == -1) then
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            option%io_buffer = 'Increase the value of ' // &
              'MAX_CELLS_SHARING_A_VERTEX within the code.'
            call printErrMsg(option)
          endif
          vert2cell_array( (ii-1)*MAX_CELLS_SHARING_A_VERTEX + kk ) = &
            aa(aaa+jj)
        enddo
      enddo
      call MatRestoreRowIJF90(mat_vert2cell_offdiag, 1, PETSC_FALSE, &
                              PETSC_FALSE, n, ia_p, ja_p, done, ierr)
      call MatGetArray(mat_vert2cell_offdiag,aa,aaa,ierr)
      !call MatGetArrayF90(mat_vert2cell_offdiag,aa,ierr)
      call MatDestroy(mat_vert2cell, ierr)
    endif
            
    !
    call VecCreateMPI(option%mycomm, (rend-rstart)*MAX_CELLS_SHARING_A_VERTEX, &
      PETSC_DECIDE, vec_vert2cell, ierr)
      
    call VecGetArrayF90(vec_vert2cell,v_loc_p,ierr)
    v_loc_p = vert2cell_array
    call VecRestoreArrayF90(vec_vert2cell,v_loc_p,ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_vert2cell.out', viewer, ierr)
    call VecView(vec_vert2cell, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    !
    ! vec_cell2facevert
    !   Contains vertices forming face. 
    !
    !  i-th cell j-th face vertex-0
    !  i-th cell j-th face vertex-1
    !  i-th cell j-th face vertex-2
    !  i-th cell j-th face vertex-3
    !  i-th cell j+1-th face vertex-0
    !  i-th cell j+1-th face vertex-1
    !  i-th cell j+1-th face vertex-2
    !  i-th cell j+1-th face vertex-3
    !  ...
    !  ...
    !  ...
    !  i+1-th cell j-th face vertex-0
    !  i+1-th cell j-th face vertex-1
    !  i+1-th cell j-th face vertex-2
    !  i+1-th cell j-th face vertex-3
    !
    call VecCreateMPI(option%mycomm, PETSC_DECIDE, &
                      ugrid%nmax*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE, &
                      vec_cell2facevert, ierr)

    ! Initialize vector
    call VecSet(vec_cell2facevert, -999.d0, ierr)
    call VecAssemblyBegin(vec_cell2facevert, ierr) ! vertex-id is 0-based
    call VecAssemblyEnd(  vec_cell2facevert, ierr) !
            
    allocate(tmp_int_array(ugrid%nlmax))
    tmp_int_array = 0
      
    do ii = 1, ugrid%max_ndual_per_cell*ugrid%ngmax
      ghosted_id = ugrid%face_to_cell_ghosted(1, ii)
      if (ghosted_id < 0 ) exit
      local_id   = grid%nG2L(ghosted_id)
      if (local_id < 1) cycle
      natural_id = grid%nG2A(ghosted_id) ! 1-based
      do jj = 1, MAX_VERT_PER_FACE
!geh        if ( ugrid%face_to_vertex_nindex(jj, ii) > 0 ) then
!gb     vertex_id_natural = ugrid%vertex_ids_natural(ugrid%face_to_vertex(jj,ii))
        if (ugrid%face_to_vertex_natural(jj,ii) > 0) then
          call VecSetValues(vec_cell2facevert,1, &
                        (natural_id - 1)*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE + &
                        tmp_int_array(local_id)*MAX_VERT_PER_FACE + jj - 1, &
                        ugrid%face_to_vertex_natural(jj, ii) - 1.d0, &
                        INSERT_VALUES, ierr)
        endif
      enddo 
      tmp_int_array(local_id) = tmp_int_array(local_id) + 1
    enddo
      
    call VecAssemblyBegin(vec_cell2facevert, ierr) ! vertex-id is 0-based
    call VecAssemblyEnd(  vec_cell2facevert, ierr) !
    deallocate(tmp_int_array)
#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_cell2facevert.out', viewer, ierr)
    call VecView(vec_cell2facevert, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
      
    allocate(tmp_int_array(region%num_verts*MAX_CELLS_SHARING_A_VERTEX))
    count = 0
    do ii = 1, region%num_verts
      do jj = 1, MAX_CELLS_SHARING_A_VERTEX
        count = count + 1
        tmp_int_array(count) = region%vertex_ids(1,ii)*MAX_CELLS_SHARING_A_VERTEX + jj
      enddo
    enddo

    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, &
                        region%num_verts*MAX_CELLS_SHARING_A_VERTEX, &
                        tmp_int_array, PETSC_COPY_VALUES, is_from, ierr)
    deallocate(tmp_int_array)
      
    call VecCreateMPI(option%mycomm, &
                region%num_verts*MAX_CELLS_SHARING_A_VERTEX, PETSC_DECIDE, &
                vec_vert2cell_reg_subset, ierr)
    call VecGetOwnershipRange(vec_vert2cell_reg_subset, istart, iend, ierr)

    allocate(tmp_int_array(region%num_verts*MAX_CELLS_SHARING_A_VERTEX))
    do ii = 1, region%num_verts*MAX_CELLS_SHARING_A_VERTEX
      tmp_int_array(ii) = ii + istart
    enddo
      
    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, &
                        region%num_verts*MAX_CELLS_SHARING_A_VERTEX, &
                        tmp_int_array, PETSC_COPY_VALUES, is_to, ierr)
    deallocate(tmp_int_array)
      
    call VecScatterCreate(vec_vert2cell, is_from, vec_vert2cell_reg_subset, &
                          is_to, vec_scat, ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)
    
    call VecScatterBegin(vec_scat, vec_vert2cell, vec_vert2cell_reg_subset, &
                          INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vec_scat, vec_vert2cell, vec_vert2cell_reg_subset, &
                          INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_vert2cell_reg_subset.out', &
                              viewer, ierr)
    call VecView(vec_vert2cell_reg_subset, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
      
    call VecGetArrayF90(vec_vert2cell_reg_subset, v_loc_p, ierr)
    allocate(cell_count(region%num_verts))
    allocate(cell_ids(region%num_verts*MAX_CELLS_SHARING_A_VERTEX))
      
    ! initialize
    cell_count    = 0
    count         = 0
      
    do ii = 1,region%num_verts
      do jj = 1,MAX_CELLS_SHARING_A_VERTEX
        if (v_loc_p( (ii-1)*MAX_CELLS_SHARING_A_VERTEX + jj) == -1) exit
        count = count + 1
        cell_count(ii) = cell_count(ii) + 1
        cell_ids(count) = v_loc_p( (ii-1)*MAX_CELLS_SHARING_A_VERTEX + jj )
      enddo
      if (cell_count(ii) < 0) then
        option%io_buffer = 'For a given vertex, no cell found. Stopping'
        call printErrMsg(option)
      endif
    enddo
    call VecRestoreArrayF90(vec_vert2cell_reg_subset, v_loc_p, ierr)
      
    allocate(tmp_int_array(count*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE))
    kk = 0
    do ii = 1,count
      do jj = 1,ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE
        kk = kk + 1
        tmp_int_array(kk) = cell_ids(ii) * ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE + jj
      enddo
    enddo
      
    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, count*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE,&
      tmp_int_array, PETSC_COPY_VALUES, is_from, ierr)
    deallocate(tmp_int_array)
      
    call VecCreateMPI(option%mycomm, count*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE, &
                      PETSC_DECIDE, vec_cell2facevert_reg_subset, ierr)
    call VecGetOwnershipRange(vec_cell2facevert_reg_subset, istart, &
                              iend, ierr)
    allocate(tmp_int_array(count*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE))
    do ii = 1, count*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE
      tmp_int_array(ii) = ii + istart
    enddo

    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, count*ugrid%max_ndual_per_cell*MAX_VERT_PER_FACE, &
                        tmp_int_array, PETSC_COPY_VALUES, is_to, ierr)
    deallocate(tmp_int_array)
      
    call VecScatterCreate(vec_cell2facevert, &
                          is_from,vec_cell2facevert_reg_subset, is_to, &
                          vec_scat, ierr)
    call ISDestroy(is_from,ierr)
    call ISDestroy(is_to,ierr)
    
    call VecScatterBegin(vec_scat, vec_cell2facevert, &
                          vec_cell2facevert_reg_subset, INSERT_VALUES, &
                          SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vec_scat, vec_cell2facevert, &
                          vec_cell2facevert_reg_subset, INSERT_VALUES, &
                          SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'vec_cell2facevert_reg_subset.out', viewer, &
                              ierr)
    call VecView(vec_cell2facevert_reg_subset, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    call VecGetArrayF90(vec_cell2facevert_reg_subset, v_loc_p, ierr)
    istart = 0
    iend   = 0
    counter1 = 0
    counter2 = 0
      
    allocate(cell_ids_for_face(region%num_verts))
    allocate(face_ids_for_face(region%num_verts))
    do jj = 1,region%num_verts
      iend   = istart + cell_count(jj) * ugrid%max_ndual_per_cell * MAX_VERT_PER_FACE
      istart = istart +  1
      counter1 = (istart-1)/MAX_VERT_PER_FACE/ugrid%max_ndual_per_cell
      counter2 = 0
      found  = PETSC_FALSE
      do ii = istart,iend,MAX_VERT_PER_FACE
        counter2 = counter2 + 1
        if (region%vertex_ids(0, jj) == 4) then
          if ((v_loc_p(ii    ) == region%vertex_ids(1,jj)).and.&
              (v_loc_p(ii + 1) == region%vertex_ids(2,jj)).and.&
              (v_loc_p(ii + 2) == region%vertex_ids(3,jj)).and.&
              (v_loc_p(ii + 3) == region%vertex_ids(4,jj)))      then
              found = PETSC_TRUE
              exit
          endif
          if ((v_loc_p(ii    ) == region%vertex_ids(2,jj)).and.&
              (v_loc_p(ii + 1) == region%vertex_ids(3,jj)).and.&
              (v_loc_p(ii + 2) == region%vertex_ids(4,jj)).and.&
              (v_loc_p(ii + 3) == region%vertex_ids(1,jj)))      then
              found = PETSC_TRUE
              exit
          endif
          if ((v_loc_p(ii    ) == region%vertex_ids(3,jj)).and.&
              (v_loc_p(ii + 1) == region%vertex_ids(4,jj)).and.&
              (v_loc_p(ii + 2) == region%vertex_ids(2,jj)).and.&
              (v_loc_p(ii + 3) == region%vertex_ids(1,jj)))      then
              found = PETSC_TRUE
              exit
          endif
          if ((v_loc_p(ii    ) == region%vertex_ids(4,jj)).and.&
              (v_loc_p(ii + 1) == region%vertex_ids(1,jj)).and.&
              (v_loc_p(ii + 2) == region%vertex_ids(2,jj)).and.&
              (v_loc_p(ii + 3) == region%vertex_ids(3,jj)))      then
              found = PETSC_TRUE
              exit
          endif
        elseif ( v_loc_p(ii + 3) < 0 ) then
          if ((v_loc_p(ii    ) == region%vertex_ids(1,jj)).and.&
              (v_loc_p(ii + 1) == region%vertex_ids(2,jj)).and.&
              (v_loc_p(ii + 2) == region%vertex_ids(3,jj)))      then
              found = PETSC_TRUE
              exit
          endif
          if ((v_loc_p(ii    ) == region%vertex_ids(2,jj)).and.&
              (v_loc_p(ii + 1) == region%vertex_ids(3,jj)).and.&
              (v_loc_p(ii + 2) == region%vertex_ids(1,jj)))      then
              found = PETSC_TRUE
              exit
          endif
          if ((v_loc_p(ii    ) == region%vertex_ids(3,jj)).and.&
              (v_loc_p(ii + 1) == region%vertex_ids(1,jj)).and.&
              (v_loc_p(ii + 2) == region%vertex_ids(2,jj)))      then
              found = PETSC_TRUE
              exit
          endif
        endif
          
        if (counter2 == ugrid%max_ndual_per_cell) then
          counter2 = 0
          counter1 = counter1 + 1
        endif
        
      enddo

      if (.not.found) then
        option%io_buffer='No cell found for vertex '
        call printErrMsg(option)
      endif
        
      cell_ids_for_face(jj) = cell_ids(counter1+1)
      face_ids_for_face(jj) = (ii-istart)/MAX_VERT_PER_FACE & 
                              - INT((ii-istart)/MAX_VERT_PER_FACE/ &
                                    ugrid%max_ndual_per_cell)* &
                                    ugrid%max_ndual_per_cell + 1
      istart = iend
    enddo
    call VecRestoreArrayF90(vec_cell2facevert_reg_subset, v_loc_p, ierr)


    call VecCreateMPI(option%mycomm, ugrid%nlmax, PETSC_DECIDE, &
                      vec_cell_ids, ierr)
    call VecCreateMPI(option%mycomm, ugrid%nlmax, PETSC_DECIDE, &
                      vec_cell_ids_loc, ierr)
    call VecCreateMPI(option%mycomm, ugrid%nlmax, PETSC_DECIDE, &
                      vec_face_ids, ierr)
    call VecCreateMPI(option%mycomm, ugrid%nlmax, PETSC_DECIDE, &
                      vec_face_ids_loc, ierr)
    
    call VecSet(vec_cell_ids, 0.d0, ierr)
    call VecAssemblyBegin(vec_cell_ids, ierr)
    call VecAssemblyEnd(  vec_cell_ids, ierr)
    
    call VecSet(vec_face_ids, 0.d0, ierr)
    call VecAssemblyBegin(vec_face_ids, ierr)
    call VecAssemblyEnd(  vec_face_ids, ierr)

    allocate(tmp_int_array(region%num_verts))
    allocate(tmp_scl_array(region%num_verts))

    do ii = 1, region%num_verts
      tmp_int_array(ii)  = cell_ids_for_face(ii)
      tmp_scl_array(ii)  = 1.d0
    enddo

    call VecSetValues(vec_cell_ids, region%num_verts, tmp_int_array, &
                      tmp_scl_array, ADD_VALUES, ierr)

    do ii = 1, region%num_verts
      tmp_int_array(ii) = cell_ids_for_face(ii)
      tmp_scl_array(ii) = face_ids_for_face(ii)
    enddo

    call VecSetValues(vec_face_ids, region%num_verts, tmp_int_array, &
                      tmp_scl_array, ADD_VALUES, ierr)
    
    deallocate(tmp_int_array)
    deallocate(tmp_scl_array)

    call VecAssemblyBegin(vec_cell_ids, ierr)
    call VecAssemblyEnd(  vec_cell_ids, ierr)
    call VecAssemblyBegin(vec_face_ids, ierr)
    call VecAssemblyEnd(  vec_face_ids, ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_cell_ids_aft_2.out', &
                              viewer, ierr)
    call VecView(vec_cell_ids, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    allocate(tmp_int_array(ugrid%nlmax))
    count = 0
    do ghosted_id=1,ugrid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (local_id < 1) cycle
      count = count + 1
      natural_id = grid%nG2A(ghosted_id)
      tmp_int_array(count) = natural_id
    enddo

    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, ugrid%nlmax, &
                      tmp_int_array, PETSC_COPY_VALUES, &
                      is_from, ierr)
    
    call VecGetOwnershipRange(vec_cell_ids_loc, istart, iend, ierr)
    do ii = 1, ugrid%nlmax
      tmp_int_array(ii) = ii + istart
    enddo

    tmp_int_array = tmp_int_array - 1
    call ISCreateBlock(option%mycomm, 1, ugrid%nlmax, &
                        tmp_int_array, PETSC_COPY_VALUES, &
                        is_to, ierr)
    deallocate(tmp_int_array)
    
    call VecScatterCreate(vec_cell_ids, is_from, vec_cell_ids_loc, &
                          is_to, vec_scat, ierr)
    call ISDestroy(is_from,ierr)
    call ISDestroy(is_to,ierr)
    
    call VecScatterBegin(vec_scat, vec_cell_ids, vec_cell_ids_loc, &
                          INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vec_scat, vec_cell_ids, vec_cell_ids_loc, &
                          INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(vec_scat, vec_face_ids, vec_face_ids_loc, &
                          INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vec_scat, vec_face_ids, vec_face_ids_loc, &
                          INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

#if GB_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, 'vec_cell_ids_loc_2.out', &
                              viewer, ierr)
    call VecView(vec_cell_ids_loc, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    call VecGetArrayF90(vec_cell_ids_loc, v_loc_p, ierr)
    call VecGetArrayF90(vec_face_ids_loc, v_loc2_p, ierr)
    count = 0
    do ii = 1, ugrid%nlmax
      if (v_loc_p(ii) == 1) count = count + 1
    enddo
    
    region%num_cells = count
    if (count > 0) then
      allocate(tmp_int_array(count))
      allocate(tmp_int_array2(count))
      count = 0
      do ii = 1, ugrid%nlmax
        if (v_loc_p(ii) == 1) then
          count = count + 1
          tmp_int_array(count) = ii
          tmp_int_array2(count) = v_loc2_p(ii)
        endif
      enddo

      !deallocate(region%cell_ids)
      allocate(region%cell_ids(region%num_cells))
      allocate(region%faces(region%num_cells))
      region%cell_ids = tmp_int_array
      region%faces    = tmp_int_array2
        
      deallocate(tmp_int_array)
      deallocate(tmp_int_array2)

    else
      nullify(region%cell_ids)
    endif
    
    call VecRestoreArrayF90(vec_cell_ids_loc, v_loc_p, ierr)


    deallocate(cell_count)
    deallocate(cell_ids)
      
    call VecDestroy(vec_cell2facevert, ierr)
    call VecDestroy(vec_vert2cell, ierr)
    deallocate(vert2cell_array)

  endif


end subroutine GridLocalizeRegionsForUGrid

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
  
  call GridVecGetArrayF90(grid, vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call GridVecRestoreArrayF90(grid, vector,vec_ptr,ierr)
  
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
  
  call GridVecGetArrayF90(grid,vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call GridVecRestoreArrayF90(grid,vector,vec_ptr,ierr)
  
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
  
  call GridVecGetArrayF90(grid,vector,vec_ptr,ierr)
  do i=1,num_values
    if (vec_ptr(i) > 0.d0) then
      array(i) = int(vec_ptr(i)+1.d-4)
    else
      array(i) = int(vec_ptr(i)-1.d-4)
    endif
  enddo
  call GridVecRestoreArrayF90(grid,vector,vec_ptr,ierr)
  
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
  
  call GridVecGetArrayF90(grid,vector,vec_ptr,ierr)
  array(1:num_values) = vec_ptr(1:num_values)
  call GridVecRestoreArrayF90(grid,vector,vec_ptr,ierr)
  
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
    if (natural_id == grid%nL2A(local_id)+1) then
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
    case(UNSTRUCTURED_GRID) 
      option%io_buffer = 'GridGetNeighbors not currently supported for ' // &
        'unstructured grids.'
      call printErrMsg(option)
  end select

end subroutine GridGetGhostedNeighbors

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
    
  if (.not.associated(grid)) return
      
  if (associated(grid%nL2G)) deallocate(grid%nL2G)
  nullify(grid%nL2G)
  if (associated(grid%nG2L)) deallocate(grid%nG2L)
  nullify(grid%nG2L)
  if (associated(grid%nL2A)) deallocate(grid%nL2A)
  nullify(grid%nL2A)
  if (associated(grid%nG2A)) deallocate(grid%nG2A)
  nullify(grid%nG2A)
  if (associated(grid%nG2P)) deallocate(grid%nG2P)
  nullify(grid%nG2P)

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
  call StructuredGridDestroy(grid%structured_grid)
                                           
  call ConnectionDestroyList(grid%internal_connection_set_list)

end subroutine GridDestroy


! ************************************************************************** !
!
! GridDestroy: Returns pointer to cell-centered vector values
! author: Bobby Philip
! date:  12/15/10
!
! ************************************************************************** !
subroutine GridVecGetMaskArrayCellF90(grid, vec, f90ptr, ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  if (associated(grid%structured_grid)) then
     call StructuredGridVecGetMaskArrayCellF90(grid%structured_grid, vec, f90ptr, ierr)
  endif

end subroutine GridVecGetMaskArrayCellF90

      
! ************************************************************************** !
!
! GridDestroy: Returns pointer to cell-centered vector values
! author: Bobby Philip
! date: 
!
! ************************************************************************** !
subroutine GridVecGetArrayCellF90(grid, vec, f90ptr, ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  if (.not.associated(grid%structured_grid)) then
     call VecGetArrayF90(vec, f90ptr, ierr)
  else
     call StructuredGridVecGetArrayF90(grid%structured_grid, vec, f90ptr, ierr)
  endif

end subroutine GridVecGetArrayCellF90

      
! ************************************************************************** !
!
! GridDestroy: Returns pointer to edge-based vector values?
! author: Bobby Philip
! date: 
!
! ************************************************************************** !
subroutine GridVecGetArraySideF90(grid, axis, vec, f90ptr, ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  PetscInt :: axis 
  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  if (.not.associated(grid%structured_grid)) then
     call VecGetArrayF90(vec, f90ptr, ierr)
  else
     call StructuredGridVecGetArrayF90(grid%structured_grid, axis, vec, f90ptr, ierr)
  endif

end subroutine GridVecGetArraySideF90

! ************************************************************************** !
!
! GridDestroy: Restores pointer to vector values
! author: Bobby Philip
! date: 
!
! ************************************************************************** !
subroutine GridVecRestoreArrayF90(grid, vec, f90ptr, ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  if (.not.associated(grid%structured_grid)) then
     call VecRestoreArrayF90(vec, f90ptr, ierr)
  else
     call StructGridVecRestoreArrayF90(grid%structured_grid, vec, f90ptr, ierr)
  endif

end subroutine GridVecRestoreArrayF90

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
      cell_id = grid%nL2A(cell_id)+1
    else if (vec_type == LOCAL) then
      cell_id = grid%nG2A(cell_id) !nG2A is 1-based
    endif
  endif
  
  call MPI_Allreduce(cell_id,GridIndexToCellID,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,PETSC_COMM_WORLD,ierr)
                     
end function GridIndexToCellID


!! ********************************************************************** !
!!
!! GridPopulateSubcontinuum: Populate subcontinuum discretization info in
!!                           the grid object
!! author: Jitendra Kumar
!! date: 11/22/2010
!!
!! *********************************************************************** !
!subroutine GridPopulateSubcontinuum(realization)
!
!  use Realization_module
!  use Discretization_module
!  use Material_module
!  use Grid_module
!  use Patch_module
!  use Level_module
!
!  implicit none
!
!  type(realization_type), pointer :: realization
!  type(grid_type), pointer :: grid
!  type(discretization_type), pointer :: discretization
!  type(patch_type), pointer :: patch 
!  type(level_type), pointer :: cur_level
!  type(patch_type), pointer :: cur_patch
!  type(subcontinuum_type), pointer :: subcontinuum_type
!
!  PetscInt :: icell, ssub
!
!  ! do this only at the last and finest level
!  cur_level => realization%level_list%last
!  if (.not.associated(cur_level)) exit
!  cur_patch => cur_level%patch_list%first
!  do 
!    grid => cur_patch%grid
!    ! Allocate storage for subcontinuum_grid_offset
!    allocate(cur_patch%grid%subcontinuum_grid_offset(cur_patch%grid%nlmax,2)
!    
!    ! Loop through all cells in the patch, copy the no. of subcontinuum 
!    ! and offset from patch%num_subcontinuum_type into
!    ! grid%subcontinuum_grid_offset.
!    ssub = 0
!    do icell=1, cur_patch%grid%nlmax
!      cur_patch%grid%subcontinuum_grid_offset(icell,1) =   &
!                          cur_patch%num_subcontinuum_type(icell,1)
!      cur_patch%grid%subcontinuum_grid_offset(icell,2) =   &
!                          cur_patch%num_subcontinuum_type(icell,2)
!      ssub = ssub + cur_patch%grid%subcontinuum_grid_offset(icell,1)
!    enddo
!
!    ! Allocate storage for grid%subcontinuum_grid and store the subgrid
!    ! information
!    allocate(cur_patch%grid%subcontinuum_grid(ssub))
!    
!    ! Loop through all the cells-> all subcontinuum and set the subgrid
!    ! information
!    do icell=1, cur_patch%grid%nlmax
!      if (associated(region)) the
!        local_id = region%cell_ids(icell)
!      else
!        local_id = icell
!      endif
!
!      ! Loop over all subcontinua
!      offset = cur_patch%grid%subcontinuum_grid_offset(icell,2)
!      do isub = 1, cur_patch%grid%subcontinuum_grid_offset(icell,1)
!        cur_patch%grid%subcontinuum_grid(offset) =  &
!        realization%subcontinuum_properties(cur_patch%subcontinuum_type_ids(offset))%num_subgrids
!        offset = offset + 1
!      enddo
!    enddo
!    cur_patch => cur_patch%next
!  enddo           
!end subroutine GridPopulateSubcontinuum
!
end module Grid_module
