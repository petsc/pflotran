module Simulation_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type,public :: simulation_aux_type

    ! Note: These are LOCAL vectors (i.e. they do not contain ghost control 
    !       volumes)

    ! Size of entire subsurface domain
    Vec :: subsurf_pres
    Vec :: subsurf_temp
    Vec :: subsurf_sat
    Vec :: subsurf_den

    ! Size of surface cells of subsurface domain
    Vec :: subsurf_pres_top_bc
    Vec :: subsurf_temp_top_bc
    Vec :: subsurf_mflux_exchange_with_surf
    Vec :: subsurf_hflux_exchange_with_surf

    ! Size of entire surface domain
    Vec :: surf_head
    Vec :: surf_temp
    Vec :: surf_mflux_exchange_with_subsurf
    Vec :: surf_hflux_exchange_with_subsurf

    VecScatter :: surf_to_subsurf
    VecScatter :: subsurf_to_surf
    VecScatter :: subsurf_to_hydrogeophyics

  end type simulation_aux_type

#ifdef SURFACE_FLOW  
  interface SimAuxCreateVecScatters
    module procedure SimAuxCreateSurfSubSurfVScats
  end interface
#endif  

  public :: SimAuxCreate, &
#ifdef SURFACE_FLOW  
            SimAuxCreateVecScatters, &
            SimAuxCreateSubSurfVecs, &
            SimAuxCreateSurfVecs, &
#endif
            SimAuxDestroy

contains

! ************************************************************************** !
!> This routine allocates auxillary object.
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
function SimAuxCreate()

  use Option_module

  implicit none

  type (simulation_aux_type),pointer :: SimAuxCreate

  type (simulation_aux_type),pointer :: aux

  allocate(aux)
  aux%subsurf_pres = 0
  aux%subsurf_temp = 0
  aux%subsurf_sat = 0
  aux%subsurf_den = 0
  aux%subsurf_pres_top_bc = 0
  aux%subsurf_temp_top_bc = 0
  aux%subsurf_mflux_exchange_with_surf = 0
  aux%subsurf_hflux_exchange_with_surf = 0

  aux%surf_head = 0
  aux%surf_temp = 0
  aux%surf_mflux_exchange_with_subsurf = 0
  aux%surf_hflux_exchange_with_subsurf = 0

  aux%surf_to_subsurf = 0
  aux%subsurf_to_surf = 0
  aux%subsurf_to_hydrogeophyics = 0

  SimAuxCreate => aux

end function SimAuxCreate

#ifdef SURFACE_FLOW
! ************************************************************************** !
!> This routine creates VecScatter between surface-subsurface grids.
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
!!
!! Algorithm:
!!  - It uses a similar logic of Matrix-Vector multiplication used in 
!!    UGridMapSideSet() subroutine. The algorithm here is extended to use 
!!    Matrix-Matrix mulitplication
!!
! ************************************************************************** !
subroutine SimAuxCreateSurfSubSurfVScats(pm_aux,realization, &
                                            surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Surface_Realization_class

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type (simulation_aux_type),pointer :: pm_aux
  type(realization_type),pointer         :: realization
  type(surface_realization_type),pointer :: surf_realization

  type(option_type),pointer           :: option
  type(unstructured_grid_type),pointer :: subsurf_grid
  type(unstructured_grid_type),pointer :: surf_grid
  type(level_type),pointer            :: cur_level
  type(patch_type),pointer            :: cur_patch 
  type(region_type),pointer           :: cur_region,top_region
  type(region_type),pointer           :: patch_region

  Mat :: Mat_vert_to_face_subsurf
  Mat :: Mat_vert_to_face_subsurf_transp
  Mat :: Mat_vert_to_face_surf
  Mat :: Mat_vert_to_face_surf_transp
  Mat :: prod
  Vec :: subsurf_petsc_ids,surf_petsc_ids

  PetscViewer :: viewer


  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4,1)
  PetscInt :: nvertices
  PetscInt :: iface
  PetscInt :: local_id,ii,jj
  PetscInt :: cell_type
  PetscInt :: ivertex,vertex_id_local
  PetscReal :: real_array4(4)
  PetscReal,pointer :: vec_ptr(:)

  PetscErrorCode :: ierr
  PetscBool :: found
  
  found = PETSC_FALSE
  
  if (.not.associated(realization)) return
  
  option => realization%option
  subsurf_grid => realization%discretization%grid%unstructured_grid
  surf_grid    => surf_realization%discretization%grid%unstructured_grid

  ! localize the regions on each patch
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      cur_region => cur_patch%regions%first
        do
          if (.not.associated(cur_region)) exit
          if (StringCompare(cur_region%name,'top')) then
            found = PETSC_TRUE
            top_region => cur_region
            exit
          endif
          cur_region => cur_region%next
        enddo
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (found.eqv.PETSC_FALSE) then
    option%io_buffer = 'When running with -DSURFACE_FLOW need to specify ' // &
      ' in the inputfile explicitly region: top '
    call printErrMsg(option)
  endif

  call MatCreateAIJ(option%mycomm, &
                       top_region%num_cells, &
                       PETSC_DECIDE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_subsurf, &
                       ierr)
  call MatCreateAIJ(option%mycomm, &
                       PETSC_DECIDE, &
                       top_region%num_cells, &
                       subsurf_grid%num_vertices_global, &
                       PETSC_DETERMINE, &
                       12, &
                       PETSC_NULL_INTEGER, &
                       12, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_subsurf_transp, &
                       ierr)
  call VecCreateMPI(option%mycomm,top_region%num_cells,PETSC_DETERMINE, &
                    subsurf_petsc_ids,ierr)  
  call MatZeroEntries(Mat_vert_to_face_subsurf,ierr)
  real_array4 = 1.d0

  call VecGetArrayF90(subsurf_petsc_ids,vec_ptr,ierr)

  offset=0
  call MPI_Exscan(top_region%num_cells,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  do ii = 1,top_region%num_cells
    local_id = top_region%cell_ids(ii)
    vec_ptr(ii) = subsurf_grid%cell_ids_petsc(local_id)
    iface    = top_region%faces(ii)
    cell_type = subsurf_grid%cell_type(local_id)
    !nfaces = UCellGetNFaces(cell_type,option)

    call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                    int_array4)
    ! For this matrix:
    !   irow = local face id
    !   icol = natural (global) vertex id
    do ivertex = 1,nvertices
      vertex_id_local = &
        subsurf_grid%cell_vertices(int_array4(ivertex),local_id)
      int_array4_0(ivertex,1) = &
        subsurf_grid%vertex_ids_natural(vertex_id_local)-1
    enddo
    call MatSetValues(Mat_vert_to_face_subsurf, &
                      1,ii-1+offset, &
                      nvertices,int_array4_0, &
                      real_array4, &
                      INSERT_VALUES,ierr)
    call MatSetValues(Mat_vert_to_face_subsurf_transp, &
                      nvertices,int_array4_0, &
                      1,ii-1+offset, &
                      real_array4, &
                      INSERT_VALUES,ierr)
  enddo

  call VecRestoreArrayF90(subsurf_petsc_ids,vec_ptr,ierr)

  call MatAssemblyBegin(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyBegin(Mat_vert_to_face_subsurf_transp,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_subsurf_transp,MAT_FINAL_ASSEMBLY,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_subsurf.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_subsurf,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  !call MatTranspose(Mat_vert_to_face_subsurf,MAT_INITIAL_MATRIX, &
  !                  Mat_vert_to_face_subsurf_transp,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_subsurf_transp.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_subsurf_transp,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  string = 'subsurf_petsc_ids.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecView(subsurf_petsc_ids,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif


  call MatCreateAIJ(option%mycomm, &
                       surf_grid%nlmax, &
                       PETSC_DECIDE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_surf, &
                       ierr)
  call MatCreateAIJ(option%mycomm, &
                       PETSC_DECIDE, &
                       surf_grid%nlmax, &
                       subsurf_grid%num_vertices_global, &
                       PETSC_DETERMINE, &
                       12, &
                       PETSC_NULL_INTEGER, &
                       12, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_surf_transp, &
                       ierr)
  call VecCreateMPI(option%mycomm,surf_grid%nlmax,PETSC_DETERMINE, &
                    surf_petsc_ids,ierr)  
  offset=0
  call MPI_Exscan(surf_grid%nlmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  call VecGetArrayF90(surf_petsc_ids,vec_ptr,ierr)

  do local_id = 1,surf_grid%nlmax
    cell_type = surf_grid%cell_type(local_id)
    vec_ptr(local_id) = surf_grid%cell_ids_petsc(local_id)
    
    int_array4_0 = 0
    nvertices = surf_grid%cell_vertices(0,local_id)
    do ivertex = 1,nvertices
      vertex_id_local = surf_grid%cell_vertices(ivertex,local_id)
      int_array4_0(ivertex,1) = &
        surf_grid%vertex_ids_natural(vertex_id_local)-1
    enddo    
    call MatSetValues(Mat_vert_to_face_surf, &
                      1,local_id-1+offset, &
                      nvertices,int_array4_0, &
                      real_array4, &
                      INSERT_VALUES,ierr)
    call MatSetValues(Mat_vert_to_face_surf_transp, &
                      nvertices,int_array4_0, &
                      1,local_id-1+offset, &
                      real_array4, &
                      INSERT_VALUES,ierr)
  enddo

  call VecRestoreArrayF90(surf_petsc_ids,vec_ptr,ierr)

  call MatAssemblyBegin(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyBegin(Mat_vert_to_face_surf_transp,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_surf_transp,MAT_FINAL_ASSEMBLY,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_surf.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_surf,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  string = 'surf_petsc_ids.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecView(surf_petsc_ids,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  !call MatTranspose(Mat_vert_to_face_surf,MAT_INITIAL_MATRIX, &
  !                  Mat_vert_to_face_surf_transp,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_surf_transp.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_surf_transp,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  
  call MatMatMult(Mat_vert_to_face_subsurf,Mat_vert_to_face_surf_transp, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,prod,ierr)

#if UGRID_DEBUG
  string = 'prod.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(prod,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  call SimulationCreateSurfSubSurfVScat(realization,surf_realization,prod, &
                                     surf_petsc_ids,pm_aux%surf_to_subsurf)
  call MatDestroy(prod,ierr)

  call MatMatMult(Mat_vert_to_face_surf,Mat_vert_to_face_subsurf_transp, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,prod,ierr)

#if UGRID_DEBUG
  string = 'prod_2.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(prod,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  call SimulationCreateSurfSubSurfVScat(realization,surf_realization,prod, &
                                        subsurf_petsc_ids,pm_aux%subsurf_to_surf)

  call MatDestroy(prod,ierr)

  call MatDestroy(Mat_vert_to_face_subsurf,ierr)
  call MatDestroy(Mat_vert_to_face_subsurf_transp,ierr)
  call MatDestroy(Mat_vert_to_face_surf,ierr)
  call MatDestroy(Mat_vert_to_face_surf_transp,ierr)

  call VecDestroy(subsurf_petsc_ids,ierr)
  call VecDestroy(surf_petsc_ids,ierr)

end subroutine SimAuxCreateSurfSubSurfVScats

! ************************************************************************** !
!> This subroutine creates a single vector scatter context
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
subroutine SimulationCreateSurfSubSurfVScat( &
              realization, &       !<
              surf_realization, &  !<
              prod_mat, &          !< Mat-Mat-Mult matrix
              source_petsc_ids, &   !< MPI-Vector containing cell ids in PETSc order
              scatter &
              )

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Field_module
  use Surface_Field_module
  use Unstructured_Grid_module
  use Discretization_module
  use Unstructured_Grid_Aux_module
  use DM_Kludge_module
  use Surface_Realization_class

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type),pointer         :: realization
  type(surface_realization_type),pointer :: surf_realization
  Mat :: prod_mat
  Vec :: source_petsc_ids
  VecScatter :: scatter

  Mat :: prod_loc_mat
  Vec :: source_loc_vec
  Vec :: corr_dest_ids_vec
  Vec :: corr_dest_ids_vec_ndof
  Vec :: source_petsc_ids_ndof
  IS  :: is_tmp1,is_tmp2
  IS  :: is_tmp3,is_tmp4
  PetscInt,pointer :: corr_v2_ids(:)
  VecScatter :: scatter_ndof

  PetscViewer :: viewer

  type(option_type),pointer :: option
  type(field_type),pointer :: field
  type(surface_field_type),pointer :: surf_field

  type(dm_ptr_type),pointer :: dm_ptr
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4,1)
  PetscReal :: real_array4(4)
  PetscInt :: ii,jj
  PetscReal,pointer :: vec_ptr(:)
  PetscInt :: ivertex,cell_id,vertex_id_local
  PetscReal :: max_value

  PetscInt,pointer               :: ia_p(:),ja_p(:)
  PetscInt                        :: nrow,rstart,rend,icol(1)
  PetscInt                        :: index
  PetscInt                        :: vertex_id
  PetscOffset                     :: iia,jja,aaa,iicol
  PetscBool                       :: done
  PetscScalar                     :: aa(1)

  PetscErrorCode :: ierr
  PetscBool :: found
  PetscInt :: nlocal

  option     => realization%option
  field      => realization%field
  surf_field => surf_realization%surf_field
  
  if (option%mycommsize > 1) then
    ! From the MPI-Matrix get the local-matrix
    call MatMPIAIJGetLocalMat(prod_mat,MAT_INITIAL_MATRIX,prod_loc_mat,ierr)
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_loc_mat,ONE_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                        nrow,ia_p,ja_p,done,ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArray(prod_loc_mat,aa,aaa,ierr)
  else
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_mat,ONE_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                        nrow,ia_p,ja_p,done,ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArray(prod_mat,aa,aaa,ierr)
  endif

  ! For each row of the local-matrix,find the column with the largest value
  allocate(corr_v2_ids(nrow))
  do ii = 1,nrow
    max_value = 0.d0
    do jj = ia_p(ii),ia_p(ii + 1) - 1
      if (aa(aaa+ jj) > max_value) then
        corr_v2_ids(ii) = ja_p(jj)
        max_value = aa(aaa+ jj)
      endif
    enddo
    if (max_value<3) then
      option%io_buffer = 'Atleast three vertices need to form a face'
      call printErrMsg(option)
    endif
  enddo

  offset = 0
  call MPI_Exscan(nrow,offset,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  allocate(int_array(nrow))
  do ii = 1,nrow
    int_array(ii) = (ii-1)+offset
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr)
  call ISCreateBlock(option%mycomm,option%nflowdof,nrow, &
                     int_array,PETSC_COPY_VALUES,is_tmp3,ierr)

  do ii = 1,nrow
    int_array(ii) = corr_v2_ids(ii)-1
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr)
  call ISCreateBlock(option%mycomm,option%nflowdof,nrow, &
                     int_array,PETSC_COPY_VALUES,is_tmp4,ierr)
  deallocate(int_array)
  
  call VecCreateMPI(option%mycomm,nrow,PETSC_DETERMINE, &
                    corr_dest_ids_vec,ierr)
  call VecScatterCreate(source_petsc_ids,is_tmp2,corr_dest_ids_vec,is_tmp1, &
                        scatter,ierr)
  call ISDestroy(is_tmp1,ierr)
  call ISDestroy(is_tmp2,ierr)
  
  call VecScatterBegin(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecDestroy(corr_dest_ids_vec,ierr)
  if (option%mycommsize>1) call MatDestroy(prod_loc_mat,ierr)

end subroutine SimulationCreateSurfSubSurfVScat

! ************************************************************************** !
!> This routine creates subsurface vectors.
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
subroutine SimAuxCreateSubSurfVecs(aux,subsurf_realization,option)

  use Realization_class
  use Coupler_module
  use Option_module
  use String_module

  implicit none

  type(simulation_aux_type),pointer :: aux
  type(realization_type),pointer :: subsurf_realization
  type(option_type),pointer :: option

  type(coupler_list_type),pointer :: coupler_list
  type(coupler_type),pointer :: coupler

  PetscInt :: num_conn
  PetscInt :: found
  PetscInt :: found_global
  PetscErrorCode :: ierr

  call VecCreate(option%mycomm,aux%subsurf_pres,ierr)
  call VecSetSizes(aux%subsurf_pres, &
                   subsurf_realization%discretization%grid%nlmax, &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(aux%subsurf_pres,ierr)
  call VecSet(aux%subsurf_pres,0.d0,ierr)

  call VecDuplicate(aux%subsurf_pres,aux%subsurf_temp,ierr)
  call VecDuplicate(aux%subsurf_pres,aux%subsurf_sat,ierr)
  call VecDuplicate(aux%subsurf_pres,aux%subsurf_den,ierr)

  found = 0
  num_conn = 0
  coupler_list => subsurf_realization%patch%source_sinks
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    if (StringCompare(coupler%name,'from_surface_ss')) then
      num_conn = coupler%connection_set%num_connections
      found = 1
    endif
    coupler => coupler%next
  enddo

  call MPI_AllReduce(found,found_global,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                     option%mycomm,ierr)

  if (found_global == 0) then
    coupler_list => subsurf_realization%patch%boundary_conditions
    coupler => coupler_list%first
    do
      if (.not.associated(coupler)) exit
      if (StringCompare(coupler%name,'from_surface_bc')) then
        num_conn = coupler%connection_set%num_connections
        found = 1
      endif
      coupler => coupler%next
    enddo
    call MPI_AllReduce(found,found_global,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                       option%mycomm,ierr)
  endif

  if (found_global > 0) then
    call VecCreate(option%mycomm,aux%subsurf_pres_top_bc,ierr)
    call VecSetSizes(aux%subsurf_pres_top_bc,num_conn,PETSC_DECIDE,ierr)
    call VecSetFromOptions(aux%subsurf_pres_top_bc,ierr)
    call VecSet(aux%subsurf_pres_top_bc,0.d0,ierr)

    if (option%iflowmode == TH_MODE) &
      call VecDuplicate(aux%subsurf_pres_top_bc,aux%subsurf_temp_top_bc,ierr)
    call VecDuplicate(aux%subsurf_pres_top_bc,aux%subsurf_mflux_exchange_with_surf,ierr)
    call VecDuplicate(aux%subsurf_pres_top_bc,aux%subsurf_hflux_exchange_with_surf,ierr)
  endif

end subroutine SimAuxCreateSubSurfVecs

! ************************************************************************** !
!> This routine creates surface vectors.
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
subroutine SimAuxCreateSurfVecs(aux,surf_realization,option)

  use Surface_Realization_class
  use Option_module

  implicit none

  type(simulation_aux_type),pointer :: aux
  type(surface_realization_type),pointer :: surf_realization
  type(option_type),pointer :: option

  PetscErrorCode :: ierr

  call VecCreate(option%mycomm,aux%surf_head,ierr)
  call VecSetSizes(aux%surf_head,surf_realization%discretization%grid%nlmax, &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(aux%surf_head,ierr)
  call VecSet(aux%surf_head,0.d0,ierr)

  call VecDuplicate(aux%surf_head,aux%surf_temp,ierr)
  call VecDuplicate(aux%surf_head,aux%surf_mflux_exchange_with_subsurf,ierr)
  call VecDuplicate(aux%surf_head,aux%surf_hflux_exchange_with_subsurf,ierr)

end subroutine SimAuxCreateSurfVecs
#endif

! ************************************************************************** !
!> This routine deallocates auxillary object.
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
subroutine SimAuxDestroy(aux)

  implicit none

  type(simulation_aux_type) :: aux

  PetscErrorCode :: ierr

  if (aux%subsurf_pres /= 0) call VecDestroy(aux%subsurf_pres,ierr)
  if (aux%subsurf_temp /= 0) call VecDestroy(aux%subsurf_temp,ierr)
  if (aux%subsurf_sat /= 0) call VecDestroy(aux%subsurf_sat,ierr)
  if (aux%subsurf_den /= 0) call VecDestroy(aux%subsurf_den,ierr)
  if (aux%subsurf_pres_top_bc /= 0) call VecDestroy(aux%subsurf_pres_top_bc,ierr)
  if (aux%subsurf_temp_top_bc /= 0) call VecDestroy(aux%subsurf_temp_top_bc,ierr)
  if (aux%subsurf_mflux_exchange_with_surf /= 0) &
    call VecDestroy(aux%subsurf_mflux_exchange_with_surf,ierr)
  if (aux%subsurf_hflux_exchange_with_surf /= 0) &
    call VecDestroy(aux%subsurf_hflux_exchange_with_surf,ierr)

  if (aux%surf_head /= 0) call VecDestroy(aux%surf_head,ierr)
  if (aux%surf_temp /= 0) call VecDestroy(aux%surf_temp,ierr)
  if (aux%surf_mflux_exchange_with_subsurf /= 0) &
    call VecDestroy(aux%surf_mflux_exchange_with_subsurf,ierr)
  if (aux%surf_hflux_exchange_with_subsurf /= 0) &
    call VecDestroy(aux%surf_hflux_exchange_with_subsurf,ierr)

  if (aux%surf_to_subsurf /= 0) call VecScatterDestroy(aux%surf_to_subsurf,ierr)
  if (aux%subsurf_to_surf /= 0) call VecScatterDestroy(aux%subsurf_to_surf,ierr)
  if (aux%subsurf_to_hydrogeophyics /= 0) &
    call VecScatterDestroy(aux%subsurf_to_hydrogeophyics,ierr)

end subroutine SimAuxDestroy

end module Simulation_Aux_module
