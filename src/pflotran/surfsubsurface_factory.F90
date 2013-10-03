#ifdef SURFACE_FLOW
module Surf_Subsurf_Factory_module

  use Surf_Subsurf_Simulation_class

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: SurfSubsurfaceInitialize

contains
! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceInitialize(simulation_base,option)

  use Option_module
  use Input_module
  use Timestepper_Base_class
  use Simulation_Base_class
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(surfsubsurface_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => SurfSubsurfaceSimulationCreate(option)
  call SurfSubsurfaceInitializePostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine SurfSubsurfaceInitialize

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceInitializePostPETSc(simulation, option)

  use Simulation_module
  use Surface_Simulation_class
  use Subsurface_Simulation_class
  use Surface_Factory_module
  use Subsurface_Factory_module
  use Option_module
  use Init_module
  use Surface_Flow_module
  use Surface_TH_module
  use Simulation_Aux_module
  use PMC_Base_class
  use PFLOTRAN_Constants_module
  
  implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(surfsubsurface_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  type(surface_simulation_type) :: surf_simulation
  type(subsurface_simulation_type) :: subsurf_simulation
  type(simulation_type), pointer :: simulation_old
  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscInt :: init_status
  VecScatter :: vscat_surf_to_subsurf
  VecScatter :: vscat_subsurf_to_surf
  Vec :: vec_subsurf_pres
  Vec :: vec_subsurf_pres_top_bc
  Vec :: vec_surf_head
  PetscErrorCode :: ierr
  
  ! process command line arguments specific to subsurface
  !call SurfSubsurfInitCommandLineSettings(option)
  
  allocate(simulation_old)
  simulation_old => SimulationCreate(option)
  call Init(simulation_old)

  call HijackSimulation(simulation_old,subsurf_simulation)
  call SubsurfaceJumpStart(subsurf_simulation)

   simulation%realization => simulation_old%realization
   simulation%flow_process_model_coupler => &
        subsurf_simulation%flow_process_model_coupler
   simulation%rt_process_model_coupler => &
        subsurf_simulation%rt_process_model_coupler
   simulation%regression => simulation_old%regression

  if(option%nsurfflowdof>0) then
    ! Both, Surface-Subsurface flow active
    call HijackSurfaceSimulation(simulation_old,surf_simulation)
    call SurfaceJumpStart(surf_simulation)

    simulation%process_model_coupler_list => &
      surf_simulation%process_model_coupler_list
    surf_simulation%process_model_coupler_list%next => &
      subsurf_simulation%process_model_coupler_list
    surf_simulation%surf_flow_process_model_coupler%subsurf_realization => &
      simulation_old%realization
    simulation%flow_process_model_coupler%realization => &
      simulation_old%realization
    simulation%process_model_coupler_list%is_master = PETSC_TRUE
    subsurf_simulation%process_model_coupler_list%is_master = PETSC_FALSE

    simulation%surf_realization => simulation_old%surf_realization
    simulation%surf_flow_process_model_coupler => &
         surf_simulation%surf_flow_process_model_coupler

    nullify(surf_simulation%process_model_coupler_list)
  
    if (option%subsurf_surf_coupling == SEQ_COUPLED .or. &
        option%subsurf_surf_coupling == SEQ_COUPLED_NEW) then
       select case(option%iflowmode)
         case (RICHARDS_MODE)
            call SurfaceFlowGetSubsurfProp(simulation%realization, &
                 simulation%surf_realization)
         case (TH_MODE)
            call SurfaceTHGetSubsurfProp(simulation%realization, &
                 simulation%surf_realization)
         end select
      end if
   else
      ! Only subsurface flow active
      simulation%process_model_coupler_list => &
           subsurf_simulation%process_model_coupler_list
      ! call printErrMsg(option,'Only subsurface-flow is active. ' // &
      !        'Check inputfile or switch -simulation_mode subsurface')
      nullify(simulation%surf_realization)
      nullify(simulation%surf_flow_process_model_coupler)
   endif

   nullify(subsurf_simulation%process_model_coupler_list)

  ! sim_aux: Create PETSc Vectors
  call SurfSubsurfCreateSubsurfVecs(simulation_old%realization, option, &
                                    vec_subsurf_pres, vec_subsurf_pres_top_bc)
  call SimAuxCopySubsurfVec(simulation%sim_aux, vec_subsurf_pres)
  call SimAuxCopySubsurfTopBCVec(simulation%sim_aux, vec_subsurf_pres_top_bc)
  call VecDestroy(vec_subsurf_pres, ierr)
  call VecDestroy(vec_subsurf_pres_top_bc, ierr)

  ! sim_aux: Create PETSc VectorScatters
  if(option%nsurfflowdof>0) &

    call SurfSubsurfCreateSurfVecs(simulation_old%surf_realization, option, &
                                   vec_surf_head)
    call SimAuxCopySurfVec(simulation%sim_aux, vec_surf_head)
    call VecDestroy(vec_surf_head, ierr)

    call SurfSubsurfCreateSurfSubSurfVScats(simulation_old%realization, &
          simulation_old%surf_realization, vscat_surf_to_subsurf, vscat_subsurf_to_surf)
    call SimAuxCopyVecScatter(simulation%sim_aux, vscat_surf_to_subsurf, SURF_TO_SUBSURF)
    call SimAuxCopyVecScatter(simulation%sim_aux, vscat_subsurf_to_surf, SUBSURF_TO_SURF)
    call VecScatterDestroy(vscat_surf_to_subsurf, ierr)
    call VecScatterDestroy(vscat_surf_to_subsurf, ierr)

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if(associated(simulation%rt_process_model_coupler)) &
    simulation%rt_process_model_coupler%sim_aux => simulation%sim_aux
  if(option%nsurfflowdof>0 .and. &
     associated(surf_simulation%surf_flow_process_model_coupler)) &
    surf_simulation%surf_flow_process_model_coupler%sim_aux => simulation%sim_aux

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  if (associated(cur_process_model_coupler%next)) then
    cur_process_model_coupler => cur_process_model_coupler%next
    call cur_process_model_coupler%SetAuxData()
  endif

  deallocate(simulation_old)

end subroutine SurfSubsurfaceInitializePostPETSc


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
subroutine SurfSubsurfCreateSurfSubSurfVScats(realization, surf_realization, &
                                              surf_to_subsurf, subsurf_to_surf)

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

  type(realization_type),pointer         :: realization
  type(surface_realization_type),pointer :: surf_realization
  VecScatter                             :: surf_to_subsurf
  VecScatter                             :: subsurf_to_surf

  type(option_type),pointer            :: option
  type(unstructured_grid_type),pointer :: subsurf_grid
  type(unstructured_grid_type),pointer :: surf_grid
  type(level_type),pointer             :: cur_level
  type(patch_type),pointer             :: cur_patch
  type(region_type),pointer            :: cur_region,top_region
  type(region_type),pointer            :: patch_region

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

  call SurfSubsurfCreateSurfSubSurfVScat(realization,surf_realization,prod, &
                                     surf_petsc_ids,surf_to_subsurf)
  call MatDestroy(prod,ierr)

  call MatMatMult(Mat_vert_to_face_surf,Mat_vert_to_face_subsurf_transp, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,prod,ierr)

#if UGRID_DEBUG
  string = 'prod_2.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(prod,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  call SurfSubsurfCreateSurfSubSurfVScat(realization,surf_realization,prod, &
                                        subsurf_petsc_ids,subsurf_to_surf)

  call MatDestroy(prod,ierr)

  call MatDestroy(Mat_vert_to_face_subsurf,ierr)
  call MatDestroy(Mat_vert_to_face_subsurf_transp,ierr)
  call MatDestroy(Mat_vert_to_face_surf,ierr)
  call MatDestroy(Mat_vert_to_face_surf_transp,ierr)

  call VecDestroy(subsurf_petsc_ids,ierr)
  call VecDestroy(surf_petsc_ids,ierr)

end subroutine SurfSubsurfCreateSurfSubSurfVScats

! ************************************************************************** !
!> This subroutine creates a single vector scatter context
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
subroutine SurfSubsurfCreateSurfSubSurfVScat( &
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

end subroutine SurfSubsurfCreateSurfSubSurfVScat

! ************************************************************************** !
!> This routine creates subsurface vectors.
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
subroutine SurfSubsurfCreateSubsurfVecs(subsurf_realization, option, &
                                        subsurf_pres, subsurf_pres_top_bc)

  use Realization_class
  use Coupler_module
  use Option_module
  use String_module

  implicit none

  type(realization_type),pointer :: subsurf_realization
  type(option_type),pointer :: option
  Vec :: subsurf_pres
  Vec :: subsurf_pres_top_bc

  type(coupler_list_type),pointer :: coupler_list
  type(coupler_type),pointer :: coupler

  PetscInt :: num_conn
  PetscInt :: found
  PetscInt :: found_global
  PetscErrorCode :: ierr

  call VecCreate(option%mycomm,subsurf_pres,ierr)
  call VecSetSizes(subsurf_pres, &
                   subsurf_realization%discretization%grid%nlmax, &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(subsurf_pres,ierr)
  call VecSet(subsurf_pres,0.d0,ierr)

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
    call VecCreate(option%mycomm,subsurf_pres_top_bc,ierr)
    call VecSetSizes(subsurf_pres_top_bc,num_conn,PETSC_DECIDE,ierr)
    call VecSetFromOptions(subsurf_pres_top_bc,ierr)
    call VecSet(subsurf_pres_top_bc,0.d0,ierr)
  endif

end subroutine SurfSubsurfCreateSubsurfVecs

! ************************************************************************** !
!> This routine creates surface vectors.
!!
!> @author
!! Gautam Bisht,LBNL
!!
!! date: 08/20/13
! ************************************************************************** !
subroutine SurfSubsurfCreateSurfVecs(surf_realization,option,surf_head)

  use Surface_Realization_class
  use Option_module

  implicit none

  type(surface_realization_type),pointer :: surf_realization
  type(option_type),pointer :: option
  Vec :: surf_head

  PetscErrorCode :: ierr

  call VecCreate(option%mycomm,surf_head,ierr)
  call VecSetSizes(surf_head,surf_realization%discretization%grid%nlmax, &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(surf_head,ierr)
  call VecSet(surf_head,0.d0,ierr)

end subroutine SurfSubsurfCreateSurfVecs

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/01/13
! ************************************************************************** !
subroutine SurfSubsurfInitCommandLineSettings(option)

  use Option_module
  use Input_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = MULTISIMULATION_SIM_TYPE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = STOCHASTIC_SIM_TYPE
  endif
  
end subroutine SurfSubsurfInitCommandLineSettings

end module Surf_Subsurf_Factory_module

#endif
