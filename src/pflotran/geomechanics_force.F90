#ifdef GEOMECH 

module Geomechanics_Force_module

  use Geomechanics_Global_Aux_module
  
  implicit none
  
  private
  
#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
!#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"
#include "finclude/petscts.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: GeomechForceSetup, &
            GeomechForceUpdateAuxVars, &
            GeomechanicsForceInitialGuess, &
            GeomechForceResidual, &
            GeomechForceJacobian
  
contains

! ************************************************************************** !
!
! GeomechForceSetup: Sets up the geomechanics calculations
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechForceSetup(geomech_realization)

  use Geomechanics_Realization_module
  
  type(geomech_realization_type) :: geomech_realization

  call GeomechForceSetPlotVariables(geomech_realization)
  
end subroutine GeomechForceSetup

! ************************************************************************** !
!
! GeomechForceSetPlotVariables: Set up of geomechanics plot variables
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechForceSetPlotVariables(geomech_realization)
  
  use Geomechanics_Realization_module
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(geomech_realization_type) :: geomech_realization
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => geomech_realization%output_option%output_variable_list
  
  name = 'disp_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GEOMECH_DISP_X)
                               
  name = 'disp_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GEOMECH_DISP_Y)
                               
  name = 'disp_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GEOMECH_DISP_Z)

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine GeomechForceSetPlotVariables

! ************************************************************************** !
!
! GeomechanicsForceInitialGuess: Sets up the inital guess for the solution
!                                The boudnary conditions are set here
! author: Satish Karra, LANL
! date: 06/19/13
!
! ************************************************************************** !
subroutine GeomechanicsForceInitialGuess(realization)

  use Geomechanics_Realization_module
  use Geomechanics_Field_module
  use Option_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Patch_module
  use Geomechanics_Coupler_module
  use Geomechanics_Region_module
  
  implicit none
  
  type(geomech_realization_type) :: realization
  
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: field
  type(geomech_patch_type), pointer :: patch
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_grid_type), pointer :: grid
  type(gm_region_type), pointer :: region
  
  PetscInt :: ghosted_id,local_id,total_verts,ivertex
  PetscReal, pointer :: xx_p(:)
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%geomech_field
  patch => realization%geomech_patch
  grid => patch%geomech_grid
  
  call GeomechGridVecGetArrayF90(grid,field%disp_xx,xx_p,ierr)
  
  boundary_condition => patch%geomech_boundary_conditions%first
  total_verts = 0
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      total_verts = total_verts + 1
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif    
      
      ! X displacement 
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_X_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Y displacement 
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Y_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Z displacement      
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Z_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
    enddo
    boundary_condition => boundary_condition%next      
  enddo
  
  call GeomechGridVecRestoreArrayF90(grid,field%disp_xx,xx_p,ierr)

end subroutine GeomechanicsForceInitialGuess

! ************************************************************************** !
!
! GeomechForceUpdateAuxVars: Updates the geomechanics variables
! author: Satish Karra, LANL
! date: 06/18/13
!
! ************************************************************************** !
subroutine GeomechForceUpdateAuxVars(geomech_realization)

  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Option_module
  use Geomechanics_Field_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Material_module
  use Geomechanics_Global_Aux_module
  use Geomechanics_Region_module

  implicit none

  type(geomech_realization_type)            :: geomech_realization
  
  type(option_type), pointer                :: option
  type(geomech_patch_type), pointer         :: patch
  type(geomech_grid_type), pointer          :: grid
  type(geomech_field_type), pointer         :: geomech_field
  type(gm_region_type), pointer             :: region
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)

  PetscInt :: ghosted_id, local_id
  PetscReal, pointer :: xx_loc_p(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  geomech_field => geomech_realization%geomech_field

  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  
  call GeomechGridVecGetArrayF90(grid,geomech_field%disp_xx_loc,xx_loc_p,ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax_node
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_X_DOF) = &
        xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)
  enddo
   
  call GeomechGridVecRestoreArrayF90(grid,geomech_field%disp_xx_loc, &
                                     xx_loc_p,ierr)

end subroutine GeomechForceUpdateAuxVars

! ************************************************************************** !
!
! GeomechForceResidual: Computes the residual equation 
! author: Satish Karra
! date: 06/21/13
!
! ************************************************************************** !
subroutine GeomechForceResidual(snes,xx,r,realization,ierr)

  use Geomechanics_Realization_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Logging_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(geomech_realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_field_type), pointer :: field
  type(option_type), pointer :: option
  
  call PetscLogEventBegin(geomech_logging%event_geomech_residual,ierr)
  
  field => realization%geomech_field
  discretization => realization%discretization
  option => realization%option

  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
  call GeomechDiscretizationGlobalToLocal(discretization,xx, &
                                          field%disp_xx_loc,NGEODOF)
  
  call GeomechForceResidualPatch(snes,xx,r,realization,ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm, &
                              'Geomech_residual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Geomech_xx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  call PetscLogEventEnd(geomech_logging%event_geomech_residual,ierr)

end subroutine GeomechForceResidual

! ************************************************************************** !
!
! GeomechForceResidualPatch: Computes the residual equation on a patch 
! author: Satish Karra
! date: 06/24/13
!
! ************************************************************************** !
subroutine GeomechForceResidualPatch(snes,xx,r,realization,ierr)

  use Geomechanics_Realization_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Logging_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(geomech_realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: res_vec(:)
  PetscInt :: ielem,ivertex 
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: petsc_id, local_id
        
                    
  field => realization%geomech_field
  discretization => realization%discretization
  patch => realization%geomech_patch
  grid => patch%geomech_grid
  option => realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars  

  call GeomechForceUpdateAuxVars(realization)
  ! Add flag for the update
  
  call VecSet(r,0.d0,ierr)
 
  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes),option%ngeomechdof))
    allocate(petsc_ids(size(elenodes)))
    allocate(ids(size(elenodes)*option%ngeomechdof))
    allocate(res_vec(size(elenodes)*option%ngeomechdof))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%EleType
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(ivertex,idof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
    enddo
    call GeomechForceLocalElemResidual(elenodes,local_coordinates, &
       local_disp,eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,res_vec,option)
    call VecSetValues(r,size(ids),ids,res_vec,ADD_VALUES,ierr) 
    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(petsc_ids)
    deallocate(ids)
    deallocate(res_vec)
  enddo
      
  call VecAssemblyBegin(r,ierr)
  call VecAssemblyEnd(r,ierr)  
  
#ifdef GEOMECH_DEBUG
  call PetscViewerASCIIOpen(realization%option%mycomm, &
                            'Geomech_residual_debug_beforeBC.out',viewer,ierr)
  call VecView(r,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  boundary_condition => patch%geomech_boundary_conditions%first
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif    
      
      ! X displacement 
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_X_DOF-1,0.d0,INSERT_VALUES,ierr)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Y displacement 
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Y_DOF-1,0.d0,INSERT_VALUES,ierr)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Z displacement      
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Z_DOF-1,0.d0,INSERT_VALUES,ierr)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
    enddo
    boundary_condition => boundary_condition%next      
  enddo

  call VecAssemblyBegin(r,ierr)
  call VecAssemblyEnd(r,ierr)  
  
#ifdef GEOMECH_DEBUG
  call PetscViewerASCIIOpen(realization%option%mycomm, &
                            'Geomech_residual_debug_afterBC.out',viewer,ierr)
  call VecView(r,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  call PetscViewerBinaryOpen(realization%option%mycomm,'Geomech_residual_debug_afterBC.bin', &
                             FILE_MODE_WRITE,viewer,ierr)
  call VecView(r,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)  
#endif   

end subroutine GeomechForceResidualPatch

! ************************************************************************** !
!
! GeomechForceLocalElemResidual: Computes the residual for a local element
! author: Satish Karra
! date: 06/24/13
!
! ************************************************************************** !
subroutine GeomechForceLocalElemResidual(elenodes,local_coordinates,local_disp, &
                                         eletype,dim,r,w,res_vec,option)
                                         
  use Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module
  
  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  
  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: zeta(:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: res_vec(:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j,d
  PetscReal :: eye_three(THREE_INTEGER)
  PetscInt :: indx(THREE_INTEGER)
  PetscInt :: dim
  PetscReal :: lambda, mu
  PetscInt :: load_type
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: den_rock
  PetscReal, allocatable :: N(:,:)
  PetscReal, pointer :: vecB_transpose(:,:), transpose_vecB_transpose(:,:)
  PetscReal, pointer :: kron_B_eye(:,:)
  PetscReal, pointer :: kron_B_transpose_eye(:,:)
  PetscReal, pointer :: Trans(:,:)
  PetscReal, pointer :: kron_eye_B_transpose(:,:)
  PetscReal, pointer :: kron_N_eye(:,:)
  PetscReal, pointer :: vec_local_disp(:,:)
  PetscReal, allocatable :: force(:), res_vec_mat(:,:)
  
  allocate(zeta(size(r,2)))
  allocate(B(size(elenodes),dim))
  allocate(Kmat(size(elenodes)*option%ngeomechdof, &
                size(elenodes)*option%ngeomechdof))  
  allocate(force(size(elenodes)*option%ngeomechdof))
  allocate(res_vec_mat(size(elenodes)*option%ngeomechdof,1))
  
  res_vec = 0.d0
  res_vec_mat = 0.d0
  Kmat = 0.d0
  force = 0.d0
  len_w = size(w)
  
  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo
  
  call Transposer(option%ngeomechdof,size(elenodes),Trans)
 
  do igpt = 1, len_w
    shapefunction%EleType = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    allocate(N(size(shapefunction%N),ONE_INTEGER))
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!' 
      call printErrMsg(option)        
    endif
    ! Find the inverse of J_map
    ! Set identity matrix
    call ludcmp(J_map,THREE_INTEGER,indx,d)
    do i = 1, THREE_INTEGER
      eye_three = 0.d0
      eye_three(i) = 1.d0
      call lubksb(J_map,THREE_INTEGER,indx,eye_three)
      inv_J_map(:,i) = eye_three
    enddo
    B = matmul(shapefunction%DN,inv_J_map)
    call GeomechGetLambdaMu(lambda,mu,x)
    call GeomechGetBodyForce(load_type,lambda,mu,x,bf) 
    call ConvertMatrixToVector(transpose(B),vecB_transpose)
    Kmat = Kmat + w(igpt)*lambda* &
      matmul(vecB_transpose,transpose(vecB_transpose))*detJ_map
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    N(:,1)= shapefunction%N    
    call Kron(N,identity,kron_N_eye)
    Kmat = Kmat + w(igpt)*mu*matmul(kron_B_eye,kron_B_transpose_eye)*detJ_map
    Kmat = Kmat + w(igpt)*mu*matmul(matmul(kron_B_eye,kron_eye_B_transpose),Trans)*detJ_map
    force = force + w(igpt)*matmul(kron_N_eye,bf)*detJ_map
    call ShapeFunctionDestroy(shapefunction)
    deallocate(N)
  enddo
  
  call ConvertMatrixToVector(transpose(local_disp),vec_local_disp)
  res_vec_mat = matmul(Kmat,vec_local_disp)
  res_vec = res_vec + res_vec_mat(:,1)
  res_vec = res_vec - force

  deallocate(zeta)
  deallocate(B)
  deallocate(force)
  deallocate(res_vec_mat)

end subroutine GeomechForceLocalElemResidual

! ************************************************************************** !
!
! GeomechForceLocalElemJacobian: Computes the Jacobian for a local element
! author: Satish Karra
! date: 06/24/13
!
! ************************************************************************** !
subroutine GeomechForceLocalElemJacobian(elenodes,local_coordinates,local_disp, &
                                         eletype,dim,r,w,Kmat,option)
                                         
  use Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module
  
  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  
  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: zeta(:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: local_disp(:)
  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j,d
  PetscReal :: eye_three(THREE_INTEGER)
  PetscInt :: indx(THREE_INTEGER)
  PetscInt :: dim
  PetscReal :: lambda, mu
  PetscInt :: load_type
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: den_rock
  PetscReal, allocatable :: N(:,:)
  PetscReal, pointer :: vecB_transpose(:,:), transpose_vecB_transpose(:,:)
  PetscReal, pointer :: kron_B_eye(:,:)
  PetscReal, pointer :: kron_B_transpose_eye(:,:)
  PetscReal, pointer :: Trans(:,:)
  PetscReal, pointer :: kron_eye_B_transpose(:,:)
  PetscReal, pointer :: kron_N_eye(:,:)
  
  allocate(zeta(size(r,2)))
  allocate(B(size(elenodes),dim))
  
  Kmat = 0.d0
  len_w = size(w)
  
  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo
  
  call Transposer(option%ngeomechdof,size(elenodes),Trans)

  do igpt = 1, len_w
    shapefunction%EleType = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    allocate(N(size(shapefunction%N),ONE_INTEGER))
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!' 
      call printErrMsg(option)        
    endif
    ! Find the inverse of J_map
    ! Set identity matrix
    call ludcmp(J_map,THREE_INTEGER,indx,d)
    do i = 1, THREE_INTEGER
      eye_three = 0.d0
      eye_three(i) = 1.d0
      call lubksb(J_map,THREE_INTEGER,indx,eye_three)
      inv_J_map(:,i) = eye_three
    enddo
    B = matmul(shapefunction%DN,inv_J_map)
    call GeomechGetLambdaMu(lambda,mu,x)
    call GeomechGetBodyForce(load_type,lambda,mu,x,bf) 
    call ConvertMatrixToVector(transpose(B),vecB_transpose)
    Kmat = Kmat + w(igpt)*lambda* &
      matmul(vecB_transpose,transpose(vecB_transpose))*detJ_map
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    N(:,1)= shapefunction%N    
    call Kron(N,identity,kron_N_eye)
    Kmat = Kmat + w(igpt)*mu*matmul(kron_B_eye,kron_B_transpose_eye)*detJ_map
    Kmat = Kmat + w(igpt)*mu*matmul(matmul(kron_B_eye,kron_eye_B_transpose),Trans)*detJ_map
    call ShapeFunctionDestroy(shapefunction)
    deallocate(N)
  enddo
  
  deallocate(zeta)
  deallocate(B)

end subroutine GeomechForceLocalElemJacobian

! ************************************************************************** !
!
! GeomechGetLambdaMu: Gets the material properties given the position 
! of the point
! author: Satish Karra
! date: 06/24/13
!
! ************************************************************************** !
subroutine GeomechGetLambdaMu(lambda,mu,coord)

  PetscReal :: lambda, mu
  PetscReal :: E, nu
  PetscReal :: coord(THREE_INTEGER)
 
  E = 1.d4
  nu = 0.0
  
  
  ! This subroutine needs major changes. For given position, it needs to give 
  ! out lambda, mu
  
  lambda = E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
  mu = E/2.d0/(1.d0+nu)


end subroutine GeomechGetLambdaMu

! ************************************************************************** !
!
! GeomechGetBodyForce: Gets the body force at a given position
! of the point
! author: Satish Karra
! date: 06/24/13
!
! ************************************************************************** !
subroutine GeomechGetBodyForce(load_type,lambda,mu,coord,bf)

  PetscInt :: load_type
  PetscReal :: lambda, mu, den_rock
  PetscReal :: coord(THREE_INTEGER)
  PetscReal :: bf(THREE_INTEGER)
 
    
  ! This subroutine needs major changes. For given position, it needs to give 
  ! out lambda, mu, also need to add density of rock
  
  select case(load_type)
    case default
      bf(GEOMECH_DISP_Z_DOF) = -9.81
  end select

end subroutine GeomechGetBodyForce

! ************************************************************************** !
!
! GeomechForceJacobian: Computes the Jacobian
! author: Satish Karra
! date: 06/21/13
!
! ************************************************************************** !
subroutine GeomechForceJacobian(snes,xx,A,B,flag,realization,ierr)

  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Logging_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(geomech_realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(geomech_grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm
  
  call PetscLogEventBegin(geomech_logging%event_geomech_jacobian,ierr)

  option => realization%option

  flag = SAME_NONZERO_PATTERN
  call MatGetType(A,mat_type,ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr)

  call GeomechForceJacobianPatch(snes,xx,J,J,flag,realization,ierr)

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Geomech_jacobian.out', &
                              viewer,ierr)
   
    call MatView(J,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_FROBENIUS,norm,ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_INFINITY,norm,ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option) 
  endif

#if 0
    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_dxx.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_dxx,viewer,ierr) 

    call PetscViewerDestroy(viewer,ierr)
 

    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_yy.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_yy,viewer,ierr) 
    call PetscViewerDestroy(viewer,ierr)
#endif

  call PetscLogEventEnd(geomech_logging%event_geomech_jacobian,ierr)
!  call printErrMsg(option)

end subroutine GeomechForceJacobian

! ************************************************************************** !
!
! GeomechForceJacobianPatch: Computes the Jacobian on a patch
! author: Satish Karra
! date: 06/21/13
!
! ************************************************************************** !
subroutine GeomechForceJacobianPatch(snes,xx,A,B,flag,realization,ierr)
       
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Field_module
  use Geomechanics_Debug_module
  use Geomechanics_Discretization_module
  use Option_module
  use Unstructured_Cell_module
  use Geomechanics_Region_module
      
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  MatStructure flag
  PetscViewer :: viewer

  PetscErrorCode :: ierr
   
  type(geomech_realization_type) :: realization 
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:)
  PetscInt, allocatable :: ids(:)
  PetscInt, allocatable :: ghosted_ids(:)
  PetscReal, allocatable :: Jac(:,:)
  PetscInt :: ielem,ivertex 
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: local_id, petsc_id
        
  field => realization%geomech_field
  discretization => realization%discretization
  patch => realization%geomech_patch
  grid => patch%geomech_grid
  option => realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars  

 
  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes)*option%ngeomechdof))
    allocate(ghosted_ids(size(elenodes)))
    allocate(ids(size(elenodes)*option%ngeomechdof))
    allocate(Jac(size(elenodes)*option%ngeomechdof,size(elenodes)*option%ngeomechdof))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%EleType
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      ghosted_ids(ivertex) = ghosted_id
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(idof + (ivertex-1)*option%ngeomechdof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (ghosted_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
    enddo
    call GeomechForceLocalElemJacobian(elenodes,local_coordinates, &
       local_disp,eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,Jac,option)
    call MatSetValuesLocal(A,size(ids),ids,size(ids),ids,Jac,ADD_VALUES,ierr) 
    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(ghosted_ids)
    deallocate(ids)
    deallocate(Jac)
  enddo
  
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)  
  
  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  boundary_condition => patch%geomech_boundary_conditions%first
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif    
      
      ! X displacement 
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            call MatZeroRows(A,1, &
              (petsc_id-1)*option%ngeomechdof + GEOMECH_DISP_X_DOF-1, &
              1.d0,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Y displacement 
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            call MatZeroRows(A,1, &
              (petsc_id-1)*option%ngeomechdof + GEOMECH_DISP_Y_DOF-1, &
              1.d0,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Z displacement      
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            call MatZeroRows(A,1, &
              (petsc_id-1)*option%ngeomechdof + GEOMECH_DISP_Z_DOF-1, &
              1.d0,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
    enddo
    boundary_condition => boundary_condition%next      
  enddo
   
#ifdef GEOMECH_DEBUG  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Geomech_jacobian_afterBC.out', &
                              viewer,ierr)
   
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)  
    call PetscViewerBinaryOpen(realization%option%mycomm,'Geomech_jacobian_afterBC.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)  
#endif    
    
end subroutine GeomechForceJacobianPatch  

end module Geomechanics_Force_module

#endif
