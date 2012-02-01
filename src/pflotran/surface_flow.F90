#ifdef SURFACE_FLOW 

module Surface_Flow_module

  use Global_Aux_module
  
  implicit none
  
  private
  
#include "definitions.h"  

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"


  public SurfaceFlowSetup, &
         SurfaceFlowResidual, &
         SurfaceFlowJacobian
  
contains

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowSetup(realization)

  use Realization_module
  
  type(realization_type) :: realization
  
end subroutine SurfaceFlowSetup


! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowResidual(snes,xx,r,realization,ierr)
!subroutine SurfaceFlowResidual(realization)

  use Realization_module
  use Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module
  use Logging_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option
  
  call PetscLogEventBegin(logging%event_r_residual,ierr)
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option

  ! pass #1 for internal and boundary flux terms
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%surf_patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call SurfaceFlowResidualPatch1(snes,xx,r,realization,ierr)
      !call SurfaceFlowResidualPatch1(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceFlowResidual

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowResidualPatch1(snes,xx,r,realization,ierr)
!subroutine SurfaceFlowResidualPatch1(realization)

  use water_eos_module

  use Connection_module
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:)
  PetscInt :: icap_up, icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  PetscViewer :: viewer
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  !richards_parameter => patch%aux%Richards%richards_parameter
  !rich_aux_vars => patch%aux%Richards%aux_vars
  !rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  !global_aux_vars => patch%aux%Global%aux_vars
  !global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  call GridVecGetArrayF90(grid,r, r_p, ierr)

  r_p = 0.d0
  
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
      
      !write(*,'(5I5)'),sum_connection,ghosted_id_dn,ghosted_id_up,local_id_dn,local_id_up

      Res(1) = 1.0d0

      if (.not.option%use_samr) then
         
        if (local_id_up>0) then
          r_p(local_id_up) = r_p(local_id_up) + Res(1)
        endif
         
        if (local_id_dn>0) then
          r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
        endif
      endif
      
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      !write(*,'(3I5)'),sum_connection,ghosted_id,local_id
  
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)

  !realization%option%io_buffer = 'stopping '
  !call printErrMsg(realization%option)

  call PetscViewerASCIIOpen(option%mycomm,'r_surfflow.out',viewer,ierr)
  call VecView(r,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

end subroutine SurfaceFlowResidualPatch1  

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Logging_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm


  call PetscLogEventBegin(logging%event_r_jacobian,ierr)

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

  ! pass #1 for internal and boundary flux terms
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%surf_patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch

      call SurfaceFlowJacobianPatch1(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceFlowJacobian


! ************************************************************************** !
!
!
! ************************************************************************** !

subroutine SurfaceFlowJacobianPatch1(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module

  use Connection_module
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscReal, pointer :: porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscInt :: icap_up,icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  !type(field_type), pointer :: field 
  !type(richards_parameter_type), pointer :: richards_parameter
  !type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:) 
  !type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  !field => realization%field
  !richards_parameter => patch%aux%Richards%richards_parameter
  !rich_aux_vars => patch%aux%Richards%aux_vars
  !rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  !global_aux_vars => patch%aux%Global%aux_vars
  !global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  !call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  !call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  !call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  !call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)

  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    
    !write(*,*),'cur_connection_set%num_connections : ',cur_connection_set%num_connections
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      
      !write(*,*),'ghosted_ids :',ghosted_id_up, ghosted_id_dn

      !if (patch%imat(ghosted_id_up) <= 0 .or. &
      !    patch%imat(ghosted_id_dn) <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      !write(*,*),'local_ids   :',ghosted_id_up, ghosted_id_dn,local_id_up,local_id_dn

      Jup = 1.0d0
      Jdn = 1.0d0

      if (local_id_up > 0) then
          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                        Jup,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                        Jdn,ADD_VALUES,ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                              Jdn,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               Jup,ADD_VALUES,ierr)
      endif
      
    enddo
    cur_connection_set => cur_connection_set%next

  enddo
   
  !realization%option%io_buffer = 'stopping '
  !call printErrMsg(realization%option)
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  call PetscViewerASCIIOpen(option%mycomm,'jacobian_surfflow.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)


end subroutine SurfaceFlowJacobianPatch1

end module Surface_Flow_module

#endif
