module Reactive_Transport_module

  use Transport_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"
  
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petsclog.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  logical :: aux_vars_up_to_date = .false.
  logical :: inactive_cells_exist = .false.
  
  public :: RTTimeCut, RTSetup, RTMaxChange, RTUpdateSolution, RTResidual, &
            RTJacobian, RTInitializeTimestep, RTGetTecplotHeader, &
            RTGetVarFromArray
  
contains

! ************************************************************************** !
!
! RTTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTTimeCut(realization)
 
  use Realization_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr

  field => realization%field
 
  call VecCopy(field%tran_xx,field%tran_yy,ierr)
 
end subroutine RTTimeCut

! ************************************************************************** !
!
! RTSetup: 
! author: Glenn Hammond
! date: 00/00/00
!
! ************************************************************************** !
subroutine RTSetup(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RTSetup

! ************************************************************************** !
!
! RTSetupPatch: 
! author: Glenn Hammond
! date: 00/00/00
!
! ************************************************************************** !
subroutine RTSetupPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Condition_module
  use Connection_module
 
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%RTAux => RTAuxCreate()
    
  ! allocate aux_var data structures for all grid cells
  allocate(patch%RTAux%aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RTAuxVarInit(patch%RTAux%aux_vars(ghosted_id),option)
  enddo
  patch%RTAux%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! aux_var data structures for them
  boundary_condition => patch%transport_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection%num_connections
    boundary_condition => boundary_condition%next
  enddo
  allocate(patch%RTAux%aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call RTAuxVarInit(patch%RTAux%aux_vars_bc(iconn),option)
  enddo
  patch%RTAux%num_aux_bc = sum_connection

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call RTCreateZeroArray(patch,option)

end subroutine RTSetupPatch

! ************************************************************************** !
!
! RTInitializeTimestep: 
! author: Glenn Hammond
! date: 00/00/00
!
! ************************************************************************** !
subroutine RTInitializeTimestep(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTInitializeTimestepPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RTInitializeTimestep

! ************************************************************************** !
!
! RTInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RTInitializeTimestepPatch(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  call RTUpdateFixedAccumulationPatch(realization)

end subroutine RTInitializeTimestepPatch
  
! ************************************************************************** !
!
! RTUpdateSolution: Updates data in module after a successful time step
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
subroutine RTUpdateSolution(realization)

  use Realization_module
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%tran_xx,realization%field%tran_yy,ierr)

end subroutine RTUpdateSolution

! ************************************************************************** !
!
! RTUpdateFixedAccumulationPatch: Computes derivative of accumulation term in 
!                           residual function 
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateFixedAccumulationPatch(realization)

  use Realization_module
  use Reactive_Transport_Aux_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: xx_p(:), porosity_loc_p(:), saturation_loc_p(:), &
                        tor_loc_p(:), volume_p(:), accum_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  aux_vars => patch%RTAux%aux_vars
  grid => patch%grid

  call VecGetArrayF90(field%tran_xx,xx_p, ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%saturation_loc,saturation_loc_p,ierr)
  call VecGetArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecGetArrayF90(grid%volume,volume_p,ierr)

  call VecGetArrayF90(field%tran_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%ncomp
    istart = iend-option%ncomp+1
    call RTAuxVarCompute(xx_p(istart:iend),aux_vars(ghosted_id),option)
    call RTAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                        saturation_loc_p(local_id), volume_p(local_id), &
                        option,accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayF90(field%tran_xx,xx_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayF90(field%saturation_loc,saturation_loc_p,ierr)
  call VecRestoreArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecRestoreArrayF90(grid%volume,volume_p,ierr)

  call VecRestoreArrayF90(field%tran_accum, accum_p, ierr)

#if 0
  call RTNumericalJacobianTest(field%tran_xx,realization)
#endif

end subroutine RTUpdateFixedAccumulationPatch

! ************************************************************************** !
!
! RTNumericalJacobianTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RTNumericalJacobianTest(xx,realization)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_type) :: realization

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  PetscReal :: derivative, perturbation
  
  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  
  PetscInt :: idof, idof2, icell

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  call VecDuplicate(xx,xx_pert,ierr)
  call VecDuplicate(xx,res,ierr)
  call VecDuplicate(xx,res_pert,ierr)
  
  call MatCreate(PETSC_COMM_WORLD,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, &
                   grid%nlmax*option%ntrandof, &
                   grid%nlmax*option%ntrandof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call RTResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr)
  do icell = 1,grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    endif
    do idof = (icell-1)*option%ntrandof+1,icell*option%ntrandof 
      call veccopy(xx,xx_pert,ierr)
      call vecgetarrayf90(xx_pert,vec_p,ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call vecrestorearrayf90(xx_pert,vec_p,ierr)
      call RTResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
      call vecgetarrayf90(res_pert,vec_p,ierr)
      do idof2 = 1, grid%nlmax*option%ntrandof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values,ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr)
    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'RTnumerical_jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call MatDestroy(A,ierr)
  
  call VecDestroy(xx_pert,ierr)
  call VecDestroy(res,ierr)
  call VecDestroy(res_pert,ierr)
  
end subroutine RTNumericalJacobianTest

! ************************************************************************** !
!
! RTAccumulationDerivative: Computes derivative of accumulation term in 
!                           residual function 
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTAccumulationDerivative(aux_var,por,sat,vol,option,Res)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  PetscReal :: por, sat, vol
  type(option_type) :: option
  PetscReal :: Res(option%ncomp)
  
  PetscInt :: icomp
  PetscReal :: psv_t
  
  psv_t = por*sat*vol/option%dt
  do icomp=1,option%ncomp
    Res(icomp) = psv_t
  enddo

end subroutine RTAccumulationDerivative

! ************************************************************************** !
!
! RTAccumulation: Computes accumulation term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTAccumulation(aux_var,por,sat,vol,option,Res)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  PetscReal :: por, sat, vol
  type(option_type) :: option
  PetscReal :: Res(option%ncomp)
  
  PetscInt :: icomp
  PetscReal :: psv_t
  
  psv_t = por*sat*vol/option%dt
  do icomp=1,option%ncomp
    Res(icomp) = psv_t*aux_var%total(icomp) 
  enddo

end subroutine RTAccumulation

! ************************************************************************** !
!
! RTResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RTResidual(snes,xx,r,realization,ierr)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  field => realization%field
  grid => realization%patch%grid
  
  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
  call GridGlobalToLocal(grid,xx,field%tran_xx_loc,NTRANDOF)
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_loc,ONEDOF)
  call GridLocalToLocal(grid,field%icap_loc,field%icap_loc,ONEDOF)

  call GridLocalToLocal(grid,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call GridLocalToLocal(grid,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTResidualPatch(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RTResidual

! ************************************************************************** !
!
! RTResidualPatch: Computes residual function for reactive transport
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTResidualPatch(snes,xx,r,realization,ierr)

  use Realization_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: xx_loc_p(:), r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), saturation_loc_p(:), tor_loc_p(:), &
                        volume_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: i, istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  PetscReal :: Res(realization%option%ncomp)
  PetscViewer :: viewer
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  aux_vars => patch%RTAux%aux_vars
  aux_vars_bc => patch%RTAux%aux_vars_bc
  
  ! Communication -----------------------------------------
  ! These 3 must be called before RTUpdateAuxVars()
  call GridGlobalToLocal(grid,xx,field%tran_xx_loc,NTRANDOF)

  call RTUpdateAuxVars(realization)
  aux_vars_up_to_date = .false. ! override flags since they will soon be out of date  

  ! Get pointer to Vector data
  call VecGetArrayF90(field%tran_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%tran_accum, accum_p, ierr)
 
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%saturation_loc, saturation_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)

  r_p = 0.d0
#if 1
  ! Accumulation terms ------------------------------------
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%ncomp
    istart = iend-option%ncomp+1
    call RTAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                        saturation_loc_p(ghosted_id),volume_p(local_id), &
                        option,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%ncomp)
  enddo
#endif

  ! Source/Sink terms -------------------------------------
  
#if 1
  ! Interior Flux Terms -----------------------------------
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      dist_up = distance*fraction_upwind
      dist_dn = distance-dist_up ! should avoid truncation error

      call TFlux(aux_vars(ghosted_id_up),porosity_loc_p(ghosted_id_up), &
                 tor_loc_p(ghosted_id_up),saturation_loc_p(ghosted_id_up), &
                 dist_up, &
                 aux_vars(ghosted_id_dn),porosity_loc_p(ghosted_id_dn), &
                 tor_loc_p(ghosted_id_dn),saturation_loc_p(ghosted_id_dn), &
                 dist_up, &
                 cur_connection_set%area(iconn),option, &
                 patch%internal_velocities(:,iconn),Res)

      if (local_id_up>0) then
        iend = local_id_up*option%ncomp
        istart = iend-option%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%ncomp)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*option%ncomp
        istart = iend-option%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%ncomp)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%transport_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      call TBCFlux(boundary_condition%condition%itype(1), &
                   aux_vars_bc(sum_connection), &
                   aux_vars(ghosted_id), &
                   porosity_loc_p(ghosted_id), &
                   tor_loc_p(ghosted_id), &
                   saturation_loc_p(ghosted_id), &
                   cur_connection_set%dist(0,iconn), &
                   cur_connection_set%area(iconn), &
                   option,patch%boundary_velocities(:,sum_connection),Res)
 
      iend = local_id*option%ncomp
      istart = iend-option%ncomp+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%ncomp)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif  

  if (inactive_cells_exist) then
    do i=1,patch%RTAux%n_zero_rows
      r_p(patch%RTAux%zero_rows_local(i)) = 0.d0
    enddo
  endif
  
  ! Restore vectors
  call VecRestoreArrayF90(field%tran_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%tran_accum, accum_p, ierr)
 
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%saturation_loc, saturation_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'RTresidual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'RTxx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
end subroutine RTResidualPatch

! ************************************************************************** !
!
! RTJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RTJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTJacobianPatch(snes,xx,A,B,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RTJacobian

! ************************************************************************** !
!
! RTJacobianPatch: Computes Jacobian for reactive transport
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTJacobianPatch(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag  
  type(realization_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: xx_loc_p(:), r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), saturation_loc_p(:), tor_loc_p(:), &
                        volume_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
    
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  PetscReal :: Jup(realization%option%ncomp,realization%option%ncomp)
  PetscReal :: Jdn(realization%option%ncomp,realization%option%ncomp)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn, rdum
  PetscViewer :: viewer
  
  option => realization%option
  field => realization%field
  patch => realization%patch  
  grid => patch%grid
  aux_vars => patch%RTAux%aux_vars
  aux_vars_bc => patch%RTAux%aux_vars_bc

  flag = SAME_NONZERO_PATTERN  
  call MatZeroEntries(A,ierr)

  ! Get pointer to Vector data
  call VecGetArrayF90(field%tran_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%tran_accum, accum_p, ierr)
 
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%saturation_loc, saturation_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
    
#if 1  
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%ncomp
    istart = iend-option%ncomp+1
    call RTAccumulationDerivative(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                                  saturation_loc_p(ghosted_id),volume_p(local_id), &
                                  option,Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)                        
  enddo
#endif

  ! Source/Sink terms -------------------------------------
  
#if 1
  ! Interior Flux Terms -----------------------------------
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      dist_up = distance*fraction_upwind
      dist_dn = distance-dist_up ! should avoid truncation error

      call TFluxDerivative(aux_vars(ghosted_id_up),porosity_loc_p(ghosted_id_up), &
                 tor_loc_p(ghosted_id_up),saturation_loc_p(ghosted_id_up), &
                 dist_up, &
                 aux_vars(ghosted_id_dn),porosity_loc_p(ghosted_id_dn), &
                 tor_loc_p(ghosted_id_dn),saturation_loc_p(ghosted_id_dn), &
                 dist_up, &
                 cur_connection_set%area(iconn),option, &
                 patch%internal_velocities(:,iconn),Jup,Jdn)

      if (local_id_up>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)        
      endif
   
      if (local_id_dn>0) then
        Jup = -1.d0*Jup
        Jdn = -1.d0*Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%transport_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      call TBCFluxDerivative(boundary_condition%condition%itype(1), &
                   aux_vars_bc(sum_connection), &
                   aux_vars(ghosted_id), &
                   porosity_loc_p(ghosted_id), &
                   tor_loc_p(ghosted_id), &
                   saturation_loc_p(ghosted_id), &
                   cur_connection_set%dist(0,iconn), &
                   cur_connection_set%area(iconn), &
                   option,patch%boundary_velocities(:,sum_connection),Jdn)
 
      Jdn = -1.d0*Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif  

  ! Restore vectors
  call VecRestoreArrayF90(field%tran_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%tran_accum, accum_p, ierr)
 
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%saturation_loc, saturation_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  
  if (inactive_cells_exist) then
    rdum = 1.d0
    call MatZeroRowsLocal(A,patch%RTAux%n_zero_rows, &
                          patch%RTAux%zero_rows_local_ghosted,rdum,ierr) 
  endif

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'RTjacobian.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
    
end subroutine RTJacobianPatch

! ************************************************************************** !
!
! RTUpdateAuxVars: Updates the auxilliary variables associated with 
!                  reactive transport
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateAuxVars(realization)

  use Realization_module
  use Reactive_Transport_Aux_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  
  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%ncomp)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  
  call VecGetArrayF90(field%tran_xx_loc,xx_loc_p, ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%ncomp
    istart = iend-option%ncomp+1
   
    call RTAuxVarCompute(xx_loc_p(istart:iend), &
                         patch%RTAux%aux_vars(ghosted_id), &
                         option)
  enddo

  boundary_condition => patch%transport_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      do idof=1,option%ncomp
        select case(boundary_condition%condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,NEUMANN_BC)
            xxbc(idof) = boundary_condition%aux_real_var(idof,iconn)
          case(ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%ncomp+idof)
        end select
      enddo
      
      call RTAuxVarCompute(xxbc, &
                           patch%RTAux%aux_vars_bc(sum_connection), &
                           option)
    enddo
    boundary_condition => boundary_condition%next
  enddo


  call VecRestoreArrayF90(field%tran_xx_loc,xx_loc_p, ierr)
  
  aux_vars_up_to_date = .true.

end subroutine RTUpdateAuxVars

! ************************************************************************** !
!
! RTCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RTCreateZeroArray(patch,option)

  use Patch_module
  use Grid_module
  use Option_module
  
  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id, icomp

  type(grid_type), pointer :: grid
  PetscInt :: flag = 0
  PetscInt :: n_zero_rows
  PetscInt, pointer :: zero_rows_local(:)
  PetscInt, pointer :: zero_rows_local_ghosted(:)
  PetscErrorCode :: ierr

  grid => patch%grid
  
  n_zero_rows = 0

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
        n_zero_rows = n_zero_rows + option%ncomp
      else
      endif
    enddo
  else
  endif

  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))

  zero_rows_local = 0
  zero_rows_local_ghosted = 0
  ncount = 0

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
        do icomp = 1, option%ncomp
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*option%ncomp+icomp
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%ncomp+icomp-1
        enddo
      else
      endif
    enddo
  endif

  patch%RTAux%zero_rows_local => zero_rows_local
  patch%RTAux%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%RTAux%n_zero_rows = n_zero_rows  

  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                     PETSC_COMM_WORLD,ierr)
  if (flag > 0) inactive_cells_exist = .true.

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif
  
end subroutine RTCreateZeroArray

! ************************************************************************** !
!
! RTGetVarFromArray: Extracts variables indexed by ivar and isubvar
!                          from RT type
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine RTGetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none

  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field

  if (.not.aux_vars_up_to_date) call RTUpdateAuxVars(realization)

  select case(ivar)
    case(FREE_ION_CONCENTRATION)
      call VecStrideGather(field%tran_xx,isubvar,vec,INSERT_VALUES,ierr)
    case(MATERIAL_ID,TOTAL_CONCENTRATION)
      call VecGetArrayF90(vec,vec_ptr,ierr)
      select case(ivar)
        case(TOTAL_CONCENTRATION)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)    
            vec_ptr(local_id) = patch%RTAux%aux_vars(ghosted_id)%total(isubvar)
          enddo
        case(MATERIAL_ID)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
          enddo
      end select
      call VecRestoreArrayF90(vec,vec_ptr,ierr)
  end select
  
end subroutine RTGetVarFromArray

! ************************************************************************** !
!
! RTGetVarFromArrayAtCell: Returns variables indexed by ivar,
!                          isubvar, local id from Reactive Transport type
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function RTGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none

  PetscReal :: RTGetVarFromArrayAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: local_id

  PetscReal :: value
  PetscInt :: ghosted_id
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field

  if (.not.aux_vars_up_to_date) call RTUpdateAuxVars(realization)

  select case(ivar)
    case(FREE_ION_CONCENTRATION)
      call VecGetArrayF90(field%tran_xx,vec_ptr,ierr)
      value = vec_ptr((local_id-1)*option%ntrandof+isubvar)
      call VecRestoreArrayF90(field%tran_xx,vec_ptr,ierr)
    case(TOTAL_CONCENTRATION)
      ghosted_id = grid%nL2G(local_id)    
      value = patch%RTAux%aux_vars(ghosted_id)%total(isubvar)
  end select
  
  RTGetVarFromArrayAtCell = value
  
end function RTGetVarFromArrayAtCell

! ************************************************************************** !
!
! RTMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTMaxChange(realization)

  use Realization_module
  use Option_module
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field

  option%dcmax=0.D0
  
  call VecWAXPY(field%tran_dxx,-1.d0,field%tran_xx,field%tran_yy,ierr)
  call VecStrideNormAll(field%tran_dxx,ZERO_INTEGER,NORM_INFINITY,option%dcmax,ierr)
    
end subroutine RTMaxChange

! ************************************************************************** !
!
! RTGetTecplotHeader: Returns Reactive Transport contribution to 
!                     Tecplot file header
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
function RTGetTecplotHeader(realization)

  use Realization_module
  use Option_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: RTGetTecplotHeader
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  PetscInt :: i
  
  option => realization%option
  
  string = '' 
  do i=1,option%ntrandof
    write(string2,'('',"COMP('',i2,'')"'')') i
    string = trim(string) // trim(string2)
  enddo
  
  RTGetTecplotHeader = string

end function RTGetTecplotHeader

end module Reactive_Transport_module
