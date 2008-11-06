module Reactive_Transport_module

  use Transport_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  
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
  
  public :: RTTimeCut, RTSetup, RTMaxChange, RTUpdateSolution, RTResidual, &
            RTJacobian, RTInitializeTimestep, RTGetTecplotHeader, &
            RTUpdateAuxVars, RTUpdateDenAndSat
  
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
  call RTInitializeTimestep(realization)  
  ! set densities and weights to t+dt
  call RTUpdateDenAndSat(realization,realization%option%tran_weight_t1)
 
end subroutine RTTimeCut

! ************************************************************************** !
!
! RTSetup: 
! author: Glenn Hammond
! date: 02/22/08
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
! date: 02/22/08
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
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction

  patch%aux%RT => RTAuxCreate()
    
  ! allocate aux_var data structures for all grid cells
  allocate(patch%aux%RT%aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RTAuxVarInit(patch%aux%RT%aux_vars(ghosted_id),reaction,option)
  enddo
  patch%aux%RT%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! aux_var data structures for them
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo
  allocate(patch%aux%RT%aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call RTAuxVarInit(patch%aux%RT%aux_vars_bc(iconn),reaction,option)
  enddo
  patch%aux%RT%num_aux_bc = sum_connection

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call RTCreateZeroArray(patch,reaction,option)
  
end subroutine RTSetupPatch

! ************************************************************************** !
!
! RTInitializeTimestep: 
! author: Glenn Hammond
! date: 02/22/08
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

  ! set densities and weights to t, as opposed to t+dt
  call RTUpdateDenAndSatPatch(realization,realization%option%tran_weight_t0)
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
  use Level_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization

  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%tran_xx,realization%field%tran_yy,ierr)
  
  ! update mineral volume fractions
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTUpdateSolutionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  

end subroutine RTUpdateSolution

! ************************************************************************** !
!
! RTUpdateSolutionPatch: 
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RTUpdateSolutionPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Reaction_module
 
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:)  

  PetscInt :: ghosted_id, imnrl
  
  option => realization%option
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid

  aux_vars => patch%aux%RT%aux_vars


  ! update activity coefficients
  if (reaction%compute_activity) then
    do ghosted_id = 1, grid%ngmax
      call RActivity(aux_vars(ghosted_id),reaction,option)
    enddo  
  endif

  ! update mineral volume fractions
  if (reaction%nkinmnrl > 0) then
    do ghosted_id = 1, grid%ngmax
      do imnrl = 1, reaction%nkinmnrl
        aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = aux_vars(ghosted_id)%mnrl_volfrac(imnrl) + &
                                                   aux_vars(ghosted_id)%mnrl_rate(imnrl)* &
                                                   reaction%kinmnrl_molar_vol(imnrl)* &
                                                   option%dt
      enddo
    enddo
  endif
  
end subroutine RTUpdateSolutionPatch

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
  use Patch_module
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
  type(reaction_type), pointer :: reaction
  PetscReal, pointer :: xx_p(:), porosity_loc_p(:), saturation_loc_p(:), &
                        tor_loc_p(:), volume_p(:), accum_p(:), density_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  aux_vars => patch%aux%RT%aux_vars
  grid => patch%grid
  reaction => realization%reaction

  call GridVecGetArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)

  call RTUpdateAuxVarsPatch(realization,PETSC_FALSE)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*reaction%ncomp
    istart = iend-reaction%ncomp+1

    call RTAuxVarCompute(xx_p(istart:iend),aux_vars(ghosted_id),reaction, &
                         option)
    call RTAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                        volume_p(local_id), &
                        reaction,option,accum_p(istart:iend)) 
  enddo

  call GridVecRestoreArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)


end subroutine RTUpdateFixedAccumulationPatch

! ************************************************************************** !
!
! RTNumericalJacobianTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RTNumericalJacobianTest(realization)

  use Realization_module
  use Patch_module
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

  call VecDuplicate(field%tran_xx,xx_pert,ierr)
  call VecDuplicate(field%tran_xx,res,ierr)
  call VecDuplicate(field%tran_xx,res_pert,ierr)
  
  call MatCreate(option%comm,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, &
                   grid%nlmax*option%ntrandof, &
                   grid%nlmax*option%ntrandof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call RTResidual(PETSC_NULL_OBJECT,field%tran_xx,res,realization,ierr)
  call GridVecGetArrayF90(grid,res,vec2_p,ierr)
  do idof = 1,grid%nlmax*option%ntrandof
    icell = (idof-1)/option%ntrandof+1
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    endif
    call veccopy(field%tran_xx,xx_pert,ierr)
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
    call GridVecRestoreArrayF90(grid,res_pert,vec_p,ierr)
  enddo
  call GridVecRestoreArrayF90(grid,res,vec2_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerASCIIOpen(option%comm,'RTnumerical_jacobian.out',viewer,ierr)
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
subroutine RTAccumulationDerivative(aux_var,por,vol,reaction,option,J)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  PetscReal :: por, vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp, iphase
  PetscReal :: psv_t, psvd_t, v_t
  
  iphase = 1
  ! units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  !         *(kg water/m^3 water) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  if (associated(aux_var%dtotal)) then ! units of dtotal = kg water/m^3 water
    if (associated(aux_var%dtotal_sorb)) then
      v_t = vol/option%dt   
      psv_t = por*aux_var%sat(iphase)*v_t
      J = aux_var%dtotal(:,:,iphase)*psv_t + aux_var%dtotal_sorb(:,:)*v_t
    else
      psv_t = por*aux_var%sat(iphase)*vol/option%dt  
      J = aux_var%dtotal(:,:,iphase)*psv_t
    endif
  else
    J = 0.d0
    psvd_t = por*aux_var%sat(iphase)*vol*aux_var%den(iphase)/option%dt ! units of den = kg water/m^3 water
    do icomp=1,reaction%ncomp
      J(icomp,icomp) = psvd_t
    enddo
  endif

end subroutine RTAccumulationDerivative

! ************************************************************************** !
!
! RTAccumulation: Computes accumulation term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTAccumulation(aux_var,por,vol,reaction,option,Res)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  PetscReal :: por, vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  
  PetscInt :: iphase
  PetscReal :: psv_t
  PetscReal :: v_t
  
  iphase = 1
  ! units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  !         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  ! 1000.d0 converts vol from m^3 -> L
  ! all residual entries should be in mol/sec
  if (associated(aux_var%total_sorb)) then
    v_t = vol*1000.d0/option%dt
    psv_t = por*aux_var%sat(iphase)*v_t
    Res(:) = psv_t*aux_var%total(:,iphase) +  v_t*aux_var%total_sorb(:)
  else
    psv_t = por*aux_var%sat(iphase)*vol*1000.d0/option%dt  
    Res(:) = psv_t*aux_var%total(:,iphase) 
  endif
  

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
  use Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module
  use Grid_module

  implicit none
  
  interface
     subroutine samrpetscobjectstateincrease(vec)
       implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
       Vec :: vec
     end subroutine samrpetscobjectstateincrease
  end interface

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscReal, pointer :: xx_p(:), log_xx_p(:)
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscViewer :: viewer  
  
  field => realization%field
  discretization => realization%discretization
  
  ! Communication -----------------------------------------
  if (realization%reaction%use_log_formulation) then
    ! have to convert the log concentration to non-log form
    cur_level => realization%level_list%first
    do
      if (.not.associated(cur_level)) exit
      cur_patch => cur_level%patch_list%first
      do
        if (.not.associated(cur_patch)) exit
        call GridVecGetArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
        call GridVecGetArrayF90(cur_patch%grid,xx,log_xx_p,ierr)
        xx_p(:) = exp(log_xx_p(:))
        call GridVecRestoreArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
        call GridVecRestoreArrayF90(cur_patch%grid,xx,log_xx_p,ierr)  
        cur_patch => cur_patch%next
      enddo
      cur_level => cur_level%next
    enddo
    call DiscretizationGlobalToLocal(discretization,field%tran_xx,field%tran_xx_loc,NTRANDOF)
  else
    call DiscretizationGlobalToLocal(discretization,xx,field%tran_xx_loc,NTRANDOF)
  endif

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

  if(discretization%itype==AMR_GRID) then
     call samrpetscobjectstateincrease(r)
  endif

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%comm,'RTresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%comm,'RTxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
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
  use Patch_module
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
  
  PetscReal, pointer :: r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), saturation_loc_p(:), tor_loc_p(:), &
                        volume_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscInt :: i, istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  PetscReal :: Res(realization%reaction%ncomp)
  PetscViewer :: viewer
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  PetscReal :: qsrc, molality
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscTruth :: volumetric

  option => realization%option
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid
  aux_vars => patch%aux%RT%aux_vars
  aux_vars_bc => patch%aux%RT%aux_vars_bc
  
  call RTUpdateAuxVarsPatch(realization,PETSC_TRUE)
  patch%aux%RT%aux_vars_up_to_date = PETSC_FALSE 

  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)

  r_p = -accum_p
#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*reaction%ncomp
    istart = iend-reaction%ncomp+1
    call RTAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                        volume_p(local_id),reaction,option,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  iphase = 1
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%pressure)) then
      qsrc = source_sink%flow_condition%pressure%dataset%cur_value(1)
      if (source_sink%flow_condition%pressure%itype == &
          VOLUMETRIC_RATE_SS) then
        volumetric = PETSC_TRUE
      else
        volumetric = PETSC_FALSE
      endif
    endif
      
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      do istart = 1, reaction%ncomp
        molality = source_sink%tran_condition%cur_constraint_coupler% &
           aqueous_species%basis_molarity(istart)/aux_vars(ghosted_id)%den(1)*1000.d0
        select case(source_sink%tran_condition%itype)
          case(EQUILIBRIUM_SS)
            ! units should be mol/sec
            Res(istart) = -1.d-6* &
                          porosity_loc_p(ghosted_id)* &
                          saturation_loc_p(ghosted_id)* &
                          volume_p(local_id)* & ! convert m^3 water -> L water
                          (molality*aux_vars(ghosted_id)%den(1) - & 
                           aux_vars(ghosted_id)%total(istart,iphase)*1000.d0) ! convert kg water/L water -> kg water/m^3 water
          case(MASS_RATE_SS)
            Res(istart) = -molality ! actually moles/sec
          case(DIRICHLET_BC)
            if (qsrc > 0) then ! injection
              if (volumetric) then ! qsrc is volumetric; must be converted to mass
                Res(istart) = -qsrc*aux_vars(ghosted_id)%den(1) * &
                              molality
              else
                 Res(istart) = -qsrc*molality
              endif
            else ! extraction
              if (volumetric) then ! qsrc is volumetric; must be converted to mass
                Res(istart) = -qsrc*aux_vars(ghosted_id)%den(1)* &
                              aux_vars(ghosted_id)%total(istart,iphase)
              else
                Res(istart) = -qsrc* &
                              aux_vars(ghosted_id)%total(istart,iphase)/aux_vars(ghosted_id)%den(1)*1000.d0 ! convert kg water/L water -> kg water/m^3 water
              endif
            endif
          case default
        end select
      enddo
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)                                  
    enddo
    source_sink => source_sink%next
  enddo
#endif
#if 1
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
        iend = local_id_up*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:reaction%ncomp)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif
#if 1
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

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      call TBCFlux(boundary_condition%tran_condition%itype, &
                   aux_vars_bc(sum_connection), &
                   aux_vars(ghosted_id), &
                   porosity_loc_p(ghosted_id), &
                   tor_loc_p(ghosted_id), &
                   saturation_loc_p(ghosted_id), &
                   cur_connection_set%dist(0,iconn), &
                   cur_connection_set%area(iconn), &
                   option,patch%boundary_velocities(:,sum_connection),Res)
 
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:reaction%ncomp)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif  
#if 1  
! Reactions
  if (associated(reaction)) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      Res = 0.d0
      Jup = 0.d0
      call RReaction(Res,Jup,PETSC_FALSE,aux_vars(ghosted_id), &
                     volume_p(local_id),reaction,option)
      r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)                    
    enddo
  endif
#endif

  if (patch%aux%RT%inactive_cells_exist) then
    do i=1,patch%aux%RT%n_zero_rows
      r_p(patch%aux%RT%zero_rows_local(i)) = 0.d0
    enddo
  endif

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

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
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module

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

#if 0
  call RTNumericalJacobianTest(realization)
#endif

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
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      grid => cur_patch%grid
      ! need to set the current patch in the Jacobian operator
      ! so that entries will be set correctly
      if(associated(grid%structured_grid) .and. &
        (.not.(grid%structured_grid%p_samr_patch.eq.0))) then
         call SAMRSetCurrentJacobianPatch(J,grid%structured_grid%p_samr_patch)
      endif

      call RTJacobianPatch(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(realization%option%comm,'RTjacobian.out', &
                              viewer,ierr)
    call MatView(J,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  if (realization%reaction%use_log_formulation) then
    call MatDiagonalScaleLocal(J,realization%field%tran_work_loc,ierr)

    if (realization%debug%matview_Jacobian) then
      call PetscViewerASCIIOpen(realization%option%comm,'RTjacobianLog.out', &
                                viewer,ierr)
      call MatView(J,viewer,ierr)
      call PetscViewerDestroy(viewer,ierr)
    endif
    
  endif
  
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
  use Patch_module
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
  
  PetscReal, pointer :: r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), saturation_loc_p(:), tor_loc_p(:), &
                        volume_p(:), work_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
      
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Jdn(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)  
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn, rdum
  PetscReal :: qsrc
  PetscTruth :: volumetric
  PetscViewer :: viewer
  
  option => realization%option
  field => realization%field
  patch => realization%patch  
  reaction => realization%reaction
  grid => patch%grid
  aux_vars => patch%aux%RT%aux_vars
  aux_vars_bc => patch%aux%RT%aux_vars_bc


  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
    
#if 1  
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*reaction%ncomp
    istart = iend-reaction%ncomp+1
    call RTAccumulationDerivative(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                                  volume_p(local_id),reaction,option,Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)                        
  enddo
#endif
#if 1
  ! Source/Sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%pressure)) then
      qsrc = source_sink%flow_condition%pressure%dataset%cur_value(1)
      if (source_sink%flow_condition%pressure%itype == &
          VOLUMETRIC_RATE_SS) then
        volumetric = PETSC_TRUE
      else
        volumetric = PETSC_FALSE
      endif
    endif
      
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      Jup = 0.d0
      do istart = 1, reaction%ncomp
        select case(source_sink%tran_condition%itype)
          case(EQUILIBRIUM_SS)
            Jup(istart,istart) = 1.d-6* &
                                 porosity_loc_p(ghosted_id)* &
                                 saturation_loc_p(ghosted_id)* &
                                 aux_vars(ghosted_id)%den(1)* &
                                 volume_p(local_id)
          case(MASS_RATE_SS)
          case(DIRICHLET_BC)
            if (qsrc < 0) then ! extraction
              if (volumetric) then ! qsrc is volumetric; must be converted to mass
                Jup(istart,istart) = -qsrc*aux_vars(ghosted_id)%den(1)
              else
                Jup(istart,istart) = -qsrc
              endif
            endif
          case default
        end select
      enddo 
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)                        
    enddo
    source_sink => source_sink%next
  enddo
#endif
#if 1
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
                 dist_dn, &
                 cur_connection_set%area(iconn),option, &
                 patch%internal_velocities(:,iconn),Jup,Jdn)

      if (local_id_up>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)        
      endif
   
      if (local_id_dn>0) then
        Jup = -Jup
        Jdn = -Jdn
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
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      call TBCFluxDerivative(boundary_condition%tran_condition%itype, &
                   aux_vars_bc(sum_connection), &
                   aux_vars(ghosted_id), &
                   porosity_loc_p(ghosted_id), &
                   tor_loc_p(ghosted_id), &
                   saturation_loc_p(ghosted_id), &
                   cur_connection_set%dist(0,iconn), &
                   cur_connection_set%area(iconn), &
                   option,patch%boundary_velocities(:,sum_connection),Jdn)
 
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif 
#if 1  
! Reactions
  if (associated(reaction)) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      Res = 0.d0
      Jup = 0.d0
      call RReactionDerivative(Res,Jup,aux_vars(ghosted_id), &
                               volume_p(local_id),reaction,option)
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1, &
                                    Jup,ADD_VALUES,ierr)                        
    enddo
  endif
#endif
 
  if (reaction%use_log_formulation) then
    call GridVecGetArrayF90(grid,field%tran_work_loc, work_loc_p, ierr)
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      iend = ghosted_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) then
          work_loc_p(istart:iend) = 1.d0
        else
          work_loc_p(istart:iend) = aux_vars(ghosted_id)%primary_molal(:)
        endif
      else
        work_loc_p(istart:iend) = aux_vars(ghosted_id)%primary_molal(:)
      endif
    enddo
    call GridVecRestoreArrayF90(grid,field%tran_work_loc, work_loc_p, ierr)
  endif

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  
  if (patch%aux%RT%inactive_cells_exist) then
    rdum = 1.d0
    call MatZeroRowsLocal(A,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,rdum,ierr) 
  endif

end subroutine RTJacobianPatch

! ************************************************************************** !
!
! RTUpdateDenAndSat: Updates the densities and saturations in auxilliary 
!                    variables associated with reactive transport
! author: Glenn Hammond
! date: 11/03/08
!
! ************************************************************************** !
subroutine RTUpdateDenAndSat(realization,weight)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: weight
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTUpdateDenAndSatPatch(realization,weight)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine RTUpdateDenAndSat

! ************************************************************************** !
!
! RTUpdateDenAndSatPatch: Updates the densities and saturations in auxilliary 
!                         variables associated with reactive transport
! author: Glenn Hammond
! date: 11/03/08
!
! ************************************************************************** !
subroutine RTUpdateDenAndSatPatch(realization,weight)

  use Realization_module
  use Patch_module
  use Reactive_Transport_Aux_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  PetscReal :: weight
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch

  PetscInt :: ghosted_id, local_id
  PetscReal, pointer :: density_loc_p(:)
  PetscReal, pointer :: density0_loc_p(:)
  PetscReal, pointer :: saturation_loc_p(:)
  PetscReal, pointer :: saturation0_loc_p(:)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  
  call GridVecGetArrayF90(grid,field%density0_loc,density0_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%density_loc, density_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%saturation0_loc,saturation0_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%saturation_loc,saturation_loc_p,ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    ! interpolate density and saturation based on weight
    patch%aux%RT%aux_vars(ghosted_id)%den(1) = &
      (weight*density_loc_p(ghosted_id)+ &
       (1.d0-weight)*density0_loc_p(ghosted_id))
    patch%aux%RT%aux_vars(ghosted_id)%sat(1) = &
      (weight*saturation_loc_p(ghosted_id)+ &
       (1.d0-weight)*saturation0_loc_p(ghosted_id))
       
  enddo    
  
  call GridVecRestoreArrayF90(grid,field%density0_loc, density0_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%density_loc, density_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%saturation0_loc, saturation0_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)  
  
end subroutine RTUpdateDenAndSatPatch

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
      call RTUpdateAuxVarsPatch(realization,PETSC_FALSE)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine RTUpdateAuxVars

! ************************************************************************** !
!
! RTUpdateAuxVarsPatch: Updates the auxilliary variables associated with 
!                       reactive transport
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateAuxVarsPatch(realization,update_bcs)

  use Realization_module
  use Patch_module
  use Reactive_Transport_Aux_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  PetscTruth :: update_bcs
  PetscTruth :: time0
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%reaction%ncomp)
  PetscReal, pointer :: basis_molarity_p(:)
  PetscReal :: weight
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  reaction => realization%reaction
  
  call GridVecGetArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*reaction%ncomp
    istart = iend-reaction%ncomp+1
    
    call RTAuxVarCompute(xx_loc_p(istart:iend), &
                         patch%aux%RT%aux_vars(ghosted_id),reaction, &
                         option)
  enddo

  if (update_bcs) then
  
    boundary_condition => patch%boundary_conditions%first
    sum_connection = 0    
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set

      basis_molarity_p => boundary_condition%tran_condition% &
        cur_constraint_coupler%aqueous_species%basis_molarity

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
        endif

        select case(boundary_condition%tran_condition%itype)
          case(CONCENTRATION_SS,DIRICHLET_BC,NEUMANN_BC)
            ! since basis_molarity is in molarity, must convert to molality
            ! by dividing by density of water (mol/L -> mol/kg)
            xxbc(1:reaction%ncomp) = basis_molarity_p(1:reaction%ncomp) / &
              patch%aux%RT%aux_vars_bc(sum_connection)%den(1) * 1000.d0
  !          xxbc(1:reaction%ncomp)
  !            boundary_condition%tran_aux_real_var(1:reaction%ncomp,iconn)
          case(ZERO_GRADIENT_BC)
            do idof=1,reaction%ncomp
              xxbc(idof) = xx_loc_p((ghosted_id-1)*reaction%ncomp+idof)
            enddo
        end select
        ! no need to update boundary fluid density since it is already set
        call RTAuxVarCompute(xxbc, &
                             patch%aux%RT%aux_vars_bc(sum_connection), &
                             reaction,option)
      enddo
      boundary_condition => boundary_condition%next
    enddo

    patch%aux%RT%aux_vars_up_to_date = PETSC_TRUE

  endif 
  
  call GridVecRestoreArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)
  
end subroutine RTUpdateAuxVarsPatch

! ************************************************************************** !
!
! RTCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RTCreateZeroArray(patch,reaction,option)

  use Patch_module
  use Grid_module
  use Option_module
  
  implicit none

  type(patch_type) :: patch
  type(reaction_type) :: reaction
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
        n_zero_rows = n_zero_rows + reaction%ncomp
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
        do icomp = 1, reaction%ncomp
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*reaction%ncomp+icomp
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*reaction%ncomp+icomp-1
        enddo
      else
      endif
    enddo
  endif

  patch%aux%RT%zero_rows_local => zero_rows_local
  patch%aux%RT%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%RT%n_zero_rows = n_zero_rows  

  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                     option%comm,ierr)

  if (flag > 0) patch%aux%RT%inactive_cells_exist = PETSC_TRUE

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif
  
end subroutine RTCreateZeroArray

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
  use Patch_module
  use Grid_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(grid_type), pointer :: grid 
  PetscReal, pointer :: dxx_ptr(:), xx_ptr(:), yy_ptr(:)
  
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  option%dcmax=0.D0
  
  call VecWAXPY(field%tran_dxx,-1.d0,field%tran_xx,field%tran_yy,ierr)
  
  call VecStrideNorm(field%tran_dxx,ZERO_INTEGER,NORM_INFINITY,option%dcmax,ierr)
      
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
  type(reaction_type), pointer :: reaction
  PetscInt :: i
  
  option => realization%option
  reaction => realization%reaction
  
  string = ''
  write(string2,'('',"'',a,''"'')') trim('pH')
  string = trim(string) // trim(string2)
  
  do i=1,option%ntrandof
    write(string2,'('',"'',a,''"'')') trim(reaction%primary_species_names(i))
    string = trim(string) // trim(string2)
  enddo
  
  do i=1,realization%reaction%nkinmnrl
    write(string2,'('',"'',a,''"'')') trim(reaction%kinmnrl_names(i))
    string = trim(string) // trim(string2)
  enddo
  
  RTGetTecplotHeader = string

end function RTGetTecplotHeader


! ************************************************************************** !
!
! RTAuxVarCompute: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine RTAuxVarCompute(x,aux_var,reaction,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: x(reaction%ncomp)
  type(reactive_transport_auxvar_type) :: aux_var
  
  PetscReal :: den ! kg water/L water

#if 0  
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp
  PetscReal :: dtotal(reaction%ncomp,reaction%ncomp)
  PetscReal :: dtotalsorb(reaction%ncomp,reaction%ncomp)
  PetscReal :: pert
  type(reactive_transport_auxvar_type) :: auxvar_pert
#endif

  ! any changes to the below must also be updated in 
  ! Reaction.F90:RReactionDerivative()
  
  aux_var%primary_molal = x
  
  den = aux_var%den(1)*1.d-3 ! convert kg water/m^3 water -> kg water/L water
  aux_var%primary_spec = aux_var%primary_molal*den
  call RTotal(aux_var,reaction,option)
  if (reaction%nsorb > 0) then
    call RTotalSorb(aux_var,reaction,option)
  endif

#if 0  
! numerical check
  Res_orig = 0.d0
  dtotal = 0.d0
  dtotalsorb = 0.d0
  call RTAuxVarInit(auxvar_pert,reaction,option)
  call RTAuxVarCopy(auxvar_pert,aux_var,option)
  do jcomp = 1, reaction%ncomp
    Res_pert = 0.d0
    call RTAuxVarCopy(auxvar_pert,aux_var,option)
    if (reaction%neqcmplx > 0) then
      aux_var%secondary_spec = 0.d0
    endif
    if (reaction%neqsurfcmplxrxn > 0) then
      auxvar_pert%eqsurfcmplx_freesite_conc = 1.d-9
      auxvar_pert%eqsurfcmplx_spec = 0.d0
    endif
    if (reaction%neqionxrxn > 0) then
      aux_var%eqionx_ref_cation_sorbed_conc = 1.d-9
    endif
    pert = auxvar_pert%primary_molal(jcomp)*perturbation_tolerance
    auxvar_pert%primary_molal(jcomp) = auxvar_pert%primary_molal(jcomp) + pert
    
    ! this is essentially what RTAuxVarCompute() performs
!      call RTAuxVarCompute(auxvar_pert%primary_molal,auxvar_pert,option)      
    auxvar_pert%primary_spec = auxvar_pert%primary_molal*den
    call RTotal(auxvar_pert,reaction,option)
    if (reaction%nsorb > 0) call RTotalSorb(auxvar_pert,reaction,option)
    dtotal(:,jcomp) = (auxvar_pert%total(:,1) - aux_var%total(:,1))/pert
    if (reaction%nsorb > 0) dtotalsorb(:,jcomp) = (auxvar_pert%total_sorb(:) - aux_var%total_sorb(:))/pert
  enddo
  do icomp = 1, reaction%ncomp
    do jcomp = 1, reaction%ncomp
      if (dabs(dtotal(icomp,jcomp)) < 1.d-16) dtotal(icomp,jcomp) = 0.d0
      if (reaction%nsorb > 0) then
        if (dabs(dtotalsorb(icomp,jcomp)) < 1.d-16) dtotalsorb(icomp,jcomp) = 0.d0
      endif
    enddo
  enddo
  dtotal = dtotal * 1000.d0
  if (reaction%nsorb > 0) dtotalsorb = dtotalsorb * 1000.d0
  aux_var%dtotal(:,:,1) = dtotal
  if (reaction%nsorb > 0) aux_var%dtotal_sorb = dtotalsorb
  call RTAuxVarDestroy(auxvar_pert)
#endif
  
end subroutine RTAuxVarCompute

end module Reactive_Transport_module
