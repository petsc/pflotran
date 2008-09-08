module Reactive_Transport_module

  use Transport_module
  use Chemistry_module
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
  call RTInitializeTimestep(realization)  
 
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
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%RT => RTAuxCreate()
    
  ! allocate aux_var data structures for all grid cells
  allocate(patch%aux%RT%aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RTAuxVarInit(patch%aux%RT%aux_vars(ghosted_id),option)
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
    call RTAuxVarInit(patch%aux%RT%aux_vars_bc(iconn),option)
  enddo
  patch%aux%RT%num_aux_bc = sum_connection

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call RTCreateZeroArray(patch,option)

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
      call RTSetupPatch(realization)
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
  use Chemistry_module
 
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(reaction_type), pointer :: chemistry
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:)  

  PetscInt :: ghosted_id, imnrl
  
  option => realization%option
  patch => realization%patch
  chemistry => realization%chemistry
  grid => patch%grid

  aux_vars => patch%aux%RT%aux_vars

  do ghosted_id = 1, grid%ngmax
    do imnrl = 1, option%nmnrl
      aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = aux_vars(ghosted_id)%mnrl_volfrac(imnrl) + &
                                                 aux_vars(ghosted_id)%mnrl_rate(imnrl)* &
                                                 chemistry%mnrl_molar_vol(imnrl)* &
                                                 option%dt
    enddo
  enddo
  
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
  type(reaction_type), pointer :: chemistry
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
  chemistry => realization%chemistry

  call GridVecGetArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%saturation0_loc,saturation_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%density0_loc,density_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%ncomp
    istart = iend-option%ncomp+1
    call RTAuxVarCompute(xx_p(istart:iend),aux_vars(ghosted_id),chemistry,option)
    call RTAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                        saturation_loc_p(ghosted_id), &
                        density_loc_p(ghosted_id), &
                        volume_p(local_id), &
                        option,accum_p(istart:iend)) 
  enddo

  call GridVecRestoreArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%saturation0_loc,saturation_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%density0_loc,density_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)

#if 0
  call RTNumericalJacobianTest(field%tran_xx,realization)
#endif

end subroutine RTUpdateFixedAccumulationPatch

! ************************************************************************** !
!
! RTAuxVarCompute: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine RTAuxVarCompute(x,aux_var,chemistry,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  PetscReal :: x(option%ncomp)  
  type(reaction_type) :: chemistry  
  type(reactive_transport_auxvar_type) :: aux_var

  aux_var%primary_spec = x
  call RTotal(aux_var,chemistry,option)
  
end subroutine RTAuxVarCompute

! ************************************************************************** !
!
! RTNumericalJacobianTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RTNumericalJacobianTest(xx,realization)

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
  call GridVecGetArrayF90(grid,res,vec2_p,ierr)
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
      call GridVecRestoreArrayF90(grid,res_pert,vec_p,ierr)
    enddo
  enddo
  call GridVecRestoreArrayF90(grid,res,vec2_p,ierr)

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
subroutine RTAccumulationDerivative(aux_var,por,sat,den,vol,option,J)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  PetscReal :: por, sat, vol, den
  type(option_type) :: option
  PetscReal :: J(option%ncomp,option%ncomp)
  
  PetscInt :: icomp, iphase
  PetscReal :: psdv_t
  
  iphase = 1
  psdv_t = por*sat*den*vol/option%dt
  if (associated(aux_var%dtotal)) then
    J = aux_var%dtotal(:,:,iphase)*psdv_t
  else
    J = 0.d0
    do icomp=1,option%ncomp
      J(icomp,icomp) = psdv_t
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
subroutine RTAccumulation(aux_var,por,sat,den,vol,option,Res)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  PetscReal :: por, sat, vol, den
  type(option_type) :: option
  PetscReal :: Res(option%ncomp)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: psdv_t
  
  iphase = 1
  psdv_t = por*sat*den*vol/option%dt
  do icomp=1,option%ncomp
    Res(icomp) = psdv_t*aux_var%total(icomp,iphase) 
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
  use Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module

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
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  field => realization%field
  discretization => realization%discretization
  
  ! Communication -----------------------------------------
  call DiscretizationGlobalToLocal(discretization,xx,field%tran_xx_loc,NTRANDOF)

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
  
  PetscReal, pointer :: xx_loc_p(:), r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), saturation_loc_p(:), tor_loc_p(:), &
                        volume_p(:), density_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscInt :: i, istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: chemistry
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  PetscReal :: Res(realization%option%ncomp)
  PetscViewer :: viewer
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  PetscReal :: qsrc
  PetscReal :: Jup(realization%option%ncomp,realization%option%ncomp)
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  chemistry => realization%chemistry
  grid => patch%grid
  aux_vars => patch%aux%RT%aux_vars
  aux_vars_bc => patch%aux%RT%aux_vars_bc
  
  call RTUpdateAuxVars(realization)
  aux_vars_up_to_date = .false. ! override flags since they will soon be out of date  

  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%tran_xx_loc, xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%density_loc, density_loc_p, ierr)
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
    iend = local_id*option%ncomp
    istart = iend-option%ncomp+1
    call RTAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                        saturation_loc_p(ghosted_id), &
                        density_loc_p(ghosted_id), &
                        volume_p(local_id),option,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%ncomp)
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  iphase = 1
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      if (associated(source_sink%flow_condition) .and. &
          associated(source_sink%flow_condition%pressure)) then
        qsrc = source_sink%flow_condition%pressure%dataset%cur_value(1)
      endif
      
      do istart = 1, option%ncomp
        select case(source_sink%tran_condition%sub_condition_ptr(istart)%ptr%itype)
          case(EQUILIBRIUM_SS)
            Res(istart) = -1.d-6* &
                          porosity_loc_p(ghosted_id)* &
                          saturation_loc_p(ghosted_id)* &
                          density_loc_p(ghosted_id)* &
                          volume_p(local_id)* &
                          (source_sink%tran_condition%sub_condition_ptr(istart)%ptr%dataset%cur_value(1)- &
                           aux_vars(ghosted_id)%total(istart,iphase))
          case(MASS_RATE_SS)
            Res(istart) = -source_sink%tran_condition%sub_condition_ptr(istart)%ptr%dataset%cur_value(1)
          case(CONCENTRATION_SS)
            if (qsrc > 0) then ! injection
              Res(istart) = -qsrc*source_sink%tran_condition%sub_condition_ptr(istart)%ptr%dataset%cur_value(1)
            else ! extraction
              Res(istart) = -qsrc*aux_vars(ghosted_id)%total(istart,iphase)
            endif
          case default
        end select
      enddo
      iend = local_id*option%ncomp
      istart = iend-option%ncomp+1
      r_p(istart:iend) = r_p(istart:iend) + Res(1:option%ncomp)                                  
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
                 density_loc_p(ghosted_id_up), &
                 dist_up, &
                 aux_vars(ghosted_id_dn),porosity_loc_p(ghosted_id_dn), &
                 tor_loc_p(ghosted_id_dn),saturation_loc_p(ghosted_id_dn), &
                 density_loc_p(ghosted_id_dn), &
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

      call TBCFlux(boundary_condition%tran_condition%itype(1), &
                   aux_vars_bc(sum_connection), &
                   aux_vars(ghosted_id), &
                   porosity_loc_p(ghosted_id), &
                   tor_loc_p(ghosted_id), &
                   saturation_loc_p(ghosted_id), &
                   density_loc_p(ghosted_id), &
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
#if 1  
! Reactions
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%ncomp
    istart = iend-option%ncomp+1
    call RKineticMineral(Res,Jup,PETSC_TRUE,aux_vars(ghosted_id),chemistry,option)
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%ncomp)                    
  enddo
#endif

  if (inactive_cells_exist) then
    do i=1,patch%aux%RT%n_zero_rows
      r_p(patch%aux%RT%zero_rows_local(i)) = 0.d0
    enddo
  endif
  
  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_xx_loc, xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%density_loc, density_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

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
  use Level_module
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
  type(grid_type),  pointer :: grid
  
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
         call SAMRSetCurrentJacobianPatch(A, grid%structured_grid%p_samr_patch)
      endif

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
  
  PetscReal, pointer :: xx_loc_p(:), r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), saturation_loc_p(:), tor_loc_p(:), &
                        volume_p(:), density_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: chemistry
      
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  PetscReal :: Jup(realization%option%ncomp,realization%option%ncomp)
  PetscReal :: Jdn(realization%option%ncomp,realization%option%ncomp)
  PetscReal :: Res(realization%option%ncomp)  
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn, rdum
  PetscReal :: qsrc
  PetscViewer :: viewer
  
  option => realization%option
  field => realization%field
  patch => realization%patch  
  chemistry => realization%chemistry
  grid => patch%grid
  aux_vars => patch%aux%RT%aux_vars
  aux_vars_bc => patch%aux%RT%aux_vars_bc

  flag = SAME_NONZERO_PATTERN  
  call MatZeroEntries(A,ierr)

  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%tran_xx_loc, xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%density_loc, density_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
    
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
                                  saturation_loc_p(ghosted_id), &
                                  density_loc_p(ghosted_id), &
                                  volume_p(local_id),option,Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)                        
  enddo
#endif
#if 1
  ! Source/Sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      if (associated(source_sink%flow_condition) .and. &
          associated(source_sink%flow_condition%pressure)) then
        qsrc = source_sink%flow_condition%pressure%dataset%cur_value(1)
      endif
      
      Jup = 0.d0
      do istart = 1, option%ncomp
        select case(source_sink%tran_condition%sub_condition_ptr(istart)%ptr%itype)
          case(EQUILIBRIUM_SS)
            Jup(istart,istart) = 1.d-6* &
                                 porosity_loc_p(ghosted_id)* &
                                 saturation_loc_p(ghosted_id)* &
                                 density_loc_p(ghosted_id)* &
                                 volume_p(local_id)
          case(MASS_RATE_SS)
          case(CONCENTRATION_SS)
            if (qsrc < 0) then ! extraction
              Jup(istart,istart) = -qsrc
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
                 density_loc_p(ghosted_id_up),dist_up, &
                 aux_vars(ghosted_id_dn),porosity_loc_p(ghosted_id_dn), &
                 tor_loc_p(ghosted_id_dn),saturation_loc_p(ghosted_id_dn), &
                 density_loc_p(ghosted_id_dn),dist_up, &
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

      call TBCFluxDerivative(boundary_condition%tran_condition%itype(1), &
                   aux_vars_bc(sum_connection), &
                   aux_vars(ghosted_id), &
                   porosity_loc_p(ghosted_id), &
                   tor_loc_p(ghosted_id), &
                   saturation_loc_p(ghosted_id), &
                   density_loc_p(ghosted_id), &
                   cur_connection_set%dist(0,iconn), &
                   cur_connection_set%area(iconn), &
                   option,patch%boundary_velocities(:,sum_connection),Jdn)
 
      Jdn = -1.d0*Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif 
#if 1  
! Reactions
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%ncomp
    istart = iend-option%ncomp+1
!    call RKineticMineral(Res,Jup,aux_vars(ghosted_id),chemistry,option)
    call RKineticMineralDerivative(Res,Jup,aux_vars(ghosted_id),chemistry,option)
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)                        
  enddo
#endif
 

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_xx_loc, xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%saturation_loc, saturation_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%density_loc, density_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  
  if (inactive_cells_exist) then
    rdum = 1.d0
    call MatZeroRowsLocal(A,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,rdum,ierr) 
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
  use Patch_module
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
  type(reaction_type), pointer :: chemistry
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%ncomp)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  chemistry => realization%chemistry
  
  call GridVecGetArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%ncomp
    istart = iend-option%ncomp+1
   
    call RTAuxVarCompute(xx_loc_p(istart:iend), &
                         patch%aux%RT%aux_vars(ghosted_id), &
                         chemistry,option)
  enddo

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

      do idof=1,option%ncomp
        select case(boundary_condition%tran_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,NEUMANN_BC)
            xxbc(idof) = boundary_condition%tran_aux_real_var(idof,iconn)
          case(ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%ncomp+idof)
        end select
      enddo
      
      call RTAuxVarCompute(xxbc, &
                           patch%aux%RT%aux_vars_bc(sum_connection), &
                           chemistry,option)
    enddo
    boundary_condition => boundary_condition%next
  enddo


  call GridVecRestoreArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)
  
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

  patch%aux%RT%zero_rows_local => zero_rows_local
  patch%aux%RT%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%RT%n_zero_rows = n_zero_rows  

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
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none

  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase

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

  iphase = 1
  select case(ivar)
    case(PRIMARY_SPEC_CONCENTRATION)
      call VecStrideGather(field%tran_xx,isubvar,vec,INSERT_VALUES,ierr)
    case(MATERIAL_ID,TOTAL_CONCENTRATION)
      call GridVecGetArrayF90(grid,vec,vec_ptr,ierr)
      select case(ivar)
        case(TOTAL_CONCENTRATION)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)    
            vec_ptr(local_id) = patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase)
          enddo
        case(MATERIAL_ID)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
          enddo
      end select
      call GridVecRestoreArrayF90(grid,vec,vec_ptr,ierr)
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
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none

  PetscReal :: RTGetVarFromArrayAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase
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

  iphase = 1
  select case(ivar)
    case(PRIMARY_SPEC_CONCENTRATION)
      call GridVecGetArrayF90(grid,field%tran_xx,vec_ptr,ierr)
      value = vec_ptr((local_id-1)*option%ntrandof+isubvar)
      call GridVecRestoreArrayF90(grid,field%tran_xx,vec_ptr,ierr)
    case(TOTAL_CONCENTRATION)
      ghosted_id = grid%nL2G(local_id)    
      value = patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase)
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
  PetscInt :: i
  
  option => realization%option
  
  string = '' 
  do i=1,option%ntrandof
    write(string2,'('',"'',a,''"'')') trim(option%comp_names(i))
    string = trim(string) // trim(string2)
  enddo
  
  RTGetTecplotHeader = string

end function RTGetTecplotHeader

! ************************************************************************** !
!
! RTotal: Computes the total component concentrations and derivative with
!         respect to free-ion
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine RTotal(auxvar,chemistry,option)

  use Option_module
  
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type) :: chemistry
  type(option_type) :: option
  
  PetscInt :: i, j, icplx, icomp, jcomp, iphase, ncomp
  PetscReal :: ln_conc(option%ncomp)
  PetscReal :: lnQK, tempreal
  PetscReal, parameter :: log_to_ln = 2.30258509299d0

  iphase = 1                         

  ln_conc = log(auxvar%primary_spec)
  auxvar%total(:,iphase) = auxvar%primary_spec
  ! initialize derivatives
  auxvar%dtotal = 0.d0
  do icomp = 1, option%ncomp
    auxvar%dtotal(icomp,icomp,iphase) = 1.d0
  enddo
  
  do icplx = 1, chemistry%neqcmplx ! for each secondary species
    ! compute secondary species concentration
    lnQK = -1.d0*chemistry%eqcmplx_K(icplx)*log_to_ln
    ncomp = chemistry%eqcmplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = chemistry%eqcmplxspecid(i,icplx)
      lnQK = lnQK + chemistry%eqcmplxstoich(i,icplx)*ln_conc(icomp)
    enddo
    auxvar%secondary_spec(icplx) = exp(lnQK)
  
    ! add contribution to primary totals
    do i = 1, ncomp
      icomp = chemistry%eqcmplxspecid(i,icplx)
      auxvar%total(icomp,iphase) = auxvar%total(icomp,iphase) + &
                                   chemistry%eqcmplxstoich(i,icplx)* &
                                   auxvar%secondary_spec(icplx)
    enddo
    
    ! add contribution to derivatives of total with respect to free
    do j = 1, ncomp
      jcomp = chemistry%eqcmplxspecid(j,icplx)
      tempreal = chemistry%eqcmplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))
      do i = 1, ncomp
        icomp = chemistry%eqcmplxspecid(i,icplx)
        auxvar%dtotal(icomp,jcomp,iphase) = auxvar%dtotal(icomp,jcomp,iphase) + &
                                      chemistry%eqcmplxstoich(i,icplx)*tempreal
      enddo
    enddo
  enddo
  
end subroutine RTotal

! ************************************************************************** !
!
! RKineticMineral: Computes the kinetic mineral precipitation/dissolution
!                  rates
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RKineticMineral(Res,Jac,derivative,auxvar,chemistry,option)

  use Option_module
  
  type(option_type) :: option
  PetscTruth :: derivative
  PetscReal :: Res(option%ncomp)
  PetscReal :: Jac(option%ncomp,option%ncomp)
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type) :: chemistry
  
  PetscInt :: i, j, k, imnrl, icomp, jcomp, kcplx, iphase, ncomp, ipref
  PetscReal :: prefactor(10), sum_prefactor_rate
  PetscReal :: dIm_dsum_prefactor_rate, dIm_dprefactor_rate
  PetscReal :: dprefactor_dcomp_numerator, dprefactor_dcomp_denominator
  PetscReal :: tempreal, tempreal2
  PetscReal :: affinity_factor, sign_
  PetscReal :: Im, Im_const, dIm_dQK
  PetscReal :: ln_conc(option%ncomp)
  PetscReal :: ln_sec(option%ncmplx)
  PetscReal :: QK, lnQK, dQK_dCj
  PetscReal, parameter :: log_to_ln = 2.30258509299d0
  PetscTruth :: prefactor_exists

  iphase = 1                         

  ln_conc = log(auxvar%primary_spec)
  ln_sec = log(auxvar%secondary_spec)

  Res = 0.d0
  Jac = 0.d0
  
  do imnrl = 1, chemistry%nkinmnrl ! for each mineral
    ! compute secondary species concentration
    lnQK = -1.d0*chemistry%kinmnrl_K(imnrl)*log_to_ln
    ncomp = chemistry%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = chemistry%kinmnrlspecid(i,imnrl)
      lnQK = lnQK + chemistry%kinmnrlstoich(i,imnrl)*ln_conc(icomp)
    enddo
    QK = exp(lnQK)
    
    if (associated(chemistry%kinmnrl_Tempkin_const)) then
      affinity_factor = 1.d0-QK**(1/chemistry%kinmnrl_Tempkin_const(imnrl))
    else
      affinity_factor = 1.d0-QK
    endif
    
    sign_ = sign(1.d0,affinity_factor)

    if (auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then
      ! compute prefactor
      if (chemistry%kinmnrl_num_prefactors(imnrl) > 0) then
        print *, 'Kinetic mineral reaction prefactor calculations have not been verified.  Ask Glenn.'
        stop
        sum_prefactor_rate = 0
        do ipref = 1, chemistry%kinmnrl_num_prefactors(imnrl)
          prefactor(ipref) = 1.d0
          do i = 1, chemistry%kinmnrl_pri_prefactor_id(0,ipref,imnrl) ! primary contribution
            icomp = chemistry%kinmnrl_pri_prefactor_id(i,ipref,imnrl)
            prefactor(ipref) = prefactor(ipref) * &
                               exp(chemistry%kinmnrl_pri_pref_alpha_stoich(i,ipref,imnrl)* &
                                   ln_conc(icomp))/ &
                               ((1.d0+chemistry%kinmnrl_pri_pref_atten_coef(i,ipref,imnrl))* &
                                 exp(chemistry%kinmnrl_pri_pref_beta_stoich(i,ipref,imnrl)* &
                                     ln_conc(icomp)))
          enddo
          do k = 1, chemistry%kinmnrl_sec_prefactor_id(0,ipref,imnrl) ! secondary contribution
            kcplx = chemistry%kinmnrl_sec_prefactor_id(k,ipref,imnrl)
            prefactor(ipref) = prefactor(ipref) * &
                               exp(chemistry%kinmnrl_sec_pref_alpha_stoich(k,ipref,imnrl)* &
                                   ln_sec(kcplx))/ &
                               ((1.d0+chemistry%kinmnrl_sec_pref_atten_coef(i,ipref,imnrl))* &
                                 exp(chemistry%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                     ln_sec(kcplx)))
          enddo
          sum_prefactor_rate = sum_prefactor_rate + prefactor(ipref)*chemistry%kinmnrl_rate(ipref,imnrl)
        enddo
      else
        sum_prefactor_rate = chemistry%kinmnrl_rate(1,imnrl)
      endif

      ! compute rate
      Im_const = -1.d0*sign_*auxvar%mnrl_area0(imnrl)
      if (associated(chemistry%kinmnrl_affinity_power)) then
        Im = Im_const*abs(affinity_factor)**chemistry%kinmnrl_affinity_power(imnrl)
      else
        Im = Im_const*abs(affinity_factor)
      endif
      if (chemistry%kinmnrl_num_prefactors(imnrl) > 0) then
        Im = Im*sum_prefactor_rate
      endif
      auxvar%mnrl_rate(imnrl) = Im
    else
      auxvar%mnrl_rate(imnrl) = 0.d0
      cycle
    endif

    ncomp = chemistry%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = chemistry%kinmnrlspecid(i,imnrl)
      Res(icomp) = Res(icomp) + chemistry%kinmnrlstoich(i,imnrl)*Im
    enddo 
    
    if (.not. derivative) cycle   

    ! calculate derivatives of rate with respect to free
    dIm_dQK = -1.d0*Im
    if (associated(chemistry%kinmnrl_affinity_power)) then
      dIm_dQK = dIm_dQK*chemistry%kinmnrl_affinity_power(imnrl)/affinity_factor
    endif
    if (associated(chemistry%kinmnrl_Tempkin_const)) then
      dIm_dQK = dIm_dQK*(1/chemistry%kinmnrl_Tempkin_const(imnrl))/QK
    endif
    
    ! derivatives with respect to primary species in reaction quotient
    do j = 1, ncomp
      jcomp = chemistry%kinmnrlspecid(j,imnrl)
      dQK_dCj = chemistry%kinmnrlstoich(j,imnrl)*exp(lnQK-ln_conc(jcomp))
      do i = 1, ncomp
        icomp = chemistry%kinmnrlspecid(i,imnrl)
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                           chemistry%kinmnrlstoich(i,imnrl)*dIm_dQK*dQK_dCj
      enddo
    enddo

    if (chemistry%kinmnrl_num_prefactors(imnrl) > 0) then ! add contribution of derivative in prefactor - messy
      print *, 'Kinetic mineral reaction prefactor calculations have not been verified.  Ask Glenn.'
      stop
      
      dIm_dsum_prefactor_rate = Im/sum_prefactor_rate
      do ipref = 1, chemistry%kinmnrl_num_prefactors(imnrl)
        dIm_dprefactor_rate = dIm_dsum_prefactor_rate*chemistry%kinmnrl_rate(ipref,imnrl)
        do j = 1, chemistry%kinmnrl_pri_prefactor_id(0,ipref,imnrl) ! primary contribution
          jcomp = chemistry%kinmnrl_pri_prefactor_id(j,ipref,imnrl)
          ! numerator
          dprefactor_dcomp_numerator = chemistry%kinmnrl_pri_pref_alpha_stoich(j,ipref,imnrl)* &
                                       prefactor(ipref)/auxvar%primary_spec(jcomp)
          ! denominator
          dprefactor_dcomp_denominator = -1.d0*prefactor(ipref)/ &
                                         chemistry%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)* &
                                         chemistry%kinmnrl_pri_pref_atten_coef(j,ipref,imnrl)* &
                                         exp((chemistry%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)-1.d0)* &
                                             ln_conc(jcomp))* &
                                         ((1.d0+chemistry%kinmnrl_pri_pref_atten_coef(j,ipref,imnrl))* &
                                           exp(chemistry%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)* &
                                               ln_conc(jcomp)))
          tempreal = dIm_dprefactor_rate*(dprefactor_dcomp_numerator+dprefactor_dcomp_denominator)
          do i = 1, ncomp
            icomp = chemistry%kinmnrlspecid(i,imnrl)
            Jac(icomp,jcomp) = Jac(icomp,jcomp) + chemistry%kinmnrlstoich(i,imnrl)*tempreal
          enddo  ! loop over col
        enddo !loop over row
        do k = 1, chemistry%kinmnrl_sec_prefactor_id(0,ipref,imnrl) ! secondary contribution
          kcplx = chemistry%kinmnrl_sec_prefactor_id(k,ipref,imnrl)
          ! numerator
          dprefactor_dcomp_numerator = chemistry%kinmnrl_sec_pref_alpha_stoich(k,ipref,imnrl)* &
                                       prefactor(ipref)/auxvar%secondary_spec(jcomp)
          ! denominator
          dprefactor_dcomp_denominator = -1.d0*prefactor(ipref)/ &
                                         (1.d0+chemistry%kinmnrl_sec_pref_atten_coef(k,ipref,imnrl)* &
                                          exp(chemistry%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                              ln_sec(kcplx)))* &
                                         chemistry%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                         chemistry%kinmnrl_sec_pref_atten_coef(k,ipref,imnrl)* &
                                         exp((chemistry%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)-1.d0)* &
                                              ln_sec(kcplx))
          tempreal = dIm_dprefactor_rate*(dprefactor_dcomp_numerator+dprefactor_dcomp_denominator)
          do j = 1, chemistry%eqcmplxstoich(0,kcplx)
            jcomp = chemistry%eqcmplxstoich(j,kcplx)
            tempreal2 = chemistry%eqcmplxstoich(j,kcplx)*exp(ln_sec(kcplx)-ln_conc(jcomp))
            do i = 1, ncomp
              icomp = chemistry%kinmnrlspecid(i,imnrl)
              Jac(icomp,jcomp) = Jac(icomp,jcomp) + chemistry%kinmnrlstoich(i,imnrl)*tempreal* &
                                                    tempreal2
            enddo  ! loop over col
          enddo  ! loop over row
        enddo  ! loop over complexes
      enddo  ! loop over prefactors
    endif
  enddo  ! loop over minerals
    
end subroutine RKineticMineral

! ************************************************************************** !
!
! RKineticMineralDerivative: Computes the kinetic mineral precipitation/dissolution
!                  rates
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RKineticMineralDerivative(Res,Jac,auxvar,chemistry,option)

  use Option_module
  
  type(option_type) :: option
  PetscReal :: Res(option%ncomp)
  PetscReal :: Jac(option%ncomp,option%ncomp)
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type) :: chemistry
  
  type(reactive_transport_auxvar_type) :: auxvar_pert
  PetscReal :: Res_orig(option%ncomp)
  PetscReal :: Res_pert(option%ncomp)
  PetscReal :: Jac_dummy(option%ncomp,option%ncomp)
  PetscReal :: pert

  PetscInt :: icomp, jcomp

  Res = 0.d0
  Jac = 0.d0

  if (option%numerical_derivatives) then
    call RTAuxVarInit(auxvar_pert,option)
    call RTAuxVarCopy(auxvar,auxvar_pert,option)
    call RKineticMineral(Res_orig,Jac_dummy,PETSC_FALSE,auxvar,chemistry,option)
    do jcomp = 1, option%ncomp
      call RTAuxVarCopy(auxvar,auxvar_pert,option)
      pert = auxvar_pert%primary_spec(jcomp)*perturbation_tolerance
      auxvar_pert%primary_spec(jcomp) = auxvar_pert%primary_spec(jcomp) + pert
      call RTotal(auxvar,chemistry,option)
      call RKineticMineral(Res_pert,Jac_dummy,PETSC_FALSE,auxvar,chemistry,option)
      do icomp = 1, option%ncomp
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + (Res_pert(icomp)-Res_orig(icomp))/pert
      enddo
    enddo
    call RTAuxVarDestroy(auxvar_pert)
  else
    call RKineticMineral(Res_orig,Jac,PETSC_TRUE,auxvar,chemistry,option)
  endif
    
end subroutine RKineticMineralDerivative

end module Reactive_Transport_module
