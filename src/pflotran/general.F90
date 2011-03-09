module General_module

  use General_Aux_module
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


! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-5

  public GeneralResidual, GeneralJacobian, &
         GeneralUpdateFixedAccum, GeneralTimeCut,&
         GeneralSetup, GeneralNumericalJacTest, &
         GeneralInitializeTimestep, GeneralUpdateAuxVars, &
         GeneralMaxChange, GeneralUpdateSolution, &
         GeneralGetTecplotHeader, GeneralComputeMassBalance, &
         GeneralDestroy

contains

! ************************************************************************** !
!
! GeneralTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralTimeCut(realization)
 
  use Realization_module
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr
  PetscInt :: local_id

  option => realization%option
  field => realization%field

  call VecCopy(field%flow_yy,field%flow_xx,ierr)
  call GeneralInitializeTimestep(realization)  
 
end subroutine GeneralTimeCut

! ************************************************************************** !
!
! GeneralSetup: 
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralSetup(realization)

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
      call GeneralSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeneralSetup

! ************************************************************************** !
!
! GeneralSetupPatch: Creates arrays for auxilliary variables
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralSetupPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Material_module
  use Material_Aux_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, idof
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:), gen_aux_vars_bc(:)
  type(material_type), pointer :: material  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%General => GeneralAuxCreate()
  patch%aux%Material => MaterialAuxCreate()

  allocate(patch%aux%Material%material_parameter% &
    sir(size(realization%saturation_function_array),option%nphase))
  patch%aux%Material%material_parameter%sir = -999.d0
  do i = 1, size(realization%saturation_function_array)
    if (associated(realization%saturation_function_array(i)%ptr)) then
      patch%aux%Material%material_parameter% &
        sir(realization%saturation_function_array(i)%ptr%id,:) = &
        realization%saturation_function_array(i)%ptr%Sr(:)
    endif
  enddo

  allocate(patch%aux%Material%material_parameter% &
    dencpr(size(realization%material_property_array)))
  patch%aux%Material%material_parameter%dencpr = -999.d0
  do i = 1, size(realization%material_property_array)
    if (associated(realization%material_property_array(i)%ptr)) then
      patch%aux%Material%material_parameter% &
        dencpr(realization%saturation_function_array(i)%ptr%id) = &
        realization%material_property_array(i)%ptr%rock_density * &
        realization%material_property_array(i)%ptr%specific_heat
    endif
  enddo

  ! allocate aux_var data structures for all grid cells  
  allocate(gen_aux_vars(option%nflowdof+1,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 1, option%nflowdof+1
      call GeneralAuxVarInit(gen_aux_vars(idof,ghosted_id),option)
    enddo
  enddo
  patch%aux%General%aux_vars => gen_aux_vars
  patch%aux%General%num_aux = grid%ngmax

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
  allocate(gen_aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call GeneralAuxVarInit(gen_aux_vars_bc(iconn),option)
  enddo
  patch%aux%General%aux_vars_bc => gen_aux_vars_bc
  patch%aux%General%num_aux_bc = sum_connection
  
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call GeneralCreateZeroArray(patch,option)

end subroutine GeneralSetupPatch

! ************************************************************************** !
!
! GeneralSetup: 
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralComputeMassBalance(realization,mass_balance)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  mass_balance = 0.d0
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GeneralComputeMassBalancePatch(realization,mass_balance)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeneralComputeMassBalance

! ************************************************************************** !
!
! GeneralComputeMassBalancePatch: Initializes mass balance
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralComputeMassBalancePatch(realization,mass_balance)
 
  use Realization_module
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  PetscReal, pointer :: volume_p(:), porosity_loc_p(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_aux_vars => patch%aux%Global%aux_vars

  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! mass = volume*saturation*density
    mass_balance = mass_balance + &
      global_aux_vars(ghosted_id)%den_kg* &
      global_aux_vars(ghosted_id)%sat* &
      porosity_loc_p(ghosted_id)*volume_p(local_id)
  enddo

  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  
end subroutine GeneralComputeMassBalancePatch

! ************************************************************************** !
!
! GeneralZeroMassBalDeltaPatch: Zeros mass balance delta array
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralZeroMassBalDeltaPatch(realization)
 
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  do iconn = 1, patch%aux%General%num_aux_bc
    global_aux_vars_bc(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine GeneralZeroMassBalDeltaPatch

! ************************************************************************** !
!
! GeneralUpdateMassBalancePatch: Updates mass balance
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralUpdateMassBalancePatch(realization)
 
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  do iconn = 1, patch%aux%General%num_aux_bc
    global_aux_vars_bc(iconn)%mass_balance = &
      global_aux_vars_bc(iconn)%mass_balance + &
      global_aux_vars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
  enddo

end subroutine GeneralUpdateMassBalancePatch

! ************************************************************************** !
!
! GeneralUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the General problem
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralUpdateAuxVars(realization)

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
      call GeneralUpdateAuxVarsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeneralUpdateAuxVars

! ************************************************************************** !
!
! GeneralUpdateAuxVarsPatch: Updates the auxilliary variables associated with 
!                        the General problem
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralUpdateAuxVarsPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:), gen_aux_vars_bc(:)  
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)  

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)  
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  gen_aux_vars => patch%aux%General%aux_vars
  gen_aux_vars_bc => patch%aux%General%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
    
  call GridVecGetArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  

  do ghosted_id = 1, grid%ngmax
     if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
   
    call GeneralAuxVarCompute(xx_loc_p(ghosted_id:ghosted_id), &
                       gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                       global_aux_vars(ghosted_id), &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
                       option)
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
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
          xxbc(1) = boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          xxbc(1) = xx_loc_p(ghosted_id)
      end select
      
      call GeneralAuxVarCompute(xxbc(1),gen_aux_vars_bc(sum_connection), &
                         global_aux_vars_bc(sum_connection), &
                         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                         
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo


  call GridVecRestoreArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  
  
  patch%aux%General%aux_vars_up_to_date = PETSC_TRUE

end subroutine GeneralUpdateAuxVarsPatch

! ************************************************************************** !
!
! GeneralInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralInitializeTimestep(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  call GeneralUpdateFixedAccum(realization)

end subroutine GeneralInitializeTimestep

! ************************************************************************** !
!
! GeneralUpdateSolution: Updates data in module after a successful time 
!                             step
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralUpdateSolution(realization)

  use Realization_module
  use Field_module
  use Level_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization

  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  
  field => realization%field
  
  call VecCopy(field%flow_xx,field%flow_yy,ierr)   

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GeneralUpdateSolutionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine GeneralUpdateSolution


! ************************************************************************** !
!
! GeneralUpdateSolutionPatch: Updates data in module after a successful time 
!                             step
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralUpdateSolutionPatch(realization)

  use Realization_module
    
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call GeneralUpdateMassBalancePatch(realization)
  endif

end subroutine GeneralUpdateSolutionPatch

! ************************************************************************** !
!
! GeneralUpdateFixedAccum: Updates the fixed portion of the 
!                                  accumulation term
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralUpdateFixedAccum(realization)

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
      call GeneralUpdateFixedAccumPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeneralUpdateFixedAccum

! ************************************************************************** !
!
! GeneralUpdateFixedAccumPatch: Updates the fixed portion of the 
!                                accumulation term
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralUpdateFixedAccumPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_module

  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          accum_p(:), perm_xx_loc_p(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  gen_aux_vars => patch%aux%General%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  material_parameter => patch%aux%Material%material_parameter
    
  call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)  

  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call GeneralAuxVarCompute(xx_p(local_id:local_id), &
                              gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                              global_aux_vars(ghosted_id), &
                              realization%saturation_function_array( &
                                int(icap_loc_p(ghosted_id)))%ptr, &
                              porosity_loc_p(ghosted_id), &
                              perm_xx_loc_p(ghosted_id), &                        
                              option)
    call GeneralAccumulation(gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                             global_aux_vars(ghosted_id), &
                             material_parameter%dencpr(imat), &
                             porosity_loc_p(ghosted_id), &
                             volume_p(local_id), &
                             option,accum_p(local_id:local_id)) 
  enddo

  call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)  

  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)

end subroutine GeneralUpdateFixedAccumPatch

! ************************************************************************** !
!
! GeneralNumericalJacTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralNumericalJacTest(xx,realization)

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
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  
  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  
  call VecDuplicate(xx,xx_pert,ierr)
  call VecDuplicate(xx,res,ierr)
  call VecDuplicate(xx,res_pert,ierr)
  
  call MatCreate(option%mycomm,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof,grid%nlmax*option%nflowdof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call GeneralResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
  call GridVecGetArrayF90(grid,res,vec2_p,ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    idof = icell
!    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call veccopy(xx,xx_pert,ierr)
      call vecgetarrayf90(xx_pert,vec_p,ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call vecrestorearrayf90(xx_pert,vec_p,ierr)
      call GeneralResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
      call vecgetarrayf90(res_pert,vec_p,ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values,ierr)
        endif
      enddo
      call GridVecRestoreArrayF90(grid,res_pert,vec_p,ierr)
!    enddo
  enddo
  call GridVecRestoreArrayF90(grid,res,vec2_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerASCIIOpen(option%mycomm,'numerical_jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call MatDestroy(A,ierr)
  
  call VecDestroy(xx_pert,ierr)
  call VecDestroy(res,ierr)
  call VecDestroy(res_pert,ierr)
  
end subroutine GeneralNumericalJacTest

! ************************************************************************** !
!
! GeneralAuxVarPerturb: Calculates auxilliary variables for perturbed system
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralAuxVarPerturb(gen_aux_var,global_aux_var, &
                                saturation_function,option)

  use Option_module
  use Saturation_Function_module

  implicit none

  type(option_type) :: option
  type(general_auxvar_type) :: gen_aux_var(0:)
  type(global_auxvar_type) :: global_aux_var
  type(saturation_function_type) :: saturation_function
     
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert(option%nflowdof)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  PetscInt :: idof

  select case(global_aux_var%istate)
    case(LIQUID_STATE)
       x(ONE_INTEGER) = gen_aux_var(ZERO_INTEGER)%pres(option%liquid_phase)
       x(TWO_INTEGER) = gen_aux_var(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
       x(THREE_INTEGER) = gen_aux_var(ZERO_INTEGER)%temp
       pert(ONE_INTEGER) = 1.d0
       pert(TWO_INTEGER) = perturbation_tolerance*x(TWO_INTEGER)
       pert(THREE_INTEGER) = perturbation_tolerance*x(THREE_INTEGER)
    case(GAS_STATE)
       x(ONE_INTEGER) = gen_aux_var(ZERO_INTEGER)%pres(option%gas_phase)
       x(TWO_INTEGER) = gen_aux_var(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(THREE_INTEGER) = gen_aux_var(ZERO_INTEGER)%temp
       pert(ONE_INTEGER) = 1.d0
       if (x(ONE_INTEGER) - x(TWO_INTEGER) > 1.d0) then 
         pert(TWO_INTEGER) = 1.d0
       else
         pert(TWO_INTEGER) = -1.d0
       endif
       pert(THREE_INTEGER) = perturbation_tolerance*x(THREE_INTEGER)
    case(TWO_PHASE_STATE)
       x(ONE_INTEGER) = gen_aux_var(ZERO_INTEGER)%pres(option%gas_phase)
       x(TWO_INTEGER) = gen_aux_var(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(THREE_INTEGER) = gen_aux_var(ZERO_INTEGER)%sat(option%gas_phase)
       pert(ONE_INTEGER) = 1.d0
       if (x(ONE_INTEGER) - x(TWO_INTEGER) > 1.d0) then 
         pert(TWO_INTEGER) = 1.d0
       else
         pert(TWO_INTEGER) = -1.d0
       endif
       if (x(THREE_INTEGER) > 0.5d0) then 
         pert(THREE_INTEGER) = -perturbation_tolerance*x(THREE_INTEGER)
       else
         pert(THREE_INTEGER) = perturbation_tolerance*x(THREE_INTEGER)
       endif
  end select
  
  do idof = 1, option%nflowdof
    gen_aux_var(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    call GeneralAuxVarCompute(x_pert,gen_aux_var(idof),global_aux_var, &
                              saturation_function,0.d0,0.d0,option)
  enddo
  
end subroutine GeneralAuxVarPerturb
  
! ************************************************************************** !
!
! GeneralAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !  
subroutine GeneralAccumulation(gen_aux_var,global_aux_var,dencpr,por,vol, &
                               option,Res)

  use Option_module
  
  implicit none

  type(general_auxvar_type) :: gen_aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  PetscReal :: vol, por, dencpr
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: pv_over_t
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  pv_over_t = por * vol / option%flow_dt
  
  ! accumulation term units = kmol/s
  Res = 0.d0
  do icomp = 1, option%nflowspec
    do iphase = 1, option%nphase
      Res(icomp) = Res(icomp) + gen_aux_var%sat(iphase) * &
                                gen_aux_var%den(iphase) * &
                                gen_aux_var%xmol(icomp,iphase)
    enddo
  enddo
  
  do icomp = 2, option%nflowspec
    Res(ONE_INTEGER) = Res(ONE_INTEGER) + Res(icomp)
  enddo
  
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * pv_over_t
  
  do iphase = 1, option%nphase
    Res(energy_id) = Res(energy_id) + gen_aux_var%sat(iphase) * &
                                      gen_aux_var%den(iphase) * &
                                      gen_aux_var%u(iphase)
  enddo

  Res(energy_id) = Res(energy_id) * pv_over_t + &
                                    (1.d0 - por) * dencpr * gen_aux_var%temp

end subroutine GeneralAccumulation
  
! ************************************************************************** !
!
! GeneralAccumDerivative: Computes derivatives of the accumulation 
!                                 term for the Jacobian
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralAccumDerivative(gen_aux_var,global_aux_var,dencpr,por,vol, &
                                  option,J)

  use Option_module
  use Saturation_Function_module
  
  implicit none

  type(general_auxvar_type) :: gen_aux_var(0:)
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: dencpr, vol, por
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  call GeneralAccumulation(gen_aux_var(ZERO_INTEGER),global_aux_var,dencpr, &
                           por,vol,option,res)
                           
  do idof = 1, option%nflowdof
    call GeneralAccumulation(gen_aux_var(idof),global_aux_var,dencpr, &
                             por,vol,option,res_pert)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/gen_aux_var(idof)%pert
    enddo !irow
  enddo ! idof

end subroutine GeneralAccumDerivative

! ************************************************************************** !
!
! GeneralFlux: Computes the internal flux terms for the residual
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** ! 
subroutine GeneralFlux(gen_aux_var_up,global_aux_var_up, &
                       por_up,sir_up,dd_up,perm_up, &
                       gen_aux_var_dn,global_aux_var_dn, &
                       por_dn,sir_dn,dd_dn,perm_dn, &
                       area,dist_gravity,upweight, &
                       option,v_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(general_auxvar_type) :: gen_aux_var_up, gen_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: upweight
  PetscReal :: Res(option%nflowdof)

  PetscReal :: upweight_adj
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: fmw_phase(option%nphase)
  PetscReal :: xmol(option%nflowspec)
  PetscReal :: density_ave
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist
  PetscReal :: delta_pressure
  PetscReal :: gravity
  PetscReal :: ukvr, mole_flux, q
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  fmw_phase(option%liquid_phase) = FMWH2O
  fmw_phase(option%gas_phase) = FMWAIR

  perm_ave_over_dist = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  Res = 0.d0
  do iphase = 1, option%nphase
    
    ! using residual saturation cannot be correct! - geh
    if (gen_aux_var_up%sat(iphase) > sir_up .or. &
        gen_aux_var_dn%sat(iphase) > sir_dn) then
      upweight_adj = upweight
      if (gen_aux_var_up%sat(iphase) < eps) then 
        upweight_adj=0.d0
      else if (gen_aux_var_dn%sat(iphase) < eps) then 
        upweight_adj=1.d0
      endif    
      density_ave = upweight_adj*gen_aux_var_up%den(iphase)+ &
                    (1.D0-upweight_adj)*gen_aux_var_dn%den(iphase)
      H_ave = upweight_adj*gen_aux_var_up%H(iphase)+ &
              (1.D0-upweight_adj)*gen_aux_var_dn%H(iphase)

      gravity = (upweight_adj*gen_aux_var_up%den(iphase) + &
                (1.D0-upweight)*gen_aux_var_dn%den(iphase)) &
                * fmw_phase(iphase) * dist_gravity 

      delta_pressure = gen_aux_var_up%pres(iphase) - &
             gen_aux_var_dn%pres(iphase) + &
             gravity

      if (delta_pressure >= 0.D0) then
        ukvr = gen_aux_var_up%kvr(iphase)
        xmol(:) = gen_aux_var_up%xmol(:,iphase)
      else
        ukvr = gen_aux_var_dn%kvr(iphase)
        xmol(:) = gen_aux_var_dn%xmol(:,iphase)
      endif      

      if (ukvr > floweps) then
        v_darcy(iphase) = perm_ave_over_dist * ukvr * delta_pressure
        q = v_darcy(iphase) * area
        mole_flux = q*density_ave       
        do icomp = 1, option%nflowspec
          Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
        enddo
        
        Res(energy_id) = Res(energy_id) + mole_flux * &
                                          gen_aux_var_dn%H(iphase)
        
      endif                   
    endif ! sat > eps
  enddo
  
  do icomp = 2, option%nflowspec
    Res(ONE_INTEGER) = Res(ONE_INTEGER) + Res(icomp)
  enddo

end subroutine GeneralFlux

! ************************************************************************** !
!
! GeneralFluxDerivative: Computes the derivatives of the internal flux terms
!                         for the Jacobian
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** ! 
subroutine GeneralFluxDerivative(gen_aux_var_up,global_aux_var_up,por_up, &
                                  sir_up,dd_up,perm_up, &
                                  gen_aux_var_dn,global_aux_var_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn, &
                                  area,dist_gravity,upweight, &
                                  option,Jup,Jdn)
  use Option_module 
  
  implicit none
  
  type(general_auxvar_type) :: gen_aux_var_up(0:), gen_aux_var_dn(0:)
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: area
  PetscReal :: dist_gravity
  PetscReal :: upweight
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  call GeneralFlux(gen_aux_var_up(ZERO_INTEGER),global_aux_var_up, &
                   por_up,sir_up,dd_up,perm_up, &
                   gen_aux_var_dn(ZERO_INTEGER),global_aux_var_dn, &
                   por_dn,sir_dn,dd_dn,perm_dn, &
                   area,dist_gravity,upweight, &
                   option,v_darcy,res)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_aux_var_up(idof),global_aux_var_up, &
                     por_up,sir_up,dd_up,perm_up, &
                     gen_aux_var_dn(ZERO_INTEGER),global_aux_var_dn, &
                     por_dn,sir_dn,dd_dn,perm_dn, &
                     area,dist_gravity,upweight, &
                     option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res(irow))/gen_aux_var_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_aux_var_up(ZERO_INTEGER),global_aux_var_up, &
                     por_up,sir_up,dd_up,perm_up, &
                     gen_aux_var_dn(idof),global_aux_var_dn, &
                     por_dn,sir_dn,dd_dn,perm_dn, &
                     area,dist_gravity,upweight, &
                     option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_aux_var_dn(idof)%pert
    enddo !irow
  enddo ! idof

end subroutine GeneralFluxDerivative

! ************************************************************************** !
!
! GeneralBCFlux: Computes the boundary flux terms for the residual
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** ! 
subroutine GeneralBCFlux(ibndtype,aux_vars, &
                         gen_aux_var_up,global_aux_var_up, &
                         gen_aux_var_dn,global_aux_var_dn, &
                         por_dn,sir_dn,dd_up,perm_dn, &
                         area,dist_gravity,option,v_darcy,Res)
  use Option_module                              
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(general_auxvar_type) :: gen_aux_var_up, gen_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn(option%nphase)
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn
  PetscReal :: v_darcy(option%nphase), area
  PetscReal :: dist_gravity
  PetscReal :: Res(1:option%nflowdof)

  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase

  PetscReal :: fmw_phase(option%nphase)
  PetscReal :: xmol(option%nflowspec)  
  PetscReal :: density_ave
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist
  PetscReal :: delta_pressure
  PetscReal :: gravity
  PetscReal :: ukvr, mole_flux, q
  
  PetscInt :: bc_type
  PetscInt :: idof
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  fmw_phase(option%liquid_phase) = FMWH2O
  fmw_phase(option%gas_phase) = FMWAIR

  Res = 0.d0

  do iphase = 1, option%nphase
  
    select case(iphase)
      case(LIQUID_PHASE)
        bc_type = ibndtype(GENERAL_LIQUID_PRESSURE_DOF)
      case(GAS_PHASE)
        bc_type = ibndtype(GENERAL_GAS_PRESSURE_DOF)
    end select

    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)

        if (bc_type == CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = GENERAL_LIQUID_CONDUCTANCE_DOF
            case(GAS_PHASE)
              idof = GENERAL_GAS_CONDUCTANCE_DOF
          end select        
          perm_ave_over_dist = aux_vars(idof)
        else
          perm_ave_over_dist = perm_dn / dd_up
        endif
        
          
        ! using residual saturation cannot be correct! - geh
        if (gen_aux_var_up%sat(iphase) > sir_dn(iphase) .or. &
            gen_aux_var_dn%sat(iphase) > sir_dn(iphase)) then
          upweight = upweight
          if (gen_aux_var_up%sat(iphase) < eps) then 
            upweight=0.d0
          else if (gen_aux_var_dn%sat(iphase) < eps) then 
            upweight=1.d0
          endif    
          density_ave = upweight*gen_aux_var_up%den(iphase)+ &
                        (1.D0-upweight)*gen_aux_var_dn%den(iphase)
          H_ave = upweight*gen_aux_var_up%H(iphase)+ &
                  (1.D0-upweight)*gen_aux_var_dn%H(iphase)

          gravity = (upweight*gen_aux_var_up%den(iphase) + &
                    (1.D0-upweight)*gen_aux_var_dn%den(iphase)) &
                    * fmw_phase(iphase) * dist_gravity 

          delta_pressure = gen_aux_var_up%pres(iphase) - &
                 gen_aux_var_dn%pres(iphase) + &
                 gravity

          if (bc_type == SEEPAGE_BC .or. &
              bc_type == CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                global_aux_var_up%pres(iphase)-option%reference_pressure < eps) then
              delta_pressure = 0.d0
            endif
          endif
            
          if (delta_pressure >= 0.D0) then
            ukvr = gen_aux_var_up%kvr(iphase)
            xmol(:) = gen_aux_var_up%xmol(:,iphase)
          else
            ukvr = gen_aux_var_dn%kvr(iphase)
            xmol(:) = gen_aux_var_dn%xmol(:,iphase)
          endif      

          if (ukvr > floweps) then
            v_darcy(iphase) = perm_ave_over_dist * ukvr * delta_pressure
          endif                   
        endif ! sat > eps

      case(NEUMANN_BC)
        select case(iphase)
          case(LIQUID_PHASE)
            idof = GENERAL_LIQUID_FLUX_DOF
          case(GAS_PHASE)
            idof = GENERAL_GAS_FLUX_DOF
        end select
      
        if (dabs(aux_vars(idof)) > floweps) then
          v_darcy(iphase) = aux_vars(idof)
          if (v_darcy(iphase) > 0.d0) then 
            density_ave = global_aux_var_up%den(iphase)
          else 
            density_ave = global_aux_var_dn%den(iphase)
          endif 
        endif
    end select

    q = v_darcy(iphase) * area
    mole_flux = q*density_ave       
    do icomp = 1, option%nflowspec
      Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
    enddo
    
    Res(energy_id) = Res(energy_id) + mole_flux * &
                                      gen_aux_var_dn%H(iphase)

  enddo
  
  do icomp = 2, option%nflowspec
    Res(ONE_INTEGER) = Res(ONE_INTEGER) + Res(icomp)
  enddo

end subroutine GeneralBCFlux

! ************************************************************************** !
!
! GeneralBCFluxDerivative: Computes the derivatives of the boundary flux terms
!                          for the Jacobian
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** ! 
subroutine GeneralBCFluxDerivative(ibndtype,aux_vars, &
                                    gen_aux_var_up, &
                                    global_aux_var_up, &
                                    gen_aux_var_dn,global_aux_var_dn, &
                                    por_dn,sir_dn,dd_up,perm_dn, &
                                    area,dist_gravity,option,Jdn)

  use Option_module 
  
  implicit none

  PetscInt :: ibndtype(:)
  PetscReal :: aux_vars(:) ! from aux_real_var array
  type(general_auxvar_type) :: gen_aux_var_up, gen_aux_var_dn(0:)
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn(option%nphase)
  PetscReal :: por_dn,perm_dn
  PetscReal :: area
  PetscReal :: dist_gravity
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  call GeneralBCFlux(ibndtype,aux_vars, &
                     gen_aux_var_up,global_aux_var_up, &
                     gen_aux_var_dn(ZERO_INTEGER),global_aux_var_dn, &
                     por_dn,sir_dn,dd_up,perm_dn, &
                     area,dist_gravity,option,v_darcy,res)                     
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralBCFlux(ibndtype,aux_vars, &
                       gen_aux_var_up,global_aux_var_up, &
                       gen_aux_var_dn(idof),global_aux_var_dn, &
                       por_dn,sir_dn,dd_up,perm_dn, &
                       area,dist_gravity,option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_aux_var_dn(idof)%pert
    enddo !irow
  enddo ! idof

end subroutine GeneralBCFluxDerivative

! ************************************************************************** !
!
! GeneralResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralResidual(snes,xx,r,realization,ierr)

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
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       Vec :: vec
     end subroutine samrpetscobjectstateincrease
     
     subroutine SAMRCoarsenFaceFluxes(p_application, vec, ierr)
       implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       PetscFortranAddr :: p_application
       Vec :: vec
       PetscErrorCode :: ierr
     end subroutine SAMRCoarsenFaceFluxes
  end interface

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
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option
  
  ! Communication -----------------------------------------
  ! These 3 must be called before GeneralUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  
  ! pass #1 for internal and boundary flux terms
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GeneralResidualPatch1(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! now coarsen all face fluxes in case we are using SAMRAI to 
  ! ensure consistent fluxes at coarse-fine interfaces
  if(option%use_samr) then
     call SAMRCoarsenFaceFluxes(discretization%amrgrid%p_application, field%flow_face_fluxes, ierr)

     cur_level => realization%level_list%first
     do
        if (.not.associated(cur_level)) exit
        cur_patch => cur_level%patch_list%first
        do
           if (.not.associated(cur_patch)) exit
           realization%patch => cur_patch
           call GeneralResidualFluxContribPatch(r,realization,ierr)
           cur_patch => cur_patch%next
        enddo
        cur_level => cur_level%next
     enddo
  endif


  ! pass #2 for everything else
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GeneralResidualPatch2(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (discretization%itype==AMR_GRID) then
     call samrpetscobjectstateincrease(r)
  endif
   
  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
end subroutine GeneralResidual

! ************************************************************************** !
!
! GeneralResidualfuxContribsPatch: should be called only for SAMR
! author: Bobby Philip
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralResidualFluxContribPatch(r,realization,ierr)
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Debug_module
  
  implicit none

  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscErrorCode :: ierr

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: face_fluxes_p(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscInt :: axis, nlx, nly, nlz
  PetscInt :: iconn, i, j, k
  PetscInt :: xup_id, xdn_id, yup_id, ydn_id, zup_id, zdn_id

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,r, r_p, ierr)

  do axis=0,2  
     call GridVecGetArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
  enddo

  nlx = grid%structured_grid%nlx  
  nly = grid%structured_grid%nly  
  nlz = grid%structured_grid%nlz 
  
  iconn=0
  do k=1,nlz
     do j=1,nly
        do i=1,nlx
           iconn=iconn+1
           xup_id = ((k-1)*nly+j-1)*(nlx+1)+i
           xdn_id = xup_id+1
           yup_id = ((k-1)*(nly+1)+(j-1))*nlx+i
           ydn_id = yup_id+nlx
           zup_id = ((k-1)*nly+(j-1))*nlx+i
           zdn_id = zup_id+nlx*nly

           r_p(iconn) = r_p(iconn)+fluxes(0)%flux_p(xdn_id)-fluxes(0)%flux_p(xup_id) &
                                  +fluxes(1)%flux_p(ydn_id)-fluxes(1)%flux_p(yup_id) &
                                  +fluxes(2)%flux_p(zdn_id)-fluxes(2)%flux_p(zup_id)

        enddo
     enddo
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
!!$
!!$  do axis=0,2  
!!$     call GridVecRestoreArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
!!$  enddo

end subroutine GeneralResidualFluxContribPatch

! ************************************************************************** !
!
! GeneralResidualPatch1: Computes the residual equation 
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralResidualPatch1(snes,xx,r,realization,ierr)

  use water_eos_module

  use Connection_module
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use Material_Aux_module
  
  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
     implicit none
     
#include "finclude/petscsysdef.h"
     
     PetscFortranAddr :: p_patch
     PetscInt :: axis,dim
   end function samr_patch_at_bc
  end interface

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
  PetscReal, pointer :: icap_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:)
  PetscInt :: icap_up, icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy(realization%option%nphase)


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:), gen_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: iphase
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  gen_aux_vars => patch%aux%General%aux_vars
  gen_aux_vars_bc => patch%aux%General%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  call GeneralUpdateAuxVarsPatch(realization)
  patch%aux%General%aux_vars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
    call GeneralZeroMassBalDeltaPatch(realization)
  endif

! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
  !print *,' Finished scattering non deriv'

  if (option%use_samr) then
     do axis=0,2  
        call GridVecGetArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
     enddo
  endif

  r_p = 0.d0

  if (option%use_samr) then
    nlx = grid%structured_grid%nlx  
    nly = grid%structured_grid%nly  
    nlz = grid%structured_grid%nlz 

    ngx = grid%structured_grid%ngx   
    ngxy = grid%structured_grid%ngxy

    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 0, 0)==1) nlx = nlx-1
    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 0, 1)==1) nlx = nlx-1
    
    max_x_conn = (nlx+1)*nly*nlz
    ! reinitialize nlx
    nlx = grid%structured_grid%nlx  

    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 1, 0)==1) nly = nly-1
    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 1, 1)==1) nly = nly-1
    
    max_y_conn = max_x_conn + nlx*(nly+1)*nlz

    ! reinitialize nly
    nly = grid%structured_grid%nly  
  endif

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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
        
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(3,iconn))

      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
   
      call GeneralFlux(gen_aux_vars(ZERO_INTEGER,ghosted_id_up), &
                        global_aux_vars(ghosted_id_up), &
                          porosity_loc_p(ghosted_id_up), &
                          material_parameter%sir(1,icap_up), &
                          dd_up,perm_up, &
                        gen_aux_vars(ZERO_INTEGER,ghosted_id_dn), &
                        global_aux_vars(ghosted_id_dn), &
                          porosity_loc_p(ghosted_id_dn), &
                          material_parameter%sir(1,icap_dn), &
                          dd_dn,perm_dn, &
                        cur_connection_set%area(iconn),distance_gravity, &
                        upweight,option,v_darcy,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy
      
      if (option%use_samr) then
        if (sum_connection <= max_x_conn) then
          direction = 0
          if(mod(mod(ghosted_id_dn,ngxy),ngx).eq.0) then
             flux_id = ((ghosted_id_dn/ngxy)-1)*(nlx+1)*nly + &
                       ((mod(ghosted_id_dn,ngxy))/ngx-1)*(nlx+1)
          else
             flux_id = ((ghosted_id_dn/ngxy)-1)*(nlx+1)*nly + &
                       ((mod(ghosted_id_dn,ngxy))/ngx-1)*(nlx+1)+ &
                       mod(mod(ghosted_id_dn,ngxy),ngx)-1
          endif

        else if (sum_connection <= max_y_conn) then
          direction = 1
          flux_id = ((ghosted_id_dn/ngxy)-1)*nlx*(nly+1) + &
                    ((mod(ghosted_id_dn,ngxy))/ngx-1)*nlx + &
                    mod(mod(ghosted_id_dn,ngxy),ngx)-1
        else
          direction = 2
          flux_id = ((ghosted_id_dn/ngxy)-1)*nlx*nly &
                   +((mod(ghosted_id_dn,ngxy))/ngx-1)*nlx &
                   +mod(mod(ghosted_id_dn,ngxy),ngx)-1
        endif
        fluxes(direction)%flux_p(flux_id) = Res(1)
      endif
      
      if(.not.option%use_samr) then
         
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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*dabs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))  

      call GeneralBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                gen_aux_vars_bc(sum_connection), &
                                global_aux_vars_bc(sum_connection), &
                                gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                material_parameter%sir(:,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn, &
                                cur_connection_set%area(iconn), &
                                distance_gravity,option, &
                                v_darcy,Res)
      patch%boundary_velocities(:,sum_connection) = v_darcy
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_aux_vars_bc(sum_connection)%mass_balance_delta(1,iphase) = &
          global_aux_vars_bc(sum_connection)%mass_balance_delta(1,iphase) - Res(1)
        ! contribution to internal 
!        global_aux_vars(ghosted_id)%mass_balance_delta(1) = &
!          global_aux_vars(ghosted_id)%mass_balance_delta(1) + Res(1)
      endif

      if (option%use_samr) then
         direction =  (boundary_condition%region%faces(iconn)-1)/2

         ! the ghosted_id gives the id of the cell. Since the
         ! flux_id is based on the ghosted_id of the downwind
         ! cell this has to be adjusted in the case of the east, 
         ! north and top faces before the flux_id is computed
         select case(boundary_condition%region%faces(iconn)) 
           case(WEST_FACE)
              flux_id = ((ghosted_id/ngxy)-1)*(nlx+1)*nly + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*(nlx+1)+ &
                        mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = Res(1)
           case(EAST_FACE)
              ghosted_id = ghosted_id+1
              flux_id = ((ghosted_id/ngxy)-1)*(nlx+1)*nly + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*(nlx+1)
              fluxes(direction)%flux_p(flux_id) = -Res(1)
           case(SOUTH_FACE)
              flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                        mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = Res(1)
           case(NORTH_FACE)
              ghosted_id = ghosted_id+ngx
              flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                        mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = -Res(1)
           case(BOTTOM_FACE)
              flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                       +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                       +mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = Res(1)
           case(TOP_FACE)
              ghosted_id = ghosted_id+ngxy
              flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                       +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                       +mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = -Res(1)
         end select

!         fluxes(direction)%flux_p(flux_id) = Res(1)

      else
         r_p(local_id)= r_p(local_id) - Res(1)
      endif

   enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)

end subroutine GeneralResidualPatch1

! ************************************************************************** !
!
! GeneralResidualPatch2: Computes the residual equation 
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralResidualPatch2(snes,xx,r,realization,ierr)

  use water_eos_module

  use Connection_module
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use Material_Aux_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscErrorCode :: ierr

  PetscInt :: i, imat
  PetscInt :: local_id, ghosted_id

  PetscReal, pointer :: accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:)

  PetscReal :: qsrc, qsrc_mol
  PetscReal :: Res(realization%option%nflowdof)


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:), gen_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(material_parameter_type), pointer :: material_parameter
  PetscInt :: iconn
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  gen_aux_vars => patch%aux%General%aux_vars
  gen_aux_vars_bc => patch%aux%General%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  material_parameter => patch%aux%Material%material_parameter

! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)

  ! Accumulation terms ------------------------------------
  if (.not.option%steady_state) then

    r_p = r_p - accum_p

    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      call GeneralAccumulation(gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                material_parameter%dencpr(imat), &
                                volume_p(local_id), &
                                option,Res) 
      r_p(local_id) = r_p(local_id) + Res(1)
    enddo

  endif

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    qsrc = source_sink%flow_condition%rate%dataset%cur_value(1)
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O ! kg/sec -> kmol/sec
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc_mol = qsrc*global_aux_vars(ghosted_id)%den(1) ! den = kmol/m^3
      end select
!      if (option%compute_mass_balance_new) then
        ! need to added global aux_var for src/sink
!        global_aux_vars_ss(ghosted_id)%mass_balance_delta(1) = &
!          global_aux_vars_ss(ghosted_id)%mass_balance_delta(1) - qsrc_kg
!      endif
      r_p(local_id) = r_p(local_id) - qsrc_mol
    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%General%inactive_cells_exist) then
    do i=1,patch%aux%General%n_zero_rows
      r_p(patch%aux%General%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

end subroutine GeneralResidualPatch2

! ************************************************************************** !
!
! GeneralJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  interface
     subroutine SAMRSetCurrentJacobianPatch(mat,patch) 
#include "finclude/petscsysdef.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
       
       Mat :: mat
       PetscFortranAddr :: patch
     end subroutine SAMRSetCurrentJacobianPatch
  end interface

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
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      grid => cur_patch%grid
      ! need to set the current patch in the Jacobian operator
      ! so that entries will be set correctly
      if(option%use_samr) then
         call SAMRSetCurrentJacobianPatch(J, grid%structured_grid%p_samr_patch)
      endif

      call GeneralJacobianPatch1(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! pass #2 for everything else
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
      if(option%use_samr) then
         call SAMRSetCurrentJacobianPatch(J, grid%structured_grid%p_samr_patch)
      endif

      call GeneralJacobianPatch2(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (realization%debug%matview_Jacobian) then
#if 1  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rjacobian.out', &
                              viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rjacobian.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
#endif    
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

end subroutine GeneralJacobian
                
! ************************************************************************** !
!
! GeneralJacobianPatch1: Computes the Jacobian
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralJacobianPatch1(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module

  use Connection_module
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscReal, pointer :: porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: icap_loc_p(:)
  PetscInt :: icap,icap_up,icap_dn
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
  type(field_type), pointer :: field 
  type(material_parameter_type), pointer :: material_parameter
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:), gen_aux_vars_bc(:)
  type(general_auxvar_type) :: test
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  gen_aux_vars => patch%aux%General%aux_vars
  gen_aux_vars_bc => patch%aux%General%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
  
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

      if (patch%imat(ghosted_id_up) <= 0 .or. &
          patch%imat(ghosted_id_dn) <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
    
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(3,iconn))
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
                              
      call GeneralFluxDerivative(gen_aux_vars(:,ghosted_id_up), &
                                  global_aux_vars(ghosted_id_up), &
                                    porosity_loc_p(ghosted_id_up), &
                                    material_parameter%sir(1,icap_up), &
                                    dd_up,perm_up, &
                                  gen_aux_vars(:,ghosted_id_dn), &
                                  global_aux_vars(ghosted_id_dn), &
                                    porosity_loc_p(ghosted_id_dn), &
                                    material_parameter%sir(1,icap_dn), &
                                    dd_dn,perm_dn, &
                                  cur_connection_set%area(iconn),distance_gravity, &
                                  upweight,option,&
                                  Jup,Jdn)
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

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_flux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*dabs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))  

      call GeneralBCFluxDerivative(boundary_condition%flow_condition%itype, &
                                  boundary_condition%flow_aux_real_var(:,iconn), &
                                  gen_aux_vars_bc(sum_connection), &
                                  global_aux_vars_bc(sum_connection), &
                                  gen_aux_vars(:,ghosted_id), &
                                  global_aux_vars(ghosted_id), &
                                  porosity_loc_p(ghosted_id), &
                                  material_parameter%sir(:,icap_dn), &
                                  cur_connection_set%dist(0,iconn),perm_dn, &
                                  cur_connection_set%area(iconn), &
                                  distance_gravity,option, &
                                  Jdn)

      Jdn = -Jdn
        call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                               ADD_VALUES,ierr) 
    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)

end subroutine GeneralJacobianPatch1

! ************************************************************************** !
!
! GeneralJacobianPatch2: Computes the Jacobian
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralJacobianPatch2(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module

  use Connection_module
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_module
    
  implicit none

  interface
     subroutine SAMRSetJacobianSourceOnPatch(which_pc, index, val, p_application, p_patch) 
#include "finclude/petscsysdef.h"

       PetscInt :: which_pc
       PetscInt :: index
       PetscReal :: val
       PetscFortranAddr :: p_application
       PetscFortranAddr :: p_patch
     end subroutine SAMRSetJacobianSourceOnPatch

     subroutine SAMRSetJacobianSrcCoeffsOnPatch(which_pc, p_application, p_patch) 
#include "finclude/petscsysdef.h"

       PetscInt :: which_pc
       PetscFortranAddr :: p_application
       PetscFortranAddr :: p_patch
     end subroutine SAMRSetJacobianSrcCoeffsOnPatch
  end interface

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), icap_loc_p(:)
  PetscReal :: qsrc
  PetscInt :: icap, imat
  PetscInt :: local_id, ghosted_id
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(material_parameter_type), pointer :: material_parameter
  PetscInt :: flow_pc
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  gen_aux_vars => patch%aux%General%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  material_parameter => patch%aux%Material%material_parameter

  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
  
  if (.not.option%steady_state) then

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    icap = int(icap_loc_p(ghosted_id))
    call GeneralAccumDerivative(gen_aux_vars(:,ghosted_id), &
                              global_aux_vars(ghosted_id), &
                              material_parameter%dencpr(imat), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option, &
                              Jup) 
      call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                             ADD_VALUES,ierr)
  enddo
  endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_accum.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    qsrc = source_sink%flow_condition%rate%dataset%cur_value(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
      
      Jup = 0.d0
      select case(source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
#ifdef GENERAL2
          Jup(1,1) = -qsrc*gen_aux_vars(ghosted_id)%dden_dp*FMWH2O
#endif
      end select
      call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)  

    enddo
    source_sink => source_sink%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
  if (patch%aux%General%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%General%n_zero_rows, &
                          patch%aux%General%zero_rows_local_ghosted, &
                          qsrc,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif

  if(option%use_samr) then
     flow_pc = 0
     call SAMRSetJacobianSrcCoeffsOnPatch(flow_pc, &
          realization%discretization%amrgrid%p_application, grid%structured_grid%p_samr_patch)
  endif

end subroutine GeneralJacobianPatch2

! ************************************************************************** !
!
! GeneralCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralCreateZeroArray(patch,option)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id

  type(grid_type), pointer :: grid
  PetscInt :: flag
  PetscInt :: n_zero_rows
  PetscInt, pointer :: zero_rows_local(:)
  PetscInt, pointer :: zero_rows_local_ghosted(:)
  PetscErrorCode :: ierr
    
  flag = 0
  grid => patch%grid
  
  n_zero_rows = 0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      n_zero_rows = n_zero_rows + option%nflowdof
    endif
  enddo

  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))

  zero_rows_local = 0
  zero_rows_local_ghosted = 0
  ncount = 0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      do idof = 1, option%nflowdof
        ncount = ncount + 1
        zero_rows_local(ncount) = (local_id-1)*option%nflowdof+idof
        zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%nflowdof+idof-1
      enddo
    endif
  enddo

  patch%aux%General%zero_rows_local => zero_rows_local
  patch%aux%General%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%General%n_zero_rows = n_zero_rows
  
  if(.not. (option%use_samr)) then
     call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX, &
                        option%mycomm,ierr)
     if (flag > 0) patch%aux%General%inactive_cells_exist = PETSC_TRUE
     
     if (ncount /= n_zero_rows) then
        print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
        stop
     endif
  endif

end subroutine GeneralCreateZeroArray

! ************************************************************************** !
!
! GeneralMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralMaxChange(realization)

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

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

end subroutine GeneralMaxChange

! ************************************************************************** !
!
! GeneralGetTecplotHeader: Returns General Lite contribution to 
!                               Tecplot file header
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
function GeneralGetTecplotHeader(realization,icolumn)
  
  use Realization_module
  use Option_module
  use Field_module
    
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: GeneralGetTecplotHeader
  type(realization_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  

  option => realization%option
  field => realization%field
  
  string = ''
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')') 
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-sl"'')') icolumn
  else
    write(string2,'('',"sl"'')') 
  endif
  string = trim(string) // trim(string2)
 
  GeneralGetTecplotHeader = string

end function GeneralGetTecplotHeader

! ************************************************************************** !
!
! GeneralDestroy: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralDestroy(realization)

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
      call GeneralDestroyPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeneralDestroy

! ************************************************************************** !
!
! GeneralDestroyPatch: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralDestroyPatch(realization)

  use Realization_module

  implicit none

  type(realization_type) :: realization
  
  ! place anything that needs to be freed here.

end subroutine GeneralDestroyPatch

end module General_module
