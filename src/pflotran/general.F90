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
         GeneralDestroy, GeneralSetPlotVariables, &
         GeneralCheckUpdatePre, GeneralCheckUpdatePost

contains

! ************************************************************************** !
!
! GeneralTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralTimeCut(realization)
 
  use Realization_class
  use Option_module
  use Field_module
  use Level_module
  use Patch_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  
  PetscErrorCode :: ierr
  PetscInt :: local_id

  option => realization%option
  field => realization%field

  call VecCopy(field%flow_yy,field%flow_xx,ierr)
  call GeneralInitializeTimestep(realization)  
  
  ! loop over patches
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GeneralTimeCutPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo  
 
end subroutine GeneralTimeCut

! ************************************************************************** !
!
! GeneralTimeCutPatch: Resets arrays for time step cut
! author: Glenn Hammond
! date: 09/26/11
!
! ************************************************************************** !
subroutine GeneralTimeCutPatch(realization)
 
  use Realization_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  
  PetscInt :: local_id, ghosted_id
  PetscReal, pointer :: iphas_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  global_aux_vars => patch%aux%Global%aux_vars
  
  ! restore stored state
  call VecGetArrayReadF90(field%iphas_loc,iphas_loc_p, ierr)
  do ghosted_id = 1, grid%ngmax
    global_aux_vars(ghosted_id)%istate = int(iphas_loc_p(ghosted_id))
  enddo
  call VecRestoreArrayReadF90(field%iphas_loc,iphas_loc_p, ierr)  
 
end subroutine GeneralTimeCutPatch

! ************************************************************************** !
!
! GeneralSetup: 
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralSetup(realization)

  use Realization_class
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
  
  call GeneralSetPlotVariables(realization)  

end subroutine GeneralSetup

! ************************************************************************** !
!
! GeneralSetupPatch: Creates arrays for auxiliary variables
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralSetupPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Material_module
  use Material_Aux_module
  use Fluid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, idof
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:) ! extra index for derivatives
  type(general_auxvar_type), pointer :: gen_aux_vars_bc(:)
  type(general_auxvar_type), pointer :: gen_aux_vars_ss(:)
  type(material_type), pointer :: material
  type(fluid_property_type), pointer :: cur_fluid_property  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%General => GeneralAuxCreate(option)
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
      ! kg rock/m^3 rock * J/kg rock-K * 1.e-6 MJ/J
      patch%aux%Material%material_parameter% &
        dencpr(realization%saturation_function_array(i)%ptr%id) = &
        realization%material_property_array(i)%ptr%rock_density * &
        realization%material_property_array(i)%ptr%specific_heat * option%scale
    endif
  enddo

  ! allocate aux_var data structures for all grid cells  
  allocate(gen_aux_vars(0:option%nflowdof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      call GeneralAuxVarInit(gen_aux_vars(idof,ghosted_id),option)
    enddo
  enddo
  patch%aux%General%aux_vars => gen_aux_vars
  patch%aux%General%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! aux_var data structures for them 
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  if (sum_connection > 0) then
    allocate(gen_aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_aux_vars_bc(iconn),option)
    enddo
    patch%aux%General%aux_vars_bc => gen_aux_vars_bc
  endif
  patch%aux%General%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! aux_var data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    allocate(gen_aux_vars_ss(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_aux_vars_ss(iconn),option)
    enddo
    patch%aux%General%aux_vars_ss => gen_aux_vars_ss
  endif
  patch%aux%General%num_aux_ss = sum_connection
    
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call GeneralCreateZeroArray(patch,option)
  
  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    patch%aux%General%general_parameter% &
      diffusion_coefficient(cur_fluid_property%phase_id) = &
        cur_fluid_property%diffusion_coefficient
    cur_fluid_property => cur_fluid_property%next
  enddo  

end subroutine GeneralSetupPatch

! ************************************************************************** !
!
! GeneralSetup: 
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralComputeMassBalance(realization,mass_balance)

  use Realization_class
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
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralComputeMassBalancePatch(realization,mass_balance)
 
  use Realization_class
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

  call VecGetArrayReadF90(field%volume,volume_p,ierr)
  call VecGetArrayReadF90(field%porosity_loc,porosity_loc_p,ierr)

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

  call VecRestoreArrayReadF90(field%volume,volume_p,ierr)
  call VecRestoreArrayReadF90(field%porosity_loc,porosity_loc_p,ierr)
  
end subroutine GeneralComputeMassBalancePatch

! ************************************************************************** !
!
! GeneralZeroMassBalDeltaPatch: Zeros mass balance delta array
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralZeroMassBalDeltaPatch(realization)
 
  use Realization_class
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
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralUpdateMassBalancePatch(realization)
 
  use Realization_class
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
! GeneralUpdateAuxVars: Updates the auxiliary variables associated with 
!                        the General problem
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralUpdateAuxVars(realization,update_state)

  use Realization_class
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscBool :: update_state
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch             
      call GeneralUpdateAuxVarsPatch(realization,update_state)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeneralUpdateAuxVars

! ************************************************************************** !
!
! GeneralUpdateAuxVarsPatch: Updates the auxiliary variables associated with 
!                        the General problem
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralUpdateAuxVarsPatch(realization,update_state)

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  
  implicit none

  type(realization_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:), gen_aux_vars_bc(:)  
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)  

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: ghosted_start, ghosted_end
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:)
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
    
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecGetArrayReadF90(field%perm_xx_loc,perm_xx_loc_p,ierr)
  call VecGetArrayReadF90(field%porosity_loc,porosity_loc_p,ierr)  

  do ghosted_id = 1, grid%ngmax
     if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    call GeneralAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                       global_aux_vars(ghosted_id), &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
                       option)
    if (update_state) then
      call GeneralAuxVarUpdateState(xx_loc_p(ghosted_start:ghosted_end), &
                                    gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                                    global_aux_vars(ghosted_id), &
                                    patch%saturation_function_array( &
                                      patch%sat_func_id(ghosted_id))%ptr, &
                                    porosity_loc_p(ghosted_id), &
                                    perm_xx_loc_p(ghosted_id), &
                                    ghosted_id, &  ! for debugging
                                    option)
    endif
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

! right now, hardware to Dirichlet or Hydrostatic
#if 0      
      do idof = 1, option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
#endif
      do idof = 1, option%nflowdof
        xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
      enddo

      ! set this based on data given 
      global_aux_vars_bc(sum_connection)%istate = &
        boundary_condition%flow_condition%iphase
      ! update state and update aux var; this could result in two update to 
      ! the aux var as update state updates if the state changes
      call GeneralAuxVarUpdateState(xxbc,gen_aux_vars_bc(sum_connection), &
                         global_aux_vars_bc(sum_connection), &
                         patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                         porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                         
                         ghosted_id,option)
      call GeneralAuxVarCompute(xxbc,gen_aux_vars_bc(sum_connection), &
                         global_aux_vars_bc(sum_connection), &
                         patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                         porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                         
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecRestoreArrayReadF90(field%perm_xx_loc,perm_xx_loc_p,ierr)
  call VecRestoreArrayReadF90(field%porosity_loc,porosity_loc_p,ierr)  

  patch%aux%General%aux_vars_up_to_date = PETSC_TRUE

end subroutine GeneralUpdateAuxVarsPatch

! ************************************************************************** !
!
! GeneralInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralInitializeTimestep(realization)

  use Realization_class
  
  implicit none
  
  type(realization_type) :: realization

  call GeneralUpdateFixedAccum(realization)

end subroutine GeneralInitializeTimestep

! ************************************************************************** !
!
! GeneralUpdateSolution: Updates data in module after a successful time 
!                             step
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralUpdateSolution(realization)

  use Realization_class
  use Field_module
  use Level_module
  use Patch_module
  use Discretization_module
  
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
  
  ! update ghosted iphase_loc values (must come after 
  ! GeneralUpdateSolutionPatch)
  call DiscretizationLocalToLocal(realization%discretization, &
                                  field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)
  
end subroutine GeneralUpdateSolution


! ************************************************************************** !
!
! GeneralUpdateSolutionPatch: Updates data in module after a successful time 
!                             step
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralUpdateSolutionPatch(realization)

  use Realization_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
    
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  
  PetscInt :: local_id, ghosted_id
  PetscReal, pointer :: iphas_loc_p(:)
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  gen_aux_vars => patch%aux%General%aux_vars  
  global_aux_vars => patch%aux%Global%aux_vars
  
  if (realization%option%compute_mass_balance_new) then
    call GeneralUpdateMassBalancePatch(realization)
  endif
  
  ! update stored state
  call VecGetArrayF90(field%iphas_loc,iphas_loc_p,ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    iphas_loc_p(ghosted_id) = global_aux_vars(ghosted_id)%istate
  enddo
  call VecRestoreArrayF90(field%iphas_loc,iphas_loc_p,ierr)

end subroutine GeneralUpdateSolutionPatch

! ************************************************************************** !
!
! GeneralUpdateFixedAccum: Updates the fixed portion of the 
!                                  accumulation term
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralUpdateFixedAccum(realization)

  use Realization_class
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
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralUpdateFixedAccumPatch(realization)

  use Realization_class
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

  PetscInt :: ghosted_id, local_id, local_start, local_end
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
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
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr)
  call VecGetArrayReadF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayReadF90(field%tortuosity_loc,tor_loc_p,ierr)
  call VecGetArrayReadF90(field%volume,volume_p,ierr)
  call VecGetArrayReadF90(field%perm_xx_loc,perm_xx_loc_p,ierr)  

  call VecGetArrayF90(field%flow_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call GeneralAuxVarCompute(xx_p(local_start:local_end), &
                              gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                              global_aux_vars(ghosted_id), &
                              patch%saturation_function_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              porosity_loc_p(ghosted_id), &
                              perm_xx_loc_p(ghosted_id), &                        
                              option)
    call GeneralAccumulation(gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                             global_aux_vars(ghosted_id), &
                             material_parameter%dencpr(imat), &
                             porosity_loc_p(ghosted_id), &
                             volume_p(local_id), &
                             option,accum_p(local_start:local_end)) 
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr)
  call VecRestoreArrayReadF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayReadF90(field%tortuosity_loc,tor_loc_p,ierr)
  call VecRestoreArrayReadF90(field%volume,volume_p,ierr)
  call VecRestoreArrayReadF90(field%perm_xx_loc,perm_xx_loc_p,ierr)  

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)

end subroutine GeneralUpdateFixedAccumPatch

! ************************************************************************** !
!
! GeneralNumericalJacTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 03/10/11
!
! ************************************************************************** !
subroutine GeneralNumericalJacTest(xx,realization)

  use Realization_class
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
  call VecGetArrayF90(res,vec2_p,ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    idof = icell
!    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call VecCopy(xx,xx_pert,ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr)
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
      call VecRestoreArrayF90(res_pert,vec_p,ierr)
!    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr)

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
! GeneralAuxVarPerturb: Calculates auxiliary variables for perturbed system
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralAuxVarPerturb(gen_aux_var,global_aux_var, &
                                saturation_function,ghosted_id, &
                                option)

  use Option_module
  use Saturation_Function_module

  implicit none

  type(option_type) :: option
  PetscInt :: ghosted_id
  type(general_auxvar_type) :: gen_aux_var(0:)
  type(global_auxvar_type) :: global_aux_var
  type(saturation_function_type) :: saturation_function
     
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert(option%nflowdof)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  PetscInt :: idof

#ifdef DEBUG_GENERAL
  character(len=MAXWORDLENGTH) :: word
  type(global_auxvar_type) :: global_aux_var_debug
  type(general_auxvar_type) :: general_aux_var_debug
  call GlobalAuxVarInit(global_aux_var_debug,option)
  call GeneralAuxVarInit(general_aux_var_debug,option)
#endif

  select case(global_aux_var%istate)
    case(LIQUID_STATE)
       x(GENERAL_LIQUID_PRESSURE_DOF) = gen_aux_var(ZERO_INTEGER)%pres(option%liquid_phase)
       x(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF) = gen_aux_var(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
       x(GENERAL_LIQUID_STATE_TEMPERATURE_DOF) = gen_aux_var(ZERO_INTEGER)%temp
       pert(GENERAL_LIQUID_PRESSURE_DOF) = 1.d0
       pert(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF) = -1.d0*perturbation_tolerance*x(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF)
       pert(GENERAL_LIQUID_STATE_TEMPERATURE_DOF) = -1.d0*perturbation_tolerance*x(GENERAL_LIQUID_STATE_TEMPERATURE_DOF)
    case(GAS_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = gen_aux_var(ZERO_INTEGER)%pres(option%gas_phase)
       x(GENERAL_AIR_PRESSURE_DOF) = gen_aux_var(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(GENERAL_GAS_STATE_TEMPERATURE_DOF) = gen_aux_var(ZERO_INTEGER)%temp
       ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a) and
       ! p(v) >= Psat (i.e. an increase in p(v)) results in two phase.
       pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0
       if (x(GENERAL_GAS_PRESSURE_DOF) - x(GENERAL_AIR_PRESSURE_DOF) > 1.d0) then 
         pert(GENERAL_AIR_PRESSURE_DOF) = 1.d0
       else
         pert(GENERAL_AIR_PRESSURE_DOF) = -1.d0
       endif
       pert(GENERAL_GAS_STATE_TEMPERATURE_DOF) = perturbation_tolerance*x(GENERAL_GAS_STATE_TEMPERATURE_DOF)
    case(TWO_PHASE_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = gen_aux_var(ZERO_INTEGER)%pres(option%gas_phase)
       x(GENERAL_AIR_PRESSURE_DOF) = gen_aux_var(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(GENERAL_GAS_SATURATION_DOF) = gen_aux_var(ZERO_INTEGER)%sat(option%gas_phase)
       pert(GENERAL_GAS_PRESSURE_DOF) = 1.d0
       if (x(GENERAL_GAS_PRESSURE_DOF) - x(GENERAL_AIR_PRESSURE_DOF) > 1.d0) then 
         pert(GENERAL_AIR_PRESSURE_DOF) = 1.d0
       else
         pert(GENERAL_AIR_PRESSURE_DOF) = -1.d0
       endif
       if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then 
         pert(GENERAL_GAS_SATURATION_DOF) = -perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)
       else
         pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)
       endif
  end select
  
  do idof = 1, option%nflowdof
    gen_aux_var(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    call GeneralAuxVarCompute(x_pert,gen_aux_var(idof),global_aux_var, &
                              saturation_function,0.d0,0.d0,option)
#ifdef DEBUG_GENERAL
    call GlobalAuxVarCopy(global_aux_var,global_aux_var_debug,option)
    call GeneralAuxVarCopy(gen_aux_var(idof),general_aux_var_debug,option)
    call GeneralAuxVarUpdateState(x_pert,general_aux_var_debug, &
                                  global_aux_var_debug,saturation_function, &
                                  0.d0,0.d0,ghosted_id,option)
    if (global_aux_var%istate /= global_aux_var_debug%istate) then
      write(option%io_buffer, &
            &'(''Change in state due to perturbation: '',i3,'' -> '',i3)') &
        global_aux_var%istate, global_aux_var_debug%istate
      call printMsg(option)
    endif
#endif

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
  
  PetscReal :: v_over_t
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  v_over_t = vol / option%flow_dt
  
  ! accumulation term units = kmol/s
  Res = 0.d0
  do icomp = 1, option%nflowspec
    do iphase = 1, option%nphase
      ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] * den[kmol phase/m^3 phase] * 
      !                          xmol[kmol comp/kmol phase]
      Res(icomp) = Res(icomp) + gen_aux_var%sat(iphase) * &
                                gen_aux_var%den(iphase) * &
                                gen_aux_var%xmol(icomp,iphase)
    enddo
  enddo
  
  ! some all components into first dof
  do icomp = 2, option%nflowspec
    ! Resk[mol total/m^3 void] = sum(Res[kmol comp/m^3 void])
    Res(ONE_INTEGER) = Res(ONE_INTEGER) + Res(icomp)
  enddo
  
  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * por * v_over_t
  
  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
    Res(energy_id) = Res(energy_id) + gen_aux_var%sat(iphase) * &
                                      gen_aux_var%den(iphase) * &
                                      gen_aux_var%U(iphase)
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * por + &
                    (1.d0 - por) * dencpr * gen_aux_var%temp) * v_over_t 

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
                       por_up,sir_up,dd_up,perm_up,tor_up, &
                       gen_aux_var_dn,global_aux_var_dn, &
                       por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                       area,dist_gravity,upweight,general_parameter, &
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
  PetscReal :: tor_up, tor_dn
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: upweight
  PetscReal :: uH
  type(general_parameter_type) :: general_parameter
  PetscReal :: Res(option%nflowdof)

  PetscReal :: upweight_adj
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: fmw_phase(option%nphase)
  PetscReal :: xmol(option%nflowspec)
  PetscReal :: den
  PetscReal :: density_ave
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: ukvr, mole_flux, q
  PetscReal :: stp_up, stp_dn
  PetscReal :: temp_ave, stp_ave, theta, v_air
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  
#ifdef DEBUG_GENERAL_LOCAL 
  PetscReal :: flux(option%nflowdof)
  call GeneralPrintAuxVars(gen_aux_var_up,global_aux_var_up,-999, &
                           'Upwind',option)  
  call GeneralPrintAuxVars(gen_aux_var_dn,global_aux_var_dn,-999, &
                           'Downwind',option)  
#endif

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  fmw_phase(option%liquid_phase) = FMWH2O
  fmw_phase(option%gas_phase) = FMWAIR

  perm_ave_over_dist = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  Res = 0.d0
  v_darcy = 0.d0
  
#ifdef DEBUG_GENERAL_LOCAL
  flux = 0.d0
#endif
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
      ! MJ/kmol
      H_ave = upweight_adj*gen_aux_var_up%H(iphase)+ &
              (1.D0-upweight_adj)*gen_aux_var_dn%H(iphase)

      !geh: dist_gravity is the distance * gravity in the direction of 
      !     gravity (negative if gravity is down)      
      gravity_term = (upweight_adj*gen_aux_var_up%den(iphase) + &
                     (1.D0-upweight)*gen_aux_var_dn%den(iphase)) &
                     * fmw_phase(iphase) * dist_gravity 

      delta_pressure = gen_aux_var_up%pres(iphase) - &
                       gen_aux_var_dn%pres(iphase) + &
                       gravity_term

      if (delta_pressure >= 0.D0) then
        ukvr = gen_aux_var_up%kvr(iphase)
        xmol(:) = gen_aux_var_up%xmol(:,iphase)
        den = density_ave
        uH = H_ave
      else
        ukvr = gen_aux_var_dn%kvr(iphase)
        xmol(:) = gen_aux_var_dn%xmol(:,iphase)
        den = density_ave
        uH = H_ave
      endif      

      if (ukvr > floweps) then
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(iphase) = perm_ave_over_dist * ukvr * delta_pressure
        ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy(iphase) * area  
        ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
        !                             density_ave[kmol phase/m^3 phase]        
        mole_flux = q*den       
        do icomp = 1, option%nflowspec
          ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
          !                      xmol[kmol comp/kmol phase]
#ifdef DEBUG_GENERAL_LOCAL 
          flux(icomp) = flux(icomp) + mole_flux * xmol(icomp)
#endif
          Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
        enddo
        ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
#ifdef DEBUG_GENERAL_LOCAL 
        flux(energy_id) = flux(energy_id) + mole_flux * uH
#endif
        Res(energy_id) = Res(energy_id) + mole_flux * uH
      endif                   
    endif ! sat > eps
  enddo
  
  do icomp = 2, option%nflowspec
    ! Res[kmol total/sec] = sum(Res[kmol comp/sec])  
    Res(ONE_INTEGER) = Res(ONE_INTEGER) + Res(icomp)
  enddo
  
#ifdef DEBUG_GENERAL_LOCAL 
  print *, 'Darcy flux: ', flux(1:3)
  flux = 0.d0
#endif  

#if 0
  ! add in gas component diffusion in gas and liquid phases
  do iphase = 1, option%nphase
    theta = 1.8d0
    if (gen_aux_var_up%sat(iphase) > eps .and. &
        gen_aux_var_dn%sat(iphase) > eps) then
      upweight_adj = upweight
      if (gen_aux_var_up%sat(iphase) < eps) then 
        upweight_adj=0.d0
      else if (gen_aux_var_dn%sat(iphase) < eps) then 
        upweight_adj=1.d0
      endif         
! not useing harmonic mean
!      stp_up = gen_aux_var_up%sat(iphase)*tor_up*por_up
!      stp_dn = gen_aux_var_dn%sat(iphase)*tor_dn*por_dn
      ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
      density_ave = upweight_adj*gen_aux_var_up%den(iphase)+ &
                    (1.D0-upweight_adj)*gen_aux_var_dn%den(iphase)
!      stp_ave = (stp_up*stp_dn)/(stp_up*dd_dn+stp_dn*dd_up)
      stp_ave = sqrt(gen_aux_var_up%sat(iphase)*gen_aux_var_dn%sat(iphase))* &
                sqrt(tor_up*tor_dn)* &
                sqrt(por_up*por_dn)
      delta_xmol = gen_aux_var_up%xmol(air_comp_id,iphase) - &
                   gen_aux_var_dn%xmol(air_comp_id,iphase)
      ! need to account for multiple phases
      ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
      if (iphase == option%liquid_phase) then
        ! Eq. 1.9a.  The water density is added below
        v_air = stp_ave * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol
      else
        temp_ave = upweight_adj*gen_aux_var_up%temp + &
                   (1.d0-upweight_adj)*gen_aux_var_dn%temp
        pressure_ave = upweight_adj*gen_aux_var_up%pres(iphase)+ &
                      (1.D0-upweight_adj)*gen_aux_var_dn%pres(iphase)
        ! Eq. 1.9b.  The gas density is added below
        v_air = stp_ave * &
                ((temp_ave+273.15)/273.15d0)**theta * &
                option%reference_pressure / pressure_ave * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol      
      endif      
      q =  v_air * area
      mole_flux = q * density_ave
#ifdef DEBUG_GENERAL_LOCAL 
      flux(air_comp_id) = flux(air_comp_id) + mole_flux
#endif
      Res(air_comp_id) = Res(air_comp_id) + mole_flux
    endif
  enddo
#endif

#ifdef DEBUG_GENERAL_LOCAL 
  print *, 'Gas diffusive flux: ', flux(air_comp_id)
  flux = 0.d0
#endif  
    
#if 0
  ! add heat conduction flux
  k_eff_up = gen_aux_var_up%sat(option%liquid_phase) * &
             general_parameter%thermal_conductivity(option%liquid_phase) + &
             gen_aux_var_up%sat(option%gas_phase) * &
             general_parameter%thermal_conductivity(option%gas_phase)
  k_eff_dn = gen_aux_var_dn%sat(option%liquid_phase) * &
             general_parameter%thermal_conductivity(option%liquid_phase) + &
             gen_aux_var_dn%sat(option%gas_phase) * &
             general_parameter%thermal_conductivity(option%gas_phase)
  if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dd_dn+k_eff_dn*dd_up)
  else
    k_eff_ave = 0.d0
  endif
  delta_temp = gen_aux_var_up%temp - gen_aux_var_dn%temp
  heat_flux = k_eff_ave * delta_temp * area
#ifdef DEBUG_GENERAL_LOCAL 
  flux(energy_id) = flux(energy_id) + heat_flux
#endif 
  Res(energy_id) = Res(energy_id) + heat_flux
    
#ifdef DEBUG_GENERAL_LOCAL 
  print *, 'Energy conductive flux: ', flux(energy_id)
#endif

! if 0
#endif
  
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
                                  sir_up,dd_up,perm_up,tor_up, &
                                  gen_aux_var_dn,global_aux_var_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn,tor_dn, &
                                  area,dist_gravity,upweight, &
                                  general_parameter, &
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
  PetscReal :: tor_up, tor_dn
  PetscReal :: area
  PetscReal :: dist_gravity
  PetscReal :: upweight
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  call GeneralFlux(gen_aux_var_up(ZERO_INTEGER),global_aux_var_up, &
                   por_up,sir_up,dd_up,perm_up,tor_up, &
                   gen_aux_var_dn(ZERO_INTEGER),global_aux_var_dn, &
                   por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                   area,dist_gravity,upweight,general_parameter, &
                   option,v_darcy,res)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_aux_var_up(idof),global_aux_var_up, &
                     por_up,sir_up,dd_up,perm_up,tor_up, &
                     gen_aux_var_dn(ZERO_INTEGER),global_aux_var_dn, &
                     por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                     area,dist_gravity,upweight,general_parameter, &
                     option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res(irow))/gen_aux_var_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_aux_var_up(ZERO_INTEGER),global_aux_var_up, &
                     por_up,sir_up,dd_up,perm_up,tor_up, &
                     gen_aux_var_dn(idof),global_aux_var_dn, &
                     por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                     area,dist_gravity,upweight,general_parameter, &
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
                         por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                         area,dist_gravity,general_parameter, &
                         option,v_darcy,Res)
  use Option_module                              
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(general_auxvar_type) :: gen_aux_var_up, gen_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_dn, sir_dn(option%nphase)
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn, perm_dn, tor_dn
  PetscReal :: v_darcy(option%nphase), area
  PetscReal :: dist_gravity
  type(general_parameter_type) :: general_parameter
  PetscReal :: Res(1:option%nflowdof)

  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase

  PetscReal :: fmw_phase(option%nphase)
  PetscReal :: xmol(option%nflowspec)  
  PetscReal :: density_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_ave_over_dist
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: gravity
  PetscReal :: ukvr, mole_flux, q
  PetscReal :: temp_ave, stp_ave, theta, v_air
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  
  PetscInt :: bc_type
  PetscInt :: idof
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  fmw_phase(option%liquid_phase) = FMWH2O
  fmw_phase(option%gas_phase) = FMWAIR

  Res = 0.d0
  v_darcy = 0.d0
  
  do iphase = 1, option%nphase
  
#if 0  
    select case(iphase)
      case(LIQUID_PHASE)
        bc_type = ibndtype(GENERAL_LIQUID_PRESSURE_DOF)
      case(GAS_PHASE)
        bc_type = ibndtype(GENERAL_GAS_PRESSURE_DOF)
    end select
#endif
    bc_type = DIRICHLET_BC

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
          perm_ave_over_dist = perm_dn / dd_dn
        endif
        
          
        ! using residual saturation cannot be correct! - geh
        if (gen_aux_var_up%sat(iphase) > sir_dn(iphase) .or. &
            gen_aux_var_dn%sat(iphase) > sir_dn(iphase)) then
          upweight = 1.d0
          if (gen_aux_var_up%sat(iphase) < eps) then 
            upweight=0.d0
          else if (gen_aux_var_dn%sat(iphase) < eps) then 
            upweight=1.d0
          endif    
          density_ave = upweight*gen_aux_var_up%den(iphase)+ &
                        (1.D0-upweight)*gen_aux_var_dn%den(iphase)
          ! MJ/kmol
!geh          H_ave = upweight*gen_aux_var_up%H(iphase)+ &
!geh                  (1.D0-upweight)*gen_aux_var_dn%H(iphase)

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
                gen_aux_var_up%pres(iphase)-option%reference_pressure < eps) then
              delta_pressure = 0.d0
            endif
          endif
            
          if (delta_pressure >= 0.D0) then
            ukvr = gen_aux_var_up%kvr(iphase)
            xmol(:) = gen_aux_var_up%xmol(:,iphase)
            uH = gen_aux_var_up%H(iphase)
          else
            ukvr = gen_aux_var_dn%kvr(iphase)
            xmol(:) = gen_aux_var_dn%xmol(:,iphase)
            uH = gen_aux_var_dn%H(iphase)
          endif      

          if (ukvr > floweps) then
            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
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
            density_ave = gen_aux_var_up%den(iphase)
          else 
            density_ave = gen_aux_var_dn%den(iphase)
          endif 
        endif
    end select

    if (dabs(v_darcy(iphase)) > 0.d0) then
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                              density_ave[kmol phase/m^3 phase]
      mole_flux = q*density_ave       
      do icomp = 1, option%nflowspec
        ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
        !                      xmol[kmol comp/mol phase]
        Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      enddo
      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH ! H_ave
    endif
  enddo
  
  do icomp = 2, option%nflowspec
    ! Res[kmol total/sec] = sum(Res[kmol comp/sec])
    Res(ONE_INTEGER) = Res(ONE_INTEGER) + Res(icomp)
  enddo

#if 1
  ! add in gas component diffusion in gas and liquid phases
  do iphase = 1, option%nphase
    theta = 1.d0
  
    if (gen_aux_var_up%sat(iphase) > eps .and. &
        gen_aux_var_dn%sat(iphase) > eps) then
      upweight = 1.d0
      if (gen_aux_var_up%sat(iphase) < eps) then 
        upweight=0.d0
      else if (gen_aux_var_dn%sat(iphase) < eps) then 
        upweight=1.d0
      endif       
      ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
      temp_ave = upweight*gen_aux_var_up%temp + &
              (1.d0-upweight)*gen_aux_var_dn%temp
      density_ave = upweight*gen_aux_var_up%den(iphase)+ &
                    (1.D0-upweight)*gen_aux_var_dn%den(iphase)             
  !    stp_ave = tor_dn*por_dn*(gen_aux_var_up%sat(iphase)*gen_aux_var_dn%sat(iphase))/ &
  !              ((gen_aux_var_up%sat(iphase)+gen_aux_var_dn%sat(iphase))*dd_dn)
      ! should saturation be distance weighted?
      stp_ave = tor_dn*por_dn*gen_aux_var_dn%sat(iphase)/dd_dn
      delta_xmol = gen_aux_var_up%xmol(air_comp_id,iphase) - &
                   gen_aux_var_dn%xmol(air_comp_id,iphase)
      ! need to account for multiple phases
      ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
      if (iphase == option%liquid_phase) then
        ! Eq. 1.9a.  The water density is added below
        v_air = stp_ave * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol
      else
        ! Eq. 1.9b.  The gas density is added below
        v_air = stp_ave * &
                (temp_ave/273.15d0)**theta * &
                option%reference_pressure / gen_aux_var_dn%pres(iphase) * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol      
      endif
      q =  v_air * area
      mole_flux = q * density_ave
      Res(air_comp_id) = Res(air_comp_id) + mole_flux
    endif
  enddo
#endif

  ! add heat conduction flux
!  k_eff_up = gen_aux_var_up%sat(option%liquid_phase) * &
!             general_parameter%diffusion_coefficient(option%liquid_phase) + &
!             gen_aux_var_up%sat(option%gas_phase) * &
!             general_parameter%diffusion_coefficient(option%gas_phase)
  k_eff_dn = gen_aux_var_dn%sat(option%liquid_phase) * &
             general_parameter%diffusion_coefficient(option%liquid_phase) + &
             gen_aux_var_dn%sat(option%gas_phase) * &
             general_parameter%diffusion_coefficient(option%gas_phase)
!  k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dd_dn+k_eff_dn*0)
  ! should k_eff be distance weighted?
  k_eff_ave = k_eff_dn/dd_dn
  delta_temp = gen_aux_var_up%temp - gen_aux_var_dn%temp
  heat_flux = k_eff_ave * delta_temp * area
  Res(energy_id) = Res(energy_id) + heat_flux

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
                                    por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                                    area,dist_gravity,general_parameter, &
                                    option,Jdn)

  use Option_module 
  
  implicit none

  PetscInt :: ibndtype(:)
  PetscReal :: aux_vars(:) ! from aux_real_var array
  type(general_auxvar_type) :: gen_aux_var_up, gen_aux_var_dn(0:)
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_dn, sir_dn(option%nphase)
  PetscReal :: por_dn,perm_dn,tor_dn
  PetscReal :: area
  PetscReal :: dist_gravity
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  call GeneralBCFlux(ibndtype,aux_vars, &
                     gen_aux_var_up,global_aux_var_up, &
                     gen_aux_var_dn(ZERO_INTEGER),global_aux_var_dn, &
                     por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                     area,dist_gravity,general_parameter, &
                     option,v_darcy,res)                     
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralBCFlux(ibndtype,aux_vars, &
                       gen_aux_var_up,global_aux_var_up, &
                       gen_aux_var_dn(idof),global_aux_var_dn, &
                       por_dn,sir_dn,dd_dn,perm_dn,tor_dn, &
                       area,dist_gravity,general_parameter, &
                       option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_aux_var_dn(idof)%pert
    enddo !irow
  enddo ! idof

end subroutine GeneralBCFluxDerivative

! ************************************************************************** !
!
! GeneralSrcSink: Computes the source/sink terms for the residual
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                          gen_aux_var,scale,res)

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(option%nphase)
  PetscInt :: flow_src_sink_type
  type(general_auxvar_type) :: gen_aux_var
  PetscReal :: scale
  PetscReal :: res(option%nflowdof)
      
  PetscReal :: fmw_phase(option%nphase)
  PetscReal :: qsrc_mol
  PetscInt :: icomp
  

  fmw_phase(option%liquid_phase) = FMWH2O
  fmw_phase(option%gas_phase) = FMWAIR

  res = 0.d0
  do icomp = 1, option%nflowspec
    select case(flow_src_sink_type)
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(icomp)/fmw_phase(icomp) ! kg/sec -> kmol/sec
      case(SCALED_MASS_RATE_SS)
        qsrc_mol = qsrc(icomp)/fmw_phase(icomp)*scale ! kg/sec -> kmol/sec
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec
        qsrc_mol = qsrc(icomp)*gen_aux_var%den(icomp) ! den = kmol/m^3
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec
        qsrc_mol = qsrc(icomp)*gen_aux_var%den(icomp)*scale ! den = kmol/m^3
    end select
    res(icomp) = qsrc_mol
    res(option%energy_id) = res(option%energy_id) ! * some sort of heat transfer
  enddo

end subroutine GeneralSrcSink
  
! ************************************************************************** !
!
! GeneralResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralResidual(snes,xx,r,realization,ierr)

  use Realization_class
  use Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module

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
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option
  
  ! Communication -----------------------------------------
  ! These 3 must be called before GeneralUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  
  option%variables_swapped = PETSC_FALSE
  call GeneralResidualPatch1(snes,xx,r,realization,ierr)
  
  if (option%variables_swapped) then
    if (option%mycommsize > 1) then
      option%io_buffer = 'Update of primary dep vars needs to be fixed for parallel.'
      call printErrMsg(option)
    endif
    call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx,NFLOWDOF)
  endif

  call GeneralResidualPatch2(snes,xx,r,realization,ierr)
   
  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Gresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Gxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
end subroutine GeneralResidual

! ************************************************************************** !
!
! GeneralResidualPatch1: Computes the residual equation 
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralResidualPatch1(snes,xx,r,realization,ierr)

  use Water_EOS_module

  use Connection_module
  use Realization_class
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

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: tort_loc_p(:)

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
  type(general_parameter_type), pointer :: general_parameter
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:), gen_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: iphase
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: local_start, local_end
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  gen_aux_vars => patch%aux%General%aux_vars
  gen_aux_vars_bc => patch%aux%General%aux_vars_bc
  general_parameter => patch%aux%General%general_parameter
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
                                             ! do update state
  call GeneralUpdateAuxVarsPatch(realization,PETSC_TRUE)
  patch%aux%General%aux_vars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
    call GeneralZeroMassBalDeltaPatch(realization)
  endif

! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayReadF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayReadF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayReadF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayReadF90(field%tortuosity_loc, tort_loc_p, ierr)
  !print *,' Finished scattering non deriv'

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

      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
   
      call GeneralFlux(gen_aux_vars(ZERO_INTEGER,ghosted_id_up), &
                        global_aux_vars(ghosted_id_up), &
                          porosity_loc_p(ghosted_id_up), &
                          material_parameter%sir(1,icap_up), &
                          dd_up,perm_up,tort_loc_p(ghosted_id_up), &
                        gen_aux_vars(ZERO_INTEGER,ghosted_id_dn), &
                        global_aux_vars(ghosted_id_dn), &
                          porosity_loc_p(ghosted_id_dn), &
                          material_parameter%sir(1,icap_dn), &
                          dd_dn,perm_dn,tort_loc_p(ghosted_id_dn), &
                        cur_connection_set%area(iconn),distance_gravity, &
                        upweight,general_parameter,option,v_darcy,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy
      
      if (local_id_up > 0) then
        local_end = local_id_up * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) + Res(:)
      endif
         
      if (local_id_dn > 0) then
        local_end = local_id_dn * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) - Res(:)
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

      icap_dn = patch%sat_func_id(ghosted_id)

      call GeneralBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                gen_aux_vars_bc(sum_connection), &
                                global_aux_vars_bc(sum_connection), &
                                gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                material_parameter%sir(:,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn, &
                                tort_loc_p(ghosted_id), &
                                cur_connection_set%area(iconn), &
                                distance_gravity,general_parameter,option, &
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

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      r_p(local_start:local_end)= r_p(local_start:local_end) - Res(:)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(r,r_p,ierr)
  call VecRestoreArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayReadF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayReadF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayReadF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayReadF90(field%tortuosity_loc, tort_loc_p, ierr)

end subroutine GeneralResidualPatch1

! ************************************************************************** !
!
! GeneralResidualPatch2: Computes the residual equation 
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralResidualPatch2(snes,xx,r,realization,ierr)

  use Water_EOS_module

  use Connection_module
  use Realization_class
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

  PetscReal :: scale
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
  PetscInt :: local_start, local_end
  
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
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr)
  call VecGetArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayReadF90(field%volume, volume_p, ierr)

  ! Accumulation terms ------------------------------------
  if (.not.option%steady_state) then

    r_p = r_p - accum_p
    
#ifdef DEBUG_GENERAL_LOCAL
    open(unit=86,file='accum.txt')
    write(86,*) accum_p(:)
    close(86)
#endif     

    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      call GeneralAccumulation(gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                material_parameter%dencpr(imat), &
                                porosity_loc_p(ghosted_id), &
                                volume_p(local_id), &
                                option,Res) 
      r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    enddo

  endif
#ifdef GENERAL
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      call GeneralSrcSink(option,source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(:), &
                        source_sink%flow_condition%rate%itype, &
                        gen_aux_vars(ZERO_INTEGER,ghosted_id), &
                        scale,Res)

      r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)

    enddo
    source_sink => source_sink%next
  enddo
#endif
  if (patch%aux%General%inactive_cells_exist) then
    do i=1,patch%aux%General%n_zero_rows
      r_p(patch%aux%General%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayReadF90(field%volume, volume_p, ierr)

end subroutine GeneralResidualPatch2

! ************************************************************************** !
!
! GeneralJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_class
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

  call GeneralJacobianPatch1(snes,xx,J,J,flag,realization,ierr)
  call GeneralJacobianPatch2(snes,xx,J,J,flag,realization,ierr)

!  norm = 1.d0
!  call MatZeroRowsLocal(A,1,1,norm,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
!  call MatZeroRowsLocal(A,1,2,norm,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  if (realization%debug%matview_Jacobian) then
#if 1  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Gjacobian.out', &
                              viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Gjacobian.bin', &
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
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralJacobianPatch1(snes,xx,A,B,flag,realization,ierr)
       
  use Water_EOS_module

  use Connection_module
  use Realization_class
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
  PetscReal, pointer :: tort_loc_p(:)
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
  type(general_parameter_type), pointer :: general_parameter
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
  general_parameter => patch%aux%General%general_parameter
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  call VecGetArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayReadF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayReadF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayReadF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayReadF90(field%tortuosity_loc, tort_loc_p, ierr)

  ! Perturb aux vars
  do ghosted_id = 1, grid%ngmax  ! For each local node do...
    if (patch%imat(ghosted_id) <= 0) cycle
    call GeneralAuxVarPerturb(gen_aux_vars(:,ghosted_id), &
                              global_aux_vars(ghosted_id), &
                              patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                              ghosted_id,option)
  enddo
  
#ifdef DEBUG_GENERAL_LOCAL
  call GeneralOutputAuxVars(gen_aux_vars,global_aux_vars,option)
#endif 
  
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
    
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
                              
      call GeneralFluxDerivative(gen_aux_vars(:,ghosted_id_up), &
                                  global_aux_vars(ghosted_id_up), &
                                    porosity_loc_p(ghosted_id_up), &
                                    material_parameter%sir(1,icap_up), &
                                    dd_up,perm_up,tort_loc_p(ghosted_id_up), &
                                  gen_aux_vars(:,ghosted_id_dn), &
                                  global_aux_vars(ghosted_id_dn), &
                                    porosity_loc_p(ghosted_id_dn), &
                                    material_parameter%sir(1,icap_dn), &
                                    dd_dn,perm_dn,tort_loc_p(ghosted_id_dn), &
                                  cur_connection_set%area(iconn),distance_gravity, &
                                  upweight,general_parameter,option,&
                                  Jup,Jdn)
      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)
      endif
      if (local_id_dn > 0) then
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

      icap_dn = patch%sat_func_id(ghosted_id)

      call GeneralBCFluxDerivative(boundary_condition%flow_condition%itype, &
                                  boundary_condition%flow_aux_real_var(:,iconn), &
                                  gen_aux_vars_bc(sum_connection), &
                                  global_aux_vars_bc(sum_connection), &
                                  gen_aux_vars(:,ghosted_id), &
                                  global_aux_vars(ghosted_id), &
                                  porosity_loc_p(ghosted_id), &
                                  material_parameter%sir(:,icap_dn), &
                                  cur_connection_set%dist(0,iconn),perm_dn, &
                                  tort_loc_p(ghosted_id), &
                                  cur_connection_set%area(iconn), &
                                  distance_gravity,general_parameter,option, &
                                  Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
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
  
  call VecRestoreArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayReadF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayReadF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayReadF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayReadF90(field%tortuosity_loc, tort_loc_p, ierr)

end subroutine GeneralJacobianPatch1

! ************************************************************************** !
!
! GeneralJacobianPatch2: Computes the Jacobian
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralJacobianPatch2(snes,xx,A,B,flag,realization,ierr)
       
  use Water_EOS_module

  use Connection_module
  use Realization_class
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

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:)
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

  call VecGetArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayReadF90(field%volume, volume_p, ierr)
  
  if (.not.option%steady_state) then

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    icap = patch%sat_func_id(ghosted_id)
    call GeneralAccumDerivative(gen_aux_vars(:,ghosted_id), &
                              global_aux_vars(ghosted_id), &
                              material_parameter%dencpr(imat), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option, &
                              Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
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
    
    qsrc = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)

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
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)  

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
  
  call VecRestoreArrayReadF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayReadF90(field%volume, volume_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
  if (patch%aux%General%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%General%n_zero_rows, &
                          patch%aux%General%zero_rows_local_ghosted, &
                          qsrc,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif

end subroutine GeneralJacobianPatch2

! ************************************************************************** !
!
! GeneralCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralCreateZeroArray(patch,option)

  use Realization_class
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
  
  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX, &
                     option%mycomm,ierr)
  if (flag > 0) patch%aux%General%inactive_cells_exist = PETSC_TRUE
  if (ncount /= n_zero_rows) then
    if (option%myrank == option%io_rank) then
      print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    endif
    stop
  endif

end subroutine GeneralCreateZeroArray

! ************************************************************************** !
!
! GeneralMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralMaxChange(realization)

  use Realization_class
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
!  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
!  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dcmax,ierr)
!  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
!  call VecStrideNorm(field%flow_dxx,TWO_INTEGER,NORM_INFINITY,option%dtmpmax,ierr)

end subroutine GeneralMaxChange

! ************************************************************************** !
!
! GeneralCheckUpdatePre: Checks update prior to update
! author: Glenn Hammond
! date: 03/06/13
!
! ************************************************************************** !
subroutine GeneralCheckUpdatePre(line_search,X,dX,changed,realization,ierr)

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  ! ignore changed flag for now.
  PetscBool :: changed
  type(realization_type) :: realization
  
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: P_R, X0, X1, delX
  PetscReal :: scale
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  gen_aux_vars => realization%patch%aux%General%aux_vars
  global_aux_vars => realization%patch%aux%Global%aux_vars

  patch => realization%patch

  call VecGetArrayF90(dX,dX_p,ierr)
  call VecGetArrayF90(X,X_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    print *, '-----------------------------------------'
    print *, 'X'
    print *, X_p
    print *, 'dX_p'
    print *, dX_p
    print *, '-----------------------------------------'  
  enddo

  call VecRestoreArrayF90(dX,dX_p,ierr)
  call VecRestoreArrayF90(X,X_p,ierr)

  if (dabs(option%pressure_dampening_factor) > 0.d0) then

    scale = option%pressure_dampening_factor

    call VecGetArrayF90(dX,dX_p,ierr)
    call VecGetArrayF90(X,X_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
    do local_id = 1, grid%nlmax
      delX = dX_p(local_id)
      X0 = X_p(local_id)
      X1 = X0 - delX
      if (X0 < P_R .and. X1 > P_R) then
        write(option%io_buffer,'("U -> S:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),X0,X1 
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_aux_vars(ghosted_id), &
                                  global_aux_vars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call printMsgAnyRank(option)
#endif
      else if (X1 < P_R .and. X0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),X0,X1
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_aux_vars(ghosted_id), &
                                  global_aux_vars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call printMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (X0 < P_R .and. X1 > P_R) then
        dX_p(local_id) = scale*delX
      endif
    enddo
    call VecRestoreArrayF90(dX,dX_p,ierr)
    call VecRestoreArrayF90(X,X_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
  endif

end subroutine GeneralCheckUpdatePre

! ************************************************************************** !
!
! GeneralCheckUpdatePost: Checks update after to update
! author: Glenn Hammond
! date: 03/06/13
!
! ************************************************************************** !
subroutine GeneralCheckUpdatePost(line_search,X0,dX,X1,dX_changed, &
                                   X1_changed,realization,ierr)

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  type(realization_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_aux_vars(:,:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: Res(1)
  PetscReal :: inf_norm
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  gen_aux_vars => realization%patch%aux%General%aux_vars
  global_aux_vars => realization%patch%aux%Global%aux_vars
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
#if 0  
  if (option%check_stomp_norm) then
    call VecGetArrayF90(dP,dP_p,ierr)
    call VecGetArrayF90(P1,P1_p,ierr)
    call VecGetArrayF90(field%volume,volume_p,ierr)
    call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
    
    inf_norm = 0.d0
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
    
      call RichardsAccumulation(rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                volume_p(local_id), &
                                option,Res)
      inf_norm = max(inf_norm,min(dabs(dP_p(local_id)/P1_p(local_id)), &
                                  dabs(r_p(local_id)/Res(1))))
    enddo
    call MPI_Allreduce(inf_norm,option%stomp_norm,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    call VecGetArrayF90(dP,dP_p,ierr)
    call VecGetArrayF90(P1,P1_p,ierr)
    call VecGetArrayF90(field%volume,volume_p,ierr)
    call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
  endif
#endif  
  
end subroutine GeneralCheckUpdatePost

! ************************************************************************** !
!
! GeneralGetTecplotHeader: Returns General Lite contribution to 
!                               Tecplot file header
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
function GeneralGetTecplotHeader(realization,icolumn)
  
  use Realization_class
  use Option_module
  use Field_module
    
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: GeneralGetTecplotHeader
  type(realization_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i

  option => realization%option
  field => realization%field
  
  string = ''
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-T [C]"'')') icolumn
  else
    write(string2,'('',"T [C]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-State"'')') icolumn
  else
    write(string2,'('',"State"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sat(l)"'')') icolumn
  else
    write(string2,'('',"Sat(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sat(g)"'')') icolumn
  else
    write(string2,'('',"Sat(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Rho(l)"'')') icolumn
  else
    write(string2,'('',"Rho(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Rho(g)"'')') icolumn
  else
    write(string2,'('',"Rho(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(l)"'')') icolumn
  else
    write(string2,'('',"U(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(g)"'')') icolumn
  else
    write(string2,'('',"U(g)"'')')
  endif
  string = trim(string) // trim(string2)
  
  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i2,'')"'')') icolumn, i
    else
      write(string2,'('',"Xl('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xg('',i2,'')"'')') icolumn, i
    else
      write(string2,'('',"Xg('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo
 
  GeneralGetTecplotHeader = string

end function GeneralGetTecplotHeader

! ************************************************************************** !
!
! GeneralSetPlotVariables: Adds variables to be printed to list
! author: Glenn Hammond
! date: 02/15/13
!
! ************************************************************************** !
subroutine GeneralSetPlotVariables(realization)
  
  use Realization_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(realization_type) :: realization
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  type(output_variable_type), pointer :: output_variable
  
  list => realization%output_option%output_variable_list

  if (associated(list%first)) then
    return
  endif
  
  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               TEMPERATURE)

  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)

  name = 'Gas Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               GAS_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)
  
  name = 'Gas Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               GAS_SATURATION)
  
  name = 'Liquid Density'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)
  
  name = 'Gas Density'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_DENSITY)
  
  name = 'X_g^l'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION, &
                               realization%option%air_id)
  
  name = 'X_l^l'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION, &
                               realization%option%water_id)
  
  name = 'X_g^g'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOLE_FRACTION, &
                               realization%option%air_id)
  
  name = 'X_l^g'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOLE_FRACTION, &
                               realization%option%water_id)
  
  name = 'Thermodynamic State'
  units = ''
  output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE,units,STATE)
  output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList( &
         realization%output_option%output_variable_list,output_variable)   
  
end subroutine GeneralSetPlotVariables

! ************************************************************************** !
!
! GeneralDestroy: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralDestroy(realization)

  use Realization_class
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
! date: 03/09/11
!
! ************************************************************************** !
subroutine GeneralDestroyPatch(realization)

  use Realization_class

  implicit none

  type(realization_type) :: realization
  
  ! place anything that needs to be freed here.

end subroutine GeneralDestroyPatch

end module General_module
