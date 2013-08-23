module THMC_module

  use THMC_Aux_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"
  
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
  PetscReal, parameter :: perturbation_tolerance = 1.d-8

  public THMCResidual,THMCJacobian, &
         THMCUpdateFixedAccumulation,THMCTimeCut,&
         THMCSetup, THMCNumericalJacobianTest, &
         THMCMaxChange, THMCUpdateSolution, &
         THMCGetTecplotHeader, THMCInitializeTimestep, &
         THMCComputeMassBalance, THMCResidualToMass, &
         THMCUpdateAuxVars, THMCDestroy, &
         THMCComputeDisplacementGradient, &
         THMCComputeStressFromDispGrad, &
         THMCComputeDisplacementGradientPert, &
         THMCComputeMinv, THMCFluxDerivativeAnalytical

  PetscInt, parameter :: jh2o = 1

contains

! ************************************************************************** !
!
! THMCTimeCut: Resets arrays for time step cut
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCTimeCut(realization)
 
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
 
  call VecCopy(field%flow_yy,field%flow_xx,ierr)
  call THMCInitializeTimestep(realization)
 
end subroutine THMCTimeCut


! ************************************************************************** !
!
! THMCSetup: Set up the THMC mode
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCSetup(realization)

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
      call THMCSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THMCSetup
  
! ************************************************************************** !
!
! THMCSetupPatch: Creates arrays for auxiliary variables
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCSetupPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(thmc_auxvar_type), pointer :: thmc_aux_vars(:), thmc_aux_vars_bc(:)
  type(fluid_property_type), pointer :: cur_fluid_property

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, iphase
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
    
  patch%aux%THMC => THMCAuxCreate(option)

  allocate(patch%aux%THMC%thmc_parameter%sir(option%nphase, &
                                  size(realization%saturation_function_array)))
  
  allocate(patch%aux%THMC%thmc_parameter%dencpr(size(realization%material_property_array)))
  allocate(patch%aux%THMC%thmc_parameter%rock_den(size(realization%material_property_array)))
  allocate(patch%aux%THMC%thmc_parameter%ckwet(size(realization%material_property_array)))
  allocate(patch%aux%THMC%thmc_parameter%ckdry(size(realization%material_property_array)))
  allocate(patch%aux%THMC%thmc_parameter%alpha(size(realization%material_property_array)))
  allocate(patch%aux%THMC%thmc_parameter%youngs_modulus(size(realization%material_property_array)))
  allocate(patch%aux%THMC%thmc_parameter%poissons_ratio(size(realization%material_property_array)))

#ifdef ICE
  allocate(patch%aux%THMC%thmc_parameter%ckfrozen(size(realization%material_property_array)))
  allocate(patch%aux%THMC%thmc_parameter%alpha_fr(size(realization%material_property_array)))
#endif

  do i = 1, size(realization%material_property_array)
    patch%aux%THMC%thmc_parameter%dencpr(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%rock_density*option%scale* &
        realization%material_property_array(i)%ptr%specific_heat
    patch%aux%THMC%thmc_parameter%rock_den(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%rock_density 
    patch%aux%THMC%thmc_parameter%ckwet(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%thermal_conductivity_wet*option%scale  
    patch%aux%THMC%thmc_parameter%ckdry(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%thermal_conductivity_dry*option%scale
    patch%aux%THMC%thmc_parameter%alpha(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%alpha
    patch%aux%THMC%thmc_parameter%youngs_modulus(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%youngs_modulus 
    patch%aux%THMC%thmc_parameter%poissons_ratio(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%poissons_ratio
#ifdef ICE
    patch%aux%THMC%thmc_parameter%ckfrozen(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%thermal_conductivity_frozen*option%scale
    patch%aux%THMC%thmc_parameter%alpha_fr(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%alpha_fr
#endif

  enddo 

  do i = 1, size(realization%saturation_function_array)
    patch%aux%THMC%thmc_parameter%sir(:,realization%saturation_function_array(i)%ptr%id) = &
      realization%saturation_function_array(i)%ptr%Sr(:)
  enddo

  ! allocate aux_var data structures for all grid cells
  allocate(thmc_aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call THMCAuxVarInit(thmc_aux_vars(ghosted_id),option)
    ! currently, hardwire to first fluid
    thmc_aux_vars(ghosted_id)%diff(1:option%nflowspec) = &
      realization%fluid_properties%diffusion_coefficient
    ! Calculate M inverse at the beginning of simulation and store it
    call THMCComputeMinv(grid, ghosted_id, thmc_aux_vars(ghosted_id)%Minv, option)   
  enddo

  patch%aux%THMC%aux_vars => thmc_aux_vars
  patch%aux%THMC%num_aux = grid%ngmax
  
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

  if (sum_connection > 0) then 
    allocate(thmc_aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call THMCAuxVarInit(thmc_aux_vars_bc(iconn),option)
      ! currently, hardwire to first fluid
      thmc_aux_vars_bc(iconn)%diff(1:option%nflowspec) = &
        realization%fluid_properties%diffusion_coefficient
    enddo
    patch%aux%THMC%aux_vars_bc => thmc_aux_vars_bc
  endif
  patch%aux%THMC%num_aux_bc = sum_connection

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call THMCCreateZeroArray(patch,option)
  
  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    patch%aux%THMC%thmc_parameter%diffusion_coefficient(iphase) = &
      cur_fluid_property%diffusion_coefficient
    patch%aux%THMC%thmc_parameter%diffusion_activation_energy(iphase) = &
      cur_fluid_property%diffusion_activation_energy
    cur_fluid_property => cur_fluid_property%next
  enddo

end subroutine THMCSetupPatch


! ************************************************************************** !
!
! THMComputeMassBalance: Compute the mass balance for a realization
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCComputeMassBalance(realization, mass_balance)

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
      call THMCComputeMassBalancePatch(realization, mass_balance)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THMCComputeMassBalance    

! ************************************************************************** !
!
! THMComputeMassBalancePatch: Compute mass balance for a patch
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCComputeMassBalancePatch(realization,mass_balance)
 
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

  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    ! mass = volume*saturation*density
    mass_balance = mass_balance + &
      global_aux_vars(ghosted_id)%den_kg* &
      global_aux_vars(ghosted_id)%sat* &
      porosity_loc_p(ghosted_id)*volume_p(local_id)
  enddo

  call VecRestoreArrayF90(field%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  
end subroutine THMCComputeMassBalancePatch

! ************************************************************************** !
!
! THMCZeroMassBalDeltaPatch: Zeros mass balance delta array
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCZeroMassBalDeltaPatch(realization)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%THMC%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%THMC%num_aux_bc > 0) then
    do iconn = 1, patch%aux%THMC%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
 
end subroutine THMCZeroMassBalDeltaPatch

! ************************************************************************** !
!
! THMCUpdateMassBalancePatch: Updates mass balance
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCUpdateMassBalancePatch(realization)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%THMC%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance = &
      patch%aux%Global%aux_vars(iconn)%mass_balance + &
      patch%aux%Global%aux_vars(iconn)%mass_balance_delta*FMWH2O* &
      option%flow_dt
  enddo
#endif

  if (patch%aux%THMC%num_aux_bc > 0) then
    do iconn = 1, patch%aux%THMC%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance = &
        global_aux_vars_bc(iconn)%mass_balance + &
        global_aux_vars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif


end subroutine THMCUpdateMassBalancePatch

! ************************************************************************** !
!
! THMCUpdateAuxVars: Updates the auxiliary variables in a realization
!                    associated with the THMC problem
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCUpdateAuxVars(realization)

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
      call THMCUpdateAuxVarsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THMCUpdateAuxVars

! ************************************************************************** !
!
! THMCUpdateAuxVarsPatch: Updates the auxiliary variables in a patch
!                         associated with the THMC problem
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCUpdateAuxVarsPatch(realization)

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
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(thmc_auxvar_type), pointer :: thmc_aux_vars(:)
  type(thmc_auxvar_type), pointer :: thmc_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(thmc_parameter_type), pointer :: thmc_parameter

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
 !  PetscReal, allocatable :: gradient(:,:,:)
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  
!  allocate(gradient(grid%ngmax,3,3))
!  gradient = 0.d0
  
  thmc_parameter => patch%aux%THMC%thmc_parameter
  thmc_aux_vars => patch%aux%THMC%aux_vars
  thmc_aux_vars_bc => patch%aux%THMC%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
   
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecGetArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
    

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))
    
   
#ifdef ICE
    call THMCAuxVarComputeIce(xx_loc_p(istart:iend), &
        thmc_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
        iphase, &
        realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
        porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
        option)
#else
    call THMCAuxVarCompute(xx_loc_p(istart:iend), &
        thmc_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
        iphase, &
        realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
        porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
        option)
        
#endif

    call THMCComputeDisplacementGradient(grid, global_aux_vars, ghosted_id, &
                       thmc_aux_vars(ghosted_id)%gradient, & 
                       thmc_aux_vars(ghosted_id)%Minv,option) 

    call THMCComputeStressFromDispGrad(thmc_aux_vars(ghosted_id)%gradient, &
              thmc_parameter%youngs_modulus(int(ithrm_loc_p(ghosted_id))), &
              thmc_parameter%poissons_ratio(int(ithrm_loc_p(ghosted_id))), &
              thmc_aux_vars(ghosted_id)%stress)

    iphase_loc_p(ghosted_id) = iphase
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

      do idof=1,option%nflowdof-option%nmechdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
      
      select case(boundary_condition%flow_condition%itype(THMC_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
          iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

!      do idof = option%nflowdof-option%nmechdof+1,option%nflowdof
!        select case(boundary_condition%flow_condition%itype(idof))
!          case(DIRICHLET_BC)
!            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
!        end select
!      enddo

#ifdef ICE
      call THMCAuxVarComputeIce(xxbc,thmc_aux_vars_bc(sum_connection), &
                      global_aux_vars_bc(sum_connection), &
                      iphasebc, &
                      realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                      porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
                      option)
#else
      call THMCAuxVarCompute(xxbc,thmc_aux_vars_bc(sum_connection), &
                      global_aux_vars_bc(sum_connection), &
                      iphasebc, &
                      realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                      porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
                      option)
#endif

    enddo
    boundary_condition => boundary_condition%next
  enddo


  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  
  patch%aux%THMC%aux_vars_up_to_date = PETSC_TRUE
    
!  deallocate(gradient)
end subroutine THMCUpdateAuxVarsPatch

! ************************************************************************** !
!
! THMCInitializeTimestep: Update data in module prior to time step
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCInitializeTimestep(realization)

  use Realization_class
  
  implicit none
  
  type(realization_type) :: realization

  call THMCUpdateFixedAccumulation(realization)

end subroutine THMCInitializeTimestep

! ************************************************************************** !
!
! THMCUpdateSolution: Updates data in module after a successful time step
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCUpdateSolution(realization)

  use Realization_class
  use Field_module
  use Level_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization

  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  field => realization%field
    
  call VecCopy(field%flow_xx,field%flow_yy,ierr)   

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call THMCUpdateSolutionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THMCUpdateSolution

! ************************************************************************** !
!
! THMCUpdateSolutionPatch: Updates data in module after a successful time 
!                          step
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCUpdateSolutionPatch(realization)

  use Realization_class
    
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call THMCUpdateMassBalancePatch(realization)
  endif

end subroutine THMCUpdateSolutionPatch

! ************************************************************************** !
!
! THMCUpdateFixedAccumulation: Updates the fixed portion of the 
!                              accumulation term for a realizations
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCUpdateFixedAccumulation(realization)

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
      call THMCUpdateFixedAccumPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THMCUpdateFixedAccumulation

! ************************************************************************** !
!
! THMCUpdateFixedAccumPatch: Updates the fixed portion of the 
!                            accumulation term for patch
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCUpdateFixedAccumPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(thmc_auxvar_type), pointer :: thmc_aux_vars(:)
  type(thmc_parameter_type), pointer :: thmc_parameter

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:), accum_p(:), perm_xx_loc_p(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  thmc_parameter => patch%aux%THMC%thmc_parameter
  thmc_aux_vars => patch%aux%THMC%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
    
  call VecGetArrayF90(field%flow_xx,xx_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%tortuosity_loc,tor_loc_p,ierr)
  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecGetArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    iend = local_id*option%nflowdof
    istart = iend - option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))

#ifdef ICE
    call THMCAuxVarComputeIce(xx_p(istart:iend), &
                       thmc_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
                       option)
#else
    call THMCAuxVarCompute(xx_p(istart:iend), &
                       thmc_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
                       option)
#endif

    iphase_loc_p(ghosted_id) = iphase
    call THMCAccumulation(thmc_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              thmc_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              thmc_parameter%rock_den(int(ithrm_loc_p(ghosted_id))), &
                              option,accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayF90(field%tortuosity_loc,tor_loc_p,ierr)
  call VecRestoreArrayF90(field%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecRestoreArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)

#if 0
   call THMCNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine THMCUpdateFixedAccumPatch

! ************************************************************************** !
!
! THMCNumericalJacobianTest: Computes the a test numerical jacobian
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCNumericalJacobianTest(xx,realization)

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
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof, &
                   grid%nlmax*option%nflowdof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call THMCResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr)
  do icell = 1,grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    endif
    do idof = (icell-1)*option%nflowdof+1, &
              icell*option%nflowdof 
      call VecCopy(xx,xx_pert,ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr)
      call THMCResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
      call VecGetArrayF90(res_pert,vec_p,ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
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
  call PetscViewerASCIIOpen(option%mycomm,'numerical_jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call MatDestroy(A,ierr)
  
  call VecDestroy(xx_pert,ierr)
  call VecDestroy(res,ierr)
  call VecDestroy(res_pert,ierr)
  
end subroutine THMCNumericalJacobianTest

! ************************************************************************** !
!
! THMCAccumDerivative: Computes derivatives of the accumulation 
!                      term for the Jacobian
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCAccumDerivative(thmc_aux_var,global_aux_var,por,vol, &
                                     rock_dencpr,rock_den,option,sat_func,J)

  use Option_module
  use Saturation_Function_module
  use Water_EOS_module
  
  implicit none

  type(thmc_auxvar_type) :: thmc_aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: vol,por,rock_dencpr,rock_den
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof, option%nflowdof)
     
  PetscInt :: ispec 
  PetscReal :: porXvol, mol(option%nflowspec), eng

  PetscInt :: iphase, ideriv
  type(thmc_auxvar_type) :: thmc_aux_var_pert
  type(global_auxvar_type) :: global_aux_var_pert
  PetscReal :: x(3), x_pert(3), pert, res(3), res_pert(3), J_pert(3,3)
  
#ifdef ICE
  PetscReal :: sat_g, p_g, den_g, p_sat, mol_g, u_g, C_g
  PetscReal :: dpsat_dt, ddeng_dt, dmolg_dt, dsatg_dp, dsatg_dt, dug_dt
  PetscReal :: sat_i, den_i, u_i
  PetscReal :: dsati_dp, dsati_dt
  PetscReal :: ddeni_dp, ddeni_dt
  PetscReal :: dui_dt
  PetscReal, parameter :: C_a = 1.86d-3 ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! in MJ/kg/K
  PetscErrorCode :: ierr
#endif

  porXvol = por*vol
  
  ! X = {p, T, x_2}; R = {R_x1, R_x2, R_T} = {R_p, R_x2, R_T}
  J = 0.d0

  J(1,1) = (global_aux_var%sat(1)*thmc_aux_var%dden_dp + &
           thmc_aux_var%dsat_dp*global_aux_var%den(1))*porXvol !*thmc_aux_var%xmol(1)
  J(1,2) = global_aux_var%sat(1)*thmc_aux_var%dden_dt*porXvol !*thmc_aux_var%xmol(1)
  J(1,3) = 0.d0 !-global_aux_var%sat(1)*global_aux_var%den(1)*porXvol
  J(2,1) = (global_aux_var%sat(1)*thmc_aux_var%dden_dp + &
           thmc_aux_var%dsat_dp*global_aux_var%den(1))*porXvol*thmc_aux_var%xmol(2)
  J(2,2) = global_aux_var%sat(1)*thmc_aux_var%dden_dt*porXvol*thmc_aux_var%xmol(2)
  J(2,3) = global_aux_var%sat(1)*global_aux_var%den(1)*porXvol
  J(3,1) = (thmc_aux_var%dsat_dp*global_aux_var%den(1)*thmc_aux_var%u + &
            global_aux_var%sat(1)*thmc_aux_var%dden_dp*thmc_aux_var%u + &
            global_aux_var%sat(1)*global_aux_var%den(1)*thmc_aux_var%du_dp)*porXvol
  J(3,2) = global_aux_var%sat(1)* &
           (thmc_aux_var%dden_dt*thmc_aux_var%u + &  ! pull %sat outside
            global_aux_var%den(1)*thmc_aux_var%du_dt)*porXvol +  &
           (1.d0 - por)*vol*rock_dencpr 
  J(3,3) = 0.d0

#ifdef ICE 
  ! SK, 11/17/11
  sat_g = thmc_aux_var%sat_gas
  sat_i = thmc_aux_var%sat_ice
  u_i = thmc_aux_var%u_ice
  den_i = thmc_aux_var%den_ice
  p_g = option%reference_pressure ! set to reference pressure
  den_g = p_g/(IDEAL_GAS_CONST*(global_aux_var%temp(1) + 273.15d0))*1.d-3 
  call PSAT(global_aux_var%temp(1), p_sat, dpsat_dt, ierr)
  mol_g = p_sat/p_g
  C_g = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR !in MJ/kmol/K, expression might be different
  u_g = C_g*(global_aux_var%temp(1) + 273.15d0)
  ddeng_dt = - p_g/(IDEAL_GAS_CONST*(global_aux_var%temp(1) + 273.15d0)**2)*1.d-3
  dmolg_dt = dpsat_dt/p_g
  dsatg_dp = thmc_aux_var%dsat_gas_dp
  dsatg_dt = thmc_aux_var%dsat_gas_dt
  dug_dt = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)* &
                 (global_aux_var%temp(1) + 273.15d0)
  dsati_dt = thmc_aux_var%dsat_ice_dt
  dsati_dp = thmc_aux_var%dsat_ice_dp
  ddeni_dt = thmc_aux_var%dden_ice_dt
  ddeni_dp = thmc_aux_var%dden_ice_dp
  dui_dt = thmc_aux_var%du_ice_dt
 
  J(1,1) = J(1,1) + (dsatg_dp*den_g*mol_g + dsati_dp*den_i + &
                    sat_i*ddeni_dp)*porXvol
  J(1,2) = J(1,2) + (thmc_aux_var%dsat_dt*global_aux_var%den(1) + &
                    dsatg_dt*den_g*mol_g + sat_g*ddeng_dt*mol_g + &
                    sat_g*den_g*dmolg_dt + dsati_dt*den_i + sat_i*ddeni_dt)* &
                    porXvol
  J(2,2) = J(2,2) + thmc_aux_var%dsat_dt*global_aux_var%den(1)* &
                    thmc_aux_var%xmol(2)*porXvol
  J(3,1) = J(3,1) + (dsatg_dp*den_g*u_g + dsati_dp*den_i*u_i + &
                    sat_i*ddeni_dp*u_i)*porXvol
  J(3,2) = J(3,2) + (thmc_aux_var%dsat_dt*global_aux_var%den(1)*thmc_aux_var%u + &
                    dsatg_dt*den_g*u_g + sat_g*ddeng_dt*u_g + &
                    sat_g*den_g*dug_dt + dsati_dt*den_i*u_i + &
                    sat_i*ddeni_dt*u_i + sat_i*den_i*dui_dt)*porXvol
#endif

  J = J/option%flow_dt 

  if (option%numerical_derivatives_flow) then
    allocate(thmc_aux_var_pert%xmol(option%nflowspec),thmc_aux_var_pert%diff(option%nflowspec))
    call GlobalAuxVarInit(global_aux_var_pert,option)  
    call THMCAuxVarCopy(thmc_aux_var,thmc_aux_var_pert,option)
    call GlobalAuxVarCopy(global_aux_var,global_aux_var_pert,option)

    x(1) = global_aux_var%pres(1)
    x(2) = global_aux_var%temp(1)
    x(3) = thmc_aux_var%xmol(2)
    

    call THMCAccumulation(thmc_aux_var,global_aux_var, &
                          por,vol,rock_dencpr,rock_den,option,res)
    do ideriv = 1,3
      pert = x(ideriv)*perturbation_tolerance
      x_pert = x
      
#ifdef ICE
      
      if (ideriv == 1) then
        if (x_pert(ideriv) < option%reference_pressure) then
          pert = - pert
        endif
          x_pert(ideriv) = x_pert(ideriv) + pert

      endif
      
      if (ideriv == 2) then
        if(x_pert(ideriv) < 0.d0) then
          pert = - 1.d-8
        else
          pert =  1.d-8
        endif
          x_pert(ideriv) = x_pert(ideriv) + pert
      endif
      
      if (ideriv == 3) then
        x_pert(ideriv) = x_pert(ideriv) + pert
      endif
      
#else
      x_pert(ideriv) = x_pert(ideriv) + pert
      
#endif


#ifdef ICE
      call THMCAuxVarComputeIce(x_pert,thmc_aux_var_pert,global_aux_var_pert,iphase,sat_func, &
                                 0.d0,0.d0,option)
#else
      call THMCAuxVarCompute(x_pert,thmc_aux_var_pert,global_aux_var_pert,iphase,sat_func, &
                                 0.d0,0.d0,option)
#endif

#if 0      
      select case(ideriv)
        case(1)
          print *, 'dsat_dp:', thmc_aux_var%dsat_dp, (global_aux_var_pert%sat(1)-global_aux_var%sat(1))/pert
          print *, 'dden_dp:', thmc_aux_var%dden_dp, (global_aux_var_pert%den(1)-global_aux_var%den(1))/pert
          print *, 'dkvr_dp:', thmc_aux_var%dkvr_dp, (thmc_aux_var_pert%kvr-thmc_aux_var%kvr)/pert
          print *, 'dh_dp:', thmc_aux_var%dh_dp, (thmc_aux_var_pert%h-thmc_aux_var%h)/pert
          print *, 'du_dp:', thmc_aux_var%du_dp, (thmc_aux_var_pert%u-thmc_aux_var%u)/pert
#ifdef ICE
          print *, 'dsati_dp:', thmc_aux_var%dsat_ice_dp, (thmc_aux_var_pert%sat_ice - thmc_aux_var%sat_ice)/pert
          print *, 'dsatg_dp:', thmc_aux_var%dsat_gas_dp, (thmc_aux_var_pert%sat_gas - thmc_aux_var%sat_gas)/pert
          print *, 'ddeni_dp:', thmc_aux_var%dden_ice_dp, (thmc_aux_var_pert%den_ice - thmc_aux_var%den_ice)/pert
#endif
          
        case(2)
          print *, 'dden_dt:', thmc_aux_var%dden_dt, (global_aux_var_pert%den(1)-global_aux_var%den(1))/pert
          print *, 'dkvr_dt:', thmc_aux_var%dkvr_dt, (thmc_aux_var_pert%kvr-thmc_aux_var%kvr)/pert
          print *, 'dh_dt:', thmc_aux_var%dh_dt, (thmc_aux_var_pert%h-thmc_aux_var%h)/pert
          print *, 'du_dt:', thmc_aux_var%du_dt, (thmc_aux_var_pert%u-thmc_aux_var%u)/pert
#ifdef ICE
          print *, 'ddeni_dt:', thmc_aux_var%dden_ice_dt, (thmc_aux_var_pert%den_ice - thmc_aux_var%den_ice)/pert
          print *, 'dsati_dt:', thmc_aux_var%dsat_ice_dt, (thmc_aux_var_pert%sat_ice - thmc_aux_var%sat_ice)/pert
          print *, 'dsatg_dt:', thmc_aux_var%dsat_gas_dt, (thmc_aux_var_pert%sat_gas - thmc_aux_var%sat_gas)/pert
          print *, 'dsat_dt:', thmc_aux_var%dsat_dt, (global_aux_var_pert%sat(1) - global_aux_var%sat(1))/pert
          print *, 'dui_dt:', thmc_aux_var%du_ice_dt, (thmc_aux_var_pert%u_ice - thmc_aux_var%u_ice)/pert
#endif
      end select     
#endif     
      call THMCAccumulation(thmc_aux_var_pert,global_aux_var_pert, &
                          por,vol,rock_dencpr,rock_den,option,res_pert)
      J_pert(:,ideriv) = (res_pert(:)-res(:))/pert
    enddo

    deallocate(thmc_aux_var_pert%xmol,thmc_aux_var_pert%diff)
    J = J_pert
    call GlobalAuxVarStrip(global_aux_var_pert)  
  endif
   
  
end subroutine THMCAccumDerivative

! ************************************************************************** !
!
! THMCAccumulation: Computes the non-fixed portion of the accumulation
!                   term for the residual
! author:
! date: 3/2/12
!
! ************************************************************************** !  
subroutine THMCAccumulation(aux_var,global_aux_var,por,vol,rock_dencpr, &
                            rock_den,option,Res)

  use Option_module
  use Water_EOS_module
  
  implicit none

  type(thmc_auxvar_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: vol,por,rock_dencpr,rock_den
     
  PetscInt :: ispec 
  PetscReal :: porXvol, mol(option%nflowspec), eng

#ifdef ICE
  PetscReal :: sat_g, p_g, den_g, p_sat, mol_g, u_g, C_g
  PetscReal :: sat_i, den_i, u_i
  PetscReal, parameter :: C_a = 1.86d-3 ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! in MJ/kg/K
  PetscErrorCode :: ierr
#endif

  Res = 0.d0
  
! TechNotes, thmc Mode: First term of Equation 8
  porXvol = por*vol
  mol(1) = global_aux_var%sat(1)*global_aux_var%den(1)*porXvol
  mol(2) = global_aux_var%sat(1)*global_aux_var%den(1)*aux_var%xmol(2)*porXvol

! TechNotes, thmc Mode: First term of Equation 9
  eng = global_aux_var%sat(1) * &
        global_aux_var%den(1) * &
        aux_var%u * porXvol + &
        (1.d0 - por) * vol * rock_dencpr * global_aux_var%temp(1)

#ifdef ICE 
  ! SK, 11/17/11
  sat_g = aux_var%sat_gas
  sat_i = aux_var%sat_ice
  u_i = aux_var%u_ice
  den_i = aux_var%den_ice
  p_g = option%reference_pressure
  den_g = p_g/(IDEAL_GAS_CONST*(global_aux_var%temp(1) + 273.15d0))*1.d-3 !in kmol/m3
  call PSAT(global_aux_var%temp(1), p_sat, ierr)
  mol_g = p_sat/p_g
  C_g = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR ! in MJ/kmol/K
  u_g = C_g*(global_aux_var%temp(1) + 273.15d0)       ! in MJ/kmol
  mol(1) = mol(1) + (sat_g*den_g*mol_g + sat_i*den_i)*porXvol
  eng = eng + (sat_g*den_g*u_g + sat_i*den_i*u_i)*porXvol
#endif

  Res(1:option%nflowspec) = mol(:)/option%flow_dt
  Res(option%nflowdof-option%nmechdof) = eng/option%flow_dt
 

end subroutine THMCAccumulation

! ************************************************************************** !
!
! THMCFluxDerivative1: Computes the derivatives of the internal flux terms
!                     for the Jacobian
! author: Satish Karra
! date: 3/2/12
!
! ************************************************************************** ! 
subroutine THMCFluxDerivative1(aux_var_up,global_aux_var_up,por_up,tor_up, &
                             sir_up,dd_up,perm_up,Dk_up, &
                             aux_var_dn,global_aux_var_dn,por_dn,tor_dn, &
                             sir_dn,dd_dn,perm_dn,Dk_dn, &
                             area,dist_gravity,upweight, &
                             option,sat_func_up,sat_func_dn, &
                             Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
                             Dk_ice_up,Dk_ice_dn, &
                             alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                             disp_grad_up,disp_grad_dn,unit_normal, &
                             youngs_modulus_up,poissons_ratio_up, &
                             youngs_modulus_dn,poissons_ratio_dn, &
                             Jup,Jdn)
                             
  use Option_module 
  use Saturation_Function_module             
  use Water_EOS_module       
  
  implicit none
  
  type(thmc_auxvar_type) :: aux_var_up, aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn
  PetscReal :: Dk_ice_up, Dk_ice_dn
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: Diff_up, Diff_dn
  PetscReal :: v_darcy, area
  PetscReal :: dist_gravity  ! distance along gravity vector
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: ddifff_dp_up, ddifff_dp_dn, ddifff_dt_up, ddifff_dt_dn
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn, dden_ave_dt_up, dden_ave_dt_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn, dphi_dt_up, dphi_dt_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn, dukvr_dt_up, dukvr_dt_dn
  PetscReal :: duh_dp_up, duh_dp_dn, duh_dt_up, duh_dt_dn
  PetscReal :: dq_dp_up, dq_dp_dn, dq_dt_up, dq_dt_dn
  PetscReal :: duxmol_dxmol_up, duxmol_dxmol_dn
  
  PetscReal :: Dk_eff_up, Dk_eff_dn
  PetscReal :: Ke_up,Ke_dn   ! unfrozen soil Kersten numbers 
  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal :: dKe_dt_up, dKe_dp_up
  PetscReal :: dKe_dt_dn, dKe_dp_dn
  PetscReal :: dDk_dt_up, dDk_dt_dn
  PetscReal :: dDk_dp_up, dDk_dp_dn

  PetscInt :: iphase, ideriv
  type(thmc_auxvar_type) :: aux_var_pert_up, aux_var_pert_dn
  type(global_auxvar_type) :: global_aux_var_pert_up, global_aux_var_pert_dn
  PetscReal :: x_up(option%nflowdof-option%nmechdof), &
               x_dn(option%nflowdof-option%nmechdof), &
               x_pert_up(option%nflowdof-option%nmechdof), &
               x_pert_dn(option%nflowdof-option%nmechdof), pert_up, pert_dn
  PetscReal :: res(option%nflowdof), res_pert_up(option%nflowdof), &
               res_pert_dn(option%nflowdof)
  PetscReal :: J_pert_up(option%nflowdof,option%nflowdof), &
               J_pert_dn(option%nflowdof,option%nflowdof)
  PetscReal :: disp_grad_up(3,3), disp_grad_dn(3,3)
  PetscReal :: unit_normal(3)  
  PetscReal :: youngs_modulus_up, youngs_modulus_dn
  PetscReal :: poissons_ratio_up, poissons_ratio_dn


#ifdef ICE  
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: ddeng_dt_up, ddeng_dt_dn
  PetscReal :: dpsat_dt_up, dpsat_dt_dn
  PetscReal :: dmolg_dt_up, dmolg_dt_dn
  PetscReal :: dDiffg_dt_up, dDiffg_dt_dn
  PetscReal :: dDiffg_dp_up, dDiffg_dp_dn
  PetscReal :: dsatg_dp_up, dsatg_dp_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: Ke_fr_up,Ke_fr_dn   ! frozen soil Kersten numbers
  PetscReal :: dKe_fr_dt_up, dKe_fr_dt_dn
  PetscReal :: dKe_fr_dp_up, dKe_fr_dp_dn
#endif
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up*tor_up*por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0 
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dt_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  ddifff_dp_up = 0.d0
  ddifff_dp_dn = 0.d0
  ddifff_dt_up = 0.d0
  ddifff_dt_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_up = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_up = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_up = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_up = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_up = 0.d0
  dq_dt_dn = 0.d0
  duxmol_dxmol_up = 0.d0
  duxmol_dxmol_dn = 0.d0
  dDk_dt_up = 0.d0
  dDk_dt_dn = 0.d0
  dDk_dp_up = 0.d0
  dDk_dp_dn = 0.d0
  
! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_aux_var_up%den(1)+(1.D0-upweight)*global_aux_var_dn%den(1)
    dden_ave_dp_up = upweight*aux_var_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*aux_var_dn%dden_dp
    dden_ave_dt_up = upweight*aux_var_up%dden_dt
    dden_ave_dt_dn = (1.D0-upweight)*aux_var_dn%dden_dt

    gravity = (upweight*global_aux_var_up%den(1)*aux_var_up%avgmw + &
              (1.D0-upweight)*global_aux_var_dn%den(1)*aux_var_dn%avgmw) &
              * dist_gravity
    dgravity_dden_up = upweight*aux_var_up%avgmw*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity

    dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*aux_var_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp
    dphi_dt_up = dgravity_dden_up*aux_var_up%dden_dt
    dphi_dt_dn = dgravity_dden_dn*aux_var_dn%dden_dt

! note uxmol only contains one phase xmol
    if (dphi>=0.D0) then
      ukvr = aux_var_up%kvr
      dukvr_dp_up = aux_var_up%dkvr_dp
      dukvr_dt_up = aux_var_up%dkvr_dt
      
      uh = aux_var_up%h
      duh_dp_up = aux_var_up%dh_dp
      duh_dt_up = aux_var_up%dh_dt
      
      uxmol(1:option%nflowspec) = aux_var_up%xmol(1:option%nflowspec)
      duxmol_dxmol_up = 1.d0
    else
      ukvr = aux_var_dn%kvr
      dukvr_dp_dn = aux_var_dn%dkvr_dp
      dukvr_dt_dn = aux_var_dn%dkvr_dt
      
      uh = aux_var_dn%h
      duh_dp_dn = aux_var_dn%dh_dp
      duh_dt_dn = aux_var_dn%dh_dt
      
      uxmol(1:option%nflowspec) = aux_var_dn%xmol(1:option%nflowspec)
      duxmol_dxmol_dn = 1.d0
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      dq_dt_up = Dq*(dukvr_dt_up*dphi+ukvr*dphi_dt_up)*area
      dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jup(1,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)

      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
      Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)

      Jup(2,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uxmol(2)
      Jup(2,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uxmol(2)
      Jup(2,3) = q*density_ave*duxmol_dxmol_up

      Jdn(2,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(2)
      Jdn(2,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(2)
      Jdn(2,3) = q*density_ave*duxmol_dxmol_dn
     
      ! based on flux = q*density_ave*uh
      Jup(option%nflowdof-option%nmechdof,1) = (dq_dp_up*density_ave+ &
                                   q*dden_ave_dp_up)*uh+q*density_ave*duh_dp_up
      Jup(option%nflowdof-option%nmechdof,2) = (dq_dt_up*density_ave+ &
                                   q*dden_ave_dt_up)*uh+q*density_ave*duh_dt_up

      Jdn(option%nflowdof-option%nmechdof,1) = (dq_dp_dn*density_ave+ &
                                   q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn
      Jdn(option%nflowdof-option%nmechdof,2) = (dq_dt_dn*density_ave+ &
                                   q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn

    endif
  endif 

    difff = diffdp*0.25D0*(global_aux_var_up%sat(1) + global_aux_var_dn%sat(1))* &
                            (global_aux_var_up%den(1) + global_aux_var_dn%den(1))
    ddifff_dp_up = diffdp*0.25D0*(aux_var_up%dsat_dp*(global_aux_var_up%den(1) + global_aux_var_dn%den(1))+ &
                  (global_aux_var_up%sat(1) + global_aux_var_dn%sat(1))*aux_var_up%dden_dp)
    ddifff_dt_up = diffdp*0.25D0*(global_aux_var_up%sat(1) + global_aux_var_dn%sat(1))*aux_var_up%dden_dt

    ddifff_dp_dn = diffdp*0.25D0*(aux_var_dn%dsat_dp*(global_aux_var_up%den(1) + global_aux_var_dn%den(1))+ &
                  (global_aux_var_up%sat(1) + global_aux_var_dn%sat(1))*aux_var_dn%dden_dp)
    ddifff_dt_dn = diffdp*0.25D0*(global_aux_var_up%sat(1) + global_aux_var_dn%sat(1))*aux_var_dn%dden_dt
                                    
! SK   

    Jup(2,1) = Jup(2,1) + ddifff_dp_up*0.5d0*(Diff_up + Diff_dn)* &
                         (aux_var_up%xmol(2) - aux_var_dn%xmol(2))
    Jup(2,2) = Jup(2,2) + ddifff_dt_up*0.5d0*(Diff_up + Diff_dn)* &
                         (aux_var_up%xmol(2) - aux_var_dn%xmol(2)) 
    Jup(2,3) = Jup(2,3) + difff*0.5d0*(Diff_up + Diff_dn)

    Jdn(2,1) = Jdn(2,1) + ddifff_dp_dn*0.5d0*(Diff_up + Diff_dn)* &
                         (aux_var_up%xmol(2) - aux_var_dn%xmol(2))
    Jdn(2,2) = Jdn(2,2) + ddifff_dt_dn*0.5d0*(Diff_up + Diff_dn)* &
                         (aux_var_up%xmol(2) - aux_var_dn%xmol(2))
    Jdn(2,3) = Jdn(2,3) + difff*0.5d0*(Diff_up + Diff_dn)*(-1.d0)


#ifdef ICE
  ! Added by Satish Karra, updated 11/11/11
  satg_up = aux_var_up%sat_gas
  satg_dn = aux_var_dn%sat_gas
  if ((satg_up > eps) .and. (satg_dn > eps)) then
  p_g = option%reference_pressure  ! set to reference pressure
  deng_up = p_g/(IDEAL_GAS_CONST*(global_aux_var_up%temp(1) + 273.15d0))*1.d-3
  deng_dn = p_g/(IDEAL_GAS_CONST*(global_aux_var_dn%temp(1) + 273.15d0))*1.d-3
    
  Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
  p_ref = 1.01325d5   ! in Pa
  T_ref = 25.d0       ! in deg C
    
  Diffg_up = Diffg_ref*(p_ref/p_g)*((global_aux_var_up%temp(1) + 273.15d0)/ &
             (T_ref + 273.15d0))**(1.8)  
  Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_aux_var_dn%temp(1) + 273.15d0)/ &
             (T_ref + 273.15d0))**(1.8)    
  Ddiffgas_up = por_up*tor_up*satg_up*deng_up*Diffg_up
  Ddiffgas_dn = por_dn*tor_dn*satg_dn*deng_dn*Diffg_dn
  call PSAT(global_aux_var_up%temp(1), psat_up, dpsat_dt_up, ierr)
  call PSAT(global_aux_var_dn%temp(1), psat_dn, dpsat_dt_dn, ierr)
  molg_up = psat_up/p_g
  molg_dn = psat_dn/p_g
  ddeng_dt_up = - p_g/(IDEAL_GAS_CONST*(global_aux_var_up%temp(1) + &
                  273.15d0)**2)*1.d-3
  dmolg_dt_up = (1/p_g)*dpsat_dt_up
  ddeng_dt_dn = - p_g/(IDEAL_GAS_CONST*(global_aux_var_dn%temp(1) + &
                  273.15d0)**2)*1.d-3
  dmolg_dt_dn = (1/p_g)*dpsat_dt_dn
  dDiffg_dt_up = 1.8*Diffg_up/(global_aux_var_up%temp(1) + 273.15d0)
  dDiffg_dt_dn = 1.8*Diffg_dn/(global_aux_var_dn%temp(1) + 273.15d0)
  dDiffg_dp_up = 0.d0
  dDiffg_dp_dn = 0.d0
  dsatg_dp_up = aux_var_up%dsat_gas_dp
  dsatg_dp_dn = aux_var_dn%dsat_gas_dp
     
  if (molg_up > molg_dn) then 
    upweight = 0.d0
  else 
    upweight = 1.d0
  endif 

  Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn 
#if 0
  Jup(1,1) = Jup(1,1) + upweight*(por_up*tor_up*Diffg_up*dsatg_dp_up* &
             deng_up + por_up*tor_up*satg_up*deng_up*dDiffg_dp_up)* &
             (molg_up - molg_dn)/(dd_up + dd_dn)*area 

  Jup(1,2) = Jup(1,2) + upweight*(por_up*tor_up*Diffg_up*satg_up* &
             ddeng_dt_up + por_up*satg_up*tor_up*deng_up*dDiffg_dt_up)* &
             (molg_up - molg_dn)/(dd_up + dd_dn)*area + Ddiffgas_avg* &
             dmolg_dt_up/(dd_up + dd_dn)*area
    
  Jdn(1,1) = Jdn(1,1) + (1.D0 - upweight)*(por_dn*tor_dn*Diffg_dn*dsatg_dp_dn* &
             deng_dn + por_dn*tor_dn*satg_dn*deng_dn*dDiffg_dp_dn)* &
             (molg_up - molg_dn)/(dd_up + dd_dn)*area
             
  Jdn(1,2) = Jdn(1,2) + (1.D0 - upweight)*(por_dn*tor_dn*Diffg_dn*satg_dn* &
             ddeng_dt_dn + por_dn*satg_dn*tor_dn*deng_dn*dDiffg_dp_dn)* &
             (molg_up - molg_dn)/(dd_up + dd_dn)*area + Ddiffgas_avg* &
             (-dmolg_dt_dn)/(dd_up + dd_dn)*area
#endif
#if 1
  Jup(1,1) = Jup(1,1) + upweight*por_up*tor_up*deng_up*(Diffg_up*dsatg_dp_up &
              + satg_up*dDiffg_dp_up)* &
             (molg_up - molg_dn)/(dd_up + dd_dn)*area 

  Jup(1,2) = Jup(1,2) + (upweight*por_up*tor_up*satg_up*(Diffg_up* &
             ddeng_dt_up + deng_up*dDiffg_dt_up)*(molg_up - molg_dn) &
             + Ddiffgas_avg*dmolg_dt_up)/(dd_up + dd_dn)*area
  
  Jdn(1,1) = Jdn(1,1) + (1.D0 - upweight)*por_dn*tor_dn*deng_dn* &
             (Diffg_dn*dsatg_dp_dn + satg_dn*dDiffg_dp_dn)* &
             (molg_up - molg_dn)/(dd_up + dd_dn)*area
              
  Jdn(1,2) = Jdn(1,2) + ((1.D0 - upweight)*por_dn*tor_dn*satg_dn*(Diffg_dn* &
             ddeng_dt_dn + deng_dn*dDiffg_dp_dn)*(molg_up - molg_dn) &
             + Ddiffgas_avg*(-dmolg_dt_dn))/(dd_up + dd_dn)*area
#endif
  endif


  ddifff_dt_up = diffdp*0.25d0*(aux_var_up%dsat_dt)*&
                 (global_aux_var_up%den(1) + global_aux_var_dn%den(1))
  ddifff_dt_dn = diffdp*0.25d0*(aux_var_dn%dsat_dt)*&
                 (global_aux_var_up%den(1) + global_aux_var_dn%den(1))

  Jup(2,2) = Jup(2,2) + ddifff_dt_up*0.5d0*(Diff_up + Diff_dn)* &
                         (aux_var_up%xmol(2) - aux_var_dn%xmol(2))
  Jdn(2,2) = Jdn(2,2) + ddifff_dt_dn*0.5d0*(Diff_up + Diff_dn)* &
                         (aux_var_up%xmol(2) - aux_var_dn%xmol(2))

#endif 

        
! conduction term  
  Ke_up = (global_aux_var_up%sat(1) + epsilon)**(alpha_up)   !unfrozen soil Kersten number
  Ke_dn = (global_aux_var_dn%sat(1) + epsilon)**(alpha_dn)
  
  dKe_dp_up = alpha_up*(global_aux_var_up%sat(1) + epsilon)**(alpha_up - 1.d0)* &
               aux_var_up%dsat_dp
  dKe_dp_dn = alpha_dn*(global_aux_var_dn%sat(1) + epsilon)**(alpha_dn - 1.d0)* &
               aux_var_dn%dsat_dp
  
#ifdef ICE
            
  Ke_fr_up = (aux_var_up%sat_ice + epsilon)**(alpha_fr_up)
  Ke_fr_dn = (aux_var_dn%sat_ice + epsilon)**(alpha_fr_dn)

  Dk_eff_up = Dk_up*Ke_up + Dk_ice_up*Ke_fr_up + &
              (1.d0 - Ke_up - Ke_fr_up)*Dk_dry_up
  Dk_eff_dn = Dk_dn*Ke_dn + Dk_ice_dn*Ke_fr_dn + &
              (1.d0 - Ke_dn - Ke_fr_dn)*Dk_dry_dn

  dKe_dt_up = alpha_up*(global_aux_var_up%sat(1) + epsilon)**(alpha_up - 1.d0)* &
               aux_var_up%dsat_dt
  dKe_dt_dn = alpha_dn*(global_aux_var_dn%sat(1) + epsilon)**(alpha_dn - 1.d0)* &
               aux_var_dn%dsat_dt
               
  dKe_fr_dt_up = alpha_fr_up*(global_aux_var_up%sat(1) + epsilon)** &
                 (alpha_fr_up - 1.d0)*aux_var_up%dsat_dt
  dKe_fr_dt_dn = alpha_fr_dn*(global_aux_var_dn%sat(1) + epsilon)** &
                 (alpha_fr_dn - 1.d0)*aux_var_dn%dsat_dt
                 
  dKe_fr_dp_up = alpha_fr_up*(global_aux_var_up%sat(1) + epsilon)** &
                 (alpha_fr_up - 1.d0)*aux_var_up%dsat_dp
  dKe_fr_dp_dn = alpha_fr_dn*(global_aux_var_dn%sat(1) + epsilon)** &
                 (alpha_fr_dn - 1.d0)*aux_var_dn%dsat_dp
                 
#else 

  Dk_eff_up = Dk_dry_up + (Dk_up - Dk_dry_up)*Ke_up
  Dk_eff_dn = Dk_dry_dn + (Dk_dn - Dk_dry_dn)*Ke_dn 
  
  dKe_dt_up = 0.d0
  dKe_dt_dn = 0.d0
       
#endif
 
  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)
  
#ifdef ICE

  dDk_dt_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up*dKe_dt_up + &
              Dk_ice_up*dKe_fr_dt_up + (- dKe_dt_up - dKe_fr_dt_up)* &
              Dk_dry_up)
  dDk_dt_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn*dKe_dt_dn + &
              Dk_ice_dn*dKe_fr_dt_dn + (- dKe_dt_dn - dKe_fr_dt_dn)* &
              Dk_dry_dn)
              
  dDk_dp_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up*dKe_dp_up + &
              Dk_ice_up*dKe_fr_dp_up + (- dKe_dp_up - dKe_fr_dp_up)* &
              Dk_dry_up)
              
  dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn*dKe_dp_dn + &
              Dk_ice_dn*dKe_fr_dp_dn + (- dKe_dp_dn - dKe_fr_dp_dn)* &
              Dk_dry_dn)  

#else
  
  dDk_dt_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dt_up
  dDk_dt_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dt_dn
  
  dDk_dp_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dp_up
  dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dp_dn

#endif  
    
  !  cond = Dk*area*(global_aux_var_up%temp(1)-global_aux_var_dn%temp(1)) 
  Jup(option%nflowdof-option%nmechdof,1) = Jup(option%nflowdof-option%nmechdof,1) + &
                           area*(global_aux_var_up%temp(1) - &
                           global_aux_var_dn%temp(1))*dDk_dp_up
  Jdn(option%nflowdof-option%nmechdof,1) = Jdn(option%nflowdof-option%nmechdof,1) + &
                           area*(global_aux_var_up%temp(1) - &
                           global_aux_var_dn%temp(1))*dDk_dp_dn
                           
  Jup(option%nflowdof-option%nmechdof,2) = Jup(option%nflowdof-option%nmechdof,2) + &
                           Dk*area + &
                           area*(global_aux_var_up%temp(1) - & 
                           global_aux_var_dn%temp(1))*dDk_dt_up 
  Jdn(option%nflowdof-option%nmechdof,2) = Jdn(option%nflowdof-option%nmechdof,2) + &
                           Dk*area*(-1.d0) + &
                           area*(global_aux_var_up%temp(1) - & 
                           global_aux_var_dn%temp(1))*dDk_dt_dn 
                           
 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives_flow) then
    allocate(aux_var_pert_up%xmol(option%nflowspec),aux_var_pert_up%diff(option%nflowspec))
    allocate(aux_var_pert_dn%xmol(option%nflowspec),aux_var_pert_dn%diff(option%nflowspec))
    call GlobalAuxVarInit(global_aux_var_pert_up,option)
    call GlobalAuxVarInit(global_aux_var_pert_dn,option)  
    call THMCAuxVarCopy(aux_var_up,aux_var_pert_up,option)
    call THMCAuxVarCopy(aux_var_dn,aux_var_pert_dn,option)
    call GlobalAuxVarCopy(global_aux_var_up,global_aux_var_pert_up,option)
    call GlobalAuxVarCopy(global_aux_var_dn,global_aux_var_pert_dn,option)
    x_up(1) = global_aux_var_up%pres(1)
    x_up(2) = global_aux_var_up%temp(1)
    x_up(3) = aux_var_up%xmol(2)
    x_dn(1) = global_aux_var_dn%pres(1)
    x_dn(2) = global_aux_var_dn%temp(1)
    x_dn(3) = aux_var_dn%xmol(2)
    
    J_pert_up = 0.d0
    J_pert_dn = 0.d0

    call THMCFlux( &
      aux_var_up,global_aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
      aux_var_dn,global_aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
      area,dist_gravity,upweight, &
      option,v_darcy,Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
      Dk_ice_up,Dk_ice_dn, &
      alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
      disp_grad_up,disp_grad_dn,unit_normal, &
      youngs_modulus_up,poissons_ratio_up, &
      youngs_modulus_dn,poissons_ratio_dn, &
      res)
    do ideriv = 1,option%nflowdof-option%nmechdof
      pert_up = x_up(ideriv)*perturbation_tolerance
      pert_dn = x_dn(ideriv)*perturbation_tolerance
      x_pert_up = x_up
      x_pert_dn = x_dn

#ifdef ICE
      
      if (ideriv == 1) then
        if (x_pert_up(ideriv) < option%reference_pressure) then
          pert_up = - pert_up
        endif
          x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
      
        if (x_pert_dn(ideriv) < option%reference_pressure) then
          pert_dn = - pert_dn
        endif
          x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
      
      if (ideriv == 2) then
        if (x_pert_up(ideriv) < 0.d0) then
          pert_up = - 1.d-5
        else
          pert_up = 1.d-5
        endif
          x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
        
        if (x_pert_dn(ideriv) < 0.d0) then
          pert_dn = - 1.d-5
        else
          pert_dn = 1.d-5
        endif
          x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
      
      if (ideriv == 3) then
        x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
        x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
      
#else
      
      x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
      x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn

#endif

#ifdef ICE
      call THMCAuxVarComputeIce(x_pert_up,aux_var_pert_up, &
                            global_aux_var_pert_up, &
                            iphase,sat_func_up, &
                            0.d0,0.d0,option)
      call THMCAuxVarComputeIce(x_pert_dn,aux_var_pert_dn, &
                            global_aux_var_pert_dn,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
#else
      call THMCAuxVarCompute(x_pert_up,aux_var_pert_up, &
                            global_aux_var_pert_up, &
                            iphase,sat_func_up, &
                            0.d0,0.d0,option)
      call THMCAuxVarCompute(x_pert_dn,aux_var_pert_dn, &
                            global_aux_var_pert_dn,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
#endif

      call THMCFlux(aux_var_pert_up,global_aux_var_pert_up, &
                   por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                   aux_var_dn,global_aux_var_dn, &
                   por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                   area,dist_gravity,upweight, &
                   option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   disp_grad_up,disp_grad_dn,unit_normal, &
                   youngs_modulus_up,poissons_ratio_up, &
                   youngs_modulus_dn,poissons_ratio_dn, &
                   res_pert_up)
      call THMCFlux(aux_var_up,global_aux_var_up, &
                   por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                   aux_var_pert_dn,global_aux_var_pert_dn, &
                   por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                   area,dist_gravity,upweight, &
                   option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   disp_grad_up,disp_grad_dn,unit_normal, &
                   youngs_modulus_up,poissons_ratio_up, &
                   youngs_modulus_dn,poissons_ratio_dn, &
                   res_pert_dn)
      J_pert_up(:,ideriv) = (res_pert_up(:) - res(:))/pert_up
      J_pert_dn(:,ideriv) = (res_pert_dn(:) - res(:))/pert_dn
    enddo
    deallocate(aux_var_pert_up%xmol,aux_var_pert_up%diff)
    deallocate(aux_var_pert_dn%xmol,aux_var_pert_dn%diff)
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_aux_var_pert_up)
    call GlobalAuxVarStrip(global_aux_var_pert_dn)    
  endif

end subroutine THMCFluxDerivative1


! ************************************************************************** !
!
! THMCFluxDerivative2: Computes the derivatives of the internal flux terms
!                     for the Jacobian
! author: Satish Karra
! date: 4/13/12
!
! ************************************************************************** ! 
subroutine THMCFluxDerivative2(aux_var_up,global_aux_var_up,por_up,tor_up, &
                             sir_up,dd_up,perm_up,Dk_up, &
                             aux_var_dn,global_aux_var_dn,por_dn,tor_dn, &
                             sir_dn,dd_dn,perm_dn,Dk_dn, &
                             area,dist_gravity,upweight, &
                             option,sat_func_up,sat_func_dn, &
                             Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
                             Dk_ice_up,Dk_ice_dn, &
                             alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                             disp_grad_up,disp_grad_dn,unit_normal, &
                             youngs_modulus_up,poissons_ratio_up, &
                             youngs_modulus_dn,poissons_ratio_dn, &
                             disp_grad_pert_up,disp_grad_pert_dn, &
                             pert_up,pert_dn, &
                             Jup,Jdn)
                             
  use Option_module 
  use Saturation_Function_module             
  use Water_EOS_module       
  
  implicit none
  
  type(thmc_auxvar_type) :: aux_var_up, aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn
  PetscReal :: Dk_ice_up, Dk_ice_dn
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: Diff_up, Diff_dn
  PetscReal :: v_darcy, area
  PetscReal :: dist_gravity  ! distance along gravity vector
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscReal :: Jup(option%nflowdof)
  PetscReal :: Jdn(option%nflowdof)
     
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
    
  PetscReal :: Dk_eff_up, Dk_eff_dn
  PetscReal :: Ke_up,Ke_dn   ! unfrozen soil Kersten numbers 
  PetscReal, parameter :: epsilon = 1.d-6

  PetscInt :: iphase, ideriv
  type(thmc_auxvar_type) :: aux_var_pert_up, aux_var_pert_dn
  type(global_auxvar_type) :: global_aux_var_pert_up, global_aux_var_pert_dn
  PetscReal :: res(option%nflowdof), res_pert_up(option%nflowdof), &
               res_pert_dn(option%nflowdof)
  PetscReal :: disp_grad_up(3,3), disp_grad_dn(3,3)
  PetscReal :: disp_grad_pert_up(3,3), disp_grad_pert_dn(3,3)
  PetscReal :: unit_normal(3)  
  PetscReal :: youngs_modulus_up, youngs_modulus_dn
  PetscReal :: poissons_ratio_up, poissons_ratio_dn
  PetscReal :: pert_up, pert_dn
  
  Jup = 0.d0
  Jdn = 0.d0 


  call THMCFlux( &
      aux_var_up,global_aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
      aux_var_dn,global_aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
      area,dist_gravity,upweight, &
      option,v_darcy,Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
      Dk_ice_up,Dk_ice_dn, &
      alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
      disp_grad_up,disp_grad_dn,unit_normal, &
      youngs_modulus_up,poissons_ratio_up, &
      youngs_modulus_dn,poissons_ratio_dn, &
      res)

  call THMCFlux( &
      aux_var_up,global_aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
      aux_var_dn,global_aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
      area,dist_gravity,upweight, &
      option,v_darcy,Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
      Dk_ice_up,Dk_ice_dn, &
      alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
      disp_grad_pert_up,disp_grad_dn,unit_normal, &
      youngs_modulus_up,poissons_ratio_up, &
      youngs_modulus_dn,poissons_ratio_dn, &
      res_pert_up)

  call THMCFlux( &
      aux_var_up,global_aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
      aux_var_dn,global_aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
      area,dist_gravity,upweight, &
      option,v_darcy,Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
      Dk_ice_up,Dk_ice_dn, &
      alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
      disp_grad_up,disp_grad_pert_dn,unit_normal, &
      youngs_modulus_up,poissons_ratio_up, &
      youngs_modulus_dn,poissons_ratio_dn, &
      res_pert_dn)

  Jup = (res_pert_up - res)/pert_up
  Jdn = (res_pert_dn - res)/pert_dn

end subroutine THMCFluxDerivative2


! ************************************************************************** !
!
! THMCFluxDerivativeAnalytical: Computes the internal flux terms for the residual
! author:
! date: 8/21/12
!
! ************************************************************************** ! 
subroutine THMCFluxDerivativeAnalytical(grid, &
                  ghosted_id_up, ghosted_id_dn, &
                  aux_var_up,global_aux_var_up, &
                  aux_var_dn,global_aux_var_dn, &
                  dd_up,dd_dn, &
                  area,option, &
                  disp_grad_up,disp_grad_dn,normal, &
                  youngs_modulus_up,poissons_ratio_up, &
                  youngs_modulus_dn,poissons_ratio_dn, &
                  Jup,Jdn)
                  
  use Option_module                              
  use Water_EOS_module
  use Grid_module

  implicit none
  
#include "finclude/petscdmda.h"

  type(grid_type), pointer :: grid
  type(thmc_auxvar_type) :: aux_var_up, aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, dd_dn
  PetscReal :: area
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof) 
  PetscReal :: disp_grad_up(3,3), disp_grad_dn(3,3)
  PetscReal :: normal(3)
  PetscReal :: disp_grad_face(3,3)
  PetscReal :: stress(3,3)
  PetscReal :: youngs_modulus_up, youngs_modulus_dn
  PetscReal :: poissons_ratio_up, poissons_ratio_dn
  PetscReal :: youngs_modulus_avg, poissons_ratio_avg
  PetscReal :: lambda_avg, mu_avg
  PetscInt :: i, j
  PetscReal :: Minv_up(3,3), Minv_dn(3,3)
  PetscInt :: ghosted_neighbors_size
  PetscInt :: ghosted_neighbors(26)
  PetscReal :: disp_vec_sum_up(3), disp_vec_sum_dn(3)
  PetscReal :: disp_up_minus_dn(3), disp_dn_minus_up(3)
  PetscInt :: ghosted_id_up, ghosted_id_dn
 
         
!  Jup = 0.d0
!  Jdn = 0.d0
  
  disp_vec_sum_up = 0.d0
  disp_vec_sum_dn = 0.d0
  
 ! Calculating the residual for the mechanical component

  youngs_modulus_avg = (youngs_modulus_up*youngs_modulus_dn)/ &
                       (dd_dn*youngs_modulus_up + dd_up*youngs_modulus_dn)* &
                       (dd_up + dd_dn)
  poissons_ratio_avg = (poissons_ratio_up*poissons_ratio_dn)/ &
                       (dd_dn*poissons_ratio_up + dd_up*poissons_ratio_dn)* &
                       (dd_up + dd_dn)

  disp_grad_face = dd_up*disp_grad_up + dd_dn*disp_grad_dn

  lambda_avg = youngs_modulus_avg*poissons_ratio_avg/((1.d0 + poissons_ratio_avg)* &
                                    (1.d0 - 2.d0*poissons_ratio_avg))
  mu_avg = youngs_modulus_avg/(2.d0*(1.d0 + poissons_ratio_avg))

  call THMCComputeMinv(grid, ghosted_id_up, Minv_up, option)
  call THMCComputeMinv(grid, ghosted_id_dn, Minv_dn, option)
    
  call GridGetGhostedNeighborsWithCorners(grid,ghosted_id_up, &
                                         DMDA_STENCIL_STAR, &
                                         ONE_INTEGER,ONE_INTEGER,ONE_INTEGER, &
                                         ghosted_neighbors_size, &
                                         ghosted_neighbors, &
                                         option)     
                                         
  do i = 1, ghosted_neighbors_size
    disp_vec_sum_up(1) = disp_vec_sum_up(1) + grid%x(ghosted_neighbors(i)) - grid%x(ghosted_id_up)
    disp_vec_sum_up(2) = disp_vec_sum_up(2) + grid%y(ghosted_neighbors(i)) - grid%y(ghosted_id_up)
    disp_vec_sum_up(3) = disp_vec_sum_up(3) + grid%z(ghosted_neighbors(i)) - grid%z(ghosted_id_up)
  enddo  
  
  call GridGetGhostedNeighborsWithCorners(grid,ghosted_id_dn, &
                                         DMDA_STENCIL_STAR, &
                                         ONE_INTEGER,ONE_INTEGER,ONE_INTEGER, &
                                         ghosted_neighbors_size, &
                                         ghosted_neighbors, &
                                         option)
                                         
  do i = 1, ghosted_neighbors_size
    disp_vec_sum_dn(1) = disp_vec_sum_dn(1) + grid%x(ghosted_neighbors(i)) - grid%x(ghosted_id_up)
    disp_vec_sum_dn(2) = disp_vec_sum_dn(2) + grid%y(ghosted_neighbors(i)) - grid%y(ghosted_id_up)
    disp_vec_sum_dn(3) = disp_vec_sum_dn(3) + grid%z(ghosted_neighbors(i)) - grid%z(ghosted_id_up)
  enddo                                      
                                        
                                         
  disp_up_minus_dn(1) = grid%x(ghosted_id_up) - grid%x(ghosted_id_dn) 
  disp_up_minus_dn(2) = grid%y(ghosted_id_up) - grid%y(ghosted_id_dn)                                         
  disp_up_minus_dn(3) = grid%z(ghosted_id_up) - grid%z(ghosted_id_dn)        
  
  disp_dn_minus_up = -disp_up_minus_dn                                                                  

! ======== Calculating mechanical components of Jup ==========================   
                                             
  do i = 1, option%nmechdof
    do j = 1, option%nmechdof
      Jup(i+option%nflowdof-option%nmechdof,j+option%nflowdof-option%nmechdof)  = &
      lambda_avg*(-dd_up*(Minv_up(1,j)*disp_vec_sum_up(1) + Minv_up(2,j)* &
      disp_vec_sum_up(2) + Minv_up(3,j)*disp_vec_sum_up(3)) + &
      dd_dn*(Minv_dn(1,j)*disp_up_minus_dn(1) + Minv_dn(2,j)* &
      disp_up_minus_dn(2) + Minv_dn(3,j)*disp_up_minus_dn(3)))*normal(i)*area
    enddo
  enddo
      
  do i = 1, option%nmechdof
    do j = 1, option%nmechdof
      Jup(i+option%nflowdof-option%nmechdof,j+option%nflowdof-option%nmechdof)  = &
      mu_avg*(-dd_up*Minv_up(i,j)*(disp_vec_sum_up(1)*normal(1) + &
      disp_vec_sum_up(2)*normal(2) + disp_vec_sum_up(3)*normal(3)) + &
      dd_dn*Minv_dn(i,j)*(disp_up_minus_dn(1)*normal(1) + &
      disp_up_minus_dn(2)*normal(2) + disp_up_minus_dn(3)*normal(3)))*area
    enddo
  enddo      
      
  do i = 1, option%nmechdof
    do j = 1, option%nmechdof
      Jup(i+option%nflowdof-option%nmechdof,j+option%nflowdof-option%nmechdof)  = &
      mu_avg*(-dd_up*disp_vec_sum_up(i)*(Minv_up(1,j)*normal(1) + &
      Minv_up(2,j)*normal(2) + Minv_up(3,j)*normal(3)) + &
      dd_dn*disp_up_minus_dn(i)*(Minv_dn(1,j)*normal(1) + &
      Minv_dn(2,j)*normal(2) + Minv_dn(3,j)*normal(3)))*area
    enddo
  enddo  

! ======== Calculating mechanical components of Jdn ==========================                                                

  do i = 1, option%nmechdof
    do j = 1, option%nmechdof
      Jdn(i+option%nflowdof-option%nmechdof,j+option%nflowdof-option%nmechdof)  = &
      lambda_avg*(-dd_dn*(Minv_dn(1,j)*disp_vec_sum_dn(1) + Minv_dn(2,j)* &
      disp_vec_sum_dn(2) + Minv_dn(3,j)*disp_vec_sum_dn(3)) + &
      dd_up*(Minv_up(1,j)*disp_dn_minus_up(1) + Minv_up(2,j)* &
      disp_dn_minus_up(2) + Minv_up(3,j)*disp_dn_minus_up(3)))*normal(i)*area
    enddo
  enddo
      
  do i = 1, option%nmechdof
    do j = 1, option%nmechdof
      Jdn(i+option%nflowdof-option%nmechdof,j+option%nflowdof-option%nmechdof)  = &
      mu_avg*(-dd_dn*Minv_dn(i,j)*(disp_vec_sum_dn(1)*normal(1) + &
      disp_vec_sum_dn(2)*normal(2) + disp_vec_sum_dn(3)*normal(3)) + &
      dd_up*Minv_up(i,j)*(disp_dn_minus_up(1)*normal(1) + &
      disp_dn_minus_up(2)*normal(2) + disp_dn_minus_up(3)*normal(3)))*area
    enddo
  enddo      
      
  do i = 1, option%nmechdof
    do j = 1, option%nmechdof
      Jdn(i+option%nflowdof-option%nmechdof,j+option%nflowdof-option%nmechdof)  = &
      mu_avg*(-dd_dn*disp_vec_sum_dn(i)*(Minv_dn(1,j)*normal(1) + &
      Minv_dn(2,j)*normal(2) + Minv_dn(3,j)*normal(3)) + &
      dd_up*disp_dn_minus_up(i)*(Minv_up(1,j)*normal(1) + &
      Minv_up(2,j)*normal(2) + Minv_up(3,j)*normal(3)))*area
    enddo
  enddo  

                                  
 
end subroutine THMCFluxDerivativeAnalytical



! ************************************************************************** !
!
! THMCFlux: Computes the internal flux terms for the residual
! author:
! date: 3/2/12
!
! ************************************************************************** ! 
subroutine THMCFlux(aux_var_up,global_aux_var_up, &
                  por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                  aux_var_dn,global_aux_var_dn, &
                  por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                  area,dist_gravity,upweight, &
                  option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                  disp_grad_up,disp_grad_dn,unit_normal, &
                  youngs_modulus_up,poissons_ratio_up, &
                  youngs_modulus_dn,poissons_ratio_dn, &
                  Res)
                  
  use Option_module                              
  use Water_EOS_module

  implicit none
  
  type(thmc_auxvar_type) :: aux_var_up, aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn
  PetscReal :: Dk_ice_up, Dk_ice_dn
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: Dk_eff_up, Dk_eff_dn
  PetscReal :: Diff_up,Diff_dn
  PetscReal :: v_darcy,area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: Ke_up,Ke_dn   ! unfrozen soil Kersten numbers 
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal :: disp_grad_up(3,3), disp_grad_dn(3,3)
  PetscReal :: unit_normal(3)
  PetscReal :: disp_grad_face(3,3)
  PetscReal :: stress(3,3)
  PetscReal :: youngs_modulus_up, youngs_modulus_dn
  PetscReal :: poissons_ratio_up, poissons_ratio_dn
  PetscReal :: youngs_modulus_avg, poissons_ratio_avg

#ifdef ICE  
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: Ke_fr_up,Ke_fr_dn   ! frozen soil Kersten numbers
#endif
     
     
  Res = 0.d0
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up*tor_up*por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0  
  
! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) < eps) then 
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) < eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_aux_var_up%den(1)+(1.D0-upweight)*global_aux_var_dn%den(1) 

    gravity = (upweight*global_aux_var_up%den(1)*aux_var_up%avgmw + &
              (1.D0-upweight)*global_aux_var_dn%den(1)*aux_var_dn%avgmw) &
              * dist_gravity

    dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity

!   note uxmol only contains one component xmol
    if (dphi >= 0.D0) then
      ukvr = aux_var_up%kvr
      uh = aux_var_up%h
      uxmol(1:option%nflowspec) = aux_var_up%xmol(1:option%nflowspec)
    else
      ukvr = aux_var_dn%kvr
      uh = aux_var_dn%h
      uxmol(1:option%nflowspec) = aux_var_dn%xmol(1:option%nflowspec)
    endif

    if (ukvr > floweps) then
      v_darcy = Dq * ukvr * dphi
   
      q = v_darcy * area
        
      fluxm(1) = fluxm(1) + q*density_ave
      fluxm(2) = fluxm(2) + q*density_ave*uxmol(2)
      fluxe = fluxe + q*density_ave*uh 
    endif
  endif 

  
  difff = diffdp * 0.25D0*(global_aux_var_up%sat(1)+global_aux_var_dn%sat(1))* &
                            (global_aux_var_up%den(1)+global_aux_var_dn%den(1))
  fluxm(2) = fluxm(2) + difff * .5D0 * (Diff_up + Diff_dn)* &
                 (aux_var_up%xmol(2) - aux_var_dn%xmol(2)) 


#ifdef ICE
  ! Added by Satish Karra, 10/24/11
  satg_up = aux_var_up%sat_gas
  satg_dn = aux_var_dn%sat_gas
 if ((satg_up > eps) .and. (satg_dn > eps)) then
  p_g = option%reference_pressure ! set to reference pressure
  deng_up = p_g/(IDEAL_GAS_CONST*(global_aux_var_up%temp(1) + 273.15d0))*1.d-3
  deng_dn = p_g/(IDEAL_GAS_CONST*(global_aux_var_dn%temp(1) + 273.15d0))*1.d-3
    
  Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
  p_ref = 1.01325d5 ! in Pa
  T_ref = 25.d0 ! in deg C
    
  Diffg_up = Diffg_ref*(p_ref/p_g)*((global_aux_var_up%temp(1) + 273.15d0)/ &
             (T_ref + 273.15d0))**(1.8)  
  Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_aux_var_dn%temp(1) + 273.15d0)/ &
             (T_ref + 273.15d0))**(1.8)
  Ddiffgas_up = por_up*tor_up*satg_up*deng_up*Diffg_up
  Ddiffgas_dn = por_dn*tor_dn*satg_dn*deng_dn*Diffg_dn
  call PSAT(global_aux_var_up%temp(1), psat_up, ierr)
  call PSAT(global_aux_var_dn%temp(1), psat_dn, ierr)
  molg_up = psat_up/p_g
  molg_dn = psat_dn/p_g
  
  if (molg_up > molg_dn) then 
    upweight = 0.d0
  else 
    upweight = 1.d0
  endif
    
  Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn 
  fluxm(1) = fluxm(1) + Ddiffgas_avg*area*(molg_up - molg_dn)/ &
             (dd_up + dd_dn)
             
 endif
#endif 

! conduction term  
  Ke_up = (global_aux_var_up%sat(1) + epsilon)**(alpha_up)   !unfrozen soil Kersten number
  Ke_dn = (global_aux_var_dn%sat(1) + epsilon)**(alpha_dn)
     
#ifdef ICE

  Ke_fr_up = (aux_var_up%sat_ice + epsilon)**(alpha_fr_up)
  Ke_fr_dn = (aux_var_dn%sat_ice + epsilon)**(alpha_fr_dn)

  Dk_eff_up = Dk_up*Ke_up + Dk_ice_up*Ke_fr_up + &
              (1.d0 - Ke_up - Ke_fr_up)*Dk_dry_up
  Dk_eff_dn = Dk_dn*Ke_dn + Dk_ice_dn*Ke_fr_dn + &
              (1.d0 - Ke_dn - Ke_fr_dn)*Dk_dry_dn
#else

  Dk_eff_up = Dk_dry_up + (Dk_up - Dk_dry_up)*Ke_up
  Dk_eff_dn = Dk_dry_dn + (Dk_dn - Dk_dry_dn)*Ke_dn      

#endif
 
  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)
  cond = Dk*area*(global_aux_var_up%temp(1) - global_aux_var_dn%temp(1)) 
  fluxe = fluxe + cond

  Res(1:option%nflowspec) = fluxm(:)
  Res(option%nflowdof-option%nmechdof) = fluxe
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

 ! Calculating the residual for the mechanical component

  youngs_modulus_avg = (youngs_modulus_up*youngs_modulus_dn)/(dd_dn*youngs_modulus_up + &
                        dd_up*youngs_modulus_dn)*(dd_up + dd_dn)
  poissons_ratio_avg = (poissons_ratio_up*poissons_ratio_dn)/(dd_dn*poissons_ratio_up + &
                        dd_up*poissons_ratio_dn)*(dd_up + dd_dn)

  disp_grad_face = dd_up*disp_grad_up + dd_dn*disp_grad_dn


  call THMCComputeStressFromDispGrad(disp_grad_face,youngs_modulus_avg, &
                                poissons_ratio_avg,stress)
                                  
  Res(option%nflowdof-option%nmechdof+1) = (stress(1,1)*unit_normal(1) + &
                                           stress(1,2)*unit_normal(2) + &
                                           stress(1,3)*unit_normal(3))*area

  Res(option%nflowdof-option%nmechdof+2) = (stress(2,1)*unit_normal(1) + &
                                           stress(2,2)*unit_normal(2) + &
                                           stress(2,3)*unit_normal(3))*area

  Res(option%nflowdof-option%nmechdof+3) = (stress(3,1)*unit_normal(1) + &
                                           stress(3,2)*unit_normal(2) + &
                                           stress(3,3)*unit_normal(3))*area
  
 
end subroutine THMCFlux

! ************************************************************************** !
!
! THMCBCFluxDerivative: Computes the derivatives of the boundary flux 
!                           terms for the Jacobian
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCBCFluxDerivative(ibndtype,aux_vars, &
                              aux_var_up,global_aux_var_up, &
                              aux_var_dn,global_aux_var_dn, &
                              por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                              area,dist_gravity,option, &
                              sat_func_dn,Jdn)
  use Option_module
  use Saturation_Function_module
  use Water_EOS_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(thmc_auxvar_type) :: aux_var_up, aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array in boundary condition
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn,Diff_dn
  PetscReal :: area
  type(saturation_function_type) :: sat_func_dn  
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: v_darcy
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: ddiff_dp_dn, ddiff_dt_dn
  PetscReal :: dden_ave_dp_dn, dden_ave_dt_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn, dphi_dt_dn
  PetscReal :: dukvr_dp_dn, dukvr_dt_dn
  PetscReal :: duh_dp_dn, duh_dt_dn
  PetscReal :: dq_dp_dn, dq_dt_dn
  PetscReal :: duxmol_dxmol_dn

  PetscInt :: iphase, ideriv
  type(thmc_auxvar_type) :: aux_var_pert_dn, aux_var_pert_up
  type(global_auxvar_type) :: global_aux_var_pert_dn, global_aux_var_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(3), x_up(3), x_pert_dn(3), x_pert_up(3), pert_dn, res(3), &
               res_pert_dn(3), J_pert_dn(3,3)

#ifdef ICE  
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: ddeng_dt_dn
  PetscReal :: dpsat_dt_dn
  PetscReal :: dmolg_dt_dn
  PetscReal :: dDiffg_dt_dn
  PetscReal :: dDiffg_dp_dn
  PetscReal :: dsatg_dp_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
#endif
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  ddiff_dp_dn = 0.d0
  ddiff_dt_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_dn = 0.d0
  duxmol_dxmol_dn = 0.d0
        
  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(THMC_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dq = perm_dn / dd_up
      ! Flow term
      if (global_aux_var_up%sat(1) > sir_dn .or. global_aux_var_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_aux_var_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_aux_var_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_aux_var_up%den(1)+(1.D0-upweight)*global_aux_var_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*aux_var_dn%dden_dp
        dden_ave_dt_dn = (1.D0-upweight)*aux_var_dn%dden_dt

        if (ibndtype(THMC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dt_dn = dden_ave_dt_dn + upweight*aux_var_up%dden_dt
        endif
        
        gravity = (upweight*global_aux_var_up%den(1)*aux_var_up%avgmw + &
                  (1.D0-upweight)*global_aux_var_dn%den(1)*aux_var_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity

        dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp
        dphi_dt_dn = dgravity_dden_dn*aux_var_dn%dden_dt

        if (ibndtype(THMC_PRESSURE_DOF) == SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_aux_var_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
            dphi_dp_dn = 0.d0
            dphi_dt_dn = 0.d0
          endif
        endif        
        
        if (ibndtype(THMC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                                   !( dgravity_dden_up                   ) (dden_dt_up)
          dphi_dt_dn = dphi_dt_dn + upweight*aux_var_up%avgmw*dist_gravity*aux_var_up%dden_dt
        endif
        
        if (dphi>=0.D0) then
          ukvr = aux_var_up%kvr
          if (ibndtype(THMC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dt_dn = aux_var_up%dkvr_dt
          endif
        else
          ukvr = aux_var_dn%kvr
          dukvr_dp_dn = aux_var_dn%dkvr_dp
          dukvr_dt_dn = aux_var_dn%dkvr_dt
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(THMC_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(THMC_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_aux_var_up%den(1)
          if (ibndtype(THMC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dden_ave_dt_dn = aux_var_up%dden_dt
          endif
        else 
          density_ave = global_aux_var_dn%den(1)
          dden_ave_dp_dn = aux_var_dn%dden_dp
          dden_ave_dt_dn = aux_var_dn%dden_dt
        endif 
        q = v_darcy * area
      endif

  end select

  if (v_darcy >= 0.D0) then
    uh = aux_var_up%h
    uxmol(:)=aux_var_up%xmol(1:option%nflowspec)
    if (ibndtype(THMC_PRESSURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dp_dn = aux_var_up%dh_dp
    endif
    if (ibndtype(THMC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dt_dn = aux_var_up%dh_dt
    endif
    if (ibndtype(THMC_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      duxmol_dxmol_dn = 1.d0
    endif
  else
    uh = aux_var_dn%h
    duh_dp_dn = aux_var_dn%dh_dp
    duh_dt_dn = aux_var_dn%dh_dt

    uxmol(:)=aux_var_dn%xmol(1:option%nflowspec)
    duxmol_dxmol_dn = 1.d0
  endif      

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
  Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)
!  Jdn(1,3:option%nflowdof) = 0.d0
  do ispec=2,option%nflowspec 
    ! based on flux = q*density_ave*uxmol
    Jdn(ispec,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(ispec)
    Jdn(ispec,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(ispec)
    Jdn(ispec,ispec+1) = q*density_ave*duxmol_dxmol_dn
  enddo
      ! based on flux = q*density_ave*uh
  Jdn(option%nflowdof-option%nmechdof,1) = (dq_dp_dn*density_ave+ &
                                   q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn
  Jdn(option%nflowdof-option%nmechdof,2) = (dq_dt_dn*density_ave+ &
                                   q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn
!  Jdn(option%nflowdof,3:option%nflowdof) = 0.d0

  ! Diffusion term   
  select case(ibndtype(THMC_CONCENTRATION_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC) 

      diff = diffdp * global_aux_var_dn%sat(1)*global_aux_var_dn%den(1)
      ddiff_dp_dn = diffdp * (aux_var_dn%dsat_dp*global_aux_var_dn%den(1)+ &
                                global_aux_var_dn%sat(1)*aux_var_dn%dden_dp)
      ddiff_dt_dn = diffdp * global_aux_var_dn%sat(1)*aux_var_dn%dden_dt
! SK
      Jdn(2,1) = Jdn(2,1) + ddiff_dp_dn*Diff_dn* &
                            (aux_var_up%xmol(2)-aux_var_dn%xmol(2))
      Jdn(2,2) = Jdn(2,2) + ddiff_dt_dn*Diff_dn* &
                            (aux_var_up%xmol(2)-aux_var_dn%xmol(2))
      Jdn(2,3) = Jdn(2,3) + diff*Diff_dn*(-1.d0)

  end select
   
  ! Conduction term
  select case(ibndtype(THMC_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dk =  Dk_dn / dd_up
      !cond = Dk*area*(global_aux_var_up%temp(1)-global_aux_var_dn%temp(1)) 
      Jdn(option%nflowdof-option%nmechdof,2) = & 
                         Jdn(option%nflowdof-option%nmechdof,2)+Dk*area*(-1.d0)
#ifdef ICE
      ! Added by Satish Karra, 11/21/11
      satg_up = aux_var_up%sat_gas
      satg_dn = aux_var_dn%sat_gas
      if ((satg_up > eps) .and. (satg_dn > eps)) then
        p_g = option%reference_pressure  ! set to reference pressure
        deng_up = p_g/(IDEAL_GAS_CONST*(global_aux_var_up%temp(1) + &
                  273.15d0))*1.d-3
        deng_dn = p_g/(IDEAL_GAS_CONST*(global_aux_var_dn%temp(1) + &
                  273.15d0))*1.d-3
        
        Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
        p_ref = 1.01325d5 ! in Pa
        T_ref = 25.d0 ! in deg C 
        
        Diffg_up = Diffg_ref*(p_ref/p_g)*((global_aux_var_up%temp(1) + &
                   273.15d0)/(T_ref + 273.15d0))**(1.8)  
        Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_aux_var_dn%temp(1) + &
                   273.15d0)/(T_ref + 273.15d0))**(1.8)  
        Ddiffgas_up = satg_up*deng_up*Diffg_up
        Ddiffgas_dn = satg_dn*deng_dn*Diffg_dn
        call PSAT(global_aux_var_up%temp(1), psat_up, ierr)
        call PSAT(global_aux_var_dn%temp(1), psat_dn, dpsat_dt_dn, ierr)
        molg_up = psat_up/p_g
        molg_dn = psat_dn/p_g
        ddeng_dt_dn = - p_g/(IDEAL_GAS_CONST*(global_aux_var_dn%temp(1) + &
                        273.15d0)**2)*1.d-3
        dmolg_dt_dn = (1/p_g)*dpsat_dt_dn
        dDiffg_dt_dn = 1.8*Diffg_dn/(global_aux_var_dn%temp(1) + 273.15d0)
        dDiffg_dp_dn = 0.d0
        dsatg_dp_dn = aux_var_dn%dsat_gas_dp
        
        if (molg_up > molg_dn) then 
          upweight = 0.d0
        else 
          upweight = 1.d0
        endif
        
        Ddiffgas_avg = upweight*Ddiffgas_up+(1.D0 - upweight)*Ddiffgas_dn 
    
        Jdn(1,1) = Jdn(1,1) + por_dn*tor_dn*(1.D0 - upweight)* &
                   Ddiffgas_dn/satg_dn*dsatg_dp_dn*(molg_up - molg_dn)/dd_up* &
                   area
        Jdn(1,2) = Jdn(1,2) + por_dn*tor_dn*(1.D0 - upweight)* &
                   (Ddiffgas_avg/deng_dn*ddeng_dt_dn + Ddiffgas_avg/Diffg_dn* &
                   dDiffg_dt_dn)*(molg_up - molg_dn)/dd_up*area + por_dn* &
                   tor_dn*Ddiffgas_avg*(-dmolg_dt_dn)/dd_up*area
      endif
#endif  

  end select

  if (option%numerical_derivatives_flow) then
    allocate(aux_var_pert_dn%xmol(option%nflowspec),aux_var_pert_dn%diff(option%nflowspec))
    allocate(aux_var_pert_up%xmol(option%nflowspec),aux_var_pert_up%diff(option%nflowspec))
    
    call GlobalAuxVarInit(global_aux_var_pert_up,option)
    call GlobalAuxVarInit(global_aux_var_pert_dn,option)  
    call THMCAuxVarCopy(aux_var_up,aux_var_pert_up,option)
    call THMCAuxVarCopy(aux_var_dn,aux_var_pert_dn,option)
    call GlobalAuxVarCopy(global_aux_var_up,global_aux_var_pert_up,option)
    call GlobalAuxVarCopy(global_aux_var_dn,global_aux_var_pert_dn,option)
    
    x_up(1) = global_aux_var_up%pres(1)
    x_up(2) = global_aux_var_up%temp(1)
    x_up(3) = aux_var_up%xmol(2)
    x_dn(1) = global_aux_var_dn%pres(1)
    x_dn(2) = global_aux_var_dn%temp(1)
    x_dn(3) = aux_var_dn%xmol(2)
    do ideriv = 1,3
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_up(ideriv) = x_dn(ideriv)
      endif
    enddo
#ifdef ICE
      call THMCAuxVarComputeIce(x_dn,aux_var_dn, &
                            global_aux_var_dn,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
      call THMCAuxVarComputeIce(x_up,aux_var_up, &
                            global_aux_var_up,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
#else
      call THMCAuxVarCompute(x_dn,aux_var_dn, &
                            global_aux_var_dn,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
      call THMCAuxVarCompute(x_up,aux_var_up, &
                            global_aux_var_up,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
#endif
    
    call THMCBCFlux(ibndtype,aux_vars,aux_var_up,global_aux_var_up, &
                  aux_var_dn,global_aux_var_dn, &
                  por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                  area,dist_gravity,option,v_darcy,Diff_dn,res)
    if (ibndtype(THMC_PRESSURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(THMC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(THMC_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    do ideriv = 1,3
      pert_dn = x_dn(ideriv)*perturbation_tolerance    
      x_pert_dn = x_dn
     
#ifdef ICE
      
      if (ideriv == 1) then
        if (x_pert_dn(ideriv) < option%reference_pressure) then
          pert_dn = - pert_dn
        endif
        x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
      
      if (ideriv == 2) then
        if (x_pert_dn(ideriv) < 0.d0) then
           pert_dn = - 1.d-5
        else
           pert_dn = 1.d-5
        endif
        x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
      
      if (ideriv == 3) then
        x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
      
#else
      
      x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn

#endif    
        
      x_pert_up = x_up
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_pert_up(ideriv) = x_pert_dn(ideriv)
      endif   

#ifdef ICE
      call THMCAuxVarComputeIce(x_pert_dn,aux_var_pert_dn, &
                            global_aux_var_pert_dn,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
      call THMCAuxVarComputeIce(x_pert_up,aux_var_pert_up, &
                            global_aux_var_pert_up,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
#else
      call THMCAuxVarCompute(x_pert_dn,aux_var_pert_dn, &
                            global_aux_var_pert_dn,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
      call THMCAuxVarCompute(x_pert_up,aux_var_pert_up, &
                            global_aux_var_pert_up,iphase,sat_func_dn, &
                            0.d0,0.d0,option)
#endif

      call THMCBCFlux(ibndtype,aux_vars,aux_var_pert_up,global_aux_var_pert_up, &
                    aux_var_pert_dn,global_aux_var_pert_dn, &
                    por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                    area,dist_gravity,option,v_darcy,Diff_dn,res_pert_dn)
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    deallocate(aux_var_pert_dn%xmol,aux_var_pert_dn%diff)
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_aux_var_pert_up)
    call GlobalAuxVarStrip(global_aux_var_pert_dn)      
  endif

end subroutine THMCBCFluxDerivative

! ************************************************************************** !
!
! THMCBCFlux: Computes the  boundary flux terms for the residual
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCBCFlux(ibndtype,aux_vars,aux_var_up,global_aux_var_up, &
                    aux_var_dn,global_aux_var_dn, &
                    por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                    area,dist_gravity,option,v_darcy,Diff_dn, & 
                    Res)
!                  disp_grad_dn,unit_normal,
!                  youngs_modulus_up,poissons_ratio_up, &
!                  youngs_modulus_dn,poissons_ratio_dn, &Res)
  use Option_module
  use Water_EOS_module 
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(thmc_auxvar_type) :: aux_var_up, aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn, Diff_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: v_darcy, area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  
#ifdef ICE  
  PetscReal :: Ddiffgas_avg, Ddiffgas_dn, Ddiffgas_up
  PetscReal :: p_g
  PetscReal :: deng_dn, deng_up
  PetscReal :: psat_dn, psat_up
  PetscReal :: molg_dn, molg_up
  PetscReal :: satg_dn, satg_up
  PetscReal :: Diffg_dn, Diffg_up
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
#endif
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  Res = 0.d0

  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(THMC_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dq = perm_dn / dd_up
      ! Flow term
      if (global_aux_var_up%sat(1) > sir_dn .or. global_aux_var_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_aux_var_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_aux_var_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_aux_var_up%den(1)+(1.D0-upweight)*global_aux_var_dn%den(1)
   
        gravity = (upweight*global_aux_var_up%den(1)*aux_var_up%avgmw + &
                  (1.D0-upweight)*global_aux_var_dn%den(1)*aux_var_dn%avgmw) &
                  * dist_gravity

        dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity

        if (ibndtype(THMC_PRESSURE_DOF) == SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_aux_var_up%pres(1) - option%reference_pressure < eps) then
            dphi = 0.d0
          endif
        endif
        
        if (dphi>=0.D0) then
          ukvr = aux_var_up%kvr
        else
          ukvr = aux_var_dn%kvr
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(THMC_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(THMC_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_aux_var_up%den(1)
        else 
          density_ave = global_aux_var_dn%den(1)
        endif 
      endif

  end select

  q = v_darcy * area

  if (v_darcy >= 0.D0) then
    uh = aux_var_up%h
    uxmol(:)=aux_var_up%xmol(1:option%nflowspec)
  else
    uh = aux_var_dn%h
    uxmol(:)=aux_var_dn%xmol(1:option%nflowspec)
  endif      
    
  do ispec=1, option%nflowspec 
    fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
  enddo
  fluxe = fluxe + q*density_ave*uh

  ! Diffusion term   
  select case(ibndtype(THMC_CONCENTRATION_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
    
!      if (global_aux_var_up%sat > eps .and. global_aux_var_dn%sat > eps) then
!        diff = diffdp * 0.25D0*(global_aux_var_up%sat(1)+global_aux_var_dn%sat(1))* &
!               (global_aux_var_up%den+global_aux_var_dn%den)

!pcl  if (global_aux_var_dn%sat(1) > eps) then
        diff = diffdp * global_aux_var_dn%sat(1)*global_aux_var_dn%den(1)
 !         fluxm(ispec) = fluxm(ispec) + diff*aux_var_dn%diff(ispec)* &
!                          (aux_var_up%xmol(ispec)-aux_var_dn%xmol(ispec))
        fluxm(2) = fluxm(2) + diff*Diff_dn* &
                           (aux_var_up%xmol(2)-aux_var_dn%xmol(2))
!pcl  endif
  end select
  


  ! Conduction term
  select case(ibndtype(THMC_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dk = Dk_dn / dd_up
      cond = Dk*area*(global_aux_var_up%temp(1)-global_aux_var_dn%temp(1)) 
      fluxe = fluxe + cond

#ifdef ICE
  ! Added by Satish Karra,
      satg_up = aux_var_up%sat_gas
      satg_dn = aux_var_dn%sat_gas
 if ((satg_up > eps) .and. (satg_dn > eps)) then
      p_g = option%reference_pressure ! set to reference pressure
      deng_up = p_g/(IDEAL_GAS_CONST*(global_aux_var_up%temp(1) + 273.15d0))*1.d-3
      deng_dn = p_g/(IDEAL_GAS_CONST*(global_aux_var_dn%temp(1) + 273.15d0))*1.d-3
        
      Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
      p_ref = 1.01325d5 ! in Pa
      T_ref = 25.d0 ! in deg C
        
      Diffg_up = Diffg_ref*(p_ref/p_g)*((global_aux_var_up%temp(1) + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
      Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_aux_var_dn%temp(1) + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
      Ddiffgas_up = satg_up*deng_up*Diffg_up
      Ddiffgas_dn = satg_dn*deng_dn*Diffg_dn
      call PSAT(global_aux_var_up%temp(1), psat_up, ierr)
      call PSAT(global_aux_var_dn%temp(1), psat_dn, ierr)
      molg_up = psat_up/p_g
      molg_dn = psat_dn/p_g
        
      if (molg_up > molg_dn) then 
        upweight = 0.d0
      else 
        upweight = 1.d0
      endif
        
      Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn 
      fluxm(1) = fluxm(1) + por_dn*tor_dn*Ddiffgas_avg*(molg_up - molg_dn)/ &
                 dd_up*area
  endif
#endif 

    case(NEUMANN_BC)
      fluxe = fluxe + aux_vars(THMC_TEMPERATURE_DOF)*area*option%scale ! added by SK 10/18/11
    case(ZERO_GRADIENT_BC)
      ! No change in fluxe
  end select

  Res(1:option%nflowspec) = fluxm(:)
  Res(option%nflowdof-option%nmechdof) = fluxe

end subroutine THMCBCFlux

! ************************************************************************** !
!
! THMCResidual: Computes the residual equation 
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCResidual(snes,xx,r,realization,ierr)

  use Realization_class
  use Level_module
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
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
  ! These 3 must be called before THMCUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  
  ! Compute internal and boundary flux terms as well as source/sink terms
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call THMCResidualPatch(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THMCResidual


! ************************************************************************** !
!
! THMCResidualPatch: Computes the residual equation at patch level
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCResidualPatch(snes,xx,r,realization,ierr)

  use Water_EOS_module

  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use Utility_module

  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: i, jn
  PetscInt :: ip1, ip2
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
                          
               
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: iphase
  PetscInt :: icap_up, icap_dn, ithrm_up, ithrm_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dd, f_up, f_dn, ff
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: Dk_dry_up, Dk_dry_dn ! dry thermal conductivities
  PetscReal :: Dk_ice_up, Dk_ice_dn ! frozen soil thermal conductivities
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: Diff_up, Diff_dn  ! "Diffusion" constants at upstream, downstream faces.
  PetscReal :: dw_kg, dw_mol
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(thmc_parameter_type), pointer :: thmc_parameter
  type(thmc_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: normal(3), unit_normal(3)
  PetscReal :: disp_grad_up(3,3), disp_grad_dn(3,3)
  PetscReal :: youngs_modulus_up, youngs_modulus_dn
  PetscReal :: poissons_ratio_up, poissons_ratio_dn
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  thmc_parameter => patch%aux%THMC%thmc_parameter
  aux_vars => patch%aux%THMC%aux_vars
  aux_vars_bc => patch%aux%THMC%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  
  call THMCUpdateAuxVarsPatch(realization)
  ! override flags since they will soon be out of date  
  patch%aux%THMC%aux_vars_up_to_date = PETSC_FALSE
  if (option%compute_mass_balance_new) then
    call THMCZeroMassBalDeltaPatch(realization)
  endif


! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90( r, r_p, ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)
 
  call VecGetArrayF90(field%flow_yy,yy_p,ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'

  r_p = 0.d0
#if 1
  ! Accumulation terms ------------------------------------

  r_p = -accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    call THMCAccumulation(aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                        porosity_loc_p(ghosted_id), &
                        volume_p(local_id), &
                        thmc_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                        thmc_parameter%rock_den(int(ithrm_loc_p(ghosted_id))), &
                        option,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)

    r_p(iend-2) =  r_p(iend-2) + &
         thmc_parameter%rock_den(int(ithrm_loc_p(ghosted_id)))* &
         option%gravity(1)*volume_p(local_id)*1.d-6 !convert to MN
    r_p(iend-1) =  r_p(iend-2) + &
         thmc_parameter%rock_den(int(ithrm_loc_p(ghosted_id)))* &
         option%gravity(2)*volume_p(local_id)*1.d-6 !convert to MN
    r_p(iend) =  r_p(iend) + &
         thmc_parameter%rock_den(int(ithrm_loc_p(ghosted_id)))* &
         option%gravity(3)*volume_p(local_id)*1.d-6 !convert to MN

  enddo 

#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%flow_condition%num_sub_conditions > THMC_CONCENTRATION_DOF) then
      enthalpy_flag = PETSC_TRUE
    else
      enthalpy_flag = PETSC_FALSE
    endif

    qsrc1 = source_sink%flow_condition%pressure%flow_dataset%time_series%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%flow_dataset%time_series%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%flow_dataset%time_series%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%flow_dataset%time_series%cur_value(1)

    qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / FMWCO2
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      if (enthalpy_flag) then
        r_p(local_id*(option%nflowdof-option%nmechdof)) = &
                         r_p(local_id*(option%nflowdof-option%nmechdof)) - &
                         hsrc1
      endif         

!      if (qsrc1 > 0.d0) then ! injection
!        call wateos_noderiv(tsrc1,global_aux_vars(ghosted_id)%pres(1), &
!                            dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
!        r_p((local_id-1)*option%nflowdof + jh2o) = r_p((local_id-1)*option%nflowdof + jh2o) &
!                                               - qsrc1
!        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - qsrc1*enth_src_h2o
!      endif  
!    
!      if (csrc1 > 0.d0) then ! injection
!        call printErrMsg(option,"concentration source not yet implemented in THMC")
!      endif
!  !  else if (qsrc1 < 0.d0) then ! withdrawal
!  !  endif
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

      normal(1) = grid%x(ghosted_id_dn) - grid%x(ghosted_id_up)
      normal(2) = grid%y(ghosted_id_dn) - grid%y(ghosted_id_up)
      normal(3) = grid%z(ghosted_id_dn) - grid%z(ghosted_id_up)

      unit_normal = normal/sqrt(DotProduct(normal,normal))

      disp_grad_up = aux_vars(ghosted_id_up)%gradient
      disp_grad_dn = aux_vars(ghosted_id_dn)%gradient
        
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
   
      D_up = thmc_parameter%ckwet(ithrm_up)
      D_dn = thmc_parameter%ckwet(ithrm_dn)
      
      Dk_dry_up = thmc_parameter%ckdry(ithrm_up)
      Dk_dry_dn = thmc_parameter%ckdry(ithrm_dn)
      
      alpha_up = thmc_parameter%alpha(ithrm_up)
      alpha_dn = thmc_parameter%alpha(ithrm_dn)

      youngs_modulus_up = thmc_parameter%youngs_modulus(ithrm_up)
      youngs_modulus_dn = thmc_parameter%youngs_modulus(ithrm_dn)

      poissons_ratio_up = thmc_parameter%poissons_ratio(ithrm_up)
      poissons_ratio_dn = thmc_parameter%poissons_ratio(ithrm_dn)

#ifdef ICE
      Dk_ice_up = thmc_parameter%ckfrozen(ithrm_up)
      DK_ice_dn = thmc_parameter%ckfrozen(ithrm_dn)
      
      alpha_fr_up = thmc_parameter%alpha_fr(ithrm_up)
      alpha_fr_dn = thmc_parameter%alpha_fr(ithrm_dn)
#else
      Dk_ice_up = Dk_dry_up
      Dk_ice_dn = Dk_dry_dn
      
      alpha_fr_up = alpha_up
      alpha_fr_dn = alpha_dn
#endif

      Diff_up = thmc_parameter%diffusion_coefficient(1)
      Diff_dn = thmc_parameter%diffusion_coefficient(1)

      call THMCFlux(aux_vars(ghosted_id_up),global_aux_vars(ghosted_id_up), &
                  porosity_loc_p(ghosted_id_up), &
                  tor_loc_p(ghosted_id_up),thmc_parameter%sir(1,icap_up), &
                  dd_up,perm_up,D_up, &
                  aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
                  porosity_loc_p(ghosted_id_dn), &
                  tor_loc_p(ghosted_id_dn),thmc_parameter%sir(1,icap_dn), &
                  dd_dn,perm_dn,D_dn, &
                  cur_connection_set%area(iconn),distance_gravity, &
                  upweight,option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                  disp_grad_up,disp_grad_dn,unit_normal, &
                  youngs_modulus_up,poissons_ratio_up, &
                  youngs_modulus_dn,poissons_ratio_dn, &
                  Res)
                  
                  

      patch%internal_velocities(1,sum_connection) = v_darcy
     
      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = thmc_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))

      call THMCBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                aux_vars_bc(sum_connection), &
                                global_aux_vars_bc(sum_connection), &
                                aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id), &
                                thmc_parameter%sir(1,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
                                cur_connection_set%area(iconn), &
                                distance_gravity,option, &
                                v_darcy,Diff_dn,Res)
      patch%boundary_velocities(1,sum_connection) = v_darcy

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo  
#endif 
 
  if (option%use_isothermal) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      istart = 3 + (local_id-1)*option%nflowdof
      r_p(istart)=xx_loc_p(2 + (ghosted_id-1)*option%nflowdof)-yy_p(istart-1)
    enddo
  endif
  
  if (patch%aux%THMC%inactive_cells_exist) then
    do i=1,patch%aux%THMC%n_zero_rows
      r_p(patch%aux%THMC%zero_rows_local(i)) = 0.d0
    enddo
  endif
    
  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(option%mycomm,'THMCresidual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(option%mycomm,'THMCxx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif 

end subroutine THMCResidualPatch

! ************************************************************************** !
!
! THMCJacobian: Computes the Jacobian for a realization
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_class
  use Patch_module
  use Level_module
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
  type(option_type),  pointer :: option
  PetscReal :: norm
  
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
      call THMCJacobianPatch(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'THMCjacobian.out', &
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
  
end subroutine THMCJacobian

! ************************************************************************** !
!
! THMCJacobianPatch: Computes the Jacobian on a patch
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCJacobianPatch(snes,xx,A,B,flag,realization,ierr)
       
  use Water_EOS_module

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Utility_module

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i
  PetscInt :: ip1, ip2 

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  PetscInt :: icap,iphas,iphas_up,iphas_dn,icap_up,icap_dn
  PetscInt :: ii, jj
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  PetscReal :: tsrc1,qsrc1,csrc1,hsrc1
  PetscReal :: dd_up, dd_dn, dd, f_up, f_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: dw_dp,dw_dt,hw_dp,hw_dt,dresT_dp,dresT_dt
  PetscReal :: D_up, D_dn  ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: Dk_dry_up, Dk_dry_dn ! dry thermal conductivities
  PetscReal :: Dk_ice_up, Dk_ice_dn ! frozen soil thermal conductivities
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: Diff_up, Diff_dn ! "Diffusion" constants upstream and downstream of a face.
  PetscReal :: zero, norm
  PetscReal :: upweight
  PetscReal :: max_dev  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
            Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(thmc_parameter_type), pointer :: thmc_parameter
  type(thmc_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 

  PetscReal :: normal(3), unit_normal(3)
  PetscReal :: disp_grad_up(3,3), disp_grad_dn(3,3)
  PetscReal :: disp_grad_pert_up(3,3), disp_grad_pert_dn(3,3)
  PetscReal :: youngs_modulus_up, youngs_modulus_dn
  PetscReal :: poissons_ratio_up, poissons_ratio_dn
  PetscReal :: pert_up, pert_dn
  
  PetscViewer :: viewer
  Vec :: debug_vec

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  thmc_parameter => patch%aux%THMC%thmc_parameter
  aux_vars => patch%aux%THMC%aux_vars
  aux_vars_bc => patch%aux%THMC%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  
#if 0
   call THMCNumericalJacobianTest(xx,realization)
#endif

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  pert_up = perturbation_tolerance
  pert_dn = perturbation_tolerance


#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    icap = int(icap_loc_p(ghosted_id))

    call THMCAccumDerivative(aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              thmc_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              thmc_parameter%rock_den(int(ithrm_loc_p(ghosted_id))), &
                              option, &
                              realization%saturation_function_array(icap)%ptr, &
                              Jup) 
        
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_accum.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%flow_condition%num_sub_conditions > THMC_CONCENTRATION_DOF) then
      enthalpy_flag = PETSC_TRUE
    else
      enthalpy_flag = PETSC_FALSE
    endif

    qsrc1 = source_sink%flow_condition%pressure%flow_dataset%time_series%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%flow_dataset%time_series%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%flow_dataset%time_series%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%flow_dataset%time_series%cur_value(1)

    qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / FMWCO2
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
!      if (enthalpy_flag) then
!        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1   
!      endif         

      if (qsrc1 > 0.d0) then ! injection
        call wateos(tsrc1,global_aux_vars(ghosted_id)%pres(1),dw_kg,dw_mol,dw_dp,dw_dt, &
              enth_src_h2o,hw_dp,hw_dt,option%scale,ierr)        
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        ! base on r_p() = r_p() - qsrc1*enth_src_h2o
        dresT_dp = -qsrc1*hw_dp
        ! dresT_dt = -qsrc1*hw_dt ! since tsrc1 is prescribed, there is no derivative
        istart = ghosted_id*option%nflowdof
        call MatSetValuesLocal(A,1,istart-1,1,istart-option%nflowdof,dresT_dp,ADD_VALUES,ierr)
        ! call MatSetValuesLocal(A,1,istart-1,1,istart-1,dresT_dt,ADD_VALUES,ierr)
      endif  
    
      if (csrc1 > 0.d0) then ! injection
        call printErrMsg(option,"concentration source not yet implemented in THMC")
      endif
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
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

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or. &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

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
      
      normal(1) = grid%x(ghosted_id_dn) - grid%x(ghosted_id_up)
      normal(2) = grid%y(ghosted_id_dn) - grid%y(ghosted_id_up)
      normal(3) = grid%z(ghosted_id_dn) - grid%z(ghosted_id_up)

      unit_normal = normal/sqrt(DotProduct(normal,normal))

      disp_grad_up = aux_vars(ghosted_id_up)%gradient
      disp_grad_dn = aux_vars(ghosted_id_dn)%gradient

      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))
    
      iphas_up = int(iphase_loc_p(ghosted_id_up))
      iphas_dn = int(iphase_loc_p(ghosted_id_dn))

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      
      D_up = thmc_parameter%ckwet(ithrm_up)
      D_dn = thmc_parameter%ckwet(ithrm_dn)
    
      Dk_dry_up = thmc_parameter%ckdry(ithrm_up)
      Dk_dry_dn = thmc_parameter%ckdry(ithrm_dn)
      
      alpha_up = thmc_parameter%alpha(ithrm_up)
      alpha_dn = thmc_parameter%alpha(ithrm_dn)

      youngs_modulus_up = thmc_parameter%youngs_modulus(ithrm_up)
      youngs_modulus_dn = thmc_parameter%youngs_modulus(ithrm_dn)

      poissons_ratio_up = thmc_parameter%poissons_ratio(ithrm_up)
      poissons_ratio_dn = thmc_parameter%poissons_ratio(ithrm_dn)

#ifdef ICE
      Dk_ice_up = thmc_parameter%ckfrozen(ithrm_up)
      DK_ice_dn = thmc_parameter%ckfrozen(ithrm_dn)
      
      alpha_fr_up = thmc_parameter%alpha_fr(ithrm_up)
      alpha_fr_dn = thmc_parameter%alpha_fr(ithrm_dn)
#else
      Dk_ice_up = Dk_dry_up
      Dk_ice_dn = Dk_dry_dn
      
      alpha_fr_up = alpha_up
      alpha_fr_dn = alpha_dn
#endif

      Diff_up = thmc_parameter%diffusion_coefficient(1)
      Diff_dn = thmc_parameter%diffusion_coefficient(1)

      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
                              
      call THMCFluxDerivative1(aux_vars(ghosted_id_up),global_aux_vars(ghosted_id_up), &
                             porosity_loc_p(ghosted_id_up), &
                             tor_loc_p(ghosted_id_up),thmc_parameter%sir(1,icap_up), &
                             dd_up,perm_up,D_up, &
                             aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
                             porosity_loc_p(ghosted_id_dn), &
                             tor_loc_p(ghosted_id_dn),thmc_parameter%sir(1,icap_dn), &
                             dd_dn,perm_dn,D_dn, &
                             cur_connection_set%area(iconn),distance_gravity, &
                             upweight,option, &
                             realization%saturation_function_array(icap_up)%ptr, &
                             realization%saturation_function_array(icap_dn)%ptr, &
                             Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
                             Dk_ice_up,Dk_ice_dn, &
                             alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                             disp_grad_up,disp_grad_dn,unit_normal, &
                             youngs_modulus_up,poissons_ratio_up, &
                             youngs_modulus_dn,poissons_ratio_dn, &
                             Jup,Jdn)


!      do i = 1,realization%option%nmechdof
!          call THMCComputeDisplacementGradientPert(grid, global_aux_vars, &
!                                               ghosted_id_up, &
!                                               disp_grad_pert_up, i, &
!                                               pert_up, option)
!          call THMCComputeDisplacementGradientPert(grid, global_aux_vars, &
!                                               ghosted_id_dn, &
!                                               disp_grad_pert_dn, i, &
!                                               pert_dn, option)
!          call THMCFluxDerivative2(aux_vars(ghosted_id_up), &
!                             global_aux_vars(ghosted_id_up), &
!                             porosity_loc_p(ghosted_id_up), &
!                             tor_loc_p(ghosted_id_up),thmc_parameter%sir(1,icap_up), &
!                             dd_up,perm_up,D_up, &
!                             aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
!                             porosity_loc_p(ghosted_id_dn), &
!                             tor_loc_p(ghosted_id_dn),thmc_parameter%sir(1,icap_dn), &
!                             dd_dn,perm_dn,D_dn, &
!                             cur_connection_set%area(iconn),distance_gravity, &
!                             upweight,option, &
!                             realization%saturation_function_array(icap_up)%ptr, &
!                             realization%saturation_function_array(icap_dn)%ptr, &
!                             Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
!                             Dk_ice_up,Dk_ice_dn, &
!                             alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
!                             disp_grad_up,disp_grad_dn,unit_normal, &
!                             youngs_modulus_up,poissons_ratio_up, &
!                             youngs_modulus_dn,poissons_ratio_dn, &
!                             disp_grad_pert_up,disp_grad_pert_dn, &
!                             pert_up, pert_dn, &
!                             Jup(:,realization%option%nflowdof-realization%option%nmechdof+i), &
!                             Jdn(:,realization%option%nflowdof-realization%option%nmechdof+i)) 
!                                        
!                             
!      enddo


     call THMCFluxDerivativeAnalytical(grid, &
                  ghosted_id_up, ghosted_id_dn, &
                  aux_vars(ghosted_id_up),global_aux_vars(ghosted_id_up), &
                  aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
                  dd_up,dd_dn, &
                  cur_connection_set%area(iconn),option, &
                  disp_grad_up,disp_grad_dn,unit_normal, &
                  youngs_modulus_up,poissons_ratio_up, &
                  youngs_modulus_dn,poissons_ratio_dn, &
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
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_flux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = thmc_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      icap_dn = int(icap_loc_p(ghosted_id))  

      call THMCBCFluxDerivative(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                aux_vars_bc(sum_connection), &
                                global_aux_vars_bc(sum_connection), &
                                aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id), &
                                thmc_parameter%sir(1,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
                                cur_connection_set%area(iconn), &
                                distance_gravity,option, &
                                realization%saturation_function_array(icap_dn)%ptr,&
                                Jdn)
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  do i=1, n_zero_rows
    ii = mod(zero_rows_local(i),option%nflowdof)
    ip1 = zero_rows_local_ghosted(i)
    if (ii == 0) then
      ip2 = ip1-1
    else if (ii == option%nflowdof-1) then
      ip2 = ip1+1
    else
      ip2 = ip1
    endif
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
#else
  if (patch%aux%THMC%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%THMC%n_zero_rows, &
                          patch%aux%THMC%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif
#endif

end subroutine THMCJacobianPatch

! ************************************************************************** !
!
! THMCCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCCreateZeroArray(patch,option)

  use Patch_module
  use Grid_module
  use Option_module
  
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

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
        n_zero_rows = n_zero_rows + option%nflowdof
      else
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
        n_zero_rows = n_zero_rows + 1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
    n_zero_rows = n_zero_rows + grid%nlmax
#endif
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
        do idof = 1, option%nflowdof
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*option%nflowdof+idof
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%nflowdof+idof-1
        enddo
      else
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
        ncount = ncount + 1
        zero_rows_local(ncount) = local_id*option%nflowdof
        zero_rows_local_ghosted(ncount) = ghosted_id*option%nflowdof-1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ncount = ncount + 1
      zero_rows_local(ncount) = local_id*option%nflowdof
      zero_rows_local_ghosted(ncount) = ghosted_id*option%nflowdof-1
    enddo
#endif
  endif

  patch%aux%THMC%zero_rows_local => zero_rows_local
  patch%aux%THMC%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%THMC%n_zero_rows = n_zero_rows

  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
  if (flag > 0) patch%aux%THMC%inactive_cells_exist = PETSC_TRUE

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine THMCCreateZeroArray

! ************************************************************************** !
!
! THMCMaxChange: Computes the maximum change in the solution vector
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCMaxChange(realization)

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

  option%dcmax=0.D0
  
  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dtmpmax,ierr)
  if (option%nflowdof > 2) &
    call VecStrideNorm(field%flow_dxx,TWO_INTEGER,NORM_INFINITY,option%dcmax,ierr)
    
end subroutine THMCMaxChange

! ************************************************************************** !
!
! THMCResidualToMass: Computes mass balance from residual equation
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCResidualToMass(realization)

  use Realization_class
  use Level_module
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module

  implicit none

  Vec :: ts_mass_balance
  type(realization_type) :: realization
  
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  PetscReal, pointer :: mass_balance_p(:)
  type(thmc_auxvar_type), pointer :: aux_vars(:) 
  type(global_auxvar_type), pointer :: global_aux_vars(:) 
  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  
  option => realization%option
  field => realization%field

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid
      aux_vars => cur_patch%aux%THMC%aux_vars

      call VecGetArrayF90(field%flow_ts_mass_balance,mass_balance_p, ierr)
  
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (cur_patch%imat(ghosted_id) <= 0) cycle

        istart = (ghosted_id-1)*option%nflowdof+1
        mass_balance_p(istart) = mass_balance_p(istart)/ &
                                 global_aux_vars(ghosted_id)%den(1)* &
                                 global_aux_vars(ghosted_id)%den_kg(1)
      enddo

      call VecRestoreArrayF90(field%flow_ts_mass_balance,mass_balance_p, ierr)

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THMCResidualToMass

! ************************************************************************** !
! 
! THMCComputeStressFromDispGrad: Computes the stress from given displacement
! gradient
! Author: Satish Karra
! Date: 3/20/12
! Stress is in MPa
!
! ************************************************************************** !


subroutine THMCComputeStressFromDispGrad(disp_grad,youngs_modulus, &
                                         poissons_ratio,stress) 


  use Utility_module

  implicit none

  PetscReal :: disp_grad(3,3)
  PetscReal :: youngs_modulus, poissons_ratio
  PetscReal :: lambda, mu
  PetscReal :: stress(3,3)
  PetscReal :: identity(3,3)
  PetscReal :: trace_disp_grad
  
  lambda = youngs_modulus*poissons_ratio/((1.d0 + poissons_ratio)* &
                                    (1.d0 - 2.d0*poissons_ratio))
  mu = youngs_modulus/(2.d0*(1.d0 + poissons_ratio))
  identity = reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))
  trace_disp_grad = DotProduct(identity(1,:),disp_grad(1,:)) + &
                    DotProduct(identity(2,:),disp_grad(2,:)) + &
                    DotProduct(identity(3,:),disp_grad(3,:))
  stress = mu*(disp_grad + transpose(disp_grad)) + lambda*trace_disp_grad* &
                                                          identity
  stress = stress*1.d-6 ! convert to MPa
  
end subroutine THMCComputeStressFromDispGrad

! ************************************************************************** !
! 
! THMCComputeDisplacementGradient: Computes the gradient of displacement using
! least square fit of values from neighboring cells
! See:I. Bijelonja, I. Demirdzic, S. Muzaferija -- A finite volume method 
! for incompressible linear elasticity, CMAME
! Author: Satish Karra
! Date: 2/20/12
!
! ************************************************************************** !


subroutine THMCComputeDisplacementGradient(grid, global_aux_vars, ghosted_id, &
                                           gradient, Minv, option) 


  use Grid_module
  use Global_Aux_module
  use Option_module
  use Utility_module

  implicit none
#include "finclude/petscdmda.h"

  type(option_type) :: option
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_aux_vars(:)

  
  PetscInt :: ghosted_neighbors_size, ghosted_id
  PetscInt :: ghosted_neighbors(26)
  PetscReal :: gradient(3,3), disp_vec(3), disp_mat(3,3)
  PetscReal :: ux_weighted(3)
  PetscReal :: uy_weighted(3)
  PetscReal :: uz_weighted(3)
  PetscReal :: Minv(3,3)
  PetscInt :: i
  
   
  call GridGetGhostedNeighborsWithCorners(grid,ghosted_id, &
                                         DMDA_STENCIL_STAR, &
                                         ONE_INTEGER,ONE_INTEGER,ONE_INTEGER, &
                                         ghosted_neighbors_size, &
                                         ghosted_neighbors, &
                                         option)   

  disp_vec = 0.d0
  ux_weighted = 0.d0
  uy_weighted = 0.d0
  uz_weighted = 0.d0
  
  
  do i = 1, ghosted_neighbors_size
    disp_vec(1) = grid%x(ghosted_neighbors(i)) - grid%x(ghosted_id)
    disp_vec(2) = grid%y(ghosted_neighbors(i)) - grid%y(ghosted_id)
    disp_vec(3) = grid%z(ghosted_neighbors(i)) - grid%z(ghosted_id)
    ux_weighted = ux_weighted + disp_vec* &
                    (global_aux_vars(ghosted_neighbors(i))%displacement(1) - &
                     global_aux_vars(ghosted_id)%displacement(1))
    uy_weighted = uy_weighted + disp_vec* &
                    (global_aux_vars(ghosted_neighbors(i))%displacement(2) - &
                     global_aux_vars(ghosted_id)%displacement(2))
    uz_weighted = uz_weighted + disp_vec* &
                    (global_aux_vars(ghosted_neighbors(i))%displacement(3) - &
                     global_aux_vars(ghosted_id)%displacement(3))

  enddo
  
  gradient = Minv*reshape((/ux_weighted(1),ux_weighted(2),ux_weighted(3), &
                           uy_weighted(1),uy_weighted(2),uy_weighted(3), &
                           uz_weighted(1),uz_weighted(2),uz_weighted(3)/), &
                           (/3,3/))  
  
end subroutine THMCComputeDisplacementGradient

! ************************************************************************** !
! 
! THMCComputeMinv: Computes the inverse of the 
! Author: Satish Karra
! Date: 2/20/12
!
! ************************************************************************** !


subroutine THMCComputeMinv(grid, ghosted_id, M_inv, option) 


  use Grid_module
  use Global_Aux_module
  use Option_module
  use Utility_module

  implicit none
#include "finclude/petscdmda.h"

  type(option_type) :: option
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  

  
  PetscInt :: ghosted_neighbors_size, ghosted_id
  PetscInt :: ghosted_neighbors(26)
  PetscReal :: M_inv(3,3), disp_vec(3,1), disp_mat(3,3)
  PetscReal :: identity(3,3)

  PetscInt :: i
  
  PetscInt :: INDX(3)
  PetscInt :: D
   
  call GridGetGhostedNeighborsWithCorners(grid,ghosted_id, &
                                         DMDA_STENCIL_STAR, &
                                         ONE_INTEGER,ONE_INTEGER,ONE_INTEGER, &
                                         ghosted_neighbors_size, &
                                         ghosted_neighbors, &
                                         option)   

  disp_vec = 0.d0
  disp_mat = 0.d0
  identity = reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))
  
  do i = 1, ghosted_neighbors_size
    disp_vec(1,1) = grid%x(ghosted_neighbors(i)) - grid%x(ghosted_id)
    disp_vec(2,1) = grid%y(ghosted_neighbors(i)) - grid%y(ghosted_id)
    disp_vec(3,1) = grid%z(ghosted_neighbors(i)) - grid%z(ghosted_id)
    disp_mat = disp_mat + matmul(disp_vec,transpose(disp_vec))
  enddo


  call ludcmp(disp_mat,THREE_INTEGER,INDX,D)
  call lubksb(disp_mat,THREE_INTEGER,INDX,identity)

  M_inv = identity
      
end subroutine THMCComputeMinv


! ************************************************************************** !
! 
! THMCComputeDisplacementGradientPert: Computes the perturbation of the 
! gradient of displacement
! Author: Satish Karra
! Date: 4/13/12
!
! ************************************************************************** !


subroutine THMCComputeDisplacementGradientPert(grid, global_aux_vars, &
                                               ghosted_id, &
                                               gradient_pert, direction_flag, &
                                               perturbation_tolerance, option) 


  use Grid_module
  use Global_Aux_module
  use Option_module
  use Utility_module

  implicit none
#include "finclude/petscdmda.h"

  type(option_type) :: option
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_aux_vars(:)

  
  PetscInt :: ghosted_neighbors_size, ghosted_id
  PetscInt :: ghosted_neighbors(26)
  PetscReal :: gradient_pert(3,3), disp_vec(3,1), disp_mat(3,3)
  PetscReal :: ux_weighted(3,1)
  PetscReal :: uy_weighted(3,1)
  PetscReal :: uz_weighted(3,1)
  PetscReal :: perturbation_tolerance

  PetscInt :: i, direction_flag, flag_x, flag_y, flag_z
  PetscInt :: INDX(3)
  PetscInt :: D
   
  call GridGetGhostedNeighborsWithCorners(grid,ghosted_id, &
                                         DMDA_STENCIL_STAR, &
                                         ONE_INTEGER,ONE_INTEGER,ONE_INTEGER, &
                                         ghosted_neighbors_size, &
                                         ghosted_neighbors, &
                                         option)   
 
  disp_vec = 0.d0
  disp_mat = 0.d0
  ux_weighted = 0.d0
  uy_weighted = 0.d0
  uz_weighted = 0.d0
  
  select case(direction_flag)
    case (ONE_INTEGER)
      flag_x = 1
      flag_y = 0
      flag_z = 0
    case (TWO_INTEGER)
      flag_x = 0
      flag_y = 1
      flag_z = 0
    case (THREE_INTEGER)
      flag_x = 0
      flag_y = 0
      flag_z = 1
   end select
 
 
  do i = 1, ghosted_neighbors_size
    disp_vec(1,1) = grid%x(ghosted_neighbors(i)) - grid%x(ghosted_id)
    disp_vec(2,1) = grid%y(ghosted_neighbors(i)) - grid%y(ghosted_id)
    disp_vec(3,1) = grid%z(ghosted_neighbors(i)) - grid%z(ghosted_id)
    disp_mat = disp_mat + matmul(disp_vec,transpose(disp_vec))
    ux_weighted = ux_weighted + disp_vec* &
                    (global_aux_vars(ghosted_neighbors(i))%displacement(1) - &
                     (global_aux_vars(ghosted_id)%displacement(1) &
                      + flag_x*perturbation_tolerance))
    uy_weighted = uy_weighted + disp_vec* &
                    (global_aux_vars(ghosted_neighbors(i))%displacement(2) - &
                     (global_aux_vars(ghosted_id)%displacement(2) &
                      + flag_y*perturbation_tolerance))
    uz_weighted = uz_weighted + disp_vec* &
                    (global_aux_vars(ghosted_neighbors(i))%displacement(3) - &
                     (global_aux_vars(ghosted_id)%displacement(3) &
                      + flag_z*perturbation_tolerance))
  enddo

  call ludcmp(disp_mat,THREE_INTEGER,INDX,D)
  call lubksb(disp_mat,THREE_INTEGER,INDX,ux_weighted)
  call lubksb(disp_mat,THREE_INTEGER,INDX,uy_weighted)
  call lubksb(disp_mat,THREE_INTEGER,INDX,uz_weighted)

  gradient_pert(:,1) = ux_weighted(:,1)
  gradient_pert(:,2) = uy_weighted(:,1)
  gradient_pert(:,3) = uz_weighted(:,1)
  
  
end subroutine THMCComputeDisplacementGradientPert

! ************************************************************************** !
!
! THMCLiteGetTecplotHeader: Returns THMC contribution to 
!                               Tecplot file header
! author:
! date: 3/2/12
!
! ************************************************************************** !
function THMCGetTecplotHeader(realization,icolumn)

  use Realization_class
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: THMCGetTecplotHeader
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
    write(string2,'('',"'',i2,''-disp_x"'')') icolumn
  else
    write(string2,'('',"disp_x"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-disp_y"'')') icolumn
  else
    write(string2,'('',"disp_y"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-disp_z"'')') icolumn
  else
    write(string2,'('',"disp_z"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sl"'')') icolumn
  else
    write(string2,'('',"Sl"'')')
  endif
  string = trim(string) // trim(string2)

#ifdef ICE
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sg"'')') icolumn
  else
    write(string2,'('',"Sg"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Si"'')') icolumn
  else
    write(string2,'('',"Si"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-deni"'')') icolumn
  else
    write(string2,'('',"deni"'')')
  endif
  string = trim(string) // trim(string2)
#endif

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-denl"'')') icolumn
  else
    write(string2,'('',"denl"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Ul"'')') icolumn
  else
    write(string2,'('',"Ul"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-visl"'')') icolumn
  else
    write(string2,'('',"visl"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-mobilityl"'')') icolumn
  else
    write(string2,'('',"mobilityl"'')')
  endif
  string = trim(string) // trim(string2)

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i2,'')"'')') icolumn,i
    else
      write(string2,'('',"Xl('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo
  
  THMCGetTecplotHeader = string

end function THMCGetTecplotHeader

! ************************************************************************** !
!
! THMCDestroy: Deallocates variables associated with Richard
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCDestroy(patch)

  use Patch_module

  implicit none
  
  type(patch_type) :: patch
  
  ! need to free array in aux vars
  call THMCAuxDestroy(patch%aux%THMC)

end subroutine THMCDestroy

end module THMC_module
