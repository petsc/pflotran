module THC_module
#include "finclude/petscmat.h"
  use petscmat
  use THC_Aux_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-8

  public THCResidual,THCJacobian, &
         THCUpdateFixedAccumulation,THCTimeCut,&
         THCSetup, THCNumericalJacobianTest, &
         THCMaxChange, THCUpdateSolution, &
         THCGetTecplotHeader, THCInitializeTimestep, &
         THCComputeMassBalance, THCResidualToMass, &
         THCSecondaryHeat, THCSecondaryHeatJacobian, & 
         THCUpdateAuxVars, THCDestroy, &
         THCCheckUpdatePost, THCCheckUpdatePre
         
  PetscInt, parameter :: jh2o = 1

contains

! ************************************************************************** !

subroutine THCTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
 
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
 
  call VecCopy(field%flow_yy,field%flow_xx,ierr);CHKERRQ(ierr)
  call THCInitializeTimestep(realization)
 
end subroutine THCTimeCut

! ************************************************************************** !

subroutine THCSetup(realization)
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THCSetupPatch(realization)
    cur_patch => cur_patch%next
  enddo

  call THCSetPlotVariables(realization)

end subroutine THCSetup

! ************************************************************************** !

subroutine THCSetupPatch(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(thc_auxvar_type), pointer :: thc_auxvars(:), thc_auxvars_bc(:)
  type(thc_auxvar_type), pointer :: thc_auxvars_ss(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  type(sec_heat_type), pointer :: thc_sec_heat_vars(:)
  type(coupler_type), pointer :: initial_condition
  PetscReal :: area_per_vol

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, iphase, local_id
  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
    
  patch%aux%THC => THCAuxCreate(option)
  patch%aux%SC_heat => SecondaryAuxHeatCreate(option)

! option%io_buffer = 'Before THC can be run, the thc_parameter object ' // &
!                    'must be initialized with the proper variables ' // &
!                    'THCAuxCreate() is called anywhere.'
! call printErrMsg(option)
  allocate(patch%aux%THC%thc_parameter%sir(option%nphase, &
                                  size(realization%saturation_function_array)))
  
  !Jitu, 08/04/2010: Check these allocations. Currently assumes only single value in the array	<modified pcl 1-13-11>
  allocate(patch%aux%THC%thc_parameter%dencpr(size(realization%material_property_array)))
  allocate(patch%aux%THC%thc_parameter%ckwet(size(realization%material_property_array)))
  allocate(patch%aux%THC%thc_parameter%ckdry(size(realization%material_property_array)))
  allocate(patch%aux%THC%thc_parameter%alpha(size(realization%material_property_array)))
  if (option%use_th_freezing) then
     allocate(patch%aux%THC%thc_parameter%ckfrozen(size(realization%material_property_array)))
     allocate(patch%aux%THC%thc_parameter%alpha_fr(size(realization%material_property_array)))
  endif

  !Copy the values in the thc_parameter from the global realization 
  do i = 1, size(realization%material_property_array)
    patch%aux%THC%thc_parameter%dencpr(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%rock_density*option%scale* &
        realization%material_property_array(i)%ptr%specific_heat
 
    patch%aux%THC%thc_parameter%ckwet(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%thermal_conductivity_wet*option%scale  
    patch%aux%THC%thc_parameter%ckdry(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%thermal_conductivity_dry*option%scale
    patch%aux%THC%thc_parameter%alpha(realization%material_property_array(i)%ptr%id) = &
      realization%material_property_array(i)%ptr%alpha
    if (option%use_th_freezing) then
       patch%aux%THC%thc_parameter%ckfrozen(realization%material_property_array(i)%ptr%id) = &
            realization%material_property_array(i)%ptr%thermal_conductivity_frozen*option%scale
       patch%aux%THC%thc_parameter%alpha_fr(realization%material_property_array(i)%ptr%id) = &
            realization%material_property_array(i)%ptr%alpha_fr
    endif

  enddo 

  do i = 1, size(realization%saturation_function_array)
    patch%aux%THC%thc_parameter%sir(:,realization%saturation_function_array(i)%ptr%id) = &
      realization%saturation_function_array(i)%ptr%Sr(:)
  enddo

  ! allocate auxvar data structures for all grid cells
  allocate(thc_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call THCAuxVarInit(thc_auxvars(ghosted_id),option)
    ! currently, hardwire to first fluid
    thc_auxvars(ghosted_id)%diff(1:option%nflowspec) = &
      realization%fluid_properties%diffusion_coefficient
  enddo

 
  if (option%use_mc) then
    initial_condition => patch%initial_conditions%first
    allocate(thc_sec_heat_vars(grid%nlmax))
  
    do local_id = 1, grid%nlmax
  
    ! Assuming the same secondary continuum for all regions (need to make it an array)
    ! S. Karra 07/18/12
      call SecondaryContinuumSetProperties( &
        thc_sec_heat_vars(local_id)%sec_continuum, &
        realization%material_property_array(1)%ptr%secondary_continuum_name, &
        realization%material_property_array(1)%ptr%secondary_continuum_length, &
        realization%material_property_array(1)%ptr%secondary_continuum_matrix_block_size, &
        realization%material_property_array(1)%ptr%secondary_continuum_fracture_spacing, &
        realization%material_property_array(1)%ptr%secondary_continuum_radius, &
        realization%material_property_array(1)%ptr%secondary_continuum_area, &
        option)
        
      thc_sec_heat_vars(local_id)%ncells = &
        realization%material_property_array(1)%ptr%secondary_continuum_ncells
      thc_sec_heat_vars(local_id)%aperture = &
        realization%material_property_array(1)%ptr%secondary_continuum_aperture
      thc_sec_heat_vars(local_id)%epsilon = &
        realization%material_property_array(1)%ptr%secondary_continuum_epsilon
      thc_sec_heat_vars(local_id)%log_spacing = &
        realization%material_property_array(1)%ptr%secondary_continuum_log_spacing
      thc_sec_heat_vars(local_id)%outer_spacing = &
        realization%material_property_array(1)%ptr%secondary_continuum_outer_spacing
                
      allocate(thc_sec_heat_vars(local_id)%area(thc_sec_heat_vars(local_id)%ncells))
      allocate(thc_sec_heat_vars(local_id)%vol(thc_sec_heat_vars(local_id)%ncells))
      allocate(thc_sec_heat_vars(local_id)%dm_minus(thc_sec_heat_vars(local_id)%ncells))
      allocate(thc_sec_heat_vars(local_id)%dm_plus(thc_sec_heat_vars(local_id)%ncells))
      allocate(thc_sec_heat_vars(local_id)%sec_continuum% &
             distance(thc_sec_heat_vars(local_id)%ncells))
                 
      call SecondaryContinuumType(thc_sec_heat_vars(local_id)%sec_continuum, &
                                  thc_sec_heat_vars(local_id)%ncells, &
                                  thc_sec_heat_vars(local_id)%area, &
                                  thc_sec_heat_vars(local_id)%vol, &
                                  thc_sec_heat_vars(local_id)%dm_minus, &
                                  thc_sec_heat_vars(local_id)%dm_plus, &
                                  thc_sec_heat_vars(local_id)%aperture, &
                                  thc_sec_heat_vars(local_id)%epsilon, &
                                  thc_sec_heat_vars(local_id)%log_spacing, &
                                  thc_sec_heat_vars(local_id)%outer_spacing, &
                                  area_per_vol,option)
                                
      thc_sec_heat_vars(local_id)%interfacial_area = area_per_vol* &
        (1.d0 - thc_sec_heat_vars(local_id)%epsilon)* &
        realization%material_property_array(1)%ptr% &
        secondary_continuum_area_scaling

    ! Setting the initial values of all secondary node temperatures same as primary node 
    ! temperatures (with initial dirichlet BC only) -- sk 06/26/12
      allocate(thc_sec_heat_vars(local_id)%sec_temp(thc_sec_heat_vars(local_id)%ncells))
      
      if (option%set_secondary_init_temp) then
        thc_sec_heat_vars(local_id)%sec_temp = &
          realization%material_property_array(1)%ptr%secondary_continuum_init_temp
      else
        thc_sec_heat_vars(local_id)%sec_temp = &
        initial_condition%flow_condition%temperature%dataset%rarray(1)
      endif
        
    enddo
      
    patch%aux%SC_heat%sec_heat_vars => thc_sec_heat_vars    

  endif

  
  patch%aux%THC%auxvars => thc_auxvars
  patch%aux%THC%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  boundary_condition => patch%boundary_conditions%first

  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

!  write(*,*)'Sum_connection', sum_connection
  if (sum_connection > 0) then 
    allocate(thc_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call THCAuxVarInit(thc_auxvars_bc(iconn),option)
      ! currently, hardwire to first fluid
      thc_auxvars_bc(iconn)%diff(1:option%nflowspec) = &
        realization%fluid_properties%diffusion_coefficient
    enddo
    patch%aux%THC%auxvars_bc => thc_auxvars_bc
  endif
  patch%aux%THC%num_aux_bc = sum_connection

  ! Create aux vars for source/sink
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    allocate(thc_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call THCAuxVarInit(thc_auxvars_ss(iconn),option)
      ! currently, hardwire to first fluid
      thc_auxvars_ss(iconn)%diff(1:option%nflowspec) = &
        realization%fluid_properties%diffusion_coefficient
    enddo
    patch%aux%THC%auxvars_ss => thc_auxvars_ss
  endif
  patch%aux%THC%num_aux_ss = sum_connection


  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call THCCreateZeroArray(patch,option)
  
  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    patch%aux%THC%thc_parameter%diffusion_coefficient(iphase) = &
      cur_fluid_property%diffusion_coefficient
    patch%aux%THC%thc_parameter%diffusion_activation_energy(iphase) = &
      cur_fluid_property%diffusion_activation_energy
    cur_fluid_property => cur_fluid_property%next
  enddo

end subroutine THCSetupPatch

! ************************************************************************** !

subroutine THCCheckUpdatePre(line_search,P,dP,changed,realization,ierr)
  ! 
  ! Checks update prior to update
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 08/02/12
  ! 

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  ! ignore changed flag for now.
  PetscBool :: changed
  type(realization_type) :: realization
  
  PetscReal, pointer :: P_p(:)
  PetscReal, pointer :: dP_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(thc_auxvar_type), pointer :: thc_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: P0, P1, P_R, delP, delP_old
  PetscReal :: scale, press_limit, temp_limit
  PetscInt :: iend, istart
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  thc_auxvars => realization%patch%aux%THC%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars

  if (dabs(option%pressure_change_limit) > 0.d0) then

    patch => realization%patch

    call VecGetArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(P,P_p,ierr);CHKERRQ(ierr)

    press_limit = dabs(option%pressure_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      P0 = P_p(istart)
      delP = dP_p(istart)
      if (press_limit < dabs(delP)) then
        write(option%io_buffer,'("dP_trunc:",1i7,2es15.7)') &         
          grid%nG2A(grid%nL2G(local_id)),press_limit,dabs(delP)
        call printMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),press_limit),delP)
      dP_p(istart) = delP
    enddo
    
    call VecRestoreArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(P,P_p,ierr);CHKERRQ(ierr)

  endif
  
  if (dabs(option%temperature_change_limit) > 0.d0) then
      
    patch => realization%patch

    call VecGetArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(P,P_p,ierr);CHKERRQ(ierr)

    temp_limit = dabs(option%temperature_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      P0 = P_p(iend)
      delP = dP_p(iend)
      if (abs(delP) > abs(temp_limit)) then
        write(option%io_buffer,'("dT_trunc:",1i7,2es15.7)') &
          grid%nG2A(grid%nL2G(local_id)),temp_limit,dabs(delP)
        call printMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),temp_limit),delP)
      dP_p(iend) = delP
    enddo
    
    call VecRestoreArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(P,P_p,ierr);CHKERRQ(ierr)
    
  endif


  if (dabs(option%pressure_dampening_factor) > 0.d0) then
    ! P^p+1 = P^p - dP^p
    P_R = option%reference_pressure
    scale = option%pressure_dampening_factor

    call VecGetArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(P,P_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      delP = dP_p(istart)
      P0 = P_p(istart)
      P1 = P0 - delP
      if (P0 < P_R .and. P1 > P_R) then
        write(option%io_buffer,'("U -> S:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1 
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(istart)
        call printMsgAnyRank(option)
#endif
      else if (P1 < P_R .and. P0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(istart)
        call printMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (P0 < P_R .and. P1 > P_R) then
        dP_p(istart) = scale*delP
      endif
    enddo
    call VecRestoreArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(P,P_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THCCheckUpdatePre

! ************************************************************************** !

subroutine THCCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                                   P1_changed,realization,ierr)
  ! 
  ! Checks update after each update
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/25/12
  ! 

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
  use Secondary_Continuum_Aux_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  type(realization_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  
  PetscReal, pointer :: P1_p(:)
  PetscReal, pointer :: dP_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(thc_auxvar_type), pointer :: thc_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  type(thc_parameter_type), pointer :: thc_parameter
  type(sec_heat_type), pointer :: thc_sec_heat_vars(:)


  PetscInt :: local_id, ghosted_id
  PetscReal :: Res(3)
  PetscReal :: inf_norm, global_inf_norm
  PetscErrorCode :: ierr
  PetscReal :: vol_frac_prim
  PetscInt :: istart, iend
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  thc_auxvars => realization%patch%aux%THC%auxvars
  thc_parameter => realization%patch%aux%THC%thc_parameter
  global_auxvars => realization%patch%aux%Global%auxvars
  thc_sec_heat_vars => realization%patch%aux%SC_heat%sec_heat_vars

  
  dP_changed = PETSC_FALSE
  P1_changed = PETSC_FALSE
  
  if (option%flow%check_post_convergence) then
    call VecGetArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(P1,P1_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%volume,volume_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    
    inf_norm = 0.d0
    vol_frac_prim = 1.d0
    
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
    
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      
      if (option%use_mc) then
        vol_frac_prim = thc_sec_heat_vars(local_id)%epsilon
      endif

      call THCAccumulation(thc_auxvars(ghosted_id), &
                           global_auxvars(ghosted_id), &
                           porosity_loc_p(ghosted_id), &
                           volume_p(local_id), &
                           thc_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                           option,vol_frac_prim,Res)
                                                        
      inf_norm = max(inf_norm,min(dabs(dP_p(istart)/P1_p(istart)), &
                                  dabs(r_p(istart)/Res(1)), &
                                  dabs(dP_p(iend)/P1_p(iend)), &
                                  dabs(r_p(iend)/Res(3))))
    enddo
    call MPI_Allreduce(inf_norm,global_inf_norm,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    option%converged = PETSC_TRUE
    if (global_inf_norm > option%flow%post_convergence_tol) &
      option%converged = PETSC_FALSE
    call VecRestoreArrayF90(dP,dP_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(P1,P1_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%volume,volume_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p, &
                            ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif
  
end subroutine THCCheckUpdatePost

! ************************************************************************** !

subroutine THCComputeMassBalance(realization, mass_balance)
  ! 
  ! THComputeMassBalance:
  ! Adapted from RichardsComputeMassBalance: need to be checked
  ! 
  ! Author: Jitendra Kumar
  ! Date: 07/21/2010
  ! 

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)
   
  type(patch_type), pointer :: cur_patch

  mass_balance = 0.d0

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THCComputeMassBalancePatch(realization, mass_balance)
    cur_patch => cur_patch%next
  enddo

end subroutine THCComputeMassBalance    

! ************************************************************************** !

subroutine THCComputeMassBalancePatch(realization,mass_balance)
  ! 
  ! THComputeMassBalancePatch:
  ! Adapted from RichardsComputeMassBalancePatch: need to be checked
  ! 
  ! Author: Jitendra Kumar
  ! Date: 07/21/2010
  ! 
 
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
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: volume_p(:), porosity_loc_p(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_auxvars => patch%aux%Global%auxvars

  call VecGetArrayF90(field%volume,volume_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    ! mass = volume*saturation*density
    mass_balance = mass_balance + &
      global_auxvars(ghosted_id)%den_kg* &
      global_auxvars(ghosted_id)%sat* &
      porosity_loc_p(ghosted_id)*volume_p(local_id)
  enddo

  call VecRestoreArrayF90(field%volume,volume_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p, &
                          ierr);CHKERRQ(ierr)
  
end subroutine THCComputeMassBalancePatch

! ************************************************************************** !

subroutine THCZeroMassBalDeltaPatch(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11
  ! 
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%THC%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%THC%num_aux_bc > 0) then
    do iconn = 1, patch%aux%THC%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
 
end subroutine THCZeroMassBalDeltaPatch

! ************************************************************************** !

subroutine THCUpdateMassBalancePatch(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: ???
  ! Date: 12/13/11
  ! 
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%THC%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance = &
      patch%aux%Global%auxvars(iconn)%mass_balance + &
      patch%aux%Global%auxvars(iconn)%mass_balance_delta*FMWH2O* &
      option%flow_dt
  enddo
#endif

  if (patch%aux%THC%num_aux_bc > 0) then
    do iconn = 1, patch%aux%THC%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance = &
        global_auxvars_bc(iconn)%mass_balance + &
        global_auxvars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif


end subroutine THCUpdateMassBalancePatch

! ************************************************************************** !

subroutine THCUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the THC problem
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THCUpdateAuxVarsPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THCUpdateAuxVars

! ************************************************************************** !

subroutine THCUpdateAuxVarsPatch(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the THC problem
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

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
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(thc_auxvar_type), pointer :: thc_auxvars(:)
  type(thc_auxvar_type), pointer :: thc_auxvars_bc(:)
  type(thc_auxvar_type), pointer :: thc_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscReal, pointer :: xx(:)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr
  
  !!
!  PetscReal, allocatable :: gradient(:,:)
  !!
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  !!
!  allocate(gradient(grid%ngmax,3))
!  gradient = 0.d0
  !!
  
  thc_auxvars => patch%aux%THC%auxvars
  thc_auxvars_bc => patch%aux%THC%auxvars_bc
  thc_auxvars_ss => patch%aux%THC%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))

    if (option%use_th_freezing) then
       call THCAuxVarComputeIce(xx_loc_p(istart:iend), &
            thc_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            iphase, &
            realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
            porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
            option)
    else
       call THCAuxVarCompute(xx_loc_p(istart:iend), &
            thc_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            iphase, &
            realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
            porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
            option)
    endif

!    call THCComputeGradient(grid, global_auxvars, ghosted_id, &
!                       gradient(ghosted_id,:), option) 


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

      do idof=1,option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
      
      select case(boundary_condition%flow_condition%itype(THC_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
          iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      if (option%use_th_freezing) then
         call THCAuxVarComputeIce(xxbc,thc_auxvars_bc(sum_connection), &
              global_auxvars_bc(sum_connection), &
              iphasebc, &
              realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
              porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
              option)
      else
         call THCAuxVarCompute(xxbc,thc_auxvars_bc(sum_connection), &
              global_auxvars_bc(sum_connection), &
              iphasebc, &
              realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
              porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
              option)
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! source/sinks
  source_sink => patch%source_sinks%first
  sum_connection = 0
  allocate(xx(option%nflowdof))
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      iend = ghosted_id*option%nflowdof
      istart = iend-option%nflowdof+1
      iphase = int(iphase_loc_p(ghosted_id))

      tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
      xx = xx_loc_p(istart:iend)
      xx(2) = tsrc1

      if (option%use_th_freezing) then
         call THCAuxVarComputeIce(xx, &
              thc_auxvars_ss(sum_connection),global_auxvars_ss(sum_connection), &
              iphase, &
              realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
              porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
              option)
      else
         call THCAuxVarCompute(xx, &
              thc_auxvars_ss(sum_connection),global_auxvars_ss(sum_connection), &
              iphase, &
              realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
              porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
              option)
      endif
    enddo
    source_sink => source_sink%next
  enddo
  deallocate(xx)

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p, &
                          ierr);CHKERRQ(ierr)
  
  patch%aux%THC%auxvars_up_to_date = PETSC_TRUE
    
!  deallocate(gradient)

end subroutine THCUpdateAuxVarsPatch

! ************************************************************************** !

subroutine THCInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: ???
  ! Date: 02/20/08
  ! 

  use Realization_class
  
  implicit none
  
  type(realization_type) :: realization

  call THCUpdateFixedAccumulation(realization)

end subroutine THCInitializeTimestep

! ************************************************************************** !

subroutine THCUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time step
  ! 
  ! Author: ???
  ! Date: 02/13/08
  ! 

  use Realization_class
  use Field_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization

  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  field => realization%field
    
  call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THCUpdateSolutionPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THCUpdateSolution

! ************************************************************************** !

subroutine THCUpdateSolutionPatch(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11, 02/28/14
  ! 


  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
    
  implicit none
  
  type(realization_type) :: realization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(THC_parameter_type), pointer :: THC_parameter
  type(THC_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(sec_heat_type), pointer :: THC_sec_heat_vars(:)

  PetscInt :: istart, iend  
  PetscInt :: local_id, ghosted_id
  ! secondary continuum variables
  PetscReal :: sec_dencpr
  PetscErrorCode :: ierr
  PetscReal, pointer :: ithrm_loc_p(:)

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  THC_parameter => patch%aux%THC%THC_parameter
  auxvars => patch%aux%THC%auxvars
  global_auxvars => patch%aux%Global%auxvars
  
  if (option%use_mc) then
    THC_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  endif

  if (realization%option%compute_mass_balance_new) then
    call THCUpdateMassBalancePatch(realization)
  endif

  if (option%use_mc) then
    call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
    
      sec_dencpr = THC_parameter%dencpr(int(ithrm_loc_p(local_id))) ! secondary rho*c_p same as primary for now
        
      call THCSecHeatAuxVarCompute(THC_sec_heat_vars(local_id), &
                            global_auxvars(ghosted_id), &
                            THC_parameter%ckwet(int(ithrm_loc_p(local_id))), &
                            sec_dencpr, &
                            option)
                            
    enddo   
    call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  endif


end subroutine THCUpdateSolutionPatch

! ************************************************************************** !

subroutine THCUpdateFixedAccumulation(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THCUpdateFixedAccumPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THCUpdateFixedAccumulation

! ************************************************************************** !

subroutine THCUpdateFixedAccumPatch(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Secondary_Continuum_Aux_module


  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(thc_auxvar_type), pointer :: thc_auxvars(:)
  type(thc_parameter_type), pointer :: thc_parameter
  type(sec_heat_type), pointer :: thc_sec_heat_vars(:)

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:), accum_p(:), perm_xx_loc_p(:)
  PetscReal :: vol_frac_prim
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  thc_parameter => patch%aux%THC%thc_parameter
  thc_auxvars => patch%aux%THC%auxvars
  global_auxvars => patch%aux%Global%auxvars
  thc_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

          
  call VecGetArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tortuosity_loc,tor_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%volume,volume_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)


  vol_frac_prim = 1.d0
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))

    if (option%use_th_freezing) then
       call THCAuxVarComputeIce(xx_p(istart:iend), &
            thc_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            iphase, &
            realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
            porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
            option)
    else
       call THCAuxVarCompute(xx_p(istart:iend), &
            thc_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            iphase, &
            realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
            porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
            option)
    endif


    if (option%use_mc) then
      vol_frac_prim = thc_sec_heat_vars(local_id)%epsilon
    endif
    
    iphase_loc_p(ghosted_id) = iphase
    call THCAccumulation(thc_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              thc_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,vol_frac_prim,accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tortuosity_loc,tor_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%volume,volume_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

#if 0
   call THCNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine THCUpdateFixedAccumPatch

! ************************************************************************** !

subroutine THCNumericalJacobianTest(xx,realization)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

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
  
  call VecDuplicate(xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res_pert,ierr);CHKERRQ(ierr)
  
  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof,grid%nlmax*option%nflowdof, &
                   ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
    
  call THCResidual(PETSC_NULL_SNES,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do icell = 1,grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    endif
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call VecCopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      call THCResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization,ierr)
      call VecGetArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values, &
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(option%mycomm,'numerical_jacobian.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  call MatDestroy(A,ierr);CHKERRQ(ierr)
  
  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)
  
end subroutine THCNumericalJacobianTest

! ************************************************************************** !

subroutine THCAccumDerivative(thc_auxvar,global_auxvar,por,vol, &
                              rock_dencpr,option,sat_func, &
                              vol_frac_prim,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module
  use Saturation_Function_module
  
  use EOS_Water_module
  
  implicit none

  type(thc_auxvar_type) :: thc_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: vol,por,rock_dencpr
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec 
  PetscReal :: porXvol, mol(option%nflowspec), eng, por1

  PetscInt :: iphase, ideriv
  type(thc_auxvar_type) :: thc_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  PetscReal :: x(3), x_pert(3), pert, res(3), res_pert(3), J_pert(3,3)
  PetscReal :: vol_frac_prim, tempreal
  
  ! variables required for ice
  PetscReal :: sat_g, p_g, den_g, p_sat, mol_g, u_g, C_g
  PetscReal :: dpsat_dt, ddeng_dt, dmolg_dt, dsatg_dp, dsatg_dt, dug_dt
  PetscReal :: sat_i, den_i, u_i
  PetscReal :: dsati_dp, dsati_dt
  PetscReal :: ddeni_dp, ddeni_dt
  PetscReal :: dui_dt
  PetscReal, parameter :: C_a = 1.86d-3 ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! in MJ/kg/K
  PetscErrorCode :: ierr

  
  ! X = {p, T, x_2}; R = {R_x1, R_x2, R_T} = {R_p, R_x2, R_T}
  
  
#ifndef USE_COMPRESSIBILITY
  porXvol = por*vol
  J(1,1) = (global_auxvar%sat(1)*thc_auxvar%dden_dp + &
           thc_auxvar%dsat_dp*global_auxvar%den(1))*porXvol !*thc_auxvar%xmol(1)
#else
  if (thc_auxvar%pc > 0.d0) then
    por1 = por
  else
    por1 = 1.d0-(1.d0-por)*exp(-1.d-10*(abs(global_auxvar%pres(1)- &
                                       option%reference_pressure)))
  endif
  
  porXvol = por1*vol
  
  if (thc_auxvar%pc > 0.d0) then
    J(1,1) = (global_auxvar%sat(1)*thc_auxvar%dden_dp + &
             thc_auxvar%dsat_dp*global_auxvar%den(1))*porXvol
  else
    tempreal = exp(-1.d-10*(abs(global_auxvar%pres(1)-option%reference_pressure)))
    J(1,1) = (global_auxvar%sat(1)*thc_auxvar%dden_dp + &
             thc_auxvar%dsat_dp*global_auxvar%den(1))*porXvol + &
             global_auxvar%sat(1)*global_auxvar%den(1)*vol*1.d-10* &
             (1.d0 - por)*tempreal*abs(global_auxvar%pres(1)- &
             option%reference_pressure)/(global_auxvar%pres(1)- &
             option%reference_pressure)
  endif
#endif

  J(1,2) = global_auxvar%sat(1)*thc_auxvar%dden_dt*porXvol !*thc_auxvar%xmol(1)
  J(1,3) = 0.d0 !-global_auxvar%sat(1)*global_auxvar%den(1)*porXvol
  J(2,1) = (global_auxvar%sat(1)*thc_auxvar%dden_dp + &
           thc_auxvar%dsat_dp*global_auxvar%den(1))*porXvol*thc_auxvar%xmol(2)
  J(2,2) = global_auxvar%sat(1)*thc_auxvar%dden_dt*porXvol*thc_auxvar%xmol(2)
  J(2,3) = global_auxvar%sat(1)*global_auxvar%den(1)*porXvol
  J(3,1) = (thc_auxvar%dsat_dp*global_auxvar%den(1)*thc_auxvar%u + &
            global_auxvar%sat(1)*thc_auxvar%dden_dp*thc_auxvar%u + &
            global_auxvar%sat(1)*global_auxvar%den(1)*thc_auxvar%du_dp)*porXvol
  J(3,2) = global_auxvar%sat(1)* &
           (thc_auxvar%dden_dt*thc_auxvar%u + &  ! pull %sat outside
            global_auxvar%den(1)*thc_auxvar%du_dt)*porXvol +  &
           (1.d0 - por)*vol*rock_dencpr 
  J(3,3) = 0.d0

  if (option%use_th_freezing) then
     ! SK, 11/17/11
     sat_g = thc_auxvar%sat_gas
     sat_i = thc_auxvar%sat_ice
     u_i = thc_auxvar%u_ice
     den_i = thc_auxvar%den_ice
     p_g = option%reference_pressure ! set to reference pressure
     den_g = p_g/(IDEAL_GAS_CONST*(global_auxvar%temp + 273.15d0))*1.d-3
     call EOSWaterSaturationPressure(global_auxvar%temp, p_sat, dpsat_dt, ierr)
     mol_g = p_sat/p_g
     C_g = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR !in MJ/kmol/K, expression might be different
     u_g = C_g*(global_auxvar%temp + 273.15d0)
     ddeng_dt = - p_g/(IDEAL_GAS_CONST*(global_auxvar%temp + 273.15d0)**2)*1.d-3
     dmolg_dt = dpsat_dt/p_g
     dsatg_dp = thc_auxvar%dsat_gas_dp
     dsatg_dt = thc_auxvar%dsat_gas_dt
     dug_dt = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)* &
          (global_auxvar%temp + 273.15d0)
     dsati_dt = thc_auxvar%dsat_ice_dt
     dsati_dp = thc_auxvar%dsat_ice_dp
     ddeni_dt = thc_auxvar%dden_ice_dt
     ddeni_dp = thc_auxvar%dden_ice_dp
     dui_dt = thc_auxvar%du_ice_dt
 
     J(1,1) = J(1,1) + (dsatg_dp*den_g*mol_g + dsati_dp*den_i + &
          sat_i*ddeni_dp)*porXvol
     J(1,2) = J(1,2) + (thc_auxvar%dsat_dt*global_auxvar%den(1) + &
          dsatg_dt*den_g*mol_g + sat_g*ddeng_dt*mol_g + &
          sat_g*den_g*dmolg_dt + dsati_dt*den_i + sat_i*ddeni_dt)* &
          porXvol
     J(2,2) = J(2,2) + thc_auxvar%dsat_dt*global_auxvar%den(1)* &
          thc_auxvar%xmol(2)*porXvol
     J(3,1) = J(3,1) + (dsatg_dp*den_g*u_g + dsati_dp*den_i*u_i + &
          sat_i*ddeni_dp*u_i)*porXvol
     J(3,2) = J(3,2) + (thc_auxvar%dsat_dt*global_auxvar%den(1)*thc_auxvar%u + &
          dsatg_dt*den_g*u_g + sat_g*ddeng_dt*u_g + &
          sat_g*den_g*dug_dt + dsati_dt*den_i*u_i + &
          sat_i*ddeni_dt*u_i + sat_i*den_i*dui_dt)*porXvol
                    
#ifdef REMOVE_SATURATION
     J(2,1) = thc_auxvar%dden_dp*porXvol*thc_auxvar%xmol(2)
     J(2,2) = thc_auxvar%dden_dt*porXvol*thc_auxvar%xmol(2)            
     J(2,3) = global_auxvar%den(1)*porXvol
#endif                    
  endif ! if (use_th_freezing)

  J = J/option%flow_dt
  
  J(option%nflowdof,:) = vol_frac_prim*J(option%nflowdof,:)

  if (option%numerical_derivatives_flow) then
    allocate(thc_auxvar_pert%xmol(option%nflowspec),thc_auxvar_pert%diff(option%nflowspec))
    call GlobalAuxVarInit(global_auxvar_pert,option)  
    call THCAuxVarCopy(thc_auxvar,thc_auxvar_pert,option)
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)

    x(1) = global_auxvar%pres(1)
    x(2) = global_auxvar%temp
    x(3) = thc_auxvar%xmol(2)
    

    call THCAccumulation(thc_auxvar,global_auxvar, &
                         por,vol,rock_dencpr,option, &
                         vol_frac_prim,res)
    
    do ideriv = 1,3
      pert = x(ideriv)*perturbation_tolerance
      x_pert = x
      
      if (option%use_th_freezing) then
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
      else
         x_pert(ideriv) = x_pert(ideriv) + pert
      endif


      if (option%use_th_freezing) then
         call THCAuxVarComputeIce(x_pert,thc_auxvar_pert,global_auxvar_pert,iphase,sat_func, &
              0.d0,0.d0,option)
      else
         call THCAuxVarCompute(x_pert,thc_auxvar_pert,global_auxvar_pert,iphase,sat_func, &
              0.d0,0.d0,option)
      endif

#if 0      
      select case(ideriv)
        case(1)
          print *, 'dsat_dp:', thc_auxvar%dsat_dp, (global_auxvar_pert%sat(1)-global_auxvar%sat(1))/pert
          print *, 'dden_dp:', thc_auxvar%dden_dp, (global_auxvar_pert%den(1)-global_auxvar%den(1))/pert
          print *, 'dkvr_dp:', thc_auxvar%dkvr_dp, (thc_auxvar_pert%kvr-thc_auxvar%kvr)/pert
          print *, 'dh_dp:', thc_auxvar%dh_dp, (thc_auxvar_pert%h-thc_auxvar%h)/pert
          print *, 'du_dp:', thc_auxvar%du_dp, (thc_auxvar_pert%u-thc_auxvar%u)/pert
          if (option%use_th_freezing) then
             print *, 'dsati_dp:', thc_auxvar%dsat_ice_dp, (thc_auxvar_pert%sat_ice - thc_auxvar%sat_ice)/pert
             print *, 'dsatg_dp:', thc_auxvar%dsat_gas_dp, (thc_auxvar_pert%sat_gas - thc_auxvar%sat_gas)/pert
             print *, 'ddeni_dp:', thc_auxvar%dden_ice_dp, (thc_auxvar_pert%den_ice - thc_auxvar%den_ice)/pert
          endif
          
        case(2)
          print *, 'dden_dt:', thc_auxvar%dden_dt, (global_auxvar_pert%den(1)-global_auxvar%den(1))/pert
          print *, 'dkvr_dt:', thc_auxvar%dkvr_dt, (thc_auxvar_pert%kvr-thc_auxvar%kvr)/pert
          print *, 'dh_dt:', thc_auxvar%dh_dt, (thc_auxvar_pert%h-thc_auxvar%h)/pert
          print *, 'du_dt:', thc_auxvar%du_dt, (thc_auxvar_pert%u-thc_auxvar%u)/pert
          if (option%use_th_freezing) then
             print *, 'ddeni_dt:', thc_auxvar%dden_ice_dt, (thc_auxvar_pert%den_ice - thc_auxvar%den_ice)/pert
             print *, 'dsati_dt:', thc_auxvar%dsat_ice_dt, (thc_auxvar_pert%sat_ice - thc_auxvar%sat_ice)/pert
             print *, 'dsatg_dt:', thc_auxvar%dsat_gas_dt, (thc_auxvar_pert%sat_gas - thc_auxvar%sat_gas)/pert
             print *, 'dsat_dt:', thc_auxvar%dsat_dt, (global_auxvar_pert%sat(1) - global_auxvar%sat(1))/pert
             print *, 'dui_dt:', thc_auxvar%du_ice_dt, (thc_auxvar_pert%u_ice - thc_auxvar%u_ice)/pert
          endif
      end select     
#endif     
      call THCAccumulation(thc_auxvar_pert,global_auxvar_pert, &
                           por,vol,rock_dencpr,option,vol_frac_prim, &
                           res_pert)
      J_pert(:,ideriv) = (res_pert(:)-res(:))/pert
    enddo

    deallocate(thc_auxvar_pert%xmol,thc_auxvar_pert%diff)
    J = J_pert
    call GlobalAuxVarStrip(global_auxvar_pert)  
  endif
   
end subroutine THCAccumDerivative

! ************************************************************************** !

subroutine THCAccumulation(auxvar,global_auxvar,por,vol, &
                           rock_dencpr,option,vol_frac_prim,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module
  
  use EOS_Water_module
  
  implicit none

  type(thc_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: vol,por,rock_dencpr,por1
     
  PetscInt :: ispec 
  PetscReal :: porXvol, mol(option%nflowspec), eng
  PetscReal :: vol_frac_prim

  ! variables for ice
  PetscReal :: sat_g, p_g, den_g, p_sat, mol_g, u_g, C_g
  PetscReal :: sat_i, den_i, u_i
  PetscReal, parameter :: C_a = 1.86d-3 ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! in MJ/kg/K
  PetscErrorCode :: ierr
  
! TechNotes, THC Mode: First term of Equation 8
  porXvol = por*vol
  
#ifndef USE_COMPRESSIBILITY  
  mol(1) = global_auxvar%sat(1)*global_auxvar%den(1)*porXvol
#else
  if (auxvar%pc > 0.d0) then
    por1 = por 
  else
    por1 = 1.d0-(1.d0-por)*exp(-1.d-10*(abs(global_auxvar%pres(1)- &
                                       option%reference_pressure)))
  endif
  mol(1) = global_auxvar%sat(1)*global_auxvar%den(1)*por1*vol
#endif
    
  mol(2) = global_auxvar%sat(1)*global_auxvar%den(1)*auxvar%xmol(2)*porXvol

! TechNotes, THC Mode: First term of Equation 9
  eng = global_auxvar%sat(1) * &
        global_auxvar%den(1) * &
        auxvar%u * porXvol + &
        (1.d0 - por) * vol * rock_dencpr * global_auxvar%temp

  if (option%use_th_freezing) then
     ! SK, 11/17/11
     sat_g = auxvar%sat_gas
     sat_i = auxvar%sat_ice
     u_i = auxvar%u_ice
     den_i = auxvar%den_ice
     p_g = option%reference_pressure
     den_g = p_g/(IDEAL_GAS_CONST*(global_auxvar%temp + 273.15d0))*1.d-3 !in kmol/m3
     call EOSWaterSaturationPressure(global_auxvar%temp, p_sat, ierr)
     mol_g = p_sat/p_g
     C_g = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR ! in MJ/kmol/K
     u_g = C_g*(global_auxvar%temp + 273.15d0)       ! in MJ/kmol
     mol(1) = mol(1) + (sat_g*den_g*mol_g + sat_i*den_i)*porXvol
     eng = eng + (sat_g*den_g*u_g + sat_i*den_i*u_i)*porXvol
#ifdef REMOVE_SATURATION
     mol(2) = global_auxvar%den(1)*auxvar%xmol(2)*porXvol
#endif
  endif

  Res(1:option%nflowdof-1) = mol(:)/option%flow_dt
  Res(option%nflowdof) = vol_frac_prim*eng/option%flow_dt

end subroutine THCAccumulation

! ************************************************************************** !

subroutine THCFluxDerivative(auxvar_up,global_auxvar_up,por_up,tor_up, &
                             sir_up,dd_up,perm_up,Dk_up, &
                             auxvar_dn,global_auxvar_dn,por_dn,tor_dn, &
                             sir_dn,dd_dn,perm_dn,Dk_dn, &
                             area,dist_gravity,upweight, &
                             option,sat_func_up,sat_func_dn, &
                             Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
                             Dk_ice_up,Dk_ice_dn, &
                             alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                             Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
                             
  use Option_module 
  use Saturation_Function_module             
  
  use EOS_Water_module
  
  implicit none
  
  type(thc_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
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
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)
     
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
  type(thc_auxvar_type) :: auxvar_pert_up, auxvar_pert_dn
  type(global_auxvar_type) :: global_auxvar_pert_up, global_auxvar_pert_dn
  PetscReal :: x_up(3), x_dn(3), x_pert_up(3), x_pert_dn(3), pert_up, pert_dn, &
            res(3), res_pert_up(3), res_pert_dn(3), J_pert_up(3,3), J_pert_dn(3,3)

  ! ice variables
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
  PetscReal, parameter :: R_gas_constant = 8.3144621 ! Gas constant in J/mol/K
  PetscReal :: fv_up, fv_dn
  PetscReal :: dfv_dt_up, dfv_dt_dn
  PetscReal :: dfv_dp_up, dfv_dp_dn
  PetscReal :: dmolg_dp_up, dmolg_dp_dn
  PetscReal :: fv_up_pert
  
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
  
  if (option%use_th_freezing) then
     dfv_dt_up = 0.d0
     dfv_dt_dn = 0.d0
     dfv_dp_up = 0.d0
     dfv_dp_dn = 0.d0
     dmolg_dp_up = 0.d0
     dmolg_dp_dn = 0.d0
     dmolg_dt_up = 0.d0
     dmolg_dt_dn = 0.d0
  endif
  
! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or. global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
    dden_ave_dp_up = upweight*auxvar_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
    dden_ave_dt_up = upweight*auxvar_up%dden_dt
    dden_ave_dt_dn = (1.D0-upweight)*auxvar_dn%dden_dt

    gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
              (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
              * dist_gravity
    dgravity_dden_up = upweight*auxvar_up%avgmw*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

    dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*auxvar_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
    dphi_dt_up = dgravity_dden_up*auxvar_up%dden_dt
    dphi_dt_dn = dgravity_dden_dn*auxvar_dn%dden_dt

! note uxmol only contains one phase xmol
    if (dphi>=0.D0) then
      ukvr = auxvar_up%kvr
      dukvr_dp_up = auxvar_up%dkvr_dp
      dukvr_dt_up = auxvar_up%dkvr_dt
      
      uh = auxvar_up%h
      duh_dp_up = auxvar_up%dh_dp
      duh_dt_up = auxvar_up%dh_dt
      
      uxmol(1:option%nflowspec) = auxvar_up%xmol(1:option%nflowspec)
      duxmol_dxmol_up = 1.d0
    else
      ukvr = auxvar_dn%kvr
      dukvr_dp_dn = auxvar_dn%dkvr_dp
      dukvr_dt_dn = auxvar_dn%dkvr_dt
      
      uh = auxvar_dn%h
      duh_dp_dn = auxvar_dn%dh_dp
      duh_dt_dn = auxvar_dn%dh_dt
      
      uxmol(1:option%nflowspec) = auxvar_dn%xmol(1:option%nflowspec)
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
      Jup(option%nflowdof,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uh+q*density_ave*duh_dp_up
      Jup(option%nflowdof,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uh+q*density_ave*duh_dt_up

      Jdn(option%nflowdof,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn
      Jdn(option%nflowdof,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn

    endif
  endif 

    difff = diffdp*0.25D0*(global_auxvar_up%sat(1) + global_auxvar_dn%sat(1))* &
                            (global_auxvar_up%den(1) + global_auxvar_dn%den(1))
    ddifff_dp_up = diffdp*0.25D0*(auxvar_up%dsat_dp*(global_auxvar_up%den(1) + global_auxvar_dn%den(1))+ &
                  (global_auxvar_up%sat(1) + global_auxvar_dn%sat(1))*auxvar_up%dden_dp)
    ddifff_dt_up = diffdp*0.25D0*(global_auxvar_up%sat(1) + global_auxvar_dn%sat(1))*auxvar_up%dden_dt

    ddifff_dp_dn = diffdp*0.25D0*(auxvar_dn%dsat_dp*(global_auxvar_up%den(1) + global_auxvar_dn%den(1))+ &
                  (global_auxvar_up%sat(1) + global_auxvar_dn%sat(1))*auxvar_dn%dden_dp)
    ddifff_dt_dn = diffdp*0.25D0*(global_auxvar_up%sat(1) + global_auxvar_dn%sat(1))*auxvar_dn%dden_dt
                                    
! SK   

    Jup(2,1) = Jup(2,1) + ddifff_dp_up*0.5d0*(Diff_up + Diff_dn)* &
                         (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
    Jup(2,2) = Jup(2,2) + ddifff_dt_up*0.5d0*(Diff_up + Diff_dn)* &
                         (auxvar_up%xmol(2) - auxvar_dn%xmol(2)) 
    Jup(2,3) = Jup(2,3) + difff*0.5d0*(Diff_up + Diff_dn)

    Jdn(2,1) = Jdn(2,1) + ddifff_dp_dn*0.5d0*(Diff_up + Diff_dn)* &
                         (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
    Jdn(2,2) = Jdn(2,2) + ddifff_dt_dn*0.5d0*(Diff_up + Diff_dn)* &
                         (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
    Jdn(2,3) = Jdn(2,3) + difff*0.5d0*(Diff_up + Diff_dn)*(-1.d0)
    
 

    if (option%use_th_freezing) then
       ! Added by Satish Karra, LANL, updated 11/11/11
       satg_up = auxvar_up%sat_gas
       satg_dn = auxvar_dn%sat_gas
       if ((satg_up > eps) .and. (satg_dn > eps)) then
          p_g = option%reference_pressure  ! set to reference pressure
          deng_up = p_g/(IDEAL_GAS_CONST*(global_auxvar_up%temp + 273.15d0))*1.d-3
          deng_dn = p_g/(IDEAL_GAS_CONST*(global_auxvar_dn%temp + 273.15d0))*1.d-3

          Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
          p_ref = 1.01325d5   ! in Pa
          T_ref = 25.d0       ! in deg C
    
          Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + 273.15d0)/ &
               (T_ref + 273.15d0))**(1.8)  
          Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + 273.15d0)/ &
               (T_ref + 273.15d0))**(1.8)    

             
          Ddiffgas_up = por_up*tor_up*satg_up*deng_up*Diffg_up
          Ddiffgas_dn = por_dn*tor_dn*satg_dn*deng_dn*Diffg_dn
          call EOSWaterSaturationPressure(global_auxvar_up%temp, psat_up, dpsat_dt_up, ierr)
          call EOSWaterSaturationPressure(global_auxvar_dn%temp, psat_dn, dpsat_dt_dn, ierr)
  
          ! vapor pressure lowering due to capillary pressure
          fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
               R_gas_constant*(global_auxvar_up%temp + 273.15d0)))
          fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
               R_gas_constant*(global_auxvar_dn%temp + 273.15d0)))
  
          molg_up = psat_up*fv_up/p_g
          molg_dn = psat_dn*fv_dn/p_g
  
          dfv_dt_up = fv_up*(auxvar_up%pc/R_gas_constant/(global_auxvar_up%den(1)* &
               (global_auxvar_up%temp + 273.15d0))**2)* &
               (auxvar_up%dden_dt*(global_auxvar_up%temp + 273.15d0) &
               + global_auxvar_up%den(1))
          dfv_dt_dn = fv_dn*(auxvar_dn%pc/R_gas_constant/(global_auxvar_dn%den(1)* &
               (global_auxvar_dn%temp + 273.15d0))**2)* &
               (auxvar_dn%dden_dt*(global_auxvar_dn%temp + 273.15d0) &
               + global_auxvar_dn%den(1))
  
          dfv_dp_up = fv_up*(auxvar_up%pc/R_gas_constant/(global_auxvar_up%den(1))**2/ &
               (global_auxvar_up%temp + 273.15d0)*auxvar_up%dden_dp &
               + 1.d0/R_gas_constant/global_auxvar_up%den(1)/ &
               (global_auxvar_up%temp + 273.15d0))
          dfv_dp_dn = fv_dn*(auxvar_dn%pc/R_gas_constant/(global_auxvar_dn%den(1))**2/ &
               (global_auxvar_dn%temp + 273.15d0)*auxvar_dn%dden_dp &
               + 1.d0/R_gas_constant/global_auxvar_dn%den(1)/ &
               (global_auxvar_dn%temp + 273.15d0))

          dmolg_dt_up = (1/p_g)*dpsat_dt_up*fv_up + psat_up/p_g*dfv_dt_up
          dmolg_dt_dn = (1/p_g)*dpsat_dt_dn*fv_dn + psat_dn/p_g*dfv_dt_dn
  
          dmolg_dp_up = psat_up/p_g*dfv_dp_up
          dmolg_dp_dn = psat_dn/p_g*dfv_dp_dn
  
          ddeng_dt_up = - p_g/(IDEAL_GAS_CONST*(global_auxvar_up%temp + &
               273.15d0)**2)*1.d-3
          ddeng_dt_dn = - p_g/(IDEAL_GAS_CONST*(global_auxvar_dn%temp + &
               273.15d0)**2)*1.d-3
  
          dDiffg_dt_up = 1.8*Diffg_up/(global_auxvar_up%temp + 273.15d0)
          dDiffg_dt_dn = 1.8*Diffg_dn/(global_auxvar_dn%temp + 273.15d0)
  
          dDiffg_dp_up = 0.d0
          dDiffg_dp_dn = 0.d0

          dsatg_dp_up = auxvar_up%dsat_gas_dp
          dsatg_dp_dn = auxvar_dn%dsat_gas_dp
     
          if (molg_up > molg_dn) then 
             upweight = 0.d0
          else 
             upweight = 1.d0
          endif

          Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn 

#ifndef NO_VAPOR_DIFFUION  
          Jup(1,1) = Jup(1,1) + (upweight*por_up*tor_up*deng_up*(Diffg_up*dsatg_dp_up &
               + satg_up*dDiffg_dp_up)* &
               (molg_up - molg_dn) + Ddiffgas_up*dmolg_dp_up)/ &
               (dd_up + dd_dn)*area 

          Jup(1,2) = Jup(1,2) + (upweight*por_up*tor_up*satg_up*(Diffg_up* &
               ddeng_dt_up + deng_up*dDiffg_dt_up)*(molg_up - molg_dn) &
               + Ddiffgas_avg*dmolg_dt_up)/(dd_up + dd_dn)*area
  
          Jdn(1,1) = Jdn(1,1) + ((1.D0 - upweight)*por_dn*tor_dn*deng_dn* &
               (Diffg_dn*dsatg_dp_dn + satg_dn*dDiffg_dp_dn)* &
               (molg_up - molg_dn) + Ddiffgas_avg*(-dmolg_dp_dn))/ &
               (dd_up + dd_dn)*area
              
          Jdn(1,2) = Jdn(1,2) + ((1.D0 - upweight)*por_dn*tor_dn*satg_dn*(Diffg_dn* &
               ddeng_dt_dn + deng_dn*dDiffg_dp_dn)*(molg_up - molg_dn) &
               + Ddiffgas_avg*(-dmolg_dt_dn))/(dd_up + dd_dn)*area
#endif
       endif


       ddifff_dt_up = diffdp*0.25d0*(auxvar_up%dsat_dt)*&
            (global_auxvar_up%den(1) + global_auxvar_dn%den(1))
       ddifff_dt_dn = diffdp*0.25d0*(auxvar_dn%dsat_dt)*&
            (global_auxvar_up%den(1) + global_auxvar_dn%den(1))

       Jup(2,2) = Jup(2,2) + ddifff_dt_up*0.5d0*(Diff_up + Diff_dn)* &
            (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
       Jdn(2,2) = Jdn(2,2) + ddifff_dt_dn*0.5d0*(Diff_up + Diff_dn)* &
            (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
                         
#ifdef REMOVE_SATURATION 

       difff = diffdp * 0.5D0* &
            (global_auxvar_up%den(1)+global_auxvar_dn%den(1))
       ddifff_dp_up = diffdp * 0.5d0 * auxvar_up%dden_dp
       ddifff_dp_dn = diffdp * 0.5d0 * auxvar_dn%dden_dp
       ddifff_dt_up = diffdp * 0.5d0 * auxvar_up%dden_dt
       ddifff_dt_dn = diffdp * 0.5d0 * auxvar_dn%dden_dt
  
       Jup(2,1) = Jup(2,1) + ddifff_dp_up*0.5d0*(Diff_up + Diff_dn)* &
            (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
       Jup(2,2) = Jup(2,2) + ddifff_dt_up*0.5d0*(Diff_up + Diff_dn)* &
            (auxvar_up%xmol(2) - auxvar_dn%xmol(2)) 
       Jup(2,3) = Jup(2,3) + difff*0.5d0*(Diff_up + Diff_dn)

       Jdn(2,1) = Jdn(2,1) + ddifff_dp_dn*0.5d0*(Diff_up + Diff_dn)* &
            (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
       Jdn(2,2) = Jdn(2,2) + ddifff_dt_dn*0.5d0*(Diff_up + Diff_dn)* &
            (auxvar_up%xmol(2) - auxvar_dn%xmol(2))
       Jdn(2,3) = Jdn(2,3) + difff*0.5d0*(Diff_up + Diff_dn)*(-1.d0)

#endif
    endif ! if (use_th_freezing)

        
! conduction term  
  Ke_up = (global_auxvar_up%sat(1) + epsilon)**(alpha_up)   !unfrozen soil Kersten number
  Ke_dn = (global_auxvar_dn%sat(1) + epsilon)**(alpha_dn)
  
  dKe_dp_up = alpha_up*(global_auxvar_up%sat(1) + epsilon)**(alpha_up - 1.d0)* &
               auxvar_up%dsat_dp
  dKe_dp_dn = alpha_dn*(global_auxvar_dn%sat(1) + epsilon)**(alpha_dn - 1.d0)* &
               auxvar_dn%dsat_dp
  
  if (option%use_th_freezing) then
     Ke_fr_up = (auxvar_up%sat_ice + epsilon)**(alpha_fr_up)
     Ke_fr_dn = (auxvar_dn%sat_ice + epsilon)**(alpha_fr_dn)

     Dk_eff_up = Dk_up*Ke_up + Dk_ice_up*Ke_fr_up + &
          (1.d0 - Ke_up - Ke_fr_up)*Dk_dry_up
     Dk_eff_dn = Dk_dn*Ke_dn + Dk_ice_dn*Ke_fr_dn + &
          (1.d0 - Ke_dn - Ke_fr_dn)*Dk_dry_dn

     dKe_dt_up = alpha_up*(global_auxvar_up%sat(1) + epsilon)**(alpha_up - 1.d0)* &
          auxvar_up%dsat_dt
     dKe_dt_dn = alpha_dn*(global_auxvar_dn%sat(1) + epsilon)**(alpha_dn - 1.d0)* &
          auxvar_dn%dsat_dt
               
     dKe_fr_dt_up = alpha_fr_up*(global_auxvar_up%sat(1) + epsilon)** &
          (alpha_fr_up - 1.d0)*auxvar_up%dsat_dt
     dKe_fr_dt_dn = alpha_fr_dn*(global_auxvar_dn%sat(1) + epsilon)** &
          (alpha_fr_dn - 1.d0)*auxvar_dn%dsat_dt
                 
     dKe_fr_dp_up = alpha_fr_up*(global_auxvar_up%sat(1) + epsilon)** &
          (alpha_fr_up - 1.d0)*auxvar_up%dsat_dp
     dKe_fr_dp_dn = alpha_fr_dn*(global_auxvar_dn%sat(1) + epsilon)** &
          (alpha_fr_dn - 1.d0)*auxvar_dn%dsat_dp           
  else 
     Dk_eff_up = Dk_dry_up + (Dk_up - Dk_dry_up)*Ke_up
     Dk_eff_dn = Dk_dry_dn + (Dk_dn - Dk_dry_dn)*Ke_dn 
  
     dKe_dt_up = 0.d0
     dKe_dt_dn = 0.d0
  endif
 
  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)
  
  if (option%use_th_freezing) then
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
  else
     dDk_dt_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dt_up
     dDk_dt_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dt_dn
  
     dDk_dp_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dp_up
     dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dp_dn
  endif
    
  !  cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)
  Jup(option%nflowdof,1) = Jup(option%nflowdof,1) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dp_up
  Jdn(option%nflowdof,1) = Jdn(option%nflowdof,1) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dp_dn
                           
  Jup(option%nflowdof,2) = Jup(option%nflowdof,2) + Dk*area + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dt_up
  Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2) + Dk*area*(-1.d0) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dt_dn

 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives_flow) then
    allocate(auxvar_pert_up%xmol(option%nflowspec),auxvar_pert_up%diff(option%nflowspec))
    allocate(auxvar_pert_dn%xmol(option%nflowspec),auxvar_pert_dn%diff(option%nflowspec))
    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    call THCAuxVarCopy(auxvar_up,auxvar_pert_up,option)
    call THCAuxVarCopy(auxvar_dn,auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
    x_up(1) = global_auxvar_up%pres(1)
    x_up(2) = global_auxvar_up%temp
    x_up(3) = auxvar_up%xmol(2)
    x_dn(1) = global_auxvar_dn%pres(1)
    x_dn(2) = global_auxvar_dn%temp
    x_dn(3) = auxvar_dn%xmol(2)
    

    call THCFlux( &
      auxvar_up,global_auxvar_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
      auxvar_dn,global_auxvar_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
      area,dist_gravity,upweight, &
      option,v_darcy,Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
      Dk_ice_up,Dk_ice_dn, &
      alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
      res)
    do ideriv = 1,3
      pert_up = x_up(ideriv)*perturbation_tolerance
      pert_dn = x_dn(ideriv)*perturbation_tolerance
      x_pert_up = x_up
      x_pert_dn = x_dn

      if (option%use_th_freezing) then      
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
      else
         x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
         x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif

      if (option%use_th_freezing) then
         call THCAuxVarComputeIce(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up, &
              iphase,sat_func_up, &
              0.d0,0.d0,option)
         call THCAuxVarComputeIce(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn,iphase,sat_func_dn, &
              0.d0,0.d0,option)
      else
         call THCAuxVarCompute(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up, &
              iphase,sat_func_up, &
              0.d0,0.d0,option)
         call THCAuxVarCompute(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn,iphase,sat_func_dn, &
              0.d0,0.d0,option)
      endif

      call THCFlux(auxvar_pert_up,global_auxvar_pert_up, &
                   por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                   auxvar_dn,global_auxvar_dn, &
                   por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                   area,dist_gravity,upweight, &
                   option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   res_pert_up)
      call THCFlux(auxvar_up,global_auxvar_up, &
                   por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                   auxvar_pert_dn,global_auxvar_pert_dn, &
                   por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                   area,dist_gravity,upweight, &
                   option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   res_pert_dn)
                                             
      J_pert_up(:,ideriv) = (res_pert_up(:)-res(:))/pert_up
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    
    deallocate(auxvar_pert_up%xmol,auxvar_pert_up%diff)
    deallocate(auxvar_pert_dn%xmol,auxvar_pert_dn%diff)
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)    
  endif

end subroutine THCFluxDerivative

! ************************************************************************** !

subroutine THCFlux(auxvar_up,global_auxvar_up, &
                  por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                  auxvar_dn,global_auxvar_dn, &
                  por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                  area,dist_gravity,upweight, &
                  option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                  Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
                  
  use Option_module                              
  
  use EOS_Water_module

  implicit none
  
  type(thc_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
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

  ! ice variables
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
  PetscReal :: fv_up, fv_dn
  PetscReal, parameter :: R_gas_constant = 8.3144621 ! Gas constant in J/mol/K
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up*tor_up*por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0  
  
! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or. global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) < eps) then 
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) < eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1) 

    gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
              (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
              * dist_gravity

    dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity

!   note uxmol only contains one component xmol
    if (dphi >= 0.D0) then
      ukvr = auxvar_up%kvr
      uh = auxvar_up%h
      uxmol(1:option%nflowspec) = auxvar_up%xmol(1:option%nflowspec)
    else
      ukvr = auxvar_dn%kvr
      uh = auxvar_dn%h
      uxmol(1:option%nflowspec) = auxvar_dn%xmol(1:option%nflowspec)
    endif

    if (ukvr > floweps) then
      v_darcy = Dq * ukvr * dphi
   
      q = v_darcy * area
        
      fluxm(1) = fluxm(1) + q*density_ave
      fluxm(2) = fluxm(2) + q*density_ave*uxmol(2)
      fluxe = fluxe + q*density_ave*uh 
    endif
  endif 

  
  difff = diffdp * 0.25D0*(global_auxvar_up%sat(1)+global_auxvar_dn%sat(1))* &
                            (global_auxvar_up%den(1)+global_auxvar_dn%den(1))
  fluxm(2) = fluxm(2) + difff * .5D0 * (Diff_up + Diff_dn)* &
                 (auxvar_up%xmol(2) - auxvar_dn%xmol(2)) 


  if (option%use_th_freezing) then
     ! Added by Satish Karra, 10/24/11
     satg_up = auxvar_up%sat_gas
     satg_dn = auxvar_dn%sat_gas
     if ((satg_up > eps) .and. (satg_dn > eps)) then
        p_g = option%reference_pressure ! set to reference pressure
        deng_up = p_g/(IDEAL_GAS_CONST*(global_auxvar_up%temp + 273.15d0))*1.d-3
        deng_dn = p_g/(IDEAL_GAS_CONST*(global_auxvar_dn%temp + 273.15d0))*1.d-3
    
        Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
        p_ref = 1.01325d5 ! in Pa
        T_ref = 25.d0 ! in deg C
    
        Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + 273.15d0)/ &
             (T_ref + 273.15d0))**(1.8)  
        Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + 273.15d0)/ &
             (T_ref + 273.15d0))**(1.8)
             
        Ddiffgas_up = por_up*tor_up*satg_up*deng_up*Diffg_up
        Ddiffgas_dn = por_dn*tor_dn*satg_dn*deng_dn*Diffg_dn
        call EOSWaterSaturationPressure(global_auxvar_up%temp, psat_up, ierr)
        call EOSWaterSaturationPressure(global_auxvar_dn%temp, psat_dn, ierr)
  
        ! vapor pressure lowering due to capillary pressure
        fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
             R_gas_constant*(global_auxvar_up%temp + 273.15d0)))
        fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
             R_gas_constant*(global_auxvar_dn%temp + 273.15d0)))
  
        molg_up = psat_up*fv_up/p_g
        molg_dn = psat_dn*fv_dn/p_g
          
        if (molg_up > molg_dn) then 
           upweight = 0.d0
        else 
           upweight = 1.d0
        endif
    
        Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn 
#ifndef NO_VAPOR_DIFFUSION
        fluxm(1) = fluxm(1) + Ddiffgas_avg*area*(molg_up - molg_dn)/ &
             (dd_up + dd_dn)
#endif

     endif
#ifdef REMOVE_SATURATION 
     difff = diffdp * 0.5D0* &
          (global_auxvar_up%den(1)+global_auxvar_dn%den(1))
     fluxm(2) = fluxm(2) + difff * .5D0 * (Diff_up + Diff_dn)* &
          (auxvar_up%xmol(2) - auxvar_dn%xmol(2)) 
#endif
  endif ! if (use_th_freezing)

! conduction term  
  Ke_up = (global_auxvar_up%sat(1) + epsilon)**(alpha_up)   !unfrozen soil Kersten number
  Ke_dn = (global_auxvar_dn%sat(1) + epsilon)**(alpha_dn)
     
  if (option%use_th_freezing) then

     Ke_fr_up = (auxvar_up%sat_ice + epsilon)**(alpha_fr_up)
     Ke_fr_dn = (auxvar_dn%sat_ice + epsilon)**(alpha_fr_dn)

     Dk_eff_up = Dk_up*Ke_up + Dk_ice_up*Ke_fr_up + &
          (1.d0 - Ke_up - Ke_fr_up)*Dk_dry_up
     Dk_eff_dn = Dk_dn*Ke_dn + Dk_ice_dn*Ke_fr_dn + &
          (1.d0 - Ke_dn - Ke_fr_dn)*Dk_dry_dn
  else

     Dk_eff_up = Dk_dry_up + (Dk_up - Dk_dry_up)*Ke_up
     Dk_eff_dn = Dk_dry_dn + (Dk_dn - Dk_dry_dn)*Ke_dn      

  endif
 
  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)
  cond = Dk*area*(global_auxvar_up%temp - global_auxvar_dn%temp)
  fluxe = fluxe + cond

  Res(1:option%nflowdof-1) = fluxm(:)
  Res(option%nflowdof) = fluxe
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  
 
  
end subroutine THCFlux

! ************************************************************************** !

subroutine THCBCFluxDerivative(ibndtype,auxvars, &
                              auxvar_up,global_auxvar_up, &
                              auxvar_dn,global_auxvar_dn, &
                              por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                              area,dist_gravity,option, &
                              sat_func_dn,Diff_dn,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux
  ! terms for the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
  use Option_module
  use Saturation_Function_module
  
  use EOS_Water_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(thc_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn
  PetscReal :: auxvars(:) ! from aux_real_var array in boundary condition
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
  type(thc_auxvar_type) :: auxvar_pert_dn, auxvar_pert_up
  type(global_auxvar_type) :: global_auxvar_pert_dn, global_auxvar_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(3), x_up(3), x_pert_dn(3), x_pert_up(3), pert_dn, res(3), &
               res_pert_dn(3), J_pert_dn(3,3)

  ! ice variables
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
  select case(ibndtype(THC_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dq = perm_dn / dd_up
      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
        dden_ave_dt_dn = (1.D0-upweight)*auxvar_dn%dden_dt

        if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dt_dn = dden_ave_dt_dn + upweight*auxvar_up%dden_dt
        endif
        
        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

        dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
        dphi_dt_dn = dgravity_dden_dn*auxvar_dn%dden_dt

        if (ibndtype(THC_PRESSURE_DOF) == SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_auxvar_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
            dphi_dp_dn = 0.d0
            dphi_dt_dn = 0.d0
          endif
        endif        
        
        if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                                   !( dgravity_dden_up                   ) (dden_dt_up)
          dphi_dt_dn = dphi_dt_dn + upweight*auxvar_up%avgmw*dist_gravity*auxvar_up%dden_dt
        endif
        
        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
          if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dt_dn = auxvar_up%dkvr_dt
          endif
        else
          ukvr = auxvar_dn%kvr
          dukvr_dp_dn = auxvar_dn%dkvr_dp
          dukvr_dt_dn = auxvar_dn%dkvr_dt
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(auxvars(THC_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(THC_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
          if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dden_ave_dt_dn = auxvar_up%dden_dt
          endif
        else 
          density_ave = global_auxvar_dn%den(1)
          dden_ave_dp_dn = auxvar_dn%dden_dp
          dden_ave_dt_dn = auxvar_dn%dden_dt
        endif 
        q = v_darcy * area
      endif

  end select

  if (v_darcy >= 0.D0) then
    uh = auxvar_up%h
    uxmol(:)=auxvar_up%xmol(1:option%nflowspec)
    if (ibndtype(THC_PRESSURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dp_dn = auxvar_up%dh_dp
    endif
    if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dt_dn = auxvar_up%dh_dt
    endif
    if (ibndtype(THC_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      duxmol_dxmol_dn = 1.d0
    endif
  else
    uh = auxvar_dn%h
    duh_dp_dn = auxvar_dn%dh_dp
    duh_dt_dn = auxvar_dn%dh_dt

    uxmol(:)=auxvar_dn%xmol(1:option%nflowspec)
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
  Jdn(option%nflowdof,1) =  &
     ((dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn)
  Jdn(option%nflowdof,2) =  &
     ((dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn)
!  Jdn(option%nflowdof,3:option%nflowdof) = 0.d0

  ! Diffusion term   
  select case(ibndtype(THC_CONCENTRATION_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC) 

      diff = diffdp * global_auxvar_dn%sat(1)*global_auxvar_dn%den(1)
      ddiff_dp_dn = diffdp * (auxvar_dn%dsat_dp*global_auxvar_dn%den(1)+ &
                                global_auxvar_dn%sat(1)*auxvar_dn%dden_dp)
      ddiff_dt_dn = diffdp * global_auxvar_dn%sat(1)*auxvar_dn%dden_dt
! SK
      Jdn(2,1) = Jdn(2,1) + ddiff_dp_dn*Diff_dn* &
                            (auxvar_up%xmol(2)-auxvar_dn%xmol(2))
      Jdn(2,2) = Jdn(2,2) + ddiff_dt_dn*Diff_dn* &
                            (auxvar_up%xmol(2)-auxvar_dn%xmol(2))
      Jdn(2,3) = Jdn(2,3) + diff*Diff_dn*(-1.d0)

      if (option%use_th_freezing) then
#ifdef REMOVE_SATURATION
         diff = diffdp * global_auxvar_dn%den(1)
      
         ddiff_dp_dn = diffdp * auxvar_dn%dden_dp
         ddiff_dt_dn = diffdp * auxvar_dn%dden_dt

         Jdn(2,1) = Jdn(2,1) + ddiff_dp_dn*Diff_dn* &
              (auxvar_up%xmol(2)-auxvar_dn%xmol(2))
         Jdn(2,2) = Jdn(2,2) + ddiff_dt_dn*Diff_dn* &
              (auxvar_up%xmol(2)-auxvar_dn%xmol(2))
         Jdn(2,3) = Jdn(2,3) + diff*Diff_dn*(-1.d0)
      
#endif
      endif
  end select
   
  ! Conduction term
  select case(ibndtype(THC_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dk =  Dk_dn / dd_up
      !cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)
      Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2)+Dk*area*(-1.d0)
      if (option%use_th_freezing) then
         ! Added by Satish Karra, 11/21/11
         satg_up = auxvar_up%sat_gas
         satg_dn = auxvar_dn%sat_gas
         if ((satg_up > eps) .and. (satg_dn > eps)) then
            p_g = option%reference_pressure  ! set to reference pressure
            deng_up = p_g/(IDEAL_GAS_CONST*(global_auxvar_up%temp + &
                 273.15d0))*1.d-3
            deng_dn = p_g/(IDEAL_GAS_CONST*(global_auxvar_dn%temp + &
                 273.15d0))*1.d-3
        
            Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
            p_ref = 1.01325d5 ! in Pa
            T_ref = 25.d0 ! in deg C 
        
            Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)  
            Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)  
            Ddiffgas_up = satg_up*deng_up*Diffg_up
            Ddiffgas_dn = satg_dn*deng_dn*Diffg_dn
            call EOSWaterSaturationPressure(global_auxvar_up%temp, psat_up, ierr)
            call EOSWaterSaturationPressure(global_auxvar_dn%temp, psat_dn, dpsat_dt_dn, ierr)
            molg_up = psat_up/p_g
            molg_dn = psat_dn/p_g
            ddeng_dt_dn = - p_g/(IDEAL_GAS_CONST*(global_auxvar_dn%temp + &
                 273.15d0)**2)*1.d-3
            dmolg_dt_dn = (1/p_g)*dpsat_dt_dn
            dDiffg_dt_dn = 1.8*Diffg_dn/(global_auxvar_dn%temp + 273.15d0)
            dDiffg_dp_dn = 0.d0
            dsatg_dp_dn = auxvar_dn%dsat_gas_dp
        
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
      endif

  end select

  if (option%numerical_derivatives_flow) then
    allocate(auxvar_pert_dn%xmol(option%nflowspec),auxvar_pert_dn%diff(option%nflowspec))
    allocate(auxvar_pert_up%xmol(option%nflowspec),auxvar_pert_up%diff(option%nflowspec))
    
    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    call THCAuxVarCopy(auxvar_up,auxvar_pert_up,option)
    call THCAuxVarCopy(auxvar_dn,auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
    
    x_up(1) = global_auxvar_up%pres(1)
    x_up(2) = global_auxvar_up%temp
    x_up(3) = auxvar_up%xmol(2)
    x_dn(1) = global_auxvar_dn%pres(1)
    x_dn(2) = global_auxvar_dn%temp
    x_dn(3) = auxvar_dn%xmol(2)
    do ideriv = 1,3
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_up(ideriv) = x_dn(ideriv)
      endif
    enddo
    if (option%use_th_freezing) then
       call THCAuxVarComputeIce(x_dn,auxvar_dn, &
            global_auxvar_dn,iphase,sat_func_dn, &
            0.d0,0.d0,option)
       call THCAuxVarComputeIce(x_up,auxvar_up, &
            global_auxvar_up,iphase,sat_func_dn, &
            0.d0,0.d0,option)
    else
       call THCAuxVarCompute(x_dn,auxvar_dn, &
            global_auxvar_dn,iphase,sat_func_dn, &
            0.d0,0.d0,option)
       call THCAuxVarCompute(x_up,auxvar_up, &
            global_auxvar_up,iphase,sat_func_dn, &
            0.d0,0.d0,option)
    endif
    
    call THCBCFlux(ibndtype,auxvars,auxvar_up,global_auxvar_up, &
                  auxvar_dn,global_auxvar_dn, &
                  por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                  area,dist_gravity,option,v_darcy,Diff_dn, &
                  res)
    if (ibndtype(THC_PRESSURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(THC_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    do ideriv = 1,3
      pert_dn = x_dn(ideriv)*perturbation_tolerance    
      x_pert_dn = x_dn
     
      if (option%use_th_freezing) then
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
      else
         x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
        
      x_pert_up = x_up
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_pert_up(ideriv) = x_pert_dn(ideriv)
      endif   

      if (option%use_th_freezing) then
         call THCAuxVarComputeIce(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn,iphase,sat_func_dn, &
              0.d0,0.d0,option)
         call THCAuxVarComputeIce(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up,iphase,sat_func_dn, &
              0.d0,0.d0,option)
      else
         call THCAuxVarCompute(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn,iphase,sat_func_dn, &
              0.d0,0.d0,option)
         call THCAuxVarCompute(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up,iphase,sat_func_dn, &
              0.d0,0.d0,option)
      endif

      call THCBCFlux(ibndtype,auxvars,auxvar_pert_up,global_auxvar_pert_up, &
                    auxvar_pert_dn,global_auxvar_pert_dn, &
                    por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                    area,dist_gravity,option,v_darcy,Diff_dn, &
                    res_pert_dn)
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    deallocate(auxvar_pert_dn%xmol,auxvar_pert_dn%diff)
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)      
  endif

end subroutine THCBCFluxDerivative

! ************************************************************************** !

subroutine THCBCFlux(ibndtype,auxvars,auxvar_up,global_auxvar_up, &
                    auxvar_dn,global_auxvar_dn, &
                    por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                    area,dist_gravity,option,v_darcy,Diff_dn, &
                    Res)
  ! 
  ! Computes the  boundary flux terms for the residual
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
  use Option_module
   
  use EOS_Water_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(thc_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn, Diff_dn
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: v_darcy, area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  
  ! ice variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_dn, Ddiffgas_up
  PetscReal :: p_g
  PetscReal :: deng_dn, deng_up
  PetscReal :: psat_dn, psat_up
  PetscReal :: molg_dn, molg_up
  PetscReal :: satg_dn, satg_up
  PetscReal :: Diffg_dn, Diffg_up
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: fv_up, fv_dn
  PetscReal, parameter :: R_gas_constant = 8.3144621 ! Gas constant in J/mol/K
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(THC_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dq = perm_dn / dd_up
      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
   
        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity

        dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity

        if (ibndtype(THC_PRESSURE_DOF) == SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_auxvar_up%pres(1) - option%reference_pressure < eps) then
            dphi = 0.d0
          endif
        endif
        
        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
        else
          ukvr = auxvar_dn%kvr
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(auxvars(THC_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(THC_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
        else 
          density_ave = global_auxvar_dn%den(1)
        endif 
      endif

  end select

  q = v_darcy * area

  if (v_darcy >= 0.D0) then
    uh = auxvar_up%h
    uxmol(:)=auxvar_up%xmol(1:option%nflowspec)
  else
    uh = auxvar_dn%h
    uxmol(:)=auxvar_dn%xmol(1:option%nflowspec)
  endif      
    
  do ispec=1, option%nflowspec 
    fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
  enddo
  fluxe = fluxe + q*density_ave*uh

  ! Diffusion term   
  select case(ibndtype(THC_CONCENTRATION_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
    
!      if (global_auxvar_up%sat > eps .and. global_auxvar_dn%sat > eps) then
!        diff = diffdp * 0.25D0*(global_auxvar_up%sat(1)+global_auxvar_dn%sat(1))* &
!               (global_auxvar_up%den+global_auxvar_dn%den)

!pcl  if (global_auxvar_dn%sat(1) > eps) then
        diff = diffdp * global_auxvar_dn%sat(1)*global_auxvar_dn%den(1)
 !         fluxm(ispec) = fluxm(ispec) + diff*auxvar_dn%diff(ispec)* &
!                          (auxvar_up%xmol(ispec)-auxvar_dn%xmol(ispec))
        fluxm(2) = fluxm(2) + diff*Diff_dn* &
                           (auxvar_up%xmol(2)-auxvar_dn%xmol(2))
!pcl  endif
        if (option%use_th_freezing) then
#ifdef REMOVE_SATURATION
           diff = diffdp*global_auxvar_dn%den(1)
           fluxm(2) = fluxm(2) + diff*Diff_dn* &
                (auxvar_up%xmol(2)-auxvar_dn%xmol(2))
#endif
        endif
  end select
  


  ! Conduction term
  select case(ibndtype(THC_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      write(*,*) 'Dk value in incorrect in THCBCFlux. Look at THBCFlux.'
      stop
      Dk = Dk_dn / dd_up
      cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)
      fluxe = fluxe + cond

      if (option%use_th_freezing) then
         ! Added by Satish Karra,
         satg_up = auxvar_up%sat_gas
         satg_dn = auxvar_dn%sat_gas
         if ((satg_up > eps) .and. (satg_dn > eps)) then
            p_g = option%reference_pressure ! set to reference pressure
            deng_up = p_g/(IDEAL_GAS_CONST*(global_auxvar_up%temp + 273.15d0))*1.d-3
            deng_dn = p_g/(IDEAL_GAS_CONST*(global_auxvar_dn%temp + 273.15d0))*1.d-3
        
            Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
            p_ref = 1.01325d5 ! in Pa
            T_ref = 25.d0 ! in deg C
        
            Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
            Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
            Ddiffgas_up = satg_up*deng_up*Diffg_up
            Ddiffgas_dn = satg_dn*deng_dn*Diffg_dn
            call EOSWaterSaturationPressure(global_auxvar_up%temp, psat_up, ierr)
            call EOSWaterSaturationPressure(global_auxvar_dn%temp, psat_dn, ierr)
      
            ! vapor pressure lowering due to capillary pressure
            fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
                 R_gas_constant*(global_auxvar_up%temp + 273.15d0)))
            fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
                 R_gas_constant*(global_auxvar_dn%temp + 273.15d0)))

            molg_up = psat_up*fv_up/p_g
            molg_dn = psat_dn*fv_dn/p_g
        
            if (molg_up > molg_dn) then 
               upweight = 0.d0
            else 
               upweight = 1.d0
            endif
        
            Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn 
            fluxm(1) = fluxm(1) + por_dn*tor_dn*Ddiffgas_avg*(molg_up - molg_dn)/ &
                 dd_up*area
         endif

      endif

    case(NEUMANN_BC)
      fluxe = fluxe + auxvars(THC_TEMPERATURE_DOF)*area*option%scale ! added by SK 10/18/11
    case(ZERO_GRADIENT_BC)
      ! No change in fluxe
  end select

  Res(1:option%nflowspec) = fluxm(:)
  Res(option%nflowdof) = fluxe

end subroutine THCBCFlux

! ************************************************************************** !

subroutine THCResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_class
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
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option
  
  ! Communication -----------------------------------------
  ! These 3 must be called before THCUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  
  ! Compute internal and boundary flux terms as well as source/sink terms
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THCResidualPatch(snes,xx,r,realization,ierr)
    cur_patch => cur_patch%next
  enddo

end subroutine THCResidual

! ************************************************************************** !

subroutine THCResidualPatch(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation at patch level
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  

  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  
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
  type(thc_parameter_type), pointer :: thc_parameter
  type(thc_auxvar_type), pointer :: auxvars(:), auxvars_bc(:)
  type(thc_auxvar_type), pointer :: auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set  
  type(sec_heat_type), pointer :: thc_sec_heat_vars(:)

  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: vol_frac_prim

  ! secondary continuum variables
  PetscReal :: sec_density
  PetscReal :: sec_dencpr
  PetscReal :: res_sec_heat
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  thc_parameter => patch%aux%THC%thc_parameter
  auxvars => patch%aux%THC%auxvars
  auxvars_bc => patch%aux%THC%auxvars_bc
  auxvars_ss => patch%aux%THC%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  thc_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  
  call THCUpdateAuxVarsPatch(realization)
  ! override flags since they will soon be out of date  
  patch%aux%THC%auxvars_up_to_date = PETSC_FALSE

  if (option%compute_mass_balance_new) then
    call THCZeroMassBalDeltaPatch(realization)
  endif


! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90( r, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
 
  call VecGetArrayF90(field%flow_yy,yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tortuosity_loc, tor_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)
  !print *,' Finished scattering non deriv'
  
  
  ! Calculating volume fractions for primary and secondary continua

  vol_frac_prim = 1.d0
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
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1


    if (option%use_mc) then
      vol_frac_prim = thc_sec_heat_vars(local_id)%epsilon
    endif

    call THCAccumulation(auxvars(ghosted_id),global_auxvars(ghosted_id), &
                        porosity_loc_p(ghosted_id), &
                        volume_p(local_id), &
                        thc_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                        option,vol_frac_prim,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res
  enddo 
#endif

! ================== Secondary continuum heat source terms =====================

#if 1
  if (option%use_mc) then
  ! Secondary continuum contribution (Added by SK 06/02/2012)
  ! only one secondary continuum for now for each primary continuum node
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
    
      sec_dencpr = thc_parameter%dencpr(int(ithrm_loc_p(local_id))) ! secondary rho*c_p same as primary for now
    
      call THCSecondaryHeat(thc_sec_heat_vars(local_id), &
                          global_auxvars(ghosted_id), &
!                         thc_parameter%ckdry(int(ithrm_loc_p(local_id))), &
                          thc_parameter%ckwet(int(ithrm_loc_p(local_id))), &
                          sec_dencpr, &
                          option,res_sec_heat)

      r_p(iend) = r_p(iend) - res_sec_heat*volume_p(local_id)
    enddo   
  endif
#endif

! ============== end secondary continuum heat source ===========================

#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%flow_condition%num_sub_conditions > THC_CONCENTRATION_DOF) then
      enthalpy_flag = PETSC_TRUE
    else
      enthalpy_flag = PETSC_FALSE
    endif

    qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)

    qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / FMWCO2
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      if (enthalpy_flag) then
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1
      endif

      if (qsrc1 > 0.d0) then ! injection
        ! units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
        ! qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        r_p((local_id-1)*option%nflowdof + jh2o) =  &
                                     r_p((local_id-1)*option%nflowdof + jh2o) &
                                     - qsrc1
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - &
          qsrc1*auxvars_ss(sum_connection)%h
      else
        ! extraction
        r_p((local_id)*option%nflowdof+jh2o) = r_p((local_id-1)*option%nflowdof+jh2o) &
                                               - qsrc1
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - &
                                        qsrc1*auxvars(ghosted_id)%h
                                        
      endif  
    
      if (csrc1 > 0.d0) then ! injection
        call printErrMsg(option,"concentration source not yet implemented in THC")
      endif
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
   
      D_up = thc_parameter%ckwet(ithrm_up)
      D_dn = thc_parameter%ckwet(ithrm_dn)
      
      Dk_dry_up = thc_parameter%ckdry(ithrm_up)
      Dk_dry_dn = thc_parameter%ckdry(ithrm_dn)
      
      alpha_up = thc_parameter%alpha(ithrm_up)
      alpha_dn = thc_parameter%alpha(ithrm_dn)

      if (option%use_th_freezing) then
         Dk_ice_up = thc_parameter%ckfrozen(ithrm_up)
         DK_ice_dn = thc_parameter%ckfrozen(ithrm_dn)
      
         alpha_fr_up = thc_parameter%alpha_fr(ithrm_up)
         alpha_fr_dn = thc_parameter%alpha_fr(ithrm_dn)
      else
         Dk_ice_up = Dk_dry_up
         Dk_ice_dn = Dk_dry_dn

         alpha_fr_up = alpha_up
         alpha_fr_dn = alpha_dn
      endif

      Diff_up = thc_parameter%diffusion_coefficient(1)
      Diff_dn = thc_parameter%diffusion_coefficient(1)
      

      call THCFlux(auxvars(ghosted_id_up),global_auxvars(ghosted_id_up), &
                  porosity_loc_p(ghosted_id_up), &
                  tor_loc_p(ghosted_id_up),thc_parameter%sir(1,icap_up), &
                  dd_up,perm_up,D_up, &
                  auxvars(ghosted_id_dn),global_auxvars(ghosted_id_dn), &
                  porosity_loc_p(ghosted_id_dn), &
                  tor_loc_p(ghosted_id_dn),thc_parameter%sir(1,icap_dn), &
                  dd_dn,perm_dn,D_dn, &
                  cur_connection_set%area(iconn),distance_gravity, &
                  upweight,option,v_darcy,Diff_up,Diff_dn,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
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
      D_dn = thc_parameter%ckwet(ithrm_dn)

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
	
      call THCBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id), &
                                thc_parameter%sir(1,icap_dn), &
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

  do local_id = 1, grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    r_p (istart:iend)= r_p(istart:iend)/volume_p(local_id)
  enddo

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

  if (patch%aux%THC%inactive_cells_exist) then
    do i=1,patch%aux%THC%n_zero_rows
      r_p(patch%aux%THC%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tor_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(option%mycomm,'THCresidual.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(option%mycomm,'THCxx.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif 

end subroutine THCResidualPatch

! ************************************************************************** !

subroutine THCJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  type(option_type),  pointer :: option
  PetscReal :: norm
  
  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THCJacobianPatch(snes,xx,J,J,realization,ierr)
    cur_patch => cur_patch%next
  enddo

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'THCjacobian.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option)
  endif
  
end subroutine THCJacobian

! ************************************************************************** !

subroutine THCJacobianPatch(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
       
  

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i
  PetscInt :: ip1, ip2 

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  PetscInt :: icap,iphas,icap_up,icap_dn
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
  type(thc_parameter_type), pointer :: thc_parameter
  type(thc_auxvar_type), pointer :: auxvars(:), auxvars_bc(:),auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 

  type(sec_heat_type), pointer :: sec_heat_vars(:)

  PetscViewer :: viewer
  Vec :: debug_vec
  PetscReal :: vol_frac_prim
  
  ! secondary continuum variables
  PetscReal :: area_prim_sec
  PetscReal :: jac_sec_heat

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  thc_parameter => patch%aux%THC%thc_parameter
  auxvars => patch%aux%THC%auxvars
  auxvars_bc => patch%aux%THC%auxvars_bc
  auxvars_ss => patch%aux%THC%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc

  sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  
#if 0
   call THCNumericalJacobianTest(xx,realization)
#endif

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tortuosity_loc, tor_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)
  
  vol_frac_prim = 1.d0

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
    
    if (option%use_mc) then    
      vol_frac_prim = sec_heat_vars(local_id)%epsilon
    endif
    
    call THCAccumDerivative(auxvars(ghosted_id),global_auxvars(ghosted_id), &
                            porosity_loc_p(ghosted_id), &
                            volume_p(local_id), &
                            thc_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                            option, &
                            realization%saturation_function_array(icap)%ptr, &
                            vol_frac_prim,Jup) 

    if (option%use_mc) then
      call THCSecondaryHeatJacobian(sec_heat_vars(local_id), &
                        thc_parameter%ckwet(int(ithrm_loc_p(local_id))), &
                        thc_parameter%dencpr(int(ithrm_loc_p(local_id))), &
                        option,jac_sec_heat)
                        
      Jup(option%nflowdof,2) = Jup(option%nflowdof,2) - &
                               jac_sec_heat*volume_p(local_id)
    endif
                            
! scale by the volume of the cell
    Jup = Jup/volume_p(local_id)
                           
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_accum.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%flow_condition%num_sub_conditions > THC_CONCENTRATION_DOF) then
      enthalpy_flag = PETSC_TRUE
    else
      enthalpy_flag = PETSC_FALSE
    endif

    qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)

    qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / FMWCO2
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif        

      if (qsrc1 > 0.d0) then ! injection
        ! units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
        ! qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        ! base on r_p() = r_p() - qsrc1*enth_src_h2o
        ! dresT_dp = -qsrc1*hw_dp
        dresT_dp = -qsrc1*auxvars_ss(sum_connection)%dh_dp
        ! dresT_dt = -qsrc1*hw_dt ! since tsrc1 is prescribed, there is no derivative
        istart = ghosted_id*option%nflowdof
        call MatSetValuesLocal(A,1,istart-1,1,istart-option%nflowdof,dresT_dp,ADD_VALUES, &
                               ierr);CHKERRQ(ierr)
        ! call MatSetValuesLocal(A,1,istart-1,1,istart-1,dresT_dt,ADD_VALUES,ierr)
      else
        ! extraction
        dresT_dp = -qsrc1*auxvars(ghosted_id)%dh_dp
        dresT_dt = -qsrc1*auxvars(ghosted_id)%dh_dt
        istart = ghosted_id*option%nflowdof
        call MatSetValuesLocal(A,1,istart-1,1,istart-option%nflowdof,dresT_dp,ADD_VALUES, &
                               ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(A,1,istart-1,1,istart-1,dresT_dt,ADD_VALUES, &
                               ierr);CHKERRQ(ierr)
      endif
    
      if (csrc1 > 0.d0) then ! injection
        call printErrMsg(option,"concentration source not yet implemented in THC")
      endif
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_srcsink.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
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
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))
    
      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      
      D_up = thc_parameter%ckwet(ithrm_up)
      D_dn = thc_parameter%ckwet(ithrm_dn)
    
      Dk_dry_up = thc_parameter%ckdry(ithrm_up)
      Dk_dry_dn = thc_parameter%ckdry(ithrm_dn)
      
      alpha_up = thc_parameter%alpha(ithrm_up)
      alpha_dn = thc_parameter%alpha(ithrm_dn)

      if (option%use_th_freezing) then
         Dk_ice_up = thc_parameter%ckfrozen(ithrm_up)
         DK_ice_dn = thc_parameter%ckfrozen(ithrm_dn)
      
         alpha_fr_up = thc_parameter%alpha_fr(ithrm_up)
         alpha_fr_dn = thc_parameter%alpha_fr(ithrm_dn)
      else
         Dk_ice_up = Dk_dry_up
         Dk_ice_dn = Dk_dry_dn
      
         alpha_fr_up = alpha_up
         alpha_fr_dn = alpha_dn
      endif

      Diff_up = thc_parameter%diffusion_coefficient(1)
      Diff_dn = thc_parameter%diffusion_coefficient(1)

      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
	  

      call THCFluxDerivative(auxvars(ghosted_id_up),global_auxvars(ghosted_id_up), &
                             porosity_loc_p(ghosted_id_up), &
                             tor_loc_p(ghosted_id_up),thc_parameter%sir(1,icap_up), &
                             dd_up,perm_up,D_up, &
                             auxvars(ghosted_id_dn),global_auxvars(ghosted_id_dn), &
                             porosity_loc_p(ghosted_id_dn), &
                             tor_loc_p(ghosted_id_dn),thc_parameter%sir(1,icap_dn), &
                             dd_dn,perm_dn,D_dn, &
                             cur_connection_set%area(iconn),distance_gravity, &
                             upweight,option, &
                             realization%saturation_function_array(icap_up)%ptr, &
                             realization%saturation_function_array(icap_dn)%ptr, &
                             Diff_up,Diff_dn,Dk_dry_up,Dk_dry_dn, &
                             Dk_ice_up,Dk_ice_dn, &
                             alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                             Jup,Jdn)
      
!  scale by the volume of the cell                      
      
      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup/volume_p(local_id_up),ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn/volume_p(local_id_up),ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn/volume_p(local_id_dn),ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup/volume_p(local_id_dn),ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_flux.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
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
      D_dn = thc_parameter%ckwet(ithrm_dn)

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
	  
      Diff_dn = thc_parameter%diffusion_coefficient(1)

      call THCBCFluxDerivative(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id), &
                                thc_parameter%sir(1,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
                                cur_connection_set%area(iconn), &
                                distance_gravity,option, &
                                realization%saturation_function_array(icap_dn)%ptr,&
                                Diff_dn,Jdn)
                                
      Jdn = -Jdn
  
!  scale by the volume of the cell
        Jdn = Jdn/volume_p(local_id)
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tor_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr);CHKERRQ(ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero, &
                        PETSC_NULL_VEC,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)
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
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
#else
  if (patch%aux%THC%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%THC%n_zero_rows, &
                          patch%aux%THC%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif
#endif

end subroutine THCJacobianPatch

! ************************************************************************** !

subroutine THCCreateZeroArray(patch,option)
  ! 
  ! Computes the zeroed rows for inactive grid cells
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

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

  patch%aux%THC%zero_rows_local => zero_rows_local
  patch%aux%THC%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%THC%n_zero_rows = n_zero_rows

  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
  if (flag > 0) patch%aux%THC%inactive_cells_exist = PETSC_TRUE

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine THCCreateZeroArray

! ************************************************************************** !

subroutine THCMaxChange(realization)
  ! 
  ! Computes the maximum change in the solution vector
  ! 
  ! Author: ???
  ! Date: 01/15/08
  ! 

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
  
  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy, &
                ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax, &
                     ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dtmpmax, &
                     ierr);CHKERRQ(ierr)
  if (option%nflowdof > 2) then
    call VecStrideNorm(field%flow_dxx,TWO_INTEGER,NORM_INFINITY,option%dcmax, &
                       ierr);CHKERRQ(ierr)
  endif
    
end subroutine THCMaxChange

! ************************************************************************** !

subroutine THCResidualToMass(realization)
  ! 
  ! Computes mass balance from residual equation
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module

  implicit none

  Vec :: ts_mass_balance
  type(realization_type) :: realization
  
  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  PetscReal, pointer :: mass_balance_p(:)
  type(thc_auxvar_type), pointer :: auxvars(:) 
  type(global_auxvar_type), pointer :: global_auxvars(:) 
  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  
  option => realization%option
  field => realization%field

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid
    auxvars => cur_patch%aux%THC%auxvars

    call VecGetArrayF90(field%flow_ts_mass_balance,mass_balance_p,  &
                        ierr);CHKERRQ(ierr)
  
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (cur_patch%imat(ghosted_id) <= 0) cycle
        
      istart = (ghosted_id-1)*option%nflowdof+1
      mass_balance_p(istart) = mass_balance_p(istart)/ &
                                global_auxvars(ghosted_id)%den(1)* &
                                global_auxvars(ghosted_id)%den_kg(1)
    enddo

    call VecRestoreArrayF90(field%flow_ts_mass_balance,mass_balance_p,  &
                            ierr);CHKERRQ(ierr)

    cur_patch => cur_patch%next
  enddo

end subroutine THCResidualToMass

! ************************************************************************** !

function THCGetTecplotHeader(realization,icolumn)
  ! 
  ! THCLiteGetTecplotHeader: Returns THC contribution to
  ! Tecplot file header
  ! 
  ! Author: ???
  ! Date: 02/13/08
  ! 

  use Realization_class
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: THCGetTecplotHeader
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
    write(string2,'('',"'',i2,''-Sl"'')') icolumn
  else
    write(string2,'('',"Sl"'')')
  endif
  string = trim(string) // trim(string2)

  if (option%use_th_freezing) then
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
  endif

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
  
  THCGetTecplotHeader = string

end function THCGetTecplotHeader

! ************************************************************************** !

subroutine THCSetPlotVariables(realization)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
  use Realization_class
  use Output_Aux_module
  use Variables_module

  implicit none

  type(realization_type) :: realization
  type(output_variable_type) :: output_variable
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => realization%output_option%output_variable_list
  
  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               TEMPERATURE)
  
  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)

  if (realization%option%use_th_freezing) then
     name = 'Gas Saturation'
     units = ''
     call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
          GAS_SATURATION)

     name = 'Ice Saturation'
     units = ''
     call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
          ICE_SATURATION)

     name = 'Ice Density'
     units = 'kg/m^3'
     call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
          ICE_DENSITY)
  endif

  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)

! name = 'Gas Density'
! units = ''
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              GAS_DENSITY)

  name = 'Liquid Energy'
  units = 'kJ/mol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)

! name = 'Gas Energy'
! units = ''
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              GAS_ENERGY)

  name = 'Liquid Viscosity'
  units = 'Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_VISCOSITY)

! name = 'Gas Viscosity'
! units = ''
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              GAS_VISCOSITY)

  name = 'Liquid Mobility'
  units = '1/Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOBILITY)

! name = 'Gas Mobility'
! units = ''
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              GAS_MOBILITY)

  name = 'Liquid Mole Fraction H2O'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION,ONE_INTEGER)

! name = 'Liquid Mole Fraction CO2'
! units = ''
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              LIQUID_MOLE_FRACTION,TWO_INTEGER)

  name = 'Gas Mole Fraction H2O'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOLE_FRACTION,ONE_INTEGER)

! name = 'Phase'
! units = ''
! output_variable%iformat = 1 ! integer
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              PHASE)

end subroutine THCSetPlotVariables

! ************************************************************************** !

subroutine THCComputeGradient(grid, global_auxvars, ghosted_id, gradient, &
                              option) 
  ! 
  ! Computes the gradient of temperature (for now) using
  ! least square fit of values from neighboring cells
  ! See:I. Bijelonja, I. Demirdzic, S. Muzaferija -- A finite volume method
  ! for incompressible linear elasticity, CMAME
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 2/20/12
  ! 

#include "finclude/petscdmda.h"
  use petscdmda
  use Grid_module
  use Global_Aux_module
  use Option_module
  use Utility_module

  implicit none

  type(option_type) :: option
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)

  
  PetscInt :: ghosted_neighbors_size, ghosted_id
  PetscInt :: ghosted_neighbors(26)
  PetscReal :: gradient(3), disp_vec(3,1), disp_mat(3,3)
  PetscReal :: temp_weighted(3,1)
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
  temp_weighted = 0.d0
  do i = 1, ghosted_neighbors_size
    disp_vec(1,1) = grid%x(ghosted_neighbors(i)) - grid%x(ghosted_id)
    disp_vec(2,1) = grid%y(ghosted_neighbors(i)) - grid%y(ghosted_id)
    disp_vec(3,1) = grid%z(ghosted_neighbors(i)) - grid%z(ghosted_id)
    disp_mat = disp_mat + matmul(disp_vec,transpose(disp_vec))
    temp_weighted = temp_weighted + disp_vec* &
                    (global_auxvars(ghosted_neighbors(i))%temp - &
                     global_auxvars(ghosted_id)%temp)
  enddo

  call ludcmp(disp_mat,THREE_INTEGER,INDX,D)
  call lubksb(disp_mat,THREE_INTEGER,INDX,temp_weighted)
  
  gradient(:) = temp_weighted(:,1)
  
end subroutine THCComputeGradient

! ************************************************************************** !

subroutine THCSecondaryHeat(sec_heat_vars,global_auxvar, &
                            therm_conductivity,dencpr, &
                            option,res_heat)
  ! 
  ! Calculates the source term contribution due to secondary
  ! continuum in the primary continuum residual
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/2/12
  ! 
                            
  use Option_module 
  use Global_Aux_module
  use Secondary_Continuum_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  PetscReal :: temp_current_N
  PetscReal :: res_heat
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_auxvar%temp
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! We only need the temperature at the outer-most node (closest to primary node)
  temp_current_N = rhs(ngcells)/coeff_diag(ngcells)
  
  ! Calculate the coupling term
  res_heat = area_fm*therm_conductivity*(temp_current_N - temp_primary_node)/ &
             dm_plus(ngcells)
                          
end subroutine THCSecondaryHeat

! ************************************************************************** !

subroutine THCSecondaryHeatJacobian(sec_heat_vars, &
                                    therm_conductivity, &
                                    dencpr, &
                                    option,jac_heat)
  ! 
  ! Calculates the source term jacobian contribution
  ! due to secondary continuum in the primary continuum residual
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/6/12
  ! 
                                    
  use Option_module 
  use Global_Aux_module
  use Secondary_Continuum_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: m
  PetscReal :: Dtemp_N_Dtemp_prim
  PetscReal :: jac_heat
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  area_fm = sec_heat_vars%interfacial_area
  dm_minus = sec_heat_vars%dm_minus
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                                        
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    ! We do not have to calculate rhs terms
  enddo

  ! We need the temperature derivative at the outer-most node (closest to primary node)
  Dtemp_N_Dtemp_prim = 1.d0/coeff_diag(ngcells)*alpha*area(ngcells)/ &
                       (dm_plus(ngcells)*vol(ngcells))
  
  ! Calculate the jacobian term
  jac_heat = area_fm*therm_conductivity*(Dtemp_N_Dtemp_prim - 1.d0)/ &
             dm_plus(ngcells)
                            
              
end subroutine THCSecondaryHeatJacobian                                  

! ************************************************************************** !

subroutine THCDestroy(patch)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  use Patch_module

  implicit none
  
  type(patch_type) :: patch
  
  ! need to free array in aux vars
  call THCAuxDestroy(patch%aux%THC)

end subroutine THCDestroy

end module THC_module
