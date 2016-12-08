module TOWG_module

  use PM_TOWG_Aux_module
  use AuxVars_TOWG_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petsclog.h"

#define CONVECTION
#define DIFFUSION
#define LIQUID_DIFFUSION
#define CONDUCTION

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  !pointing to null() function
  procedure(TOWGUpdateAuxVarsDummy), pointer :: TOWGUpdateAuxVars => null()
  procedure(TOWGAccumulationDummy), pointer :: TOWGAccumulation => null()

  abstract interface
    subroutine TOWGUpdateAuxVarsDummy(realization,update_state)
      use Realization_Subsurface_class  
      implicit none

      type(realization_subsurface_type) :: realization
      PetscBool :: update_state

    end subroutine TOWGUpdateAuxVarsDummy

    subroutine TOWGAccumulationDummy(auxvar,global_auxvar,material_auxvar, &
                                     soil_heat_capacity,option,Res,debug_cell)
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Option_module
      use Material_module
      use Material_Aux_class
      implicit none

      class(auxvar_towg_type) :: auxvar
      type(global_auxvar_type) :: global_auxvar
      class(material_auxvar_type) :: material_auxvar
      PetscReal :: soil_heat_capacity
      type(option_type) :: option
      PetscReal :: Res(option%nflowdof) 
      PetscBool :: debug_cell
    end subroutine TOWGAccumulationDummy
  end interface

#ifdef TOWG_DEBUG
  PetscInt, parameter :: debug_unit = 87
  PetscInt, parameter :: debug_info_unit = 86
  character(len=MAXWORDLENGTH) :: debug_filename
  PetscInt :: debug_flag = 0
  PetscInt :: debug_iteration_count
  PetscInt :: debug_timestep_cut_count
  PetscInt :: debug_timestep_count
#endif


  public :: TOWGSetup,&
            TOWGUpdateAuxVars, &
            TOWGUpdateSolution, &
            TOWGMapBCAuxVarsToGlobal, &
            TOWGInitializeTimestep

contains

! ************************************************************************** !

subroutine TOWGSetup(realization)
  ! 
  ! Creates arrays for TOWG auxiliary variables
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/07/16
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  !use Fluid_module
  use Material_Aux_class
  use Output_Aux_module

  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter
  type(output_variable_list_type), pointer :: list

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: num_bc_connection, num_ss_connection
  PetscInt :: i, idof, count
  PetscBool :: error_found
  PetscInt :: flag(10)

  class(material_auxvar_type), pointer :: material_auxvars(:)
  !type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%TOWG => TOWGAuxCreate(option)

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE

  if (minval(material_parameter%soil_residual_saturation(:,:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil residual saturation.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_heat_capacity(:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil heat capacity.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_thermal_conductivity(:,:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil thermal conductivity.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  
  material_auxvars => patch%aux%Material%auxvars
  flag = 0

  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'Non-initialized cell volume.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'Non-initialized porosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'Non-initialized tortuosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%soil_particle_density < 0.d0 .and. &
        flag(4) == 0) then
      flag(4) = 1
      option%io_buffer = 'Non-initialized soil particle density.'
      call printMsg(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'Non-initialized permeability.'
      call printMsg(option)
      flag(5) = 1
    endif
  enddo

  if (error_found .or. maxval(flag) > 0) then
    option%io_buffer = 'Material property errors found in TOWGSetup.'
    call printErrMsg(option)
  endif

  num_bc_connection = &
               CouplerGetNumConnectionsInList(patch%boundary_condition_list)

  num_ss_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)

  call patch%aux%TOWG%Init(grid,num_bc_connection,num_ss_connection,option)

  !initialise here diffusion coefficient when supporting the oil and gas
  ! diffusione terms
  ! initialize parameters
  !cur_fluid_property => realization%fluid_properties
  !do 
  !  if (.not.associated(cur_fluid_property)) exit
  !  patch%aux%TOWGl%parameter% &
  !    diffusion_coefficient(cur_fluid_property%phase_id) = &
  !      cur_fluid_property%diffusion_coefficient
  !  cur_fluid_property => cur_fluid_property%next
  !enddo  
  ! check whether diffusion coefficients are initialized.
  ! check initialisation of diffusiont coefficient phase by phase as for 
  ! example below
  !if (Uninitialized(patch%aux%TOWG%parameter% &
  !    diffusion_coefficient(option%liquid_phase))) then
  !  option%io_buffer = &
  !    UninitializedMessage('Liquid phase diffusion coefficient','')
  !  call printErrMsg(option)
  !endif
  !

  list => realization%output_option%output_snap_variable_list
  call TOWGSetPlotVariables(list)
  list => realization%output_option%output_obs_variable_list
  call TOWGSetPlotVariables(list) 

  ! covergence creteria to be chosen (can use TOUGH or TOWG type) 
  ! set up here tough convergnce creteria if needed.

  !set up TOWG functions that can vary with the miscibility model
  select case(towg_miscibility_model)
    case(TOWG_IMMISCIBLE)
      TOWGUpdateAuxVars => TOWGImsUpdateAuxVars
      TOWGAccumulation => TOWGImsTLAccumulation
      call TOWGImsAuxVarComputeSetup()
    case default
      option%io_buffer = 'TOWGSetup: only TOWG_IMMISCIBLE is supported.'
      call printErrMsg(option)
  end select

end subroutine TOWGSetup

! ************************************************************************** !

subroutine TOWGImsUpdateAuxVars(realization,update_state)
  ! 
  ! Updates the auxiliary variables associated with the TOWG_IMS problem
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/30/16
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_class
  use EOS_Water_module
  use Saturation_Function_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)  

  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end
  PetscInt :: iphasebc, iphase
  PetscInt :: offset
  PetscInt :: istate
  PetscReal :: gas_pressure, capillary_pressure, liquid_saturation
  PetscReal :: saturation_pressure, temperature
  PetscInt :: real_index, variable
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)

!#define DEBUG_AUXVARS
#ifdef DEBUG_AUXVARS
  character(len=MAXWORDLENGTH) :: word
  PetscInt, save :: icall = 0
#endif

  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  towg => patch%aux%TOWG
  !gen_auxvars => patch%aux%General%auxvars
  !gen_auxvars_bc => patch%aux%General%auxvars_bc

  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

#ifdef DEBUG_AUXVARS
  icall = icall + 1
  write(word,*) icall
  word = 'genaux' // trim(adjustl(word))
#endif

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    !TOWG_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = TOWG_UPDATE_FOR_FIXED_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
 
    !to replace with TOWGImsAuxVarCompute (there is no need for pointer here)
    call TOWGAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                        towg%auxvars(ZERO_INTEGER,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        patch%characteristic_curves_array( &
                        patch%sat_func_id(ghosted_id))%ptr, &
                       natural_id, &
                       option)
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      !geh: negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id) 
      offset = (ghosted_id-1)*option%nflowdof
      if (patch%imat(ghosted_id) <= 0) cycle

      xxbc(:) = xx_loc_p(offset+1:offset+option%nflowdof)
      istate = boundary_condition%flow_aux_int_var(TOWG_STATE_INDEX,iconn)
      if (istate == TOWG_NULL_STATE) then !this is applied to flux (Neumann) conditions
        istate = global_auxvars(ghosted_id)%istate !look into the state of down cell
        select case(istate)
          !only two states are possible for TOWGIms: TOWG_THREE_PHASE_STATE, TOWG_NULL_STATE
          case(TOWG_THREE_PHASE_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(HYDROSTATIC_BC)
                  real_index = &
                    boundary_condition% &
                    flow_aux_mapping(towg_dof_to_primary_variable(idof,istate))
                  xxbc(idof) = &
                    boundary_condition%flow_aux_real_var(real_index,iconn)
                case(DIRICHLET_BC)
                  variable = towg_dof_to_primary_variable(idof,istate)
                  select case(variable)
                    ! for oil pressure dof
                    case(TOWG_OIL_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'TOWG Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs oil pressure defined.'
                        call printErrMsg(option)
                      endif
                    ! for oil saturation dof
                    case(TOWG_OIL_SATURATION_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        ! should be able to use oil saturation in DN cell
                        !option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                        !  trim(boundary_condition%flow_condition%name) // &
                        !  '" needs oil saturation defined.'
                        !call printErrMsg(option)
                      endif
                    case(TOWG_GAS_SATURATION_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        ! should be able to use oil saturation in DN cell
                        !option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                        !  trim(boundary_condition%flow_condition%name) // &
                        !  '" needs gas saturation defined.'
                        !call printErrMsg(option)
                      endif
                    case(TOWG_TEMPERATURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'TOWG Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs temperature defined.'
                        call printErrMsg(option)
                      endif
                  end select
                case(NEUMANN_BC)
                case default
                  option%io_buffer = 'Unknown BC type in TOWGUpdateAuxVars().'
                  call printErrMsg(option)
              end select
            enddo  
        end select
      else
        ! we do this for all BCs; Neumann bcs will be set later
        do idof = 1, option%nflowdof
          real_index = boundary_condition% &
              flow_aux_mapping(towg_dof_to_primary_variable(idof,istate))
          if (real_index > 0) then
            xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
          else
            option%io_buffer = 'Error setting up boundary ' // &
                                'condition in TOWGUpdateAuxVars'
            call printErrMsg(option)
          endif
        enddo
      endif
          
      ! set this based on data given 
      global_auxvars_bc(sum_connection)%istate = istate
      ! TOWG_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = TOWG_UPDATE_FOR_BOUNDARY
      call TOWGAuxVarCompute(xxbc,towg%auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%TOWG%auxvars_up_to_date = PETSC_TRUE

end subroutine TOWGImsUpdateAuxVars

! ************************************************************************** !

subroutine TOWGInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/07/16 
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call TOWGUpdateFixedAccum(realization)
  
#ifdef TOWG_DEBUG
  debug_flag = 0
  if (.true.) then
    debug_iteration_count = 0
    debug_flag = 1
  endif
  debug_iteration_count = 0
#endif

end subroutine TOWGInitializeTimestep

! ************************************************************************** !

subroutine TOWGUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  towg => patch%aux%TOWG
  global_auxvars => patch%aux%Global%auxvars
  
  if (realization%option%compute_mass_balance_new) then
    call TOWGUpdateMassBalance(realization)
  endif
  
  ! update stored state - currently not needed 
  !  - if required define istate_store in the auxvar_towg extension 
  !do ghosted_id = 1, grid%ngmax
  !  towg%auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS) = &
  !    global_auxvars(ghosted_id)%istate
  !enddo
    
end subroutine TOWGUpdateSolution

! ************************************************************************** !

subroutine TOWGZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array  
  ! PO: identical for many flow modes (Genral, Toil_Ims, TOWG), where can it 
  !     be located to be shared?? flow_mode_common.F90 ?
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%TOWG%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%TOWG%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine TOWGZeroMassBalanceDelta

! ************************************************************************** !

subroutine TOWGUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  
  PetscInt :: iconn
  PetscInt :: icomp

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%TOWG%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        towg_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%TOWG%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        towg_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine TOWGUpdateMassBalance

! ************************************************************************** !

subroutine TOWGMapBCAuxVarsToGlobal(realization)
  !
  ! Map BC auxvars for TOWG problem coupled with reactive transport
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update
  
  towg => patch%aux%TOWG
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        towg%auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        towg%auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        towg%auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine TOWGMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine TOWGSetPlotVariables(list)
  ! 
  ! Adds variables to be printed to list for TOWG module
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/29/16
  ! 
  
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(output_variable_list_type), pointer :: list
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_type), pointer :: output_variable

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

  name = 'Oil Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               OIL_PRESSURE)

  name = 'Gas Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               GAS_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)
  
  name = 'Oil Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               OIL_SATURATION)

  name = 'Gas Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               GAS_SATURATION)

  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)
  
  name = 'Oil Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               OIL_DENSITY)

  name = 'Gas Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_DENSITY)
  
  name = 'Liquid Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)

  name = 'Oil Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               OIL_ENERGY)

  name = 'Gas Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_ENERGY)

  !to be added later for BLACK oil model and solvent model
  if ( towg_miscibility_model == TOWG_BLACK_OIL .or. &
       towg_miscibility_model == TOWG_SOLVENT_TL ) then
     write(*,*) "error: TOWGSetPlotVariables: TOWG_BLACK_OIL and " // &
                "TOWG_SOLVENT_TL not currently supported"
     stop
     ! add the following output variables to print:
     ! - gas mol fraction in oil (dissolved gas)
     ! - oil mol fraction in oil
     ! - gas mol fraction in gas
     ! - oil mol fraction in gas (oil vaour)
  end if

  if ( towg_miscibility_model == TOWG_SOLVENT_TL ) then
    ! add solvent saturation 
     write(*,*) "error: TOWGSetPlotVariables: TOWG_SOLVENT_TL " // &
                "not currently supported"
     stop
  end if

 ! to switch on when TOWG_BLACK_OIL is implemented 
 ! name = 'Thermodynamic State'
 ! units = ''
 ! output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE,units,STATE)
 ! output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
 ! output_variable%iformat = 1 ! integer
 ! call OutputVariableAddToList(list,output_variable)   
  
end subroutine TOWGSetPlotVariables

! ************************************************************************** !

subroutine TOWGUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/07/16 
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  towg => patch%aux%TOWG
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  !IF towg convergence criteria required: 
  !   initialize dynamic accumulation term for every p iteration step
  !if (towg_tough2_conv_criteria) then
  !  call VecGetArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! TOWG_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = TOWG_UPDATE_FOR_FIXED_ACCUM
    call TOWGAuxVarCompute(xx_p(local_start:local_end), &
                           towg%auxvars(ZERO_INTEGER,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           patch%characteristic_curves_array( &
                           patch%sat_func_id(ghosted_id))%ptr, &
                           natural_id, &
                           option)
    call TOWGAccumulation(towg%auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          material_parameter%soil_heat_capacity(imat), &
                          option,accum_p(local_start:local_end), &
                          local_id == towg_debug_cell_id) 
  enddo
  
  !for tough2 convergence criteria
  !if (towg_tough2_conv_criteria) then
  !  accum_p2 = accum_p
  !endif
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  !tough2 convergence criteria:
  ! initialize dynamic accumulation term for every p iteration step
  !if (towg_tough2_conv_criteria) then
  !  call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif
  
end subroutine TOWGUpdateFixedAccum

! ************************************************************************** !

subroutine TOWGImsTLAccumulation(auxvar,global_auxvar,material_auxvar, &
                                 soil_heat_capacity,option,Res,debug_cell)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual - TOWG_IMS and TOWG_TL models
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/07/16 
  ! 

  use Option_module
  use Material_module
  use Material_Aux_class
  
  implicit none

  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  PetscBool :: debug_cell
  
  PetscInt :: iphase, energy_id
  
  PetscReal :: porosity
  PetscReal :: volume_over_dt
 
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use gen_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = auxvar%effective_porosity
  
  ! accumulation term units = kmol/s 
  ! not for TOWG IMS and TL (kmol phase) = (kmol comp)
  ! and nphase = nflowspec
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol phase/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] * 
    Res(iphase) = Res(iphase) + auxvar%sat(iphase) * &
                                auxvar%den(iphase) 
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * volume_over_dt

  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
    Res(energy_id) = Res(energy_id) + auxvar%sat(iphase) * &
                                      auxvar%den(iphase) * &
                                      auxvar%U(iphase)
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * auxvar%temp) * volume_over_dt
  
#ifdef TOWG_DEBUG
  if (debug_flag > 0) then
    write(debug_unit,'(a,7es24.15)') 'accum:', Res
  endif
#endif                    

end subroutine TOWGImsTLAccumulation

! ************************************************************************** !



end module TOWG_module
