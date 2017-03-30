module TOWG_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_TOWG_Aux_module
  use AuxVars_TOWG_module
  use Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#define CONVECTION
#define DIFFUSION
#define LIQUID_DIFFUSION
#define CONDUCTION

!#define DEBUG_TOWG_FILEOUTPUT
!#define DEBUG_TOWG_FLUXES  

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

#ifdef DEBUG_TOWG_FILEOUTPUT
  PetscInt, parameter :: debug_unit = 87
  PetscInt, parameter :: debug_info_unit = 86
  character(len=MAXWORDLENGTH) :: debug_filename
  PetscInt :: debug_flag = 0
  PetscInt :: debug_iteration_count
  PetscInt :: debug_timestep_cut_count
  PetscInt :: debug_timestep_count
#endif

  !pointing to null() function
  procedure(TOWGUpdateAuxVarsDummy), pointer :: TOWGUpdateAuxVars => null()
  procedure(TOWGAccumulationDummy), pointer :: TOWGAccumulation => null()
  procedure(TOWGComputeMassBalanceDummy), pointer :: &
                                           TOWGComputeMassBalance => null()
  procedure(TOWGFluxDummy), pointer :: TOWGFlux => null()
  procedure(TOWGBCFluxDummy), pointer :: TOWGBCFlux => null()
  procedure(TOWGSrcSinkDummy), pointer :: TOWGSrcSink => null()
  procedure(TOWGCheckUpdatePreDummy), pointer :: TOWGCheckUpdatePre => null()
  procedure(TOWGMaxChangeDummy), pointer :: TOWGMaxChange => null()

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

    subroutine TOWGComputeMassBalanceDummy(realization,mass_balance)
      use Realization_Subsurface_class 
      implicit none
      type(realization_subsurface_type) :: realization
      PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)
    end subroutine TOWGComputeMassBalanceDummy

    subroutine TOWGFluxDummy(auxvar_up,global_auxvar_up, &
                             material_auxvar_up, &
                             sir_up, &
                             thermal_conductivity_up, &
                             auxvar_dn,global_auxvar_dn, &
                             material_auxvar_dn, &
                             sir_dn, &
                             thermal_conductivity_dn, &
                             area, dist, towg_parameter, &
                             option,v_darcy,Res, &
                             debug_connection)
      use PM_TOWG_Aux_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_class
      implicit none
      class(auxvar_towg_type) :: auxvar_up, auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: sir_up(:), sir_dn(:)
      PetscReal :: v_darcy(option%nphase)
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(towg_parameter_type) :: towg_parameter
      PetscReal :: thermal_conductivity_dn(2)
      PetscReal :: thermal_conductivity_up(2)
      PetscReal :: Res(option%nflowdof)
      PetscBool :: debug_connection
    end subroutine TOWGFluxDummy

    subroutine TOWGBCFluxDummy(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                               auxvar_up,global_auxvar_up, &
                               auxvar_dn,global_auxvar_dn, &
                               material_auxvar_dn, &
                               sir_dn, &
                               thermal_conductivity_dn, &
                               area,dist,towg_parameter, &
                               option,v_darcy,Res,debug_connection)
      use PM_TOWG_Aux_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Option_module                              
      use Material_Aux_class
      implicit none
      class(auxvar_towg_type) :: auxvar_up, auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_dn
      type(option_type) :: option
      PetscReal :: sir_dn(:)
      PetscReal :: bc_auxvars(:)
      PetscReal :: v_darcy(option%nphase), area
      type(towg_parameter_type) :: towg_parameter
      PetscReal :: dist(-1:3)
      PetscReal :: Res(1:option%nflowdof)
      PetscInt :: ibndtype(1:option%nflowdof)
      PetscInt :: bc_auxvar_mapping(TOWG_MAX_INDEX)
      PetscReal :: thermal_conductivity_dn(2)
      PetscBool :: debug_connection
    end subroutine TOWGBCFluxDummy

    subroutine TOWGSrcSinkDummy(option,src_sink_condition, auxvar, &
                            global_auxvar,ss_flow_vol_flux,scale,Res)
      use Option_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Condition_module  
      implicit none
      type(option_type) :: option
      type(flow_towg_condition_type), pointer :: src_sink_condition
      class(auxvar_towg_type) :: auxvar
      type(global_auxvar_type) :: global_auxvar
      PetscReal :: ss_flow_vol_flux(option%nphase)
      PetscReal :: scale  
      PetscReal :: Res(option%nflowdof)
    end subroutine TOWGSrcSinkDummy

    subroutine TOWGCheckUpdatePreDummy(line_search,X,dX,changed,realization, &
                                       max_it_before_damping,damping_factor, &
                                       max_pressure_change,ierr)
      use Realization_Subsurface_class
#include "petsc/finclude/petscsnes.h"
      use petscsnes
      implicit none
      SNESLineSearch :: line_search
      Vec :: X
      Vec :: dX
      PetscBool :: changed
      PetscErrorCode :: ierr
      type(realization_subsurface_type) :: realization
      PetscInt :: max_it_before_damping
      PetscReal :: damping_factor
      PetscReal :: max_pressure_change
    end subroutine TOWGCheckUpdatePreDummy

    subroutine TOWGMaxChangeDummy(realization,max_change_ivar, &
                                  max_change_isubvar,max_pressure_change, &
                                  max_xmol_change,max_saturation_change, &
                                  max_temperature_change)
      use Realization_Subsurface_class
      implicit none
      class(realization_subsurface_type), pointer :: realization
      PetscInt :: max_change_ivar(:)
      PetscInt :: max_change_isubvar(:)
      PetscReal :: max_pressure_change
      PetscReal :: max_xmol_change
      PetscReal :: max_saturation_change
      PetscReal :: max_temperature_change
    end subroutine TOWGMaxChangeDummy

  end interface 

  public :: TOWGSetup,&
            TOWGUpdateAuxVars, &
            TOWGUpdateSolution, &
            TOWGMapBCAuxVarsToGlobal, &
            TOWGInitializeTimestep, &
            TOWGComputeMassBalance, &
            TOWGResidual, &
            TOWGJacobian, &
            TOWGCheckUpdatePre, &
            TOWGTimeCut, &
            TOWGMaxChange, &
            TOWGDestroy

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
      TOWGComputeMassBalance => TOWGImsTLComputeMassBalance
      TOWGFlux => TOWGImsTLFlux
      TOWGBCFlux => TOWGImsTLBCFlux
      TOWGSrcSink => TOWGImsTLSrcSink
      TOWGCheckUpdatePre => TOWGImsTLCheckUpdatePre
      TOWGMaxChange => TOWGImsTLMaxChange
      call TOWGImsAuxVarComputeSetup()
    case default
      option%io_buffer = 'TOWGSetup: only TOWG_IMMISCIBLE is supported.'
      call printErrMsg(option)
  end select

#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_flag = 0
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = 0
  ! create new file
  open(debug_info_unit, file='debug_towg_info.txt', action="write", &
       status="unknown")
  write(debug_info_unit,*) 'type timestep cut iteration'
  close(debug_info_unit)
#endif  

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
  word = 'towgaux' // trim(adjustl(word))
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
  
#ifdef DEBUG_TOWG_FILEOUTPUT
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
  
  do ghosted_id = 1, grid%ngmax
    towg%auxvars(ZERO_INTEGER,ghosted_id)%istate_store(TOWG_PREV_TS) = &
      global_auxvars(ghosted_id)%istate
  enddo

#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = debug_timestep_count + 1
#endif 
    
end subroutine TOWGUpdateSolution

! ************************************************************************** !

subroutine TOWGTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/16
  ! 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(pm_towg_aux_type), pointer :: towg
  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  towg => patch%aux%TOWG

  ! restore stored state
  do ghosted_id = 1, grid%ngmax
    global_auxvars(ghosted_id)%istate = &
      towg%auxvars(ZERO_INTEGER,ghosted_id)%istate_store(TOWG_PREV_TS)
  enddo

#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_timestep_cut_count = debug_timestep_cut_count + 1
#endif 

  call TOWGInitializeTimestep(realization)  

end subroutine TOWGTimeCut

! ************************************************************************** !

subroutine TOWGZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array  
  ! PO: identical for many flow modes (Genral, Toil_Ims, TOWG), where can it 
  !     be located to be shared?? flow_mode_common.F90 ?
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/08/16 
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

subroutine TOWGImsTLComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/21/16
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  class(pm_towg_aux_type), pointer :: towg
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase, icomp
  PetscReal :: vol_phase

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  towg => patch%aux%TOWG
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    !not for both TOWG_IMS and TL phases and components coincides
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = &
        towg%auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        towg%auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density 
      mass_balance(iphase,1) = mass_balance(iphase,1) + &
        towg%auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)* &
        towg_fmw_comp(iphase) * vol_phase
    enddo
  enddo

end subroutine TOWGImsTLComputeMassBalance

! ************************************************************************** !

subroutine TOWGImsTLMaxChange(realization,max_change_ivar,max_change_isubvar,&
                              max_pressure_change,max_xmol_change, &
                              max_saturation_change,max_temperature_change)
  ! 
  ! Compute primary variable max changes
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/16
  ! 

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  !use Global_Aux_module
  !use General_Aux_module
  !use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
  !                             TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                             GAS_SATURATION
  implicit none
  
  class(realization_subsurface_type), pointer :: realization
  PetscInt :: max_change_ivar(:)
  PetscInt :: max_change_isubvar(:)
  PetscReal :: max_pressure_change
  PetscReal :: max_xmol_change
  PetscReal :: max_saturation_change
  PetscReal :: max_temperature_change

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(4)
  PetscReal :: max_change_global(4)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global = 0.d0
  max_change_local = 0.d0
  
  ! 'TOWG_IMMISCIBLE','TODD_LONGOSTAFF'
  ! max change variables = [OIL_PRESSURE, OIL_SATURATION, &
  !                         GAS_SATURATION,TEMPERATURE]
  do i = 1, 4
    call RealizationGetVariable(realization,field%work, &
                                max_change_ivar(i),max_change_isubvar(i))
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo
    max_change_local(i) = max_change
    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(max_change_local,max_change_global,FOUR_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4, &
      & " dsg= ",1pe12.4," dt= ",1pe12.4)') &
      max_change_global(1:4)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4,&
      & " dsg= ",1pe12.4,/,15x,"  dt= ",1pe12.4)') &
      max_change_global(1:4)
  endif
 
  ! 'TOWG_IMMISCIBLE','TODD_LONGOSTAFF'
  ! max change variables = [OIL_PRESSURE, OIL_SATURATION, &
  !                         GAS_SATURATION,TEMPERATURE]
  max_pressure_change = max_change_global(1)
  max_xmol_change = 0.0d0
  max_saturation_change = maxval(max_change_global(2:3))
  max_temperature_change = max_change_global(4)
  
end subroutine TOWGImsTLMaxChange

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
  
#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,7es24.15)') 'accum:', Res
  endif
#endif                    

end subroutine TOWGImsTLAccumulation

! ************************************************************************** !

subroutine TOWGImsTLFlux(auxvar_up,global_auxvar_up, &
                         material_auxvar_up, &
                         sir_up, &
                         thermal_conductivity_up, &
                         auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         sir_dn, &
                         thermal_conductivity_dn, &
                         area, dist, towg_parameter, &
                         option,v_darcy,Res, &
                         debug_connection)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/08/16 
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  !use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  class(auxvar_towg_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscInt :: energy_id
  PetscInt :: iphase
  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure, delta_temp
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_liquid
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  !for debugging
  PetscReal :: adv_flux(option%nflowdof)
  PetscReal :: debug_flux(3), debug_dphi(3)
  
  PetscReal :: dummy_dperm_up, dummy_dperm_dn
  PetscReal :: temp_perm_up, temp_perm_dn

  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  ! Fracture permeability change only available for structured grid (Heeho)
  !if (associated(material_auxvar_up%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_up,perm_up,temp_perm_up, &
  !                            dummy_dperm_up,dist)
  !  perm_up = temp_perm_up
  !endif
  !if (associated(material_auxvar_dn%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
  !                            dummy_dperm_dn,dist)
  !  perm_dn = temp_perm_dn
  !endif
  
  !if (associated(klinkenberg)) then
  !  perm_ave_over_dist(1) = (perm_up * perm_dn) / &
  !                          (dist_up*perm_dn + dist_dn*perm_up)
  !  temp_perm_up = klinkenberg%Evaluate(perm_up, &
  !                                       auxvar_up%pres(option%gas_phase))
  !  temp_perm_dn = klinkenberg%Evaluate(perm_dn, &
  !                                       auxvar_dn%pres(option%gas_phase))
  !  perm_ave_over_dist(2) = (temp_perm_up * temp_perm_dn) / &
  !                          (dist_up*temp_perm_dn + dist_dn*temp_perm_up)
  !else
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
  !endif
      
  Res = 0.d0
  
  v_darcy = 0.d0
#ifdef DEBUG_FLUXES  
  adv_flux = 0.d0
#endif
#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

#ifdef CONVECTION
  do iphase = 1, option%nphase
 
    if (auxvar_up%mobility(iphase) + &
        auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif

    density_kg_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                             auxvar_dn%sat(iphase), &
                                             auxvar_up%den_kg(iphase), &
                                             auxvar_dn%den_kg(iphase) )

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = auxvar_up%pres(iphase) - &
                     auxvar_dn%pres(iphase) + &
                     gravity_term

#ifdef TOWG_DEBUG
      debug_dphi(iphase) = delta_pressure
#endif

    if (delta_pressure >= 0.D0) then
      mobility = auxvar_up%mobility(iphase)
      H_ave = auxvar_up%H(iphase)
      uH = H_ave
    else
      mobility = auxvar_dn%mobility(iphase)
      H_ave = auxvar_dn%H(iphase)
      uH = H_ave
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                            auxvar_dn%sat(iphase), &
                                            auxvar_up%den(iphase), &
                                            auxvar_dn%den(iphase) )
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave

      ! Res[kmol total/sec] = mole_flux[kmol phase/sec]
      Res(iphase) = mole_flux

      !do icomp = 1, option%nflowspec
      !  ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
      !  !                      xmol[kmol comp/kmol phase]
      !  Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      !enddo

#ifdef DEBUG_FLUXES  
      adv_flux(iphase) = mole_flux
      !do icomp = 1, option%nflowspec
      !  adv_flux(icomp) = adv_flux(icomp) + mole_flux * xmol(icomp)
      !enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
#endif

      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH

#ifdef DEBUG_FLUXES  
      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
#endif

#ifdef DEBUG_TOWG_FILEOUTPUT
      debug_dphi(iphase) = delta_pressure
      debug_flux(iphase) = mole_flux * uH
#endif
    endif                   

  enddo
! CONVECTION
#endif

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'delta pressure :', debug_dphi(:)
    write(debug_unit,'(a,7es24.15)') 'adv flux (energy):', debug_flux(:)
  endif
  debug_flux = 0.d0
#endif                    


#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry) 
  ! Assuming that oil and water have same conductivity:
  ! s_l = S_water + S_oil 
  sat_liquid = auxvar_up%sat(option%liquid_phase) + &
               auxvar_up%sat(option%oil_phase)
  k_eff_up = thermal_conductivity_up(1) + &
             sqrt(sat_liquid) * &
             (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  sat_liquid = auxvar_dn%sat(option%liquid_phase) + &
               auxvar_dn%sat(option%oil_phase)
  k_eff_dn = thermal_conductivity_dn(1) + &
             sqrt(sat_liquid) * &
             (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  if (k_eff_up > 0.d0 .or. k_eff_dn > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  else
    k_eff_ave = 0.d0
  endif
  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = auxvar_up%temp - auxvar_dn%temp
  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux
! CONDUCTION
#endif
  
#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
!    write(*,'(a,7es12.4)') 'in: ', adv_flux(:)*dist(1), diff_flux(:)*dist(1)
    write(*,'('' phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(1), auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(1), auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(1) * 1.d6
    write(*,'('' phase: oil'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(2), auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(2), auxvar_dn%sat(2)
    write(*,'(''  oil --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(2) * 1.d6
    write(*,'('' phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(3), auxvar_dn%pres(3)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(3), auxvar_dn%sat(3)
    write(*,'(''  gas --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(3)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(3) * 1.d6
    write(*,'(''  energy --'')')
    write(*,'(''   advective heat flux:'',es12.4)') adv_flux(4) * 1.d6
    write(*,'(''   conductive heat flux:'',es12.4)') heat_flux * 1.d6
    write(*,'(''   total heat flux:'',es12.4)') (heat_flux + adv_flux(4))*1.d6

  endif
#endif

end subroutine TOWGImsTLFlux

! ************************************************************************** !

subroutine TOWGImsTLBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                           auxvar_up,global_auxvar_up, &
                           auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           sir_dn, &
                           thermal_conductivity_dn, &
                           area,dist,towg_parameter, &
                           option,v_darcy,Res,debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/22/16
  ! 
  use Option_module                              
  use Material_Aux_class
  !use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  class(auxvar_towg_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: bc_auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy(option%nphase), area
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: bc_auxvar_mapping(TOWG_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)
  PetscBool :: debug_connection
  
  PetscInt :: energy_id
  PetscInt :: iphase !, icomp
  PetscInt :: bc_type
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_temp !, delta_xmol
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux

  PetscReal :: adv_flux(option%nflowdof)
  PetscReal :: debug_flux(3), debug_dphi(3)

  PetscReal :: boundary_pressure
  PetscReal :: sat_liquid

  PetscReal :: dden_dn, dden_up

  PetscInt :: idof
  PetscBool :: neumann_bc_present
  
  PetscReal :: temp_perm_dn
  PetscReal :: dummy_dperm_dn
  
  energy_id = option%energy_id

  Res = 0.d0
  v_darcy = 0.d0  

#ifdef DEBUG_FLUXES  
  adv_flux = 0.d0
#endif
#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

  neumann_bc_present = PETSC_FALSE
  
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  ! Fracture permeability change only available for structured grid (Heeho)
  !if (associated(material_auxvar_dn%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
  !                            dummy_dperm_dn,dist)
  !  perm_dn = temp_perm_dn
  !endif  
  
  !if (associated(klinkenberg)) then
  !  perm_dn_adj(1) = perm_dn
  !                                        
  !  perm_dn_adj(2) = klinkenberg%Evaluate(perm_dn, &
  !                                        gen_auxvar_dn%pres(option%gas_phase))
  !else
    perm_dn_adj(:) = perm_dn
  !endif
  
#ifdef CONVECTION  
  do iphase = 1, option%nphase
 
    bc_type = ibndtype(iphase)
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then
          select case(option%phase_map(iphase)) 
            case(LIQUID_PHASE)
              idof = bc_auxvar_mapping(TOWG_LIQ_CONDUCTANCE_INDEX)
            case(OIL_PHASE)
              idof = bc_auxvar_mapping(TOWG_OIL_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = bc_auxvar_mapping(TOWG_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = bc_auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
          
        ! using residual saturation cannot be correct! - geh
        ! reusing sir_dn for bounary auxvar
#define BAD_MOVE1 ! this works
#ifndef BAD_MOVE1       
        if (auxvar_up%sat(iphase) > sir_dn(iphase) .or. &
            auxvar_dn%sat(iphase) > sir_dn(iphase)) then
#endif
          boundary_pressure = auxvar_up%pres(iphase)

          !PO: no free surfce boundaries considered           
          !if (iphase == LIQUID_PHASE .and. &
          !    global_auxvar_up%istate == GAS_STATE) then
          !  ! the idea here is to accommodate a free surface boundary
          !  ! face.  this will not work for an interior grid cell as
          !  ! there should be capillary pressure in force.
          !  boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
          !endif

          density_kg_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                                   auxvar_dn%sat(iphase), &
                                                   auxvar_up%den_kg(iphase), &
                                                   auxvar_dn%den_kg(iphase) )


          gravity_term = density_kg_ave * dist_gravity
          delta_pressure = boundary_pressure - &
                           auxvar_dn%pres(iphase) + &
                           gravity_term

#ifdef DEBUG_TOWG_FILEOUTPUT
          debug_dphi(iphase) = delta_pressure
#endif

          ! PO CONDUCTANCE_BC and SEEPAGE_BC to be implemented/tested
          if (bc_type == SEEPAGE_BC .or. &
              bc_type == CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                auxvar_up%pres(iphase) - &
                 option%reference_pressure < eps) then
              delta_pressure = 0.d0
            endif
          endif
          
          !upwinding mobility and enthalpy  
          if (delta_pressure >= 0.D0) then
            mobility = auxvar_up%mobility(iphase)
            uH = auxvar_up%H(iphase)
          else
            mobility = auxvar_dn%mobility(iphase)
            uH = auxvar_dn%H(iphase)
          endif      

          if (mobility > floweps) then
            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
            v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
            ! only need average density if velocity > 0.
            density_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                                  auxvar_dn%sat(iphase), &
                                                  auxvar_up%den(iphase), &
                                                  auxvar_dn%den(iphase) )
          endif
#ifndef BAD_MOVE1        
        endif ! sat > eps
#endif

      case(NEUMANN_BC)
        select case(option%phase_map(iphase))
          case(LIQUID_PHASE)
            idof = bc_auxvar_mapping(TOWG_LIQUID_FLUX_INDEX)
          case(OIL_PHASE)
            idof = bc_auxvar_mapping(TOWG_OIL_FLUX_INDEX)
          case(GAS_PHASE)
            idof = bc_auxvar_mapping(TOWG_GAS_FLUX_INDEX)
        end select
        
        neumann_bc_present = PETSC_TRUE
        !xmol = 0.d0
        !xmol(iphase) = 1.d0
        if (dabs(bc_auxvars(idof)) > floweps) then
          v_darcy(iphase) = bc_auxvars(idof)
          !upwinding based on given BC flux sign
          if (v_darcy(iphase) > 0.d0) then 
            density_ave = auxvar_up%den(iphase)
            uH = auxvar_up%H(iphase)
          else 
            density_ave = auxvar_dn%den(iphase)
            uH = auxvar_dn%H(iphase)
          endif 
        endif
      case default
        option%io_buffer = &
         'Boundary condition type not recognized in TOWGImsTLBCFlux phase loop'
        call printErrMsg(option)
    end select

    if (dabs(v_darcy(iphase)) > 0.d0) then
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area
      if (density_ave < 1.d-40) then
        option%io_buffer = 'Zero density in TOWGImsTLBCFlux()'
        call printErrMsgByRank(option)
      endif
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                              density_ave[kmol phase/m^3 phase]
      mole_flux = q*density_ave
      ! Res[kmol phase/sec] 
      Res(iphase) = mole_flux 

      !do icomp = 1, option%nflowspec
      !  ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
      !  !                      xmol[kmol comp/mol phase]
      !  Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      !enddo
#ifdef DEBUG_FLUXES  
      adv_flux(iphase) = mole_flux 
      !do icomp = 1, option%nflowspec
      !  adv_flux(icomp,iphase) = adv_flux(icomp,iphase) + mole_flux * xmol(icomp)
      !enddo
#endif
!#ifdef DEBUG_TOWG_FILEOUTPUT
!      do icomp = 1, option%nflowspec
!        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo
!#endif
      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH ! H_ave

#ifdef DEBUG_FLUXES  
      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
#endif
#ifdef DEBUG_TOWG_FILEOUTPUT
      debug_flux(iphase) = mole_flux * uH
#endif
    endif
  enddo
! CONVECTION
#endif

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'bc delta pressure :', debug_dphi(:)
    write(debug_unit,'(a,7es24.15)') 'bc adv flux (energy):', debug_flux(:)
  endif
  debug_flux = 0.d0
#endif                    

#ifdef CONDUCTION
  ! add heat conduction flux
  heat_flux = 0.d0
  select case (ibndtype(towg_energy_eq_idx))
    case (DIRICHLET_BC)
      ! based on Somerton et al., 1974:
      ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
      ! Assuming that oil and water have same conductivity:
      ! s_l = S_water + S_oil 
      sat_liquid = auxvar_dn%sat(option%liquid_phase) + &
                   auxvar_dn%sat(option%oil_phase)

      k_eff_dn = thermal_conductivity_dn(1) + &
                 sqrt(sat_liquid) * &
                 (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
      ! units:
      ! k_eff = W/K/m/m = J/s/K/m/m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = J/s
      k_eff_ave = k_eff_dn / dist(0)
      delta_temp = auxvar_up%temp - auxvar_dn%temp
      heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! convert W -> MW
    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = bc_auxvars(bc_auxvar_mapping(TOWG_ENERGY_FLUX_INDEX)) * area

    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'TOWGImsTLBCFlux heat conduction loop.'
      call printErrMsg(option)
  end select
  Res(energy_id) = Res(energy_id) + heat_flux ! MW
! CONDUCTION
#endif


#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
    write(*,'('' bc phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(1), auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(1), auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(1) * 1.d6
    write(*,'('' bc phase: oil'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(2), auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(2), auxvar_dn%sat(2)
    write(*,'(''  oil --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(2) * 1.d6
    write(*,'('' bc phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(3), auxvar_dn%pres(3)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(3), auxvar_dn%sat(3)
    write(*,'(''  gas --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(3)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(3) * 1.d6
    write(*,'(''  bc energy --'')')
    write(*,'(''   advective heat flux:'',es12.4)') adv_flux(4) * 1.d6
    write(*,'(''   conductive heat flux:'',es12.4)') heat_flux * 1.d6
    write(*,'(''   total heat flux:'',es12.4)') (heat_flux + adv_flux(4))*1.d6

  endif
#endif
  
end subroutine TOWGImsTLBCFlux

! ************************************************************************** !

subroutine TOWGImsTLSrcSink(option,src_sink_condition, auxvar, &
                            global_auxvar,ss_flow_vol_flux,scale,Res)
  ! 
  ! Computes the source/sink terms for the residual 
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  ! 

  use Option_module
  use Condition_module  

  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  type(flow_towg_condition_type), pointer :: src_sink_condition
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar !keep global_auxvar for salinity
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale  
  PetscReal :: Res(option%nflowdof)

  ! local parameter
  PetscInt, parameter :: SRC_TEMPERATURE = 1
  PetscInt, parameter :: SRC_ENTHALPY = 2 
  ! local variables
  PetscReal, pointer :: qsrc(:)
  PetscInt :: flow_src_sink_type    
  PetscReal :: qsrc_mol
  PetscReal :: den, den_kg, enthalpy, internal_energy_dummy, temperature
  PetscReal :: cell_pressure
  PetscInt :: iphase
  PetscInt :: energy_var
  PetscErrorCode :: ierr

  ! this can be removed if extending to pressure condition
  if (.not.associated(src_sink_condition%rate) ) then
    option%io_buffer = 'TOWGImsTLSrcSink fow condition rate not defined ' // &
    'rate is needed for a valid src/sink term'
    call printErrMsg(option)  
  end if

  qsrc => src_sink_condition%rate%dataset%rarray

  energy_var = 0
  if ( associated(src_sink_condition%temperature) ) then
    energy_var = SRC_TEMPERATURE 
  else if ( associated(src_sink_condition%enthalpy) ) then
    energy_var = SRC_ENTHALPY
  end if

  flow_src_sink_type = src_sink_condition%rate%itype

  ! checks that qsrc(liquid_phase), qsrc(oil_phase), qsrc(gas_phase) 
  ! do not have different signs
  ! if ( (qsrc(option%liquid_phase)>0.0d0 .and. qsrc(option%oil_phase)<0.d0).or.&
  !     (qsrc(option%liquid_phase)<0.0d0 .and. qsrc(option%oil_phase)>0.d0)  & 
  !   ) then
  !   option%io_buffer = "TOilImsSrcSink error: " // &
  !     "src(wat) and src(oil) with opposite sign"
  !   call printErrMsg(option)
  ! end if

  ! if not given, approximates BHP with pressure of perforated cell
  if ( associated(src_sink_condition%bhp_pressure) ) then
    cell_pressure = src_sink_condition%bhp_pressure%dataset%rarray(1)
  else
    cell_pressure = &
        maxval(auxvar%pres(option%liquid_phase:option%gas_phase))
  end if

  ! if enthalpy is used to define enthalpy or energy rate is used  
  ! approximate bottom hole temperature (BHT) with local temp
  if ( associated(src_sink_condition%temperature) ) then
    temperature = src_sink_condition%temperature%dataset%rarray(1)
  else   
    temperature = auxvar%temp
  end if

  Res = 0.d0
  do iphase = 1, option%nphase
    qsrc_mol = 0.d0
    if ( qsrc(iphase) > 0.d0) then 
      select case(option%phase_map(iphase))
        case(LIQUID_PHASE)
          call EOSWaterDensity(temperature,cell_pressure,den_kg,den,ierr)
        case(OIL_PHASE)
          call EOSOilDensity(temperature,cell_pressure,den,ierr)
        case(GAS_PHASE)
          call EOSGasDensity(temperature,cell_pressure,den,ierr)
      end select 
    else
      den = auxvar%den(iphase)
    end if

    select case(flow_src_sink_type)
      ! injection and production 
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(iphase)/towg_fmw_comp(iphase) ! kg/sec -> kmol/sec
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_mol = qsrc(iphase)/towg_fmw_comp(iphase)*scale 
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now 
                  ! qsrc(iphase) = m^3/sec  
        qsrc_mol = qsrc(iphase)*den ! den = kmol/m^3 
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_mol = qsrc(iphase)* den * scale
    end select
    ss_flow_vol_flux(iphase) = qsrc_mol/ den
    Res(iphase) = qsrc_mol
  enddo

  ! when using scaled src/sinks, the rates (mass or vol) scaling 
  ! at this point the scale factor is already included in Res(iphase)

  ! Res(option%energy_id), energy units: MJ/sec

  !if ( associated(src_sink_condition%temperature) .or. &
  !    associated(src_sink_condition%enthalpy) &
  !   ) then
  !if the energy rate is not given, use either temperature or enthalpy
  if ( dabs(qsrc(FOUR_INTEGER)) < 1.0d-40 ) then
    ! water injection 
    if (qsrc(option%liquid_phase) > 0.d0) then !implies qsrc(option%oil_phase)>=0
      if ( energy_var == SRC_ENTHALPY ) then
        !input as J/kg
        enthalpy = src_sink_condition%enthalpy% &
                       dataset%rarray(option%liquid_phase)
                     ! J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * towg_fmw_comp(option%liquid_phase)
      else !note: temp can either be input or taken as the one of perf. block
      !else if ( energy_var == SRC_TEMPERATURE ) then
        call EOSWaterEnthalpy(temperature, cell_pressure,enthalpy,ierr)
        ! enthalpy = [J/kmol]
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! water component mass    
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * enthalpy
    end if
    ! oil injection 
    if (qsrc(option%oil_phase) > 0.d0) then !implies qsrc(option%liquid_phase)>=0
      if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%oil_phase)
                      !J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * towg_fmw_comp(option%oil_phase)
      else !note: temp can either be input or taken as the one of perf. block
      !if ( energy_var == SRC_TEMPERATURE ) then
        call EOSOilEnthalpy(temperature,cell_pressure,enthalpy,ierr)
        ! enthalpy = [J/kmol]
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! oil component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * enthalpy
    end if
    ! gas injection 
    if (qsrc(option%gas_phase) > 0.d0) then
      if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%gas_phase)
                      !J/kg * kg/kmol = J/kmol
        enthalpy = enthalpy * towg_fmw_comp(option%gas_phase)
      else !note: temp can either be input or taken as the one of perf. block
      !if ( energy_var == SRC_TEMPERATURE ) then
        call EOSGasEnergy(temperature,cell_pressure,enthalpy, &
                              internal_energy_dummy,ierr)
        ! enthalpy = [J/kmol]
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! oil component mass           
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%gas_phase) * enthalpy
    end if
    ! water energy extraction due to water production
    if (qsrc(option%liquid_phase) < 0.d0) then !implies qsrc(option%oil_phase)<=0
      ! auxvar enthalpy units: MJ/kmol ! water component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * &
                              auxvar%H(option%liquid_phase)
    end if
    !oil energy extraction due to oil production
    if (qsrc(option%oil_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol ! water component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * &
                              auxvar%H(option%oil_phase)
    end if
    if (qsrc(option%gas_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol ! water component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%gas_phase) * &
                              auxvar%H(option%gas_phase)
    end if
  else !if the energy rate is given, it overwrites both temp and enthalpy
    ! if energy rate is given, loaded in qsrc(4) in MJ/sec 
    Res(option%energy_id) = qsrc(FOUR_INTEGER)* scale ! MJ/s
  end if

  nullify(qsrc)      
  
end subroutine TOWGImsTLSrcSink

! ************************************************************************** !

subroutine TOWGAccumDerivative(auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  ! 

  use Option_module
  use Saturation_Function_module
  use Material_Aux_class
  
  implicit none

  type(auxvar_towg_type) :: auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: jac(option%nflowdof,option%nflowdof)
  PetscReal :: jac_pert(option%nflowdof,option%nflowdof)
  PetscInt :: idof, irow

  call TOWGAccumulation(auxvar(ZERO_INTEGER),global_auxvar,material_auxvar, &
                        soil_heat_capacity,option,res,PETSC_FALSE)

  do idof = 1, option%nflowdof
    call TOWGAccumulation(auxvar(idof),global_auxvar,material_auxvar, &
                          soil_heat_capacity,option,res_pert,PETSC_FALSE)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/auxvar(idof)%pert
    enddo !irow
  enddo ! idof

  if (towg_isothermal) then
    J(towg_energy_eq_idx,:) = 0.d0
    J(:,towg_energy_eq_idx) = 0.d0
  endif

  if (towg_no_oil) then
    J(TOWG_OIL_EQ_IDX,:) = 0.d0
    J(:,TOWG_OIL_EQ_IDX) = 0.d0
  endif

  if (towg_no_gas) then
    J(TOWG_GAS_EQ_IDX,:) = 0.d0
    J(:,TOWG_GAS_EQ_IDX) = 0.d0
  endif
  

end subroutine TOWGAccumDerivative

! ************************************************************************** !

subroutine TOWGFluxDerivative(auxvar_up,global_auxvar_up, &
                              material_auxvar_up, &
                              sir_up, &
                              thermal_conductivity_up, &
                              auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              sir_dn, &
                              thermal_conductivity_dn, &
                              area, dist, &
                              towg_parameter, &
                              option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  ! 
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(auxvar_towg_type) :: auxvar_up(0:), auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0
  
  !print *, 'TOWGFluxDerivative'
  option%iflag = -2
  call TOWGFlux(auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                material_auxvar_up,sir_up, &
                thermal_conductivity_up, &
                auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                material_auxvar_dn,sir_dn, &
                thermal_conductivity_dn, &
                area,dist,towg_parameter, &
                option,v_darcy,res,PETSC_FALSE)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call TOWGFlux(auxvar_up(idof),global_auxvar_up, &
                  material_auxvar_up,sir_up, &
                  thermal_conductivity_up, &
                  auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn,sir_dn, &
                  thermal_conductivity_dn, &
                  area,dist,towg_parameter, &
                  option,v_darcy,res_pert,PETSC_FALSE)
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res(irow))/auxvar_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call TOWGFlux(auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up,sir_up, &
                  thermal_conductivity_up, &
                  auxvar_dn(idof),global_auxvar_dn, &
                  material_auxvar_dn,sir_dn, &
                  thermal_conductivity_dn, &
                  area,dist,towg_parameter, &
                  option,v_darcy,res_pert,PETSC_FALSE)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (towg_isothermal) then
    Jup(towg_energy_eq_idx,:) = 0.d0
    Jup(:,towg_energy_eq_idx) = 0.d0
    Jdn(towg_energy_eq_idx,:) = 0.d0
    Jdn(:,towg_energy_eq_idx) = 0.d0
  endif
  
  if (towg_no_oil) then
    Jup(TOWG_OIL_EQ_IDX,:) = 0.d0
    Jup(:,TOWG_OIL_EQ_IDX) = 0.d0
    Jdn(TOWG_OIL_EQ_IDX,:) = 0.d0
    Jdn(:,TOWG_OIL_EQ_IDX) = 0.d0
  endif  

  if (towg_no_gas) then
    Jup(TOWG_GAS_EQ_IDX,:) = 0.d0
    Jup(:,TOWG_GAS_EQ_IDX) = 0.d0
    Jdn(TOWG_GAS_EQ_IDX,:) = 0.d0
    Jdn(:,TOWG_GAS_EQ_IDX) = 0.d0
  endif  

end subroutine TOWGFluxDerivative

! ************************************************************************** !

subroutine TOWGBCFluxDerivative(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                                auxvar_up, global_auxvar_up, &
                                auxvar_dn,global_auxvar_dn, &
                                material_auxvar_dn,sir_dn, &
                                thermal_conductivity_dn, &
                                area,dist,towg_parameter, &
                                option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  ! 

  use Option_module 
  use Material_Aux_class
  
  implicit none

  PetscReal :: bc_auxvars(:) ! from aux_real_var array
  type(auxvar_towg_type) :: auxvar_up, auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: bc_auxvar_mapping(TOWG_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jdn = 0.d0
  !print *, 'TOWGBCFluxDerivative'

  option%iflag = -2
  call TOWGBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                  auxvar_up,global_auxvar_up, &
                  auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn,sir_dn, &
                  thermal_conductivity_dn, &
                  area,dist,towg_parameter, &
                  option,v_darcy,res,PETSC_FALSE)
                    
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call TOWGBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                    auxvar_up,global_auxvar_up, &
                    auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn,sir_dn, &
                    thermal_conductivity_dn, &
                    area,dist,towg_parameter, &
                    option,v_darcy,res_pert,PETSC_FALSE)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (towg_isothermal) then
    Jdn(towg_energy_eq_idx,:) = 0.d0
    Jdn(:,towg_energy_eq_idx) = 0.d0
  endif
  
  if (towg_no_oil) then
    Jdn(TOWG_OIL_EQ_IDX,:) = 0.d0
    Jdn(:,TOWG_OIL_EQ_IDX) = 0.d0
  endif  

  if (towg_no_gas) then
    Jdn(TOWG_GAS_EQ_IDX,:) = 0.d0
    Jdn(:,TOWG_GAS_EQ_IDX) = 0.d0
  endif  
  
end subroutine TOWGBCFluxDerivative

! ************************************************************************** !

subroutine TOWGSrcSinkDerivative(option,src_sink_condition,auxvars, &
                                 global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/28/16
  ! 

  use Option_module
  use Condition_module

  implicit none

  type(option_type) :: option
  type(flow_towg_condition_type), pointer :: src_sink_condition
  type(auxvar_towg_type) :: auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  option%iflag = -3
  call TOWGSrcSink(option,src_sink_condition,auxvars(ZERO_INTEGER), &
                   global_auxvar,dummy_real,scale,Res)

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call TOWGSrcSink(option,src_sink_condition,auxvars(idof), &
                     global_auxvar,dummy_real,scale,res_pert)
    do irow = 1, option%nflowdof
      Jac(irow,idof) = (res_pert(irow)-res(irow))/auxvars(idof)%pert
    enddo !irow
  enddo ! idof
 
  if (towg_isothermal) then
    Jac(towg_energy_eq_idx,:) = 0.d0
    Jac(:,towg_energy_eq_idx) = 0.d0
  endif
  
  if (towg_no_oil) then
    Jac(TOWG_OIL_EQ_IDX,:) = 0.d0
    Jac(:,TOWG_OIL_EQ_IDX) = 0.d0
  endif  

  if (towg_no_gas) then
    Jac(TOWG_GAS_EQ_IDX,:) = 0.d0
    Jac(:,TOWG_GAS_EQ_IDX) = 0.d0
  endif  

end subroutine TOWGSrcSinkDerivative

! ************************************************************************** !

subroutine TOWGResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/28/16
  ! 
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscmat.h"
  use petscsnes
  use petscmat
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module

  use Connection_module
  use Grid_module
  use Coupler_module  
  use Debug_module
  use Material_Aux_class

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  Mat, parameter :: null_mat = PETSC_NULL_MAT
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  class(pm_towg_aux_type), pointer :: towg
  type(towg_parameter_type), pointer :: towg_parameter
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: iphase
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)
  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn
  PetscInt, save :: iplot = 0

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word

  PetscInt :: icap_up, icap_dn
  PetscReal :: Res(realization%option%nflowdof)
  !PetscReal :: Jac_dummy(realization%option%nflowdof, &
  !                       realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)
  
  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  towg => patch%aux%TOWG
  towg_parameter => towg%parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  
#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    debug_iteration_count = debug_iteration_count + 1
    write(word,*) debug_timestep_count
    string = 'residual_debug_data_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    debug_filename = trim(string) // '_' // trim(adjustl(word)) // '.txt'
    open(debug_unit, file=debug_filename, action="write", status="unknown")
    open(debug_info_unit, file='debug_info.txt', action="write", &
         position="append", status="unknown")
    write(debug_info_unit,*) 'residual ', debug_timestep_count, &
      debug_timestep_cut_count, debug_iteration_count
    close(debug_info_unit)
  endif
#endif

  ! Communication -----------------------------------------
  ! These 3 must be called before TOWGUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  
                                     ! do update state
  call TOWGUpdateAuxVars(realization,PETSC_TRUE)

  ! override flags since they will soon be out of date
  towg%auxvars_up_to_date = PETSC_FALSE 

  ! always assume variables have been swapped; therefore, must copy back
  call VecLockPop(xx,ierr); CHKERRQ(ierr)
  call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                   NFLOWDOF)
  call VecLockPush(xx,ierr); CHKERRQ(ierr)

  if (option%compute_mass_balance_new) then
    call TOWGZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  r_p = -accum_p

  
  !Heeho dynamically update p+1 accumulation term
  !if (towg_tough2_conv_criteria) then
  !  call VecGetArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif
  
  ! accumulation at t(k+1)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call TOWGAccumulation(towg%auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          material_parameter%soil_heat_capacity(imat), &
                          option,Res,local_id == towg_debug_cell_id)


    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    
    !Heeho dynamically update p+1 accumulation term
    !if (towg_tough2_conv_criteria) then
    !  accum_p2(local_start:local_end) = Res(:)
    !endif
    
  enddo

  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  !Heeho dynamically update p+1 accumulation term
  !if (towg_tough2_conv_criteria) then
  !  call VecRestoreArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif

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

      imat_up = patch%imat(ghosted_id_up) 
      imat_dn = patch%imat(ghosted_id_dn) 
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)

      call TOWGFlux(towg%auxvars(ZERO_INTEGER,ghosted_id_up), &
                    global_auxvars(ghosted_id_up), &
                    material_auxvars(ghosted_id_up), & 
                    material_parameter%soil_residual_saturation(:,icap_up), &
                    material_parameter%soil_thermal_conductivity(:,imat_up), &
                    towg%auxvars(ZERO_INTEGER,ghosted_id_dn), &
                    global_auxvars(ghosted_id_dn), &
                    material_auxvars(ghosted_id_dn), &
                    material_parameter%soil_residual_saturation(:,icap_dn), &
                    material_parameter%soil_thermal_conductivity(:,imat_dn), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    towg_parameter,option,v_darcy,Res, &
                    (local_id_up == towg_debug_cell_id .or. &
                     local_id_dn == towg_debug_cell_id))

      patch%internal_velocities(:,sum_connection) = v_darcy
      if (associated(patch%internal_flow_fluxes)) then
        patch%internal_flow_fluxes(:,sum_connection) = Res(:)
      endif
      
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
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id)

      call TOWGBCFlux(boundary_condition%flow_bc_type, &
                    boundary_condition%flow_aux_mapping, & 
                    boundary_condition%flow_aux_real_var(:,iconn), &
                    towg%auxvars_bc(sum_connection), &
                    global_auxvars_bc(sum_connection), &
                    towg%auxvars(ZERO_INTEGER,ghosted_id), &
                    global_auxvars(ghosted_id), &
                    material_auxvars(ghosted_id), &
                    material_parameter%soil_residual_saturation(:,icap_dn), &
                    material_parameter%soil_thermal_conductivity(:,imat_dn), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    towg_parameter,option,v_darcy,Res, &
                    local_id == towg_debug_cell_id)

      patch%boundary_velocities(:,sum_connection) = v_darcy
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      r_p(local_start:local_end)= r_p(local_start:local_end) - Res(:)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
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

      call TOWGSrcSink(option,source_sink%flow_condition%towg, &
                       towg%auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), ss_flow_vol_flux, &
                       scale,Res)

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif      
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif      
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

    enddo
    source_sink => source_sink%next
  enddo

  if (towg%inactive_cells_exist) then
    do i=1,towg%n_inactive_rows
      r_p(towg%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  if (Initialized(towg_debug_cell_id)) then
    call VecGetArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
    do local_id = towg_debug_cell_id-1, towg_debug_cell_id+1
      write(*,'(''  residual   : '',i2,10es12.4)') local_id, &
        r_p((local_id-1)*option%nflowdof+1:(local_id-1)*option%nflowdof+2), &
        r_p(local_id*option%nflowdof)*1.d6
    enddo
    call VecRestoreArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  
  if (towg_isothermal) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+towg_energy_eq_idx) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  if (towg_no_oil) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+TOWG_OIL_EQ_IDX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif  
  if (towg_no_gas) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+TOWG_GAS_EQ_IDX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif  

#ifdef DEBUG_TOWG_FILEOUTPUT
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'fixed residual:', local_id, &
      accum_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'residual:', local_id, &
      r_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
#endif
  
  if (realization%debug%vecview_residual) then
    string = 'TOWGresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'TOWGxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    close(debug_unit)
  endif
#endif
  
end subroutine TOWGResidual

! ************************************************************************** !

subroutine TOWGJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/28/16
  ! 

#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscmat.h"
  use petscsnes
  use petscmat
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscInt :: icap_up,icap_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: irow
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = PETSC_NULL_VEC
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  PetscInt, pointer :: zeros(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(material_parameter_type), pointer :: material_parameter
  class(pm_towg_aux_type), pointer :: towg
  type(towg_parameter_type), pointer :: towg_parameter
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  towg => patch%aux%TOWG
  towg_parameter => towg%parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    write(word,*) debug_timestep_count
    string = 'jacobian_debug_data_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    debug_filename = trim(string) // '_' // trim(adjustl(word)) // '.txt'
    open(debug_unit, file=debug_filename, action="write", status="unknown")
    open(debug_info_unit, file='debug_info.txt', action="write", &
         position="append", status="unknown")
    write(debug_info_unit,*) 'jacobian ', debug_timestep_count, &
      debug_timestep_cut_count, debug_iteration_count
    close(debug_info_unit)
  endif
#endif

  ! Perturb aux vars
  do ghosted_id = 1, grid%ngmax  ! For each local node do...
    if (patch%imat(ghosted_id) <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    call TOWGAuxVarPerturb(towg%auxvars(:,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           patch%characteristic_curves_array( &
                           patch%sat_func_id(ghosted_id))%ptr, &
                           natural_id,option)
  enddo

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call TOWGAccumDerivative(towg%auxvars(:,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option, &
                             Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_accum'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
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

      imat_up = patch%imat(ghosted_id_up)
      imat_dn = patch%imat(ghosted_id_dn)
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
                              
      call TOWGFluxDerivative(towg%auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     material_parameter%soil_residual_saturation(:,icap_up), &
                     material_parameter%soil_thermal_conductivity(:,imat_up), &
                     towg%auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     material_parameter%soil_residual_saturation(:,icap_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     towg_parameter,option,&
                     Jup,Jdn)

      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_flux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id)

      call TOWGBCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      towg%auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      towg%auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      material_parameter%soil_residual_saturation(:,icap_dn), &
                      material_parameter%soil_thermal_conductivity(:,imat_dn), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      towg_parameter,option, &
                      Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_bcflux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  ! Source/sinks
  source_sink => patch%source_sink_list%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      Jup = 0.d0
      call TOWGSrcSinkDerivative(option, &
                        source_sink%flow_condition%towg, &
                        towg%auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo
  
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_srcsink'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out isothermal and inactive cells
  if (towg%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,towg%n_inactive_rows, &
                          towg%inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif

  if (towg_isothermal) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => towg%row_zeroing_array
    ! zero energy residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        towg_energy_eq_idx - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_VEC, &
                          PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  endif

  if (towg_no_oil) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => towg%row_zeroing_array
    ! zero gas component mass balance residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        TOWG_OIL_EQ_IDX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_VEC, &
                          PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  endif

  if (towg_no_gas) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => towg%row_zeroing_array
    ! zero gas component mass balance residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        TOWG_GAS_EQ_IDX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_VEC, &
                          PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  endif
  
  if (realization%debug%matview_Jacobian) then
    string = 'TOWGjacobian'
    call DebugCreateViewer(realization%debug,string,option,viewer)
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

!  call MatView(J,PETSC_VIEWER_STDOUT_WORLD,ierr)

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    write(word,*) debug_timestep_count
    string = 'jacobian_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    string = trim(string) // '_' // trim(adjustl(word)) // '.out'
    call PetscViewerASCIIOpen(realization%option%mycomm,trim(string), &
                              viewer,ierr);CHKERRQ(ierr)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    close(debug_unit)
  endif
#endif

end subroutine TOWGJacobian

! ************************************************************************** !

subroutine TOWGImsTLCheckUpdatePre(line_search,X,dX,changed,realization, &
                                   max_it_before_damping,damping_factor, &
                                   max_pressure_change,ierr)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/30/16
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module

  implicit none

  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  type(realization_subsurface_type) :: realization
  PetscInt :: max_it_before_damping
  PetscReal :: damping_factor
  PetscReal :: max_pressure_change

  PetscReal, pointer :: X_p(:), dX_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field

  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)
  !type(global_auxvar_type), pointer :: global_auxvars(:)  

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset

  PetscInt :: pressure_index, saturation_index, temperature_index

  PetscReal :: pressure0, pressure1, del_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation

  PetscReal :: max_saturation_change = 0.125d0
  PetscReal :: max_temperature_change = 10.d0
  PetscReal :: scale, temp_scale, temp_real
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  SNES :: snes
  PetscInt :: newton_iteration

  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  patch => realization%patch

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE

  !print *, "TOWGImsTLCheckUpdatePre"
  ! truncation
  ! Oil and Gas Saturations must be truncated.  We do not use scaling
  ! here because of the very small values.  just truncation.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*option%nflowdof
    saturation_index = offset + TOWG_OIL_SATURATION_DOF
    if ( (X_p(saturation_index) - dX_p(saturation_index)) < 0.d0 ) then
      ! we use 1.d-6 since cancelation can occur with smaller values
      ! this threshold is imposed in the initial condition
      dX_p(saturation_index) = X_p(saturation_index)
    end if
    saturation_index = offset + TOWG_GAS_SATURATION_3PH_DOF
    if ( (X_p(saturation_index) - dX_p(saturation_index)) < 0.d0 ) then
      ! we use 1.d-6 since cancelation can occur with smaller values
      ! this threshold is imposed in the initial condition
      dX_p(saturation_index) = X_p(saturation_index)
    end if
  enddo

  scale = initial_scale
  if (max_it_before_damping > 0 .and. &
      newton_iteration > max_it_before_damping) then
    scale = damping_factor
  endif

#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
!!#define LIMIT_MAX_TEMPERATURE_CHANGE
!! TRUNCATE_PRESSURE is needed for times when the solve wants
!! to pull them negative.
!#define TRUNCATE_PRESSURE

  ! scaling
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    offset = (local_id-1)*option%nflowdof
    temp_scale = 1.d0
    pressure_index = offset + TOWG_OIL_PRESSURE_DOF
    dX_p(pressure_index) = dX_p(pressure_index) * towg_pressure_scale
    del_pressure = dX_p(pressure_index)
    pressure0 = X_p(pressure_index)
    pressure1 = pressure0 - del_pressure
#ifdef LIMIT_MAX_PRESSURE_CHANGE
    if (dabs(del_pressure) > max_pressure_change) then
      temp_real = dabs(max_pressure_change/del_pressure)
      temp_scale = min(temp_scale,temp_real)
     endif
#endif
#ifdef TRUNCATE_PRESSURE
    if (pressure1 <= 0.d0) then
      if (dabs(del_pressure) > 1.d-40) then
        temp_real = tolerance * dabs(pressure0 / del_pressure)
        temp_scale = min(temp_scale,temp_real)
      endif
    endif
#endif 
!TRUNCATE_PRESSURE
#ifdef LIMIT_MAX_SATURATION_CHANGE
    !oil saturation
    saturation_index = offset + TOWG_OIL_SATURATION_DOF
    del_saturation = dX_p(saturation_index)
    !saturation0 = X_p(saturation_index)
    !saturation1 = saturation0 - del_saturation
    if (dabs(del_saturation) > max_saturation_change) then
       temp_real = dabs(max_saturation_change/del_saturation)
       temp_scale = min(temp_scale,temp_real)
    endif
    !gas saturation
    saturation_index = offset + TOWG_GAS_SATURATION_3PH_DOF
    del_saturation = dX_p(saturation_index)
    !saturation0 = X_p(saturation_index)
    !saturation1 = saturation0 - del_saturation
    if (dabs(del_saturation) > max_saturation_change) then
       temp_real = dabs(max_saturation_change/del_saturation)
       temp_scale = min(temp_scale,temp_real)
    endif
#endif 
!LIMIT_MAX_SATURATION_CHANGE
#ifdef LIMIT_MAX_TEMPERATURE_CHANGE        
    temperature_index  = offset + towg_energy_dof
    del_temperature = dX_p(temperature_index)
    if (dabs(del_temperature) > max_temperature_change) then
       temp_real = dabs(max_temperature_change/del_temperature)
       temp_scale = min(temp_scale,temp_real)
    endif
#endif 
!LIMIT_MAX_TEMPERATURE_CHANGE
    scale = min(scale,temp_scale) 
  enddo

  temp_scale = scale
  call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MIN,option%mycomm,ierr)

  ! it performs an homogenous scaling using the smallest scaling factor
  ! over all subdomains domains
  if (scale < 0.9999d0) then
    dX_p = scale*dX_p
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine TOWGImsTLCheckUpdatePre

! ************************************************************************** !

function TOWGImsTLAverageDensity(sat_up,sat_dn,density_up,density_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/18/16
  ! 

  implicit none

  PetscReal :: sat_up, sat_dn
  PetscReal :: density_up, density_dn

  PetscReal :: TOWGImsTLAverageDensity

  if (sat_up < eps ) then
    TOWGImsTLAverageDensity = density_dn
  else if (sat_dn < eps ) then 
    TOWGImsTLAverageDensity = density_up
  else ! in here we could use an armonic average, 
       ! other idea sat weighted average but it needs truncation
    TOWGImsTLAverageDensity = 0.5d0*(density_up+density_dn)
  end if

end function TOWGImsTLAverageDensity

! ************************************************************************** !

subroutine TOWGDestroy(realization)
  ! 
  ! Deallocates variables associated with TOWG
  ! 
  ! Author: Paolo Orsini
  ! Date: 01/19/17
  ! 

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! place anything that needs to be freed here.

end subroutine TOWGDestroy

! ************************************************************************** !

end module TOWG_module
