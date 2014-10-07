module General_module

  use General_Aux_module
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

#define CONVECTION
#define DIFFUSION
#define LIQUID_DIFFUSION
#define CONDUCTION
  
!#define DEBUG_GENERAL_FILEOUTPUT
!#define DEBUG_FLUXES  

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

#ifdef DEBUG_GENERAL_FILEOUTPUT
  PetscInt, parameter :: debug_unit = 87
  PetscInt, parameter :: debug_info_unit = 86
  character(len=MAXWORDLENGTH) :: debug_filename
  PetscInt :: debug_flag = 0
  PetscInt :: debug_iteration_count
  PetscInt :: debug_timestep_cut_count
  PetscInt :: debug_timestep_count
#endif

  public GeneralRead, &
         GeneralSetup, &
         GeneralInitializeTimestep, &
         GeneralUpdateSolution, &
         GeneralTimeCut,&
         GeneralUpdateAuxVars, &
         GeneralUpdateFixedAccum, &
         GeneralComputeMassBalance, &
         GeneralResidual, &
         GeneralJacobian, &
         GeneralGetTecplotHeader, &
         GeneralSetPlotVariables, &
         GeneralCheckUpdatePre, &
         GeneralCheckUpdatePost, &
         GeneralMapBCAuxvarsToGlobal, &
         GeneralDestroy

contains

! ************************************************************************** !

subroutine GeneralRead(input,option)
  ! 
  ! Reads parameters for general phase
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  call InputReadWord(input,option,keyword,PETSC_TRUE)
  if (input%ierr /= 0) then
    return
  endif
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GENERAL_MODE')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('WINDOW_EPSILON') 
        call InputReadDouble(input,option,window_epsilon)
        call InputErrorMsg(input,option,'diffusion coefficient','GENERAL_MODE')
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        !geh: assuming gas component is index 2
        call InputReadDouble(input,option,fmw_comp(2))
        call InputErrorMsg(input,option,'gas component formula wt.', &
                           'GENERAL_MODE')
      case('TWO_PHASE_ENERGY_DOF')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'two_phase_energy_dof','GENERAL_MODE')
        call GeneralAuxSetEnergyDOF(word,option)
      case('ISOTHERMAL')
        general_isothermal = PETSC_TRUE
      case('NO_AIR')
        general_no_air = PETSC_TRUE
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,general_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           'GENERAL_MODE')
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,general_max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           'GENERAL_MODE')
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,general_damping_factor)
        call InputErrorMsg(input,option,'damping factor','GENERAL_MODE')
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in General Mode'    
        call printErrMsg(option)
    end select
    
  enddo  
  
  if (general_isothermal .and. &
      general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
    option%io_buffer = 'Isothermal GENERAL mode may only be run with ' // &
                       'temperature as the two phase energy dof.'
    call printErrMsg(option)
  endif

end subroutine GeneralRead

! ************************************************************************** !

subroutine GeneralSetup(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Fluid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: i, idof, count
  PetscBool :: error_found
  PetscInt :: flag(10)
                                                ! extra index for derivatives
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(general_auxvar_type), pointer :: gen_auxvars_bc(:)
  type(general_auxvar_type), pointer :: gen_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%General => GeneralAuxCreate(option)

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
    option%io_buffer = 'Material property errors found in GeneralSetup.'
    call printErrMsg(option)
  endif
  
  ! allocate auxvar data structures for all grid cells  
  allocate(gen_auxvars(0:option%nflowdof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      call GeneralAuxVarInit(gen_auxvars(idof,ghosted_id),option)
    enddo
  enddo
  patch%aux%General%auxvars => gen_auxvars
  patch%aux%General%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them 
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  if (sum_connection > 0) then
    allocate(gen_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_auxvars_bc(iconn),option)
    enddo
    patch%aux%General%auxvars_bc => gen_auxvars_bc
  endif
  patch%aux%General%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    allocate(gen_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_auxvars_ss(iconn),option)
    enddo
    patch%aux%General%auxvars_ss => gen_auxvars_ss
  endif
  patch%aux%General%num_aux_ss = sum_connection
    
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call GeneralCreateZeroArray(patch,option)

  ! create array for zeroing Jacobian entries if isothermal and/or no air
  allocate(patch%aux%General%row_zeroing_array(grid%nlmax))
  patch%aux%General%row_zeroing_array = 0
  
  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    patch%aux%General%general_parameter% &
      diffusion_coefficient(cur_fluid_property%phase_id) = &
        cur_fluid_property%diffusion_coefficient
    cur_fluid_property => cur_fluid_property%next
  enddo  

  if (realization%output_option%print_fluxes) then
    allocate(patch%internal_fluxes(option%nflowdof,1, &
           ConnectionGetNumberInList(patch%grid%internal_connection_set_list)))
    patch%internal_fluxes = 0.d0
  endif

  call GeneralSetPlotVariables(realization) 
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flag = 0
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = 0
  ! create new file
  open(debug_info_unit, file='debug_info.txt', action="write", &
       status="unknown")
  write(debug_info_unit,*) 'type timestep cut iteration'
  close(debug_info_unit)
#endif  

end subroutine GeneralSetup

! ************************************************************************** !

subroutine GeneralInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_class
  
  implicit none
  
  type(realization_type) :: realization

  call GeneralUpdateFixedAccum(realization)
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flag = 0
!  if (realization%option%time >= 35.6d0*3600d0*24.d0*365.d0 - 1.d-40) then
!  if (.false.) then
  if (.true.) then
    debug_iteration_count = 0
    debug_flag = 1
  endif
  debug_iteration_count = 0
#endif

end subroutine GeneralInitializeTimestep

! ************************************************************************** !

subroutine GeneralUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal, pointer :: iphas_loc_p(:)
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  gen_auxvars => patch%aux%General%auxvars  
  global_auxvars => patch%aux%Global%auxvars
  
  call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)

  if (realization%option%compute_mass_balance_new) then
    call GeneralUpdateMassBalance(realization)
  endif
  
  ! update stored state
  call VecGetArrayF90(field%iphas_loc,iphas_loc_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    iphas_loc_p(ghosted_id) = global_auxvars(ghosted_id)%istate
    gen_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS) = &
      global_auxvars(ghosted_id)%istate
  enddo
  call VecRestoreArrayF90(field%iphas_loc,iphas_loc_p,ierr);CHKERRQ(ierr)
  
  ! update ghosted iphas_loc values (must come after 
  ! GeneralUpdateSolutionPatch)
  call DiscretizationLocalToLocal(realization%discretization, &
                                  field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)
  
  ! Set states of ghosted cells
  call VecGetArrayF90(field%iphas_loc,iphas_loc_p,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, realization%patch%grid%ngmax
    realization%patch%aux%Global%auxvars(ghosted_id)%istate = &
      int(iphas_loc_p(ghosted_id))
  enddo
  call VecRestoreArrayF90(field%iphas_loc,iphas_loc_p,ierr);CHKERRQ(ierr)

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = debug_timestep_count + 1
#endif   
  
end subroutine GeneralUpdateSolution

! ************************************************************************** !

subroutine GeneralTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 
  use Realization_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  
  PetscInt :: local_id, ghosted_id
  PetscReal, pointer :: iphas_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars

  call VecCopy(field%flow_yy,field%flow_xx,ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(realization%discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)
  
  ! restore stored state
  call VecGetArrayReadF90(field%iphas_loc,iphas_loc_p, ierr);CHKERRQ(ierr)
  do ghosted_id = 1, grid%ngmax
    global_auxvars(ghosted_id)%istate = int(iphas_loc_p(ghosted_id))
  enddo
  call VecRestoreArrayReadF90(field%iphas_loc,iphas_loc_p, ierr);CHKERRQ(ierr)

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_timestep_cut_count = debug_timestep_cut_count + 1
#endif 

  call GeneralInitializeTimestep(realization)  

end subroutine GeneralTimeCut

! ************************************************************************** !

subroutine GeneralComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 
 
  use Realization_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(general_auxvar_type), pointer :: general_auxvars(:,:)
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

  general_auxvars => patch%aux%General%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = &
        general_auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        material_auxvars(ghosted_id)%porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density
      do icomp = 1, option%nflowspec
        mass_balance(icomp,iphase) = mass_balance(icomp,iphase) + &
          general_auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)* &
          general_auxvars(ZERO_INTEGER,ghosted_id)%xmol(icomp,iphase) * &
          fmw_comp(icomp)*vol_phase
      enddo
    enddo
  enddo

end subroutine GeneralComputeMassBalance

! ************************************************************************** !

subroutine GeneralZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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

  do iconn = 1, patch%aux%General%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%General%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine GeneralZeroMassBalanceDelta

! ************************************************************************** !

subroutine GeneralUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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
  PetscInt :: icomp

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%General%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%General%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine GeneralUpdateMassBalance

! ************************************************************************** !

subroutine GeneralUpdateAuxVars(realization,update_state)
  ! 
  ! Updates the auxiliary variables associated with the General problem
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_class
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

  type(realization_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(general_auxvar_type), pointer :: gen_auxvars(:,:), gen_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: ghosted_start, ghosted_end
  PetscInt :: iphasebc, iphase
  PetscInt :: offset
  PetscInt :: istate
  PetscReal :: gas_pressure, capillary_pressure, liquid_saturation
  PetscReal :: saturation_pressure, temperature
  PetscInt :: real_index, variable
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  gen_auxvars => patch%aux%General%auxvars
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
     if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    ! flag(1) indicates call from non-perturbation
    option%iflag = 1
    call GeneralAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       gen_auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%characteristic_curves_array( &
                         patch%sat_func_id(ghosted_id))%ptr, &
                       ghosted_id, &
                       option)
    if (update_state) then
      call GeneralAuxVarUpdateState(xx_loc_p(ghosted_start:ghosted_end), &
                                    gen_auxvars(ZERO_INTEGER,ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%sat_func_id(ghosted_id))%ptr, &
                                    ghosted_id, &  ! for debugging
                                    option)
    endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,i5,i3,7es24.15)') 'auxvar:', ghosted_id, &
                        global_auxvars(ghosted_id)%istate, &
                        xx_loc_p(ghosted_start:ghosted_end)
  endif
#endif
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
      offset = (ghosted_id-1)*option%nflowdof
      if (patch%imat(ghosted_id) <= 0) cycle

      xxbc(:) = xx_loc_p(offset+1:offset+option%nflowdof)
      istate = boundary_condition%flow_aux_int_var(GENERAL_STATE_INDEX,iconn)
      if (istate == ANY_STATE) then
        istate = global_auxvars(ghosted_id)%istate
        select case(istate)
          case(LIQUID_STATE,GAS_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(DIRICHLET_BC,HYDROSTATIC_BC)
                  real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
              end select   
            enddo
          case(TWO_PHASE_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(HYDROSTATIC_BC)
                  real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                case(DIRICHLET_BC)
                  variable = dof_to_primary_variable(idof,istate)
                  select case(variable)
                    ! for gas pressure dof
                    case(GENERAL_GAS_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs gas pressure defined.'
                        call printErrMsg(option)
                      endif
                    ! for air pressure dof
                    case(GENERAL_AIR_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index == 0) then ! air pressure not found
                        ! if air pressure is not available, let's try temperature 
                        real_index = boundary_condition%flow_aux_mapping(GENERAL_TEMPERATURE_INDEX)
                        if (real_index /= 0) then
                          temperature = boundary_condition%flow_aux_real_var(real_index,iconn)
                          call EOSWaterSaturationPressure(temperature,saturation_pressure,ierr)
                          ! now verify whether gas pressure is provided through BC
                          if (boundary_condition%flow_bc_type(ONE_INTEGER) == NEUMANN_BC) then
                            gas_pressure = xxbc(ONE_INTEGER)
                          else
                            real_index = boundary_condition%flow_aux_mapping(GENERAL_GAS_PRESSURE_INDEX)
                            if (real_index /= 0) then
                              gas_pressure = boundary_condition%flow_aux_real_var(real_index,iconn)
                            else
                              option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                                trim(boundary_condition%flow_condition%name) // &
                                '" needs gas pressure defined to calculate air ' // &
                                'pressure from temperature.'
                              call printErrMsg(option)
                            endif
                          endif
                          xxbc(idof) = gas_pressure - saturation_pressure
                        else
                          option%io_buffer = 'Cannot find boundary constraint for air pressure.'
                          call printErrMsg(option)
                        endif
                      else
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      endif
                    ! for gas saturation dof
                    case(GENERAL_GAS_SATURATION_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
!geh: should be able to use the saturation within the cell
!                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
!                          trim(boundary_condition%flow_condition%name) // &
!                          '" needs saturation defined.'
!                        call printErrMsg(option)
                      endif
                    case(GENERAL_TEMPERATURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs temperature defined.'
                        call printErrMsg(option)
                      endif
                  end select
                case(NEUMANN_BC)
                case default
                  option%io_buffer = 'Unknown BC type in GeneralUpdateAuxVars().'
                  call printErrMsg(option)
              end select
            enddo  
        end select
      else
        ! we do this for all BCs; Neumann bcs will be set later
        do idof = 1, option%nflowdof
          real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
          if (real_index > 0) then
            xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
          else
            option%io_buffer = 'Error setting up boundary condition in GeneralUpdateAuxVars'
            call printErrMsg(option)
          endif
        enddo
      endif
          
      ! set this based on data given 
      global_auxvars_bc(sum_connection)%istate = istate
      ! flag(2) indicates call from non-perturbation
      option%iflag = 2
      call GeneralAuxVarCompute(xxbc,gen_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                ghosted_id, &
                                option)
      ! update state and update aux var; this could result in two update to 
      ! the aux var as update state updates if the state changes
      call GeneralAuxVarUpdateState(xxbc,gen_auxvars_bc(sum_connection), &
                                    global_auxvars_bc(sum_connection), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%sat_func_id(ghosted_id))%ptr, &
                                    ghosted_id,option)
#ifdef DEBUG_GENERAL_FILEOUTPUT
      if (debug_flag > 0) then
        write(debug_unit,'(a,i5,i3,7es24.15)') 'bc_auxvar:', ghosted_id, &
                           global_auxvars_bc(ghosted_id)%istate, &
                            xxbc(:)
      endif
#endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%General%auxvars_up_to_date = PETSC_TRUE

end subroutine GeneralUpdateAuxVars

! ************************************************************************** !

subroutine GeneralUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  gen_auxvars => patch%aux%General%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! flag(0) indicates call from non-perturbation
    option%iflag = 0
    call GeneralAuxVarCompute(xx_p(local_start:local_end), &
                              gen_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              ghosted_id, &
                              option)
    call GeneralAccumulation(gen_auxvars(ZERO_INTEGER,ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,accum_p(local_start:local_end)) 
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

end subroutine GeneralUpdateFixedAccum

! ************************************************************************** !

subroutine GeneralAccumulation(gen_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Material_module
  use Material_Aux_class
  
  implicit none

  type(general_auxvar_type) :: gen_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: porosity, compressed_porosity, dcompressed_porosity_dp
  PetscReal :: v_over_t
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  v_over_t = material_auxvar%volume / option%flow_dt

  porosity = material_auxvar%porosity
  if (soil_compressibility_index > 0) then
    call MaterialCompressSoil(material_auxvar, &
                              maxval(gen_auxvar%pres(1:2)), &
                              compressed_porosity,dcompressed_porosity_dp)
    porosity = compressed_porosity
  endif
  
  ! accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] * 
    !                           xmol[kmol comp/kmol phase]
    do icomp = 1, option%nflowspec
#ifdef DEBUG_GENERAL
      ! for debug version, aux var entries are initialized to NaNs.  even if
      ! saturation is zero, density may be a NaN.  So the conditional prevents
      ! this calculation.  For non-debug, aux var entries are initialized to
      ! 0.d0
      if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
      Res(icomp) = Res(icomp) + gen_auxvar%sat(iphase) * &
                                gen_auxvar%den(iphase) * &
                                gen_auxvar%xmol(icomp,iphase)
#ifdef DEBUG_GENERAL
      endif
#endif
    enddo
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * v_over_t

  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
#ifdef DEBUG_GENERAL
    ! for debug version, aux var entries are initialized to NaNs.  even if
    ! saturation is zero, density may be a NaN.  So the conditional prevents
    ! this calculation.  For non-debug, aux var entries are initialized to
    ! 0.d0
    if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
    Res(energy_id) = Res(energy_id) + gen_auxvar%sat(iphase) * &
                                      gen_auxvar%den(iphase) * &
                                      gen_auxvar%U(iphase)
#ifdef DEBUG_GENERAL
    endif
#endif
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * gen_auxvar%temp) * v_over_t
                    
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,7es24.15)') 'accum:', Res
  endif
#endif                    

end subroutine GeneralAccumulation

! ************************************************************************** !

subroutine GeneralFlux(gen_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       sir_up, &
                       thermal_conductivity_up, &
                       gen_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area, dist, general_parameter, &
                       option,v_darcy,Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
    
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: xmol(option%nflowspec)
  PetscReal :: perm_up, perm_dn
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: stp_up, stp_dn
  PetscReal :: sat_up, sat_dn
  PetscReal :: temp_ave, stp_ave_over_dist, v_air
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: adv_flux(3), diff_flux(3)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
   
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  perm_ave_over_dist = (perm_up * perm_dn)/(dist_up*perm_dn + dist_dn*perm_up)

  Res = 0.d0
  v_darcy = 0.d0
#ifdef DEBUG_FLUXES  
  adv_flux = 0.d0
  diff_flux = 0.d0
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

#ifdef CONVECTION
  do iphase = 1, option%nphase
 
    if (gen_auxvar_up%mobility(iphase) + &
        gen_auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif

    density_kg_ave = GeneralAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           gen_auxvar_up%den_kg, &
                                           gen_auxvar_dn%den_kg)
    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = gen_auxvar_up%pres(iphase) - &
                     gen_auxvar_dn%pres(iphase) + &
                     gravity_term

#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_dphi(iphase) = delta_pressure
#endif

    if (delta_pressure >= 0.D0) then
      mobility = gen_auxvar_up%mobility(iphase)
      xmol(:) = gen_auxvar_up%xmol(:,iphase)
      H_ave = gen_auxvar_up%H(iphase)
      uH = H_ave
    else
      mobility = gen_auxvar_dn%mobility(iphase)
      xmol(:) = gen_auxvar_dn%xmol(:,iphase)
      H_ave = gen_auxvar_dn%H(iphase)
      uH = H_ave
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
      density_ave = GeneralAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          gen_auxvar_up%den, &
                                          gen_auxvar_dn%den)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave
      ! Res[kmol total/sec]
      do icomp = 1, option%nflowspec
        ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
        !                      xmol[kmol comp/kmol phase]
        Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      enddo
#ifdef DEBUG_FLUXES  
      do icomp = 1, option%nflowspec
        adv_flux(icomp) = adv_flux(icomp) + mole_flux * xmol(icomp)
      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      do icomp = 1, option%nflowspec
        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
#endif
      Res(energy_id) = Res(energy_id) + mole_flux * uH
#ifdef DEBUG_FLUXES  
      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_dphi(iphase) = delta_pressure
      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
#endif
    endif                   

  enddo
! CONVECTION
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'delta pressure :', debug_dphi(:)
    write(debug_unit,'(a,7es24.15)') 'adv flux (liquid):', debug_flux(:,1)
    write(debug_unit,'(a,7es24.15)') 'adv flux (gas):', debug_flux(:,2)
  endif
  debug_flux = 0.d0
#endif                    

#ifdef DIFFUSION
  ! add in gas component diffusion in gas and liquid phases
  do iphase = 1, option%nphase
    
#ifndef LIQUID_DIFFUSION  
    if (iphase == LIQUID_PHASE) cycle
#endif    
    
    sat_up = gen_auxvar_up%sat(iphase)
    sat_dn = gen_auxvar_dn%sat(iphase)
    !geh: changed to .and. -> .or.
    if (sqrt(sat_up*sat_dn) < eps) cycle
    if (sat_up > eps .or. sat_dn > eps) then
      ! for now, if liquid state neighboring gas, we allow for minute
      ! diffusion in liquid phase.
      if (iphase == option%liquid_phase) then
        if ((sat_up > eps .or. sat_dn > eps)) then
          sat_up = max(sat_up,eps)
          sat_dn = max(sat_dn,eps)
        endif
      endif
      density_ave = GeneralAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          gen_auxvar_up%den, &
                                          gen_auxvar_dn%den)
      stp_up = sat_up*material_auxvar_up%tortuosity*material_auxvar_up%porosity
      stp_dn = sat_dn*material_auxvar_dn%tortuosity*material_auxvar_dn%porosity
      stp_ave_over_dist = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
      delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                   gen_auxvar_dn%xmol(air_comp_id,iphase)
      ! need to account for multiple phases
      ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
      if (iphase == option%liquid_phase) then
        ! Eq. 1.9a.  The water density is added below
        v_air = stp_ave_over_dist * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol
      else
        temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
        pressure_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+gen_auxvar_dn%pres(iphase))
        ! Eq. 1.9b.  The gas density is added below
        v_air = stp_ave_over_dist * &
                ((temp_ave+273.15)/273.15d0)**1.8d0 * &
                101325.d0 / pressure_ave * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol      
      endif      
      q =  v_air * area
      mole_flux = q * density_ave
      Res(wat_comp_id) = Res(wat_comp_id) - mole_flux
      Res(air_comp_id) = Res(air_comp_id) + mole_flux
#ifdef DEBUG_FLUXES  
      diff_flux(wat_comp_id) = diff_flux(wat_comp_id) - mole_flux
      diff_flux(air_comp_id) = diff_flux(air_comp_id) + mole_flux      
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_flux(wat_comp_id,iphase) = debug_flux(wat_comp_id,iphase) - mole_flux 
      debug_flux(air_comp_id,iphase) = debug_flux(air_comp_id,iphase) + mole_flux 
#endif
    endif
  enddo
! DIFFUSION
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
  k_eff_up = thermal_conductivity_up(1) + &
             sqrt(gen_auxvar_up%sat(option%liquid_phase)) * &
             (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  k_eff_dn = thermal_conductivity_dn(1) + &
             sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
             (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  else
    k_eff_ave = 0.d0
  endif
  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = gen_auxvar_up%temp - gen_auxvar_dn%temp
  heat_flux = k_eff_ave * delta_temp * area
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux * 1.d-6 ! J/s -> MJ/s
! CONDUCTION
#endif
#ifdef DEBUG_FLUXES  
  diff_flux(3) = diff_flux(3) + heat_flux * 1.d-6

  if (option%iflag == 1) then  
    write(*,'(a,7es12.4)') 'in: ', adv_flux(:)*dist(3), diff_flux(:)*dist(3)
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux * 1.d-6
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'dif flux (liquid):', debug_flux(:,1)
    write(debug_unit,'(a,7es24.15)') 'dif flux (gas):', debug_flux(:,2)
  endif
#endif

end subroutine GeneralFlux

! ************************************************************************** !

subroutine GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         gen_auxvar_up,global_auxvar_up, &
                         gen_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         sir_dn, &
                         thermal_conductivity_dn, &
                         area,dist,general_parameter, &
                         option,v_darcy,Res)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module                              
  use Material_Aux_class
  
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy(option%nphase), area
  type(general_parameter_type) :: general_parameter
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  PetscInt :: bc_type
  PetscReal :: xmol(option%nflowspec)  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_ave_over_dist, dist_gravity
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_dn, perm_dn
  PetscReal :: temp_ave, stp_ave_over_dist, v_air, pres_ave
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: adv_flux(3), diff_flux(3)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  PetscReal :: boundary_pressure
  
  PetscInt :: idof
  PetscBool :: neumann_bc_present
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  Res = 0.d0
  v_darcy = 0.d0  
#ifdef DEBUG_FLUXES    
  adv_flux = 0.d0
  diff_flux = 0.d0
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

  neumann_bc_present = PETSC_FALSE
  
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
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
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(GENERAL_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(GENERAL_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn / dist(0)
        endif
        
          
        ! using residual saturation cannot be correct! - geh
        ! reusing sir_dn for bounary auxvar
#define BAD_MOVE1 ! this works
#ifndef BAD_MOVE1       
        if (gen_auxvar_up%sat(iphase) > sir_dn(iphase) .or. &
            gen_auxvar_dn%sat(iphase) > sir_dn(iphase)) then
#endif
          boundary_pressure = gen_auxvar_up%pres(iphase)
          if (iphase == LIQUID_PHASE .and. &
              global_auxvar_up%istate == GAS_STATE) then
            ! the idea here is to accommodate a free surface boundary
            ! face.  this will not work for an interior grid cell as
            ! there should be capillary pressure in force.
            boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
          endif
          density_kg_ave = GeneralAverageDensity(iphase, &
                                                 global_auxvar_up%istate, &
                                                 global_auxvar_dn%istate, &
                                                 gen_auxvar_up%den_kg, &
                                                 gen_auxvar_dn%den_kg)
          gravity_term = density_kg_ave * dist_gravity
          delta_pressure = boundary_pressure - &
                           gen_auxvar_dn%pres(iphase) + &
                           gravity_term

#ifdef DEBUG_GENERAL_FILEOUTPUT
          debug_dphi(iphase) = delta_pressure
#endif

          if (bc_type == SEEPAGE_BC .or. &
              bc_type == CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                gen_auxvar_up%pres(iphase) - &
                 option%reference_pressure < eps) then
              delta_pressure = 0.d0
            endif
          endif
            
          if (delta_pressure >= 0.D0) then
            mobility = gen_auxvar_up%mobility(iphase)
            xmol(:) = gen_auxvar_up%xmol(:,iphase)
            uH = gen_auxvar_up%H(iphase)
          else
            mobility = gen_auxvar_dn%mobility(iphase)
            xmol(:) = gen_auxvar_dn%xmol(:,iphase)
            uH = gen_auxvar_dn%H(iphase)
          endif      

          if (mobility > floweps) then
            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
            v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
            ! only need average density if velocity > 0.
            density_ave = GeneralAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                gen_auxvar_up%den, &
                                                gen_auxvar_dn%den)
          endif
#ifndef BAD_MOVE1        
        endif ! sat > eps
#endif

      case(NEUMANN_BC)
        select case(iphase)
          case(LIQUID_PHASE)
            idof = auxvar_mapping(GENERAL_LIQUID_FLUX_INDEX)
          case(GAS_PHASE)
            idof = auxvar_mapping(GENERAL_GAS_FLUX_INDEX)
        end select
        
        neumann_bc_present = PETSC_TRUE
        xmol = 0.d0
        xmol(iphase) = 1.d0
        if (dabs(auxvars(idof)) > floweps) then
          v_darcy(iphase) = auxvars(idof)
          if (v_darcy(iphase) > 0.d0) then 
            density_ave = gen_auxvar_up%den(iphase)
            uH = gen_auxvar_up%H(iphase)
          else 
            density_ave = gen_auxvar_dn%den(iphase)
            uH = gen_auxvar_dn%H(iphase)
          endif 
        endif
      case default
        option%io_buffer = &
          'Boundary condition type not recognized in GeneralBCFlux phase loop.'
        call printErrMsg(option)
    end select

    if (dabs(v_darcy(iphase)) > 0.d0) then
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                              density_ave[kmol phase/m^3 phase]
      mole_flux = q*density_ave       
      ! Res[kmol total/sec]
      do icomp = 1, option%nflowspec
        ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
        !                      xmol[kmol comp/mol phase]
        Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      enddo
#ifdef DEBUG_FLUXES  
      do icomp = 1, option%nflowspec
        adv_flux(icomp) = adv_flux(icomp) + mole_flux * xmol(icomp)
      enddo
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      do icomp = 1, option%nflowspec
        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
      enddo
#endif
      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH ! H_ave
#ifdef DEBUG_FLUXES  
      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
#endif
    endif
  enddo
! CONVECTION
#endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then 
    write(debug_unit,'(a,7es24.15)') 'bc delta pressure :', debug_dphi(:)  
    write(debug_unit,'(a,7es24.15)') 'bc adv flux (liquid):', debug_flux(:,1)
    write(debug_unit,'(a,7es24.15)') 'bc adv flux (gas):', debug_flux(:,2)
  endif
  debug_flux = 0.d0
#endif  

#ifdef DIFFUSION
  ! add in gas component diffusion in gas and liquid phases
  do iphase = 1, option%nphase
  
#ifdef LIQUID_DIFFUSION    
!    if (neumann_bc_present) cycle
    if (ibndtype(iphase) == NEUMANN_BC) cycle
#else
    if (iphase == LIQUID_PHASE) cycle
#endif
    
    ! diffusion all depends upon the downwind cell.  phase diffusion only
    ! occurs if a phase exists in both auxvars (boundary and internal) or
    ! a liquid phase exists in the internal cell. so, one could say that
    ! liquid diffusion always exists as the internal cell has a liquid phase,
    ! but gas phase diffusion only occurs if the internal cell has a gas
    ! phase.
    if (gen_auxvar_dn%sat(iphase) > eps) then
      sat_dn = gen_auxvar_dn%sat(iphase)
      density_ave = GeneralAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          gen_auxvar_up%den, &
                                          gen_auxvar_dn%den)
      ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) 
      !       = m^3 water/m^4 bulk 
      !    stp_ave = tor_dn*por_dn*(sat_up*sat_dn)/ &
      !              ((sat_up+sat_dn)*dist_dn)
      ! should saturation be distance weighted?
      stp_ave_over_dist = material_auxvar_dn%tortuosity * &
                          material_auxvar_dn%porosity * &
                          sat_dn / dist(0)
      delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                   gen_auxvar_dn%xmol(air_comp_id,iphase)
      ! need to account for multiple phases
      ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
      if (iphase == option%liquid_phase) then
        ! Eq. 1.9a.  The water density is added below
        v_air = stp_ave_over_dist * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol
      else
        temp_ave = 0.5d0*(gen_auxvar_up%temp + gen_auxvar_dn%temp)
        pres_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+gen_auxvar_dn%pres(iphase))
        ! Eq. 1.9b.  The gas density is added below
        v_air = stp_ave_over_dist * &
                ((temp_ave+273.15)/273.15d0)**1.8d0 * &
                101325.d0 / &
                pres_ave * &
                general_parameter%diffusion_coefficient(iphase) * delta_xmol      
      endif
      q =  v_air * area
      mole_flux = q * density_ave
      Res(wat_comp_id) = Res(wat_comp_id) - mole_flux
      Res(air_comp_id) = Res(air_comp_id) + mole_flux
#ifdef DEBUG_FLUXES  
      ! equal but opposite
      diff_flux(wat_comp_id) = diff_flux(wat_comp_id) - mole_flux
      diff_flux(air_comp_id) = diff_flux(air_comp_id) + mole_flux
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_flux(wat_comp_id,iphase) = debug_flux(wat_comp_id,iphase) - mole_flux
      debug_flux(air_comp_id,iphase) = debug_flux(air_comp_id,iphase) + mole_flux
#endif
    endif
  enddo
! DIFFUSION
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  select case (ibndtype(GENERAL_ENERGY_EQUATION_INDEX))
    case (DIRICHLET_BC)
      ! based on Somerton et al., 1974:
      ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
      k_eff_dn = thermal_conductivity_dn(1) + &
                 sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
                 (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
      ! units:
      ! k_eff = W/K/m/m = J/s/K/m/m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = J/s
      k_eff_ave = k_eff_dn / dist(0)
      delta_temp = gen_auxvar_up%temp - gen_auxvar_dn%temp
      heat_flux = k_eff_ave * delta_temp * area
    case(NEUMANN_BC)
                  ! flux prescribed as W/m^2
      heat_flux = auxvars(auxvar_mapping(GENERAL_ENERGY_FLUX_INDEX)) * area

    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'GeneralBCFlux heat conduction loop.'
      call printErrMsg(option)
  end select
  Res(energy_id) = Res(energy_id) + heat_flux * 1.d-6 ! J/s -> MJ/s
! CONDUCTION
#endif

#ifdef DEBUG_FLUXES  
  diff_flux(3) = diff_flux(3) + heat_flux * 1.d-6
  if (option%iflag == 1) then
    write(*,'(a,7es12.4)') 'bc: ', adv_flux(:)*dist(3), diff_flux(:)*dist(3)
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux * 1.d-6
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'bc dif flux (liquid):', debug_flux(:,1)*dist(3)
    write(debug_unit,'(a,7es24.15)') 'bc dif flux (gas):', debug_flux(:,2)*dist(3)
  endif
#endif
  
end subroutine GeneralBCFlux

! ************************************************************************** !

subroutine GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                          gen_auxvar,global_auxvar,scale,Res)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
      
  PetscReal :: qsrc_mol(option%nphase)
  PetscReal :: den, den_kg, enthalpy, internal_energy
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: icomp, ierr

  Res = 0.d0
  do icomp = 1, option%nflowspec
    select case(flow_src_sink_type)
      case(MASS_RATE_SS)
        qsrc_mol(icomp) = qsrc(icomp)/fmw_comp(icomp) ! kg/sec -> kmol/sec
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_mol(icomp) = qsrc(icomp)/fmw_comp(icomp)*scale 
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec
        qsrc_mol(icomp) = qsrc(icomp)*gen_auxvar%den(icomp) ! den = kmol/m^3
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_mol(icomp) = qsrc(icomp)*gen_auxvar%den(icomp)*scale 
    end select
    Res(icomp) = qsrc_mol(icomp)
  enddo
  ! energy units: MJ/sec
  if (size(qsrc) == THREE_INTEGER) then
    if (dabs(qsrc(THREE_INTEGER)) < 1.d-40) then
      cell_pressure = &
        maxval(gen_auxvar%pres(option%liquid_phase:option%gas_phase))
      if (dabs(qsrc(ONE_INTEGER)) > 0.d0) then
        call EOSWaterDensityEnthalpy(gen_auxvar%temp,cell_pressure, &
                                     den_kg,den,enthalpy,ierr)
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
        ! enthalpy units: MJ/kmol
        Res(option%energy_id) = Res(option%energy_id) + &
                                qsrc_mol(ONE_INTEGER) * enthalpy
      endif
      if (dabs(qsrc(TWO_INTEGER)) > 0.d0) then
        ! this is pure air, we use the enthalpy of air, NOT the air/water
        ! mixture in gas
        ! air enthalpy is only a function of temperature and the 
        dummy_pressure = 0.d0
        call EOSGasEnergy(gen_auxvar%temp,dummy_pressure, &
                          enthalpy,internal_energy,ierr)
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> MJ/kmol                                  
        ! enthalpy units: MJ/kmol
        Res(option%energy_id) = Res(option%energy_id) + &
          qsrc_mol(TWO_INTEGER) * enthalpy
      endif
    else
      Res(option%energy_id) = qsrc(THREE_INTEGER) ! MJ/s
    endif
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'src/sink:', Res(1)-Res(2),Res(12:3)
  endif
#endif   
  
end subroutine GeneralSrcSink

! ************************************************************************** !

subroutine GeneralAccumDerivative(gen_auxvar,material_auxvar, &
                                  soil_heat_capacity,option,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Saturation_Function_module
  use Material_Aux_class
  
  implicit none

  type(general_auxvar_type) :: gen_auxvar(0:)
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

!geh:print *, 'GeneralAccumDerivative'

  call GeneralAccumulation(gen_auxvar(ZERO_INTEGER), &
                           material_auxvar,soil_heat_capacity,option,res)
                           
  do idof = 1, option%nflowdof
    call GeneralAccumulation(gen_auxvar(idof), &
                             material_auxvar,soil_heat_capacity, &
                             option,res_pert)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar(idof)%pert
!geh:print *, irow, idof, J(irow,idof), gen_auxvar(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_isothermal) then
    J(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    J(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    J(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    J(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,10es24.15)') 'accum deriv:', J
  endif
#endif

end subroutine GeneralAccumDerivative

! ************************************************************************** !

subroutine GeneralFluxDerivative(gen_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 sir_up, &
                                 thermal_conductivity_up, &
                                 gen_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 sir_dn, &
                                 thermal_conductivity_dn, &
                                 area, dist, &
                                 general_parameter, &
                                 option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up(0:), gen_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0
  
!geh:print *, 'GeneralFluxDerivative'
  option%iflag = -2
  call GeneralFlux(gen_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up,sir_up, &
                   thermal_conductivity_up, &
                   gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn,sir_dn, &
                   thermal_conductivity_dn, &
                   area,dist,general_parameter, &
                   option,v_darcy,res)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_auxvar_up(idof),global_auxvar_up, &
                     material_auxvar_up,sir_up, &
                     thermal_conductivity_up, &
                     gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn,sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_up(idof)%pert
!geh:print *, 'up: ', irow, idof, Jup(irow,idof), gen_auxvar_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                     material_auxvar_up,sir_up, &
                     thermal_conductivity_up, &
                     gen_auxvar_dn(idof),global_auxvar_dn, &
                     material_auxvar_dn,sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_dn(idof)%pert
!geh:print *, 'dn: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_isothermal) then
    Jup(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jup(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
    Jdn(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jup(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jup(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
    Jdn(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,20es24.15)') 'flux deriv:', Jup, Jdn
  endif
#endif
  
end subroutine GeneralFluxDerivative

! ************************************************************************** !

subroutine GeneralBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   gen_auxvar_up, &
                                   global_auxvar_up, &
                                   gen_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   sir_dn, &
                                   thermal_conductivity_dn, &
                                   area,dist,general_parameter, &
                                   option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module 
  use Material_Aux_class
  
  implicit none

  PetscReal :: auxvars(:) ! from aux_real_var array
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jdn = 0.d0
!geh:print *, 'GeneralBCFluxDerivative'

  option%iflag = -2
  call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     gen_auxvar_up,global_auxvar_up, &
                     gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res)                     
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                       gen_auxvar_up,global_auxvar_up, &
                       gen_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area,dist,general_parameter, &
                       option,v_darcy,res_pert)   
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_dn(idof)%pert
!print *, 'bc: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_isothermal) then
    Jdn(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jdn(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,10es24.15)') 'bc flux deriv:', Jdn
  endif
#endif
  
end subroutine GeneralBCFluxDerivative

! ************************************************************************** !

subroutine GeneralSrcSinkDerivative(option,qsrc,flow_src_sink_type, &
                                    gen_auxvars,global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(general_auxvar_type) :: gen_auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  option%iflag = -3
  call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                      gen_auxvars(ZERO_INTEGER),global_auxvar,scale,res)
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                        gen_auxvars(idof),global_auxvar,scale,res_pert)            
    do irow = 1, option%nflowdof
      Jac(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvars(idof)%pert
    enddo !irow
  enddo ! idof
  
  if (general_isothermal) then
    Jac(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jac(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jac(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jac(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,20es24.15)') 'src/sink deriv:', Jac
  endif
#endif

end subroutine GeneralSrcSinkDerivative

! ************************************************************************** !

subroutine GeneralResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Realization_class
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
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  Mat, parameter :: null_mat = 0
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  type(general_parameter_type), pointer :: general_parameter
  type(general_auxvar_type), pointer :: gen_auxvars(:,:), gen_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: iphase
  PetscReal :: scale
  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word

  PetscInt :: icap_up, icap_dn
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)
  
  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  gen_auxvars => patch%aux%General%auxvars
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  general_parameter => patch%aux%General%general_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
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
  ! These 3 must be called before GeneralUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)
  
  option%variables_swapped = PETSC_FALSE
                                             ! do update state
  call GeneralUpdateAuxVars(realization,PETSC_TRUE)

! for debugging a single grid cell
!  i = 90
!  call GeneralOutputAuxVars(gen_auxvars(0,i),global_auxvars(i),i,'genaux', &
!                            PETSC_TRUE,option)

  ! override flags since they will soon be out of date
  patch%aux%General%auxvars_up_to_date = PETSC_FALSE 
  if (option%variables_swapped) then
    !geh: since this operation is not collective (i.e. all processors may
    !     not swap), this operation may fail....
    call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                     NFLOWDOF)
  endif

  if (option%compute_mass_balance_new) then
    call GeneralZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  r_p = -accum_p
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  ! accumulation at t(k+1)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call GeneralAccumulation(gen_auxvars(ZERO_INTEGER,ghosted_id), &
                              material_auxvars(ghosted_id), &
                              material_parameter%soil_heat_capacity(imat), &
                              option,Res) 
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
  enddo

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
   
      call GeneralFlux(gen_auxvars(ZERO_INTEGER,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       material_parameter%soil_residual_saturation(:,icap_up), &
                       material_parameter%soil_thermal_conductivity(:,imat_up), &
                       gen_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       material_parameter%soil_residual_saturation(:,icap_dn), &
                       material_parameter%soil_thermal_conductivity(:,imat_dn), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       general_parameter,option,v_darcy,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy
      if (associated(patch%internal_fluxes)) then
        patch%internal_fluxes(:,1,sum_connection) = Res(:)
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
  boundary_condition => patch%boundary_conditions%first
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

      call GeneralBCFlux(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     gen_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     gen_auxvars(ZERO_INTEGER,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     material_parameter%soil_residual_saturation(:,icap_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     general_parameter,option, &
                     v_darcy,Res)
      patch%boundary_velocities(:,sum_connection) = v_darcy
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
  source_sink => patch%source_sinks%first 
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
      
      call GeneralSrcSink(option,source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                        source_sink%flow_condition%general%rate%itype, &
                        gen_auxvars(ZERO_INTEGER,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Res)

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)
      
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%General%inactive_cells_exist) then
    do i=1,patch%aux%General%n_inactive_rows
      r_p(patch%aux%General%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  call GeneralSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                        gen_auxvars,option)
  
  if (general_isothermal) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+GENERAL_ENERGY_EQUATION_INDEX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  if (general_no_air) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+GENERAL_GAS_EQUATION_INDEX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif  

#ifdef DEBUG_GENERAL_FILEOUTPUT
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'fixed residual:', local_id, &
      accum_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'residual:', local_id, &
      r_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
#endif
  
  
  if (realization%debug%vecview_residual) then
    string = 'Gresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'Gxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    close(debug_unit)
  endif
#endif
  
end subroutine GeneralResidual

! ************************************************************************** !

subroutine GeneralJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Realization_class
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
  type(realization_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscInt :: icap_up,icap_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id
  PetscInt :: irow
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = 0
  
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
  type(general_parameter_type), pointer :: general_parameter
  type(general_auxvar_type), pointer :: gen_auxvars(:,:), gen_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  gen_auxvars => patch%aux%General%auxvars
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  general_parameter => patch%aux%General%general_parameter
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

#ifdef DEBUG_GENERAL_FILEOUTPUT
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
    call GeneralAuxVarPerturb(gen_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              ghosted_id,option)
  enddo
  
#ifdef DEBUG_GENERAL_LOCAL
  call GeneralOutputAuxVars(gen_auxvars,global_auxvars,option)
#endif 

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call GeneralAccumDerivative(gen_auxvars(:,ghosted_id), &
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
                              
      call GeneralFluxDerivative(gen_auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     material_parameter%soil_residual_saturation(:,icap_up), &
                     material_parameter%soil_thermal_conductivity(:,imat_up), &
                     gen_auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     material_parameter%soil_residual_saturation(:,icap_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     general_parameter,option,&
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
  boundary_condition => patch%boundary_conditions%first
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

      call GeneralBCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      gen_auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      gen_auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      material_parameter%soil_residual_saturation(:,icap_dn), &
                      material_parameter%soil_thermal_conductivity(:,imat_dn), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      general_parameter,option, &
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
  source_sink => patch%source_sinks%first 
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
      call GeneralSrcSinkDerivative(option, &
                        source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                        source_sink%flow_condition%general%rate%itype, &
                        gen_auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo
  
  call GeneralSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
                        gen_auxvars,option)

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
  if (patch%aux%General%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%General%n_inactive_rows, &
                          patch%aux%General%inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          ierr);CHKERRQ(ierr)
  endif

  if (general_isothermal) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => patch%aux%General%row_zeroing_array
    ! zero energy residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        GENERAL_ENERGY_EQUATION_INDEX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
  endif

  if (general_no_air) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => patch%aux%General%row_zeroing_array
    ! zero gas component mass balance residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        GENERAL_GAS_EQUATION_INDEX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
  endif
  
  if (realization%debug%matview_Jacobian) then
    string = 'Gjacobian'
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

#if 0
  imat = 1
  if (imat == 1) then
    call GeneralNumericalJacTest(xx,realization) 
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
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

end subroutine GeneralJacobian

! ************************************************************************** !

subroutine GeneralCreateZeroArray(patch,option)
  ! 
  ! Computes the zeroed rows for inactive grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

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
  PetscInt :: n_inactive_rows
  PetscInt, pointer :: inactive_rows_local(:)
  PetscInt, pointer :: inactive_rows_local_ghosted(:)
  PetscErrorCode :: ierr
    
  flag = 0
  grid => patch%grid
  
  n_inactive_rows = 0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      n_inactive_rows = n_inactive_rows + option%nflowdof
    endif
  enddo

  allocate(inactive_rows_local(n_inactive_rows))
  allocate(inactive_rows_local_ghosted(n_inactive_rows))

  inactive_rows_local = 0
  inactive_rows_local_ghosted = 0
  ncount = 0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      do idof = 1, option%nflowdof
        ncount = ncount + 1
        inactive_rows_local(ncount) = (local_id-1)*option%nflowdof+idof
        inactive_rows_local_ghosted(ncount) = (ghosted_id-1)*option%nflowdof + &
                                              idof-1
      enddo
    endif
  enddo

  patch%aux%General%inactive_rows_local => inactive_rows_local
  patch%aux%General%inactive_rows_local_ghosted => inactive_rows_local_ghosted
  patch%aux%General%n_inactive_rows = n_inactive_rows
  
  call MPI_Allreduce(n_inactive_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
  if (flag > 0) patch%aux%General%inactive_cells_exist = PETSC_TRUE
  if (ncount /= n_inactive_rows) then
    if (option%myrank == option%io_rank) then
      print *, 'Error:  Mismatch in non-zero row count!', ncount, &
        n_inactive_rows
    endif
    stop
  endif

end subroutine GeneralCreateZeroArray

! ************************************************************************** !

subroutine GeneralCheckUpdatePre(line_search,X,dX,changed,realization,ierr)
  ! 
  ! Checks update prior to update
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/06/13
  ! 

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
  PetscBool :: changed
  type(realization_type) :: realization
  
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
#ifdef DEBUG_GENERAL_INFO
  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  character(len=MAXWORDLENGTH) :: cell_id_word
  PetscInt, parameter :: max_cell_id = 10
  PetscInt :: cell_id, cell_locator(0:max_cell_id)
  PetscInt :: i, ii
#endif
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: liquid_pressure_index, gas_pressure_index, air_pressure_index
  PetscInt :: saturation_index, temperature_index
  PetscInt :: lid, gid, apid, cpid, vpid, spid
  PetscReal :: liquid_pressure0, liquid_pressure1, del_liquid_pressure
  PetscReal :: gas_pressure0, gas_pressure1, del_gas_pressure
  PetscReal :: air_pressure0, air_pressure1, del_air_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation
  PetscReal :: xmol_air_in_water0, xmol_air_in_water1, del_xmol_air_in_water
  PetscReal :: max_saturation_change = 0.125d0
  PetscReal :: max_temperature_change = 10.d0
  PetscReal :: min_pressure
  PetscReal :: scale, temp_scale, temp_real
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  SNES :: snes
  PetscInt :: newton_iteration
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  gen_auxvars => realization%patch%aux%General%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars

  patch => realization%patch

  spid = option%saturation_pressure_id
  apid = option%air_pressure_id

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  scale = initial_scale
  if (general_max_it_before_damping > 0 .and. &
      newton_iteration > general_max_it_before_damping) then
    scale = general_damping_factor
  endif

  changed = PETSC_TRUE
  
#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
!!#define LIMIT_MAX_TEMPERATURE_CHANGE
!#define TRUNCATE_LIQUID_PRESSURE
!! TRUNCATE_GAS/AIR_PRESSURE is needed for times when the solve wants
!! to pull them negative.
!#define TRUNCATE_GAS_PRESSURE
!#define TRUNCATE_AIR_PRESSURE

#ifdef DEBUG_GENERAL_INFO
  cell_locator = 0
#endif

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    offset = (local_id-1)*option%nflowdof
    temp_scale = 1.d0
#ifdef DEBUG_GENERAL_INFO
    cell_id = grid%nG2A(ghosted_id)
    write(cell_id_word,*) cell_id
    cell_id_word = '(Cell ' // trim(adjustl(cell_id_word)) // '): '
#endif
    select case(global_auxvars(ghosted_id)%istate)
      case(LIQUID_STATE)
        liquid_pressure_index  = offset + GENERAL_LIQUID_PRESSURE_DOF
        temperature_index  = offset + GENERAL_ENERGY_DOF
        dX_p(liquid_pressure_index) = dX_p(liquid_pressure_index) * &
                                      general_pressure_scale
        temp_scale = 1.d0
        del_liquid_pressure = dX_p(liquid_pressure_index)
        liquid_pressure0 = X_p(liquid_pressure_index)
        liquid_pressure1 = liquid_pressure0 - del_liquid_pressure
        del_temperature = dX_p(temperature_index)
        temperature0 = X_p(temperature_index)
        temperature1 = temperature0 - del_temperature
#ifdef LIMIT_MAX_PRESSURE_CHANGE
        if (dabs(del_liquid_pressure) > general_max_pressure_change) then
          temp_real = dabs(general_max_pressure_change/del_liquid_pressure)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Liquid pressure change scaled to truncate at max_pressure_change: '
          call printMsg(option,string)
          write(string2,*) liquid_pressure0
          string = '  Liquid Pressure 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) liquid_pressure1
          string = '  Liquid Pressure 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_liquid_pressure
          string = 'Liquid Pressure change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif !LIMIT_MAX_PRESSURE_CHANGE
#ifdef TRUNCATE_LIQUID_PRESSURE
        ! truncate liquid pressure change to prevent liquid pressure from 
        ! dropping below the air pressure while in the liquid state
        min_pressure = gen_auxvars(ZERO_INTEGER,ghosted_id)%pres(apid) + &
                       gen_auxvars(ZERO_INTEGER,ghosted_id)%pres(spid)
        if (liquid_pressure1 < min_pressure) then
          temp_real = tolerance * (liquid_pressure0 - min_pressure)
          if (dabs(del_liquid_pressure) > 1.d-40) then
            temp_real = dabs(temp_real / del_liquid_pressure)
          else
            option%io_buffer = 'Something is seriously wrong in ' // &
              'GeneralCheckUpdatePre(liquid<min).  Contact Glenn through ' // &
              'pflotran-dev@googlegroups.com.'
            call printErrMsg(option)
          endif
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Liquid pressure change scaled to prevent liquid ' // &
            'pressure from dropping below air pressure: '
          call printMsg(option,string)
          write(string2,*) liquid_pressure0
          string = '  Liquid pressure 0: ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) liquid_pressure1
          string = '  Liquid pressure 1: ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_liquid_pressure
          string = '  pressure change  : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif !TRUNCATE_LIQUID_PRESSURE  
#ifdef LIMIT_MAX_TEMPERATURE_CHANGE
        if (dabs(del_temperature) > max_temperature_change) then
          temp_real = dabs(max_temperature_change/del_temperature)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Temperature change scaled to truncate at max_temperature_change: '
          call printMsg(option,string)
          write(string2,*) temperature0
          string = '  Temperature 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temperature1
          string = '  Temperature 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_temperature
          string = 'Temperature change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif !LIMIT_MAX_TEMPERATURE_CHANGE        
      case(TWO_PHASE_STATE)
        gas_pressure_index = offset + GENERAL_GAS_PRESSURE_DOF
!        air_pressure_index = offset + 2
        saturation_index = offset + GENERAL_GAS_SATURATION_DOF
        temperature_index  = offset + GENERAL_ENERGY_DOF
        dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                   general_pressure_scale
        if (general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then                                   
          air_pressure_index = offset + GENERAL_ENERGY_DOF
          dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                     general_pressure_scale
          del_air_pressure = dX_p(air_pressure_index)
          air_pressure0 = X_p(air_pressure_index)
          air_pressure1 = air_pressure0 - del_air_pressure
        endif
        temp_scale = 1.d0
        del_gas_pressure = dX_p(gas_pressure_index)
        gas_pressure0 = X_p(gas_pressure_index)
        gas_pressure1 = gas_pressure0 - del_gas_pressure
        del_saturation = dX_p(saturation_index)
        saturation0 = X_p(saturation_index)
        saturation1 = saturation0 - del_saturation
#ifdef LIMIT_MAX_PRESSURE_CHANGE
        if (dabs(del_gas_pressure) > general_max_pressure_change) then
          temp_real = dabs(general_max_pressure_change/del_gas_pressure)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Gas pressure change scaled to truncate at max_pressure_change: '
          call printMsg(option,string)
          write(string2,*) gas_pressure0
          string = '  Gas Pressure 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) gas_pressure1
          string = '  Gas Pressure 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_gas_pressure
          string = 'Gas Pressure change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif
#ifdef TRUNCATE_GAS_PRESSURE
        if (gas_pressure1 <= 0.d0) then
          if (dabs(del_gas_pressure) > 1.d-40) then
            temp_real = tolerance * dabs(gas_pressure0 / del_gas_pressure)
#ifdef DEBUG_GENERAL_INFO
            if (cell_locator(0) < max_cell_id) then
              cell_locator(0) = cell_locator(0) + 1
              cell_locator(cell_locator(0)) = ghosted_id
            endif
            string = trim(cell_id_word) // &
              'Gas pressure change scaled to prevent gas ' // &
              'pressure from dropping below zero: '
            call printMsg(option,string)
            write(string2,*) gas_pressure0
            string = '  Gas pressure 0   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) gas_pressure1
            string = '  Gas pressure 1   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) -1.d0*del_gas_pressure
            string = '  pressure change  : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) temp_real
            string = '          scaling  : ' // adjustl(string2)
            call printMsg(option,string)
#endif
            temp_scale = min(temp_scale,temp_real)
          endif
        endif
#endif !TRUNCATE_GAS_PRESSURE
#ifdef TRUNCATE_AIR_PRESSURE
        if (air_pressure1 <= 0.d0) then
          if (dabs(del_air_pressure) > 1.d-40) then
            temp_real = tolerance * dabs(air_pressure0 / del_air_pressure)
#ifdef DEBUG_GENERAL_INFO
            if (cell_locator(0) < max_cell_id) then
              cell_locator(0) = cell_locator(0) + 1
              cell_locator(cell_locator(0)) = ghosted_id
            endif
            string = trim(cell_id_word) // &
              'Air pressure change scaled to prevent air ' // &
              'pressure from dropping below zero: '
            call printMsg(option,string)
            write(string2,*) air_pressure0
            string = '  Air pressure 0   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) air_pressure1
            string = '  Air pressure 1   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) -1.d0*del_air_pressure
            string = '  pressure change  : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) temp_real
            string = '          scaling  : ' // adjustl(string2)
            call printMsg(option,string)
#endif
            temp_scale = min(temp_scale,temp_real)
          endif
        endif
#endif !TRUNCATE_AIR_PRESSURE
#if defined(TRUNCATE_GAS_PRESSURE) && defined(TRUNCATE_AIR_PRESSURE)
        ! have to factor in scaled update from previous conditionals
        gas_pressure1 = gas_pressure0 - temp_scale * del_gas_pressure
        air_pressure1 = air_pressure0 - temp_scale * del_air_pressure
        if (gas_pressure1 <= air_pressure1) then
!          temp_real = (air_pressure0 - gas_pressure0) / &
!                      (temp_scale * (del_air_pressure - del_gas_pressure))
!          temp_real = temp_real * tolerance * temp_scale
          ! refactored to prevent divide by 0
          temp_real = temp_scale * (del_air_pressure - del_gas_pressure)
          if (temp_real > 0.d0) then
            temp_real = (air_pressure0 - gas_pressure0) / temp_real
            temp_real = temp_real * tolerance * temp_scale
          endif
          if (temp_real <= 0.d0) then
            ! add info on pressures here
            option%io_buffer = 'Something is seriously wrong in ' // &
              'GeneralCheckUpdatePre(gas<=air).  Contact Glenn through ' // &
              'pflotran-dev@googlegroups.com.'
            call printErrMsg(option)
          endif
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Gas/Air pressure change scaled again to prevent gas ' // &
            'pressure from dropping below air pressure: '
          call printMsg(option,string)
          write(string2,*) gas_pressure0
          string = '  Gas pressure 0       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) gas_pressure1
          string = '  Gas pressure 1       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*temp_real*del_gas_pressure
          string = '  Gas pressure change  : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) air_pressure0
          string = '  Air pressure 0       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) air_pressure1
          string = '  Air pressure 1       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*temp_real*del_air_pressure
          string = '  Air pressure change  : ' // adjustl(string2)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif !TRUNCATE_GAS_PRESSURE && TRUNCATE_AIR_PRESSURE
#ifdef LIMIT_MAX_SATURATION_CHANGE
        if (dabs(del_saturation) > max_saturation_change) then
          temp_real = dabs(max_saturation_change/del_saturation)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Gas saturation change scaled to truncate at ' // &
            'max_saturation_change: '
          call printMsg(option,string)
          write(string2,*) saturation0
          string = '  Saturation 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) saturation1
          string = '  Saturation 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_saturation
          string = 'Saturation change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif !LIMIT_MAX_SATURATION_CHANGE        
      case(GAS_STATE) 
        gas_pressure_index = offset + GENERAL_GAS_PRESSURE_DOF
        air_pressure_index = offset + GENERAL_GAS_STATE_AIR_PRESSURE_DOF
        dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                   general_pressure_scale
        dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                   general_pressure_scale
    end select
    scale = min(scale,temp_scale) 
  enddo

  temp_scale = scale
  call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MIN,option%mycomm,ierr)

  if (scale < 0.9999d0) then
#ifdef DEBUG_GENERAL_INFO
    string  = '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    call printMsg(option,string)
    write(string2,*) scale, (grid%nG2A(cell_locator(i)),i=1,cell_locator(0))
    string = 'Final scaling: : ' // adjustl(string2)
    call printMsg(option,string)
    do i = 1, cell_locator(0)
      ghosted_id = cell_locator(i)
      offset = (ghosted_id-1)*option%nflowdof
      write(string2,*) grid%nG2A(ghosted_id)
      string = 'Cell ' // trim(adjustl(string2))
      write(string2,*) global_auxvars(ghosted_id)%istate
      string = trim(string) // ' (State = ' // trim(adjustl(string2)) // ') '
      call printMsg(option,string)
      ! for some reason cannot perform array operation on dX_p(:)
      write(string2,*) (X_p(offset+ii),ii=1,3)
      string = '   Orig. Solution: ' // trim(adjustl(string2))
      call printMsg(option,string)
      write(string2,*) (X_p(offset+ii)-dX_p(offset+ii),ii=1,3)
      string = '  Solution before: ' // trim(adjustl(string2))
      call printMsg(option,string)
      write(string2,*) (X_p(offset+ii)-scale*dX_p(offset+ii),ii=1,3)
      string = '   Solution after: ' // trim(adjustl(string2))
      call printMsg(option,string)
    enddo
    string  = '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    call printMsg(option,string)
#endif
    dX_p = scale*dX_p
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine GeneralCheckUpdatePre

! ************************************************************************** !

subroutine GeneralCheckUpdatePost(line_search,X0,dX,X1,dX_changed, &
                                   X1_changed,realization,ierr)
  ! 
  ! Checks update after to update
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/06/13
  ! 

  use Realization_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_class
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  type(realization_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  
  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(general_auxvar_type), pointer :: general_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  type(material_parameter_type), pointer :: material_parameter
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset , ival, idof
#ifdef DEBUG_GENERAL_INFO
  PetscInt :: icell_max_rel_update(3), icell_max_scaled_residual(3)
  PetscInt :: istate_max_rel_update(3), istate_max_scaled_residual(3)
  character(len=2) :: state_char
#endif
  PetscReal :: dX_X0, R_A
  PetscReal :: inf_norm_rel_update(3,3), global_inf_norm_rel_update(3,3)
  PetscReal :: inf_norm_scaled_residual(3,3), global_inf_norm_scaled_residual(3,3)
  PetscReal :: inf_norm_update(3,3), global_inf_norm_update(3,3)
  PetscReal, parameter :: inf_pres_tol = 1.d-1
  PetscReal, parameter :: inf_temp_tol = 1.d-5
  PetscReal, parameter :: inf_sat_tol = 1.d-6
  PetscReal, parameter :: inf_xmol_tol = 1.d-6
  PetscReal, parameter :: inf_norm_update_tol(3,3) = &
    reshape([inf_pres_tol,inf_xmol_tol,inf_temp_tol, &
             inf_pres_tol,inf_pres_tol,inf_temp_tol, &
             inf_pres_tol,inf_pres_tol,inf_sat_tol], &
            shape(inf_norm_update_tol)) * &
            0.d0
  PetscReal :: temp(3,9), global_temp(3,9)
  PetscMPIInt :: mpi_int
  PetscBool :: converged_abs_update
  PetscBool :: converged_rel_update
  PetscBool :: converged_scaled_residual
  PetscInt :: istate
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  patch => realization%patch
  general_auxvars => patch%aux%General%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  option%converged = PETSC_FALSE
  if (option%flow%check_post_convergence) then
    call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
#ifdef DEBUG_GENERAL_INFO
    icell_max_rel_update = 0
    istate_max_rel_update = 0
    icell_max_scaled_residual = 0
    istate_max_scaled_residual = 0
#endif
    inf_norm_update(:,:) = -1.d20
    inf_norm_rel_update(:,:) = -1.d20
    inf_norm_scaled_residual(:,:) = -1.d20
    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%nflowdof
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      istate = global_auxvars(ghosted_id)%istate
      do idof = 1, option%nflowdof
        ival = offset+idof
        R_A = dabs(r_p(ival)/accum_p(ival))
        dX_X0 = dabs(dX_p(ival)/X0_p(ival))
        inf_norm_update(idof,istate) = max(inf_norm_update(idof,istate), &
                                           dabs(dX_p(ival)))
        if (inf_norm_rel_update(idof,istate) < dX_X0) then
#ifdef DEBUG_GENERAL_INFO
          if (maxval(inf_norm_rel_update(idof,:)) < dX_X0) then
            icell_max_rel_update(idof) = grid%nG2A(ghosted_id)
            istate_max_rel_update(idof) = global_auxvars(ghosted_id)%istate
          endif
#endif
          inf_norm_rel_update(idof,istate) = dX_X0
        endif
        if (inf_norm_scaled_residual(idof,istate) < R_A) then
#ifdef DEBUG_GENERAL_INFO
          if (maxval(inf_norm_scaled_residual(idof,:)) < R_A) then
            icell_max_scaled_residual(idof) = grid%nG2A(ghosted_id)
            istate_max_scaled_residual(idof) = global_auxvars(ghosted_id)%istate
          endif
#endif
          inf_norm_scaled_residual(idof,istate) = R_A
        endif
      enddo
    enddo
    temp(1:3,1:3) = inf_norm_update(:,:)
    temp(1:3,4:6) = inf_norm_rel_update(:,:)
    temp(1:3,7:9) = inf_norm_scaled_residual(:,:)
    mpi_int = 27
    call MPI_Allreduce(temp,global_temp,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    global_inf_norm_update(:,:) = global_temp(1:3,1:3)
    global_inf_norm_rel_update(:,:) = global_temp(1:3,4:6)
    global_inf_norm_scaled_residual(:,:) = global_temp(1:3,7:9)
    converged_abs_update = PETSC_TRUE
    do istate = 1, 3
      do idof = 1, option%nflowdof
        if (global_inf_norm_update(idof,istate) > &
          inf_norm_update_tol(idof,istate)) then
          converged_abs_update = PETSC_FALSE
        endif
      enddo  
    enddo  
    converged_rel_update = maxval(global_inf_norm_rel_update) < &
                           option%flow%inf_rel_update_tol
    converged_scaled_residual = maxval(global_inf_norm_scaled_residual) < &
                                option%flow%inf_scaled_res_tol
#if 0
    do idof = 1, option%nflowdof
      if (global_inf_norm(idof) > option%flow%post_convergence_tol) then
        converged_rel_update = PETSC_FALSE
#ifdef DEBUG_GENERAL_INFO
        select case(istate_max(idof))
          case(1)
            state_char = 'L'
          case(2)
            state_char = 'G'
          case(3)
            state_char = '2P'
        end select
        write(*,'(''-+ '',a3,i2,''('',i5,''):'',es12.4, &
                 &'' dX_X/dX/X:'',3es12.4, &
                 &'' R_A/R/A:'',3es12.4)') state_char,idof, &
           icell_max(idof),global_inf_norm(idof), &
           dX_X1_max(idof), dX_max(idof),  X1_max(idof), &
           R_A_max(idof), R_max(idof), A_max(idof)
#endif
      endif
    enddo
#endif
#ifdef DEBUG_GENERAL_INFO
    write(*,'(4x,''-+  dpl:'',es12.4,''  dxa:'',es12.4,''  dt:'',es12.4)') &
      (max(global_inf_norm_update(idof,1),0.d0),idof=1,3)
    write(*,'(4x,''-+  dpg:'',es12.4,''  dpa:'',es12.4,''  dt:'',es12.4)') &
      (max(global_inf_norm_update(idof,2),0.d0),idof=1,3)
    if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
      write(*,'(4x,''-+  dpg:'',es12.4,''  dsg:'',es12.4,''  dt:'',es12.4)') &
        (max(global_inf_norm_update(idof,3),0.d0),idof=1,3)
    else
      write(*,'(4x,''-+  dpg:'',es12.4,''  dsg:'',es12.4,'' dpa:'',es12.4)') &
        (max(global_inf_norm_update(idof,3),0.d0),idof=1,3)
    endif
    write(*,'(4x,''-+ rupl:'',es12.4,'' ruxa:'',es12.4,'' rut:'',es12.4)') &
      (max(global_inf_norm_rel_update(idof,1),0.d0),idof=1,3)
    write(*,'(4x,''-+ rupg:'',es12.4,'' rupa:'',es12.4,'' rut:'',es12.4)') &
      (max(global_inf_norm_rel_update(idof,2),0.d0),idof=1,3)
    write(*,'(4x,''-+ rupg:'',es12.4,'' rusg:'',es12.4,'' rut:'',es12.4)') &
      (max(global_inf_norm_rel_update(idof,3),0.d0),idof=1,3)
    write(*,'(4x,''-+  srl:'',es12.4,''  srg:'',es12.4,'' sre:'',es12.4)') &
      (max(global_inf_norm_scaled_residual(idof,1),0.d0),idof=1,3)
    write(*,'(4x,''-+  srl:'',es12.4,''  srg:'',es12.4,'' sre:'',es12.4)') &
      (max(global_inf_norm_scaled_residual(idof,2),0.d0),idof=1,3)
    write(*,'(4x,''-+  srl:'',es12.4,''  srg:'',es12.4,'' sre:'',es12.4)') &
      (max(global_inf_norm_scaled_residual(idof,3),0.d0),idof=1,3)
    write(*,'(4x,''-+ ru1 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_rel_update(1), istate_max_rel_update(1), &
      X0_p((icell_max_rel_update(1)-1)*3+1), &
      -1.d0*dX_p((icell_max_rel_update(1)-1)*3+1), &
      r_p((icell_max_rel_update(1)-1)*3+1)
    write(*,'(4x,''-+ ru2 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_rel_update(2), istate_max_rel_update(2), &
      X0_p((icell_max_rel_update(2)-1)*3+2), &
      -1.d0*dX_p((icell_max_rel_update(2)-1)*3+2), &
      r_p((icell_max_rel_update(2)-1)*3+2)
    write(*,'(4x,''-+ ru3 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_rel_update(3), istate_max_rel_update(3), &
      X0_p((icell_max_rel_update(3)-1)*3+3), &
      -1.d0*dX_p((icell_max_rel_update(3)-1)*3+3), &
      r_p((icell_max_rel_update(3)-1)*3+3)
    write(*,'(4x,''-+ sr1 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_scaled_residual(1), istate_max_scaled_residual(1), &
      X0_p((icell_max_scaled_residual(1)-1)*3+1), &
      -1.d0*dX_p((icell_max_scaled_residual(1)-1)*3+1), &
      r_p((icell_max_scaled_residual(1)-1)*3+1)
    write(*,'(4x,''-+ sr2 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_scaled_residual(2), istate_max_scaled_residual(2), &
      X0_p((icell_max_scaled_residual(2)-1)*3+2), &
      -1.d0*dX_p((icell_max_scaled_residual(2)-1)*3+2), &
      r_p((icell_max_scaled_residual(2)-1)*3+2)
    write(*,'(4x,''-+ sr3 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_scaled_residual(3), istate_max_scaled_residual(3), &
      X0_p((icell_max_scaled_residual(3)-1)*3+3), &
      -1.d0*dX_p((icell_max_scaled_residual(3)-1)*3+3), &
      r_p((icell_max_scaled_residual(3)-1)*3+3)
#endif
    option%converged = PETSC_FALSE
    if (converged_abs_update .or. converged_rel_update .or. &
        converged_scaled_residual) then
      option%converged = PETSC_TRUE
    endif
    call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  endif

  call VecGetArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    istate = global_auxvars(ghosted_id)%istate
    if (istate == LIQUID_STATE .and. X1_p(offset+2) < 0.d0) then
        X1_p(offset+2) = 1.d-40
        X1_changed = PETSC_TRUE
    endif
  enddo
  call VecRestoreArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
  
end subroutine GeneralCheckUpdatePost

! ************************************************************************** !

function GeneralGetTecplotHeader(realization,icolumn)
  ! 
  ! Returns General Lite contribution to
  ! Tecplot file header
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  
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

subroutine GeneralSetPlotVariables(realization)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/13
  ! 
  
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
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)
  
  name = 'Gas Density'
  units = 'kg/m^3'
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
  
  name = 'Liquid Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)
  
  name = 'Gas Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_ENERGY)
  
  name = 'Thermodynamic State'
  units = ''
  output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE,units,STATE)
  output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList( &
         realization%output_option%output_variable_list,output_variable)   
  
end subroutine GeneralSetPlotVariables

! ************************************************************************** !

function GeneralAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/14
  ! 

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)

  PetscReal :: GeneralAverageDensity

  if (iphase == LIQUID_PHASE) then
    if (istate_up == GAS_STATE) then
      GeneralAverageDensity = density_dn(iphase)
    else if (istate_dn == GAS_STATE) then
      GeneralAverageDensity = density_up(iphase)
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == LIQUID_STATE) then
      GeneralAverageDensity = density_dn(iphase)
    else if (istate_dn == LIQUID_STATE) then
      GeneralAverageDensity = density_up(iphase)
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
    endif
  endif

end function GeneralAverageDensity

! ************************************************************************** !

subroutine GeneralSSSandbox(residual,Jacobian,compute_derivative, &
                            grid,material_auxvars,general_auxvars,option)
  ! 
  ! Evaluates source/sink term storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Grid_module
  use Material_Aux_class, only: material_auxvar_type
  use SrcSink_Sandbox_module
  use SrcSink_Sandbox_Base_class
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  PetscBool :: compute_derivative
  Vec :: residual
  Mat :: Jacobian
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(general_auxvar_type), pointer :: general_auxvars(:,:)
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscReal, pointer :: r_p(:)
  PetscReal :: res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  PetscInt :: i, local_id, ghosted_id, istart, iend, idof, irow
  PetscReal :: res_pert(option%nflowdof)
  PetscReal :: aux_real(10)
  PetscErrorCode :: ierr
  
  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif
  
  cur_srcsink => sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
      aux_real = 0.d0

      do i = 1, size(cur_srcsink%region%cell_ids)
        local_id = cur_srcsink%region%cell_ids(i)
        ghosted_id = grid%nL2G(local_id)
        res = 0.d0
        Jac = 0.d0
        call GeneralSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                          general_auxvars(ZERO_INTEGER,ghosted_id),option)
        call cur_srcsink%Evaluate(res,Jac,PETSC_FALSE, &
                                  material_auxvars(ghosted_id), &
                                  aux_real,option)
        if (compute_derivative) then
          do idof = 1, option%nflowdof
            res_pert = 0.d0
            call GeneralSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                                       general_auxvars(idof,ghosted_id),option)
            call cur_srcsink%Evaluate(res_pert,Jac,PETSC_FALSE, &
                                      material_auxvars(ghosted_id), &
                                      aux_real,option)
            do irow = 1, option%nflowdof
              Jac(irow,idof) = (res_pert(irow)-res(irow)) / &
                               general_auxvars(idof,ghosted_id)%pert
            enddo
          enddo
          if (general_isothermal) then
            Jac(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
            Jac(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
          endif         
          if (general_no_air) then
            Jac(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
            Jac(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
          endif          
          call MatSetValuesBlockedLocal(Jacobian,1,ghosted_id-1,1, &
                                        ghosted_id-1,Jac,ADD_VALUES, &
                                        ierr);CHKERRQ(ierr)
        else
          iend = local_id*option%nflowdof
          istart = iend - option%nflowdof + 1
          r_p(istart:iend) = r_p(istart:iend) - res
        endif
      enddo
    cur_srcsink => cur_srcsink%next
  enddo
  
  if (.not.compute_derivative) then
    call VecRestoreArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine GeneralSSSandbox

! ************************************************************************** !

subroutine GeneralSSSandboxLoadAuxReal(srcsink,aux_real,gen_auxvar,option)

  use Option_module
  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_WIPP_Gas_class
  use SrcSink_Sandbox_WIPP_Well_class

  implicit none

  class(srcsink_sandbox_base_type) :: srcsink
  PetscReal :: aux_real(:)
  type(general_auxvar_type) gen_auxvar
  type(option_type) :: option
  
  aux_real = 0.d0
  select type(srcsink)
    class is(srcsink_sandbox_wipp_gas_type)
      aux_real(WIPP_GAS_WATER_SATURATION_INDEX) = &
        gen_auxvar%sat(option%liquid_phase)
      aux_real(WIPP_GAS_TEMPERATURE_INDEX) = &
        gen_auxvar%temp
    class is(srcsink_sandbox_wipp_well_type)
      aux_real(WIPP_WELL_LIQUID_MOBILITY) = &
        gen_auxvar%mobility(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_MOBILITY) = &
        gen_auxvar%mobility(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_PRESSURE) = &
        gen_auxvar%pres(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_PRESSURE) = &
        gen_auxvar%pres(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_ENTHALPY) = &
        gen_auxvar%H(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_ENTHALPY) = &
        gen_auxvar%H(option%gas_phase)
      aux_real(WIPP_WELL_XMOL_AIR_IN_LIQUID) = &
        gen_auxvar%xmol(option%air_id,option%liquid_phase)
      aux_real(WIPP_WELL_XMOL_WATER_IN_GAS) = &
        gen_auxvar%xmol(option%water_id,option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_DENSITY) = &
        gen_auxvar%den(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_DENSITY) = &
        gen_auxvar%den(option%gas_phase)
  end select
  
end subroutine GeneralSSSandboxLoadAuxReal

! ************************************************************************** !

subroutine GeneralMapBCAuxvarsToGlobal(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Realization_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(general_auxvar_type), pointer :: gen_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update
  
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        gen_auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        gen_auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        gen_auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine GeneralMapBCAuxvarsToGlobal

! ************************************************************************** !

subroutine GeneralDestroy(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Realization_class

  implicit none

  type(realization_type) :: realization
  
  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine GeneralDestroy

end module General_module
