module TOilIms_module
! Brief description for the module
! Pimary variables for ToilIms: oil_pressure, oil_saturation, temperature

  use TOilIms_Aux_module
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

! Cutoff parameters - no public
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  public :: TOilImsSetup, &
            TOilImsUpdateAuxVars, &
            TOilImsInitializeTimestep, &
            TOilImsCheckUpdatePre, &
            TOilImsUpdateSolution, &
            TOilImsMapBCAuxVarsToGlobal

contains

! ************************************************************************** !

subroutine TOilImsSetup(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/20/15
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  !use Fluid_module
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
  type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)
  type(toil_ims_auxvar_type), pointer :: toil_auxvars_bc(:)
  type(toil_ims_auxvar_type), pointer :: toil_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  !type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%TOil_ims => TOilImsAuxCreate(option)

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
    option%io_buffer = 'Material property errors found in GeneralSetup.'
    call printErrMsg(option)
  endif

  ! allocate auxvar data structures for all grid cells  
  allocate(toil_auxvars(0:option%nflowdof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      !call GeneralAuxVarInit(gen_auxvars(idof,ghosted_id),option)
      call TOilImsAuxVarInit(toil_auxvars(idof,ghosted_id),option)
    enddo
  enddo
  patch%aux%TOil_ims%auxvars => toil_auxvars
  patch%aux%TOil_ims%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them 
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(toil_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call TOilImsAuxVarInit(toil_auxvars_bc(iconn),option)
    enddo
    patch%aux%TOil_ims%auxvars_bc => toil_auxvars_bc
  endif
  patch%aux%TOil_ims%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(toil_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call TOilImsAuxVarInit(toil_auxvars_ss(iconn),option)
    enddo
    patch%aux%TOil_ims%auxvars_ss => toil_auxvars_ss
  endif
  patch%aux%TOil_ims%num_aux_ss = sum_connection

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call TOilImsCreateZeroArray(patch,option)

  ! create array for zeroing Jacobian entries if isothermal
  allocate(patch%aux%TOil_ims%row_zeroing_array(grid%nlmax))
  patch%aux%TOil_ims%row_zeroing_array = 0

  call TOilImsSetPlotVariables(realization)
 
  ! covergence creteria to be chosen (can use TOUGH or general type) 
  !if (general_tough2_conv_criteria .and. &
  !    Initialized(option%flow%inf_scaled_res_tol)) then
  !  ! override what was set in OPTION block of GENERAL process model
  !  general_tough2_itol_scaled_res_e1 = option%flow%inf_scaled_res_tol
  !endif

end subroutine TOilImsSetup

! ************************************************************************** !

subroutine TOilImsInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_class
  
  implicit none
  
  type(realization_type) :: realization

  call TOilImsUpdateFixedAccum(realization)
  

end subroutine TOilImsInitializeTimestep

! ************************************************************************** !

subroutine TOilImsCreateZeroArray(patch,option)
  ! 
  ! Computes the zeroed rows for inactive grid cells
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/20/15
  ! 

  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  !use Field_module
  
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

  patch%aux%TOil_ims%inactive_rows_local => inactive_rows_local
  patch%aux%TOil_ims%inactive_rows_local_ghosted => inactive_rows_local_ghosted
  patch%aux%TOil_ims%n_inactive_rows = n_inactive_rows
  
  call MPI_Allreduce(n_inactive_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
  if (flag > 0) patch%aux%TOil_ims%inactive_cells_exist = PETSC_TRUE
  if (ncount /= n_inactive_rows) then
    if (option%myrank == option%io_rank) then
      print *, 'Error:  Mismatch in non-zero row count!', ncount, &
        n_inactive_rows
    endif
    stop
  endif

end subroutine TOilImsCreateZeroArray

! ************************************************************************** !

subroutine TOilImsCheckUpdatePre(line_search,X,dX,changed,realization,ierr)
  ! 
  ! Checks update prior to update
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/22/15
  ! 

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
  !use Saturation_Function_module
  use Patch_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  type(realization_type) :: realization
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field

  type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  

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
  toil_auxvars => realization%patch%aux%TOil_ims%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars

  patch => realization%patch

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE

  scale = initial_scale
  if (toil_ims_max_it_before_damping > 0 .and. &
      newton_iteration > toil_ims_max_it_before_damping) then
    scale = toil_ims_damping_factor
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

    pressure_index = offset + TOIL_IMS_PRESSURE_DOF
    saturation_index = offset + TOIL_IMS_SATURATION_DOF
    temperature_index  = offset + TOIL_IMS_ENERGY_DOF
    dX_p(pressure_index) = dX_p(pressure_index) * &
                             toil_ims_pressure_scale
    temp_scale = 1.d0
    del_pressure = dX_p(pressure_index)
    pressure0 = X_p(pressure_index)
    pressure1 = pressure0 - del_pressure
    del_saturation = dX_p(saturation_index)
    saturation0 = X_p(saturation_index)
    saturation1 = saturation0 - del_saturation
#ifdef LIMIT_MAX_PRESSURE_CHANGE
    if (dabs(del_pressure) > toil_ims_max_pressure_change) then
      temp_real = dabs(toil_ims_max_pressure_change/del_pressure)
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
#endif !TRUNCATE_PRESSURE

#ifdef LIMIT_MAX_SATURATION_CHANGE
    if (dabs(del_saturation) > max_saturation_change) then
       temp_real = dabs(max_saturation_change/del_saturation)
       temp_scale = min(temp_scale,temp_real)
    endif
#endif !LIMIT_MAX_SATURATION_CHANGE        
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

end subroutine TOilImsCheckUpdatePre

! ************************************************************************** !

subroutine TOilImsSetPlotVariables(realization)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/20/15
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

  name = 'Oil Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               OIL_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)
  
  name = 'Oil Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               OIL_SATURATION)
  
  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)
  
  name = 'Oil Density'
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
  
 !name = 'Thermodynamic State'
 ! units = ''
 ! output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE,units,STATE)
 ! output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
 ! output_variable%iformat = 1 ! integer
 ! call OutputVariableAddToList( &
 !        realization%output_option%output_variable_list,output_variable)   
  
end subroutine TOilImsSetPlotVariables

! ************************************************************************** !

subroutine TOilImsUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with the TOilIms problem
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/21/15
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
  !use EOS_Water_module 
  !use Saturation_Function_module
  
  implicit none

  type(realization_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:), toil_auxvars_bc(:)  

  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)  

  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: ghosted_start, ghosted_end
  !PetscInt :: iphasebc, iphase
  PetscInt :: offset
  PetscInt :: istate

  !PetscReal :: gas_pressure, capillary_pressure, liquid_saturation
  !PetscReal :: saturation_pressure, temperature

  PetscInt :: real_index, variable
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)

  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  toil_auxvars => patch%aux%TOil_ims%auxvars
  toil_auxvars_bc => patch%aux%TOil_ims%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    ! TOIL_IMS_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = TOIL_IMS_UPDATE_FOR_ACCUM

    call TOilImsAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       toil_auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%characteristic_curves_array( &
                         patch%sat_func_id(ghosted_id))%ptr, &
                       ghosted_id, &
                       option)

  enddo
  
  ! compute auxiliary variables for boundary cells
  boundary_condition => patch%boundary_condition_list%first
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
      !istate = boundary_condition%flow_aux_int_var(GENERAL_STATE_INDEX,iconn)

      ! we do this for all BCs; Neumann bcs will be set later
      do idof = 1, option%nflowdof
        real_index = boundary_condition% &
                       flow_aux_mapping(toil_ims_dof_to_primary_vars(idof))
        if (real_index > 0) then
          xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
        else
          option%io_buffer = 'Error setting up boundary condition' // &
                             ' in TOilImsUpdateAuxVars'
          call printErrMsg(option)
        endif
      enddo
        
      ! state not required 
      !toil_auxvars_bc(sum_connection)%istate = istate
      ! TOIL_IMS_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = TOIL_IMS_UPDATE_FOR_BOUNDARY
      call TOilImsAuxVarCompute(xxbc,toil_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                ghosted_id, &
                                option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%TOil_ims%auxvars_up_to_date = PETSC_TRUE

end subroutine TOilImsUpdateAuxVars

! ************************************************************************** !

subroutine TOilImsUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
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
  type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)

  type(global_auxvar_type), pointer :: global_auxvars(:)

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  toil_auxvars => patch%aux%TOil_ims%auxvars

  global_auxvars => patch%aux%Global%auxvars

  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  !Tough2 conv. creteria: initialize accumulation term for every iteration
  if (toil_ims_tough2_conv_criteria) then
    call VecGetArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! TOIL_IMS_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = TOIL_IMS_UPDATE_FOR_FIXED_ACCUM ! not currently used
    call TOilImsAuxVarCompute(xx_p(local_start:local_end), &
                              toil_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              ghosted_id, &
                              option)
    call TOilImsAccumulation(toil_auxvars(ZERO_INTEGER,ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,accum_p(local_start:local_end) )
  enddo
  
  !Tough2 conv. creteria: initialize accumulation term for every iteration
  if (toil_ims_tough2_conv_criteria) then
    accum_p2 = accum_p
  endif
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  !Tough2 conv. creteria: initialize accumulation term for every iteration
  if (toil_ims_tough2_conv_criteria) then
    call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
end subroutine TOilImsUpdateFixedAccum

! ************************************************************************** !

subroutine TOilImsUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step: currently it updates only mass balance
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 

  use Realization_class
  !use Field_module
  !use Patch_module
  !use Discretization_module
  !use Option_module
  !use Grid_module
  
  implicit none
  
  type(realization_type) :: realization

  !type(option_type), pointer :: option
  !type(patch_type), pointer :: patch
  !type(grid_type), pointer :: grid
  !type(field_type), pointer :: field
  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)
  !type(global_auxvar_type), pointer :: global_auxvars(:)  
  !PetscInt :: local_id, ghosted_id
  !PetscErrorCode :: ierr
  
  !option => realization%option
  !field => realization%field
  !patch => realization%patch
  !grid => patch%grid
  !gen_auxvars => patch%aux%General%auxvars  
  !global_auxvars => patch%aux%Global%auxvars
  
  if (realization%option%compute_mass_balance_new) then
    call TOilImsUpdateMassBalance(realization)
  endif
  
  
  
end subroutine TOilImsUpdateSolution

! ************************************************************************** !

subroutine TOilImsUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! Using existing data structure for two phase compositional
  ! For memory efficiency define new data structure 
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use EOS_Oil_module
 
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

  ! option%nflowspec = 2,
  ! two species (H2O,OIL): each present only in its own rich phase

  ! updating with mass balance, assuming molar quantity loaded in:
  ! mass_balance and mass_balance_delta

  !write(*,*) "toil fmw", toil_ims_fmw_comp(1), toil_ims_fmw_comp(2)

  do iconn = 1, patch%aux%TOil_ims%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        toil_ims_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%TOil_ims%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        toil_ims_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine TOilImsUpdateMassBalance

! ************************************************************************** !

subroutine TOilImsZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
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

end subroutine TOilImsZeroMassBalanceDelta

! ************************************************************************** !

subroutine TOilImsMapBCAuxVarsToGlobal(realization)
  ! 
  ! Maps toil_ims BC Auxvars to global auxvars
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
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
  type(toil_ims_auxvar_type), pointer :: toil_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  !NOTE: this is called only if cpoupling to RT
  if (option%ntrandof == 0) return ! no need to update
  
  toil_auxvars_bc => patch%aux%TOil_ims%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        toil_auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        toil_auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        toil_auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine TOilImsMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine TOilImsAccumulation(toil_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 

  use Option_module
  use Material_module
  use Material_Aux_class
  
  implicit none

  type(toil_ims_auxvar_type) :: toil_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  !PetscBool :: debug_cell
  
  PetscInt :: iphase, energy_id
  
  PetscReal :: porosity
  PetscReal :: v_over_t
  
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  v_over_t = material_auxvar%volume / option%flow_dt
  ! must use toil_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = toil_auxvar%effective_porosity
  
  ! flow accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] 
    ! molar balance formulation (kmol)
    Res(iphase) = Res(iphase) + toil_auxvar%sat(iphase) * &
                                toil_auxvar%den(iphase)  
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nphase) = Res(1:option%nphase) * &
                            porosity * v_over_t

  ! energy accumulation term units = MJ/s
  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
    Res(energy_id) = Res(energy_id) + toil_auxvar%sat(iphase) * &
                                      toil_auxvar%den(iphase) * &
                                      toil_auxvar%U(iphase)
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * toil_auxvar%temp) * v_over_t
                    
end subroutine TOilImsAccumulation

! ************************************************************************** !

end module TOilIms_module

