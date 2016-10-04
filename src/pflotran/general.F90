module General_module

  use General_Aux_module
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

  public :: GeneralSetup, &
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
            GeneralMapBCAuxVarsToGlobal, &
            GeneralDestroy

contains

! ************************************************************************** !

subroutine GeneralSetup(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Fluid_module
  use Material_Aux_class
  use Output_Aux_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
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
    option%io_buffer = 'ERROR: Non-initialized soil residual saturation.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_heat_capacity(:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil heat capacity.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_thermal_conductivity(:,:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil thermal conductivity.'
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
      option%io_buffer = 'ERROR: Non-initialized cell volume.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'ERROR: Non-initialized tortuosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%soil_particle_density < 0.d0 .and. &
        flag(4) == 0) then
      flag(4) = 1
      option%io_buffer = 'ERROR: Non-initialized soil particle density.'
      call printMsg(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'ERROR: Non-initialized permeability.'
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
      call GeneralAuxVarInit(gen_auxvars(idof,ghosted_id), &
                             (general_analytical_derivatives .and. idof==0), &
                             option)
    enddo
  enddo
  patch%aux%General%auxvars => gen_auxvars
  patch%aux%General%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them 
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(gen_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_auxvars_bc(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%General%auxvars_bc => gen_auxvars_bc
  endif
  patch%aux%General%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(gen_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_auxvars_ss(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%General%auxvars_ss => gen_auxvars_ss
  endif
  patch%aux%General%num_aux_ss = sum_connection

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
  ! check whether diffusion coefficients are initialized.
  if (Uninitialized(patch%aux%General%general_parameter% &
      diffusion_coefficient(LIQUID_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Liquid phase diffusion coefficient','')
    call printErrMsg(option)
  endif
  if (Uninitialized(patch%aux%General%general_parameter% &
      diffusion_coefficient(GAS_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Gas phase diffusion coefficient','')
    call printErrMsg(option)
  endif

  list => realization%output_option%output_snap_variable_list
  call GeneralSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call GeneralSetPlotVariables(realization,list)
  
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

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

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
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  gen_auxvars => patch%aux%General%auxvars  
  global_auxvars => patch%aux%Global%auxvars
  
  if (realization%option%compute_mass_balance_new) then
    call GeneralUpdateMassBalance(realization)
  endif
  
  ! update stored state
  do ghosted_id = 1, grid%ngmax
    gen_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS) = &
      global_auxvars(ghosted_id)%istate
  enddo
  
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
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  gen_auxvars => patch%aux%General%auxvars

  ! restore stored state
  do ghosted_id = 1, grid%ngmax
    global_auxvars(ghosted_id)%istate = &
      gen_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS)
  enddo

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_timestep_cut_count = debug_timestep_cut_count + 1
#endif 

  call GeneralInitializeTimestep(realization)  

end subroutine GeneralTimeCut

! ************************************************************************** !

subroutine GeneralNumericalJacobianTest(xx,realization,B)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/03/15
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_subsurface_type) :: realization
  Mat :: B

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  PetscReal :: derivative, perturbation
  PetscReal :: perturbation_tolerance = 1.d-6
  PetscInt, save :: icall = 0
  character(len=MAXWORDLENGTH) :: word, word2

  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  icall = icall + 1
  call VecDuplicate(xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res_pert,ierr);CHKERRQ(ierr)

  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof, &
                   grid%nlmax*option%nflowdof, &
                   ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(A,27,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)

  call VecZeroEntries(res,ierr);CHKERRQ(ierr)
  call GeneralResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
#if 0
  word  = 'num_0.dat'
  call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
  call VecView(res,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call VecCopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      call VecZeroEntries(res_pert,ierr);CHKERRQ(ierr)
      call GeneralResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
#if 0
      write(word,*) idof
      word  = 'num_' // trim(adjustl(word)) // '.dat'
      call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
      call VecView(res_pert,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
      call VecGetArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call MatSetValue(A,idof2-1,idof-1,derivative,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if 1
  write(word,*) icall
  word = 'numerical_jacobian-' // trim(adjustl(word)) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

!geh: uncomment to overwrite numerical Jacobian
!  call MatCopy(A,B,DIFFERENT_NONZERO_PATTERN,ierr)
  call MatDestroy(A,ierr);CHKERRQ(ierr)

  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)

end subroutine GeneralNumericalJacobianTest

! ************************************************************************** !

subroutine GeneralComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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
        general_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
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
  type(general_auxvar_type), pointer :: gen_auxvars(:,:), gen_auxvars_bc(:)  
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

  gen_auxvars => patch%aux%General%auxvars
  gen_auxvars_bc => patch%aux%General%auxvars_bc
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
    ! GENERAL_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = GENERAL_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call GeneralAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       gen_auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%characteristic_curves_array( &
                         patch%sat_func_id(ghosted_id))%ptr, &
                       natural_id, &
                       option)
    if (update_state) then
      call GeneralAuxVarUpdateState(xx_loc_p(ghosted_start:ghosted_end), &
                                    gen_auxvars(ZERO_INTEGER,ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%sat_func_id(ghosted_id))%ptr, &
                                    natural_id, &  ! for debugging
                                    option)
    endif
#ifdef DEBUG_AUXVARS
!geh: for debugging
    call GeneralOutputAuxVars(gen_auxvars(0,ghosted_id), &
                              global_auxvars(ghosted_id),natural_id,word, &
                              PETSC_TRUE,option)
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,i5,i3,7es24.15)') 'auxvar:', natural_id, &
                        global_auxvars(ghosted_id)%istate, &
                        xx_loc_p(ghosted_start:ghosted_end)
  endif
#endif
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
      ! GENERAL_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = GENERAL_UPDATE_FOR_BOUNDARY
      call GeneralAuxVarCompute(xxbc,gen_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
      ! update state and update aux var; this could result in two update to 
      ! the aux var as update state updates if the state changes
      call GeneralAuxVarUpdateState(xxbc,gen_auxvars_bc(sum_connection), &
                                    global_auxvars_bc(sum_connection), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%sat_func_id(ghosted_id))%ptr, &
                                    natural_id,option)
#ifdef DEBUG_GENERAL_FILEOUTPUT
      if (debug_flag > 0) then
        write(debug_unit,'(a,i5,i3,7es24.15)') 'bc_auxvar:', natural_id, &
                           global_auxvars_bc(ghosted_id)%istate, &
                            xxbc(:)
      endif
#endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

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
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
                          
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

  !Heeho initialize dynamic accumulation term for every p iteration step
  if (general_tough2_conv_criteria) then
    call VecGetArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! GENERAL_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = GENERAL_UPDATE_FOR_FIXED_ACCUM
    call GeneralAuxVarCompute(xx_p(local_start:local_end), &
                              gen_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              natural_id, &
                              option)
    call GeneralAccumulation(gen_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,accum_p(local_start:local_end), &
                             Jac_dummy,PETSC_FALSE, &
                             local_id == general_debug_cell_id) 
  enddo
  
  if (general_tough2_conv_criteria) then
    accum_p2 = accum_p
  endif
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  !Heeho initialize dynamic accumulation term for every p iteration step
  if (general_tough2_conv_criteria) then
    call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
end subroutine GeneralUpdateFixedAccum

! ************************************************************************** !

subroutine GeneralAccumulation(gen_auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res,Jac, &
                               analytical_derivatives,debug_cell)
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
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_cell
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: porosity
  PetscReal :: volume_over_dt
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use gen_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = gen_auxvar%effective_porosity
  
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
                            porosity * volume_over_dt

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
                    soil_heat_capacity * gen_auxvar%temp) * volume_over_dt
  
  if (analytical_derivatives) then
    Jac = 0.d0
    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        ! satl = 1
        ! ----------
        ! Water Equation
        ! por * satl * denl * Xwl
        ! ---
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Xwl + 
        ! por * ddenl_dpl * Xwl
        Jac(1,1) = &
          gen_auxvar%d%por_pl * gen_auxvar%den(1) * gen_auxvar%xmol(1,1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%xmol(1,1)
        ! w/respect to air mole fraction
        ! liquid phase density is indepenent of air mole fraction
        ! por * denl * dXwl_dXal
        ! Xwl = 1. - Xal
        ! dXwl_dXal = -1.
        Jac(1,2) = porosity * gen_auxvar%den(1) * (-1.d0)
        ! w/repect to temperature
        ! por * ddenl_dT * Xwl
        Jac(1,3) = porosity * gen_auxvar%d%denl_T * gen_auxvar%xmol(1,1)
        ! ----------
        ! Air Equation
        ! por * satl * denl * Xal
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Xal + 
        ! por * ddenl_dpl * Xal
        Jac(2,1) = &
          gen_auxvar%d%por_pl * gen_auxvar%den(1) * gen_auxvar%xmol(2,1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%xmol(2,1)
        ! w/respect to air mole fraction
        Jac(2,2) = porosity * gen_auxvar%den(1)
        ! w/repect to temperature
        ! por * ddenl_dT * Xwl
        Jac(2,3) = porosity * gen_auxvar%d%denl_T * gen_auxvar%xmol(2,1)
        ! ----------
        ! Energy Equation
        ! por * satl * denl * Ul + (1-por) * dens * Cp * T
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Ul + 
        ! por * ddenl_dpl * Ul + 
        ! por * denl * dUl_dpl + 
        ! -dpor_dpl * dens * Cp * T
        Jac(3,1) = &
          gen_auxvar%d%por_pl * gen_auxvar%den(1) * gen_auxvar%U(1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%U(1) + &
          porosity * gen_auxvar%den(1) * gen_auxvar%d%Ul_pl + &
          (-1.d0) * gen_auxvar%d%por_pl * &
            material_auxvar%soil_particle_density * &
            soil_heat_capacity * gen_auxvar%temp
        ! w/respect to air mole fraction
        Jac(3,2) = 0.d0
        ! w/respect to temperature
        Jac(3,3) = &
          porosity * gen_auxvar%den(1) * gen_auxvar%d%Ul_T + &
          (1.d0 - porosity) * material_auxvar%soil_particle_density * &
            soil_heat_capacity
      case(GAS_STATE)
      case(TWO_PHASE_STATE)
!        if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
    end select
    Jac = Jac * volume_over_dt
  endif
  
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
                       option,v_darcy,Res,Jup,Jdn, &
                       analytical_derivatives, &
                       debug_connection)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  use Fracture_module
  use Klinkenberg_module
  
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
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: xmol(option%nflowspec)
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass
  PetscReal :: delta_X_whatever
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  
  PetscReal :: dummy_dperm_up, dummy_dperm_dn
  PetscReal :: temp_perm_up, temp_perm_dn

  PetscReal :: dden_up, dden_dn
  PetscReal :: dden_dden_kg_up, dden_dden_kg_dn
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  PetscReal :: dmobility_dpup, dmobility_dpdn
  PetscReal :: dmobility_dsatup, dmobility_dsatdn
  PetscReal :: dmobility_dTup, dmobility_dTdn
  PetscReal :: dmole_flux_dpup, dmole_flux_dpdn
  PetscReal :: dmole_flux_dTup, dmole_flux_dTdn
  PetscReal :: dv_darcy_dpup, dv_darcy_dpdn
  PetscReal :: dv_darcy_dTup, dv_darcy_dTdn
  PetscReal :: duH_dpup, duH_dpdn
  PetscReal :: duH_dTup, duH_dTdn
  PetscReal :: dxmol_up, dxmol_dn
   
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_up%fracture)) then
    call FracturePermEvaluate(material_auxvar_up,perm_up,temp_perm_up, &
                              dummy_dperm_up,dist)
    perm_up = temp_perm_up
  endif
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif
  
  if (associated(klinkenberg)) then
    perm_ave_over_dist(1) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
    temp_perm_up = klinkenberg%Evaluate(perm_up, &
                                         gen_auxvar_up%pres(option%gas_phase))
    temp_perm_dn = klinkenberg%Evaluate(perm_dn, &
                                         gen_auxvar_dn%pres(option%gas_phase))
    perm_ave_over_dist(2) = (temp_perm_up * temp_perm_dn) / &
                            (dist_up*temp_perm_dn + dist_dn*temp_perm_up)
  else
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
  endif
      
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
                                           gen_auxvar_dn%den_kg, &
                                           dden_up,dden_dn)
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
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = GeneralAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          gen_auxvar_up%den, &
                                          gen_auxvar_dn%den, &
                                          dden_up,dden_dn)
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
      if (general_harmonic_diff_density) then
        den_up = gen_auxvar_up%den(iphase)
        den_dn = gen_auxvar_dn%den(iphase)
      else
        ! we use upstream weighting when iphase is not equal, otherwise
        ! arithmetic with 50/50 weighting
        den_up = GeneralAverageDensity(iphase, &
                                       global_auxvar_up%istate, &
                                       global_auxvar_dn%istate, &
                                       gen_auxvar_up%den, &
                                       gen_auxvar_dn%den, &
                                       dden_up,dden_dn)
        ! by setting both equal, we avoid the harmonic weighting below
        den_dn = den_up
      endif
      stpd_up = sat_up*material_auxvar_up%tortuosity* &
                gen_auxvar_up%effective_porosity*den_up
      stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
                gen_auxvar_dn%effective_porosity*den_dn
      if (general_diffuse_xmol) then ! delta of mole fraction
        delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                     gen_auxvar_dn%xmol(air_comp_id,iphase)
        delta_X_whatever = delta_xmol
      else ! delta of mass fraction
        xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
        xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
        xmass_air_up = xmol_air_up*fmw_comp(2) / &
                   (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
        xmass_air_dn = xmol_air_dn*fmw_comp(2) / &
                   (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
        delta_xmass = xmass_air_up - xmass_air_dn
        delta_X_whatever = delta_xmass
      endif
      ! units = [mole/m^4 bulk]
      stpd_ave_over_dist = (stpd_up*stpd_dn)/(stpd_up*dist_dn+stpd_dn*dist_up)
      ! need to account for multiple phases
      tempreal = 1.d0
      ! Eq. 1.9b.  The gas density is added below
      if (general_temp_dep_gas_air_diff .and. &
          iphase == option%gas_phase) then
        temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
        pressure_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                              gen_auxvar_dn%pres(iphase))
        tempreal = ((temp_ave+273.15d0)/273.15d0)**1.8d0 * &
                    101325.d0 / pressure_ave
      endif
      ! units = mole/sec
      mole_flux = stpd_ave_over_dist * tempreal * &
                  general_parameter%diffusion_coefficient(iphase) * &
                  delta_X_whatever * area
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
  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux
! CONDUCTION
#endif

  if (analytical_derivatives) then
    Jup = 0.d0
    Jdn = 0.d0
    
    do iphase = 1, option%nphase
 
      if (gen_auxvar_up%mobility(iphase) + &
          gen_auxvar_dn%mobility(iphase) < eps) then
        cycle
      endif

      density_kg_ave = GeneralAverageDensity(iphase, &
                                             global_auxvar_up%istate, &
                                             global_auxvar_dn%istate, &
                                             gen_auxvar_up%den_kg, &
                                             gen_auxvar_dn%den_kg, &
                                             dden_dden_kg_up,dden_dden_kg_dn)
      gravity_term = density_kg_ave * dist_gravity
      delta_pressure = gen_auxvar_up%pres(iphase) - &
                       gen_auxvar_dn%pres(iphase) + &
                       gravity_term
                       
      ddelta_pressure_dpup = 1.d0 + dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_pl*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_pl*FMWH2O)
      ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_pl*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_pl*FMWH2O)
      ddelta_pressure_dTup = dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_T*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_T*FMWH2O)
      ddelta_pressure_dTdn = dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_T*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_T*FMWH2O)

      if (delta_pressure >= 0.D0) then
        mobility = gen_auxvar_up%mobility(iphase)
        xmol(:) = gen_auxvar_up%xmol(:,iphase)
        H_ave = gen_auxvar_up%H(iphase)
        uH = H_ave
        
        duH_dpup = gen_auxvar_up%d%Hl_pl
        duH_dpdn = 0.d0
        duH_dTup = gen_auxvar_up%d%Hl_T
        duH_dTdn = 0.d0
        dxmol_up = 1.d0
        dxmol_dn = 0.d0
        dmobility_dpup = gen_auxvar_up%d%mobilityl_pl
        dmobility_dsatup = gen_auxvar_up%d%mobilityl_sat
        dmobility_dTup = gen_auxvar_up%d%mobilityl_T
        dmobility_dpdn = 0.d0
        dmobility_dsatdn = 0.d0
        dmobility_dTdn = 0.d0        
      else
        mobility = gen_auxvar_dn%mobility(iphase)
        xmol(:) = gen_auxvar_dn%xmol(:,iphase)
        H_ave = gen_auxvar_dn%H(iphase)
        uH = H_ave

        duH_dpup = 0.d0
        duH_dTup = 0.d0
        duH_dTdn = gen_auxvar_dn%d%Hl_T
        dxmol_up = 0.d0
        dxmol_dn = 1.d0
        dmobility_dpup = 0.d0
        dmobility_dsatup = 0.d0
        dmobility_dTup = 0.d0
        dmobility_dpdn = gen_auxvar_dn%d%mobilityl_pl
        dmobility_dsatdn = gen_auxvar_dn%d%mobilityl_sat
        dmobility_dTdn = gen_auxvar_dn%d%mobilityl_T
      endif      

      if (mobility > floweps) then
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
        
        tempreal = perm_ave_over_dist(iphase)
        dv_darcy_dpup = tempreal * &
          (dmobility_dpup * delta_pressure + mobility * ddelta_pressure_dpup)
        dv_darcy_dTup = tempreal * &
          (dmobility_dTup * delta_pressure + mobility * ddelta_pressure_dTup)
        dv_darcy_dpdn = tempreal * &
          (dmobility_dpdn * delta_pressure + mobility * ddelta_pressure_dpdn)
        dv_darcy_dTdn = tempreal * &
          (dmobility_dTdn * delta_pressure + mobility * ddelta_pressure_dTdn)
        
        density_ave = GeneralAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            gen_auxvar_up%den, &
                                            gen_auxvar_dn%den, &
                                            dden_up,dden_dn)
        ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy(iphase) * area  
        ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
        !                             density_ave[kmol phase/m^3 phase]        
        mole_flux = q*density_ave
        ! Res[kmol total/sec]
!        do icomp = 1, option%nflowspec
!          ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
!          !                      xmol[kmol comp/kmol phase]
!          Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
!        enddo
!        Res(energy_id) = Res(energy_id) + mole_flux * uH

        select case(global_auxvar_up%istate)
          case(LIQUID_STATE)
            dmole_flux_dpup = &
              (dv_darcy_dpup * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_pl + &
                 dden_dn * gen_auxvar_dn%d%denl_pl))
            dmole_flux_dTup = &
              (dv_darcy_dTup * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_T + &
                 dden_dn * gen_auxvar_dn%d%denl_T))
            do icomp = 1, option%nflowspec
              Jup(icomp,1) = Jup(icomp,1) + dmole_flux_dpup * xmol(icomp)
              Jup(icomp,2) = Jup(icomp,2) + mole_flux * dxmol_up
              Jup(icomp,3) = Jup(icomp,3) + dmole_flux_dTup * xmol(icomp)
            enddo
            Jup(energy_id,1) = Jup(energy_id,1) + &
              (dmole_flux_dpup * uH + mole_flux * duH_dpup)
            Jup(energy_id,3) = Jup(energy_id,3) + &
              (dmole_flux_dTup * uH + mole_flux * duH_dTup)
          case(GAS_STATE)
          case(TWO_PHASE_STATE)
        end select
        select case(global_auxvar_dn%istate)
          case(LIQUID_STATE)
            dmole_flux_dpdn = &
              (dv_darcy_dpdn * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_pl + &
                 dden_dn * gen_auxvar_dn%d%denl_pl)) 
            dmole_flux_dTdn = &
              (dv_darcy_dTdn * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_T + &
                 dden_dn * gen_auxvar_dn%d%denl_T))
            do icomp = 1, option%nflowspec
              Jdn(icomp,1) = Jdn(icomp,1) + dmole_flux_dpdn * xmol(icomp)
              Jdn(icomp,2) = Jdn(icomp,2) + mole_flux * dxmol_dn
              Jdn(icomp,3) = Jdn(icomp,3) + dmole_flux_dTdn * xmol(icomp)
            enddo
            Jdn(energy_id,1) = Jdn(energy_id,1) + &
              (dmole_flux_dpdn * uH + mole_flux * duH_dpdn)
            Jdn(energy_id,3) = Jdn(energy_id,3) + &
              (dmole_flux_dTdn * uH + mole_flux * duH_dTdn)
          case(GAS_STATE)
          case(TWO_PHASE_STATE)
        end select        
      endif                   
    enddo  
    Jup = Jup * area
    Jdn = Jdn * area

    ! add in gas component diffusion in gas and liquid phases
    do iphase = 1, option%nphase
    
      if (iphase == LIQUID_PHASE) cycle
    
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
        if (general_harmonic_diff_density) then
          den_up = gen_auxvar_up%den(iphase)
          den_dn = gen_auxvar_dn%den(iphase)
        else
          ! we use upstream weighting when iphase is not equal, otherwise
          ! arithmetic with 50/50 weighting
          den_up = GeneralAverageDensity(iphase, &
                                         global_auxvar_up%istate, &
                                         global_auxvar_dn%istate, &
                                         gen_auxvar_up%den, &
                                         gen_auxvar_dn%den, &
                                         dden_up,dden_dn)
          ! by setting both equal, we avoid the harmonic weighting below
          den_dn = den_up
        endif
        stpd_up = sat_up*material_auxvar_up%tortuosity* &
                  gen_auxvar_up%effective_porosity*den_up
        stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
                  gen_auxvar_dn%effective_porosity*den_dn
        if (general_diffuse_xmol) then
          delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                       gen_auxvar_dn%xmol(air_comp_id,iphase)
          delta_X_whatever = delta_xmol
        else
          xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
          xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
          xmass_air_up = xmol_air_up*fmw_comp(2) / &
                     (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
          xmass_air_dn = xmol_air_dn*fmw_comp(2) / &
                     (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
          delta_xmass = xmass_air_up - xmass_air_dn
          delta_X_whatever = delta_xmass
        endif
        ! units = [mole/m^4 bulk]
        stpd_ave_over_dist = (stpd_up*stpd_dn)/(stpd_up*dist_dn+stpd_dn*dist_up)
        ! need to account for multiple phases
        tempreal = 1.d0
        ! Eq. 1.9b.  The gas density is added below
        if (general_temp_dep_gas_air_diff .and. &
            iphase == option%gas_phase) then
          temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
          pressure_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                                gen_auxvar_dn%pres(iphase))
          tempreal = ((temp_ave+273.15d0)/273.15d0)**1.8d0 * &
                      101325.d0 / pressure_ave
        endif
        ! units = mole/sec
        mole_flux = stpd_ave_over_dist * tempreal * &
                    general_parameter%diffusion_coefficient(iphase) * &
                    delta_X_whatever * area
        Res(wat_comp_id) = Res(wat_comp_id) - mole_flux
        Res(air_comp_id) = Res(air_comp_id) + mole_flux
      endif
    enddo
  ! DIFFUSION
  endif
  
#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
!    write(*,'(a,7es12.4)') 'in: ', adv_flux(:)*dist(1), diff_flux(:)*dist(1)
    write(*,'('' phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(2), gen_auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(2), gen_auxvar_dn%sat(2)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,2), gen_auxvar_dn%xmol(1,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,2)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,2), gen_auxvar_dn%xmol(2,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,2)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,2) + heat_flux)*1.d6
    write(*,'('' phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(1), gen_auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(1), gen_auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,1), gen_auxvar_dn%xmol(1,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,1)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,1), gen_auxvar_dn%xmol(2,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,1)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,1) + heat_flux)*1.d6
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux
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
                         option,v_darcy,Res,debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module                              
  use Material_Aux_class
  use Fracture_module
  use Klinkenberg_module
  
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
  PetscBool :: debug_connection
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  PetscInt :: bc_type
  PetscReal :: xmol(option%nflowspec)  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  PetscReal :: boundary_pressure
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass  
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: tempreal
  PetscReal :: delta_X_whatever

  PetscReal :: dden_dn, dden_up

  PetscInt :: idof
  PetscBool :: neumann_bc_present
  
  PetscReal :: temp_perm_dn
  PetscReal :: dummy_dperm_dn
  
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

  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif  
  
  if (associated(klinkenberg)) then
    perm_dn_adj(1) = perm_dn
                                          
    perm_dn_adj(2) = klinkenberg%Evaluate(perm_dn, &
                                          gen_auxvar_dn%pres(option%gas_phase))
  else
    perm_dn_adj(:) = perm_dn
  endif
  
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
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
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
                                                 gen_auxvar_dn%den_kg, &
                                                 dden_up, dden_dn)
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
                                                gen_auxvar_dn%den, &
                                                dden_up,dden_dn)
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
      if (density_ave < 1.d-40) then
        option%io_buffer = 'Zero density in GeneralBCFlux()'
        call printErrMsgByRank(option)
      endif
      mole_flux = q*density_ave       
      ! Res[kmol total/sec]
      do icomp = 1, option%nflowspec
        ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
        !                      xmol[kmol comp/mol phase]
        Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      enddo
#ifdef DEBUG_FLUXES  
      do icomp = 1, option%nflowspec
        adv_flux(icomp,iphase) = adv_flux(icomp,iphase) + mole_flux * xmol(icomp)
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
      adv_flux(energy_id,iphase) = adv_flux(energy_id,iphase) + mole_flux * uH
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
      if (general_harmonic_diff_density) then
        den_dn = gen_auxvar_dn%den(iphase)
      else
        ! we use upstream weighting when iphase is not equal, otherwise
        ! arithmetic with 50/50 weighting
        den_dn = GeneralAverageDensity(iphase, &
                                       global_auxvar_up%istate, &
                                       global_auxvar_dn%istate, &
                                       gen_auxvar_up%den, &
                                       gen_auxvar_dn%den, &
                                       dden_up,dden_dn)
      endif
      ! units = [mole/m^4 bulk]
      stpd_ave_over_dist = sat_dn*material_auxvar_dn%tortuosity * &
                           gen_auxvar_dn%effective_porosity * &
                           den_dn / dist(0)
      if (general_diffuse_xmol) then
        delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                     gen_auxvar_dn%xmol(air_comp_id,iphase)
        delta_X_whatever = delta_xmol
      else
        xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
        xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
        xmass_air_up = xmol_air_up*fmw_comp(2) / &
                   (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
        xmass_air_dn = xmol_air_dn*fmw_comp(2) / &
                   (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
        delta_xmass = xmass_air_up - xmass_air_dn
        delta_X_whatever = delta_xmass
      endif
      ! need to account for multiple phases
      ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
      tempreal = 1.d0
      ! Eq. 1.9b.  The gas density is added below
      if (general_temp_dep_gas_air_diff .and. &
          iphase == option%gas_phase) then
        temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
        pres_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                          gen_auxvar_dn%pres(iphase))
        tempreal = ((temp_ave+273.15d0)/273.15d0)**1.8d0 * &
                    101325.d0 / pres_ave
      endif
      ! units = mole/sec
      mole_flux = stpd_ave_over_dist * tempreal * &
                  general_parameter%diffusion_coefficient(iphase) * &
                  delta_X_whatever * area
      Res(wat_comp_id) = Res(wat_comp_id) - mole_flux
      Res(air_comp_id) = Res(air_comp_id) + mole_flux
#ifdef DEBUG_FLUXES  
      ! equal but opposite
      diff_flux(wat_comp_id,iphase) = diff_flux(wat_comp_id,iphase) - mole_flux
      diff_flux(air_comp_id,iphase) = diff_flux(air_comp_id,iphase) + mole_flux
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
  heat_flux = 0.d0
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
      heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! convert W -> MW
    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = auxvars(auxvar_mapping(GENERAL_ENERGY_FLUX_INDEX)) * area

    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'GeneralBCFlux heat conduction loop.'
      call printErrMsg(option)
  end select
  Res(energy_id) = Res(energy_id) + heat_flux ! MW
! CONDUCTION
#endif

#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
!    write(*,'(a,7es12.4)') 'in: ', adv_flux(:)*dist(1), diff_flux(:)*dist(1)
    write(*,'('' phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(2), gen_auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(2), gen_auxvar_dn%sat(2)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,2), gen_auxvar_dn%xmol(1,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,2)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,2), gen_auxvar_dn%xmol(2,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,2)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,2) + heat_flux)*1.d6
    write(*,'('' phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(1), gen_auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(1), gen_auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,1), gen_auxvar_dn%xmol(1,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,1)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,1), gen_auxvar_dn%xmol(2,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,1)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,1) + heat_flux)*1.d6
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'bc dif flux (liquid):', debug_flux(:,1)*dist(3)
    write(debug_unit,'(a,7es24.15)') 'bc dif flux (gas):', debug_flux(:,2)*dist(3)
  endif
#endif
  
end subroutine GeneralBCFlux

! ************************************************************************** !

subroutine GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                          gen_auxvar,global_auxvar,ss_flow_vol_flux, &
                          scale,Res,debug_cell)
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
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
  PetscBool :: debug_cell
      
  PetscReal :: qsrc_mol
  PetscReal :: enthalpy, internal_energy
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: icomp
  PetscErrorCode :: ierr

  Res = 0.d0
  do icomp = 1, option%nflowspec
    qsrc_mol = 0.d0
    select case(flow_src_sink_type)
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(icomp)/fmw_comp(icomp) ! kg/sec -> kmol/sec
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_mol = qsrc(icomp)/fmw_comp(icomp)*scale 
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec
        qsrc_mol = qsrc(icomp)*gen_auxvar%den(icomp) ! den = kmol/m^3
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_mol = qsrc(icomp)*gen_auxvar%den(icomp)*scale 
    end select
    ! icomp here is really iphase
    ss_flow_vol_flux(icomp) = qsrc_mol/gen_auxvar%den(icomp)
    Res(icomp) = qsrc_mol
  enddo
  if (dabs(qsrc(TWO_INTEGER)) < 1.d-40 .and. &
      qsrc(ONE_INTEGER) < 0.d0) then ! extraction only
    ! Res(1) holds qsrc_mol for water.  If the src/sink value for air is zero,
    ! remove/add the equivalent mole fraction of air in the liquid phase.
    qsrc_mol = Res(ONE_INTEGER)*gen_auxvar%xmol(TWO_INTEGER,ONE_INTEGER)
    Res(TWO_INTEGER) = qsrc_mol
    ss_flow_vol_flux(TWO_INTEGER) = qsrc_mol/gen_auxvar%den(TWO_INTEGER)
  endif
  ! energy units: MJ/sec
  if (size(qsrc) == THREE_INTEGER) then
    if (dabs(qsrc(THREE_INTEGER)) < 1.d-40) then
      cell_pressure = &
        maxval(gen_auxvar%pres(option%liquid_phase:option%gas_phase))
      if (dabs(qsrc(ONE_INTEGER)) > 0.d0) then
        call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,enthalpy,ierr)
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
        ! enthalpy units: MJ/kmol                       ! water component mass
        Res(option%energy_id) = Res(option%energy_id) + Res(ONE_INTEGER) * &
                                                        enthalpy
      endif
      if (dabs(qsrc(TWO_INTEGER)) > 0.d0) then
        ! this is pure air, we use the enthalpy of air, NOT the air/water
        ! mixture in gas
        ! air enthalpy is only a function of temperature and the 
        dummy_pressure = 0.d0
        call EOSGasEnergy(gen_auxvar%temp,dummy_pressure, &
                          enthalpy,internal_energy,ierr)
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> MJ/kmol                                  
        ! enthalpy units: MJ/kmol                       ! air component mass
        Res(option%energy_id) = Res(option%energy_id) + Res(TWO_INTEGER) * &
                                                        enthalpy
      endif
    else
      Res(option%energy_id) = qsrc(THREE_INTEGER)*scale ! MJ/s
    endif
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'src/sink:', Res(1)-Res(2),Res(12:3)
  endif
#endif   
  
end subroutine GeneralSrcSink

! ************************************************************************** !

subroutine GeneralAccumDerivative(gen_auxvar,global_auxvar,material_auxvar, &
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
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: jac(option%nflowdof,option%nflowdof)
  PetscReal :: jac_pert(option%nflowdof,option%nflowdof)
  PetscInt :: idof, irow

!geh:print *, 'GeneralAccumDerivative'

  call GeneralAccumulation(gen_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar,soil_heat_capacity,option, &
                           res,jac,general_analytical_derivatives, &
                           PETSC_FALSE)
                           
  do idof = 1, option%nflowdof
    call GeneralAccumulation(gen_auxvar(idof), &
                             global_auxvar, &
                             material_auxvar,soil_heat_capacity, &
                             option,res_pert,jac_pert,PETSC_FALSE,PETSC_FALSE)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar(idof)%pert
!geh:print *, irow, idof, J(irow,idof), gen_auxvar(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_analytical_derivatives) then
    J = jac
  endif

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
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_up(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_dn(option%nflowdof,option%nflowdof)
  PetscReal :: Jdummy(option%nflowdof,option%nflowdof)

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
                   option,v_darcy,res,Janal_up,Janal_dn,&
                   general_analytical_derivatives,PETSC_FALSE)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_auxvar_up(idof),global_auxvar_up, &
                     material_auxvar_up,sir_up, &
                     thermal_conductivity_up, &
                     gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn,sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res_pert,Jdummy,Jdummy, &
                     PETSC_FALSE,PETSC_FALSE)
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
                     option,v_darcy,res_pert,Jdummy,Jdummy, &
                     PETSC_FALSE,PETSC_FALSE)
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
                     option,v_darcy,res,PETSC_FALSE)                     
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                       gen_auxvar_up,global_auxvar_up, &
                       gen_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area,dist,general_parameter, &
                       option,v_darcy,res_pert,PETSC_FALSE)   
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
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  option%iflag = -3
  call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                      gen_auxvars(ZERO_INTEGER),global_auxvar,dummy_real, &
                      scale,res,PETSC_FALSE)
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                        gen_auxvars(idof),global_auxvar,dummy_real, &
                        scale,res_pert,PETSC_FALSE)            
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

!#define DEBUG_WITH_TECPLOT
#ifdef DEBUG_WITH_TECPLOT
  use Output_Tecplot_module
#endif

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
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
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
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
  
                                             ! do update state
  call GeneralUpdateAuxVars(realization,PETSC_TRUE)

! for debugging a single grid cell
!  i = 6
!  call GeneralOutputAuxVars(gen_auxvars(0,i),global_auxvars(i),i,'genaux', &
!                            PETSC_TRUE,option)
#ifdef DEBUG_WITH_TECPLOT
! for debugging entire solution over a single SNES solve
  write(word,*) iplot
  iplot = iplot + 1
  realization%output_option%plot_name = 'general-ni-' // trim(adjustl(word))
  call OutputTecplotPoint(realization)
#endif

  ! override flags since they will soon be out of date
  patch%aux%General%auxvars_up_to_date = PETSC_FALSE 

  ! always assume variables have been swapped; therefore, must copy back
  call VecLockPop(xx,ierr); CHKERRQ(ierr)
  call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                   NFLOWDOF)
  call VecLockPush(xx,ierr); CHKERRQ(ierr)

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

  
  !Heeho dynamically update p+1 accumulation term
  if (general_tough2_conv_criteria) then
    call VecGetArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
  ! accumulation at t(k+1)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call GeneralAccumulation(gen_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,Res,Jac_dummy, &
                             general_analytical_derivatives, &
                             local_id == general_debug_cell_id) 
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    
    !Heeho dynamically update p+1 accumulation term
    if (general_tough2_conv_criteria) then
      accum_p2(local_start:local_end) = Res(:)
    endif
    
  enddo

  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  !Heeho dynamically update p+1 accumulation term
  if (general_tough2_conv_criteria) then
    call VecRestoreArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
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
                       general_parameter,option,v_darcy,Res, &
                       Jac_dummy,Jac_dummy, &
                       general_analytical_derivatives, &
                       (local_id_up == general_debug_cell_id .or. &
                        local_id_dn == general_debug_cell_id))

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
                     v_darcy,Res, &
                     local_id == general_debug_cell_id)
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
      
      call GeneralSrcSink(option,source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                        source_sink%flow_condition%general%rate%itype, &
                        gen_auxvars(ZERO_INTEGER,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        ss_flow_vol_flux, &
                        scale,Res, &
                        local_id == general_debug_cell_id)

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

  if (patch%aux%General%inactive_cells_exist) then
    do i=1,patch%aux%General%n_inactive_rows
      r_p(patch%aux%General%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  call GeneralSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                        gen_auxvars,option)

  if (Initialized(general_debug_cell_id)) then
    call VecGetArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
    do local_id = general_debug_cell_id-1, general_debug_cell_id+1
      write(*,'(''  residual   : '',i2,10es12.4)') local_id, &
        r_p((local_id-1)*option%nflowdof+1:(local_id-1)*option%nflowdof+2), &
        r_p(local_id*option%nflowdof)*1.d6
    enddo
    call VecRestoreArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  
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
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'residual:', local_id, &
      r_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
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
    natural_id = grid%nG2A(ghosted_id)
    call GeneralAuxVarPerturb(gen_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              natural_id,option)
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
    call GeneralNumericalJacobianTest(xx,realization,J) 
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

function GeneralGetTecplotHeader(realization,icolumn)
  ! 
  ! Returns General Lite contribution to
  ! Tecplot file header
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  
  use Realization_Subsurface_class
  use Option_module
  use Field_module
    
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: GeneralGetTecplotHeader
  type(realization_subsurface_type) :: realization
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

subroutine GeneralSetPlotVariables(realization,list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/13
  ! 
  
  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(realization_subsurface_type) :: realization
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
  call OutputVariableAddToList(list,output_variable)   
  
end subroutine GeneralSetPlotVariables

! ************************************************************************** !

function GeneralAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn,dden_up,dden_dn)
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
  PetscReal :: dden_up, dden_dn

  PetscReal :: GeneralAverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0
  if (iphase == LIQUID_PHASE) then
    if (istate_up == GAS_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0
    else if (istate_dn == GAS_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == LIQUID_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0      
    else if (istate_dn == LIQUID_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0      
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0      
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
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"

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
  PetscInt :: local_id, ghosted_id, istart, iend, irow, idof
  PetscReal :: res_pert(option%nflowdof)
  PetscReal :: aux_real(10)
  PetscErrorCode :: ierr
  
  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif
  
  cur_srcsink => ss_sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    aux_real = 0.d0
    local_id = cur_srcsink%local_cell_id
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

subroutine GeneralMapBCAuxVarsToGlobal(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
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
  type(general_auxvar_type), pointer :: gen_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update
  
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
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
  
end subroutine GeneralMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine GeneralDestroy(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine GeneralDestroy

end module General_module
