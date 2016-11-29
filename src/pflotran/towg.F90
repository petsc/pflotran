module TOWG_module

  use PM_TOWG_Aux_module
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

  public :: TOWGSetup

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

  ! use TOWGSetPlotVariables pointer function, pointing to different functions
  ! based on the miscibility model in use:
  !  TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL,TOWG_SOLVENT_TL
  !list => realization%output_option%output_snap_variable_list
  !call TOilImsSetPlotVariables(list)
  !list => realization%output_option%output_obs_variable_list
  !call TOilImsSetPlotVariables(list) 

  ! covergence creteria to be chosen (can use TOUGH or TOWG type) 
  ! set up here tough convergnce creteria if needed.

end subroutine TOWGSetup

! ************************************************************************** !




end module TOWG_module
