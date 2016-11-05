module PM_TOWG_Aux_module

  use PFLOTRAN_Constants_module

  use PM_Base_Aux_module
  use AuxVars_TOWG_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: towg_window_epsilon = 1.d-4
  PetscReal, public :: towg_fmw_comp(3) = & ! initialised after EOSread
             [UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE] 
  
  PetscInt, public :: towg_debug_cell_id = UNINITIALIZED_INTEGER
  PetscReal, parameter, public :: towg_pressure_scale = 1.d0
  ! these variables, which are global to general, can be modified
  !PetscInt, public :: toil_ims_dof_to_primary_vars(3) ! only one 2ph state 
  PetscBool, public :: towg_isothermal = PETSC_FALSE
  PetscInt, public :: towg_miscibility_model = UNINITIALIZED_INTEGER

  !list of TOWG paramters
  !PetscInt, parameter, public :: TOWG_PREV_TS = 1
  !PetscInt, parameter, public :: TOWG_PREV_IT = 2

  !available miscibility models
  PetscInt, parameter, public :: TOWG_IMMISCIBLE = 1
  PetscInt, parameter, public :: TOWG_TODD_LONGSTAFF = 2
  PetscInt, parameter, public :: TOWG_BLACK_OIL = 3
  PetscInt, parameter, public :: TOWG_SOLVENT_TL = 4

  ! thermodynamic state of fluid ids - for BLACK OIL 
  PetscInt, parameter, public :: TOWG_NULL_STATE = 0
  PetscInt, parameter, public :: TOWG_LIQ_OIL_STATE = 1
  PetscInt, parameter, public :: TOWG_LIQ_GAS_STATE = 2
  PetscInt, parameter, public :: TOWG_THREE_PHASE_STATE = 3
  PetscInt, parameter, public :: TOWG_ANY_STATE = 4

  ! Primary DOF indices - 
  PetscInt, parameter, public :: TOWG_OIL_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TOWG_OIL_SATURATION_DOF = 2
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_DOF = 2
  PetscInt, parameter, public :: TOWG_X_GAS_IN_OIL_DOF = 3
  PetscInt, parameter, public :: TOWG_X_OIL_IN_GAS_DOF = 3
  PetscInt, parameter, public :: TOWG_SOLV_SATURATION = 4
  PetscInt, parameter, public :: TOWG_3CMPS_ENERGY_DOF = 4
  PetscInt, parameter, public :: TOWG_SOLV_TL_ENERGY_DOF = 5 
  !towg_energy_dof assigned TOWG_3CMPS_ENERGY_DOF or TOWG_SOLV_TL_ENERGY_DOF
  !in PMTOWGCreate 
  PetscInt, public :: towg_energy_dof = UNINITIALIZED_INTEGER 
  

  ! Equation indices -   
  PetscInt, parameter, public :: TOWG_LIQ_EQ_IDX = 1
  PetscInt, parameter, public :: TOWG_OIL_EQ_IDX = 2
  PetscInt, parameter, public :: TOWG_GAS_EQ_IDX = 3
  PetscInt, parameter, public :: TOWG_SOLV_EQ_IDX = 4
  PetscInt, parameter, public :: TOWG_3CMPS_ENERGY_EQ_IDX = 4
  PetscInt, parameter, public :: TOWG_SOLV_TL_ENERGY_EQ_IDX = 5
  PetscInt, public :: towg_energy_eq_idx = UNINITIALIZED_INTEGER

  !Indices used to map aux_real for condition values 
  PetscInt, parameter, public :: TOWG_OIL_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: TOWG_OIL_SATURATION_INDEX = 2 
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_INDEX = 3  
  PetscInt, parameter, public :: TOWG_SOLV_SATURATION_INDEX = 4
  PetscInt, parameter, public :: TOWG_X_GAS_IN_OIL_INDEX = 5
  PetscInt, parameter, public :: TOWG_X_OIL_IN_GAS_INDEX = 6
  PetscInt, parameter, public :: TOWG_TEMPERATURE_INDEX = 7
  PetscInt, parameter, public :: TOWG_LIQUID_FLUX_INDEX = 8
  PetscInt, parameter, public :: TOWG_OIL_FLUX_INDEX = 9
  PetscInt, parameter, public :: TOWG_GAS_FLUX_INDEX = 10
  PetscInt, parameter, public :: TOWG_SOLV_FLUX_INDEX = 11
  PetscInt, parameter, public :: TOWG_ENERGY_FLUX_INDEX = 12
  PetscInt, parameter, public :: TOWG_LIQ_CONDUCTANCE_INDEX = 13
  PetscInt, parameter, public :: TOWG_OIL_CONDUCTANCE_INDEX = 14
  PetscInt, parameter, public :: TOWG_GAS_CONDUCTANCE_INDEX = 15
  PetscInt, parameter, public :: TOWG_MAX_INDEX = 15

  !Indices used to map aux_int_var for condition values
  PetscInt, parameter, public :: TOWG_STATE_INDEX = 1

  !PetscInt, parameter, public :: TOWG_UPDATE_FOR_DERIVATIVE = -1
  !PetscInt, parameter, public :: TOWG_UPDATE_FOR_FIXED_ACCUM = 0
  !PetscInt, parameter, public :: TOWG_UPDATE_FOR_ACCUM = 1
  !PetscInt, parameter, public :: TOWG_UPDATE_FOR_BOUNDARY = 2

  ! it might be required for thermal diffusion terms and tough conv criteria
  ! consider to pput here all 
  type, public :: towg_parameter_type
     !  PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
     !  PetscReal :: newton_inf_scaled_res_tol
     !  PetscBool :: check_post_converged
  end type towg_parameter_type

  !if required, could add other intermediate classes:
  type, public, extends(pm_base_aux_type) :: pm_towg_aux_type 
    type(towg_parameter_type), pointer :: parameter
    class(auxvar_towg_type), pointer :: auxvars(:,:)
    class(auxvar_towg_type), pointer :: auxvars_bc(:)
    class(auxvar_towg_type), pointer :: auxvars_ss(:)
  contains
    !add bound-procedure
    !procedure, public :: Init => InitTOilImsAuxVars
    !procedure, public :: Perturb => PerturbTOilIms
  end type pm_towg_aux_type

  !interface TOilImsAuxVarStrip
  !  module procedure TOilImsAuxVarArray1Strip
  !  module procedure TOilImsAuxVarArray2Strip
  !end interface TOilImsAuxVarStrip
 
  !public :: TOilImsAuxCreate, TOilImsAuxVarCompute, &
  !          TOilImsAuxVarPerturb, TOilImsAuxDestroy, &
  !          TOilImsAuxVarStrip

contains

! ************************************************************************** !


end module PM_TOWG_Aux_module
