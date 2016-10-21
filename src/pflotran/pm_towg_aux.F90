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
  !PetscInt, parameter, public :: TOIL_IMS_PRESSURE_DOF = 1
  !PetscInt, parameter, public :: TOIL_IMS_SATURATION_DOF = 2
  !PetscInt, parameter, public :: TOIL_IMS_ENERGY_DOF = 3

  ! Equation indices -   
  !PetscInt, parameter, public :: TOIL_IMS_LIQUID_EQUATION_INDEX = 1
  !PetscInt, parameter, public :: TOIL_IMS_OIL_EQUATION_INDEX = 2
  !PetscInt, parameter, public :: TOIL_IMS_ENERGY_EQUATION_INDEX = 3

  ! Indices used to map aux_real for condition values 
  !PetscInt, parameter, public :: TOIL_IMS_PRESSURE_INDEX = 1
  !PetscInt, parameter, public :: TOIL_IMS_OIL_SATURATION_INDEX = 2
  !PetscInt, parameter, public :: TOIL_IMS_TEMPERATURE_INDEX = 3
  !PetscInt, parameter, public :: TOIL_IMS_LIQUID_FLUX_INDEX = 4
  !PetscInt, parameter, public :: TOIL_IMS_OIL_FLUX_INDEX = 5
  !PetscInt, parameter, public :: TOIL_IMS_ENERGY_FLUX_INDEX = 6
  !PetscInt, parameter, public :: TOIL_IMS_LIQ_CONDUCTANCE_INDEX = 7
  !PetscInt, parameter, public :: TOIL_IMS_OIL_CONDUCTANCE_INDEX = 8
  !PetscInt, parameter, public :: TOIL_IMS_MAX_INDEX = 9

  !PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_DERIVATIVE = -1
  !PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_FIXED_ACCUM = 0
  !PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_ACCUM = 1
  !PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_BOUNDARY = 2

  !phase mapping:
  !While LIQUID_PHASE = 1 alway, OIL_PHASE can be:
  ! 2 for brine/oil
  ! 3 for brine/gas/oil (can be 2 again, but then GAS_PHASE requires
  !                      a different index). Need mapping in any case
  !PetscInt, parameter, public :: TOIL_IMS_OIL_PHASE = 2

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
