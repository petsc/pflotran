module TOilIms_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  PetscReal, public :: toil_ims_window_epsilon = 1.d-4
  PetscReal, public :: toil_ims_fmw_comp(2) = [FMWH2O,FMWOIL]
  PetscReal, public :: toil_ims_max_pressure_change = 5.d4
  PetscInt, public :: toil_ims_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: toil_ims_damping_factor = 0.6d0
  PetscReal, public :: toil_ims_tgh2_itol_scld_res_e1 = 1.d-5
  PetscReal, public :: toil_ims_tgh2_itol_scld_res_e2 = 1.d0
  PetscBool, public :: toil_ims_tgh2_conv_criteria = PETSC_FALSE
  PetscInt, public :: toil_ims_debug_cell_id = UNINITIALIZED_INTEGER

  ! if needed must be specifi to toil_ims (e.g. TOIL_IMS_PREV_TS)
  !PetscInt, parameter, public :: TOIL_IMS_PREV_TS = 1
  !PetscInt, parameter, public :: TOIL_IMS_PREV_IT = 2

  ! Primary DOF indices 
  PetscInt, parameter, public :: TOIL_IMS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TOIL_IMS_SATURATION = 2
  PetscInt, parameter, public :: TOIL_IMS_TEMPERATURE = 3

  ! Equation indices   
  PetscInt, parameter, public :: TOIL_IMS_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: TOIL_IMS_OIL_EQUATION_INDEX = 2
  PetscInt, parameter, public :: TOIL_IMS_ENERGY_EQUATION_INDEX = 3

  ! Indices used to map aux_real for condition values 
  PetscInt, parameter, public :: TOIL_IMS_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: TOIL_IMS_OIL_SATURATION_INDEX = 2
  PetscInt, parameter, public :: TOIL_IMS_TEMPERATURE_INDEX = 3
  PetscInt, parameter, public :: TOIL_IMS_LIQUID_FLUX_INDEX = 4
  PetscInt, parameter, public :: TOIL_IMS_OIL_FLUX_INDEX = 5
  PetscInt, parameter, public :: TOIL_IMS_ENERGY_FLUX_INDEX = 6
  PetscInt, parameter, public :: TOIL_IMS_LIQUID_CONDUCTANCE_INDEX = 7
  PetscInt, parameter, public :: TOIL_IMS_OIL_CONDUCTANCE_INDEX = 8
  PetscInt, parameter, public :: TOIL_IMS_MAX_INDEX = 9

  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_BOUNDARY = 2

  PetscReal, parameter, public :: toil_ims_pressure_scale = 1.d0
  
  ! these variables, which are global to general, can be modified
  PetscInt, public :: toil_ims_dof_to_primary_vars(3) ! only one 2ph state 
  PetscBool, public :: toil_ims_isothermal = PETSC_FALSE


  type, public :: toil_ims_auxvar_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal :: temp
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
!    PetscReal, pointer :: dsat_dp(:,:)
!    PetscReal, pointer :: dden_dp(:,:)
!    PetscReal, pointer :: dsat_dt(:)
!    PetscReal, pointer :: dden_dt(:)
    PetscReal, pointer :: mobility(:) ! relative perm / kinematic viscosity
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: pert
!    PetscReal, pointer :: dmobility_dp(:)
  end type toil_ims_auxvar_type

  ! it might be required for thermal diffusion terms and tough conv criteria
  !type, public :: general_parameter_type
  !  PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
  !  PetscReal :: newton_inf_scaled_res_tol
  !  PetscBool :: check_post_converged
  !end type general_parameter_type

  type, public :: toil_ims_type
    PetscInt :: n_inactive_rows
    PetscInt, pointer :: inactive_rows_local(:), inactive_rows_local_ghosted(:)
    PetscInt, pointer :: row_zeroing_array(:)

    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    !type(general_parameter_type), pointer :: general_parameter
    type(toil_ims_auxvar_type), pointer :: auxvars(:,:)
    type(toil_ims_auxvar_type), pointer :: auxvars_bc(:)
    type(toil_ims_auxvar_type), pointer :: auxvars_ss(:)
  end type toil_ims_type

! uncomment as needed
!  interface TOilImsAuxVarDestroy
!    module procedure TOilImsAuxVarSingleDestroy
!    module procedure TOilImsAuxVarArray1Destroy
!    module procedure TOilImsAuxVarArray2Destroy
!  end interface TOilImsAuxVarDestroy
  
!  interface TOilImsOutputAuxVars
!    module procedure TOilImsOutputAuxVars1
!    module procedure TOilImsOutputAuxVars2
!  end interface TOilImsOutputAuxVars


contains

end module TOilIms_Aux_module

