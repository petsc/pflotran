module PM_TOWG_Aux_module

  use PFLOTRAN_Constants_module

  use PM_Base_Aux_module
  use AuxVars_TOWG_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  !global variable to TOWG
  PetscReal, public :: towg_window_epsilon = 1.d-4
  !PetscReal, public :: towg_fmw_comp(3) = & ! initialised after EOSread
  !           [UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE]
  PetscReal, pointer, public :: towg_fmw_comp(:)
  PetscInt, pointer, public :: towg_dof_to_primary_variable(:,:)

  PetscInt, public :: towg_debug_cell_id = UNINITIALIZED_INTEGER
  PetscReal, parameter, public :: towg_pressure_scale = 1.d0

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
  PetscInt, parameter, public :: TOWG_X_GAS_IN_GAS_INDEX = 6
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
    procedure, public :: Init => InitTOWGAuxVars
    !procedure, public :: Perturb => PerturbTOilIms
  end type pm_towg_aux_type

  interface TOWGAuxVarStrip
    module procedure TOWGAuxVarArray1Strip
    module procedure TOWGAuxVarArray2Strip
  end interface TOWGAuxVarStrip
 
  public :: TOWGAuxCreate, TOWGAuxDestroy, TOWGAuxVarStrip
  !           , TOilImsAuxVarCompute, &
  !          TOilImsAuxVarPerturb, TOilImsAuxDestroy, &
  !          TOilImsAuxVarStrip

contains

! ************************************************************************** !

function TOWGAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object for TOWG
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/05/16
  ! 

  use Option_module
  use EOS_Oil_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
    
  class(pm_towg_aux_type), pointer :: TOWGAuxCreate

  class(pm_towg_aux_type), pointer :: aux

  allocate(towg_fmw_comp(option%nflowspec)) 
  
  towg_fmw_comp(1) = FMWH2O
  towg_fmw_comp(2) = EOSOilGetFMW()
  towg_fmw_comp(3) = EOSGasGetFMW()
  !in TOWG the gas FMW must be defined in the input deck
  if ( Initialized(EOSGasGetFMW()) ) then
    option%io_buffer = 'TOWG: gas FMW not initialised. ' // &
                       'Define its value in the the input deck' 
    call printErrMsg(option)    
  endif
  !need to add an EOS for solvent, before adding solvent 
  !if ( towg_miscibility_model == TOWG_SOLVENT_TL ) then
  !  towg_fmw_comp(4) = 
  !end if   

  allocate( towg_dof_to_primary_variable(1:option%nflowdof,1:3) )

  select case(towg_miscibility_model)
    case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL)
      towg_dof_to_primary_variable(1:option%nflowdof,1:3) = &
        reshape([TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_X_GAS_IN_OIL_INDEX, TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_GAS_SATURATION_INDEX, &
                TOWG_X_GAS_IN_GAS_INDEX, TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_GAS_SATURATION_INDEX,TOWG_TEMPERATURE_INDEX], &
                shape(towg_dof_to_primary_variable))
    case(TOWG_SOLVENT_TL)
      towg_dof_to_primary_variable(1:option%nflowdof,1:3) = &
        reshape([TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_X_GAS_IN_OIL_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_GAS_SATURATION_INDEX, &
                TOWG_X_GAS_IN_GAS_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_GAS_SATURATION_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX], &
                shape(towg_dof_to_primary_variable))
  end select

  !allocate here to define this is a pm_toil_ims_aux_type
  allocate(aux) 

  call PMBaseAuxInit(aux) 

  !nullify here and not in the parent class, because auxvars are mode dependent
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  allocate(aux%parameter)
  
  !PO - to allocate when supporting diffusion the oil and gas phase
  !allocate(aux%parameter%diffusion_coefficient(option%nphase)) 
  !aux%parameter%diffusion_coefficient(option%liquid_phase) = & 
  !                                                 UNINITIALIZED_DOUBLE
  !aux%parameter%diffusion_coefficient(option%oil_phase) = &
  !                                                 UNINITIALIZED_DOUBLE
  !aux%parameter%diffusion_coefficient(option%gas_phase) = 2.13d-5
  !aux%parameter%newton_inf_scaled_res_tol = 1.d-50
  !aux%parameter%check_post_converged = PETSC_FALSE


  TOWGAuxCreate => aux
  
end function TOWGAuxCreate

! ************************************************************************** !

subroutine InitTOWGAuxVars(this,grid,num_bc_connection, &
                              num_ss_connection,option)
  ! 
  ! Initialize pm_towg_auxvars  
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/07/16
  ! 

  use Option_module
  use Grid_module 

  implicit none

  class(pm_towg_aux_type) :: this
  PetscInt :: num_bc_connection
  PetscInt :: num_ss_connection  
  type(grid_type) :: grid
  type(option_type) :: option

  PetscInt :: ghosted_id, iconn, local_id
  PetscInt :: idof 
 
  allocate(this%auxvars(0:option%nflowdof,grid%ngmax)) 
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      call this%auxvars(idof,ghosted_id)%Init(option)
    enddo
  enddo

  this%num_aux = grid%ngmax

  if (num_bc_connection > 0) then
    allocate(this%auxvars_bc(num_bc_connection))
    do iconn = 1, num_bc_connection
      call this%auxvars_bc(iconn)%Init(option)
    enddo
  endif
  this%num_aux_bc = num_bc_connection

  if (num_ss_connection > 0) then
    allocate(this%auxvars_ss(num_ss_connection))
    do iconn = 1, num_ss_connection
      call this%auxvars_ss(iconn)%Init(option)
    enddo
  endif
  this%num_aux_ss = num_ss_connection

  call PMBaseAuxSetup(this,grid,option)

end subroutine InitTOWGAuxVars

! ************************************************************************** !

subroutine TOWGAuxDestroy(aux)
  ! 
  ! Deallocates a towg auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_towg_aux_type), pointer :: aux
  PetscInt :: iaux, idof
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars) ) then
    call TOWGAuxVarStrip(aux%auxvars)
    deallocate(aux%auxvars)
  end if 
  nullify(aux%auxvars) 

  if (associated(aux%auxvars_bc) ) then
    call TOWGAuxVarStrip(aux%auxvars_bc)
    deallocate(aux%auxvars_bc)
  end if
  nullify(aux%auxvars_bc)

  if ( associated(aux%auxvars_ss) ) then
    call TOWGAuxVarStrip(aux%auxvars_ss)
    deallocate(aux%auxvars_ss)
  end if
  nullify(aux%auxvars_ss)

  call PMBaseAuxStrip(aux)

  if (associated(aux%parameter)) then
    deallocate(aux%parameter) 
    !to add paramter strip when introducing diff in oil and gas phases
  end if
  nullify(aux%parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine TOWGAuxDestroy

! ************************************************************************** !

! ************************************************************************** !

subroutine  TOWGAuxVarArray1Strip(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes 
  ! using class(*) (unlimited polymorphic)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  ! 

  use AuxVars_TOWG_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !bug fixed in gfortran 6.2.0
  !here we can pass by pointer, we could destroy the array within the routine
  !but we don't to be consistent with TOilImsAuxVarArray2Strip 
  !class(auxvar_towg_type), pointer :: auxvars(:)
  type(auxvar_towg_type) :: auxvars(:)

  PetscInt :: iaux

  !print *, "den oil bc/ss pass = ", auxvars(1)%den(2)
  
  do iaux = 1, size(auxvars)
    call auxvars(iaux)%Strip
  enddo  

end subroutine TOWGAuxVarArray1Strip

! ************************************************************************** !

subroutine TOWGAuxVarArray2Strip(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes 
  ! using class(*) (unlimited polymorphic)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  ! 

  use AuxVars_TOWG_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !bug fixed in gfortran 6.2.0
  !cannot use type(...) with pointer attribute, therefore we deallocate and 
  !nullify pointer outide this routine 
  !because the compiler does not allow to specify lower 0-bound in auxvar
  !type(auxvar_towg_type), pointer :: auxvars(:,:)
  !class(auxvar_towg_type) :: auxvars(0:,:)
  type(auxvar_towg_type) :: auxvars(0:,:)

  PetscInt :: iaux, idof

  do iaux = 1, size(auxvars,2)
    do idof = 1, size(auxvars,1)
      call auxvars(idof-1,iaux)%Strip()
    enddo
  enddo  

end subroutine TOWGAuxVarArray2Strip

! ************************************************************************** !

end module PM_TOWG_Aux_module
