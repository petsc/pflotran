module PM_TOWG_class

  use PM_Base_class
  use PM_Subsurface_Flow_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"

  type, public, extends(pm_subsurface_flow_type) :: pm_towg_type
    !this two vectors could be moved to pm_subsurface_flow_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
    ! A) and B) could be moved to pm_subsurface
    !A) used for truncation within CheckUpdatePre
    PetscReal :: trunc_max_pressure_change = 5.d4
    PetscInt ::  max_it_before_damping = UNINITIALIZED_INTEGER
    PetscReal :: damping_factor = 0.6d0
    !B) used for convergence criteria within CheckUpdatePost
    PetscReal :: itol_rel_update = UNINITIALIZED_DOUBLE
    PetscReal :: itol_scaled_res = 1.d-5
    !check might need different dimension 
    PetscReal :: tgh2_itol_scld_res_e1(3,3) = 1.d-5
    PetscReal :: tgh2_itol_scld_res_e2 = 1.d0
    PetscBool :: tough2_conv_criteria = PETSC_FALSE
    !PetscInt :: miscibility_model = UNINITIALIZED_INTEGER 
  contains
    procedure, public :: Read => PMTOWGRead
    !procedure, public :: InitializeRun => PMGeneralInitializeRun
    !procedure, public :: InitializeTimestep => PMGeneralInitializeTimestep
    !procedure, public :: Residual => PMGeneralResidual
    !procedure, public :: Jacobian => PMGeneralJacobian
    !procedure, public :: UpdateTimestep => PMGeneralUpdateTimestep
    !procedure, public :: PreSolve => PMGeneralPreSolve
    !procedure, public :: PostSolve => PMGeneralPostSolve
    !procedure, public :: CheckUpdatePre => PMGeneralCheckUpdatePre
    !procedure, public :: CheckUpdatePost => PMGeneralCheckUpdatePost
    !procedure, public :: TimeCut => PMGeneralTimeCut
    !procedure, public :: UpdateSolution => PMGeneralUpdateSolution
    !procedure, public :: UpdateAuxVars => PMGeneralUpdateAuxVars
    !procedure, public :: MaxChange => PMGeneralMaxChange
    !procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    !procedure, public :: InputRecord => PMGeneralInputRecord
    !procedure, public :: CheckpointBinary => PMGeneralCheckpointBinary
    !procedure, public :: RestartBinary => PMGeneralRestartBinary
    !procedure, public :: Destroy => PMGeneralDestroy
  end type pm_towg_type
  
  public :: PMTOWGCreate
  
contains

! ************************************************************************** !
function PMTOWGCreate(miscibility_model,option)
  ! 
  ! Creates TOWG process models shell
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/15/16
  ! 
  use Option_module
  use Input_Aux_module
  use PM_TOWG_Aux_module
  use Variables_module, only : OIL_PRESSURE, GAS_PRESSURE, &
                               OIL_SATURATION, GAS_SATURATION, &
                               OIL_MOLE_FRACTION, GAS_MOLE_FRACTION, &  
                               TEMPERATURE, SOLVENT_SATURATION
  implicit none
  
  class(pm_towg_type), pointer :: PMTOWGCreate

  character(len=MAXWORDLENGTH) :: miscibility_model
  type(option_type), pointer :: option

  class(pm_towg_type), pointer :: towg_pm
  
#ifdef TOWG_DEBUG  
  print *, 'PMTOWGCreate()'
#endif  

  allocate(towg_pm)
  allocate(towg_pm%max_change_ivar(8))
  towg_pm%max_change_ivar = [OIL_PRESSURE, GAS_PRESSURE, &
                             OIL_SATURATION, GAS_SATURATION, &
                             OIL_MOLE_FRACTION, GAS_MOLE_FRACTION, & 
                             TEMPERATURE, SOLVENT_SATURATION]

  allocate(towg_pm%max_change_isubvar(8))
                                       !2 = gas in xmol(gas_c,oil_phase)
  towg_pm%max_change_isubvar = [0,0,0,0,2,1,0,0]
                                         !1 = gas in xmol(gas_c,gas_phase)

  towg_pm%trunc_max_pressure_change = 5.d4
  towg_pm%max_it_before_damping = UNINITIALIZED_INTEGER
  towg_pm%damping_factor = 0.6d0
  towg_pm%itol_rel_update = UNINITIALIZED_DOUBLE
  towg_pm%itol_scaled_res = 1.d-5
  towg_pm%tgh2_itol_scld_res_e1(3,3) = 1.d-5
  towg_pm%tgh2_itol_scld_res_e2 = 1.d0
  towg_pm%tough2_conv_criteria = PETSC_FALSE

  select case(trim(miscibility_model)) 
    case('TOWG_IMMISCIBLE')
      towg_miscibility_model = TOWG_IMMISCIBLE
      towg_energy_dof = TOWG_3CMPS_ENERGY_DOF
      towg_energy_eq_idx = TOWG_3CMPS_ENERGY_EQ_IDX
    case('TODD_LONGOSTAFF','TOWG_MISCIBLE')
      towg_miscibility_model = TOWG_TODD_LONGSTAFF
      towg_energy_dof = TOWG_3CMPS_ENERGY_DOF
      towg_energy_eq_idx = TOWG_3CMPS_ENERGY_EQ_IDX
    case('BLACK_OIL')
      towg_miscibility_model = TOWG_BLACK_OIL
      towg_energy_dof = TOWG_3CMPS_ENERGY_DOF
      towg_energy_eq_idx = TOWG_3CMPS_ENERGY_EQ_IDX
    case('SOLVENT_TL')
      towg_miscibility_model = TOWG_SOLVENT_TL
      towg_energy_dof = TOWG_SOLV_TL_ENERGY_DOF
      towg_energy_eq_idx = TOWG_SOLV_TL_ENERGY_EQ_IDX
    case default
      call InputKeywordUnrecognized(miscibility_model, &
                         'TOWG MISCIBILITY_MODEL',option)
      !option%io_buffer = 'TOWG - miscibility_model not recognised' 
      !call printErrMsg(option)
  end select 


  !towg_pm%miscibility_model = UNINITIALIZED_INTEGER

  if (Uninitialized(towg_energy_dof)) then 
    option%io_buffer = 'towg_energy_dof not set up'
    call printErrMsg(option)
  end if  

  if (Uninitialized(towg_energy_eq_idx)) then 
    option%io_buffer = 'towg_energy_eq_idx not set up'
    call printErrMsg(option)
  end if  


  call PMSubsurfaceFlowCreate(towg_pm)
  towg_pm%name = 'PMTOWG'

  PMTOWGCreate => towg_pm
  
end function PMTOWGCreate

! ************************************************************************** !

subroutine PMTOWGRead(this,input)
  ! 
  ! Read TOWG options and set up miscibility, and functions to be use in 
  ! the Residual and Jacobian 
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/15/16
  !
  !use TOWG_module !to set up the functions to be used in R and J
  use PM_TOWG_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none
  
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword 
  class(pm_towg_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'TOWG Options'

  !this%miscibility_model = UNINITIALIZED_INTEGER
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found,option)    
    if (found) cycle
    
    select case(trim(keyword))
      !case('MISCIBILITY_MODEL')
      !  ! read a word, add a new select case, and detect the case
      !  call InputReadWord(input,option,word,PETSC_TRUE)
      !  call InputErrorMsg(input,option,'MODEL type','MISCIBILITY_MODEL')  
      !  call StringToUpper(word)
      !  select case(trim(word)) 
      !    case('IMMISICIBLE')
      !      towg_miscibility_model = TOWG_IMMISCIBLE
      !    case('TODD_LONGOSTAFF','MISCIBLE')
      !      towg_miscibility_model = TOWG_TODD_LONGSTAFF
      !    case('BLACK_OIL')
      !      towg_miscibility_model = TOWG_BLACK_OIL
      !    case('SOLVENT_TL')
      !      towg_miscibility_model = TOWG_SOLVENT_TL
      !    case default
      !      call InputKeywordUnrecognized(keyword,'MISCIBILITY_MODEL',option)
      !  end select 
      case('ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,this%itol_scaled_res)
        call InputDefaultMsg(input,option,'towg itol_scaled_res')
        this%check_post_convergence = PETSC_TRUE
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,this%itol_rel_update)
        call InputDefaultMsg(input,option,'towg itol_rel_update')
        this%check_post_convergence = PETSC_TRUE        
      case('TOUGH2_ITOL_SCALED_RESIDUAL')
        ! tgh2_itol_scld_res_e1 is an array: assign same value to all entries
        tempreal = UNINITIALIZED_DOUBLE
        call InputReadDouble(input,option,tempreal)
        ! tempreal will remain uninitialized if the read fails.
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e1')
        if (Initialized(tempreal)) then
          this%tgh2_itol_scld_res_e1 = tempreal
        endif
        call InputReadDouble(input,option,this%tgh2_itol_scld_res_e2)
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e2')
        this%tough2_conv_criteria = PETSC_TRUE
        this%check_post_convergence = PETSC_TRUE
      case('T2_ITOL_SCALED_RESIDUAL_TEMP')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option, &
                           'tough_itol_scaled_residual_e1 for temperature', &
                           error_string)
        this%tgh2_itol_scld_res_e1(3,:) = tempreal
      case('WINDOW_EPSILON') 
        call InputReadDouble(input,option,towg_window_epsilon)
        call InputErrorMsg(input,option,'towg window epsilon',error_string)
      ! read this in the gas EOS
      !case('GAS_COMPONENT_FORMULA_WEIGHT')
      !  !geh: assuming gas component is index 2
      !  call InputReadDouble(input,option,fmw_comp(2))
      !  call InputErrorMsg(input,option,'gas component formula wt.', &
      !                     error_string)
      case('ISOTHERMAL')
        towg_isothermal = PETSC_TRUE
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,this%trunc_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           error_string)
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,this%max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           error_string)
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,this%damping_factor)
        call InputErrorMsg(input,option,'damping factor',error_string)
      case('DEBUG_CELL')
        call InputReadInt(input,option,towg_debug_cell_id)
        call InputErrorMsg(input,option,'debug cell id',error_string)
      !case('DIFFUSE_XMASS')
      !  general_diffuse_xmol = PETSC_FALSE
      !case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
      !  general_harmonic_diff_density = PETSC_TRUE
      !case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
      !  general_harmonic_diff_density = PETSC_FALSE
      !case('ANALYTICAL_DERIVATIVES')
      !  general_analytical_derivatives = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,'TOWG Mode',option)
    end select
    
  enddo  

  !if (Uninitialized(towg_miscibility_model)) then 
  !  option%io_buffer = 'TOWG MISCIBILITY_MODEL not set up'
  !  call printErrMsg(option)
  !end if  

  !here set up functions in TOWG and pm_TOWG_aux based on miscibility model 
  !select case(towg_miscibility_model)
  !  case(TOWG_IMMISCIBLE) 
  !    !set up towg and pm_towg_aux functions   
  !  case default
  !    option%io_buffer = 'only immiscible TOWG currently implemented' 
  !    call printErrMsg(option)
  !end select

end subroutine PMTOWGRead

! ************************************************************************** !



end module PM_TOWG_class

