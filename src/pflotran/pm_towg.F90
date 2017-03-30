module PM_TOWG_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

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
    !procedure(TOWGMaxChangeDummy), pointer :: MaxChange => null()
    !procedure(MaxChange), pointer :: MaxChange => null()
  contains
    procedure, public :: Read => PMTOWGRead
    procedure, public :: InitializeRun => PMTOWGInitializeRun
    procedure, public :: InitializeTimestep => PMTOWGInitializeTimestep
    procedure, public :: Residual => PMTOWGResidual
    procedure, public :: Jacobian => PMTOWGJacobian
    procedure, public :: UpdateTimestep => PMTOWGUpdateTimestep
    procedure, public :: PreSolve => PMTOWGPreSolve
    !procedure, public :: PostSolve => PMTOWGPostSolve
    procedure, public :: CheckUpdatePre => PMTOWGCheckUpdatePre
    !procedure, public :: CheckUpdatePost => PMTOWGCheckUpdatePost
    procedure, public :: TimeCut => PMTOWGTimeCut
    procedure, public :: UpdateSolution => PMTOWGUpdateSolution
    procedure, public :: UpdateAuxVars => PMTOWGUpdateAuxVars
    procedure, public :: MaxChange => PMTOWGMaxChange
    procedure, public :: ComputeMassBalance => PMTOWGComputeMassBalance
    !procedure, public :: InputRecord => PMTOWGInputRecord
    procedure, public :: CheckpointBinary => PMTOWGCheckpointBinary
    procedure, public :: RestartBinary => PMTOWGRestartBinary
    procedure, public :: Destroy => PMTOWGDestroy
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

  select case(trim(miscibility_model))
    case('TOWG_IMMISCIBLE','TODD_LONGOSTAFF','TOWG_MISCIBLE')
      allocate(towg_pm%max_change_ivar(4))
      towg_pm%max_change_ivar = [OIL_PRESSURE, OIL_SATURATION, &
                                 GAS_SATURATION,TEMPERATURE]
      allocate(towg_pm%max_change_isubvar(4))
      towg_pm%max_change_isubvar = [0,0,0,0]
    case('BLACK_OIL')
      allocate(towg_pm%max_change_ivar(7))
      towg_pm%max_change_ivar = [OIL_PRESSURE, GAS_PRESSURE, &
                                 OIL_SATURATION, GAS_SATURATION, &
                                 OIL_MOLE_FRACTION, GAS_MOLE_FRACTION, &
                                 TEMPERATURE]
      allocate(towg_pm%max_change_isubvar(8))
                                           !2 = gas in xmol(gas_c,oil_phase)
      towg_pm%max_change_isubvar = [0,0,0,0,2,1,0]
                                             !1 = gas in xmol(gas_c,gas_phase)
    case('SOLVENT_TL')
      allocate(towg_pm%max_change_ivar(8))
      towg_pm%max_change_ivar = [OIL_PRESSURE, GAS_PRESSURE, &
                                 OIL_SATURATION, GAS_SATURATION, &
                                 OIL_MOLE_FRACTION, GAS_MOLE_FRACTION, &
                                 TEMPERATURE, SOLVENT_SATURATION]
      allocate(towg_pm%max_change_isubvar(8))
                                           !2 = gas in xmol(gas_c,oil_phase)
      towg_pm%max_change_isubvar = [0,0,0,0,2,1,0,0]
                                             !1 = gas in xmol(gas_c,gas_phase)
    case default
      call InputKeywordUnrecognized(miscibility_model, &
                         'TOWG MISCIBILITY_MODEL',option)
  end select

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
  towg_pm%name = 'TOWG Flow'

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
      case('NO_OIL')
        towg_no_oil = PETSC_TRUE
      case('NO_GAS')
        towg_no_gas = PETSC_TRUE
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

recursive subroutine PMTOWGInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 

  use Realization_Base_class
  
  implicit none
  
  class(pm_towg_type) :: this
  
  PetscInt :: num_var,i
  PetscErrorCode :: ierr

  num_var = size(this%max_change_ivar(:))

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,num_var, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, num_var
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMTOWGInitializeRun

! ************************************************************************** !

subroutine PMTOWGInitializeTimestep(this)
  ! 
  ! To be replaced by PreSolve? Which is currently empty...
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/07/16 
  ! 

  use TOWG_module, only : TOWGInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_towg_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)                                 

  !PO: To be removed? And common to all flow modes?
!geh:remove   everywhere                                
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY, &
                                 ZERO_INTEGER)
                                 
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," TOWG FLOW ",64("="))')
  endif
  
  call TOWGInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)                                 
  
end subroutine PMTOWGInitializeTimestep

! ************************************************************************** !

subroutine PMTOWGUpdateAuxVars(this)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 

  use TOWG_module, only : TOWGUpdateAuxVars

  implicit none
  
  class(pm_towg_type) :: this

  call TOWGUpdateAuxVars(this%realization,PETSC_FALSE)

end subroutine PMTOWGUpdateAuxVars   

! ************************************************************************** !

subroutine PMTOWGResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/28/16
  ! 

  use TOWG_module, only : TOWGResidual

  implicit none
  
  class(pm_towg_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call TOWGResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTOWGResidual

! ************************************************************************** !

subroutine PMTOWGJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/28/16
  ! 

  use TOWG_module, only : TOWGJacobian 

  implicit none
  
  class(pm_towg_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call TOWGJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMTOWGJacobian

! ************************************************************************** !

subroutine PMTOWGPreSolve(this)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/28/16

  implicit none

  class(pm_towg_type) :: this

end subroutine PMTOWGPreSolve

! ************************************************************************** !

subroutine PMTOWGCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/22/15
  ! 

  use TOWG_module, only : TOWGCheckUpdatePre

  implicit none
  
  class(pm_towg_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  call TOWGCheckUpdatePre(line_search,X,dX,changed,this%realization, &
                          this%max_it_before_damping,this%damping_factor, &
                          this%trunc_max_pressure_change,ierr)

end subroutine PMTOWGCheckUpdatePre

! ************************************************************************** !


subroutine PMTOWGUpdateSolution(this)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/06/16 
  ! 

  use TOWG_module, only : TOWGUpdateSolution, &
                          TOWGMapBCAuxVarsToGlobal

  implicit none
  
  class(pm_towg_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call TOWGUpdateSolution(this%realization)
  call TOWGMapBCAuxVarsToGlobal(this%realization)

end subroutine PMTOWGUpdateSolution     

! ************************************************************************** !

subroutine PMTOWGTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use TOWG_module, only : TOWGTimeCut

  implicit none
  
  class(pm_towg_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call TOWGTimeCut(this%realization)

end subroutine PMTOWGTimeCut

! ************************************************************************** !

subroutine PMTOWGMaxChange(this)
  ! 
  ! Compute primary variable max changes
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/16
  ! 
  use TOWG_module, only : TOWGMaxChange

  implicit none
  
  class(pm_towg_type) :: this

  call TOWGMaxChange(this%realization, this%max_change_ivar, &
                     this%max_change_isubvar, this%max_pressure_change, &
                     this%max_xmol_change,this%max_saturation_change, &
                     this%max_temperature_change)

end subroutine PMTOWGMaxChange

! ************************************************************************** !

subroutine PMTOWGUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/16
  ! 

  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION

  implicit none
  
  class(pm_towg_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ut, ux, us, umin
  PetscReal :: dtt
  type(field_type), pointer :: field
  
#ifdef PM_TOWG_DEBUG  
  call printMsg(this%option,'PMTOWG%UpdateTimestep()')
#endif
  
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    ut = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
    ux = this%xmol_change_governor/(this%max_xmol_change+1.d-5)
    us = this%saturation_change_governor/(this%max_saturation_change+1.d-5)
    umin = min(up,ut,ux,us)
  endif
  ifac = max(min(num_newton_iterations,size(tfac)),1)
  dtt = fac * dt * (1.d0 + umin)
  dt = min(dtt,tfac(ifac)*dt,dt_max)
  dt = max(dt,dt_min)

  !if (Initialized(this%cfl_governor)) then
  !  ! Since saturations are not stored in global_auxvar for general mode, we
  !  ! must copy them over for the CFL check
  !  ! liquid saturation
  !  field => this%realization%field
  !  call RealizationGetVariable(this%realization,field%work, &
  !                              LIQUID_SATURATION,ZERO_INTEGER)
  !  call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
  !  call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
  !                             LIQUID_SATURATION,TIME_NULL)
  !  call RealizationGetVariable(this%realization,field%work, &
  !                              GAS_SATURATION,ZERO_INTEGER)
  !  call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
  !  call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
  !                             GAS_SATURATION,TIME_NULL)
  !  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  !endif

end subroutine PMTOWGUpdateTimestep

! ************************************************************************** !

subroutine PMTOWGComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/30/16
  ! 

  use TOWG_module, only : TOWGComputeMassBalance

  implicit none
  
  class(pm_towg_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call TOWGComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMTOWGComputeMassBalance

! ************************************************************************** !

subroutine PMTOWGCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with TOWG PM
  ! If saves istate for all TOWG submodels, but required only for black oil
  ! and solvent models
  ! 
  ! Author: Paolo Orsini
  ! Date: 01/19/17

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_towg_type) :: this
  PetscViewer :: viewer
  
  !call GlobalGetAuxVarVecLoc(this%realization, &
  !                           this%realization%field%iphas_loc, &
  !                           STATE,ZERO_INTEGER)
  call PMSubsurfaceFlowCheckpointBinary(this,viewer)
  
end subroutine PMTOWGCheckpointBinary

! ************************************************************************** !

subroutine PMTOWGRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with TOWG PM
  ! If loads istate for all TOWG submodels, but required only for black oil
  ! and solvent models
  !
  ! Author: Paolo Orsini
  ! Date: 01/19/17

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_towg_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowRestartBinary(this,viewer)
  call GlobalSetAuxVarVecLoc(this%realization, &
                             this%realization%field%iphas_loc, &
                             STATE,ZERO_INTEGER)
  
end subroutine PMTOWGRestartBinary

! ************************************************************************** !

subroutine PMTOWGDestroy(this)
  ! 
  ! Destroys TOWG process model
  ! 
  ! Author: Paolo Orsini
  ! Date: 01/19/17
  ! 

  use TOWG_module, only : TOWGDestroy

  implicit none
  
  class(pm_towg_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  ! preserve this ordering
  call TOWGDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMTOWGDestroy

! ************************************************************************** !

end module PM_TOWG_class

