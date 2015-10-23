module PM_TOilIms_class

  use PM_Base_class
  use PM_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(pm_subsurface_type) :: pm_toil_ims_type
    PetscReal :: dPmax
    PetscReal :: dTmax
    PetscReal :: dSmax
    PetscReal :: dPmax_allowable
    PetscReal :: dTmax_allowable
    PetscReal :: dSmax_allowable
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
  contains
    ! all the routines below needs to be replaced, uncomment as I develop them
    procedure, public :: Read => PMTOilImsRead
    !procedure, public :: SetupSolvers => PMGeneralSetupSolvers
    procedure, public :: InitializeRun => PMTOilImsInitializeRun
    procedure, public :: InitializeTimestep => PMTOilImsInitializeTimestep
    !procedure, public :: Residual => PMGeneralResidual
    !procedure, public :: Jacobian => PMGeneralJacobian
    !procedure, public :: UpdateTimestep => PMGeneralUpdateTimestep
    procedure, public :: PreSolve => PMTOilImsPreSolve
    !procedure, public :: PostSolve => PMGeneralPostSolve
    procedure, public :: CheckUpdatePre => PMTOilImsCheckUpdatePre
    !procedure, public :: CheckUpdatePost => PMGeneralCheckUpdatePost
    !procedure, public :: TimeCut => PMGeneralTimeCut
    procedure, public :: UpdateSolution => PMTOilImsUpdateSolution
    procedure, public :: UpdateAuxVars => PMTOilImsUpdateAuxVars
    !procedure, public :: MaxChange => PMGeneralMaxChange
    !procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    !procedure, public :: CheckpointBinary => PMGeneralCheckpointBinary
    !procedure, public :: RestartBinary => PMGeneralRestartBinary
    !procedure, public :: Destroy => PMGeneralDestroy
  end type pm_toil_ims_type
  
  public :: PMToilImsCreate
  
contains

! ************************************************************************** !

function PMTOilImsCreate()
  ! 
  ! Creates TOilIms process models shell
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 9/8/15
  ! 
  use Variables_module, only : LIQUID_PRESSURE, OIL_PRESSURE, OIL_SATURATION, &
                               TEMPERATURE
  implicit none
  
  class(pm_toil_ims_type), pointer :: PMToilImsCreate

  class(pm_toil_ims_type), pointer :: toil_ims_pm
  
!#ifdef PM_TOIL_IMS_DEBUG  
  print *, 'PMTOilImsCreate()'
!#endif  

  allocate(toil_ims_pm)

  toil_ims_pm%dPmax = 0.d0
  toil_ims_pm%dTmax = 0.d0
  toil_ims_pm%dSmax = 0.d0
  toil_ims_pm%dPmax_allowable = 5.d5 !Pa
  toil_ims_pm%dTmax_allowable = 5.d0
  toil_ims_pm%dSmax_allowable = 1.d0
  allocate(toil_ims_pm%max_change_ivar(4))
  toil_ims_pm%max_change_ivar = [LIQUID_PRESSURE, OIL_PRESSURE, &
                                OIL_SATURATION, TEMPERATURE]
  allocate(toil_ims_pm%max_change_isubvar(4))
  toil_ims_pm%max_change_isubvar = [0,0,0,0]
  
  call PMSubsurfaceCreate(toil_ims_pm)
  toil_ims_pm%name = 'PMTOilIms'

  PMTOilImsCreate => toil_ims_pm
  
end function PMTOilImsCreate

! ************************************************************************** !

subroutine PMTOilImsRead(this,input)
  ! 
  ! Reads input specific to pm_toil_Ims.
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: Date: 9/9/15
  !
  ! use TOilIms_module ! shouldn't need this... 
  use TOilIms_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_toil_ims_type) :: this  
  type(input_type) :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  
  type(option_type), pointer :: option

  option => this%option

  call InputReadWord(input,option,keyword,PETSC_TRUE)
  if (input%ierr /= 0) then
    return
  endif
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','TOIL_IMS_MODE')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('TOUGH2_ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,toil_ims_tgh2_itol_scld_res_e1)
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e1')
        call InputReadDouble(input,option,toil_ims_tgh2_itol_scld_res_e2)
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e2')
        toil_ims_tough2_conv_criteria = PETSC_TRUE
      case('WINDOW_EPSILON') 
        call InputReadDouble(input,option,toil_ims_window_epsilon)
        call InputErrorMsg(input,option,'window epsilon','TOIL_IMS_MODE')
      ! consider to move in this in eos_oil, since this is an eos property
      !case('OIL_COMPONENT_FORMULA_WEIGHT')
      !  !assuming oil component is index 2, H2O ois index 1
      !  call InputReadDouble(input,option,toil_ims_fmw_comp(2))
      !  call InputErrorMsg(input,option,'oil component formula wt.', &
      !                     'TOIL_IMS_MODE')
      case('ISOTHERMAL')
        toil_ims_isothermal = PETSC_TRUE
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,toil_ims_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           'TOIL_IMS_MODE')
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,toil_ims_max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           'TOIL_IMS_MODE')
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,toil_ims_damping_factor)
        call InputErrorMsg(input,option,'damping factor','TOIL_IMS_MODE')
      case('GOVERN_MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,this%dPmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable pressure change', &
                           'TOIL_IMS_MODE')
      case('GOVERN_MAXIMUM_TEMPERATURE_CHANGE')
        call InputReadDouble(input,option,this%dTmax_allowable)
        call InputErrorMsg(input,option, &
                           'maximum allowable temperature change', &
                           'TOIL_IMS_MODE')
      case('GOVERN_MAXIMUM_SATURATION_CHANGE')
        call InputReadDouble(input,option,this%dSmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable saturation change', &
                           'TOIL_IMS_MODE')
      case('DEBUG_CELL')
        call InputReadInt(input,option,toil_ims_debug_cell_id)
        call InputErrorMsg(input,option,'debug cell id','TOIL_IMS_MODE')
      ! might need some input here for the thermal diffusion model
      !case('NO_TEMP_DEPENDENT_DIFFUSION')
      !  general_temp_dep_gas_air_diff = PETSC_FALSE
      !case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
      !  general_harmonic_diff_density = PETSC_TRUE
      !case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
      !  general_harmonic_diff_density = PETSC_FALSE
      case default
        call InputKeywordUnrecognized(keyword,'TOIL_IMS Mode',option)
    end select
    
  enddo  
  
end subroutine PMTOilImsRead

! ************************************************************************** !

recursive subroutine PMTOilImsInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15

  use Realization_Base_class
  
  implicit none
  
  class(pm_toil_ims_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,FOUR_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 4
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo


  ! call parent implementation
  call PMSubsurfaceInitializeRun(this)

end subroutine PMTOilImsInitializeRun

! ************************************************************************** !

subroutine PMTOilImsInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use TOilIms_module, only : TOilImsInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_toil_ims_type) :: this

  call PMSubsurfaceInitializeTimestepA(this)                                 
!geh:remove   everywhere                                
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY,0)
                                 
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," TOIL_IMS FLOW ",64("="))')
  endif
  
  call TOilImsInitializeTimestep(this%realization)

  call PMSubsurfaceInitializeTimestepB(this)                                 
  
end subroutine PMTOilImsInitializeTimestep

! ************************************************************************** !

subroutine PMTOilImsPreSolve(this)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15

  implicit none

  class(pm_toil_ims_type) :: this

  ! currently does nothing - could add here explicit iitialization
  ! for highly het. problems

end subroutine PMTOilImsPreSolve

! ************************************************************************** !

subroutine PMTOilImsUpdateAuxVars(this)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15

  use TOilIms_module, only : TOilImsUpdateAuxVars

  implicit none
  
  class(pm_toil_ims_type) :: this

  call TOilImsUpdateAuxVars(this%realization)

end subroutine PMTOilImsUpdateAuxVars   

! ************************************************************************** !

subroutine PMTOilImsUpdateSolution(this)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15
  ! 

  use TOilIms_module, only : TOilImsUpdateSolution, &
                             TOilImsMapBCAuxVarsToGlobal 

  implicit none
  
  class(pm_toil_ims_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call TOilImsUpdateSolution(this%realization)
  call TOilImsMapBCAuxVarsToGlobal(this%realization)

end subroutine PMTOilImsUpdateSolution     

! ************************************************************************** !

subroutine PMTOilImsCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/22/15
  ! 

  use TOilIms_module, only : TOilImsCheckUpdatePre

  implicit none
  
  class(pm_toil_ims_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call TOilImsCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

end subroutine PMTOilImsCheckUpdatePre

! ************************************************************************** !

end module PM_TOilIms_class

