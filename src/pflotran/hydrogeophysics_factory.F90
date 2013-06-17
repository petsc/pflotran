module Hydrogeophysics_Factory_module

  use Hydrogeophysics_Simulation_class
!  use PMC_Hydrogeophysics_class
  
  implicit none

  private

#include "definitions.h"

  public :: HydrogeophysicsInitialize, &
            HydrogeophysicsRun, &
            HydrogeophysicsFinalize

contains

! ************************************************************************** !
!
! HydrogeophysicsInitialize: Sets up hydrogeophysics simulation 
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsInitialize(simulation_base,option)

  use Option_module
  use Input_module
  use Timestepper_module
  use Simulation_Base_class  
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(hydrogeophysics_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => HydrogeophysicsCreate(option)
  call HydrogeophysicsInitPostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine HydrogeophysicsInitialize

! ************************************************************************** !
!
! HydrogeophysicsInitializePostPETSc: Sets up hydrogeophysics simulation 
!                                     framework after to PETSc initialization
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsInitPostPETSc(simulation, option)

  use Simulation_module
  use Subsurface_Factory_module
  use Option_module
  use Init_module
  
  implicit none
  
  class(hydrogeophysics_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  call SubsurfaceInitializePostPETSc(simulation, option)
  
end subroutine HydrogeophysicsInitPostPETSc

! ************************************************************************** !
!
! HydrogeophysicsRun: Runs the hydrogeophysics simulation
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsRun(simulation, master_stepper, init_status)

  use Simulation_module
  use Timestepper_module
  
  implicit none
  
  type(hydrogeophysics_simulation_type) :: simulation
  type(stepper_type), pointer :: master_stepper
  PetscInt :: init_status

end subroutine HydrogeophysicsRun

! ************************************************************************** !
!
! HydrogeophysicsFinalize: Destroys hydrogeophysics simulation framework
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsFinalize(simulation,option)

  use Option_module
  use Simulation_Base_class
  
  implicit none
  
  class(simulation_base_type) :: simulation
  type(option_type) :: option
  
  call simulation%FinalizeRun()
!  call HydrogeophysicsSimulationDestroy(simulation)
  call simulation%Strip()

end subroutine HydrogeophysicsFinalize

! ************************************************************************** !
!
! HydrogeoInitCommandLineSettings: Initializes hydrogeophysics settings
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeoInitCommandLineSettings(option)

  use Option_module
  use Input_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
!  string = '-dummy'
!  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
!  if (option_found) then
!    option%subsurface_simulation_type = dummy
!  endif
  
end subroutine HydrogeoInitCommandLineSettings

end module Hydrogeophysics_Factory_module
