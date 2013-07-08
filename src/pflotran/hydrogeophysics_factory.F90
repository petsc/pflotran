module Hydrogeophysics_Factory_module

  use Hydrogeophysics_Simulation_class
  
  implicit none

  private

#include "definitions.h"

  public :: HydrogeophysicsInitialize

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
  use Simulation_Base_class  
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(hydrogeophysics_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => HydrogeophysicsCreate(option)
  call HydrogeophysicsInitPostPetsc(simulation,option)
  
  simulation_base => simulation

end subroutine HydrogeophysicsInitialize

! ************************************************************************** !
!
! HydrogeophysicsInitializePostPetsc: Sets up hydrogeophysics simulation 
!                                     framework after to PETSc initialization
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsInitPostPetsc(simulation, option)

  use Simulation_module
  use Subsurface_Factory_module
  use Hydrogeophysics_Wrapper_module
  use PMC_Hydrogeophysics_class
  use Option_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  class(hydrogeophysics_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  class(pmc_hydrogeophysics_type), pointer :: hydrogeophysics_coupler
  PetscErrorCode :: ierr
  
  ! Init() is called in SubsurfaceInitializePostPetsc
  call SubsurfaceInitializePostPetsc(simulation, option)
  call VecDuplicate(simulation%realization%field%work,simulation%sigma,ierr)
  
  ! add hydrogeophysics coupler to list
  hydrogeophysics_coupler => PMCHydrogeophysicsCreate()
  hydrogeophysics_coupler%option => simulation%option
  hydrogeophysics_coupler%realization => simulation%realization
  ! Petsc vectors are integer quantities, not pointers
  hydrogeophysics_coupler%sigma = simulation%sigma
  simulation%hydrogeophysics_coupler => hydrogeophysics_coupler
  simulation%process_model_coupler_list%below%below => hydrogeophysics_coupler
  
  call HydrogeophysicsWrapperInit(option)
  
end subroutine HydrogeophysicsInitPostPetsc

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
