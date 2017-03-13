module PM_UFD_Biosphere_class

  use PM_Base_class
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(pm_base_type) :: pm_ufd_biosphere_type
    class(realization_subsurface_type), pointer :: realization
  contains
    procedure, public :: PMUFDBSetRealization
    procedure, public :: Setup => PMUFDBSetup
    procedure, public :: Read => PMUFDBRead
    procedure, public :: InitializeRun => PMUFDBInitializeRun
    procedure, public :: InitializeTimestep => PMUFDBInitializeTimestep
    procedure, public :: FinalizeTimestep => PMUFDBFinalizeTimestep
    procedure, public :: Solve => PMUFDBSolve
    procedure, public :: Output => PMUFDBOutput
    procedure, public :: InputRecord => PMUFDBInputRecord
    procedure, public :: Destroy => PMUFDBDestroy
  end type pm_ufd_biosphere_type

  public :: PMUFDBCreate
  
contains

! *************************************************************************** !

function PMUFDBCreate()
  !
  ! Creates and initializes the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none
  
  class(pm_ufd_biosphere_type), pointer :: PMUFDBCreate
  
  allocate(PMUFDBCreate)
  nullify(PMUFDBCreate%realization)
  PMUFDBCreate%name = 'ufd biosphere'

  call PMBaseInit(PMUFDBCreate)
  
end function PMUFDBCreate

! *************************************************************************** !

subroutine PMUFDBSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  use Realization_Subsurface_class

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMUFDBSetRealization

! *************************************************************************** !

subroutine PMUFDBRead(this,input)
  !
  ! Reads input file parameters for the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(input_type), pointer :: input
  
end subroutine PMUFDBRead

! *************************************************************************** !

subroutine PMUFDBSetup(this)
  !
  ! 
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBSetup

! ************************************************************************** !

subroutine PMUFDBInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBInitializeRun

! *************************************************************************** !

subroutine PMUFDBInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
  PetscReal :: dt

  dt = this%option%flow_dt
  
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," UFD BIOSPHERE MODEL ",45("="))')
  endif

  
end subroutine PMUFDBInitializeTimestep

! *************************************************************************** !

 subroutine PMUFDBSolve(this,time,ierr)
  ! 
  ! 
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

end subroutine PMUFDBSolve

! ************************************************************************** !

subroutine PMUFDBFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBFinalizeTimestep

! *************************************************************************** !

 subroutine PMUFDBOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBOutput

! *************************************************************************** !

subroutine PMUFDBInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  ! 
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMUFDBInputRecord

! *************************************************************************** !

subroutine PMUFDBDestroy(this)
  ! 
  ! Strips and destroys the UFD Biosphere process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
  call PMUFDBStrip(this)
  
end subroutine PMUFDBDestroy

! ************************************************************************** !

subroutine PMUFDBStrip(this)
  ! 
  ! Strips the UFD Biosphere process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this


end subroutine PMUFDBStrip

! ************************************************************************** !

end module PM_UFD_Biosphere_class