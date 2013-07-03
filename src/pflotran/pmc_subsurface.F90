module PMC_Subsurface_class

  use PMC_Base_class
  use Realization_class

  implicit none

#include "definitions.h"
  
  private

  type, public, extends(pmc_base_type) :: pmc_subsurface_type
    class(realization_type), pointer :: realization
  contains
    procedure, public :: Init => PMCSubsurfaceInit
  end type pmc_subsurface_type
  
  public :: PMCSubsurfaceCreate
  
contains

! ************************************************************************** !
!
! PMCSubsurfaceCreate: Allocates and initializes a new process_model_coupler 
!                      object.
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function PMCSubsurfaceCreate()

  implicit none
  
  class(pmc_subsurface_type), pointer :: PMCSubsurfaceCreate
  
  class(pmc_subsurface_type), pointer :: pmc

  print *, 'PMCSubsurface%Create()'
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSubsurfaceCreate => pmc  
  
end function PMCSubsurfaceCreate

! ************************************************************************** !
!
! PMCSubsurfaceInit: Initializes a new process model coupler object.
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
subroutine PMCSubsurfaceInit(this)

  implicit none
  
  class(pmc_subsurface_type) :: this
  
  print *, 'PMCSubsurface%Init()'
  
  call PMCBaseInit(this)
  this%name = 'PMCSubsurface'
  nullify(this%realization)

end subroutine PMCSubsurfaceInit

! ************************************************************************** !
!
! PMCSubsurfaceFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCSubsurfaceFinalizeRun(this)

  use Option_module
  
  implicit none
  
  class(pmc_subsurface_type) :: this
  
  call printMsg(this%option,'PMCSubsurface%FinalizeRun()')
  
  nullify(this%realization)
  
end subroutine PMCSubsurfaceFinalizeRun

! ************************************************************************** !
!
! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
recursive subroutine Destroy(this)

  use Utility_module, only: DeallocateArray
  use Option_module

  implicit none
  
  class(pmc_subsurface_type) :: this
  
  call printMsg(this%option,'PMCSubsurface%Destroy()')
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 
  
end subroutine Destroy
  
end module PMC_Subsurface_class
