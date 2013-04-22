module Hydrogeophysics_Coupler_module
 
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
 
  type, public :: hydrogeophysic_coupler_type
    Vec :: sigma
  end type hydrogeophysic_coupler_type
  
  public :: HydrogeophysicsCouplerCreate, &
            HydrogeophysicsCouplerDestroy

contains

! ************************************************************************** !
!
! HydrogeophysicsCouplerCreate: Creates a hydrogeophysics object
! author: Glenn Hammond
! date: 04/22/13
!
! ************************************************************************** !
function HydrogeophysicsCouplerCreate()
  
  implicit none

  type(hydrogeophysic_coupler_type), pointer :: HydrogeophysicsCouplerCreate
  
  type(hydrogeophysic_coupler_type), pointer :: coupler
  
  allocate(coupler)
  coupler%sigma = 0

  HydrogeophysicsCouplerCreate => coupler

end function HydrogeophysicsCouplerCreate

! ************************************************************************** !
!
! HydrogeophysCouplerSetMapping: Sets up mapping between PFLOTRAN and E4D
! author: Glenn Hammond
! date: 04/22/13
!
! ************************************************************************** !
subroutine HydrogeophysCouplerSetMapping(coupler)

  implicit none
  
  type(hydrogeophysic_coupler_type), pointer :: coupler  

end subroutine HydrogeophysCouplerSetMapping

! ************************************************************************** !
!
! HydrogeophysicsCouplerMap: Maps data from PFLOTRAN domain to E4D domain
! author: Glenn Hammond
! date: 04/22/13
!
! ************************************************************************** !
subroutine HydrogeophysicsCouplerMap(coupler)

  implicit none
  
  type(hydrogeophysic_coupler_type), pointer :: coupler  

end subroutine HydrogeophysicsCouplerMap

! ************************************************************************** !
!
! HydrogeophysicsCouplerDestroy: hydrogeophysics object
! author: Glenn Hammond
! date: 04/22/13
!
! ************************************************************************** !
subroutine HydrogeophysicsCouplerDestroy(coupler)

  implicit none
  
  type(hydrogeophysic_coupler_type), pointer :: coupler
  
  PetscErrorCode :: ierr
  
  if (.not.associated(coupler)) return
  
  if (coupler%sigma /= 0) &
    call VecDestroy(coupler%sigma,ierr)
  
  deallocate(coupler)
  nullify(coupler)
  
end subroutine HydrogeophysicsCouplerDestroy

end module Hydrogeophysics_Coupler_module
