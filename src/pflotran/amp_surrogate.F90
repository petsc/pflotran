module Argonne_Mixed_Potential_module
  
  implicit none
  
  private

#include "finclude/petscsys.h"

  public :: AMPInit, &
            AMPRun, &
            AMPDestroy
  
contains

! ************************************************************************** !

subroutine AMPInit()
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/15
  ! 

  implicit none
  
end subroutine AMPInit

! ************************************************************************** !

subroutine AMPRun(time,conc,flux,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/15
  !

  implicit none

  interface
    subroutine AMP_step ( sTme, conc, initialRun, flux, status )
      real ( kind = 8), intent( in )  :: sTme   
      real ( kind = 8), intent( inout ),  dimension (:,:) :: conc
      logical ( kind = 4), intent( in ) :: initialRun
      real ( kind = 8), intent(out), dimension (:,:) :: flux
      integer ( kind = 4), intent(out) :: status
    end subroutine
  end interface  
  
  PetscReal :: time
  PetscReal :: conc(:,:)
  PetscReal :: flux(:,:)
  PetscErrorCode :: ierr
  
  integer ( kind = 4) :: status
  logical ( kind = 4) :: initialRun = PETSC_FALSE
  
  ierr = 0
  call AMP_step(time, conc, PETSC_FALSE, flux, status)
  if (status /= 1) ierr = 1
  
end subroutine AMPRun

! ************************************************************************** !

subroutine AMPDestroy()
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/15
  ! 

  implicit none
  
end subroutine AMPDestroy

end module Argonne_Mixed_Potential_module
