program test

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscbag.h"
#include "finclude/petscviewer.h"

  type :: header_type
    sequence
    integer*8 :: value1
    integer*8 :: value2
  end type header_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: header_type
      implicit none
#include "finclude/petscbag.h"
      PetscBag :: bag
      type(header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  PetscBag :: bag
  PetscSizeT, parameter :: bagsize = 16
  type(header_type), pointer :: header
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
  call PetscBagCreate(PETSC_COMM_WORLD,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%value1,-999,"value1","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%value2,-999,"value2","", &
                           ierr);CHKERRQ(ierr)
  header%value1 = 98
  header%value2 = 79

  print *, 'header%value1 = ', header%value1
  print *, 'header%value2 = ', header%value2
  print *, ''

  call PetscBagView(bag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  call PetscFinalize(ierr);CHKERRQ(ierr)
  
end program test
