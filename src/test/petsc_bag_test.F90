program test

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscbag.h"
#include "finclude/petscviewer.h"

  type :: header_type
    PetscInt :: int1
    PetscReal :: real1
    PetscInt :: int2
    PetscReal :: real2
    PetscInt :: int3
  end type header_type

  type, extends(header_type) :: header2_type
    PetscInt :: int21 
    PetscInt :: int22
    PetscReal :: real21
  end type header2_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: header2_type
      implicit none
#include "finclude/petscbag.h"
      PetscBag :: bag
      type(header2_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  PetscBag :: bag
  PetscSizeT, parameter :: bagsize = 48
  class(header2_type), pointer :: header
  class(header2_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscErrorCode :: ierr


  call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)

  print *, size(transfer(dummy_header,dummy_char))

#if 1
  call PetscBagCreate(PETSC_COMM_WORLD,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%int1,-999,"int1","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%int2,-999,"int2","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%int3,-999,"int3","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%real1,-999.0,"real1","", &
                            ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%real2,-999.0,"real2","", &
                            ierr);CHKERRQ(ierr)

  call PetscBagRegisterInt(bag,header%int21,-999,"int21","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%int22,-999,"int22","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%real21,-999.0,"real21","", &
                            ierr);CHKERRQ(ierr)
  header%int1 = 98
  header%int2 = 79
  header%int3 = 3
  header%real1 = 99.d0
  header%real2 = 0.75d0

  header%int21 = 56
  header%int22 = 86
  header%real21 = 69.69d0

  print *, 'header%int1 = ', header%int1
  print *, 'header%int2 = ', header%int2
  print *, 'header%int3 = ', header%int3
  print *, 'header%int21 = ', header%int21
  print *, 'header%int22 = ', header%int22
  print *, 'header%real1 = ', header%real1
  print *, 'header%real2 = ', header%real2
  print *, 'header%real21 = ', header%real21
  print *, ''

  call PetscBagView(bag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
#endif

  call PetscFinalize(ierr);CHKERRQ(ierr)
  
end program test
