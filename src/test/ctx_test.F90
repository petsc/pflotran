module Base_module
  implicit none
  private
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"
  type, public :: base_type
    PetscInt :: A
    PetscReal :: I
  contains
    procedure, public :: Print => BasePrint
  end type base_type
  public :: BaseFunction
contains
subroutine BasePrint(this)
  implicit none
  class(base_type) :: this
  print *, 'Base printout'
end subroutine BasePrint
subroutine BaseFunction(snes,xx,r,ctx,ierr)
  implicit none
  SNES :: snes
  Vec :: xx
  Vec :: r
  class(base_type) :: ctx
  PetscErrorCode :: ierr
  call ctx%Print()
end subroutine BaseFunction
end module Base_module

module Extended_module
  use Base_module
  implicit none
  private
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"
  type, public, extends(base_type) :: extended_type
    PetscInt :: B
    PetscReal :: J
  contains
    procedure, public :: Print =>  ExtendedPrint
  end type extended_type
  public :: ExtendedFunction
contains
subroutine ExtendedPrint(this)
  implicit none
  class(extended_type) :: this
  print *, 'Extended printout'
end subroutine ExtendedPrint
subroutine ExtendedFunction(snes,xx,r,ctx,ierr)
  implicit none
  SNES :: snes
  Vec :: xx
  Vec :: r
!  class(extended_type) :: ctx
  class(base_type) :: ctx
  PetscErrorCode :: ierr
  call ctx%Print()
end subroutine ExtendedFunction
end module Extended_module

program test

  use Base_module
  use Extended_module

  implicit none
  
#include "finclude/petscsys.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  PetscInt :: ndof
  PetscInt :: my_ndof
  
  SNES :: snes_base, snes_extended
  Vec :: x
  class(base_type), pointer :: base
  class(extended_type), pointer :: extended
  PetscErrorCode :: ierr

  print *, 'Start of Fortran90 test program'

  nullify(base)
  nullify(extended)
!  allocate(base)
  allocate(extended)
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  call VecCreate(PETSC_COMM_WORLD,x,ierr);CHKERRQ(ierr)
  call SNESCreate(PETSC_COMM_WORLD,snes_base,ierr);CHKERRQ(ierr)
 
  ! when I use the base class as the context
  print *, 'the base class will succeed by printing out "Base printout"'
  call SNESCreate(PETSC_COMM_WORLD,snes_base,ierr);CHKERRQ(ierr)
  call SNESSetFunction(snes_base,x,BaseFunction,base);CHKERRQ(ierr)
  call SNESComputeFunction(snes_base,x,x,ierr);CHKERRQ(ierr)
  call SNESDestroy(snes_base,ierr);CHKERRQ(ierr)

  ! when I use the extended class as the context
  print *, 'the extended class will succeed by printing out "Extended printout"'
  call SNESCreate(PETSC_COMM_WORLD,snes_extended,ierr);CHKERRQ(ierr)
  call SNESSetFunction(snes_extended,x,ExtendedFunction,extended);CHKERRQ(ierr)
  call SNESComputeFunction(snes_extended,x,x,ierr);CHKERRQ(ierr)
  call VecDestroy(x,ierr);CHKERRQ(ierr)
  call SNESDestroy(snes_extended,ierr);CHKERRQ(ierr)
  if (associated(base)) deallocate(base)
  if (associated(extended)) deallocate(extended)
  call PetscFinalize(ierr)

  print *, 'End of Fortran90 test program'
 
end program test

