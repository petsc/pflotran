module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

  type, public :: field_type 
  

    
!geh material id
    PetscInt, pointer :: imat(:)
    
    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)

    ! 1 degree of freedom
    Vec :: porosity0, porosity_loc
    Vec :: tor_loc
    Vec :: ithrm_loc
    Vec :: icap_loc
    Vec :: iphas_loc, iphas_old_loc

    Vec :: perm_xx_loc, perm_yy_loc, perm_zz_loc
    Vec :: perm0_xx, perm0_yy, perm0_zz, perm_pow
    
    ! NDOF degree of freedom
    ! residual vector
    Vec :: r            
    ! Solution vectors
    Vec :: xx, xx_loc, dxx, yy, accum
   
  end type 

  public :: FieldCreate, &
            FieldDestroy

contains

! ************************************************************************** !
!
! FieldCreate: Allocates and initializes a new Field object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function FieldCreate()

  implicit none
  
  type(field_type), pointer :: FieldCreate
  
  type(field_type), pointer :: field
  
  allocate(field)
  
  ! nullify PetscVecs
  field%porosity0 = 0
  field%porosity_loc = 0
  field%tor_loc = 0
  field%ithrm_loc = 0
  field%icap_loc = 0
  field%iphas_loc = 0
  field%iphas_old_loc = 0

  field%perm_xx_loc = 0
  field%perm_yy_loc = 0
  field%perm_zz_loc = 0
  field%perm0_xx = 0
  field%perm0_yy = 0
  field%perm0_zz = 0
  field%perm_pow = 0
  
  field%r = 0
  field%xx = 0
  field%xx_loc = 0
  field%dxx = 0
  field%yy = 0
  field%accum = 0
  
  nullify(field%imat)
  nullify(field%internal_velocities)
  nullify(field%boundary_velocities)
  
  FieldCreate => field
  
end function FieldCreate

! ************************************************************************** !
!
! FieldDestroy: Deallocates a field object
! author: Glenn Hammond
! date: 11/15/07
!
! ************************************************************************** !
subroutine FieldDestroy(field)

  

  implicit none
  
  type(field_type), pointer :: field
  
  PetscMPIInt :: myrank, ierr
  
  call MPI_Comm_Rank(PETSC_COMM_WORLD,myrank,ierr)
  if (myrank == 0) then
    print *, 'Need to implement FieldDestroy'
  endif
  
  ! all kinds of stuff needs to be added here.
  
  
  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
