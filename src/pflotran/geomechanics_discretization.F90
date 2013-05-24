#ifdef GEOMECH
module Geomechanics_Discretization_module

  use Geomech_Grid_module
  use Geomech_Grid_Aux_module
  
  implicit none

  private
 
#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmshell.h90"

  type, public :: dm_ptr_type
    DM :: dm  ! PETSc DM
    type(gmdm_type), pointer :: gmdm
  end type dm_ptr_type

  type, public :: geomech_discretization_type
    PetscInt :: itype                          ! type of discretization (e.g. structured, unstructured, etc.)
    character(len=MAXWORDLENGTH) :: ctype      ! name of discretization
    PetscReal :: origin(3)                     ! origin of global domain
    type(geomech_grid_type), pointer :: grid   ! pointer to a grid object
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: dm_index_to_ndof(3)            ! mapping between a dm_ptr to the number of degrees of freedom
    type(dm_ptr_type), pointer :: dm_1dof
    type(dm_ptr_type), pointer :: dm_ngeodof    
  end type geomech_discretization_type

  public :: GeomechDiscretizationCreate, &
            GeomechDiscretizationDestroy  
            
contains

! ************************************************************************** !
!
! GeomechDiscretizationCreate: Creates a geomechanics discretization
! author: Satish Karra, LANL
! date: 05/23/2013
!
! ************************************************************************** !
function GeomechDiscretizationCreate()

  implicit none
  
  type(geomech_discretization_type), pointer :: GeomechDiscretizationCreate
  type(geomech_discretization_type), pointer :: discretization
  
  allocate(discretization)
  discretization%ctype = ''
  discretization%itype = 0
  discretization%origin = 0.d0
  discretization%filename = ''

  ! nullify DM pointers
  allocate(discretization%dm_1dof)
  allocate(discretization%dm_ngeodof)
  discretization%dm_1dof%dm = 0
  discretization%dm_ngeodof%dm = 0
  nullify(discretization%dm_1dof%gmdm)
  nullify(discretization%dm_ngeodof%gmdm)  
  nullify(discretization%grid)
  
  GeomechDiscretizationCreate => discretization

end function GeomechDiscretizationCreate

! ************************************************************************** !
!
! GeomechDiscretizationDestroy: Deallocates a geomechanics discretization
! author: Satish Karra, LANL
! date: 05/23/2013
!
! ************************************************************************** !
subroutine GeomechDiscretizationDestroy(discretization)

  implicit none
  
  type(geomech_discretization_type), pointer :: discretization
  
  PetscErrorCode :: ierr
  PetscInt :: i
    
  if (.not.associated(discretization)) return
      
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      if (discretization%dm_1dof%dm /= 0) &
        call DMDestroy(discretization%dm_1dof%dm,ierr)
      discretization%dm_1dof%dm = 0
      if (discretization%dm_ngeodof%dm /= 0) &
        call DMDestroy(discretization%dm_ngeodof%dm,ierr)
      discretization%dm_ngeodof%dm = 0
    case(UNSTRUCTURED_GRID)
      if (associated(discretization%dm_1dof%gmdm)) &
        call GMDMDestroy(discretization%dm_1dof%gmdm)
      if (associated(discretization%dm_ngeodof%gmdm)) &
        call GMDMDestroy(discretization%dm_ngeodof%gmdm)
  end select
  if (associated(discretization%dm_1dof)) &
    deallocate(discretization%dm_1dof)
  nullify(discretization%dm_1dof)
  if (associated(discretization%dm_ngeodof)) &
    deallocate(discretization%dm_ngeodof)
  nullify(discretization%dm_ngeodof)

end subroutine GeomechDiscretizationDestroy

end module Geomechanics_Discretization_module
#endif