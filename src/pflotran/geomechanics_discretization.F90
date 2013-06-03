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

  type, public :: gmdm_ptr_type
    DM :: dm  ! PETSc DM
    type(gmdm_type), pointer :: gmdm
  end type gmdm_ptr_type

  type, public :: geomech_discretization_type
    PetscInt :: itype                          ! type of discretization (e.g. structured, unstructured, etc.)
    character(len=MAXWORDLENGTH) :: ctype      ! name of discretization
    PetscReal :: origin(3)                     ! origin of global domain
    type(geomech_grid_type), pointer :: grid   ! pointer to a grid object
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: dm_index_to_ndof(3)            ! mapping between a dm_ptr to the number of degrees of freedom
    type(gmdm_ptr_type), pointer :: dm_1dof
    type(gmdm_ptr_type), pointer :: dm_ngeodof    
  end type geomech_discretization_type

  public :: GeomechDiscretizationCreate, &
            GeomechDiscretizationDestroy, &
            GeomechDiscretizationCreateDMs
            
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
! GeomechDiscretizationCreateDMs: creates distributed, parallel meshes/grids
! If there are multiple degrees of freedom per grid cell, this will call 
! GeomechDiscretizationCreateDM() multiple times to create the DMs corresponding 
! to one degree of freedom grid cell and those corresponding to multiple 
! degrees of freedom per cell for geomechanics.
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationCreateDMs(discretization,option)
      
  use Option_module    
      
  implicit none
  
  type(geomech_discretization_type)            :: discretization
  type(option_type)                            :: option
      
  PetscInt                                     :: ndof
  PetscErrorCode                               :: ierr
  type(geomech_grid_type), pointer             :: geomech_grid

  !-----------------------------------------------------------------------
  ! Generate the DM objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call GeomechDiscretizationCreateDM(discretization,discretization%dm_1dof, &
                                     ndof,option)
  
  if (option%ngeomechdof > 0) then
    ndof = option%ngeomechdof
    call GeomechDiscretizationCreateDM(discretization,discretization%dm_ngeodof, &
                                       ndof,option)
  endif


end subroutine GeomechDiscretizationCreateDMs

! ************************************************************************** !
!
! GeomechDiscretizationCreateDM: creates a distributed, parallel mesh/grid
! for geomechanics
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationCreateDM(discretization,dm_ptr,ndof,option)

  use Option_module
  
  implicit none
  
  type(geomech_discretization_type) :: discretization
  type(gmdm_ptr_type), pointer                         :: dm_ptr
  PetscInt                                             :: ndof
  type(option_type)                                    :: option
  PetscErrorCode                                       :: ierr

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      option%io_buffer = &
        'Geomechanics currently works only with unstructured grid.'
      call printErrMsg(option)
    case(UNSTRUCTURED_GRID)
      call GMCreateGMDM(discretization%grid, &
                        dm_ptr%gmdm,ndof,option)
      call DMShellCreate(option%mycomm,dm_ptr%dm,ierr)
 !     call DMShellSetGlobalToLocalVecScatter(dm_ptr%dm,dm_ptr%gmdm%scatter_gtol,ierr)
  end select

end subroutine GeomechDiscretizationCreateDM


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