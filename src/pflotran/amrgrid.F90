module AMR_Grid_module

  use Grid_module

  implicit none
 
  private
 
#include "definitions.h"

#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"

  ! the next type encapsulates a pointer to a grid info object
  type, public:: gridptr_type
     type(grid_type), pointer :: grid_ptr
  end type gridptr_type

  ! the next type encapsulates a pointer to a patch info object
  type, public:: gridlevelptr_type
     type(gridptr_type), dimension(:), pointer :: grids
  end type gridlevelptr_type

  type, public :: amrgrid_type

  type(gridlevelptr_type), dimension(:), pointer :: gridlevel ! pointer to the grid levels

  PetscFortranAddr p_samr_hierarchy ! pointer to the SAMRAI hierarchy

  end type amrgrid_type

  public:: AMRGridCreate, &
           AMRGridDestroy, &
           AMRGridCreateLevelPatchLists, &
           AMRGridInitialize, &
           AMRGridCreateVector, &
           AMRGridCreateJacobian, &
           AMRGridComputeGeometryInformation, &
           AMRGridGlobalToLocal, &
           AMRGridLocalToGlobal, &
           AMRGridLocalToLocal

contains

! ************************************************************************** !
!
! AMRGridCreate: Creates a structured AMR grid
! author: Bobby Philip
! date: 03/10/08
!
! ************************************************************************** !
function AMRGridCreate()

  implicit none
  
  type(amrgrid_type), pointer :: AMRGridCreate
  type(amrgrid_type), pointer :: amrgrid

  allocate(amrgrid)
  nullify(amrgrid%gridlevel)
  amrgrid%p_samr_hierarchy = 0

  AMRGridCreate => amrgrid

end function AMRGridCreate

function AMRGridCreateLevelPatchLists(amrgrid)
  use Level_module
  use Patch_module

  implicit none

  type(level_list_type), pointer :: AMRGridCreateLevelPatchLists
  type(level_list_type), pointer :: level_list
  type(level_type), pointer :: level
  type(patch_type), pointer :: patch 
  type(amrgrid_type), pointer :: amrgrid

  interface
     integer function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   integer function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
   end function level_number_patches

   logical function is_local_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
     integer, intent(in) :: pn
   end function is_local_patch

   PetscFortranAddr function hierarchy_get_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
     integer, intent(in) :: pn
   end function hierarchy_get_patch
   
   subroutine samr_physical_dimensions(p_hierarchy, nx, ny, nz)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(inout) :: nx
     integer, intent(inout) :: ny
     integer, intent(inout) :: nz
   end subroutine samr_physical_dimensions
  end interface

#include "include/finclude/petsc.h"
  PetscFortranAddr :: p_samr_hierarchy
  integer :: nlevels
  integer :: npatches
  integer :: ln
  integer :: pn
  logical :: islocal

  p_samr_hierarchy = amrgrid%p_samr_hierarchy

  nlevels =  hierarchy_number_levels(p_samr_hierarchy)

  level_list => LevelCreateList()

  do ln=0,nlevels-1
     level => LevelCreate()
     call LevelAddToList(level,level_list)
     npatches = level_number_patches(p_samr_hierarchy, ln )
     do pn=0,npatches-1
        islocal = is_local_patch(p_samr_hierarchy, ln, pn);
        if(islocal) then
           patch => PatchCreate()
           patch%grid => amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr
        else
           nullify(patch)
        endif
        call PatchAddToList(patch,level%patch_list)
     end do
  end do

  AMRGridCreateLevelPatchLists => level_list
  
end function AMRGridCreateLevelPatchLists

! ************************************************************************** !
!
! AMRGridInitialize: initialize data associated with the AMR grid
! author: Bobby Philip
! date: 03/10/08
!
! ************************************************************************** !
subroutine AMRGridInitialize(amrgrid)
  use Grid_module
  use Structured_Grid_module

  implicit none

  type(amrgrid_type), pointer :: amrgrid

  interface
     integer function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   integer function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
   end function level_number_patches

   logical function is_local_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
     integer, intent(in) :: pn
   end function is_local_patch

   PetscFortranAddr function hierarchy_get_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
     integer, intent(in) :: pn
   end function hierarchy_get_patch
   
   subroutine samr_physical_dimensions(p_hierarchy, nx, ny, nz)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(inout) :: nx
     integer, intent(inout) :: ny
     integer, intent(inout) :: nz
   end subroutine samr_physical_dimensions

   subroutine samr_get_origin(p_hierarchy, x0, y0, z0)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscReal, intent(inout) :: x0
     PetscReal, intent(inout) :: y0
     PetscReal, intent(inout) :: z0
   end subroutine samr_get_origin
  end interface

#include "include/finclude/petsc.h"
  PetscFortranAddr :: p_samr_hierarchy
  PetscInt :: nx, ny, nz
  PetscReal:: x0, y0, z0

  integer :: nlevels
  integer :: npatches
  integer :: ln
  integer :: pn
  logical :: islocal

  type(gridlevelptr_type), dimension(:), pointer :: gridlevel
  type(structured_grid_type), pointer :: struct_grid

  if(associated(amrgrid) .and. amrgrid%p_samr_hierarchy/=0 ) then
     p_samr_hierarchy = amrgrid%p_samr_hierarchy

     nlevels =  hierarchy_number_levels(p_samr_hierarchy)
 
     allocate(amrgrid%gridlevel(nlevels))
     
     gridlevel => amrgrid%gridlevel

     call samr_get_origin(p_samr_hierarchy, x0, y0, z0)

     do ln=0,nlevels-1
        npatches = level_number_patches(p_samr_hierarchy, ln )
        allocate(gridlevel(ln+1)%grids(npatches))
        do pn=0,npatches-1
           islocal = is_local_patch(p_samr_hierarchy, ln, pn);
           if(islocal) then
              gridlevel(ln+1)%grids(pn+1)%grid_ptr => GridCreate()
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%itype = STRUCTURED_GRID 
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%ctype = 'structured'
              struct_grid=>StructuredGridCreate()
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%structured_grid=>struct_grid
              struct_grid%p_samr_patch=hierarchy_get_patch(p_samr_hierarchy, ln, pn)
              call samr_physical_dimensions(p_samr_hierarchy, nx,ny,nz)
              struct_grid%nx = nx
              struct_grid%ny = ny
              struct_grid%nz = nz
              struct_grid%nxy = nx*ny
              struct_grid%nmax = nx*ny*nz
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%nmax = struct_grid%nmax
              struct_grid%origin(X_DIRECTION)=x0
              struct_grid%origin(Y_DIRECTION)=y0
              struct_grid%origin(Z_DIRECTION)=z0
           else              
              nullify(gridlevel(ln+1)%grids(pn+1)%grid_ptr)
           endif
        end do
     end do
     
  endif

end subroutine AMRGridInitialize

! ************************************************************************** !
!
! AMRGridCreateVector: Creates a global PETSc vector
! author: Bobby Philip
! date: 03/10/08
!
! ************************************************************************** !
subroutine AMRGridCreateVector(amrgrid, dof,vector,vector_type)

  implicit none
  
  type(amrgrid_type), pointer:: amrgrid
  DM :: dm_ptr
  Vec :: vector
  PetscInt :: vector_type
  PetscErrorCode :: ierr
  PetscInt :: dof
  PetscTruth:: use_ghost

  select case (vector_type)
    case(GLOBAL)
      use_ghost=PETSC_FALSE
      call create_samrai_vec(amrgrid%p_samr_hierarchy, dof, use_ghost, vector)
    case(LOCAL)
       use_ghost=PETSC_TRUE
      call create_samrai_vec(amrgrid%p_samr_hierarchy, dof, use_ghost, vector)
    case(NATURAL)
       print *, 'ERROR::SAMRAI will not create PETSc Natural Vecs!!'
       stop
  end select
    
end subroutine AMRGridCreateVector

! ************************************************************************** !
!
! AMRGridCreateJacobian: Creates Jacobian matrix associated with the AMR grid
! author: Bobby Philip
! date: 03/10/08
!
! ************************************************************************** !
subroutine AMRGridCreateJacobian(dm_ptr,Jacobian,option)
   use Option_module
 
  implicit none
  
  DM :: dm_ptr
  PetscErrorCode :: ierr
  Mat :: Jacobian
  type(option_type) :: option

end subroutine AMRGridCreateJacobian

subroutine AMRGridComputeGeometryInformation(amrgrid, field, option)
  use Option_module
  use Field_module
  
  implicit none

  interface
     integer function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   integer function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
   end function level_number_patches
  end interface

  type(amrgrid_type), pointer:: amrgrid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
#include "include/finclude/petsc.h"
  PetscFortranAddr :: p_samr_hierarchy

  type(grid_type), pointer :: grid

  integer :: nlevels
  integer :: npatches
  integer :: ln
  integer :: pn

  p_samr_hierarchy = amrgrid%p_samr_hierarchy

  nlevels =  hierarchy_number_levels(p_samr_hierarchy)

  do ln=0,nlevels-1
     npatches = level_number_patches(p_samr_hierarchy, ln )
     do pn=0,npatches-1
        if (associated(amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr)) then
           grid => amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr
           call GridMapIndices(grid)
           call GridComputeSpacing(grid)
           call GridComputeCoordinates(grid,option)
           call GridComputeVolumes(grid,field%volume,option)
           ! set up internal connectivity, distance, etc.
           call GridComputeInternalConnect(grid,option)
        endif
     end do
  end do

end subroutine AMRGridComputeGeometryInformation

! ************************************************************************** !
!
! AMRGridGlobalToLocal: Performs global to local communication with DM
! author: Bobby Philip
! date: 03/10/08
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine AMRGridGlobalToLocal(dm_ptr,global_vec,local_vec)

  implicit none

  DM :: dm_ptr
  Vec :: global_vec
  Vec :: local_vec
  PetscErrorCode :: ierr

end subroutine AMRGridGlobalToLocal

  
! ************************************************************************** !
!
! AMRGridLocalToGlobal: Performs local to global communication with DM
! author: Bobby Philip
! date: 03/10/08
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine AMRGridLocalToGlobal(dm_ptr,local_vec,global_vec)

  implicit none

  DM :: dm_ptr
  Vec :: local_vec
  Vec :: global_vec
  PetscErrorCode :: ierr

end subroutine AMRGridLocalToGlobal

  
! ************************************************************************** !
!
! AMRGridLocalToLocal: Performs local to local communication with DM
! author: Bobby Philip
! date: 03/10/08
!
! ************************************************************************** !
subroutine AMRGridLocalToLocal(dm_ptr,local_vec1,local_vec2)

  implicit none

  DM :: dm_ptr
  Vec :: local_vec1
  Vec :: local_vec2
  PetscErrorCode :: ierr

end subroutine AMRGridLocalToLocal


! ************************************************************************** !
!
! AMRGridDestroy: Deallocates a discretization
! author: Bobby Philip
! date: 03/10/08
!
! ************************************************************************** !
subroutine AMRGridDestroy(amrgrid)

  implicit none

  type(amrgrid_type), pointer :: amrgrid
  
  PetscErrorCode :: ierr

end subroutine AMRGridDestroy

end module AMR_Grid_module
