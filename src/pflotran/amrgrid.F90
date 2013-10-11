module AMR_Grid_module

  use Grid_module
  use Structured_Grid_module

  use PFLOTRAN_Constants_module

  implicit none
 
  private
 
#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"

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

  PetscFortranAddr p_application ! pointer to the SAMRAI hierarchy

  PetscInt :: nlevels
  
  end type amrgrid_type

  public:: AMRGridCreate, &
           AMRGridDestroy, &
           AMRGridCreateLevelPatchLists, &
           AMRGridInitialize, &
           AMRGridComputeLocalBounds, &
           AMRGridComputeGridSpacing, &
           AMRGridCreateVector, &
           AMRGridCreateJacobian, &
           AMRGridComputeGeometryInformation, &
           AMRGridReadDXYZ, &
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
  amrgrid%p_application = 0

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
   PetscInt function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   PetscInt function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
   end function level_number_patches

   PetscInt function is_local_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
     PetscInt, intent(in) :: pn
   end function is_local_patch

   PetscFortranAddr function hierarchy_get_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
     PetscInt, intent(in) :: pn
   end function hierarchy_get_patch
   
subroutine samr_physical_dimensions(p_hierarchy, nx, ny, nz)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(inout) :: nx
     PetscInt, intent(inout) :: ny
     PetscInt, intent(inout) :: nz
end subroutine samr_physical_dimensions
end interface

#include "finclude/petscsysdef.h"
  PetscFortranAddr :: p_application
  PetscInt :: nlevels
  PetscInt :: npatches
  PetscInt :: ln
  PetscInt :: pn
  PetscInt :: islocal

  p_application = amrgrid%p_application

  nlevels =  hierarchy_number_levels(p_application)

  level_list => LevelCreateList()

  do ln=0,nlevels-1
     level => LevelCreate()
     call LevelAddToList(level,level_list)
     level%patch_list => PatchCreateList()
     npatches = level_number_patches(p_application, ln )
     do pn=0,npatches-1
        islocal = is_local_patch(p_application, ln, pn);
        if(islocal == 1) then
           patch => PatchCreate()
           patch%grid => amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr
           call PatchAddToList(patch,level%patch_list)
! this has to be called after the call to PatchAddToList since that routine
! sets the patch id according to the size of the list
           patch%id = pn
        endif
     end do
  end do

  AMRGridCreateLevelPatchLists => level_list
  
end function AMRGridCreateLevelPatchLists

subroutine AMRGridComputeLocalBounds(amrgrid)
  use Grid_module
  use Structured_Grid_module

  implicit none

  type(amrgrid_type), pointer :: amrgrid
  type(structured_grid_type), pointer :: structured_grid

  interface
   PetscInt function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   PetscInt function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
   end function level_number_patches

   PetscInt function is_local_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
     PetscInt, intent(in) :: pn
   end function is_local_patch
  end interface

#include "finclude/petscsysdef.h"
  PetscFortranAddr :: p_application

  PetscInt :: nlevels
  PetscInt :: npatches
  PetscInt :: ln
  PetscInt :: pn
  PetscInt :: islocal
  DM :: da

  p_application = amrgrid%p_application

  nlevels =  hierarchy_number_levels(p_application)

  do ln=0,nlevels-1
     npatches = level_number_patches(p_application, ln )
     do pn=0,npatches-1
        islocal = is_local_patch(p_application, ln, pn);
        if(islocal == 1) then
           structured_grid=>amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr%structured_grid
           call StructGridComputeLocalBounds(structured_grid, da)
           amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr%nlmax = structured_grid%nlmax
           amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr%ngmax = structured_grid%ngmax
        endif
     end do
  end do

end subroutine AMRGridComputeLocalBounds

subroutine AMRGridComputeGridSpacing(amrgrid)
  implicit none
  type(amrgrid_type), pointer :: amrgrid
  type(structured_grid_type), pointer :: structured_grid

  interface   
     PetscInt function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   PetscInt function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
   end function level_number_patches

   PetscInt function is_local_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
     PetscInt, intent(in) :: pn
   end function is_local_patch
  end interface

#include "finclude/petscsysdef.h"
  PetscFortranAddr :: p_application
  PetscFortranAddr :: p_samr_patch
  PetscReal :: dx, dy, dz

  PetscInt :: nlevels
  PetscInt :: npatches
  PetscInt :: ln
  PetscInt :: pn
  PetscInt :: islocal

  
  p_application = amrgrid%p_application

  nlevels =  hierarchy_number_levels(p_application)

  do ln=0,nlevels-1
     npatches = level_number_patches(p_application, ln )
     do pn=0,npatches-1
        islocal = is_local_patch(p_application, ln, pn);
        if(islocal == 1) then
           structured_grid =>amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr%structured_grid
           p_samr_patch = structured_grid%p_samr_patch
           call samr_patch_get_spacing(p_samr_patch, dx, dy, dz)
! for now we will assume that the SAMR grids will use constant grid spacing in each direction
           allocate(structured_grid%dxg_local(structured_grid%ngx))
           structured_grid%dxg_local = dx
           allocate(structured_grid%dyg_local(structured_grid%ngy))
           structured_grid%dyg_local = dy
           allocate(structured_grid%dzg_local(structured_grid%ngz))
           structured_grid%dzg_local = dz
           allocate(structured_grid%dx(structured_grid%ngmax))
           structured_grid%dx = dx
           allocate(structured_grid%dy(structured_grid%ngmax))
           structured_grid%dy = dy
           allocate(structured_grid%dz(structured_grid%ngmax))
           structured_grid%dz = dz
        endif
     end do
  end do
   
end subroutine AMRGridComputeGridSpacing


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
   PetscInt function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   PetscInt function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
   end function level_number_patches

   PetscInt function is_local_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
     PetscInt, intent(in) :: pn
   end function is_local_patch

   PetscFortranAddr function hierarchy_get_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
     PetscInt, intent(in) :: pn
   end function hierarchy_get_patch
   
subroutine samr_physical_dimensions(p_hierarchy, nx, ny, nz)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(inout) :: nx
     PetscInt, intent(inout) :: ny
     PetscInt, intent(inout) :: nz
end subroutine samr_physical_dimensions

subroutine samr_get_origin(p_hierarchy, x0, y0, z0)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscReal, intent(inout) :: x0
     PetscReal, intent(inout) :: y0
     PetscReal, intent(inout) :: z0
end subroutine samr_get_origin

subroutine samr_get_upper_corner(p_hierarchy, x1, y1, z1)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscReal, intent(inout) :: x1
     PetscReal, intent(inout) :: y1
     PetscReal, intent(inout) :: z1
end subroutine samr_get_upper_corner

end interface

#include "finclude/petscsysdef.h"
  PetscFortranAddr :: p_application
  PetscInt :: nx, ny, nz
  PetscReal:: x0, y0, z0
  PetscReal:: x1, y1, z1

  PetscInt :: nlevels
  PetscInt :: npatches
  PetscInt :: ln
  PetscInt :: pn
  PetscInt :: islocal

  type(gridlevelptr_type), dimension(:), pointer :: gridlevel
  type(structured_grid_type), pointer :: struct_grid

  if(associated(amrgrid) .and. amrgrid%p_application/=0 ) then
     p_application = amrgrid%p_application

     nlevels =  hierarchy_number_levels(p_application)
 
     allocate(amrgrid%gridlevel(nlevels))
     
     gridlevel => amrgrid%gridlevel

     call samr_get_origin(p_application, x0, y0, z0)
     call samr_get_upper_corner(p_application, x1, y1, z1)

     do ln=0,nlevels-1
        npatches = level_number_patches(p_application, ln )
        allocate(gridlevel(ln+1)%grids(npatches))
        do pn=0,npatches-1
           islocal = is_local_patch(p_application, ln, pn);
           if(islocal == 1) then
              gridlevel(ln+1)%grids(pn+1)%grid_ptr => GridCreate()
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%itype = STRUCTURED_GRID 
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%ctype = 'structured'
              struct_grid=>StructGridCreate()
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%structured_grid=>struct_grid
              struct_grid%p_samr_patch=hierarchy_get_patch(p_application, ln, pn)
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%x_min_global = x0
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%y_min_global = y0
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%z_min_global = z0
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%x_max_global = x1
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%y_max_global = y1
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%z_max_global = z1
              call samr_physical_dimensions(p_application, nx,ny,nz)
              struct_grid%nx = nx
              struct_grid%ny = ny
              struct_grid%nz = nz
              struct_grid%nxy = nx*ny
              struct_grid%nmax = nx*ny*nz
              gridlevel(ln+1)%grids(pn+1)%grid_ptr%nmax = struct_grid%nmax
              struct_grid%origin(X_DIRECTION)=x0
              struct_grid%origin(Y_DIRECTION)=y0
              struct_grid%origin(Z_DIRECTION)=z0
! for now all AMR grids are going to be cartesian though we can generalize later
              struct_grid%itype = CARTESIAN_GRID
              struct_grid%ctype = 'cartesian'
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
subroutine AMRGridCreateVector(amrgrid, dof, vector,vector_type, &
                               use_components, option)

  use Option_module

  implicit none

interface
subroutine create_samrai_vec(p_application, dof, centering, use_ghost, use_components, vec)
       use Option_module
       implicit none
       
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       
       PetscFortranAddr :: p_application
       PetscInt :: dof
       PetscInt :: centering
       PetscBool :: use_ghost
       PetscBool :: use_components
       Vec :: vec
       type(option_type) :: option
     end subroutine create_samrai_vec
  end interface

  type(amrgrid_type), pointer:: amrgrid
  DM :: dm_ptr
  Vec :: vector
  PetscInt :: vector_type
  PetscBool :: use_components
  type(option_type) :: option
  PetscErrorCode :: ierr
  PetscInt :: dof
  PetscBool:: use_ghost
  PetscInt :: var_centering

  var_centering = option%ivar_centering

  select case (vector_type)
    case(GLOBAL)
      use_ghost=PETSC_FALSE
      call create_samrai_vec(amrgrid%p_application, dof, var_centering, use_ghost, use_components, vector)
    case(LOCAL)
      use_ghost=PETSC_TRUE
      call create_samrai_vec(amrgrid%p_application, dof, var_centering, use_ghost, use_components, vector)
    case(NATURAL)
      call printErrMsg(option,'SAMRAI will not create PETSc Natural Vecs!!')
  end select

  if(dof>1) then
      call VecSetBlockSize(vector, dof)
  endif
      
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

subroutine AMRGridComputeGeometryInformation(amrgrid, origin_global, field, &
                                             option)
  use Option_module
  use Field_module
  
  implicit none

  interface
     PetscInt function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   PetscInt function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     PetscInt, intent(in) :: ln
   end function level_number_patches
  end interface

  type(amrgrid_type), pointer:: amrgrid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscReal :: origin_global(3)
  
#include "finclude/petscsysdef.h"
  PetscFortranAddr :: p_application

  type(grid_type), pointer :: grid
  DM :: dummy_ptr

  PetscInt :: nlevels
  PetscInt :: npatches
  PetscInt :: ln
  PetscInt :: pn

  dummy_ptr = PETSC_NULL

  call AMRGridComputeGridSpacing(amrgrid)
  
  p_application = amrgrid%p_application

  nlevels =  hierarchy_number_levels(p_application)

  do ln=0,nlevels-1
     npatches = level_number_patches(p_application, ln )
     do pn=0,npatches-1
        if (associated(amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr)) then
           grid => amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr
           call GridMapIndices(grid, dummy_ptr)
           call GridComputeCoordinates(grid,origin_global,option)
           call GridComputeVolumes(grid,field%volume,option)
           ! set up internal connectivity, distance, etc.
           call GridComputeInternalConnect(grid,option)
        endif
     end do
  end do

end subroutine AMRGridComputeGeometryInformation

subroutine AMRGridReadDXYZ(amrgrid, input, option)
  use Option_module
  use Input_Aux_module

  implicit none

  interface
     PetscInt function hierarchy_number_levels(p_hierarchy)
        PetscFortranAddr, intent(inout) :: p_hierarchy
     end function hierarchy_number_levels
   
     PetscInt function level_number_patches(p_hierarchy, ln)
       PetscFortranAddr, intent(inout) :: p_hierarchy
       PetscInt, intent(in) :: ln
     end function level_number_patches
 
     PetscInt function is_local_patch(p_hierarchy, ln, pn)
       PetscFortranAddr, intent(inout) :: p_hierarchy
       PetscInt, intent(in) :: ln
       PetscInt, intent(in) :: pn
     end function is_local_patch

  end interface

  type(amrgrid_type), pointer:: amrgrid
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(structured_grid_type),pointer :: structured_grid

#include "finclude/petscsysdef.h"
  PetscFortranAddr :: p_application
  PetscInt :: nlevels
  PetscInt :: npatches
  PetscInt :: ln
  PetscInt :: pn
  PetscInt :: islocal
  logical :: readdxyz
  PetscInt :: ndim    
  PetscReal :: dxdydz(3)
 
  p_application = amrgrid%p_application
  ! for now we assume uniform grid spacings 
  ndim=3    
  call AMRGridReadDXDYDZ(dxdydz,ndim, input,option)
  
  nlevels =  hierarchy_number_levels(p_application)
  do ln=0,nlevels-1
     npatches = level_number_patches(p_application, ln )
     do pn=0,npatches-1
        islocal = is_local_patch(p_application, ln, pn);
        if(islocal == 1) then
           grid => amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr
           structured_grid => grid%structured_grid           
           allocate(structured_grid%dx_global(structured_grid%nx))
           structured_grid%dx_global = dxdydz(1)
           allocate(structured_grid%dy_global(structured_grid%ny))
           structured_grid%dy_global = dxdydz(2)
           allocate(structured_grid%dz_global(structured_grid%nz))
           structured_grid%dz_global = dxdydz(3)
        endif
     end do
     ! we make the assumption of a refinement ratio of 2 in each direction   
     dxdydz(:) = dxdydz(:)/2.0   
  end do

end subroutine AMRGridReadDXYZ

! ************************************************************************** !
!
! StructGridReadArrayNew: Reads structured grid spacing along an axis from  
!                         input file
! author: Glenn Hammond
! date: 05/21/09
!
! ************************************************************************** !
subroutine AMRGridReadDXDYDZ(array,ndim,input,option)

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(input_type) :: input
  PetscInt :: ndim
  PetscReal :: array(ndim)
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: value
  PetscErrorCode :: ierr

  do i=1,ndim

      call InputReadFlotranString(input,option)
      call InputReadStringErrorMsg(input,option,'DXYZ')
      call InputReadWord(input,option,word,PETSC_TRUE)
      if (InputError(input)) exit
      string = word
      call InputReadDouble(string,option,value,input%ierr)
      call InputErrorMsg(input,option,'value','AMRGridReadDXDYDZ')
      array(i) = value
  enddo

end subroutine AMRGridReadDXDYDZ

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
