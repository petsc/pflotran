module Structured_Grid_module

  implicit none
 
  private
 
#include "definitions.h"

  type, public :: structured_grid_type

    character(len=MAXWORDLENGTH) :: ctype
    PetscInt :: itype  ! type of grid (e.g. structured, unstructured, etc.)
    PetscInt :: nx, ny, nz    ! Global domain dimensions of the grid.
    PetscInt :: nxy, nmax     ! nx * ny, nx * ny * nz
    PetscInt :: npx, npy, npz ! Processor partition in each direction.
    PetscInt :: npx_final, npy_final, npz_final ! actual decomposition
    PetscInt :: nlx, nly, nlz ! Local grid dimension w/o ghost nodes.
    PetscInt :: ngx, ngy, ngz ! Local grid dimension with ghost nodes.
    PetscInt :: nxs, nys, nzs 
      ! Global indices of non-ghosted corner (starting) of local domain.
    PetscInt :: ngxs, ngys, ngzs
      ! Global indices of ghosted starting corner of local domain.
    PetscInt :: nxe, nye, nze, ngxe, ngye, ngze
      ! Global indices of non-ghosted/ghosted ending corner of local domain.
    PetscInt :: nlxy, nlxz, nlyz
    PetscInt :: ngxy, ngxz, ngyz
    
    PetscInt :: istart, jstart, kstart, iend, jend, kend
      ! istart gives the ghosted local x-index of the non-ghosted starting 
      ! (lower left)corner. iend gives the local x-index of the non-ghosted 
      ! ending corner. jstart, jend correspond to y-index, kstart, kend to 
      ! z-index.  These are all zero-based indexing    

    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt :: nlmax_faces  ! Total number of non-ghosted faces in local domain.
    PetscInt :: ngmax_faces  ! Number of ghosted & non-ghosted faces in local domain.

    PetscReal :: origin(3) ! local origin of non-ghosted grid
    PetscReal :: bounds(3,2)

    ! grid spacing for each direction for global domain
    PetscReal, pointer :: dx_global(:), dy_global(:), dz_global(:)
    ! grid spacing for each direction for local, ghosted domain
    PetscReal, pointer :: dxg_local(:), dyg_local(:), dzg_local(:)
    
    PetscBool :: invert_z_axis
    
    PetscReal, pointer :: dx(:), dy(:), dz(:)  ! ghosted grid spacings for each grid cell
    
    PetscFortranAddr :: p_samr_patch ! pointer to a SAMRAI patch object

  end type structured_grid_type

  interface StructuredGridVecGetArrayF90
     module procedure StructuredGridVecGetArrayCellF90
     module procedure StructuredGridVecGetArraySideF90
  end interface

  public :: StructuredGridCreate, &
            StructuredGridDestroy, &
            StructuredGridCreateDM, &
            StructGridComputeLocalBounds, &
            StructGridComputeInternConnect, &
            StructGridComputeBoundConnect, &
            StructuredGridCreateVecFromDM, &
            StructuredGridMapIndices, &
            StructuredGridComputeSpacing, &
            StructuredGridComputeCoord, &
            StructuredGridReadDXYZ, &
            StructuredGridComputeVolumes, &
            StructGridPopulateConnection, &
            StructGridGetIJKFromCoordinate, &
            StructGridGetIJKFromLocalID, &
            StructGridGetIJKFromGhostedID, &
            StructuredGridVecGetMaskArrayCellF90, &
            StructuredGridVecGetArrayF90, &
            StructGridVecRestoreArrayF90, &
            StructGridGetGhostedNeighbors
contains

! ************************************************************************** !
!
! StructuredGridCreate: Creates a structured grid object
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
function StructuredGridCreate()

  implicit none
  
  type(structured_grid_type), pointer :: StructuredGridCreate

  type(structured_grid_type), pointer :: structured_grid

  allocate(structured_grid)
  
  structured_grid%ctype = ''
  structured_grid%itype = 0
  structured_grid%nx = 0
  structured_grid%ny = 0
  structured_grid%nz = 0
  structured_grid%npx = PETSC_DECIDE
  structured_grid%npy = PETSC_DECIDE
  structured_grid%npz = PETSC_DECIDE
  
  structured_grid%npx_final = 0
  structured_grid%npy_final = 0
  structured_grid%npz_final = 0

  structured_grid%nx = 0
  structured_grid%ny = 0
  structured_grid%nz = 0
  
  structured_grid%nxy = 0
  structured_grid%nmax = 0
  
  structured_grid%nlx = 0
  structured_grid%nly = 0
  structured_grid%nlz = 0
  structured_grid%nlxy = 0
  structured_grid%nlxz = 0
  structured_grid%nlyz = 0
  structured_grid%nlmax = 0
  
  structured_grid%ngx = 0
  structured_grid%ngy = 0
  structured_grid%ngz = 0
  structured_grid%ngxy = 0
  structured_grid%ngxz = 0
  structured_grid%ngyz = 0
  structured_grid%ngmax = 0

  structured_grid%nxs = 0
  structured_grid%nys = 0
  structured_grid%nzs = 0

  structured_grid%ngxs = 0
  structured_grid%ngys = 0
  structured_grid%ngzs = 0

  structured_grid%nxe = 0
  structured_grid%nye = 0
  structured_grid%nze = 0

  structured_grid%ngxe = 0
  structured_grid%ngye = 0
  structured_grid%ngze = 0

  structured_grid%istart = 0
  structured_grid%jstart = 0
  structured_grid%kstart = 0
  structured_grid%iend = 0
  structured_grid%jend = 0
  structured_grid%kend = 0
  
  nullify(structured_grid%dx_global)
  nullify(structured_grid%dy_global)
  nullify(structured_grid%dz_global)

  nullify(structured_grid%dxg_local)
  nullify(structured_grid%dyg_local)
  nullify(structured_grid%dzg_local)

  nullify(structured_grid%dx)
  nullify(structured_grid%dy)
  nullify(structured_grid%dz)
  
  
  structured_grid%origin = -1.d20
  structured_grid%bounds = -1.d20
  
  structured_grid%invert_z_axis = PETSC_FALSE
  
  structured_grid%p_samr_patch=0

  StructuredGridCreate => structured_grid
  
end function StructuredGridCreate
  
! ************************************************************************** !
!
! StructuredGridCreateDMs: Creates structured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine StructuredGridCreateDM(structured_grid,da,ndof,stencil_width, &
                                  option)

  use Option_module
        
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#ifndef DMDA_OLD
! For PETSc versions >= 3.2
#include "finclude/petscdmda.h"
#endif

  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  DM :: da
  PetscInt :: ndof
  PetscInt :: stencil_width

  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DM object that will manage communication.
  !-----------------------------------------------------------------------
#ifndef DMDA_OLD
  ! This code is for the DMDACreate3D() interface in PETSc versions >= 3.2 --RTM
  call DMDACreate3D(option%mycomm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, &
                  DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  da,ierr)
  call DMDAGetInfo(da,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,structured_grid%npx_final, &
                 structured_grid%npy_final,structured_grid%npz_final, &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,ierr)
#else
  ! This code is for the DMDACreate3D() interface in versions of PETSc
  ! prior to release 3.2.  This should be removed, eventually. --RTM
  call DMDACreate3D(option%mycomm,DMDA_NONPERIODIC,DMDA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  da,ierr)
  call DMDAGetInfo(da,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,structured_grid%npx_final, &
                 structured_grid%npy_final,structured_grid%npz_final, &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,ierr)
#endif


end subroutine StructuredGridCreateDM

! ************************************************************************** !
!
! StructGridComputeLocalBounds: Computes the corners for the local portion
!                               of the structured grid
! author: Glenn Hammond
! date: 03/13/08
!
! ************************************************************************** !
subroutine StructGridComputeLocalBounds(structured_grid,da)

implicit none

interface
subroutine samr_patch_get_corners(p_patch, nxs, nys, nzs, nlx, nly, nlz)
       implicit none
       
#include "finclude/petscsysdef.h"

       PetscFortranAddr :: p_patch
       PetscInt :: nxs, nys, nzs, nlx, nly, nlz

end subroutine samr_patch_get_corners
     
subroutine samr_patch_get_ghostcorners(p_patch, nxs, nys, nzs, nlx, nly, nlz)
       implicit none
       
#include "finclude/petscsysdef.h"
       
       PetscFortranAddr :: p_patch
       PetscInt :: nxs, nys, nzs, nlx, nly, nlz

     end subroutine samr_patch_get_ghostcorners
  end interface
     
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"

  type(structured_grid_type) :: structured_grid
  DM :: da

  PetscErrorCode :: ierr

  if(structured_grid%p_samr_patch==0) then
      ! get corner information
     call DMDAGetCorners(da, structured_grid%nxs, &
          structured_grid%nys, structured_grid%nzs, structured_grid%nlx, &
          structured_grid%nly, structured_grid%nlz, ierr)
     
     structured_grid%nxe = structured_grid%nxs + structured_grid%nlx
     structured_grid%nye = structured_grid%nys + structured_grid%nly
     structured_grid%nze = structured_grid%nzs + structured_grid%nlz
     structured_grid%nlxy = structured_grid%nlx * structured_grid%nly
     structured_grid%nlxz = structured_grid%nlx * structured_grid%nlz
     structured_grid%nlyz = structured_grid%nly * structured_grid%nlz
     structured_grid%nlmax = structured_grid%nlx * structured_grid%nly * structured_grid%nlz
     
     ! get ghosted corner information
     call DMDAGetGhostCorners(da, structured_grid%ngxs, &
          structured_grid%ngys, structured_grid%ngzs, structured_grid%ngx, &
          structured_grid%ngy, structured_grid%ngz, ierr)
     
     structured_grid%ngxe = structured_grid%ngxs + structured_grid%ngx
     structured_grid%ngye = structured_grid%ngys + structured_grid%ngy
     structured_grid%ngze = structured_grid%ngzs + structured_grid%ngz
     structured_grid%ngxy = structured_grid%ngx * structured_grid%ngy
     structured_grid%ngxz = structured_grid%ngx * structured_grid%ngz
     structured_grid%ngyz = structured_grid%ngy * structured_grid%ngz
     structured_grid%ngmax = structured_grid%ngx * structured_grid%ngy * structured_grid%ngz
  else
     ! get corner information
     call samr_patch_get_corners(structured_grid%p_samr_patch, &
          structured_grid%nxs, structured_grid%nys, structured_grid%nzs, &
          structured_grid%nlx, structured_grid%nly, structured_grid%nlz)
     
     structured_grid%nxe = structured_grid%nxs + structured_grid%nlx
     structured_grid%nye = structured_grid%nys + structured_grid%nly
     structured_grid%nze = structured_grid%nzs + structured_grid%nlz
     structured_grid%nlxy = structured_grid%nlx * structured_grid%nly
     structured_grid%nlxz = structured_grid%nlx * structured_grid%nlz
     structured_grid%nlyz = structured_grid%nly * structured_grid%nlz
     structured_grid%nlmax = structured_grid%nlx * structured_grid%nly * structured_grid%nlz
     
     ! get ghosted corner information
     call samr_patch_get_ghostcorners(structured_grid%p_samr_patch, &
          structured_grid%ngxs, structured_grid%ngys, structured_grid%ngzs, &
          structured_grid%ngx, structured_grid%ngy, structured_grid%ngz)
     structured_grid%ngxe = structured_grid%ngxs + structured_grid%ngx
     structured_grid%ngye = structured_grid%ngys + structured_grid%ngy
     structured_grid%ngze = structured_grid%ngzs + structured_grid%ngz
     structured_grid%ngxy = structured_grid%ngx * structured_grid%ngy
     structured_grid%ngxz = structured_grid%ngx * structured_grid%ngz
     structured_grid%ngyz = structured_grid%ngy * structured_grid%ngz
     structured_grid%ngmax = structured_grid%ngx * structured_grid%ngy * structured_grid%ngz
   endif
   
end subroutine StructGridComputeLocalBounds

! ************************************************************************** !
!
! StructuredGridCreateVecFromDM: Creates a PETSc vector from a DM
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine StructuredGridCreateVecFromDM(da,vector,vector_type)

  implicit none

  DM :: da
  Vec :: vector
  PetscInt :: vector_type
  
  PetscErrorCode :: ierr

  select case (vector_type)
    case(GLOBAL)
      call DMCreateGlobalVector(da,vector,ierr)
    case(LOCAL)
      call DMCreateLocalVector(da,vector,ierr)
    case(NATURAL)
      call DMDACreateNaturalVector(da,vector,ierr)
  end select

end subroutine StructuredGridCreateVecFromDM

! ************************************************************************** !
!
! StructuredGridReadDXYZ: Reads structured grid spacing from input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine StructuredGridReadDXYZ(structured_grid,input,option)

  use Option_module
  use Input_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt :: i

  allocate(structured_grid%dx_global(structured_grid%nx))
  structured_grid%dx_global = 0.d0
  allocate(structured_grid%dy_global(structured_grid%ny))
  structured_grid%dy_global = 0.d0
  allocate(structured_grid%dz_global(structured_grid%nz))
  structured_grid%dz_global = 0.d0

  word = 'X'
  call StructuredGridReadArrayNew(structured_grid%dx_global, &
                               structured_grid%nx,word,input,option)
  word = 'Y'
  call StructuredGridReadArrayNew(structured_grid%dy_global, &
                               structured_grid%ny,word,input,option)
  word = 'Z'
  call StructuredGridReadArrayNew(structured_grid%dz_global, &
                               structured_grid%nz,word,input,option)
    
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'(/," *DXYZ ")')
    write(option%fid_out,'("  dx  ",/,(1p10e12.4))') &
      (structured_grid%dx_global(i),i=1,structured_grid%nx)
    write(option%fid_out,'("  dy  ",/,(1p10e12.4))') &
      (structured_grid%dy_global(i),i=1,structured_grid%ny)
    write(option%fid_out,'("  dz  ",/,(1p10e12.4))') &
      (structured_grid%dz_global(i),i=1,structured_grid%nz)
  endif

end subroutine StructuredGridReadDXYZ

! ************************************************************************** !
!
! StructuredGridReadArray: Reads structured grid spacing along an axis from  
!                         input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine StructuredGridReadArray(a,n,input,option)

  use Input_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(input_type) :: input
  PetscInt :: fid
  PetscInt :: n
  PetscInt :: i, i1, i2, m
  PetscInt, parameter ::  nvalue=10
  PetscReal, intent(inout) :: a(*)
  character(len=MAXSTRINGLENGTH) :: string 

!  call fiReadStringErrorMsg('DXYZ',ierr)

!  call fiReadDouble(string,grid%radius_0,ierr)
!  call fiDefaultMsg(option%myrank,'radius_0',ierr)

  i2 = 0
  do
    i1 = i2+1
    i2 = i2+nvalue
    if (i2.gt.n) i2 = n
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'DXYZ')
    do i = i1, i2
      call InputReadDouble(input,option,a(i))
!geh  ierr, which comes from iostat, will not necessarily be 0 and 1, 
!geh  it could be another number
!geh      if (ierr == 1) a(i) = 0.d0
      if (input%ierr /= 0) a(i) = 0.d0
!     print *,i,i1,i2,nvalue,a(i),n,ierr
!     call fiDefaultMsg("Error reading grid spacing", ierr)
    enddo
    do i = i1,i2
      if (a(i).eq.0.d0) then

!---------if less than nx non-zero values are read, set all the zero
!         values to the last non zero value read. Only for cartesian 
!         system

        do m = i,n
          a(m) = a(i-1)
        enddo
        return
      endif
    enddo
    if (i2 >= n) exit
  enddo
    
end subroutine StructuredGridReadArray

! ************************************************************************** !
!
! StructuredGridReadArrayNew: Reads structured grid spacing along an axis from  
!                         input file
! author: Glenn Hammond
! date: 05/21/09
!
! ************************************************************************** !
subroutine StructuredGridReadArrayNew(array,array_size,axis,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(input_type) :: input
  character(len=MAXWORDLENGTH) :: axis
  PetscInt :: array_size
  PetscReal :: array(array_size)
  
  PetscInt :: i, num_values, count
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word, word2, word3
  character(len=1) :: backslash
  PetscBool :: continuation_flag
  PetscReal :: value
  PetscErrorCode :: ierr

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  count = 0
  continuation_flag = PETSC_TRUE
  do
    
    if (count >= array_size) exit
    
    if (.not.continuation_flag .and. count /= 1) then
      option%io_buffer = 'Insufficient values read for ' // &
                         trim(axis) // &
                         ' direction in DXYZ card'
      call printErrMsg(option)
    else if (count == 1) then
      array = array(count)
      exit
    endif
    
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'DXYZ')

    continuation_flag = PETSC_FALSE
    if (index(input%buf,backslash) > 0) &
      continuation_flag = PETSC_TRUE

    do 
      call InputReadWord(input,option,word,PETSC_TRUE)
      if (InputError(input) .or. StringCompare(word,backslash,ONE_INTEGER)) exit
      i = index(word,'*')
      if (i == 0) i = index(word,'@')
      if (i /= 0) then
        word2 = word(1:i-1)
        word3 = word(i+1:len_trim(word))
        string2 = word2
        call InputReadInt(string2,option,num_values,input%ierr)
        call InputErrorMsg(input,option,'# values','StructuredGridReadArrayNew')
        string2 = word3
        call InputReadDouble(string2,option,value,input%ierr)
        call InputErrorMsg(input,option,'value','StructuredGridReadArrayNew')
        do i=1, num_values
          count = count + 1
          array(count) = value
        enddo
      else
        string2 = word
        call InputReadDouble(string2,option,value,input%ierr)
        call InputErrorMsg(input,option,'value','StructuredGridReadDXYZ')
        count = count + 1
        array(count) = value
      endif
    enddo
  enddo

end subroutine StructuredGridReadArrayNew

! ************************************************************************** !
!
! StructuredGridComputeSpacing: Computes structured grid spacing
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine StructuredGridComputeSpacing(structured_grid,nG2A,nG2L,option)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  PetscInt :: nG2A(:)
  PetscInt :: nG2L(:)
  type(option_type) :: option
  
  PetscInt :: i, j, k, ghosted_id
  PetscErrorCode :: ierr

  allocate(structured_grid%dxg_local(structured_grid%ngx))
  structured_grid%dxg_local = 0.d0
  allocate(structured_grid%dyg_local(structured_grid%ngy))
  structured_grid%dyg_local = 0.d0
  allocate(structured_grid%dzg_local(structured_grid%ngz))
  structured_grid%dzg_local = 0.d0
  
  if (structured_grid%p_samr_patch == 0) then
    if (.not.associated(structured_grid%dx_global)) then
      ! indicates that the grid spacings still need to be computed
      if (structured_grid%bounds(1,1) < -1.d19) then ! bounds have not been initialized
        call printErrMsg(option,'Bounds have not been set for grid and DXYZ does not exist')
      endif
      allocate(structured_grid%dx_global(structured_grid%nx))
      allocate(structured_grid%dy_global(structured_grid%ny))
      allocate(structured_grid%dz_global(structured_grid%nz))
      
      select case(structured_grid%itype)
        case(CARTESIAN_GRID)
          structured_grid%dx_global = (structured_grid%bounds(X_DIRECTION,UPPER)- &
                                       structured_grid%bounds(X_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nx)
          structured_grid%dy_global = (structured_grid%bounds(Y_DIRECTION,UPPER)- &
                                       structured_grid%bounds(Y_DIRECTION,LOWER)) / &
                                       dble(structured_grid%ny)
          structured_grid%dz_global = (structured_grid%bounds(Z_DIRECTION,UPPER)- &
                                       structured_grid%bounds(Z_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nz)
        case(CYLINDRICAL_GRID)
          structured_grid%dx_global = (structured_grid%bounds(X_DIRECTION,UPPER)- &
                                       structured_grid%bounds(X_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nx)
          structured_grid%dy_global = 1.d0
          structured_grid%dz_global = (structured_grid%bounds(Z_DIRECTION,UPPER)- &
                                       structured_grid%bounds(Z_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nz)
        case(SPHERICAL_GRID)
          structured_grid%dx_global = (structured_grid%bounds(X_DIRECTION,UPPER)- &
                                       structured_grid%bounds(X_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nx)
          structured_grid%dy_global = 1.d0
          structured_grid%dz_global = 1.d0
      end select
      
    endif
    structured_grid%dxg_local(1:structured_grid%ngx) = &
      structured_grid%dx_global(structured_grid%ngxs+1:structured_grid%ngxe)
    structured_grid%dyg_local(1:structured_grid%ngy) = &
      structured_grid%dy_global(structured_grid%ngys+1:structured_grid%ngye)
    structured_grid%dzg_local(1:structured_grid%ngz) = &
      structured_grid%dz_global(structured_grid%ngzs+1:structured_grid%ngze)
  endif    
        
  allocate(structured_grid%dx(structured_grid%ngmax))
  structured_grid%dx = 0.d0
  allocate(structured_grid%dy(structured_grid%ngmax))
  structured_grid%dy = 0.d0
  allocate(structured_grid%dz(structured_grid%ngmax))
  structured_grid%dz = 0.d0
 
  do k = 1, structured_grid%ngz
    do j = 1, structured_grid%ngy
      do i = 1, structured_grid%ngx
        ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
        structured_grid%dx(ghosted_id) = structured_grid%dxg_local(i)
        structured_grid%dy(ghosted_id) = structured_grid%dyg_local(j)
        structured_grid%dz(ghosted_id) = structured_grid%dzg_local(k)
      enddo
    enddo
  enddo
  
end subroutine StructuredGridComputeSpacing

! ************************************************************************** !
!
! StructuredGridComputeCoord: Computes structured coordinates in x,y,z
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridComputeCoord(structured_grid,option,origin_global, &
                                      grid_x,grid_y,grid_z, &
                                      x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module
  
implicit none

interface
subroutine samr_patch_get_origin(p_patch, xs, ys, zs)
       implicit none
#include "finclude/petscsysdef.h"
       PetscFortranAddr, intent(inout) :: p_patch
       PetscReal, intent(inout) :: xs
       PetscReal, intent(inout) :: ys
       PetscReal, intent(inout) :: zs
     end subroutine samr_patch_get_origin
  end interface
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscReal :: origin_global(3)
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

! PetscErrorCode :: ierr
  PetscInt :: i, j, k, ghosted_id
  PetscReal :: x, y, z
  PetscInt :: prevnode

  if (structured_grid%p_samr_patch == 0) then
    x_min = origin_global(X_DIRECTION)
    y_min = origin_global(Y_DIRECTION)
    z_min = origin_global(Z_DIRECTION)
    
    do i=1,structured_grid%nxs
      x_min = x_min + structured_grid%dx_global(i)
    enddo
    do j=1,structured_grid%nys
      y_min = y_min + structured_grid%dy_global(j)
    enddo
    do k=1,structured_grid%nzs
      z_min = z_min + structured_grid%dz_global(k)
    enddo
  else
    call samr_patch_get_origin(structured_grid%p_samr_patch, x_min, y_min, z_min)
  endif
  
  ! set min and max bounds of domain in coordinate directions
  structured_grid%origin(X_DIRECTION) = x_min
  structured_grid%origin(Y_DIRECTION) = y_min
  structured_grid%origin(Z_DIRECTION) = z_min
  x_max = x_min
  y_max = y_min
  z_max = z_min
  
  do i=structured_grid%istart,structured_grid%iend
    x_max = x_max + structured_grid%dxg_local(i+1)
  enddo
  do j=structured_grid%jstart,structured_grid%jend
    y_max = y_max + structured_grid%dyg_local(j+1)
  enddo
  do k=structured_grid%kstart,structured_grid%kend
    z_max = z_max + structured_grid%dzg_local(k+1)
  enddo

! fill in grid cell coordinates
  ghosted_id = 0
  if (structured_grid%kstart > 0 .or. structured_grid%p_samr_patch /= 0) then
    z = -0.5d0*structured_grid%dzg_local(1)+structured_grid%origin(Z_DIRECTION)
  else
    z = 0.5d0*structured_grid%dzg_local(1)+structured_grid%origin(Z_DIRECTION)
  endif
  do k=1, structured_grid%ngz
    if (structured_grid%jstart > 0 .or. structured_grid%p_samr_patch /= 0) then
      y = -0.5d0*structured_grid%dyg_local(1)+structured_grid%origin(Y_DIRECTION)
    else
      y = 0.5d0*structured_grid%dyg_local(1)+structured_grid%origin(Y_DIRECTION)
    endif
    do j=1, structured_grid%ngy
      if (structured_grid%istart > 0 .or. structured_grid%p_samr_patch /= 0) then
        x = -0.5d0*structured_grid%dxg_local(1)+structured_grid%origin(X_DIRECTION)
      else
        x = 0.5d0*structured_grid%dxg_local(1)+structured_grid%origin(X_DIRECTION)
      endif
      do i=1, structured_grid%ngx
        ghosted_id = ghosted_id + 1
        grid_x(ghosted_id) = x
        grid_y(ghosted_id) = y
        grid_z(ghosted_id) = z
        if (i < structured_grid%ngx) &
          x = x + 0.5d0*(structured_grid%dxg_local(i)+structured_grid%dxg_local(i+1))
      enddo
      if (j < structured_grid%ngy) &
        y = y + 0.5d0*(structured_grid%dyg_local(j)+structured_grid%dyg_local(j+1))
    enddo
    if (k < structured_grid%ngz) &
      z = z + 0.5d0*(structured_grid%dzg_local(k)+structured_grid%dzg_local(k+1))
  enddo
    
end subroutine StructuredGridComputeCoord

! ************************************************************************** !
!
! StructGridGetIJKFromCoordinate: Finds local, non-ghosted i,j,k indices for  
!                                 grid cell encompassing coordinate
! author: Glenn Hammond
! date: 08/27/08
!
! ************************************************************************** !
subroutine StructGridGetIJKFromCoordinate(structured_grid,x,y,z,i,j,k)

  use Option_module
  
  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
#include "finclude/petscsysdef.h"
       PetscFortranAddr :: p_patch
       PetscInt :: axis,dim
     end function samr_patch_at_bc
  end interface
    
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: i, j, k
  PetscInt :: i_local, j_local, k_local
  PetscInt :: i_ghosted, j_ghosted, k_ghosted
  PetscReal :: x, y, z
  
  PetscReal :: x_upper_face, y_upper_face, z_upper_face

  i = -1
  j = -1
  k = -1

  x_upper_face = structured_grid%origin(X_DIRECTION)
  i_local = 1
  do i_ghosted=structured_grid%istart,structured_grid%iend
    if (x >= x_upper_face .and. &                   ! since i_ghosted is zero-based
        x <= x_upper_face+structured_grid%dxg_local(i_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in x-dir on proc
      if (i_ghosted == structured_grid%istart) then
        ! located on upwind boundary and ghosted
        if (structured_grid%p_samr_patch == 0) then
          if (x == x_upper_face .and. &
              structured_grid%nxs /= structured_grid%ngxs) exit
! geh - the below should not be necessary for AMR since patches do not span processors              
!        else if (samr_patch_at_bc(structured_grid%p_samr_patch, 0, 0) == 0) then
!          exit
        endif
      endif
      i = i_local 
      exit
    endif
    i_local = i_local + 1
    x_upper_face = x_upper_face + structured_grid%dxg_local(i_ghosted+1)
  enddo
  y_upper_face = structured_grid%origin(Y_DIRECTION)
  j_local = 1
  do j_ghosted=structured_grid%jstart,structured_grid%jend
    if (y >= y_upper_face .and. &
        y <= y_upper_face+structured_grid%dyg_local(j_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in y-dir on proc
      if (j_ghosted == structured_grid%jstart) then
        ! located on upwind boundary and ghosted
        if (structured_grid%p_samr_patch == 0) then
          if (y == y_upper_face .and. &
              structured_grid%nys /= structured_grid%ngys) exit
!        else if (samr_patch_at_bc(structured_grid%p_samr_patch, 1, 0) == 0) then
!          exit
        endif
      endif
      j = j_local
      exit
    endif
    j_local = j_local + 1
    y_upper_face = y_upper_face + structured_grid%dyg_local(j_ghosted+1)
  enddo
  z_upper_face = structured_grid%origin(Z_DIRECTION)
  k_local = 1
  do k_ghosted=structured_grid%kstart,structured_grid%kend
    if (z >= z_upper_face .and. &
        z <= z_upper_face+structured_grid%dzg_local(k_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in z-dir on proc
      if (k_ghosted == structured_grid%kstart) then
        ! if located on upwind boundary and ghosted, skip
        if (structured_grid%p_samr_patch == 0) then
          if (z == z_upper_face .and. &
            structured_grid%nzs /= structured_grid%ngzs) exit
!        else if (samr_patch_at_bc(structured_grid%p_samr_patch, 2, 0) == 0) then
!          exit
        endif
      endif          
      k = k_local
      exit
    endif
    k_local = k_local + 1
    z_upper_face = z_upper_face + structured_grid%dzg_local(k_ghosted+1)
  enddo
    
end subroutine StructGridGetIJKFromCoordinate

! ************************************************************************** !
!
! StructGridGetIJKFromLocalID: Finds i,j,k indices for grid cell defined by 
!                              local_id
! author: Glenn Hammond
! date: 04/11/08
!
! ************************************************************************** !
subroutine StructGridGetIJKFromLocalID(structured_grid,local_id,i,j,k)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: local_id
  PetscReal :: i, j, k
  
  k= int((local_id-1)/structured_grid%nlxy) + 1
  j= int(mod((local_id-1),structured_grid%nlxy)/structured_grid%nlx) + 1
  i= mod(mod((local_id-1),structured_grid%nlxy),structured_grid%nlx) + 1  
  
end subroutine StructGridGetIJKFromLocalID

! ************************************************************************** !
!
! StructGridGetIJKFromGhostedID: Finds i,j,k indices for grid cell defined by 
!                                a ghosted id
! author: Glenn Hammond
! date: 04/11/08
!
! ************************************************************************** !
subroutine StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscInt :: i, j, k
  
  k= int((ghosted_id-1)/structured_grid%ngxy) + 1
  j= int(mod((ghosted_id-1),structured_grid%ngxy)/structured_grid%ngx) + 1
  i= mod(mod((ghosted_id-1),structured_grid%ngxy),structured_grid%ngx) + 1  
  
end subroutine StructGridGetIJKFromGhostedID

! ************************************************************************** !
!
! StructGridGetLocalIDFromIJK: Finds local_id for grid cell defined by 
!                              i,j,k indices
! author: Glenn Hammond
! date: 01/28/11
!
! ************************************************************************** !
function StructGridGetLocalIDFromIJK(structured_grid,i,j,k)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: i, j, k
  
  PetscInt :: StructGridGetLocalIDFromIJK
  
  StructGridGetLocalIDFromIJK = &
    i+(j-1)*structured_grid%nlx+(k-1)*structured_grid%nlxy
  
end function StructGridGetLocalIDFromIJK

! ************************************************************************** !
!
! StructGridGetGhostedIDFromIJK: Finds ghosted_id for grid cell defined by 
!                              i,j,k indices
! author: Glenn Hammond
! date: 01/28/11
!
! ************************************************************************** !
function StructGridGetGhostedIDFromIJK(structured_grid,i,j,k)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: i, j, k
  
  PetscInt :: StructGridGetGhostedIDFromIJK
  
  StructGridGetGhostedIDFromIJK = &
    i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
  
end function StructGridGetGhostedIDFromIJK

! ************************************************************************** !
!
! StructGridComputeInternConnect: computes internal connectivity of a  
!                               structured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function StructGridComputeInternConnect(structured_grid, xc, yc, zc, option)

  use Connection_module
  use Option_module
  
  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
#include "finclude/petscsysdef.h"
       PetscFortranAddr :: p_patch
       PetscInt :: axis,dim
     end function samr_patch_at_bc
  end interface

!  PetscReal :: radius(:)
  type(connection_set_type), pointer :: StructGridComputeInternConnect
  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  PetscReal, pointer :: xc(:),yc(:),zc(:)
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  PetscInt :: i, j, k, iconn, id_up, id_dn
  PetscInt :: samr_ofx, samr_ofy, samr_ofz
  PetscInt :: nconn
  PetscInt :: lenx, leny, lenz
  PetscReal :: dist_up, dist_dn
  PetscReal :: r1, r2
  type(connection_set_type), pointer :: connections, connections_2

  PetscErrorCode :: ierr
  PetscReal, pointer :: radius(:)
  PetscInt, pointer :: int_array1(:), int_array2(:),int_array3(:),int_array4(:),int_array5(:),index(:)
  PetscInt :: count

  radius => xc
  
  samr_ofx = 0
  samr_ofy = 0
  samr_ofz = 0
  ! the adjustments in the case of AMR are based on the PIMS code adjustments by LC
  nconn = (structured_grid%ngx-1)*structured_grid%nly*structured_grid%nlz+ &
          structured_grid%nlx*(structured_grid%ngy-1)*structured_grid%nlz+ &
          structured_grid%nlx*structured_grid%nly*(structured_grid%ngz-1)
  

  structured_grid%nlmax_faces = 0
  structured_grid%ngmax_faces = 0

  lenx = structured_grid%ngx - 1
  leny = structured_grid%ngy - 1
  lenz = structured_grid%ngz - 1

  if(.not.(structured_grid%p_samr_patch == 0)) then
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ZERO_INTEGER, ZERO_INTEGER) ==1) then
        nconn = nconn - structured_grid%nlyz
        lenx = lenx-1
        samr_ofx = 1
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ZERO_INTEGER, ONE_INTEGER) ==1) then 
        nconn = nconn - structured_grid%nlyz
        lenx = lenx-1
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ONE_INTEGER, ZERO_INTEGER) ==1) then
        nconn = nconn - structured_grid%nlxz
        leny=leny-1
        samr_ofy = structured_grid%ngx
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ONE_INTEGER, ONE_INTEGER) ==1) then
        nconn = nconn - structured_grid%nlxz
        leny=leny-1
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, TWO_INTEGER, ZERO_INTEGER) ==1) then
        nconn = nconn - structured_grid%nlxy
        lenz=lenz-1
        samr_ofz = structured_grid%ngxy
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, TWO_INTEGER, ONE_INTEGER) ==1) then
        nconn = nconn - structured_grid%nlxy
        lenz=lenz-1
     endif  
  endif



  connections => ConnectionCreate(nconn, &
                                  option%nphase,INTERNAL_CONNECTION_TYPE)

  iconn = 0
  ! x-connections
  if (structured_grid%ngx > 1) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
        
              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
               
 
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy+samr_ofx
              id_dn = id_up + 1


              
              if (i==1) then
                if (structured_grid%nxs==structured_grid%ngxs) then
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
                end if
              else 
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
              end if


              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dx(id_up)
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = structured_grid%dy(id_up)* &
                                        structured_grid%dz(id_up)
              connections%cntr(1,iconn) = xc(id_up) + &
                  connections%dist(-1,iconn)*(xc(id_dn) - xc(id_up))
              connections%cntr(2,iconn) = yc(id_up)
              connections%cntr(3,iconn) = zc(id_up)
              
            enddo
          enddo
        enddo
      case(CYLINDRICAL_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy+samr_ofx
              id_dn = id_up + 1
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dx(id_up)
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = 2.d0 * pi * (radius(id_up)+0.5d0*structured_grid%dx(id_up))* &
                                        structured_grid%dz(id_up)
              connections%cntr(1,iconn) = xc(id_up) + connections%dist(-1,iconn)*(xc(id_dn) - xc(id_up))
              connections%cntr(2,iconn) = yc(id_up)
              connections%cntr(3,iconn) = zc(id_up)
            enddo
          enddo
        enddo
      case(SPHERICAL_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy+samr_ofx
              id_dn = id_up + 1
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dx(id_up)
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = 4.d0 * pi * (radius(id_up)+0.5d0*structured_grid%dx(id_up))
              connections%cntr(1,iconn) = xc(id_up) + connections%dist(-1,iconn)*(xc(id_dn) - xc(id_up))
              connections%cntr(2,iconn) = yc(id_up)
              connections%cntr(3,iconn) = zc(id_up)
            enddo
          enddo
        enddo
  end select
    
  endif

  ! y-connections
  if (structured_grid%ngy > 1) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do i = structured_grid%istart, structured_grid%iend
            do j = 1, leny
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              if (j==1) then
                if (structured_grid%nys==structured_grid%ngys) then
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
                end if
              else 
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
              end if

              id_up = i + 1 + (j-1) * structured_grid%ngx + k * structured_grid%ngxy &
                  +samr_ofy
              id_dn = id_up + structured_grid%ngx
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dy(id_up)
              dist_dn = 0.5d0*structured_grid%dy(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(2,iconn) = 1.d0  ! y component of unit vector
              connections%area(iconn) = structured_grid%dx(id_up)* &
                                    structured_grid%dz(id_up)
              connections%cntr(1,iconn) = xc(id_up) 
              connections%cntr(2,iconn) = yc(id_up) + connections%dist(-1,iconn)*(yc(id_dn) - yc(id_up)) 
              connections%cntr(3,iconn) = zc(id_up)
            enddo
          enddo
        enddo
      case(CYLINDRICAL_GRID)
        option%io_buffer = 'For cylindrical coordinates, NY must be equal to 1.'
        call printErrMsg(option)
      case(SPHERICAL_GRID)
        option%io_buffer = 'For spherical coordinates, NY must be equal to 1.'
        call printErrMsg(option)

    end select
  endif
      
  ! z-connections
  if (structured_grid%ngz > 1) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        do j = structured_grid%jstart, structured_grid%jend
          do i = structured_grid%istart, structured_grid%iend
            do k = 1, lenz
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              if (k==1) then
                if (structured_grid%nzs==structured_grid%ngzs) then
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
                end if
              else 
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
              end if

              id_up = i + 1 + j * structured_grid%ngx + (k-1) * &
                  structured_grid%ngxy + samr_ofz
              id_dn = id_up + structured_grid%ngxy
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dz(id_up)
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(3,iconn) = 1.d0  ! z component of unit vector
              connections%area(iconn) = structured_grid%dx(id_up) * &
                                        structured_grid%dy(id_up)
              connections%cntr(1,iconn) = xc(id_up) 
              connections%cntr(2,iconn) = yc(id_up) 
              connections%cntr(3,iconn) = zc(id_up) + connections%dist(-1,iconn)*(zc(id_dn) - zc(id_up)) 
            enddo
          enddo
        enddo
      case(CYLINDRICAL_GRID)
        do j = structured_grid%jstart, structured_grid%jend
          do i = structured_grid%istart, structured_grid%iend
            do k = 1, lenz
              iconn = iconn+1
              id_up = i + 1 + j * structured_grid%ngx + (k-1) * &
                  structured_grid%ngxy + samr_ofz
              id_dn = id_up + structured_grid%ngxy
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dz(id_up)
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(3,iconn) = 1.d0  ! z component of unit vector
              ! pi*(r2^2-r1^2)
              r2 = xc(id_up) + 0.5d0*structured_grid%dx(id_up)
              r1 = xc(id_up) - 0.5d0*structured_grid%dx(id_up)
              connections%area(iconn) = pi * dabs(r2*r2 - r1*r1)
              connections%cntr(1,iconn) = xc(id_up) 
              connections%cntr(2,iconn) = yc(id_up) 
              connections%cntr(3,iconn) = zc(id_up) + &
                                          connections%dist(-1,iconn)* &
                                          (zc(id_dn) - zc(id_up)) 
            enddo
          enddo
        enddo
      case(SPHERICAL_GRID)
        option%io_buffer = 'For spherical coordinates, NZ must be equal to 1.'
        call printErrMsg(option)
  end select
  endif

#ifdef REARRANGE_CONN
  allocate(int_array1(1:iconn))
  allocate(int_array2(1:iconn))
  allocate(int_array3(1:iconn))
  allocate(int_array4(1:iconn))
  allocate(int_array5(1:iconn))
  allocate(index(1:iconn))

  do i = 1,iconn
    int_array1(i) = i
    int_array2(i) = connections%id_up(i)
  enddo
  
  int_array1 = int_array1 - 1
  call PetscSortIntWithPermutation(iconn, int_array2, int_array1,ierr)
  int_array1 = int_array1 + 1

  count = 0
  i = 1
  count = count + 1
  int_array3(count) = int_array2(int_array1(i))
  int_array4(count) = connections%id_dn(int_array1(i))

  do i=2,iconn
    if( int_array3(count).ne.int_array2(int_array1(i) )) then
      do k = 1,count
        int_array5(k) = k
      enddo
      int_array5 = int_array5 - 1
      call PetscSortIntWithPermutation(count,int_array4,int_array5,ierr)
      int_array5 = int_array5 + 1
      do k = 1,count
        index(i -count +k -1) = int_array1(i -count -1 + int_array5(k) )
      enddo
      count = 1
      int_array3(count) = int_array2(int_array1(i))
      int_array4(count) = connections%id_dn(int_array1(i))
    else
      count = count + 1
      int_array3(count) = int_array2(int_array1(i))
      int_array4(count) = connections%id_dn(int_array1(i))
    endif
  enddo
    do k = 1,count
      int_array5(k) = k
    enddo
    int_array5 = int_array5 - 1
    call PetscSortIntWithPermutation(count,int_array4,int_array5,ierr)
    int_array5 = int_array5 + 1
    do k = 1,count
      index(i -count +k -1) = int_array1(i -count -1 + int_array5(k) )
    enddo
  connections_2=> ConnectionCreate(nconn, &
                                  option%nphase,INTERNAL_CONNECTION_TYPE)
  do i=1,iconn
    connections_2%local(i)     = connections%local(index(i))
    connections_2%id_up(i)     = connections%id_up(index(i))
    connections_2%id_dn(i)     = connections%id_dn(index(i))
    connections_2%dist(-1:3,i) = connections%dist(-1:3,index(i))
    connections_2%area(i)      = connections%area(index(i))
    connections_2%cntr(1:3,i)  = connections%cntr(1:3,index(i))
  enddo
  StructGridComputeInternConnect => connections_2
#else
  StructGridComputeInternConnect => connections
#endif

end function StructGridComputeInternConnect


! ************************************************************************** !
!
! StructGridComputeBoundConnect: computes boundary connectivity of a  
!                               structured grid
! author: Daniil Svyatskiy
! date: 02/04/10
!
! ************************************************************************** !
function StructGridComputeBoundConnect(structured_grid, xc, yc, zc, option)

  use Connection_module
  use Option_module
  
  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
#include "finclude/petscsysdef.h"
       PetscFortranAddr :: p_patch
       PetscInt :: axis,dim
     end function samr_patch_at_bc
  end interface

  PetscReal, pointer :: xc(:), yc(:), zc(:)
  type(connection_set_type), pointer :: StructGridComputeBoundConnect
  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  PetscInt :: i, j, k, iconn, id_up, id_dn
  PetscInt :: samr_ofx, samr_ofy, samr_ofz
  PetscInt :: nconn
  PetscInt :: lenx, leny, lenz
  PetscReal :: dist_up, dist_dn
  type(connection_set_type), pointer :: connections
  PetscErrorCode :: ierr
  PetscReal, pointer :: radius(:)

  radius => xc
  
  samr_ofx = 0
  samr_ofy = 0
  samr_ofz = 0
  ! the adjustments in the case of AMR are based on the PIMS code adjustments by LC

  nconn = structured_grid%nly*structured_grid%nlz*(structured_grid%nlx + 2 - structured_grid%ngx)&
         +structured_grid%nlx*structured_grid%nlz*(structured_grid%nly + 2 - structured_grid%ngy)&  
         +structured_grid%nlx*structured_grid%nly*(structured_grid%nlz + 2 - structured_grid%ngz)

  lenx = structured_grid%nlx 
  leny = structured_grid%nly 
  lenz = structured_grid%nlz 


!  write(*,*) option%myrank, 'nxs=',structured_grid%nxs,' nys=',structured_grid%nys, ' nzs=',structured_grid%nzs 
!  write(*,*) option%myrank, 'nxe=',structured_grid%nxe,' nye=',structured_grid%nye, ' nze=',structured_grid%nze
!  write(*,*) option%myrank, 'nx=',structured_grid%nlx,' ny=',structured_grid%ny, ' nz=',structured_grid%nz  
!  write(*,*) option%myrank, 'istart=',structured_grid%istart,' iend=',structured_grid%iend
!  write(*,*) option%myrank, 'jstart=',structured_grid%jstart,' jend=',structured_grid%jend
!  write(*,*) option%myrank, 'kstart=',structured_grid%kstart,' kend=',structured_grid%kend
!  write(*,*) option%myrank, 'boundary nconn=',nconn

!  stop


  connections => ConnectionCreate(nconn, &
                                  option%nphase,BOUNDARY_CONNECTION_TYPE)

!  StructGridComputeBoundConnect => connections

  iconn = 0
  ! x-connections
  if (structured_grid%nlx + 2 - structured_grid%ngx > 0) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        if (structured_grid%nxs == structured_grid%ngxs) then
          i = structured_grid%istart 
          do k = structured_grid%kstart, structured_grid%kend
            do j = structured_grid%jstart, structured_grid%jend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy+samr_ofx

              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = 0.
              connections%dist(0,iconn) = dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = structured_grid%dy(id_dn)* &
                                        structured_grid%dz(id_dn)
              connections%cntr(1,iconn) = xc(id_dn) - dist_dn
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) 
            enddo
          enddo
        endif
        if (structured_grid%nxe == structured_grid%ngxe) then
          i = structured_grid%iend 
          do k = structured_grid%kstart, structured_grid%kend
            do j = structured_grid%jstart, structured_grid%jend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy+samr_ofx
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(1,iconn) = -1.d0  ! x component of unit vector
              connections%area(iconn) = structured_grid%dy(id_dn)* &
                                        structured_grid%dz(id_dn)
              connections%cntr(1,iconn) = xc(id_dn) + dist_dn
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) 
            enddo
          enddo
        endif
      case(CYLINDRICAL_GRID)
        print *, 'Boundary connection - Cylindrical coordinates are not implemented.'
        stop
      case(SPHERICAL_GRID)
        print *, 'Boundary connection -  Spherical coordinates are not implemented.'
        stop
  end select
    
  endif


  ! y-connections
  if (structured_grid%nly + 2  -  structured_grid%ngy > 0) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        if (structured_grid%nys == structured_grid%ngys) then 
          j = structured_grid%jstart
          do k = structured_grid%kstart, structured_grid%kend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy &
                  +samr_ofy
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dy(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(2,iconn) = 1.d0  ! y component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn)* &
                                    structured_grid%dz(id_dn)
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) - dist_dn
              connections%cntr(3,iconn) = zc(id_dn) 
            enddo
          enddo
        endif
        if (structured_grid%nye == structured_grid%ngye) then 
          j = structured_grid%jend
          do k = structured_grid%kstart, structured_grid%kend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy &
                  +samr_ofy
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dy(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(2,iconn) = -1.d0  ! y component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn)* &
                                    structured_grid%dz(id_dn)
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) + dist_dn
              connections%cntr(3,iconn) = zc(id_dn) 
            enddo
          enddo
        endif
      case(CYLINDRICAL_GRID)
        print *, 'Boundary connection - Cylindrical coordinates are not implemented.'
        stop
      case(SPHERICAL_GRID)
        print *, 'Boundary connection -  Spherical coordinates are not implemented.'
        stop
    end select
  endif

      
  ! z-connections
  if (structured_grid%nlz + 2 - structured_grid%ngz > 0) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        if (structured_grid%nzs == structured_grid%ngzs) then
          k = structured_grid%kstart
          do j = structured_grid%jstart, structured_grid%jend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * &
                  structured_grid%ngxy + samr_ofz
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(3,iconn) = 1.d0  ! z component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn) * &
                                        structured_grid%dy(id_dn)
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) - dist_dn 
            enddo
          enddo
        endif    
        if (structured_grid%nze == structured_grid%ngze ) then
          k = structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * &
                  structured_grid%ngxy + samr_ofz
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(3,iconn) = -1.d0  ! z component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn) * &
                                        structured_grid%dy(id_dn)
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) + dist_dn 
            enddo
          enddo
        endif    
      case(CYLINDRICAL_GRID)
        print *, 'Boundary connection - Cylindrical coordinates are not implemented.'
        stop
      case(SPHERICAL_GRID)
        print *, 'Boundary connection -  Spherical coordinates are not implemented.'
        stop
  end select
  endif
 
 
  StructGridComputeBoundConnect => connections
!    stop

end function StructGridComputeBoundConnect





! ************************************************************************** !
!
! StructGridPopulateConnection: Computes details of connection (area, dist, etc)
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine StructGridPopulateConnection(radius,structured_grid,connection,iface, &
                                        iconn,ghosted_id)

  use Connection_module
  
  implicit none
 
  type(structured_grid_type) :: structured_grid
  type(connection_set_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: ghosted_id
  PetscReal :: radius(:)
  
  PetscErrorCode :: ierr
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
      select case(iface)

        case(WEST_FACE,EAST_FACE)

          select case(structured_grid%itype)
            case(CARTESIAN_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dx(ghosted_id)
              connection%area(iconn) = structured_grid%dy(ghosted_id)* &
                                   structured_grid%dz(ghosted_id)


              if (iface ==  WEST_FACE) then
                connection%dist(1,iconn) = 1.d0
              else
                connection%dist(1,iconn) = -1.d0
              endif
            case(CYLINDRICAL_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dx(ghosted_id)
              if (iface ==  WEST_FACE) then
                connection%dist(1,iconn) = 1.d0
                connection%area(iconn) = 2.d0 * pi * (radius(ghosted_id)-0.5d0*structured_grid%dx(ghosted_id)) * &
                                        structured_grid%dz(ghosted_id)
              else
                connection%dist(1,iconn) = -1.d0
                connection%area(iconn) = 2.d0 * pi * (radius(ghosted_id)+0.5d0*structured_grid%dx(ghosted_id)) * &
                                        structured_grid%dz(ghosted_id)
              endif
            case(SPHERICAL_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dx(ghosted_id)
              if (iface ==  WEST_FACE) then
                connection%dist(1,iconn) = 1.d0
                connection%area(iconn) = 0.d0
              else
                connection%dist(1,iconn) = -1.d0
                connection%area(iconn) = 4.d0 * pi * (radius(ghosted_id)+0.5d0*structured_grid%dx(ghosted_id))
              endif
          end select

        case(SOUTH_FACE,NORTH_FACE)

          select case(structured_grid%itype)
            case(CARTESIAN_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dy(ghosted_id)
              connection%area(iconn) = structured_grid%dx(ghosted_id)* &
                                   structured_grid%dz(ghosted_id)
              if (iface ==  SOUTH_FACE) then
                connection%dist(2,iconn) = 1.d0
              else
                connection%dist(2,iconn) = -1.d0
              endif
            case(CYLINDRICAL_GRID)
              print *, 'Cylindrical coordinates not applicable.'
              stop
            case(SPHERICAL_GRID)
              print *, 'Spherical coordinates not applicable.'
              stop
          end select

        case(BOTTOM_FACE,TOP_FACE)

          select case(structured_grid%itype)
            case(CARTESIAN_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dz(ghosted_id)
              connection%area(iconn) = structured_grid%dx(ghosted_id)* &
                                   structured_grid%dy(ghosted_id)
              if (structured_grid%invert_z_axis) then
                if (iface ==  TOP_FACE) then 
                  connection%dist(3,iconn) = 1.d0
                else
                  connection%dist(3,iconn) = -1.d0
                endif
              else
                if (iface ==  TOP_FACE) then 
                  connection%dist(3,iconn) = -1.d0
                else
                  connection%dist(3,iconn) = 1.d0
                endif
              endif
            case(CYLINDRICAL_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dz(ghosted_id)
              connection%area(iconn) = 2.d0 * pi * radius(ghosted_id) * &
                                        structured_grid%dx(ghosted_id)
              if (structured_grid%invert_z_axis) then
                if (iface ==  TOP_FACE) then 
                  connection%dist(3,iconn) = 1.d0
                else
                  connection%dist(3,iconn) = -1.d0
                endif
              else
                if (iface ==  TOP_FACE) then 
                  connection%dist(3,iconn) = -1.d0
                else
                  connection%dist(3,iconn) = 1.d0
                endif
              endif
            case(SPHERICAL_GRID)
              print *, 'Areas for spherical coordinates for z-axis not applicable.'
              stop
          end select
      end select
    case(INITIAL_CONNECTION_TYPE)
    case(SRC_SINK_CONNECTION_TYPE)
  end select
  
end subroutine StructGridPopulateConnection

! ************************************************************************** !
!
! StructuredGridComputeVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine StructuredGridComputeVolumes(radius,structured_grid,option,nL2G,volume)

  use Option_module
  
  implicit none

! These includes are needed for VecRestoreArrayF90() - geh
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: nL2G(:)
  PetscReal :: radius(:)
  Vec :: volume
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  PetscInt :: local_id, ghosted_id
  PetscReal, pointer :: volume_p(:)
  PetscReal :: r_up, r_down
  PetscErrorCode :: ierr
  
  call StructuredGridVecGetArrayF90(structured_grid, volume,volume_p, ierr)

  select case(structured_grid%itype)
    case(CARTESIAN_GRID)
      do local_id=1, structured_grid%nlmax
        ghosted_id = nL2G(local_id)
        volume_p(local_id) = structured_grid%dx(ghosted_id) * &
                             structured_grid%dy(ghosted_id) * &
                             structured_grid%dz(ghosted_id)
      enddo
    case(CYLINDRICAL_GRID)
      do local_id=1, structured_grid%nlmax
        ghosted_id = nL2G(local_id)
        volume_p(local_id) = 2.d0 * pi * radius(ghosted_id) * structured_grid%dx(ghosted_id) * &
                                structured_grid%dz(ghosted_id)
      enddo
    case(SPHERICAL_GRID)
      do local_id=1, structured_grid%nlmax
        ghosted_id = nL2G(local_id)
        r_up = radius(ghosted_id) + 0.5d0*structured_grid%dx(ghosted_id)
        r_down = radius(ghosted_id) - 0.5d0*structured_grid%dx(ghosted_id)
        volume_p(local_id) = 4.d0/3.d0 * pi * structured_grid%dx(ghosted_id) &
        * (r_up*r_up + r_up*r_down + r_down*r_down)
      enddo
  end select
  
  call StructGridVecRestoreArrayF90(structured_grid,volume,volume_p, ierr)
  
  if (OptionPrintToScreen(option) .and. &
      option%mycommsize > 1 .and. option%mycommsize <= 16) then
    write(*,'(" rank= ",i3,", nlmax= ",i6,", nlx,y,z= ",3i4, &
      & ", nxs,e = ",2i4,", nys,e = ",2i4,", nzs,e = ",2i4)') &
      option%myrank,structured_grid%nlmax,structured_grid%nlx, &
        structured_grid%nly,structured_grid%nlz,structured_grid%nxs, &
        structured_grid%nxe,structured_grid%nys,structured_grid%nye, &
        structured_grid%nzs,structured_grid%nze

    write(*,'(" rank= ",i3,", ngmax= ",i6,", ngx,y,z= ",3i4, &
      & ", ngxs,e= ",2i4,", ngys,e= ",2i4,", ngzs,e= ",2i4)') &
      option%myrank,structured_grid%ngmax,structured_grid%ngx, &
        structured_grid%ngy,structured_grid%ngz,structured_grid%ngxs, &
        structured_grid%ngxe,structured_grid%ngys,structured_grid%ngye, &
        structured_grid%ngzs,structured_grid%ngze
  endif

end subroutine StructuredGridComputeVolumes

! ************************************************************************** !
!
! StructuredGridMapIndices: maps global, local and natural indices of cells 
!                          to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridMapIndices(structured_grid,nG2L,nL2G,nL2A,nG2A)

  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
#include "finclude/petscsysdef.h"
       PetscFortranAddr :: p_patch
       PetscInt :: axis,dim
     end function samr_patch_at_bc
  end interface

  type(structured_grid_type) :: structured_grid
  PetscInt, pointer :: nG2L(:), nL2G(:), nL2A(:), nG2A(:)

  PetscInt :: i, j, k, local_id, ghosted_id, natural_id, count1
  PetscErrorCode :: ierr
  
  allocate(nL2G(structured_grid%nlmax))
  allocate(nG2L(structured_grid%ngmax))
! only allocate space for the next two arrays if the current grid is not
! not part of an AMR grid hierarchy
  if(structured_grid%p_samr_patch == 0) then
     allocate(nL2A(structured_grid%nlmax))
     allocate(nG2A(structured_grid%ngmax))
  endif

  structured_grid%istart = structured_grid%nxs-structured_grid%ngxs
  structured_grid%jstart = structured_grid%nys-structured_grid%ngys
  structured_grid%kstart = structured_grid%nzs-structured_grid%ngzs
  structured_grid%iend = structured_grid%istart+structured_grid%nlx-1
  structured_grid%jend = structured_grid%jstart+structured_grid%nly-1
  structured_grid%kend = structured_grid%kstart+structured_grid%nlz-1

  ! Local <-> Ghosted Transformation
  nG2L = 0  ! Must initialize this to zero!
  nL2G = 0

  if(structured_grid%p_samr_patch == 0)then
     nG2A = 0
     nL2A = 0
  endif

  local_id = 0
  do k=structured_grid%kstart,structured_grid%kend
    do j=structured_grid%jstart,structured_grid%jend
      do i=structured_grid%istart,structured_grid%iend
        local_id = local_id + 1
        ghosted_id = i+j*structured_grid%ngx+k*structured_grid%ngxy+1
        nL2G(local_id) = ghosted_id
        nG2L(ghosted_id) = local_id
      enddo
    enddo
  enddo

  do i=1,structured_grid%ngmax
    j = nG2L(i)
    if (j > 0) then
      k = nL2G(j)
      if (i /= k) then
        print *,'Error in ghost-local node numbering for ghost node =', i
        print *,'node_id_gtol(i) =', j
        print *,'node_id_ltog(node_id_gtol(i)) =', k
        stop
      endif
    endif
  enddo
  ! Local(non ghosted)->Natural(natural order starts from 0)

  !geh - set corner ghosted nodes to -1
  do k=1,structured_grid%ngz
    do j=1,structured_grid%ngy
      do i=1,structured_grid%ngx
        count1 = 0
        if (i == 1 .and. &
            abs(structured_grid%nxs-structured_grid%ngxs) > 0) &
          count1 = count1 + 1
        if (i == structured_grid%ngx .and. &
            abs(structured_grid%ngxe-structured_grid%nxe) > 0) &
          count1 = count1 + 1
        if (j == 1 .and. &
            abs(structured_grid%nys-structured_grid%ngys) > 0) &
          count1 = count1 + 1
        if (j == structured_grid%ngy .and. &
            abs(structured_grid%ngye-structured_grid%nye) > 0) &
          count1 = count1 + 1
        if (k == 1 .and. &
            abs(structured_grid%nzs-structured_grid%ngzs) > 0) &
          count1 = count1 + 1
        if (k == structured_grid%ngz .and. &
            abs(structured_grid%ngze-structured_grid%nze) > 0) &
          count1 = count1 + 1
        if (count1 > 1) then
          ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
          nG2L(ghosted_id) = -1
        endif
      enddo
    enddo
  enddo

  if(.not.(structured_grid%p_samr_patch == 0)) then
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ZERO_INTEGER, ZERO_INTEGER) ==1) then
        do k=1,structured_grid%ngz
           do j=1,structured_grid%ngy
              ghosted_id = 1+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ZERO_INTEGER, ONE_INTEGER) ==1) then 
        i=structured_grid%ngx
        do k=1,structured_grid%ngz
           do j=1,structured_grid%ngy
              ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ONE_INTEGER, ZERO_INTEGER) ==1) then
        do k=1,structured_grid%ngz
           do i=1,structured_grid%ngx
              ghosted_id = i+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, ONE_INTEGER, ONE_INTEGER) ==1) then
        j=structured_grid%ngy
        do k=1,structured_grid%ngz
           do i=1,structured_grid%ngx
              ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, TWO_INTEGER, ZERO_INTEGER) ==1) then
        do j=1,structured_grid%ngy
           do i=1,structured_grid%ngx
              ghosted_id = i+(j-1)*structured_grid%ngx
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, TWO_INTEGER, ONE_INTEGER) ==1) then
        k=structured_grid%ngz
        do j=1,structured_grid%ngy
           do i=1,structured_grid%ngx
              ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
  endif

  if(structured_grid%p_samr_patch == 0) then
     local_id=0
     do k=1,structured_grid%nlz
        do j=1,structured_grid%nly
           do i=1,structured_grid%nlx
              local_id = local_id + 1
              natural_id = i-1+structured_grid%nxs+(j-1+structured_grid%nys)*structured_grid%nx+ &
                   (k-1+structured_grid%nzs)*structured_grid%nxy
              if (natural_id>(structured_grid%nmax-1)) print *,'Wrong Nature order....'
              nL2A(local_id) = natural_id
           enddo
        enddo
     enddo
     ! Local(ghosted)->Natural(natural order starts from 0)
     local_id=0
     do k=1,structured_grid%ngz
        do j=1,structured_grid%ngy
           do i=1,structured_grid%ngx
              local_id = local_id + 1
              natural_id = i-1+structured_grid%ngxs+(j-1+structured_grid%ngys)*structured_grid%nx+ &
                   (k-1+structured_grid%ngzs)*structured_grid%nxy
              if (natural_id>(structured_grid%nmax-1)) print *,'Wrong Nature order....'
              nG2A(local_id) = natural_id
           enddo
        enddo
     enddo
  endif

end subroutine StructuredGridMapIndices

! ************************************************************************** !
!
! StructGridGetGhostedNeighbors: Returns an array of neighboring cells
! author: Glenn Hammond
! date: 01/28/11
!
! ************************************************************************** !
subroutine StructGridGetGhostedNeighbors(structured_grid,ghosted_id, &
                                         stencil_type, &
                                         stencil_width_i,stencil_width_j, &
                                         stencil_width_k,x_count,y_count, &
                                         z_count,ghosted_neighbors, &
                                         option)

  use Option_module

  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscInt :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: x_count
  PetscInt :: y_count
  PetscInt :: z_count
  PetscInt :: ghosted_neighbors(*)

  PetscInt :: i, j, k
  PetscInt :: icount
  PetscInt :: ii, jj, kk

  call StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)
  
  x_count = 0
  y_count = 0
  z_count = 0
  icount = 0
  select case(stencil_type)
    case(STAR_STENCIL)
      do ii = max(i-stencil_width_i,1), min(i+stencil_width_i,structured_grid%ngx)
        if (ii /= i) then
          icount = icount + 1
          x_count = x_count + 1
          ghosted_neighbors(icount) = &
            StructGridGetGhostedIDFromIJK(structured_grid,ii,j,k)
        endif
      enddo
      do jj = max(j-stencil_width_j,1), min(j+stencil_width_j,structured_grid%ngy)
        if (jj /= j) then
          icount = icount + 1
          y_count = y_count + 1
          ghosted_neighbors(icount) = &
            StructGridGetGhostedIDFromIJK(structured_grid,i,jj,k)
        endif
      enddo
      do kk = max(k-stencil_width_k,1), min(k+stencil_width_k,structured_grid%ngz)
        if (kk /= k) then
          icount = icount + 1
          z_count = z_count + 1
          ghosted_neighbors(icount) = &
            StructGridGetGhostedIDFromIJK(structured_grid,i,j,kk)
        endif
      enddo
    case(BOX_STENCIL)
      option%io_buffer = 'BOX_STENCIL not yet supported in ' // &
        'StructGridGetNeighbors.'
      call printErrMsg(option)
  end select

end subroutine StructGridGetGhostedNeighbors

! ************************************************************************** !
!
! StructuredGridDestroy: Deallocates a structured grid
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine StructuredGridDestroy(structured_grid)

  implicit none
  
  type(structured_grid_type), pointer :: structured_grid
  
  PetscErrorCode :: ierr
    
  if (.not.associated(structured_grid)) return
  

  if (associated(structured_grid%dx_global)) deallocate(structured_grid%dx_global)
  nullify(structured_grid%dx_global)
  if (associated(structured_grid%dy_global)) deallocate(structured_grid%dy_global)
  nullify(structured_grid%dy_global)
  if (associated(structured_grid%dz_global)) deallocate(structured_grid%dz_global)
  nullify(structured_grid%dz_global)

  if (associated(structured_grid%dxg_local)) deallocate(structured_grid%dxg_local)
  nullify(structured_grid%dxg_local)
  if (associated(structured_grid%dyg_local)) deallocate(structured_grid%dyg_local)
  nullify(structured_grid%dyg_local)
  if (associated(structured_grid%dzg_local)) deallocate(structured_grid%dzg_local)
  nullify(structured_grid%dzg_local)

  if (associated(structured_grid%dx)) deallocate(structured_grid%dx)
  nullify(structured_grid%dx)
  if (associated(structured_grid%dy)) deallocate(structured_grid%dy)
  nullify(structured_grid%dy)
  if (associated(structured_grid%dz)) deallocate(structured_grid%dz)
  nullify(structured_grid%dz)
  
  deallocate(structured_grid)
  nullify(structured_grid)

end subroutine StructuredGridDestroy

                          
! ************************************************************************** !
!
! StructuredGridVecGetArrayCellF90: Interface for SAMRAI AMR
! author: Bobby Philip
! date: 12/15/10
!
! ************************************************************************** !
subroutine StructuredGridVecGetMaskArrayCellF90(structured_grid, vec, f90ptr, ierr)

 use cf90interface_module

 implicit none 

interface
subroutine samrvecgetmaskarraycellf90(patch, petscvec, f90wrap)
      implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
      PetscFortranAddr, intent(inout):: patch
      Vec:: petscvec
      PetscFortranAddr :: f90wrap
    end subroutine samrvecgetmaskarraycellf90
 end interface

#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

 type(structured_grid_type) :: structured_grid
 Vec:: vec
 PetscReal, pointer :: f90ptr(:)
 PetscErrorCode :: ierr
 
 type(f90ptrwrap), pointer :: ptr
 PetscFortranAddr :: cptr
 
 if(structured_grid%p_samr_patch == 0) then
! we'll have to throw an error here      
 else
    ierr=0
    allocate(ptr)
    nullify(ptr%f90ptr)
    call assign_c_array_ptr(cptr, ptr)
    call samrvecgetmaskarraycellf90(structured_grid%p_samr_patch, vec, cptr)
    f90ptr => ptr%f90ptr
    deallocate(ptr)
 endif
 
end subroutine StructuredGridVecGetMaskArrayCellF90
      
! ************************************************************************** !
!
! StructuredGridVecGetArrayCellF90: Interface for SAMRAI AMR
! author: Bobby Philip
! date: 06/09/08
!
! ************************************************************************** !
subroutine StructuredGridVecGetArrayCellF90(structured_grid, vec, f90ptr, ierr)

 use cf90interface_module

 implicit none 

interface
subroutine samr_vecgetarraycellf90(patch, petscvec, f90wrap)
      implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
      PetscFortranAddr, intent(inout):: patch
      Vec:: petscvec
      PetscFortranAddr :: f90wrap
    end subroutine samr_vecgetarraycellf90
 end interface

#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

 type(structured_grid_type) :: structured_grid
 Vec:: vec
 PetscReal, pointer :: f90ptr(:)
 PetscErrorCode :: ierr
 
 type(f90ptrwrap), pointer :: ptr
 PetscFortranAddr :: cptr
 
 if(structured_grid%p_samr_patch == 0) then
    call VecGetArrayF90(vec, f90ptr, ierr)
 else
    ierr=0
    allocate(ptr)
    nullify(ptr%f90ptr)
    call assign_c_array_ptr(cptr, ptr)
    call samr_vecgetarraycellf90(structured_grid%p_samr_patch, vec, cptr)
    f90ptr => ptr%f90ptr
    deallocate(ptr)
 endif
 
end subroutine StructuredGridVecGetArrayCellF90
                          
! ************************************************************************** !
!
! StructuredGridVecGetArray2F90: Interface for SAMRAI AMR
! author: Bobby Philip
! date: 06/09/08
!
! ************************************************************************** !
subroutine StructuredGridVecGetArraySideF90(structured_grid, axis, vec, f90ptr, ierr)

 use cf90interface_module

 implicit none 

interface
subroutine samr_vecgetarraysidef90(patch, axis, petscvec, f90wrap)
      implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
      PetscFortranAddr, intent(inout):: patch
      PetscInt :: axis
      Vec:: petscvec
      PetscFortranAddr :: f90wrap
    end subroutine samr_vecgetarraysidef90
 end interface

#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

 type(structured_grid_type) :: structured_grid
 PetscInt :: axis
 Vec:: vec
 PetscReal, pointer :: f90ptr(:)
 PetscErrorCode :: ierr
 
 type(f90ptrwrap), pointer :: ptr
 PetscFortranAddr :: cptr
 
 if(structured_grid%p_samr_patch == 0) then
    call VecGetArrayF90(vec, f90ptr, ierr)
 else
    ierr=0
    allocate(ptr)
    nullify(ptr%f90ptr)
    call assign_c_array_ptr(cptr, ptr)
    call samr_vecgetarraysidef90(structured_grid%p_samr_patch, axis, vec, cptr)
    f90ptr => ptr%f90ptr
    deallocate(ptr)
 endif
 
end subroutine StructuredGridVecGetArraySideF90

! ************************************************************************** !
!
! StructuredGridVecGetArrayF90: Interface for SAMRAI AMR
! author: Bobby Philip
! date: 06/09/08
!
! ************************************************************************** !
subroutine StructGridVecRestoreArrayF90(structured_grid, vec, f90ptr, ierr)

 use cf90interface_module

 implicit none 

#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

 type(structured_grid_type) :: structured_grid
 Vec:: vec
 PetscReal, pointer :: f90ptr(:)
 PetscErrorCode :: ierr
 
 type(f90ptrwrap), pointer :: ptr
 PetscFortranAddr :: cptr
 
 if(structured_grid%p_samr_patch == 0) then
    call VecRestoreArrayF90(vec, f90ptr, ierr)
 else
    ierr = 0
 endif
 
end subroutine StructGridVecRestoreArrayF90

end module Structured_Grid_module
