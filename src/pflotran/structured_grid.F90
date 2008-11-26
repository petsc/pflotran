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

    PetscReal :: origin(3) ! local origin of non-ghosted grid
    PetscReal :: bounds(3,3)

    ! grid spacing for each direction for global domain
    PetscReal, pointer :: dx_global(:), dy_global(:), dz_global(:)
    ! grid spacing for each direction for local, ghosted domain
    PetscReal, pointer :: dxg_local(:), dyg_local(:), dzg_local(:)
    
    PetscTruth :: invert_z_axis
    
    PetscReal, pointer :: dx(:), dy(:), dz(:)  ! ghosted grid spacings for each grid cell
    
    PetscFortranAddr :: p_samr_patch ! pointer to a SAMRAI patch object

  end type structured_grid_type

  public :: StructuredGridCreate, &
            StructuredGridDestroy, &
            StructuredGridCreateDA, &
            StructGridComputeLocalBounds, &
            StructGridComputeInternConnect, &
            StructuredGridCreateVecFromDA, &
            StructuredGridMapIndices, &
            StructuredGridComputeSpacing, &
            StructuredGridComputeCoord, &
            StructuredGridReadDXYZ, &
            StructuredGridComputeVolumes, &
            StructGridPopulateConnection, &
            StructGridGetIJKFromCoordinate, &
            StructGridGetIJKFromLocalID, &
            StructGridGetIJKFromGhostedID, &
            StructuredGridVecGetArrayF90, &
            StructGridVecRestoreArrayF90
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
subroutine StructuredGridCreateDA(structured_grid,da,ndof,stencil_width, &
                                  option)

  use Option_module
        
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"

  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  DA :: da
  PetscInt :: ndof
  PetscInt :: stencil_width

  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DA object that will manage communication.
  !-----------------------------------------------------------------------
  call DACreate3D(option%comm,DA_NONPERIODIC,DA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  da,ierr)

end subroutine StructuredGridCreateDA

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
       
#include "finclude/petsc.h"

       PetscFortranAddr :: p_patch
       PetscInt :: nxs, nys, nzs, nlx, nly, nlz

     end subroutine samr_patch_get_corners
     
     subroutine samr_patch_get_ghostcorners(p_patch, nxs, nys, nzs, nlx, nly, nlz)
       implicit none
       
#include "finclude/petsc.h"
       
       PetscFortranAddr :: p_patch
       PetscInt :: nxs, nys, nzs, nlx, nly, nlz

     end subroutine samr_patch_get_ghostcorners
  end interface
     
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"

  type(structured_grid_type) :: structured_grid
  DA :: da

  PetscErrorCode :: ierr

  if(structured_grid%p_samr_patch==0) then
      ! get corner information
     call DAGetCorners(da, structured_grid%nxs, &
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
     call DAGetGhostCorners(da, structured_grid%ngxs, &
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
! StructuredGridCreateVecFromDA: Creates a PETSc vector from a DA
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine StructuredGridCreateVecFromDA(da,vector,vector_type)

  implicit none

  DA :: da
  Vec :: vector
  PetscInt :: vector_type
  
  PetscErrorCode :: ierr

  select case (vector_type)
    case(GLOBAL)
      call DACreateGlobalVector(da,vector,ierr)
    case(LOCAL)
      call DACreateLocalVector(da,vector,ierr)
    case(NATURAL)
      call DACreateNaturalVector(da,vector,ierr)
  end select

end subroutine StructuredGridCreateVecFromDA

! ************************************************************************** !
!
! StructuredGridReadDXYZ: Reads structured grid spacing from input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine StructuredGridReadDXYZ(structured_grid,fid,option)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  PetscInt :: fid
  type(option_type) :: option
  
  PetscInt :: i

  allocate(structured_grid%dx_global(structured_grid%nx))
  structured_grid%dx_global = 0.d0
  allocate(structured_grid%dy_global(structured_grid%ny))
  structured_grid%dy_global = 0.d0
  allocate(structured_grid%dz_global(structured_grid%nz))
  structured_grid%dz_global = 0.d0

  call StructuredGridReadArray(structured_grid%dx_global, &
                               structured_grid%nx,fid,option)
  call StructuredGridReadArray(structured_grid%dy_global, &
                               structured_grid%ny,fid,option)
  call StructuredGridReadArray(structured_grid%dz_global, &
                               structured_grid%nz,fid,option)
    
  if (option%myrank==0) then
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
subroutine StructuredGridReadArray(a,n,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: fid
  PetscInt :: n
  PetscInt :: i, i1, i2, m
  PetscInt ::  nvalue=10
  PetscReal, intent(inout) :: a(*)
  character(len=MAXSTRINGLENGTH) :: string 
  PetscErrorCode :: ierr

  save nvalue

!  call fiReadStringErrorMsg('DXYZ',ierr)

!  call fiReadDouble(string,grid%radius_0,ierr)
!  call fiDefaultMsg(option%myrank,'radius_0',ierr)

  i2 = 0
  do
    i1 = i2+1
    i2 = i2+nvalue
    if (i2.gt.n) i2 = n
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(option%myrank,'DXYZ',ierr)
    do i = i1, i2
      call fiReadDouble(string, a(i), ierr)
!geh  ierr, which comes from iostat, will not necessarily be 0 and 1, 
!geh  it could be another number
!geh      if (ierr .eq. 1) a(i) = 0.d0
      if (ierr /= 0) a(i) = 0.d0
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
    if (i2.ge.n) exit
  enddo
    
end subroutine StructuredGridReadArray

! ************************************************************************** !
!
! StructuredGridComputeSpacing: Computes structured grid spacing
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine StructuredGridComputeSpacing(structured_grid,nG2A,nG2L)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  PetscInt :: nG2A(:)
  PetscInt :: nG2L(:)
  
  PetscInt :: i, j, k, ghosted_id
  PetscErrorCode :: ierr

  allocate(structured_grid%dxg_local(structured_grid%ngx))
  structured_grid%dxg_local = 0.d0
  allocate(structured_grid%dyg_local(structured_grid%ngy))
  structured_grid%dyg_local = 0.d0
  allocate(structured_grid%dzg_local(structured_grid%ngz))
  structured_grid%dzg_local = 0.d0
  
  if (structured_grid%p_samr_patch .eq. 0) then
    if (.not.associated(structured_grid%dx_global)) then
      ! indicates that the grid spacings still need to be computed
      if (structured_grid%bounds(1,1) < -1.d19) then ! bounds have not been initialized
        print *, 'ERROR: Bounds have not been set for grid and DXYZ does not exist'
        stop
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
!         print *, 'CYLINDRICAL grids still need dr to be set up'
!         stop
          structured_grid%dx_global = (structured_grid%bounds(X_DIRECTION,UPPER)- &
                                       structured_grid%bounds(X_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nx)
          structured_grid%dy = 1.d0
          structured_grid%dz_global = (structured_grid%bounds(Z_DIRECTION,UPPER)- &
                                       structured_grid%bounds(Z_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nz)
        case(SPHERICAL_GRID)
!         print *, 'CYLINDRICAL grids still need dr to be set up'
!         stop
          structured_grid%dx_global = (structured_grid%bounds(X_DIRECTION,UPPER)- &
                                       structured_grid%bounds(X_DIRECTION,LOWER)) / &
                                       dble(structured_grid%nx)
          structured_grid%dy = 1.d0
          structured_grid%dz = 1.d0
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
#include "finclude/petsc.h"
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

! PetscInt :: ierr
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
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: i, j, k
  PetscInt :: i_local, j_local, k_local
  PetscInt :: i_ghosted, j_ghosted, k_ghosted
  PetscReal :: x, y, z
  
  PetscReal :: x_lower_face, y_lower_face, z_lower_face

  i = -1
  j = -1
  k = -1

  x_lower_face = structured_grid%origin(X_DIRECTION)
  i_local = 1
  do i_ghosted=structured_grid%istart,structured_grid%iend
    if (x >= x_lower_face .and. &                   ! since i_ghosted is zero-based
        x <= x_lower_face+structured_grid%dxg_local(i_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in x-dir on proc
      if (i_ghosted == structured_grid%istart) then
        ! located on upwind boundary and ghosted
        if (x == x_lower_face .and. &
            structured_grid%nxs /= structured_grid%ngxs) exit
      endif
      i = i_local 
      exit
    endif
    i_local = i_local + 1
    x_lower_face = x_lower_face + structured_grid%dxg_local(i_ghosted+1)
  enddo
  y_lower_face = structured_grid%origin(Y_DIRECTION)
  j_local = 1
  do j_ghosted=structured_grid%jstart,structured_grid%jend
    if (y >= y_lower_face .and. &
        y <= y_lower_face+structured_grid%dyg_local(j_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in y-dir on proc
      if (j_ghosted == structured_grid%jstart) then
        ! located on upwind boundary and ghosted
        if (y == y_lower_face .and. &
            structured_grid%nys /= structured_grid%ngys) exit
      endif      
      j = j_local
      exit
    endif
    j_local = j_local + 1
    y_lower_face = y_lower_face + structured_grid%dyg_local(j_ghosted+1)
  enddo
  z_lower_face = structured_grid%origin(Z_DIRECTION)
  k_local = 1
  do k_ghosted=structured_grid%kstart,structured_grid%kend
    if (z >= z_lower_face .and. &
        z <= z_lower_face+structured_grid%dzg_local(k_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in z-dir on proc
      if (k_ghosted == structured_grid%kstart) then
        ! if located on upwind boundary and ghosted, skip
        if (z == z_lower_face .and. &
            structured_grid%nzs /= structured_grid%ngzs) exit
      endif          
      k = k_local
      exit
    endif
    k_local = k_local + 1
    z_lower_face = z_lower_face + structured_grid%dzg_local(k_ghosted+1)
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
! StructGridComputeInternConnect: computes internal connectivity of a  
!                               structured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function StructGridComputeInternConnect(radius,structured_grid,option)

  use Connection_module
  use Option_module
  
  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
#include "finclude/petsc.h"
       PetscFortranAddr :: p_patch
       PetscInt :: axis,dim
     end function samr_patch_at_bc
  end interface

  PetscReal :: radius(:)
  type(connection_set_type), pointer :: StructGridComputeInternConnect
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
  
  call ConnectionAllocateLists()
  
  samr_ofx = 0
  samr_ofy = 0
  samr_ofz = 0
  ! the adjustments in the case of AMR are based on the PIMS code adjustments by LC
  nconn = (structured_grid%ngx-1)*structured_grid%nly*structured_grid%nlz+ &
          structured_grid%nlx*(structured_grid%ngy-1)*structured_grid%nlz+ &
          structured_grid%nlx*structured_grid%nly*(structured_grid%ngz-1)

  lenx = structured_grid%ngx - 1
  leny = structured_grid%ngy - 1
  lenz = structured_grid%ngz - 1

  if(.not.(structured_grid%p_samr_patch.eq.0)) then
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 0, 0) ==1) then
        nconn = nconn - structured_grid%nlyz
        lenx = lenx-1
        samr_ofx = 1
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 0, 1) ==1) then 
        nconn = nconn - structured_grid%nlyz
        lenx = lenx-1
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 1, 0) ==1) then
        nconn = nconn - structured_grid%nlxz
        leny=leny-1
        samr_ofy = structured_grid%ngx
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 1, 1) ==1) then
        nconn = nconn - structured_grid%nlxz
        leny=leny-1
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 2, 0) ==1) then
        nconn = nconn - structured_grid%nlxy
        lenz=lenz-1
        samr_ofz = structured_grid%ngxy
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 2, 1) ==1) then
        nconn = nconn - structured_grid%nlxy
        lenz=lenz-1
     endif  
  endif

  connections => &
       ConnectionCreate(nconn, &
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
              connections%area(iconn) = structured_grid%dy(id_up)* &
                                        structured_grid%dz(id_up)
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
            enddo
          enddo
        enddo
      case(CYLINDRICAL_GRID)
        print *, 'Cylindrical coordinates not applicable.'
        stop
      case(SPHERICAL_GRID)
        print *, 'Spherical coordinates not applicable.'
        stop
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
              connections%area(iconn) = 2.d0 * pi * radius(id_up) * &
                                        structured_grid%dx(id_up)
            enddo
          enddo
        enddo
      case(SPHERICAL_GRID)
        print *, 'Areas for spherical coordinates for z-axis not applicable.'
        stop
  end select
  endif
  
  StructGridComputeInternConnect => connections

end function StructGridComputeInternConnect

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
                connection%area(iconn) = 0.d0
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
  
  if (option%commsize <= 16) then
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
#include "finclude/petsc.h"
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
  if(structured_grid%p_samr_patch.eq.0) then
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

  if(structured_grid%p_samr_patch.eq.0)then
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

  if(.not.(structured_grid%p_samr_patch.eq.0)) then
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 0, 0) ==1) then
        do k=1,structured_grid%ngz
           do j=1,structured_grid%ngy
              ghosted_id = 1+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 0, 1) ==1) then 
        i=structured_grid%ngx
        do k=1,structured_grid%ngz
           do j=1,structured_grid%ngy
              ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 1, 0) ==1) then
        do k=1,structured_grid%ngz
           do i=1,structured_grid%ngx
              ghosted_id = i+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 1, 1) ==1) then
        j=structured_grid%ngy
        do k=1,structured_grid%ngz
           do i=1,structured_grid%ngx
              ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 2, 0) ==1) then
        do j=1,structured_grid%ngy
           do i=1,structured_grid%ngx
              ghosted_id = i+(j-1)*structured_grid%ngx
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
     if(samr_patch_at_bc(structured_grid%p_samr_patch, 2, 1) ==1) then
        k=structured_grid%ngz
        do j=1,structured_grid%ngy
           do i=1,structured_grid%ngx
              ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
              nG2L(ghosted_id) = -1
           enddo
        enddo
     endif  
  endif

  if(structured_grid%p_samr_patch.eq.0) then
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
! StructuredGridVecGetArrayF90: Interface for SAMRAI AMR
! author: Bobby Philip
! date: 06/09/08
!
! ************************************************************************** !
subroutine StructuredGridVecGetArrayF90(structured_grid, vec, f90ptr, ierr)

 use cf90interface_module

 implicit none 

 interface
    subroutine samr_vecgetarrayf90(patch, petscvec, f90wrap)
      implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
      PetscFortranAddr, intent(inout):: patch
      Vec:: petscvec
      PetscFortranAddr :: f90wrap
    end subroutine samr_vecgetarrayf90
 end interface

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

 type(structured_grid_type) :: structured_grid
 Vec:: vec
 PetscReal, pointer :: f90ptr(:)
 PetscInt :: ierr
 
 type(f90ptrwrap), pointer :: ptr
 PetscFortranAddr :: cptr
 
 if(structured_grid%p_samr_patch .eq. 0) then
    call VecGetArrayF90(vec, f90ptr, ierr)
 else
    ierr=0
    allocate(ptr)
    nullify(ptr%f90ptr)
    call assign_c_array_ptr(cptr, ptr)
    call samr_vecgetarrayf90(structured_grid%p_samr_patch, vec, cptr)
    f90ptr => ptr%f90ptr
    deallocate(ptr)
 endif
 
end subroutine StructuredGridVecGetArrayF90

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

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

 type(structured_grid_type) :: structured_grid
 Vec:: vec
 PetscReal, pointer :: f90ptr(:)
 PetscInt :: ierr
 
 type(f90ptrwrap), pointer :: ptr
 PetscFortranAddr :: cptr
 
 if(structured_grid%p_samr_patch .eq. 0) then
    call VecRestoreArrayF90(vec, f90ptr, ierr)
 else
    ierr = 0
 endif
 
end subroutine StructGridVecRestoreArrayF90

end module Structured_Grid_module
