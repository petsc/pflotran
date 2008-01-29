module Structured_Grid_module

  implicit none
 
  private
 
#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"

  PetscInt, parameter :: WEST = 1
  PetscInt, parameter :: EAST = 2
  PetscInt, parameter :: SOUTH = 3
  PetscInt, parameter :: NORTH = 4
  PetscInt, parameter :: BOTTOM = 5
  PetscInt, parameter :: TOP = 6

  type, public :: structured_grid_type

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
      ! istart gives the local x-index of the non-ghosted starting (lower left)
      ! corner. iend gives the local x-index of the non-ghosted ending 
      ! corner. jstart, jend correspond to y-index, kstart, kend to z-index.    

    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.

    PetscReal, pointer :: dx0(:), dy0(:), dz0(:), rd(:)
    
    PetscInt :: igeom
    
    PetscReal :: radius_0
    
    logical :: invert_z_axis
    
    Vec :: dx, dy, dz, dx_loc, dy_loc, dz_loc  ! Grid spacings

    DA :: da_1_dof, da_nphase_dof, da_3np_dof, da_ndof, da_nphancomp_dof, &
          da_nphanspec_dof, da_nphanspecncomp_dof, da_var_dof 
    
  end type

  public :: StructuredGridInit, &
            StructuredGridDestroy, &
            StructuredGridCreateDMs, &
            StructGridComputeInternConnect, &
            StructuredGridCreateVecFromDA, &
            StructuredGridMapIndices, &
            StructuredGridComputeSpacing, &
            StructuredGridComputeCoord, &
            StructuredGridCreateJacobian, &
            StructuredGridCreateColoring, &
            StructureGridGlobalToLocal, &
            StructureGridLocalToGlobal, &
            StructureGridGlobalToNatural, &
            StructureGridNaturalToGlobal, &
            StructureGridLocalToLocalBegin, &
            StructureGridLocalToLocalEnd, &
            StructureGridGlobalToLocalBegin, &
            StructureGridGlobalToLocalEnd, &
            StructureGridGlobalToNaturBegin, &
            StructureGridGlobalToNaturEnd, &
            StructureGridNaturToGlobalBegin, &
            StructureGridNaturToGlobalEnd, &
            StructureGridLocalToLocal, &
            StructuredGridReadDXYZ, &
            StructuredGridComputeVolumes, &
            StructGridPopulateConnection

contains

! ************************************************************************** !
!
! StructuredGridInit: Initializes a structured grid object
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine StructuredGridInit(structured_grid)

  implicit none
  
  type(structured_grid_type) :: structured_grid

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
  
  nullify(structured_grid%dx0)
  nullify(structured_grid%dy0)
  nullify(structured_grid%dz0)
  
  nullify(structured_grid%rd)
  
  structured_grid%igeom = 0
  structured_grid%radius_0 = 0.d0
  
  structured_grid%invert_z_axis = .false.
  
  ! nullify Vec pointers
  structured_grid%dx = 0
  structured_grid%dy = 0
  structured_grid%dz = 0
  structured_grid%dx_loc = 0
  structured_grid%dy_loc = 0
  structured_grid%dz_loc = 0
  
  ! nullify DA pointers
  structured_grid%da_1_dof = 0
  structured_grid%da_nphase_dof = 0
  structured_grid%da_3np_dof = 0
  structured_grid%da_ndof = 0
  structured_grid%da_nphancomp_dof = 0
  structured_grid%da_nphanspec_dof = 0
  structured_grid%da_nphanspecncomp_dof = 0
  structured_grid%da_var_dof = 0
  structured_grid%da_1_dof = 0
  structured_grid%da_1_dof = 0  
  
end subroutine StructuredGridInit
  
! ************************************************************************** !
!
! StructuredGridCreateDMs: Creates structured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine StructuredGridCreateDMs(structured_grid,option)
      
  use Option_module
      
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option

  PetscInt :: ndof
  PetscInt, parameter :: stencil_width = 1
  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DA objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       structured_grid%da_1_dof,ierr)

  ndof = option%nphase
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  structured_grid%da_nphase_dof,ierr)

#if 0
  ndof = 3*option%nphase
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  structured_grid%da_3np_dof,ierr)                
#endif

  ndof = option%ndof
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  option%ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  structured_grid%da_ndof,ierr)

  select case(option%imode) 
    case(TWOPH_MODE)
      ndof = option%nphase*option%npricomp
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                      structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                      structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                      ndof,stencil_width, &
                      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                      structured_grid%da_nphancomp_dof,ierr)
      ndof = option%nphase*option%nspec
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                      structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                      structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                      ndof,stencil_width, &
                      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                      structured_grid%da_nphanspec_dof,ierr)
      ndof = option%nphase*option%nspec*option%npricomp
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                      structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                      structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                      ndof,stencil_width, &
                      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                      structured_grid%da_nphanspecncomp_dof,ierr)
    case(MPH_MODE,VADOSE_MODE,FLASH_MODE,OWG_MODE)
      ndof = (option%ndof+1)*(2+7*option%nphase + 2*option%nspec*option%nphase)
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                      structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                      structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                      ndof,stencil_width, &
                      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                      structured_grid%da_var_dof,ierr)
  end select

 ! get corner information
  call DAGetCorners(structured_grid%da_nphase_dof, structured_grid%nxs, &
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
  call DAGetGhostCorners(structured_grid%da_nphase_dof, structured_grid%ngxs, &
                         structured_grid%ngys, structured_grid%ngzs, structured_grid%ngx, &
                         structured_grid%ngy, structured_grid%ngz, ierr)

  structured_grid%ngxe = structured_grid%ngxs + structured_grid%ngx
  structured_grid%ngye = structured_grid%ngys + structured_grid%ngy
  structured_grid%ngze = structured_grid%ngzs + structured_grid%ngz
  structured_grid%ngxy = structured_grid%ngx * structured_grid%ngy
  structured_grid%ngxz = structured_grid%ngx * structured_grid%ngz
  structured_grid%ngyz = structured_grid%ngy * structured_grid%ngz
  structured_grid%ngmax = structured_grid%ngx * structured_grid%ngy * structured_grid%ngz

end subroutine StructuredGridCreateDMs

! ************************************************************************** !
!
! StructuredGridCreateVecFromDA: Creates a global PETSc vector
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridCreateVecFromDA(structured_grid,da_index,vector, &
                                   vector_type)

  implicit none
  type(structured_grid_type) :: structured_grid
  PetscInt :: da_index
  Vec :: vector
  PetscInt :: vector_type
  
  PetscErrorCode :: ierr

  select case (vector_type)
    case(GLOBAL)
      select case (da_index)
        case(ONEDOF)
          call DACreateGlobalVector(structured_grid%da_1_dof,vector,ierr)
        case(NPHASEDOF)
          call DACreateGlobalVector(structured_grid%da_nphase_dof,vector,ierr)
        case(THREENPDOF)
          call DACreateGlobalVector(structured_grid%da_3np_dof,vector,ierr)
        case(NDOF)
          call DACreateGlobalVector(structured_grid%da_ndof,vector,ierr)
        case(NPHANCOMPDOF)
          call DACreateGlobalVector(structured_grid%da_nphancomp_dof,vector,ierr)
        case(NPHANSPECDOF)
          call DACreateGlobalVector(structured_grid%da_nphanspec_dof,vector,ierr)
        case(NPHANSPECNCOMPDOF)
          call DACreateGlobalVector(structured_grid%da_nphanspecncomp_dof,vector,ierr)
        case(VARDOF)
          call DACreateGlobalVector(structured_grid%da_var_dof,vector,ierr)
      end select    
    case(LOCAL)
      select case (da_index)
        case(ONEDOF)
          call DACreateLocalVector(structured_grid%da_1_dof,vector,ierr)
        case(NPHASEDOF)
          call DACreateLocalVector(structured_grid%da_nphase_dof,vector,ierr)
        case(THREENPDOF)
          call DACreateLocalVector(structured_grid%da_3np_dof,vector,ierr)
        case(NDOF)
          call DACreateLocalVector(structured_grid%da_ndof,vector,ierr)
        case(NPHANCOMPDOF)
          call DACreateLocalVector(structured_grid%da_nphancomp_dof,vector,ierr)
        case(NPHANSPECDOF)
          call DACreateLocalVector(structured_grid%da_nphanspec_dof,vector,ierr)
        case(NPHANSPECNCOMPDOF)
          call DACreateLocalVector(structured_grid%da_nphanspecncomp_dof,vector,ierr)
        case(VARDOF)
          call DACreateLocalVector(structured_grid%da_var_dof,vector,ierr)
      end select    
    case(NATURAL)
      select case (da_index)
        case(ONEDOF)
          call DACreateNaturalVector(structured_grid%da_1_dof,vector,ierr)
        case(NPHASEDOF)
          call DACreateNaturalVector(structured_grid%da_nphase_dof,vector,ierr)
        case(THREENPDOF)
          call DACreateNaturalVector(structured_grid%da_3np_dof,vector,ierr)
        case(NDOF)
          call DACreateNaturalVector(structured_grid%da_ndof,vector,ierr)
        case(NPHANCOMPDOF)
          call DACreateNaturalVector(structured_grid%da_nphancomp_dof,vector,ierr)
        case(NPHANSPECDOF)
          call DACreateNaturalVector(structured_grid%da_nphanspec_dof,vector,ierr)
        case(NPHANSPECNCOMPDOF)
          call DACreateNaturalVector(structured_grid%da_nphanspecncomp_dof,vector,ierr)
        case(VARDOF)
          call DACreateNaturalVector(structured_grid%da_var_dof,vector,ierr)
      end select    
  end select

end subroutine StructuredGridCreateVecFromDA

! ************************************************************************** !
!
! StructuredGridReadDXYZ: Reads structured grid spacing from input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine StructuredGridReadDXYZ(structured_grid,option)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  
  PetscInt :: i

  allocate(structured_grid%dx0(structured_grid%nx))
  allocate(structured_grid%dy0(structured_grid%ny))
  allocate(structured_grid%dz0(structured_grid%nz))
        
  call StructuredGridReadArray(structured_grid%dx0,structured_grid%nx,option)
  call StructuredGridReadArray(structured_grid%dy0,structured_grid%ny,option)
  call StructuredGridReadArray(structured_grid%dz0,structured_grid%nz,option)
    
  if (option%myrank==0) then
    write(IUNIT2,'(/," *DXYZ ")')
    write(IUNIT2,'("  dx  ",/,(1p10e12.4))') (structured_grid%dx0(i),i=1,structured_grid%nx)
    write(IUNIT2,'("  dy  ",/,(1p10e12.4))') (structured_grid%dy0(i),i=1,structured_grid%ny)
    write(IUNIT2,'("  dz  ",/,(1p10e12.4))') (structured_grid%dz0(i),i=1,structured_grid%nz)
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
subroutine StructuredGridReadArray(a,n,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: n
  PetscInt :: i, i1, i2, m
  PetscInt ::  ierr, nvalue=10
  PetscReal, intent(inout) :: a(*)
  character(len=MAXSTRINGLENGTH) :: string 

  save nvalue

!  call fiReadStringErrorMsg('DXYZ',ierr)

!  call fiReadDouble(string,grid%radius_0,ierr)
!  call fiDefaultMsg(option%myrank,'radius_0',ierr)

  i2 = 0
  do
    i1 = i2+1
    i2 = i2+nvalue
    if (i2.gt.n) i2 = n
    call fiReadFlotranString(IUNIT1,string,ierr)
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
subroutine StructuredGridComputeSpacing(structured_grid,nL2A)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  PetscInt :: nL2A(:)
  
  PetscInt :: i, j, k, n, na
  PetscReal, pointer :: dx_p(:), dy_p(:), dz_p(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(structured_grid%dx,dx_p,ierr)
  call VecGetArrayF90(structured_grid%dy,dy_p,ierr)
  call VecGetArrayF90(structured_grid%dz,dz_p,ierr)
  do n = 1,structured_grid%nlmax
    na = nL2A(n)
    k= int(na/structured_grid%nxy) + 1
    j= int(mod(na,structured_grid%nxy)/structured_grid%nx) + 1
    i= mod(mod(na,structured_grid%nxy),structured_grid%nx) + 1
    dx_p(n) = structured_grid%dx0(i)
    dy_p(n) = structured_grid%dy0(j)
    dz_p(n) = structured_grid%dz0(k)
  enddo
  call VecRestoreArrayF90(structured_grid%dx,dx_p,ierr)
  call VecRestoreArrayF90(structured_grid%dy,dy_p,ierr)
  call VecRestoreArrayF90(structured_grid%dz,dz_p,ierr)
  
  call DAGlobalToLocalBegin(structured_grid%da_1_dof, structured_grid%dx, INSERT_VALUES, &
                            structured_grid%dx_loc, ierr)
  call DAGlobalToLocalEnd(structured_grid%da_1_dof, structured_grid%dx, INSERT_VALUES, &
                          structured_grid%dx_loc, ierr)
  call DAGlobalToLocalBegin(structured_grid%da_1_dof, structured_grid%dy, INSERT_VALUES, &
                            structured_grid%dy_loc, ierr)
  call DAGlobalToLocalEnd(structured_grid%da_1_dof, structured_grid%dy, INSERT_VALUES, &
                          structured_grid%dy_loc,ierr)
  call DAGlobalToLocalBegin(structured_grid%da_1_dof, structured_grid%dz, INSERT_VALUES, &
                            structured_grid%dz_loc, ierr)
  call DAGlobalToLocalEnd(structured_grid%da_1_dof, structured_grid%dz, INSERT_VALUES, &
                          structured_grid%dz_loc,ierr)
  
end subroutine StructuredGridComputeSpacing

! ************************************************************************** !
!
! StructuredGridComputeCoord: Computes structured coordinates in x,y,z
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridComputeCoord(structured_grid,option, &
                                      origin,grid_x,grid_y,grid_z, &
                                      x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscReal :: origin(3)
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

! PetscInt :: ierr
  PetscInt :: i, j, k, n
  PetscReal :: x, y, z
  PetscInt :: prevnode

! set min and max bounds of domain in coordinate directions
  x_min = origin(X_DIRECTION)
  y_min = origin(Y_DIRECTION)
  z_min = origin(Z_DIRECTION)
  x_max = x_min
  y_max = y_min
  z_max = z_min
  do i=1,structured_grid%nx
    x_max = x_max + structured_grid%dx0(i)
  enddo
  do j=1,structured_grid%ny
    y_max = y_max + structured_grid%dy0(j)
  enddo
  do k=1,structured_grid%nz
    z_max = z_max + structured_grid%dz0(k)
  enddo

! set min and max bounds of domain in coordinate directions
  n = 0
  z = 0.5d0*structured_grid%dz0(1)
  do k=1, structured_grid%nz
    y = 0.5d0*structured_grid%dy0(1)
    do j=1, structured_grid%ny
      x = 0.5d0*structured_grid%dx0(1)
      do i=1, structured_grid%nx
        if (i > structured_grid%ngxs .and. i <= structured_grid%ngxe .and. &
            j > structured_grid%ngys .and. j <= structured_grid%ngye .and. &
            k > structured_grid%ngzs .and. k <= structured_grid%ngze) then
          n = n + 1
          grid_x(n) = x
          grid_y(n) = y
          grid_z(n) = z
        endif
        if (i < structured_grid%nx) x = x + 0.5d0*(structured_grid%dx0(i)+structured_grid%dx0(i+1))
      enddo
      if (j < structured_grid%ny) y = y + 0.5d0*(structured_grid%dy0(j)+structured_grid%dy0(j+1))
    enddo
    if (k < structured_grid%nz) z = z + 0.5d0*(structured_grid%dz0(k)+structured_grid%dz0(k+1))
  enddo
  if (n /= structured_grid%ngmax .and. option%myrank == 0) &
      print *, 'ERROR: Number of coordinates (',n, ') ', &
             'does not match number of ghosted cells (', structured_grid%ngmax, ')'             
    
end subroutine StructuredGridComputeCoord

! ************************************************************************** !
!
! StructGridComputeInternConnect: computes internal connectivity of a  
!                               structured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function StructGridComputeInternConnect(structured_grid,option)

  use Connection_module
  use Option_module
  
  implicit none
  
  type(connection_type), pointer :: StructGridComputeInternConnect
  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  
  
  PetscInt :: i, j, k, iconn, id_up, id_dn
  PetscReal :: dist_up, dist_dn
  type(connection_type), pointer :: connections
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  
  call VecGetArrayF90(structured_grid%dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dz_loc, dz_loc_p, ierr)
  
  call ConnectionAllocateLists()
  
  connections => &
       ConnectionCreate((structured_grid%ngx-1)*structured_grid%nly* &
                          structured_grid%nlz+ &
                        structured_grid%nlx*(structured_grid%ngy-1)* &
                          structured_grid%nlz+ &
                        structured_grid%nlx*structured_grid%nly* &
                          (structured_grid%ngz-1), &
                        option%nphase,INTERNAL_CONNECTION_TYPE)

  iconn = 0
  ! x-connections
  if (structured_grid%ngx > 1) then
    do k = structured_grid%kstart, structured_grid%kend
      do j = structured_grid%jstart, structured_grid%jend
        do i = 1, structured_grid%ngx - 1
          iconn = iconn+1
          id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy
          id_dn = id_up + 1
          connections%id_up(iconn) = id_up
          connections%id_dn(iconn) = id_dn
          connections%dist(-1:3,iconn) = 0.d0
          dist_up = 0.5d0*dx_loc_p(id_up)
          dist_dn = 0.5d0*dx_loc_p(id_dn)
          connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
          connections%dist(0,iconn) = dist_up+dist_dn
          connections%dist(1,iconn) = 1.d0  ! x component of unit vector
          connections%area(iconn) = dy_loc_p(id_up)*dz_loc_p(id_up)
        enddo
      enddo
    enddo
  endif

  ! y-connections
  if (structured_grid%ngy > 1) then
    do k = structured_grid%kstart, structured_grid%kend
      do i = structured_grid%istart, structured_grid%iend
        do j = 1, structured_grid%ngy - 1
          iconn = iconn+1
          id_up = i + 1 + (j-1) * structured_grid%ngx + k * structured_grid%ngxy
          id_dn = id_up + structured_grid%ngx
          connections%id_up(iconn) = id_up
          connections%id_dn(iconn) = id_dn
          connections%dist(-1:3,iconn) = 0.d0
          dist_up = 0.5d0*dy_loc_p(id_up)
          dist_dn = 0.5d0*dy_loc_p(id_dn)
          connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
          connections%dist(0,iconn) = dist_up+dist_dn
          connections%dist(2,iconn) = 1.d0  ! y component of unit vector
          connections%area(iconn) = dx_loc_p(id_up)*dz_loc_p(id_up)
        enddo
      enddo
    enddo
  endif
      
  ! z-connections
  if (structured_grid%ngz > 1) then
    do j = structured_grid%jstart, structured_grid%jend
      do i = structured_grid%istart, structured_grid%iend
        do k = 1, structured_grid%ngz - 1
          iconn = iconn+1
          id_up = i + 1 + j * structured_grid%ngx + (k-1) * structured_grid%ngxy
          id_dn = id_up + structured_grid%ngxy
          connections%id_up(iconn) = id_up
          connections%id_dn(iconn) = id_dn
          connections%dist(-1:3,iconn) = 0.d0
          dist_up = 0.5d0*dz_loc_p(id_up)
          dist_dn = 0.5d0*dz_loc_p(id_dn)
          connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
          connections%dist(0,iconn) = dist_up+dist_dn
          connections%dist(3,iconn) = 1.d0  ! z component of unit vector
          connections%area(iconn) = dx_loc_p(id_up)*dy_loc_p(id_up)
        enddo
      enddo
    enddo
  endif
  
  if (structured_grid%igeom == STRUCTURED_CYLINDRICAL) then
    allocate(structured_grid%rd(0:structured_grid%nx))
    structured_grid%rd = 0.D0
    structured_grid%rd(0) = structured_grid%Radius_0 
    do i = 1, structured_grid%nx
      structured_grid%rd(i) = structured_grid%rd(i-1) + structured_grid%dx0(i)
    enddo
  endif 

  call VecRestoreArrayF90(structured_grid%dx_loc, dx_loc_p, ierr)
  call VecRestoreArrayF90(structured_grid%dy_loc, dy_loc_p, ierr)
  call VecRestoreArrayF90(structured_grid%dz_loc, dz_loc_p, ierr)
  
  StructGridComputeInternConnect => connections

end function StructGridComputeInternConnect

! ************************************************************************** !
!
! StructGridPopulateConnection: Computes details of connection (area, dist, etc)
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine StructGridPopulateConnection(structured_grid,connection,iface, &
                                        iconn,cell_id_ghosted)

  use Connection_module
  use Option_module
  use Coupler_module
  use Region_module
  
  implicit none
 
  type(structured_grid_type) :: structured_grid
  type(connection_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: cell_id_ghosted
  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  
  call VecGetArrayF90(structured_grid%dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dz_loc, dz_loc_p, ierr)
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
      select case(structured_grid%igeom)
        case(STRUCTURED_CARTESIAN) ! cartesian
          select case(iface)
            case(WEST,EAST)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*dx_loc_p(cell_id_ghosted)
              connection%area(iconn) = dy_loc_p(cell_id_ghosted)* &
                                        dz_loc_p(cell_id_ghosted)
              if (iface ==  WEST) then
                connection%dist(1,iconn) = 1.d0
              else
                connection%dist(1,iconn) = -1.d0
              endif
            case(SOUTH,NORTH)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*dy_loc_p(cell_id_ghosted)
              connection%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                        dz_loc_p(cell_id_ghosted)
              if (iface ==  SOUTH) then
                connection%dist(2,iconn) = 1.d0
              else
                connection%dist(2,iconn) = -1.d0
              endif
            case(BOTTOM,TOP)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*dz_loc_p(cell_id_ghosted)
              connection%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                        dy_loc_p(cell_id_ghosted)
              if (structured_grid%invert_z_axis) then
                if (iface ==  TOP) then 
                  connection%dist(3,iconn) = 1.d0
                else
                  connection%dist(3,iconn) = -1.d0
                endif
              else
                if (iface ==  TOP) then 
                  connection%dist(3,iconn) = -1.d0
                else
                  connection%dist(3,iconn) = 1.d0
                endif
              endif
          end select
        case(STRUCTURED_CYLINDRICAL) ! cylindrical
        case(STRUCTURED_SPHERICAL) ! spherical
      end select
    case(INITIAL_CONNECTION_TYPE)
    case(SRC_SINK_CONNECTION_TYPE)
  end select
  
  call VecRestoreArrayF90(structured_grid%dx_loc, dx_loc_p, ierr)
  call VecRestoreArrayF90(structured_grid%dy_loc, dy_loc_p, ierr)
  call VecRestoreArrayF90(structured_grid%dz_loc, dz_loc_p, ierr)
    
end subroutine StructGridPopulateConnection

! ************************************************************************** !
!
! StructuredGridComputeVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine StructuredGridComputeVolumes(structured_grid,option,nL2G,volume)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: nL2G(:)
  Vec :: volume
  
  PetscReal, parameter :: Pi=3.1415926d0
  
  PetscInt :: i, n, ng
  PetscReal, pointer :: volume_p(:), dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(volume,volume_p, ierr)
  call VecGetArrayF90(structured_grid%dx_loc,dx_loc_p,ierr)
  call VecGetArrayF90(structured_grid%dy_loc,dy_loc_p,ierr)
  call VecGetArrayF90(structured_grid%dz_loc,dz_loc_p,ierr)
  do n=1, structured_grid%nlmax
    ng = nL2G(n)
    if (structured_grid%igeom == STRUCTURED_CARTESIAN) then
      volume_p(n) = dx_loc_p(ng) * dy_loc_p(ng) * dz_loc_p(ng)
    else if (structured_grid%igeom == STRUCTURED_CYLINDRICAL) then
      i = mod(mod((n),structured_grid%nlxy),structured_grid%nlx)!+(grid%ngxs-grid%nxs)
      if (i==0) i = structured_grid%nlx
      volume_p(n) = Pi * (structured_grid%rd(i+structured_grid%nxs) + structured_grid%rd(i-1+structured_grid%nxs))*&
      (structured_grid%rd(i+structured_grid%nxs) - structured_grid%rd(i-1+structured_grid%nxs)) * dz_loc_p(ng)
  !   print *, 'setup: Vol ', grid%myrank, n,ng,i, dz_loc_p(ng),grid%rd(i+grid%nxs),volume_p(n)
    else if (structured_grid%igeom == STRUCTURED_SPHERICAL) then
    endif
  enddo
  call VecRestoreArrayF90(volume,volume_p, ierr)
  call VecRestoreArrayF90(structured_grid%dx_loc,dx_loc_p,ierr)
  call VecRestoreArrayF90(structured_grid%dy_loc,dy_loc_p,ierr)
  call VecRestoreArrayF90(structured_grid%dz_loc,dz_loc_p,ierr)
  
  write(*,'(" rank= ",i3,", nlmax= ",i6,", nlx,y,z= ",3i4, &
    & ", nxs,e = ",2i4,", nys,e = ",2i4,", nzs,e = ",2i4)') &
    option%myrank,structured_grid%nlmax,structured_grid%nlx,structured_grid%nly,structured_grid%nlz, &
    structured_grid%nxs,structured_grid%nxe,structured_grid%nys,structured_grid%nye,structured_grid%nzs,structured_grid%nze

  write(*,'(" rank= ",i3,", ngmax= ",i6,", ngx,y,z= ",3i4, &
    & ", ngxs,e= ",2i4,", ngys,e= ",2i4,", ngzs,e= ",2i4)') &
    option%myrank,structured_grid%ngmax,structured_grid%ngx,structured_grid%ngy,structured_grid%ngz, &
    structured_grid%ngxs,structured_grid%ngxe,structured_grid%ngys,structured_grid%ngye,structured_grid%ngzs,structured_grid%ngze

end subroutine StructuredGridComputeVolumes

! ************************************************************************** !
!
! StructuredGridMapIndices: maps global, local and natural indices of cells 
!                          to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridMapIndices(structured_grid,nG2L,nL2G,nL2A,nG2A,nG2N)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  PetscInt, pointer :: nG2L(:), nL2G(:), nL2A(:), nG2A(:), nG2N(:)

  PetscInt :: i, j, k, n, ng, na, count1
  PetscErrorCode :: ierr
  
  allocate(nL2G(structured_grid%nlmax))
  allocate(nG2L(structured_grid%ngmax))
  allocate(nL2A(structured_grid%nlmax))
  allocate(nG2N(structured_grid%ngmax))
  allocate(nG2A(structured_grid%ngmax))
  
  structured_grid%istart = structured_grid%nxs-structured_grid%ngxs
  structured_grid%jstart = structured_grid%nys-structured_grid%ngys
  structured_grid%kstart = structured_grid%nzs-structured_grid%ngzs
  structured_grid%iend = structured_grid%istart+structured_grid%nlx-1
  structured_grid%jend = structured_grid%jstart+structured_grid%nly-1
  structured_grid%kend = structured_grid%kstart+structured_grid%nlz-1

  ! Local <-> Ghosted Transformation
  nG2L = 0  ! Must initialize this to zero!
  nL2G = 0
  nG2A = 0
  nL2A = 0
  nG2N = 0
  
  n = 0
  do k=structured_grid%kstart,structured_grid%kend
    do j=structured_grid%jstart,structured_grid%jend
      do i=structured_grid%istart,structured_grid%iend
        n = n + 1
        ng = i+j*structured_grid%ngx+k*structured_grid%ngxy+1
        nL2G(n) = ng
        nG2L(ng) = n
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
          ng = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
          nG2L(ng) = -1
        endif
      enddo
    enddo
  enddo

  n=0
  do k=1,structured_grid%nlz
    do j=1,structured_grid%nly
      do i=1,structured_grid%nlx
        n = n + 1
        na = i-1+structured_grid%nxs+(j-1+structured_grid%nys)*structured_grid%nx+(k-1+structured_grid%nzs)*structured_grid%nxy
        if (na>(structured_grid%nmax-1)) print *,'Wrong Nature order....'
        nL2A(n) = na
        !print *,grid%myrank, k,j,i,n,na
        !grid%nG2N(ng) = na
      enddo
    enddo
  enddo
  ! Local(ghosted)->Natural(natural order starts from 0)
  n=0
  do k=1,structured_grid%ngz
    do j=1,structured_grid%ngy
      do i=1,structured_grid%ngx
        n = n + 1
        na = i-1+structured_grid%ngxs+(j-1+structured_grid%ngys)*structured_grid%nx+(k-1+structured_grid%ngzs)*structured_grid%nxy
        if (na>(structured_grid%nmax-1)) print *,'Wrong Nature order....'
        nG2A(n) = na
!        print *,grid%myrank, k,j,i,n,na
        !grid%nG2N(ng) = na
      enddo
    enddo
  enddo
   
  call DAGetGlobalIndicesF90(structured_grid%da_1_dof,structured_grid%ngmax,nG2N, ierr)

end subroutine StructuredGridMapIndices

! ************************************************************************** !
!
! StructuredGridCreateJacobian: Creates Jacobian matrix associated with grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridCreateJacobian(structured_grid,solver,option)

  use Option_module
  use Solver_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(solver_type) :: solver
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  if (option%iblkfmt == 0) then
    call DAGetMatrix(structured_grid%da_ndof, MATAIJ, solver%J, ierr)
  else
    call DAGetMatrix(structured_grid%da_ndof, MATBAIJ, solver%J, ierr)
  endif
 ! call  MatSetBlocksize(grid%J,grid%ndof,ierr)
#if (PETSC_VERSION_RELEASE == 1)
  call MatSetOption(solver%J,MAT_KEEP_ZEROED_ROWS,ierr)
  call MatSetOption(solver%J,MAT_COLUMN_ORIENTED,ierr)
#else
  call MatSetOption(solver%J,MAT_KEEP_ZEROED_ROWS,PETSC_FALSE,ierr)
  call MatSetOption(solver%J,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
#endif
  
end subroutine StructuredGridCreateJacobian

! ************************************************************************** !
!
! StructuredGridCreateColoring: Creates ISColoring for grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridCreateColoring(structured_grid,option,coloring)

  use Option_module
  
  implicit none

#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  ISColoring :: coloring
  PetscErrorCode :: ierr
  
  call DAGetColoring(structured_grid%da_ndof,IS_COLORING_GLOBAL,coloring,ierr)

end subroutine StructuredGridCreateColoring

! ************************************************************************** !
!
! StructureGridGlobalToLocal: Performs global to local communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridGlobalToLocal(structured_grid,global_vec,local_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToLocalBegin(da_ptr,global_vec,INSERT_VALUES, &
                            local_vec,ierr)
  call DAGlobalToLocalEnd(da_ptr,global_vec, INSERT_VALUES, &
                          local_vec, ierr)
                          
end subroutine StructureGridGlobalToLocal

! ************************************************************************** !
!
! StructureGridLocalToGlobal: Performs local to global communication with DA
! author: Glenn Hamm8ond
! date: 1/2/07
!
! ************************************************************************** !
subroutine StructureGridLocalToGlobal(structured_grid,local_vec,global_vec, &
                                      da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DALocalToGlobal(da_ptr,local_vec,INSERT_VALUES, &
                       global_vec,ierr)
                          
end subroutine StructureGridLocalToGlobal

! ************************************************************************** !
!
! StructureGridLocalToLocal: Performs local to local communication with DA
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine StructureGridLocalToLocal(structured_grid,local_vec1,local_vec2, &
                                     da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DALocalToLocalBegin(da_ptr,local_vec1,INSERT_VALUES, &
                           local_vec2,ierr)
  call DALocalToLocalEnd(da_ptr,local_vec1, INSERT_VALUES, &
                         local_vec2, ierr)
                          
end subroutine StructureGridLocalToLocal

! ************************************************************************** !
!
! StructureGridGlobalToNatural: Performs global to natural communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridGlobalToNatural(structured_grid,global_vec,natural_vec, &
                                   da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToNaturalBegin(da_ptr,global_vec,INSERT_VALUES, &
                            natural_vec,ierr)
  call DAGlobalToNaturalEnd(da_ptr,global_vec,INSERT_VALUES, &
                          natural_vec, ierr)
                          
end subroutine StructureGridGlobalToNatural

! ************************************************************************** !
!
! StructureGridNaturalToGlobal: Performs natural to global communication with DA
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine StructureGridNaturalToGlobal(structured_grid,natural_vec, &
                                        global_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DANaturalToGlobalBegin(da_ptr,natural_vec,INSERT_VALUES, &
                              global_vec,ierr)
  call DANaturalToGlobalEnd(da_ptr,natural_vec,INSERT_VALUES, &
                            global_vec,ierr)
                          
end subroutine StructureGridNaturalToGlobal

! ************************************************************************** !
!
! StructureGridGlobalToLocalBegin: Begins global to local communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridGlobalToLocalBegin(structured_grid,global_vec, &
                                           local_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToLocalBegin(da_ptr,global_vec,INSERT_VALUES, &
                            local_vec,ierr)
                          
end subroutine StructureGridGlobalToLocalBegin

! ************************************************************************** !
!
! StructureGridLocalToLocalBegin: Begins local to local communication with DA
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine StructureGridLocalToLocalBegin(structured_grid,local_vec1, &
                                          local_vec2,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DALocalToLocalBegin(da_ptr,local_vec1,INSERT_VALUES, &
                           local_vec2,ierr)
                          
end subroutine StructureGridLocalToLocalBegin

! ************************************************************************** !
!
! StructureGridGlobalToNaturBegin: Begins global to natural communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridGlobalToNaturBegin(structured_grid,global_vec, &
                                           natural_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToNaturalBegin(da_ptr,global_vec,INSERT_VALUES, &
                            natural_vec,ierr)
                          
end subroutine StructureGridGlobalToNaturBegin

! ************************************************************************** !
!
! StructureGridGlobalToLocalEnd: Ends global to local communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridGlobalToLocalEnd(structured_grid,global_vec,local_vec, &
                                         da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToLocalEnd(da_ptr,global_vec,INSERT_VALUES, &
                          local_vec, ierr)
                          
end subroutine StructureGridGlobalToLocalEnd

! ************************************************************************** !
!
! StructureGridNaturToGlobalBegin: Begins global to natural communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridNaturToGlobalBegin(structured_grid,natural_vec, &
                                           global_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DANaturalToGlobalBegin(da_ptr,natural_vec,INSERT_VALUES, &
                              global_vec,ierr)
                          
end subroutine StructureGridNaturToGlobalBegin

! ************************************************************************** !
!
! StructureGridLocalToLocalEnd: Ends local to local communication with DA
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine StructureGridLocalToLocalEnd(structured_grid,local_vec1,local_vec2, &
                                     da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DALocalToLocalEnd(da_ptr,local_vec1,INSERT_VALUES, &
                         local_vec2, ierr)
                          
end subroutine StructureGridLocalToLocalEnd

! ************************************************************************** !
!
! StructureGridGlobalToNaturEnd: Ends global to natural communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridGlobalToNaturEnd(structured_grid,global_vec, &
                                         natural_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToNaturalEnd(da_ptr,global_vec,INSERT_VALUES, &
                          natural_vec, ierr)
                          
end subroutine StructureGridGlobalToNaturEnd

! ************************************************************************** !
!
! StructureGridNaturToGlobalEnd: Ends global to natural communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructureGridNaturToGlobalEnd(structured_grid,natural_vec, &
                                         global_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = StructGridGetDAPtrFromIndex(structured_grid,da_index)

  call DANaturalToGlobalEnd(da_ptr,natural_vec,INSERT_VALUES, &
                            global_vec,ierr)
                          
end subroutine StructureGridNaturToGlobalEnd

! ************************************************************************** !
!
! StructGridGetDAPtrFromIndex: Returns the integer pointer for the DA referenced
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
function StructGridGetDAPtrFromIndex(structured_grid,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  PetscInt :: da_index
  
  DA :: StructGridGetDAPtrFromIndex
  
  select case (da_index)
    case(ONEDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_1_dof
    case(NPHASEDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_nphase_dof
    case(THREENPDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_3np_dof
    case(NDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_ndof
    case(NPHANCOMPDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_nphancomp_dof
    case(NPHANSPECDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_nphanspec_dof
    case(NPHANSPECNCOMPDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_nphanspecncomp_dof
    case(VARDOF)
      StructGridGetDAPtrFromIndex = structured_grid%da_var_dof
  end select  
  
end function StructGridGetDAPtrFromIndex

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
  

  if (associated(structured_grid%dx0)) deallocate(structured_grid%dx0)
  nullify(structured_grid%dx0)
  if (associated(structured_grid%dy0)) deallocate(structured_grid%dy0)
  nullify(structured_grid%dy0)
  if (associated(structured_grid%dz0)) deallocate(structured_grid%dz0)
  nullify(structured_grid%dz0)
  if (associated(structured_grid%rd)) deallocate(structured_grid%rd)
  nullify(structured_grid%rd)

  if (structured_grid%dx /= 0) call VecDestroy(structured_grid%dx,ierr)
  if (structured_grid%dy /= 0) call VecDestroy(structured_grid%dy,ierr)
  if (structured_grid%dz /= 0) call VecDestroy(structured_grid%dz,ierr)
  if (structured_grid%dx_loc /= 0) call VecDestroy(structured_grid%dx_loc,ierr)
  if (structured_grid%dy_loc /= 0) call VecDestroy(structured_grid%dy_loc,ierr)
  if (structured_grid%dz_loc /= 0) call VecDestroy(structured_grid%dz_loc,ierr)

  if (structured_grid%da_1_dof /= 0) &
    call DADestroy(structured_grid%da_1_dof,ierr)
  if (structured_grid%da_nphase_dof /= 0) &
    call DADestroy(structured_grid%da_nphase_dof,ierr)
  if (structured_grid%da_3np_dof /= 0) &
    call DADestroy(structured_grid%da_3np_dof,ierr)
  if (structured_grid%da_ndof /= 0) &
    call DADestroy(structured_grid%da_ndof,ierr)
  if (structured_grid%da_nphancomp_dof /= 0) &
    call DADestroy(structured_grid%da_nphancomp_dof,ierr)
  if (structured_grid%da_nphanspec_dof /= 0) &
    call DADestroy(structured_grid%da_nphanspec_dof,ierr)
  if (structured_grid%da_nphanspecncomp_dof /= 0) &
    call DADestroy(structured_grid%da_nphanspecncomp_dof,ierr)
  if (structured_grid%da_var_dof /= 0) &
    call DADestroy(structured_grid%da_var_dof,ierr)

  deallocate(structured_grid)
  nullify(structured_grid)

end subroutine StructuredGridDestroy
                          
end module Structured_Grid_module
