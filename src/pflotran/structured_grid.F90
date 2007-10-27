module Structured_Grid_module

  implicit none
 
  private
 
#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
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


  type, public :: structured_grid_type

    integer :: nx, ny, nz    ! Global domain dimensions of the grid.
    integer :: nxy, nmax     ! nx * ny, nx * ny * nz
    integer :: npx, npy, npz ! Processor partition in each direction.
    integer :: nlx, nly, nlz ! Local grid dimension w/o ghost nodes.
    integer :: ngx, ngy, ngz ! Local grid dimension with ghost nodes.
    integer :: nxs, nys, nzs 
      ! Global indices of non-ghosted corner (starting) of local domain.
    integer :: ngxs, ngys, ngzs
      ! Global indices of ghosted starting corner of local domain.
    integer :: nxe, nye, nze, ngxe, ngye, ngze
      ! Global indices of non-ghosted/ghosted ending corner of local domain.
    integer :: nlxy, nlxz, nlyz
    integer :: ngxy, ngxz, ngyz
    
    integer :: istart, jstart, kstart, iend, jend, kend
      ! istart gives the local x-index of the non-ghosted starting (lower left)
      ! corner. iend gives the local x-index of the non-ghosted ending 
      ! corner. jstart, jend correspond to y-index, kstart, kend to z-index.    

    integer :: nlmax  ! Total number of non-ghosted nodes in local domain.
    integer :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    
    real*8 :: x_max, x_min, y_max, y_min, z_max, z_min

    real*8, pointer :: dx0(:), dy0(:), dz0(:), rd(:)
    
    integer :: igeom
    
    real*8 :: radius_0
    
    Vec :: dx, dy, dz, dx_loc, dy_loc, dz_loc  ! Grid spacings

    DA :: da_1_dof, da_nphase_dof, da_3np_dof, da_ndof, da_nphancomp_dof, &
          da_nphanspec_dof, da_nphanspecncomp_dof, da_var_dof 
    
  end type

  public :: initStructuredGrid, &
            createStructuredDMs, &
            computeStructInternalConnect, &
            computeStructBoundaryConnect, &
            createPetscVectorFromDA, &
            mapStructuredGridIndices, &
            computeStructuredGridSpacing, &
            computeStructuredGridCoordinates, &
            createStructuredGridJacobian, &
            createStructuredGridColoring, &
            DMStructGlobalToLocal, &
            DMStructGlobalToNatural, &
            readStructuredDXYZ, &
            computeStructuredCellVolumes

contains

! ************************************************************************** !
!
! initStructuredGrid: Initializes a structured grid object
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine initStructuredGrid(structured_grid)

  implicit none
  
  type(structured_grid_type) :: structured_grid

  structured_grid%nx = 0
  structured_grid%ny = 0
  structured_grid%nz = 0
  structured_grid%npx = PETSC_DECIDE
  structured_grid%npy = PETSC_DECIDE
  structured_grid%npz = PETSC_DECIDE
  
end subroutine initStructuredGrid
  
! ************************************************************************** !
!
! createStructuredDMs: Creates structured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine createStructuredDMs(structured_grid,option)
      
  use Option_module
      
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option

  integer :: ndof
  integer, parameter :: stencil_width = 1
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

  ndof = 3*option%nphase
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  structured_grid%da_3np_dof,ierr)

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
    case(MPH_MODE,RICHARDS_MODE,VADOSE_MODE,FLASH_MODE,OWG_MODE)
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

end subroutine createStructuredDMs

! ************************************************************************** !
!
! createPetscVectorFromDA: Creates a global PETSc vector
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine createPetscVectorFromDA(structured_grid,da_index,vector, &
                                   vector_type)

  implicit none
  type(structured_grid_type) :: structured_grid
  integer :: da_index
  Vec :: vector
  integer :: vector_type
  
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

end subroutine createPetscVectorFromDA

! ************************************************************************** !
!
! readStructuredDXYZ: Reads structured grid spacing from input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readStructuredDXYZ(structured_grid,option)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  
  integer :: i

  allocate(structured_grid%dx0(structured_grid%nx))
  allocate(structured_grid%dy0(structured_grid%ny))
  allocate(structured_grid%dz0(structured_grid%nz))
        
  call readDXYZ(structured_grid%dx0,structured_grid%nx)
  call readDXYZ(structured_grid%dy0,structured_grid%ny)
  call readDXYZ(structured_grid%dz0,structured_grid%nz)
    
  if (option%myrank==0) then
    write(IUNIT2,'(/," *DXYZ ")')
    write(IUNIT2,'("  dx  ",/,(1p10e12.4))') (structured_grid%dx0(i),i=1,structured_grid%nx)
    write(IUNIT2,'("  dy  ",/,(1p10e12.4))') (structured_grid%dy0(i),i=1,structured_grid%ny)
    write(IUNIT2,'("  dz  ",/,(1p10e12.4))') (structured_grid%dz0(i),i=1,structured_grid%nz)
  endif

end subroutine readStructuredDXYZ

! ************************************************************************** !
!
! readStructuredDXYZ: Reads structured grid spacing along an axis from input 
!                     file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readDXYZ(a,n)

  use fileio_module
  
  implicit none
  
  integer*4, intent(in) :: n
  integer*4 :: i, i1, i2, m
  integer ::  ierr, nvalue=10
  real*8, intent(inout) :: a(*)
  character(len=MAXSTRINGLENGTH) :: string 

  save nvalue

!  call fiReadStringErrorMsg('DXYZ',ierr)

!  call fiReadDouble(string,grid%radius_0,ierr)
!  call fiDefaultMsg('radius_0',ierr)

  i2 = 0
  do
    i1 = i2+1
    i2 = i2+nvalue
    if (i2.gt.n) i2 = n
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg('DXYZ',ierr)
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
    
end subroutine readDXYZ

! ************************************************************************** !
!
! computeStructuredGridSpacing: Computes structured grid spacing
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine computeStructuredGridSpacing(structured_grid,nL2A)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  integer :: nL2A(:)
  
  integer :: i, j, k, n, na
  PetscScalar, pointer :: dx_p(:), dy_p(:), dz_p(:)
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
  
end subroutine computeStructuredGridSpacing

! ************************************************************************** !
!
! computeStructuredGridCoordinates: Computes structured coordinates in x,y,z
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine computeStructuredGridCoordinates(structured_grid,option, &
                                            grid_x,grid_y,grid_z)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  real*8 :: grid_x(:), grid_y(:), grid_z(:)
  real*8 :: x_min, y_min, z_min

! integer :: ierr
  integer*4 :: i, j, k, n
  real*8 :: x, y, z
  integer :: prevnode

! set min and max bounds of domain in coordinate directions
  structured_grid%x_min = 0.d0
  structured_grid%y_min = 0.d0
  structured_grid%z_min = 0.d0
  structured_grid%x_max = 0.d0
  structured_grid%y_max = 0.d0
  structured_grid%z_max = 0.d0
  do i=1,structured_grid%nx
    structured_grid%x_max = structured_grid%x_max + structured_grid%dx0(i)
  enddo
  do j=1,structured_grid%ny
    structured_grid%y_max = structured_grid%y_max + structured_grid%dy0(j)
  enddo
  do k=1,structured_grid%nz
    structured_grid%z_max = structured_grid%z_max + structured_grid%dz0(k)
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
    
end subroutine computeStructuredGridCoordinates

! ************************************************************************** !
!
! computeStructInternalConnect: computes internal connectivity of a  
!                               structured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function computeStructInternalConnect(structured_grid,option)

  use Connection_module
  use option_module
  
  implicit none
  
  type(connection_type), pointer :: computeStructInternalConnect
  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  
  
  integer :: i, j, k, iconn, id_up, id_dn
  real*8 :: dist_up, dist_dn
  type(connection_type), pointer :: connections
  PetscErrorCode :: ierr
  
  PetscScalar, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  
  call VecGetArrayF90(structured_grid%dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dz_loc, dz_loc_p, ierr)
  
  call allocateConnectionLists()
  
  connections => &
       createConnection((structured_grid%ngx-1)*structured_grid%nly* &
                          structured_grid%nlz+ &
                        structured_grid%nlx*(structured_grid%ngy-1)* &
                          structured_grid%nlz+ &
                        structured_grid%nlx*structured_grid%nly* &
                          (structured_grid%ngz-1), &
                        option%nphase)

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
          connections%dist(2,iconn) = 1.d0  ! x component of unit vector
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
          connections%dist(3,iconn) = 1.d0  ! x component of unit vector
          connections%area(iconn) = dx_loc_p(id_up)*dy_loc_p(id_up)
        enddo
      enddo
    enddo
  endif
  
  if (structured_grid%igeom == 2) then
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
  
  computeStructInternalConnect => connections

end function computeStructInternalConnect

! ************************************************************************** !
!
! computeStructBoundaryConnect: computes boundary connectivity of a 
!                               structured grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function computeStructBoundaryConnect(structured_grid,option,ibconn,nL2G)

  use Connection_module
  use Option_module
  
  implicit none

  type(connection_type), pointer :: computeStructBoundaryConnect  
  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  integer, pointer :: ibconn(:)
  integer :: nL2G(:)
  
  integer :: num_conn_hypothetically
  integer :: i, j, k, iconn
  integer :: cell_id_local, cell_id_ghosted
  integer :: ii1, ii2, jj1, jj2, kk1, kk2
  integer :: ibc, ir
  type(connection_type), pointer :: connections
  PetscErrorCode :: ierr
  
  PetscScalar, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  
  call VecGetArrayF90(structured_grid%dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(structured_grid%dz_loc, dz_loc_p, ierr)
  
  iconn = 0
  if (structured_grid%nx > 1 .and. structured_grid%ny == 1 .and. &
      structured_grid%nz == 1) then
    if (structured_grid%nxs == structured_grid%ngxs) iconn = iconn + 1
    if (structured_grid%nxe == structured_grid%ngxe) iconn = iconn + 1
  else if (structured_grid%nx == 1 .and. structured_grid%ny == 1 .and. &
           structured_grid%nz > 1) then
    if (structured_grid%nzs == structured_grid%ngzs) iconn = iconn + 1
    if (structured_grid%nze == structured_grid%ngze) iconn = iconn + 1
  else
    if (structured_grid%nx > 1) then
      if (structured_grid%nxs == structured_grid%ngxs) &
        iconn = iconn + structured_grid%nlyz
      if (structured_grid%nxe == structured_grid%ngxe) &
        iconn = iconn + structured_grid%nlyz
    endif
    if (structured_grid%ny > 1) then
      if (structured_grid%nys == structured_grid%ngys) &
        iconn = iconn + structured_grid%nlxz
      if (structured_grid%nye == structured_grid%ngye) &
        iconn = iconn + structured_grid%nlxz
    endif
    if (structured_grid%nz > 1) then
      if (structured_grid%nzs == structured_grid%ngzs) &
        iconn = iconn + structured_grid%nlxy
      if (structured_grid%nze == structured_grid%ngze) &
        iconn = iconn + structured_grid%nlxy
    endif
  endif

  num_conn_hypothetically = iconn
  connections => createConnection(iconn,option%nphase)

  allocate(ibconn(iconn)) 

  iconn = 0
  if (structured_grid%nxs == structured_grid%ngxs .or. &
      structured_grid%nxe == structured_grid%ngxe .or. &
      structured_grid%nys == structured_grid%ngys .or. &
      structured_grid%nye == structured_grid%ngye .or. &
      structured_grid%nzs == structured_grid%ngzs .or. &
      structured_grid%nze == structured_grid%ngze) then

    ! calculate boundary conditions locally on only those processors which 
    ! contain a boundary!

    do ibc = 1, option%nblkbc
      do ir = option%iregbc1(ibc), option%iregbc2(ibc)
        kk1 = option%k1bc(ir) - structured_grid%nzs
        kk2 = option%k2bc(ir) - structured_grid%nzs
        jj1 = option%j1bc(ir) - structured_grid%nys
        jj2 = option%j2bc(ir) - structured_grid%nys
        ii1 = option%i1bc(ir) - structured_grid%nxs
        ii2 = option%i2bc(ir) - structured_grid%nxs

        kk1 = max(1,kk1)
        kk2 = min(structured_grid%nlz,kk2)
        jj1 = max(1,jj1)
        jj2 = min(structured_grid%nly,jj2)
        ii1 = max(1,ii1)
        ii2 = min(structured_grid%nlx,ii2)

        if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle 

        do k = kk1,kk2
          do j = jj1,jj2
            do i = ii1,ii2
              iconn = iconn + 1
              cell_id_local = i+(j-1)*structured_grid%nlx+(k-1)*structured_grid%nlxy

              connections%id_dn(iconn) = cell_id_local  ! m is a local index
              
              ! old way
              !grid%ibconn(nc) = ibc
              ibconn(iconn) = ir
              cell_id_ghosted = nL2G(cell_id_local)
            ! Use ghosted index to access dx, dy, dz because we have
            ! already done a global-to-local scatter for computing the
            ! interior node connections.
      
!               print *,'pflowgrid_mod: ',nc,ibc,ir,m,ng,ii1,ii2,kk1,kk2, &
!                        grid%nblkbc,grid%igeom
        
              select case(structured_grid%igeom)
                case(1) ! cartesian
                  if (option%iface(ibc) == 1) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dx_loc_p(cell_id_ghosted)
                    connections%dist(1,iconn) = 1.d0
                    connections%area(iconn) = dy_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (option%iface(ibc) == 2) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dx_loc_p(cell_id_ghosted)
                    connections%dist(1,iconn) = 1.d0
                    connections%area(iconn) = dy_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (option%iface(ibc) == 3) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dz_loc_p(cell_id_ghosted)
                    connections%dist(3,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dy_loc_p(cell_id_ghosted)
                  else if (option%iface(ibc) == 4) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dz_loc_p(cell_id_ghosted)
                    connections%dist(3,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dy_loc_p(cell_id_ghosted)
                  else if (option%iface(ibc) == 5) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dy_loc_p(cell_id_ghosted)
                    connections%dist(2,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (option%iface(ibc) == 6) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dy_loc_p(cell_id_ghosted)
                    connections%dist(2,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  endif
                case(2) ! cylindrical
                case(3) ! spherical
              end select
            enddo ! i
          enddo ! j
        enddo ! k
      enddo ! ir
    enddo ! ibc
  endif
  
  if (num_conn_hypothetically /= iconn) then
    print *, 'ERROR: Number of actual connections does not match hypothetical value'
    stop
  endif

  call VecRestoreArrayF90(structured_grid%dx_loc, dx_loc_p, ierr)
  call VecRestoreArrayF90(structured_grid%dy_loc, dy_loc_p, ierr)
  call VecRestoreArrayF90(structured_grid%dz_loc, dz_loc_p, ierr)
  
  computeStructBoundaryConnect => connections
  
end function computeStructBoundaryConnect

! ************************************************************************** !
!
! computeStructuredCellVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine computeStructuredCellVolumes(structured_grid,option,nL2G)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  integer :: nL2G(:)
  
  real*8, parameter :: Pi=3.1415926d0
  
  integer :: i, n, ng
  PetscScalar, pointer :: volume_p(:), dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(option%volume,volume_p, ierr)
  call VecGetArrayF90(structured_grid%dx_loc,dx_loc_p,ierr)
  call VecGetArrayF90(structured_grid%dy_loc,dy_loc_p,ierr)
  call VecGetArrayF90(structured_grid%dz_loc,dz_loc_p,ierr)
  do n=1, structured_grid%nlmax
    ng = nL2G(n)
    if (structured_grid%igeom == 1) then
      volume_p(n) = dx_loc_p(ng) * dy_loc_p(ng) * dz_loc_p(ng)
    else if (structured_grid%igeom == 2) then
      i = mod(mod((n),structured_grid%nlxy),structured_grid%nlx)!+(grid%ngxs-grid%nxs)
      if (i==0) i = structured_grid%nlx
      volume_p(n) = Pi * (structured_grid%rd(i+structured_grid%nxs) + structured_grid%rd(i-1+structured_grid%nxs))*&
      (structured_grid%rd(i+structured_grid%nxs) - structured_grid%rd(i-1+structured_grid%nxs)) * dz_loc_p(ng)
  !   print *, 'setup: Vol ', grid%myrank, n,ng,i, dz_loc_p(ng),grid%rd(i+grid%nxs),volume_p(n)
    else if (structured_grid%igeom == 3) then
    endif
  enddo
  call VecRestoreArrayF90(option%volume,volume_p, ierr)
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

end subroutine computeStructuredCellVolumes

! ************************************************************************** !
!
! mapStructuredGridIndices: maps global, local and natural indices of cells 
!                          to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine mapStructuredGridIndices(structured_grid,nG2L,nL2G,nL2A,nG2A,nG2N)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  integer, pointer :: nG2L(:), nL2G(:), nL2A(:), nG2A(:), nG2N(:)

  integer :: i, j, k, n, ng, na
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

end subroutine mapStructuredGridIndices

! ************************************************************************** !
!
! createStructuredGridJacobian: Creates Jacobian matrix associated with grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine createStructuredGridJacobian(structured_grid,option)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  if (option%iblkfmt == 0) then
    call DAGetMatrix(structured_grid%da_ndof, MATAIJ, option%J, ierr)
  else
    call DAGetMatrix(structured_grid%da_ndof, MATBAIJ, option%J, ierr)
  endif
 ! call  MatSetBlocksize(grid%J,grid%ndof,ierr)
  call MatSetOption(option%J,MAT_KEEP_ZEROED_ROWS,ierr)
  call MatSetOption(option%J,MAT_COLUMN_ORIENTED,ierr)
  
end subroutine createStructuredGridJacobian

! ************************************************************************** !
!
! createStructuredGridColoring: Creates ISColoring for grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine createStructuredGridColoring(structured_grid,option,coloring)

  use Option_module
  
  implicit none

#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  ISColoring :: coloring
  PetscErrorCode :: ierr
  
  call DAGetColoring(structured_grid%da_ndof,IS_COLORING_GLOBAL,coloring,ierr)

end subroutine createStructuredGridColoring

! ************************************************************************** !
!
! DMStructGlobalToLocal: Performs global to local communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DMStructGlobalToLocal(structured_grid,global_vec,local_vec,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: local_vec
  integer :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = getDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToLocalBegin(da_ptr,global_vec,INSERT_VALUES, &
                            local_vec,ierr)
  call DAGlobalToLocalEnd(da_ptr,global_vec, INSERT_VALUES, &
                          local_vec, ierr)
                          
end subroutine DMStructGlobalToLocal

! ************************************************************************** !
!
! DMStructGlobalToNatural: Performs global to natural communication with DA
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DMStructGlobalToNatural(structured_grid,global_vec,natural_vec, &
                                   da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  Vec :: global_vec
  Vec :: natural_vec
  integer :: da_index
  
  DA :: da_ptr
  PetscErrorCode :: ierr

  da_ptr = getDAPtrFromIndex(structured_grid,da_index)

  call DAGlobalToNaturalBegin(da_ptr,global_vec,INSERT_VALUES, &
                            natural_vec,ierr)
  call DAGlobalToNaturalEnd(da_ptr,global_vec, INSERT_VALUES, &
                          natural_vec, ierr)
                          
end subroutine DMStructGlobalToNatural

! ************************************************************************** !
!
! getDAPtrFromIndex: Returns the integer pointer for the DA referenced
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
function getDAPtrFromIndex(structured_grid,da_index)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  integer :: da_index
  
  DA :: getDAPtrFromIndex
  
  select case (da_index)
    case(ONEDOF)
      getDAPtrFromIndex = structured_grid%da_1_dof
    case(NPHASEDOF)
      getDAPtrFromIndex = structured_grid%da_nphase_dof
    case(THREENPDOF)
      getDAPtrFromIndex = structured_grid%da_3np_dof
    case(NDOF)
      getDAPtrFromIndex = structured_grid%da_ndof
    case(NPHANCOMPDOF)
      getDAPtrFromIndex = structured_grid%da_nphancomp_dof
    case(NPHANSPECDOF)
      getDAPtrFromIndex = structured_grid%da_nphanspec_dof
    case(NPHANSPECNCOMPDOF)
      getDAPtrFromIndex = structured_grid%da_nphanspecncomp_dof
    case(VARDOF)
      getDAPtrFromIndex = structured_grid%da_var_dof
  end select  
  
end function getDAPtrFromIndex
                          
end module Structured_Grid_module
