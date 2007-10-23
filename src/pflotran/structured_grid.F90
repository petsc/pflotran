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

    integer*4 :: nx, ny, nz    ! Global domain dimensions of the grid.
    integer*4 :: nxy, nmax     ! nx * ny, nx * ny * nz
    integer*4 :: npx, npy, npz ! Processor partition in each direction.
    integer*4 :: nlx, nly, nlz ! Local grid dimension w/o ghost nodes.
    integer*4 :: ngx, ngy, ngz ! Local grid dimension with ghost nodes.
    integer*4 :: nxs, nys, nzs 
      ! Global indices of non-ghosted corner (starting) of local domain.
    integer*4 :: ngxs, ngys, ngzs
      ! Global indices of ghosted starting corner of local domain.
    integer*4 :: nxe, nye, nze, ngxe, ngye, ngze
      ! Global indices of non-ghosted/ghosted ending corner of local domain.
    integer*4 :: nlxy, nlxz, nlyz
    integer*4 :: ngxy, ngxz, ngyz
    
    integer*4 :: istart, jstart, kstart, iend, jend, kend
      ! istart gives the local x-index of the non-ghosted starting (lower left)
      ! corner. iend gives the local x-index of the non-ghosted ending 
      ! corner. jstart, jend correspond to y-index, kstart, kend to z-index.    

    real*8, pointer :: dx0(:), dy0(:), dz0(:), rd(:)
    
    integer :: igeom
    
    real*8 :: radius_0
    
    
    integer*4 :: nblkbc
      ! The number of "blocks" of boundary conditions that are defined.
      ! Such a block is a specification of a set of boundary conditions.
      ! This set of boundary conditions can apply to any number of regions,
      ! so nblkbc does NOT equal the number of boundary condition regions.
    integer*4 :: nconnbc  ! The number of interfaces along boundaries.
!GEH - Structured Grid Dependence - Begin
    integer*4, pointer :: i1bc(:), i2bc(:), j1bc(:), j2bc(:), k1bc(:), k2bc(:)
!GEH - Structured Grid Dependence - End
#if 0
    integer*4, pointer :: ibconn(:)
#endif
      ! ibconn(nc) specifies the id of the of boundary condition block that
      ! applies at boundary interface nc.  
    integer*4, pointer :: ibndtyp(:)
      ! ibndtyp(ibc) specifies the type of boundary condition that applies
      ! for boundary condition block ibc.
    integer*4, pointer :: iface(:)
      ! iface(ibc) specifies the face (left, right, top, bottom, etc.) on
      ! which BC block ibc lies.
    integer*4, pointer :: mblkbc(:)
      ! mblkbc(nc) gives the local, non-ghosted index of the cell that has
      ! boundary connection nc.
    integer*4, pointer :: iregbc1(:), iregbc2(:)
      ! iregbc1(ibc) and iregbc2(ibc) give the id of the first region and 
      ! last region, respectively, that utilizes the boundary conditions in 
      ! boundary condition block ibc.
          
    Vec :: dx, dy, dz, dx_loc, dy_loc, dz_loc  ! Grid spacings

    DA :: da_1_dof, da_nphase_dof, da_3np_dof, da_ndof, da_nphancomp_dof, &
          da_nphanspec_dof, da_nphanspecncomp_dof, da_var_dof 
    
  end type


  public :: createStructuredDMs, &
            computeStructInternalConnect, &
            computeStructBoundaryConnect
  
contains

! ************************************************************************** !
!
! initStructuredGrid: Initializes a structured grid object
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine initStructuredGrid(grid)

  implicit none
  
  type(structured_grid_type) :: grid

  grid%nx = 0
  grid%ny = 0
  grid%nz = 0
  grid%npx = PETSC_DECIDE
  grid%npy = PETSC_DECIDE
  grid%npz = PETSC_DECIDE
  
end subroutine initStructuredGrid
  
! ************************************************************************** !
!
! createStructuredDMs: Creates structured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine createStructuredDMs(solution,structured_grid)
      
  use Solution_module
      
  implicit none
  
  type(solution_type) :: solution
  type(structured_grid_type) :: structured_grid
  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DA objects that will manage communication.
  !-----------------------------------------------------------------------
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,1,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       structured_grid%da_1_dof,ierr)

! call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
!      structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,3,1, &
!      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
!      structured_grid%da_3_dof,ierr)
 
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,solution%nphase,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       structured_grid%da_nphase_dof,ierr)

  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,3*solution%nphase,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       structured_grid%da_3np_dof,ierr)

  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,solution%ndof,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       structured_grid%da_ndof,ierr)

  if (solution%use_2ph == PETSC_TRUE) then
    print *,' 2ph create DA'
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,solution%nphase*solution%npricomp,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_nphancomp_dof,ierr)

    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,solution%nphase*solution%nspec,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_nphanspec_dof,ierr)

    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                    solution%nphase*solution%nspec*solution%npricomp,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_nphanspecncomp_dof,ierr)
  endif
       
  if (solution%use_mph == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,(solution%ndof+1)*(2+7*solution%nphase + 2* &
                    solution%nspec*solution%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_var_dof,ierr)

  endif
  
 if (solution%use_richards == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,(solution%ndof+1)*(2+7*solution%nphase + 2* &
                    solution%nspec*solution%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_var_dof,ierr)

  endif
    
  if (solution%use_vadose == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,(solution%ndof+1)*(2+7*solution%nphase + 2* &
                    solution%nspec*solution%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_var_dof,ierr)
  endif

  if (solution%use_flash == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,(solution%ndof+1)*(2+7*solution%nphase + 2* &
                    solution%nspec*solution%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_var_dof,ierr)
  endif
     
  if (solution%use_owg == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz,structured_grid%npx,structured_grid%npy,structured_grid%npz,(solution%ndof+1)*(2+7*solution%nphase + 2* &
                    solution%nspec*solution%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    structured_grid%da_var_dof,ierr)
  endif

  !print DA info for each processor
! call DAView(structured_grid%da_ndof,PETSC_VIEWER_STDOUT_WORLD,ierr)

end subroutine createStructuredDMs

! ************************************************************************** !
!
! readStructuredDXYZ: Reads structured grid spacing from input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readStructuredDXYZ(solution,structured_grid)

  use Solution_module
  
  implicit none
  
  type(solution_type) :: solution
  type(structured_grid_type) :: structured_grid
  
  integer :: i

  allocate(structured_grid%dx0(structured_grid%nx))
  allocate(structured_grid%dy0(structured_grid%ny))
  allocate(structured_grid%dz0(structured_grid%nz))
        
  call readDXYZ(structured_grid%dx0,structured_grid%nx)
  call readDXYZ(structured_grid%dy0,structured_grid%ny)
  call readDXYZ(structured_grid%dz0,structured_grid%nz)
    
  if (solution%myrank==0) then
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
#include "definitions.h"
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
! computeStructInternalConnect: computes internal connectivity of a  
!                               structured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function computeStructInternalConnect(solution,structured_grid)

  use Connection_module
  use Solution_module
  
  implicit none
  
  type(connection_type), pointer :: computeStructInternalConnect
  type(solution_type) :: solution
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
                        solution%nphase)

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
function computeStructBoundaryConnect(solution,structured_grid,ibconn,nL2G)

  use Connection_module
  use Solution_module
  
  implicit none

  type(connection_type), pointer :: computeStructBoundaryConnect  
  type(solution_type) :: solution
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
  connections => createConnection(iconn,solution%nphase)

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

    do ibc = 1, structured_grid%nblkbc
      do ir = structured_grid%iregbc1(ibc), structured_grid%iregbc2(ibc)
        kk1 = structured_grid%k1bc(ir) - structured_grid%nzs
        kk2 = structured_grid%k2bc(ir) - structured_grid%nzs
        jj1 = structured_grid%j1bc(ir) - structured_grid%nys
        jj2 = structured_grid%j2bc(ir) - structured_grid%nys
        ii1 = structured_grid%i1bc(ir) - structured_grid%nxs
        ii2 = structured_grid%i2bc(ir) - structured_grid%nxs

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
                  if (structured_grid%iface(ibc) == 1) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dx_loc_p(cell_id_ghosted)
                    connections%dist(1,iconn) = 1.d0
                    connections%area(iconn) = dy_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (structured_grid%iface(ibc) == 2) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dx_loc_p(cell_id_ghosted)
                    connections%dist(1,iconn) = 1.d0
                    connections%area(iconn) = dy_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (structured_grid%iface(ibc) == 3) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dz_loc_p(cell_id_ghosted)
                    connections%dist(3,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dy_loc_p(cell_id_ghosted)
                  else if (structured_grid%iface(ibc) == 4) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dz_loc_p(cell_id_ghosted)
                    connections%dist(3,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dy_loc_p(cell_id_ghosted)
                  else if (structured_grid%iface(ibc) == 5) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dy_loc_p(cell_id_ghosted)
                    connections%dist(2,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (structured_grid%iface(ibc) == 6) then
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

end module Structured_Grid_module
