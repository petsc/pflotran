module Structured_Grid_module

  use pflow_gridtype_module
 
  implicit none
 
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

  public :: createStructuredDMs, computeInternalConnectivity, &
            computeBoundaryConnectivity
  
contains

subroutine createStructuredDMs(grid)
      
  implicit none
  
  type(pflowGrid) :: grid
  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DA objects that will manage communication.
  !-----------------------------------------------------------------------
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,1,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_1_dof,ierr)

! call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
!      grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,3,1, &
!      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
!      grid%da_3_dof,ierr)
 
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,grid%nphase,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_nphase_dof,ierr)

  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,3*grid%nphase,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_3np_dof,ierr)

  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,grid%ndof,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_ndof,ierr)

  if (grid%use_2ph == PETSC_TRUE) then
    print *,' 2ph create DA'
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,grid%nphase*grid%npricomp,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_nphancomp_dof,ierr)

    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,grid%nphase*grid%nspec,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_nphanspec_dof,ierr)

    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz, &
                    grid%nphase*grid%nspec*grid%npricomp,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_nphanspecncomp_dof,ierr)
  endif
       
  if (grid%use_mph == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)

  endif
  
 if (grid%use_richards == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)

  endif
    
  if (grid%use_vadose == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)
  endif

  if (grid%use_flash == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)
  endif
     
  if (grid%use_owg == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)
  endif

  !print DA info for each processor
! call DAView(grid%da_ndof,PETSC_VIEWER_STDOUT_WORLD,ierr)

end subroutine createStructuredDMs

subroutine readStructuredDXYZ(grid)

  type(pflowGrid) :: grid
  
  integer :: i

  allocate(grid%dx0(grid%nx))
  allocate(grid%dy0(grid%ny))
  allocate(grid%dz0(grid%nz))
        
  call readDXYZ (grid%dx0,grid%nx)
  call readDXYZ (grid%dy0,grid%ny)
  call readDXYZ (grid%dz0,grid%nz)
    
  if (grid%myrank==0) then
    write(IUNIT2,'(/," *DXYZ ")')
    write(IUNIT2,'("  dx  ",/,(1p10e12.4))') (grid%dx0(i),i=1,grid%nx)
    write(IUNIT2,'("  dy  ",/,(1p10e12.4))') (grid%dy0(i),i=1,grid%ny)
    write(IUNIT2,'("  dz  ",/,(1p10e12.4))') (grid%dz0(i),i=1,grid%nz)
  endif

end subroutine readStructuredDXYZ

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
! computeInternalConnectivity: computes internal connectivity of a structured 
!                              grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
subroutine computeInternalConnectivity(grid)

  use Connection_module
  use pflow_gridtype_module
  
  implicit none
  
  type(pflowGrid) :: grid
  
  integer :: i, j, k, iconn, id_up, id_dn
  real*8 :: dist_up, dist_dn
  type(connection_type), pointer :: connections
  PetscErrorCode :: ierr
  
  PetscScalar, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  
  call VecGetArrayF90(grid%dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(grid%dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(grid%dz_loc, dz_loc_p, ierr)
  
  call allocateConnectionLists()
  
  connections => createConnection((grid%ngx-1)*grid%nly*grid%nlz+ &
                                  grid%nlx*(grid%ngy-1)*grid%nlz+ &
                                  grid%nlx*grid%nly*(grid%ngz-1))

  iconn = 0
  ! x-connections
  if (grid%ngx > 1) then
    do k = grid%kstart, grid%kend
      do j = grid%jstart, grid%jend
        do i = 1, grid%ngx - 1
          iconn = iconn+1
          id_up = i + j * grid%ngx + k * grid%ngxy
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
  if (grid%ngy > 1) then
    do k = grid%kstart, grid%kend
      do i = grid%istart, grid%iend
        do j = 1, grid%ngy - 1
          iconn = iconn+1
          id_up = i + 1 + (j-1) * grid%ngx + k * grid%ngxy
          id_dn = id_up + grid%ngx
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
  if (grid%ngz > 1) then
    do j = grid%jstart, grid%jend
      do i = grid%istart, grid%iend
        do k = 1, grid%ngz - 1
          iconn = iconn+1
          id_up = i + 1 + j * grid%ngx + (k-1) * grid%ngxy
          id_dn = id_up + grid%ngxy
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

  call VecRestoreArrayF90(grid%dx_loc, dx_loc_p, ierr)
  call VecRestoreArrayF90(grid%dy_loc, dy_loc_p, ierr)
  call VecRestoreArrayF90(grid%dz_loc, dz_loc_p, ierr)
  
  call addConnectionToList(connections,getInternalConnectionList())

end subroutine computeInternalConnectivity

! ************************************************************************** !
!
! computeBoundaryConnectivity: computes boundary connectivity of a structured 
!                              grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine computeBoundaryConnectivity(grid)

  use Connection_module
  use pflow_gridtype_module
  
  implicit none
  
  type(pflowGrid) :: grid
  
  integer :: num_conn_hypothetically
  integer :: i, j, k, iconn
  integer :: cell_id_local, cell_id_ghosted
  integer :: ii1, ii2, jj1, jj2, kk1, kk2
  integer :: ibc, ir
  type(connection_type), pointer :: connections
  PetscErrorCode :: ierr
  
  PetscScalar, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  
  call VecGetArrayF90(grid%dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(grid%dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(grid%dz_loc, dz_loc_p, ierr)
  
  iconn = 0
  if (grid%nx > 1 .and. grid%ny == 1 .and. grid%nz == 1) then
    if (grid%nxs == grid%ngxs) iconn = iconn + 1
    if (grid%nxe == grid%ngxe) iconn = iconn + 1
  else if (grid%nx == 1 .and. grid%ny == 1 .and. grid%nz > 1) then
    if (grid%nzs == grid%ngzs) iconn = iconn + 1
    if (grid%nze == grid%ngze) iconn = iconn + 1
  else
    if (grid%nx > 1) then
      if (grid%nxs == grid%ngxs) iconn = iconn + grid%nlyz
      if (grid%nxe == grid%ngxe) iconn = iconn + grid%nlyz
    endif
    if (grid%ny > 1) then
      if (grid%nys == grid%ngys) iconn = iconn + grid%nlxz
      if (grid%nye == grid%ngye) iconn = iconn + grid%nlxz
    endif
    if (grid%nz > 1) then
      if (grid%nzs == grid%ngzs) iconn = iconn + grid%nlxy
      if (grid%nze == grid%ngze) iconn = iconn + grid%nlxy
    endif
  endif

  num_conn_hypothetically = iconn
  connections => createConnection(iconn)

  allocate(grid%ibconn(iconn)) 

  iconn = 0
  if (grid%nxs == grid%ngxs .or. grid%nxe == grid%ngxe &
        .or. grid%nys == grid%ngys .or. grid%nye == grid%ngye &
        .or. grid%nzs == grid%ngzs .or. grid%nze == grid%ngze) then

    ! calculate boundary conditions locally on only those processors which 
    ! contain a boundary!

    do ibc = 1, grid%nblkbc
      do ir = grid%iregbc1(ibc), grid%iregbc2(ibc)
        kk1 = grid%k1bc(ir) - grid%nzs
        kk2 = grid%k2bc(ir) - grid%nzs
        jj1 = grid%j1bc(ir) - grid%nys
        jj2 = grid%j2bc(ir) - grid%nys
        ii1 = grid%i1bc(ir) - grid%nxs
        ii2 = grid%i2bc(ir) - grid%nxs

        kk1 = max(1,kk1)
        kk2 = min(grid%nlz,kk2)
        jj1 = max(1,jj1)
        jj2 = min(grid%nly,jj2)
        ii1 = max(1,ii1)
        ii2 = min(grid%nlx,ii2)

        if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle 

        do k = kk1,kk2
          do j = jj1,jj2
            do i = ii1,ii2
              iconn = iconn + 1
              cell_id_local = i+(j-1)*grid%nlx+(k-1)*grid%nlxy
!GEH - Structured Grid Dependence - End

              connections%id_dn(iconn) = cell_id_local  ! m is a local index
              
              ! old way
              !grid%ibconn(nc) = ibc
              grid%ibconn(iconn) = ir
              cell_id_ghosted = grid%nL2G(cell_id_local)
            ! Use ghosted index to access dx, dy, dz because we have
            ! already done a global-to-local scatter for computing the
            ! interior node connections.
      
!               print *,'pflowgrid_mod: ',nc,ibc,ir,m,ng,ii1,ii2,kk1,kk2, &
!                        grid%nblkbc,grid%igeom
        
              select case(grid%igeom)
                case(1) ! cartesian
                  if (grid%iface(ibc) == 1) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dx_loc_p(cell_id_ghosted)
                    connections%dist(1,iconn) = 1.d0
                    connections%area(iconn) = dy_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (grid%iface(ibc) == 2) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dx_loc_p(cell_id_ghosted)
                    connections%dist(1,iconn) = 1.d0
                    connections%area(iconn) = dy_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (grid%iface(ibc) == 3) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dz_loc_p(cell_id_ghosted)
                    connections%dist(3,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dy_loc_p(cell_id_ghosted)
                  else if (grid%iface(ibc) == 4) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dz_loc_p(cell_id_ghosted)
                    connections%dist(3,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dy_loc_p(cell_id_ghosted)
                  else if (grid%iface(ibc) == 5) then
                    connections%dist(:,iconn) = 0.d0
                    connections%dist(0,iconn) = 0.5d0*dy_loc_p(cell_id_ghosted)
                    connections%dist(2,iconn) = 1.d0
                    connections%area(iconn) = dx_loc_p(cell_id_ghosted)* &
                                              dz_loc_p(cell_id_ghosted)
                  else if (grid%iface(ibc) == 6) then
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

  call VecRestoreArrayF90(grid%dx_loc, dx_loc_p, ierr)
  call VecRestoreArrayF90(grid%dy_loc, dy_loc_p, ierr)
  call VecRestoreArrayF90(grid%dz_loc, dz_loc_p, ierr)
  
  call addConnectionToList(connections,getBoundaryConnectionList())

end subroutine computeBoundaryConnectivity

end module Structured_Grid_module
