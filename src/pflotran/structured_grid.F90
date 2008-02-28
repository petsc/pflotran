module Structured_Grid_module

  implicit none
 
  private
 
#include "definitions.h"

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

    PetscReal :: origin(3)

    PetscReal, pointer :: dx0(:), dy0(:), dz0(:)
    
    logical :: invert_z_axis
    
    PetscReal, pointer :: dx(:), dy(:), dz(:), dxg(:), dyg(:), dzg(:)  ! Grid spacings
    
  end type structured_grid_type

  public :: StructuredGridCreate, &
            StructuredGridDestroy, &
            StructuredGridCreateDA, &
            StructGridComputeInternConnect, &
            StructuredGridCreateVecFromDA, &
            StructuredGridMapIndices, &
            StructuredGridComputeSpacing, &
            StructuredGridComputeCoord, &
            StructuredGridReadDXYZ, &
            StructuredGridComputeVolumes, &
            StructGridPopulateConnection

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
  nullify(structured_grid%dx)
  nullify(structured_grid%dy)
  nullify(structured_grid%dz)
  nullify(structured_grid%dxg)
  nullify(structured_grid%dyg)
  nullify(structured_grid%dzg)
  
  
  structured_grid%origin = 0.d0
  
  structured_grid%invert_z_axis = .false.
  
  StructuredGridCreate => structured_grid
  
end function StructuredGridCreate
  
! ************************************************************************** !
!
! StructuredGridCreateDMs: Creates structured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine StructuredGridCreateDA(structured_grid,da,ndof,stencil_width)
      
  use Option_module
      
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"

  type(structured_grid_type) :: structured_grid
  DA :: da
  PetscInt :: ndof
  PetscInt :: stencil_width

  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DA object that will manage communication.
  !-----------------------------------------------------------------------
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                  structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                  structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                  ndof,stencil_width, &
                  PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                  da,ierr)

  if (structured_grid%nlx+structured_grid%nly+structured_grid%nlz == 0) then

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
  endif
  
end subroutine StructuredGridCreateDA

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
subroutine StructuredGridReadDXYZ(structured_grid,option)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  
  PetscInt :: i

  allocate(structured_grid%dx0(structured_grid%nx))
  structured_grid%dx0 = 0.d0
  allocate(structured_grid%dy0(structured_grid%ny))
  structured_grid%dy0 = 0.d0
  allocate(structured_grid%dz0(structured_grid%nz))
  structured_grid%dz0 = 0.d0

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
subroutine StructuredGridComputeSpacing(structured_grid,nG2A,nG2L)

  implicit none
  
  type(structured_grid_type) :: structured_grid
  PetscInt :: nG2A(:)
  PetscInt :: nG2L(:)
  
  PetscInt :: i, j, k, nl, ng, na
  PetscErrorCode :: ierr

  allocate(structured_grid%dx(structured_grid%nlmax))
  structured_grid%dx = 0.d0
  allocate(structured_grid%dy(structured_grid%nlmax))
  structured_grid%dy = 0.d0
  allocate(structured_grid%dz(structured_grid%nlmax))
  structured_grid%dz = 0.d0
  
  allocate(structured_grid%dxg(structured_grid%ngmax))
  structured_grid%dxg = 0.d0
  allocate(structured_grid%dyg(structured_grid%ngmax))
  structured_grid%dyg = 0.d0
  allocate(structured_grid%dzg(structured_grid%ngmax))
  structured_grid%dzg = 0.d0
        
  do ng = 1,structured_grid%ngmax
    na = nG2A(ng)
    nl = nG2L(ng)
    k= int(na/structured_grid%ngxy) + 1
    j= int(mod(na,structured_grid%ngxy)/structured_grid%ngx) + 1
    i= mod(mod(na,structured_grid%ngxy),structured_grid%ngx) + 1
    structured_grid%dxg(ng) = structured_grid%dx0(i)
    structured_grid%dyg(ng) = structured_grid%dy0(j)
    structured_grid%dzg(ng) = structured_grid%dz0(k)
    if (nl > 0) then
      structured_grid%dxg(nl) = structured_grid%dx0(i)
      structured_grid%dyg(nl) = structured_grid%dy0(j)
      structured_grid%dzg(nl) = structured_grid%dz0(k)
    endif
  enddo
  
end subroutine StructuredGridComputeSpacing

! ************************************************************************** !
!
! StructuredGridComputeCoord: Computes structured coordinates in x,y,z
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructuredGridComputeCoord(structured_grid,option, &
                                      grid_x,grid_y,grid_z, &
                                      x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

! PetscInt :: ierr
  PetscInt :: i, j, k, n
  PetscReal :: x, y, z
  PetscInt :: prevnode

! set min and max bounds of domain in coordinate directions
  x_min = structured_grid%origin(X_DIRECTION)
  y_min = structured_grid%origin(Y_DIRECTION)
  z_min = structured_grid%origin(Z_DIRECTION)
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
  z = 0.5d0*structured_grid%dz0(1)+structured_grid%origin(Z_DIRECTION)
  do k=1, structured_grid%nz
    y = 0.5d0*structured_grid%dy0(1)+structured_grid%origin(Y_DIRECTION)
    do j=1, structured_grid%ny
      x = 0.5d0*structured_grid%dx0(1)+structured_grid%origin(X_DIRECTION)
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
          dist_up = 0.5d0*structured_grid%dxg(id_up)
          dist_dn = 0.5d0*structured_grid%dxg(id_dn)
          connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
          connections%dist(0,iconn) = dist_up+dist_dn
          connections%dist(1,iconn) = 1.d0  ! x component of unit vector
          connections%area(iconn) = structured_grid%dyg(id_up)*structured_grid%dzg(id_up)
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
          dist_up = 0.5d0*structured_grid%dyg(id_up)
          dist_dn = 0.5d0*structured_grid%dyg(id_dn)
          connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
          connections%dist(0,iconn) = dist_up+dist_dn
          connections%dist(2,iconn) = 1.d0  ! y component of unit vector
          connections%area(iconn) = structured_grid%dxg(id_up)*structured_grid%dzg(id_up)
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
          dist_up = 0.5d0*structured_grid%dzg(id_up)
          dist_dn = 0.5d0*structured_grid%dzg(id_dn)
          connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
          connections%dist(0,iconn) = dist_up+dist_dn
          connections%dist(3,iconn) = 1.d0  ! z component of unit vector
          connections%area(iconn) = structured_grid%dxg(id_up)*structured_grid%dyg(id_up)
        enddo
      enddo
    enddo
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
subroutine StructGridPopulateConnection(structured_grid,connection,iface, &
                                        iconn,cell_id_ghosted)

  use Connection_module
  
  implicit none
 
  type(structured_grid_type) :: structured_grid
  type(connection_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: cell_id_ghosted
  
  PetscErrorCode :: ierr
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
      select case(iface)
        case(WEST_FACE,EAST_FACE)
          connection%dist(:,iconn) = 0.d0
          connection%dist(0,iconn) = 0.5d0*structured_grid%dxg(cell_id_ghosted)
          connection%area(iconn) = structured_grid%dyg(cell_id_ghosted)* &
                                   structured_grid%dzg(cell_id_ghosted)
          if (iface ==  WEST_FACE) then
            connection%dist(1,iconn) = 1.d0
          else
            connection%dist(1,iconn) = -1.d0
          endif
        case(SOUTH_FACE,NORTH_FACE)
          connection%dist(:,iconn) = 0.d0
          connection%dist(0,iconn) = 0.5d0*structured_grid%dyg(cell_id_ghosted)
          connection%area(iconn) = structured_grid%dxg(cell_id_ghosted)* &
                                   structured_grid%dzg(cell_id_ghosted)
          if (iface ==  SOUTH_FACE) then
            connection%dist(2,iconn) = 1.d0
          else
            connection%dist(2,iconn) = -1.d0
          endif
        case(BOTTOM_FACE,TOP_FACE)
          connection%dist(:,iconn) = 0.d0
          connection%dist(0,iconn) = 0.5d0*structured_grid%dzg(cell_id_ghosted)
          connection%area(iconn) = structured_grid%dxg(cell_id_ghosted)* &
                                   structured_grid%dyg(cell_id_ghosted)
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
subroutine StructuredGridComputeVolumes(structured_grid,option,nL2G,volume)

  use Option_module
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: nL2G(:)
  Vec :: volume
  
  PetscReal, parameter :: Pi=3.1415926d0
  
  PetscInt :: i, n, ng
  PetscReal, pointer :: volume_p(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(volume,volume_p, ierr)
  do n=1, structured_grid%nlmax
    ng = nL2G(n)
    volume_p(n) = structured_grid%dxg(ng) * structured_grid%dyg(ng) * &
                  structured_grid%dzg(ng)
  enddo
  call VecRestoreArrayF90(volume,volume_p, ierr)
  
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
subroutine StructuredGridMapIndices(structured_grid,nG2L,nL2G,nL2A,nG2A)

  implicit none

  type(structured_grid_type) :: structured_grid
  PetscInt, pointer :: nG2L(:), nL2G(:), nL2A(:), nG2A(:)

  PetscInt :: i, j, k, n, ng, na, count1
  PetscErrorCode :: ierr
  
  allocate(nL2G(structured_grid%nlmax))
  allocate(nG2L(structured_grid%ngmax))
  allocate(nL2A(structured_grid%nlmax))
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
        na = i-1+structured_grid%nxs+(j-1+structured_grid%nys)*structured_grid%nx+ &
             (k-1+structured_grid%nzs)*structured_grid%nxy
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
        na = i-1+structured_grid%ngxs+(j-1+structured_grid%ngys)*structured_grid%nx+ &
             (k-1+structured_grid%ngzs)*structured_grid%nxy
        if (na>(structured_grid%nmax-1)) print *,'Wrong Nature order....'
        nG2A(n) = na
!        print *,grid%myrank, k,j,i,n,na
        !grid%nG2N(ng) = na
      enddo
    enddo
  enddo
 
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
  

  if (associated(structured_grid%dx0)) deallocate(structured_grid%dx0)
  nullify(structured_grid%dx0)
  if (associated(structured_grid%dy0)) deallocate(structured_grid%dy0)
  nullify(structured_grid%dy0)
  if (associated(structured_grid%dz0)) deallocate(structured_grid%dz0)
  nullify(structured_grid%dz0)

  if (associated(structured_grid%dx)) deallocate(structured_grid%dx)
  nullify(structured_grid%dx)
  if (associated(structured_grid%dy)) deallocate(structured_grid%dy)
  nullify(structured_grid%dy)
  if (associated(structured_grid%dz)) deallocate(structured_grid%dz)
  nullify(structured_grid%dz)
  if (associated(structured_grid%dxg)) deallocate(structured_grid%dxg)
  nullify(structured_grid%dxg)
  if (associated(structured_grid%dyg)) deallocate(structured_grid%dyg)
  nullify(structured_grid%dyg)
  if (associated(structured_grid%dzg)) deallocate(structured_grid%dzg)
  nullify(structured_grid%dzg)
  
  deallocate(structured_grid)
  nullify(structured_grid)

end subroutine StructuredGridDestroy
                          
end module Structured_Grid_module
