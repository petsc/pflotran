module Structured_Grid_module

  use PFLOTRAN_Constants_module

  implicit none
 
  private
 
#include "finclude/petscsys.h"

! structured grid faces
  PetscInt, parameter, public :: NULL_FACE = 0
  PetscInt, parameter, public :: WEST_FACE = 1
  PetscInt, parameter, public :: EAST_FACE = 2
  PetscInt, parameter, public :: SOUTH_FACE = 3
  PetscInt, parameter, public :: NORTH_FACE = 4
  PetscInt, parameter, public :: BOTTOM_FACE = 5
  PetscInt, parameter, public :: TOP_FACE = 6

  PetscInt, parameter, public :: CARTESIAN_GRID = 3
  PetscInt, parameter, public :: CYLINDRICAL_GRID = 4
  PetscInt, parameter, public :: SPHERICAL_GRID = 5

  type, public :: structured_grid_type

    character(len=MAXWORDLENGTH) :: ctype
    PetscInt :: itype  ! type of grid (e.g. structured, unstructured, etc.)
    PetscInt :: global_offset
    PetscInt :: nx, ny, nz    ! Global domain dimensions of the grid.
    PetscInt :: nxy, nmax     ! nx * ny, nx * ny * nz
    PetscInt :: npx, npy, npz ! Processor partition in each direction.
    PetscInt :: npx_final, npy_final, npz_final ! actual decomposition
    PetscInt :: nlx, nly, nlz ! Local grid dimension w/o ghost nodes.
    PetscInt :: ngx, ngy, ngz ! Local grid dimension with ghost nodes.
    PetscInt :: lxs, lys, lzs 
      ! Global indices of non-ghosted corner (starting) of local domain.
    PetscInt :: gxs, gys, gzs
      ! Global indices of ghosted starting corner of local domain.
    PetscInt :: lxe, lye, lze, gxe, gye, gze
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
    
    PetscInt, pointer :: cell_neighbors(:,:)
    
  end type structured_grid_type

  public :: StructGridCreate, &
            StructGridDestroy, &
            StructGridCreateDM, &
            StructGridComputeLocalBounds, &
            StructGridComputeInternConnect, &
            StructGridComputeBoundConnect, &
            StructGridCreateVecFromDM, &
            StructGridMapIndices, &
            StructGridComputeSpacing, &
            StructGridComputeCoord, &
            StructGridReadDXYZ, &
            StructGridComputeVolumes, &
            StructGridPopulateConnection, &
            StructGridGetIJKFromCoordinate, &
            StructGridGetIJKFromLocalID, &
            StructGridGetIJKFromGhostedID, &
            StructGridGetGhostedNeighbors, &
            StructGridCreateTVDGhosts, &
            StructGridGetGhostedNeighborsCorners, &
            StructGridComputeNeighbors
  
contains

! ************************************************************************** !
!
! StructGridCreate: Creates a structured grid object
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
function StructGridCreate()

  implicit none
  
  type(structured_grid_type), pointer :: StructGridCreate

  type(structured_grid_type), pointer :: structured_grid

  allocate(structured_grid)
  
  structured_grid%ctype = ''
  structured_grid%itype = 0
  structured_grid%global_offset = 0
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

  structured_grid%lxs = 0
  structured_grid%lys = 0
  structured_grid%lzs = 0

  structured_grid%gxs = 0
  structured_grid%gys = 0
  structured_grid%gzs = 0

  structured_grid%lxe = 0
  structured_grid%lye = 0
  structured_grid%lze = 0

  structured_grid%gxe = 0
  structured_grid%gye = 0
  structured_grid%gze = 0

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
  
  nullify(structured_grid%cell_neighbors)
 
  structured_grid%origin = -1.d20
  structured_grid%bounds = -1.d20
  
  structured_grid%invert_z_axis = PETSC_FALSE
  
  StructGridCreate => structured_grid
  
end function StructGridCreate
  
! ************************************************************************** !
!
! StructGridCreateDMs: Creates structured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine StructGridCreateDM(structured_grid,da,ndof,stencil_width, &
                              stencil_type,option)

  use Option_module
        
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"

  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  DM :: da
  PetscInt :: ndof
  PetscInt :: stencil_width,stencil_type

  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DM object that will manage communication.
  !-----------------------------------------------------------------------
  ! This code is for the DMDACreate3D() interface in PETSc versions >= 3.2 --RTM
  call DMDACreate3D(option%mycomm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, &
                  DMDA_BOUNDARY_NONE,stencil_type, &
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

end subroutine StructGridCreateDM

! ************************************************************************** !
!
! StructGridComputeLocalBounds: Computes the corners for the local portion
!                               of the structured grid
! author: Glenn Hammond
! date: 03/13/08
!
! ************************************************************************** !
subroutine StructGridComputeLocalBounds(structured_grid,da,option)

  use Option_module
  
  implicit none
     
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"

  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  DM :: da

  PetscErrorCode :: ierr

  ! get corner information
  call DMDAGetCorners(da, structured_grid%lxs, &
      structured_grid%lys, structured_grid%lzs, structured_grid%nlx, &
      structured_grid%nly, structured_grid%nlz, ierr)
     
  structured_grid%lxe = structured_grid%lxs + structured_grid%nlx
  structured_grid%lye = structured_grid%lys + structured_grid%nly
  structured_grid%lze = structured_grid%lzs + structured_grid%nlz
  structured_grid%nlxy = structured_grid%nlx * structured_grid%nly
  structured_grid%nlxz = structured_grid%nlx * structured_grid%nlz
  structured_grid%nlyz = structured_grid%nly * structured_grid%nlz
  structured_grid%nlmax = structured_grid%nlx * structured_grid%nly * structured_grid%nlz
     
  ! get ghosted corner information
  call DMDAGetGhostCorners(da, structured_grid%gxs, &
      structured_grid%gys, structured_grid%gzs, structured_grid%ngx, &
      structured_grid%ngy, structured_grid%ngz, ierr)
     
  structured_grid%gxe = structured_grid%gxs + structured_grid%ngx
  structured_grid%gye = structured_grid%gys + structured_grid%ngy
  structured_grid%gze = structured_grid%gzs + structured_grid%ngz
  structured_grid%ngxy = structured_grid%ngx * structured_grid%ngy
  structured_grid%ngxz = structured_grid%ngx * structured_grid%ngz
  structured_grid%ngyz = structured_grid%ngy * structured_grid%ngz
  structured_grid%ngmax = structured_grid%ngx * structured_grid%ngy * structured_grid%ngz
  
  structured_grid%global_offset = 0
  call MPI_Exscan(structured_grid%nlmax,structured_grid%global_offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)  
   
end subroutine StructGridComputeLocalBounds

! ************************************************************************** !
!
! StructGridCreateVecFromDM: Creates a PETSc vector from a DM
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine StructGridCreateVecFromDM(da,vector,vector_type)

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

end subroutine StructGridCreateVecFromDM

! ************************************************************************** !
!
! StructGridReadDXYZ: Reads structured grid spacing from input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine StructGridReadDXYZ(structured_grid,input,option)

  use Option_module
  use Input_Aux_module
  
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
  call StructGridReadArrayNew(structured_grid%dx_global, &
                               structured_grid%nx,word,input,option)
  word = 'Y'
  call StructGridReadArrayNew(structured_grid%dy_global, &
                               structured_grid%ny,word,input,option)
  word = 'Z'
  call StructGridReadArrayNew(structured_grid%dz_global, &
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

end subroutine StructGridReadDXYZ

! ************************************************************************** !
!
! StructGridReadArray: Reads structured grid spacing along an axis from  
!                         input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine StructGridReadArray(a,n,input,option)

  use Input_Aux_module
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
    call InputReadPflotranString(input,option)
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
      if (a(i) == 0.d0) then

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
    
end subroutine StructGridReadArray

! ************************************************************************** !
!
! StructGridReadArrayNew: Reads structured grid spacing along an axis from  
!                         input file
! author: Glenn Hammond
! date: 05/21/09
!
! ************************************************************************** !
subroutine StructGridReadArrayNew(array,array_size,axis,input,option)

  use Input_Aux_module
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
    
    call InputReadPflotranString(input,option)
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
        call InputErrorMsg(input,option,'# values','StructGridReadArrayNew')
        string2 = word3
        call InputReadDouble(string2,option,value,input%ierr)
        call InputErrorMsg(input,option,'value','StructGridReadArrayNew')
        do i=1, num_values
          count = count + 1
          array(count) = value
        enddo
      else
        string2 = word
        call InputReadDouble(string2,option,value,input%ierr)
        call InputErrorMsg(input,option,'value','StructGridReadDXYZ')
        count = count + 1
        array(count) = value
      endif
    enddo
  enddo

end subroutine StructGridReadArrayNew

! ************************************************************************** !
!
! StructGridComputeSpacing: Computes structured grid spacing
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine StructGridComputeSpacing(structured_grid,option)

  use Option_module
  
  implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  
  PetscInt :: i, j, k, ghosted_id
  PetscErrorCode :: ierr

  allocate(structured_grid%dxg_local(structured_grid%ngx))
  structured_grid%dxg_local = 0.d0
  allocate(structured_grid%dyg_local(structured_grid%ngy))
  structured_grid%dyg_local = 0.d0
  allocate(structured_grid%dzg_local(structured_grid%ngz))
  structured_grid%dzg_local = 0.d0
  
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
    structured_grid%dx_global(structured_grid%gxs+1:structured_grid%gxe)
  structured_grid%dyg_local(1:structured_grid%ngy) = &
    structured_grid%dy_global(structured_grid%gys+1:structured_grid%gye)
  structured_grid%dzg_local(1:structured_grid%ngz) = &
    structured_grid%dz_global(structured_grid%gzs+1:structured_grid%gze)
        
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
  
end subroutine StructGridComputeSpacing

! ************************************************************************** !
!
! StructGridComputeCoord: Computes structured coordinates in x,y,z
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructGridComputeCoord(structured_grid,option,origin_global, &
                                      grid_x,grid_y,grid_z, &
                                      x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module
  
implicit none
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscReal :: origin_global(3)
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

! PetscErrorCode :: ierr
  PetscInt :: i, j, k, ghosted_id
  PetscReal :: x, y, z
  PetscInt :: prevnode

  x_min = origin_global(X_DIRECTION)
  y_min = origin_global(Y_DIRECTION)
  z_min = origin_global(Z_DIRECTION)
    
  do i=1,structured_grid%lxs
    x_min = x_min + structured_grid%dx_global(i)
  enddo
  do j=1,structured_grid%lys
    y_min = y_min + structured_grid%dy_global(j)
  enddo
  do k=1,structured_grid%lzs
    z_min = z_min + structured_grid%dz_global(k)
  enddo
   
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
  if (structured_grid%kstart > 0) then
    z = -0.5d0*structured_grid%dzg_local(1)+structured_grid%origin(Z_DIRECTION)
  else
    z = 0.5d0*structured_grid%dzg_local(1)+structured_grid%origin(Z_DIRECTION)
  endif
  do k=1, structured_grid%ngz
    if (structured_grid%jstart > 0) then
      y = -0.5d0*structured_grid%dyg_local(1)+structured_grid%origin(Y_DIRECTION)
    else
      y = 0.5d0*structured_grid%dyg_local(1)+structured_grid%origin(Y_DIRECTION)
    endif
    do j=1, structured_grid%ngy
      if (structured_grid%istart > 0) then
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
    
end subroutine StructGridComputeCoord

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
          if (x == x_upper_face .and. &
            structured_grid%lxs /= structured_grid%gxs) exit
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
        if (y == y_upper_face .and. &
            structured_grid%lys /= structured_grid%gys) exit
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
        if (z == z_upper_face .and. &
          structured_grid%lzs /= structured_grid%gzs) exit
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

!  PetscReal :: radius(:)
  type(connection_set_type), pointer :: StructGridComputeInternConnect
  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  PetscReal, pointer :: xc(:),yc(:),zc(:)
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  PetscInt :: i, j, k, iconn, id_up, id_dn, id_up2, id_dn2
  PetscInt :: nconn
  PetscInt :: lenx, leny, lenz
  PetscInt :: tvd_ghost_offset, ghost_count
  PetscReal :: dist_up, dist_dn
  PetscReal :: r1, r2
  type(connection_set_type), pointer :: connections, connections_2

  PetscErrorCode :: ierr
  PetscReal, pointer :: radius(:)
  PetscInt, pointer :: int_array1(:), int_array2(:),int_array3(:),int_array4(:),int_array5(:),index(:)
  PetscInt :: count

  radius => xc

  ! the adjustments in the case of AMR are based on the PIMS code adjustments by LC
  nconn = (structured_grid%ngx-1)*structured_grid%nly*structured_grid%nlz+ &
          structured_grid%nlx*(structured_grid%ngy-1)*structured_grid%nlz+ &
          structured_grid%nlx*structured_grid%nly*(structured_grid%ngz-1)
  
  structured_grid%nlmax_faces = 0
  structured_grid%ngmax_faces = 0

  lenx = structured_grid%ngx - 1
  leny = structured_grid%ngy - 1
  lenz = structured_grid%ngz - 1

  connections => ConnectionCreate(nconn,INTERNAL_CONNECTION_TYPE)
  
  ! if using higher order advection, allocate associated arrays
  if (option%itranmode == EXPLICIT_ADVECTION .and. &
      option%tvd_flux_limiter /= 1) then  ! 1 = upwind
    allocate(connections%id_up2(size(connections%id_up)))
    allocate(connections%id_dn2(size(connections%id_dn)))
    connections%id_up2 = 0
    connections%id_dn2 = 0
  endif

  iconn = 0
  tvd_ghost_offset = 0
  ghost_count = 0
  ! x-connections
  if (structured_grid%ngx > 1) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy
              id_dn = id_up + 1
#ifdef DASVYAT
              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              if (i == 1) then
                if (structured_grid%lxs==structured_grid%gxs) then
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
                end if
              else 
                structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                connections%local(iconn) = 1
              end if
#endif
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              
              if (associated(connections%id_up2)) then
                if (i == 1) then
                  ! id_up indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_up2 = 1 + j + k*structured_grid%nly
                  ghost_count = ghost_count + 1
                  id_up2 = ghost_count
                  connections%id_up2(iconn) = -id_up2
                else
                  connections%id_up2(iconn) = id_up - 1
                endif
                if (i == lenx) then
                  ! id_dn indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_dn2 = 1 + j + k*structured_grid%nly + structured_grid%nlyz
                  id_dn2 = ghost_count + structured_grid%nlyz
                  connections%id_dn2(iconn) = -id_dn2
                else
                  connections%id_dn2(iconn) = id_dn + 1
                endif
              endif
              
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dx(id_up)
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = structured_grid%dy(id_up)* &
                                        structured_grid%dz(id_up)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_up) + &
                  connections%dist(-1,iconn)*(xc(id_dn) - xc(id_up))
              connections%cntr(2,iconn) = yc(id_up)
              connections%cntr(3,iconn) = zc(id_up)
#endif              
            enddo
          enddo
        enddo
        tvd_ghost_offset = 2*structured_grid%nlyz ! west & east
        ghost_count = tvd_ghost_offset
      case(CYLINDRICAL_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy
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
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_up) + connections%dist(-1,iconn)*(xc(id_dn) - xc(id_up))
              connections%cntr(2,iconn) = yc(id_up)
              connections%cntr(3,iconn) = zc(id_up)
#endif
            enddo
          enddo
        enddo
      case(SPHERICAL_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy
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
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_up) + connections%dist(-1,iconn)*(xc(id_dn) - xc(id_up))
              connections%cntr(2,iconn) = yc(id_up)
              connections%cntr(3,iconn) = zc(id_up)
#endif
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

#ifdef DASVYAT
              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              if (j == 1) then
                if (structured_grid%lys==structured_grid%gys) then
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
                end if
              else 
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
              end if
#endif
              id_up = i + 1 + (j-1) * structured_grid%ngx + k * structured_grid%ngxy
              id_dn = id_up + structured_grid%ngx
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              
              if (associated(connections%id_up2)) then
                if (j == 1) then
                  ! id_up indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_up2 = 1 + i + k*structured_grid%nlx + tvd_ghost_offset
                  ghost_count = ghost_count + 1
                  id_up2 = ghost_count
                  connections%id_up2(iconn) = -id_up2
                else
                  connections%id_up2(iconn) = id_up - structured_grid%ngx
                endif
                if (j == leny) then
                  ! id_dn indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_dn2 = 1 + i + k*structured_grid%nlx + &
!                           structured_grid%nlxz + tvd_ghost_offset
                  id_dn2 = ghost_count + structured_grid%nlxz
                  connections%id_dn2(iconn) = -id_dn2
                else
                  connections%id_dn2(iconn) = id_dn + structured_grid%ngx
                endif
              endif
                            
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dy(id_up)
              dist_dn = 0.5d0*structured_grid%dy(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(2,iconn) = 1.d0  ! y component of unit vector
              connections%area(iconn) = structured_grid%dx(id_up)* &
                                    structured_grid%dz(id_up)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_up) 
              connections%cntr(2,iconn) = yc(id_up) + connections%dist(-1,iconn)*(yc(id_dn) - yc(id_up)) 
              connections%cntr(3,iconn) = zc(id_up)
#endif
            enddo
          enddo
        enddo
        tvd_ghost_offset = tvd_ghost_offset + &
          2*structured_grid%nlxz ! south & north
        ghost_count = tvd_ghost_offset
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

#ifdef DASVYAT
              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              if (k == 1) then
                if (structured_grid%lzs==structured_grid%gzs) then
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
                end if
              else 
                   structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1
                   connections%local(iconn) = 1
              end if
#endif
              id_up = i + 1 + j * structured_grid%ngx + (k-1) * &
                  structured_grid%ngxy 
              id_dn = id_up + structured_grid%ngxy
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              
              if (associated(connections%id_up2)) then
                if (k == 1) then
!                  id_up2 = 1 + i + j*structured_grid%nlx + tvd_ghost_offset
                  ghost_count = ghost_count + 1
                  id_up2 = ghost_count
                  connections%id_up2(iconn) = -id_up2
                else
                  connections%id_up2(iconn) = id_up - structured_grid%ngxy
                endif
                if (k == lenz) then
                  ! id_dn indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_dn2 = 1 + i + j*structured_grid%nlx + &
!                           structured_grid%nlxy + tvd_ghost_offset
                  id_dn2 = ghost_count + structured_grid%nlxy
                  connections%id_dn2(iconn) = -id_dn2
                else
                  connections%id_dn2(iconn) = id_dn + structured_grid%ngxy
                endif
              endif
                                 
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dz(id_up)
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(3,iconn) = 1.d0  ! z component of unit vector
              connections%area(iconn) = structured_grid%dx(id_up) * &
                                        structured_grid%dy(id_up)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_up) 
              connections%cntr(2,iconn) = yc(id_up) 
              connections%cntr(3,iconn) = zc(id_up) + connections%dist(-1,iconn)*(zc(id_dn) - zc(id_up)) 
#endif
            enddo
          enddo
        enddo
      case(CYLINDRICAL_GRID)
        do j = structured_grid%jstart, structured_grid%jend
          do i = structured_grid%istart, structured_grid%iend
            do k = 1, lenz
              iconn = iconn+1
              id_up = i + 1 + j * structured_grid%ngx + (k-1) * &
                  structured_grid%ngxy
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
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_up) 
              connections%cntr(2,iconn) = yc(id_up) 
              connections%cntr(3,iconn) = zc(id_up) + &
                                          connections%dist(-1,iconn)* &
                                          (zc(id_dn) - zc(id_up)) 
#endif
            enddo
          enddo
        enddo
      case(SPHERICAL_GRID)
        option%io_buffer = 'For spherical coordinates, NZ must be equal to 1.'
        call printErrMsg(option)
  end select
  endif

  StructGridComputeInternConnect => connections

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

  PetscReal, pointer :: xc(:), yc(:), zc(:)
  type(connection_set_type), pointer :: StructGridComputeBoundConnect
  type(option_type) :: option
  type(structured_grid_type) :: structured_grid
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  PetscInt :: i, j, k, iconn, id_up, id_dn, id_dn2
  PetscInt :: nconn
  PetscInt :: lenx, leny, lenz
  PetscReal :: dist_up, dist_dn
  type(connection_set_type), pointer :: connections
  PetscErrorCode :: ierr
  PetscReal, pointer :: radius(:)

  radius => xc
  
  ! the adjustments in the case of AMR are based on the PIMS code adjustments by LC

  nconn = structured_grid%nly*structured_grid%nlz*(structured_grid%nlx + 2 - structured_grid%ngx)&
         +structured_grid%nlx*structured_grid%nlz*(structured_grid%nly + 2 - structured_grid%ngy)&  
         +structured_grid%nlx*structured_grid%nly*(structured_grid%nlz + 2 - structured_grid%ngz)

  lenx = structured_grid%nlx 
  leny = structured_grid%nly 
  lenz = structured_grid%nlz 


!  write(*,*) option%myrank, 'lxs=',structured_grid%lxs,' lys=',structured_grid%lys, ' lzs=',structured_grid%lzs 
!  write(*,*) option%myrank, 'lxe=',structured_grid%lxe,' lye=',structured_grid%lye, ' lze=',structured_grid%lze
!  write(*,*) option%myrank, 'nx=',structured_grid%nlx,' ny=',structured_grid%ny, ' nz=',structured_grid%nz  
!  write(*,*) option%myrank, 'istart=',structured_grid%istart,' iend=',structured_grid%iend
!  write(*,*) option%myrank, 'jstart=',structured_grid%jstart,' jend=',structured_grid%jend
!  write(*,*) option%myrank, 'kstart=',structured_grid%kstart,' kend=',structured_grid%kend
!  write(*,*) option%myrank, 'boundary nconn=',nconn

!  stop


  connections => ConnectionCreate(nconn,BOUNDARY_CONNECTION_TYPE)

!  StructGridComputeBoundConnect => connections

  iconn = 0
  ! x-connections
  if (structured_grid%nlx + 2 - structured_grid%ngx > 0) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        if (structured_grid%lxs == structured_grid%gxs) then
          i = structured_grid%istart 
          do k = structured_grid%kstart, structured_grid%kend
            do j = structured_grid%jstart, structured_grid%jend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy
              
              connections%id_dn(iconn) = id_dn

              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = 0.
              connections%dist(0,iconn) = dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = structured_grid%dy(id_dn)* &
                                        structured_grid%dz(id_dn)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_dn) - dist_dn
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) 
#endif
            enddo
          enddo
        endif
        if (structured_grid%lxe == structured_grid%gxe) then
          i = structured_grid%iend 
          do k = structured_grid%kstart, structured_grid%kend
            do j = structured_grid%jstart, structured_grid%jend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy
              
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(1,iconn) = -1.d0  ! x component of unit vector
              connections%area(iconn) = structured_grid%dy(id_dn)* &
                                        structured_grid%dz(id_dn)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_dn) + dist_dn
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) 
#endif
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
        if (structured_grid%lys == structured_grid%gys) then 
          j = structured_grid%jstart
          do k = structured_grid%kstart, structured_grid%kend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy
              
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dy(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(2,iconn) = 1.d0  ! y component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn)* &
                                    structured_grid%dz(id_dn)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) - dist_dn
              connections%cntr(3,iconn) = zc(id_dn) 
#endif
            enddo
          enddo
        endif
        if (structured_grid%lye == structured_grid%gye) then 
          j = structured_grid%jend
          do k = structured_grid%kstart, structured_grid%kend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * structured_grid%ngxy
              
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dy(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(2,iconn) = -1.d0  ! y component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn)* &
                                    structured_grid%dz(id_dn)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) + dist_dn
              connections%cntr(3,iconn) = zc(id_dn) 
#endif
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
        if (structured_grid%lzs == structured_grid%gzs) then
          k = structured_grid%kstart
          do j = structured_grid%jstart, structured_grid%jend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * &
                  structured_grid%ngxy
              
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(3,iconn) = 1.d0  ! z component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn) * &
                                        structured_grid%dy(id_dn)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) - dist_dn 
#endif
            enddo
          enddo
        endif    
        if (structured_grid%lze == structured_grid%gze ) then
          k = structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = structured_grid%istart, structured_grid%iend
              iconn = iconn+1

              structured_grid%ngmax_faces = structured_grid%ngmax_faces + 1
              structured_grid%nlmax_faces = structured_grid%nlmax_faces + 1

              id_dn = i + 1 + j * structured_grid%ngx + k * &
                  structured_grid%ngxy
              
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(0,iconn) = dist_dn
              connections%dist(3,iconn) = -1.d0  ! z component of unit vector
              connections%area(iconn) = structured_grid%dx(id_dn) * &
                                        structured_grid%dy(id_dn)
#ifdef DASVYAT
              connections%cntr(1,iconn) = xc(id_dn) 
              connections%cntr(2,iconn) = yc(id_dn) 
              connections%cntr(3,iconn) = zc(id_dn) + dist_dn 
#endif
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
                                        iconn,ghosted_id,option)

  use Connection_module
  use Option_module
  
  implicit none
 
  type(structured_grid_type) :: structured_grid
  type(connection_set_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: ghosted_id
  type(option_type) :: option
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
              if (associated(connection%id_dn2)) then
                  connection%id_dn2(iconn) = &
                    StructGetTVDGhostConnection(ghosted_id, &
                                                structured_grid, &
                                                iface,option)
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
              if (iface == SOUTH_FACE) then
                connection%dist(2,iconn) = 1.d0
              else
                connection%dist(2,iconn) = -1.d0
              endif
              if (associated(connection%id_dn2)) then
                  connection%id_dn2(iconn) = &
                    StructGetTVDGhostConnection(ghosted_id, &
                                                structured_grid, &
                                                iface,option)
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
                if (iface == TOP_FACE) then 
                  option%io_buffer = 'Need to ensure that direction of ' // &
                    'inverted z is correct in StructGridPopulateConnection()'
                  call printErrMsg(option)
                  connection%dist(3,iconn) = -1.d0
                else
                  connection%dist(3,iconn) = 1.d0
                endif
              else
                if (iface == BOTTOM_FACE) then 
                  connection%dist(3,iconn) = 1.d0
                else
                  connection%dist(3,iconn) = -1.d0
                endif
              endif
              if (associated(connection%id_dn2)) then
                  connection%id_dn2(iconn) = &
                    StructGetTVDGhostConnection(ghosted_id, &
                                                structured_grid, &
                                                iface,option)
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
! StructGridComputeVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine StructGridComputeVolumes(radius,structured_grid,option,nL2G,volume)

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
  PetscReal :: r1, r2
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(volume,volume_p, ierr)

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
        r1 = radius(ghosted_id) - 0.5d0 * structured_grid%dx(ghosted_id)
        r2 = r1 + structured_grid%dx(ghosted_id)
        volume_p(local_id) = pi * (r2*r2 - r1*r1) * structured_grid%dz(ghosted_id)
      enddo
    case(SPHERICAL_GRID)
      do local_id=1, structured_grid%nlmax
        ghosted_id = nL2G(local_id)
        r1 = radius(ghosted_id) - 0.5d0 * structured_grid%dx(ghosted_id)
        r2 = r1 + structured_grid%dx(ghosted_id)
        volume_p(local_id) = 4.d0/3.d0 * pi * (r2*r2*r2 - r1*r1*r1)
      enddo
  end select
  
  call VecRestoreArrayF90(volume,volume_p, ierr)
  
  if (OptionPrintToScreen(option) .and. &
      option%mycommsize > 1 .and. option%mycommsize <= 16) then
    write(*,'(" rank= ",i3,", nlmax= ",i6,", nlx,y,z= ",3i4, &
      & ", nxs,e = ",2i4,", nys,e = ",2i4,", nzs,e = ",2i4)') &
      option%myrank,structured_grid%nlmax,structured_grid%nlx, &
        structured_grid%nly,structured_grid%nlz,structured_grid%lxs, &
        structured_grid%lxe,structured_grid%lys,structured_grid%lye, &
        structured_grid%lzs,structured_grid%lze

    write(*,'(" rank= ",i3,", ngmax= ",i6,", ngx,y,z= ",3i4, &
      & ", ngxs,e= ",2i4,", ngys,e= ",2i4,", ngzs,e= ",2i4)') &
      option%myrank,structured_grid%ngmax,structured_grid%ngx, &
        structured_grid%ngy,structured_grid%ngz,structured_grid%gxs, &
        structured_grid%gxe,structured_grid%gys,structured_grid%gye, &
        structured_grid%gzs,structured_grid%gze
  endif

end subroutine StructGridComputeVolumes

! ************************************************************************** !
!
! StructGridMapIndices: maps global, local and natural indices of cells 
!                          to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine StructGridMapIndices(structured_grid,stencil_type,lsm_flux_method, &
                                    nG2L,nL2G, &
                                    nG2A,ghosted_level,option)

  use Option_module

  implicit none
  
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"  

  type(structured_grid_type) :: structured_grid
  PetscInt :: stencil_type
  PetscBool :: lsm_flux_method
  PetscInt, pointer :: nG2L(:), nL2G(:), nG2A(:), ghosted_level(:)
  type(option_type) :: option

  PetscInt :: i, j, k, local_id, ghosted_id, natural_id, count1
  PetscErrorCode :: ierr
  
  allocate(nL2G(structured_grid%nlmax))
  allocate(nG2L(structured_grid%ngmax))
  allocate(nG2A(structured_grid%ngmax))

  structured_grid%istart = structured_grid%lxs-structured_grid%gxs
  structured_grid%jstart = structured_grid%lys-structured_grid%gys
  structured_grid%kstart = structured_grid%lzs-structured_grid%gzs
  structured_grid%iend = structured_grid%istart+structured_grid%nlx-1
  structured_grid%jend = structured_grid%jstart+structured_grid%nly-1
  structured_grid%kend = structured_grid%kstart+structured_grid%nlz-1

  ! Local <-> Ghosted Transformation
  nG2L = 0  ! Must initialize this to zero!
  nL2G = 0

  nG2A = 0

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

  ! if STAR stencil, need to set corner ghosted cells to -1
  if (stencil_type == DMDA_STENCIL_STAR) then
    !geh - set corner ghosted nodes to -1
    do k=1,structured_grid%ngz
      do j=1,structured_grid%ngy
        do i=1,structured_grid%ngx
          count1 = 0
          if (i == 1 .and. &
              abs(structured_grid%lxs-structured_grid%gxs) > 0) &
            count1 = count1 + 1
          if (i == structured_grid%ngx .and. &
              abs(structured_grid%gxe-structured_grid%lxe) > 0) &
            count1 = count1 + 1
          if (j == 1 .and. &
              abs(structured_grid%lys-structured_grid%gys) > 0) &
            count1 = count1 + 1
          if (j == structured_grid%ngy .and. &
              abs(structured_grid%gye-structured_grid%lye) > 0) &
            count1 = count1 + 1
          if (k == 1 .and. &
              abs(structured_grid%lzs-structured_grid%gzs) > 0) &
            count1 = count1 + 1
          if (k == structured_grid%ngz .and. &
              abs(structured_grid%gze-structured_grid%lze) > 0) &
            count1 = count1 + 1
          if (count1 > 1) then
            ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
            nG2L(ghosted_id) = -1
          endif
        enddo
      enddo
    enddo
  endif

  ! local ghosted -> natural (1-based)
  ghosted_id = 0
  do k=1,structured_grid%ngz
    do j=1,structured_grid%ngy
      do i=1,structured_grid%ngx
        ghosted_id = ghosted_id + 1
        natural_id = i + structured_grid%gxs + & ! 1-based
                      (j-1+structured_grid%gys)*structured_grid%nx+ &
                      (k-1+structured_grid%gzs)*structured_grid%nxy
        nG2A(ghosted_id) = natural_id
      enddo
    enddo
  enddo

  if (lsm_flux_method) then
  
    allocate(ghosted_level(structured_grid%ngmax))
    
    ! Save information about ghost cell belong in which level of SNES stencil
    ghosted_level = 0
    do k=structured_grid%gzs,structured_grid%gze-1
      do j=structured_grid%gys,structured_grid%gye-1
        do i=structured_grid%gxs,structured_grid%gxe-1
          if(i<structured_grid%lxs) then
            ghosted_id = StructGridGetGhostedIDFromIJK( structured_grid, &
                                                        i-structured_grid%gxs+1, &
                                                        j-structured_grid%gys+1, &
                                                        k-structured_grid%gzs+1)
            ghosted_level(ghosted_id) = structured_grid%lxs-i
          endif
          if(i>structured_grid%lxe-1) then
            ghosted_id = StructGridGetGhostedIDFromIJK( structured_grid, &
                                                        i-structured_grid%gxs+1, &
                                                        j-structured_grid%gys+1, &
                                                        k-structured_grid%gzs+1)
            ghosted_level(ghosted_id) = -(structured_grid%lxe-1-i)
          endif

          if(j<structured_grid%lys) then
            ghosted_id = StructGridGetGhostedIDFromIJK( structured_grid, &
                                                        i-structured_grid%gxs+1, &
                                                        j-structured_grid%gys+1, &
                                                        k-structured_grid%gzs+1)
            ghosted_level(ghosted_id) = structured_grid%lys-j
          endif

          if(j>structured_grid%lye-1) then
            ghosted_id = StructGridGetGhostedIDFromIJK( structured_grid, &
                                                        i-structured_grid%gxs+1, &
                                                        j-structured_grid%gys+1, &
                                                        k-structured_grid%gzs+1)
            ghosted_level(ghosted_id) = -(structured_grid%lye-1-j)
          endif

          if(k<structured_grid%lzs) then
            ghosted_id = StructGridGetGhostedIDFromIJK( structured_grid, &
                                                        i-structured_grid%gxs+1, &
                                                        j-structured_grid%gys+1, &
                                                        k-structured_grid%gzs+1)
            ghosted_level(ghosted_id) = structured_grid%lzs-k
          endif

          if(k>structured_grid%lze-1) then
            ghosted_id = StructGridGetGhostedIDFromIJK( structured_grid, &
                                                        i-structured_grid%gxs+1, &
                                                        j-structured_grid%gys+1, &
                                                        k-structured_grid%gzs+1)
            ghosted_level(ghosted_id) = -(structured_grid%lze-1-k)
          endif
          
        enddo ! i-loop
      enddo ! j-loop
    enddo ! k-loop
  endif ! if LSM_FLUX

end subroutine StructGridMapIndices

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
#include "finclude/petscdmda.h"
  
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
    case(DMDA_STENCIL_STAR)
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
    case(DMDA_STENCIL_BOX)
      option%io_buffer = 'DMDA_STENCIL_BOX not yet supported in ' // &
        'StructGridGetNeighbors.'
      call printErrMsg(option)
  end select

end subroutine StructGridGetGhostedNeighbors


! ************************************************************************** !
!
! StructGridGetGhostedNeighborsCorners: Returns an array of neighboring cells
! including the corner nodes
! Note that the previous subroutine does not return the corner nodes
! author: Satish Karra, LANL
! date: 02/19/12
!
! ************************************************************************** !
subroutine StructGridGetGhostedNeighborsCorners(structured_grid,ghosted_id, &
                                         stencil_type, &
                                         stencil_width_i,stencil_width_j, &
                                         stencil_width_k, icount, &
                                         ghosted_neighbors, &
                                         option)

  use Option_module

  implicit none
#include "finclude/petscdmda.h"
  
  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscInt :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: ghosted_neighbors(*)

  PetscInt :: i, j, k
  PetscInt :: icount
  PetscInt :: ii, jj, kk

  call StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)

  icount = 0
  
  ! gb:08/08/13 Dependence on stencil_type is not necessary.
  !select case(stencil_type)
  !  case(DMDA_STENCIL_STAR)
      do kk = max(k-stencil_width_k,1), &
                min(k+stencil_width_k,structured_grid%ngz)
        do jj = max(j-stencil_width_j,1), &
                  min(j+stencil_width_j,structured_grid%ngy)
          do ii = max(i-stencil_width_i,1), &
                    min(i+stencil_width_i,structured_grid%ngx)
            if (ii == i .and. jj == j .and. kk == k) then
            ! do nothing
            else
              icount = icount + 1
              ghosted_neighbors(icount) = &
              StructGridGetGhostedIDFromIJK(structured_grid,ii,jj,kk)
            endif
          enddo
        enddo          
      enddo
  !  case(DMDA_STENCIL_BOX)
  !    option%io_buffer = 'DMDA_STENCIL_BOX not yet supported in ' // &
  !      'StructGridGetNeighbors.'
  !    call printErrMsg(option)
  !end select

end subroutine StructGridGetGhostedNeighborsCorners

! ************************************************************************** !
!
! StructGridDestroy: Deallocates a structured grid
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine StructGridDestroy(structured_grid)

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
  
  if(associated(structured_grid%cell_neighbors)) deallocate(structured_grid%cell_neighbors)
  nullify(structured_grid%cell_neighbors)
  
  deallocate(structured_grid)
  nullify(structured_grid)

end subroutine StructGridDestroy

! ************************************************************************** !
!
! StructGridCreateTVDGhosts: Calculates the TVD ghost vector and the 
!                            associated scatter context
! author: Glenn Hammond
! date: 01/28/11
!
! ************************************************************************** !
subroutine StructGridCreateTVDGhosts(structured_grid,ndof,global_vec, &
                                     dm_1dof, &
                                     ghost_vec,scatter_ctx,option)

  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

  type(structured_grid_type) :: structured_grid
  PetscInt :: ndof
  DM :: dm_1dof
  Vec :: global_vec
  Vec :: ghost_vec
  VecScatter :: scatter_ctx
  type(option_type) :: option

  PetscInt :: vector_size
  IS :: is_petsc
  IS :: is_ghost
  PetscInt :: icount, index, offset
  PetscInt :: increment
  PetscInt :: i, j, k
  PetscInt, allocatable :: global_indices_of_local_ghosted(:)
  PetscInt, allocatable :: global_indices_from(:)
  PetscInt, allocatable :: tvd_ghost_indices_to(:)
  ISLocalToGlobalMapping :: mapping_ltog
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  ! structured grid has 6 sides to it
  vector_size = 0
  ! east-west
  if (structured_grid%nx > 1) then
    increment = structured_grid%nly*structured_grid%nlz
    vector_size = vector_size + 2*increment
  endif
  ! north-south
  if (structured_grid%ny > 1) then
    increment = structured_grid%nlx*structured_grid%nlz
    vector_size = vector_size + 2*increment
  endif
  ! top-bottom
  if (structured_grid%nz > 1) then
    increment = structured_grid%nlx*structured_grid%nly
    vector_size = vector_size + 2*increment
  endif
  
  if (vector_size == 0) then
    option%io_buffer = 'TVD does not handle a single grid cell.'
    call printErrMsg(option)
  endif
  
  call VecCreateSeq(PETSC_COMM_SELF,vector_size*ndof,ghost_vec,ierr)
  call VecSetBlockSize(ghost_vec,ndof,ierr)
  
  ! Create an IS composed of the petsc indexing of the ghost cells
  allocate(global_indices_from(vector_size))
  global_indices_from = -999 ! to catch bugs
  allocate(tvd_ghost_indices_to(vector_size))
  do i = 1, vector_size
    tvd_ghost_indices_to(i) = i-1
  enddo
  
  ! GEH: I'm going to play a trick here.  If I know the global index
  ! of my ghost cells, I can calculate the global index of the next
  ! layer of ghost cells in each direction since the block are 
  ! consistent through each dimension of the grid
  allocate(global_indices_of_local_ghosted(structured_grid%ngmax))
  do i = 1, structured_grid%ngmax
    global_indices_of_local_ghosted(i) = i-1
  enddo
  call DMGetLocalToGlobalMapping(dm_1dof,mapping_ltog,ierr)
  ! in and out integer arrays can be the same
  call ISLocalToGlobalMappingApply(mapping_ltog,structured_grid%ngmax, &
                                   global_indices_of_local_ghosted, &
                                   global_indices_of_local_ghosted,ierr)
  ! leave global_indices_of_local_ghosted() in zero-based for the below
  
  ! Need to make a list of all indices that will receive updates through
  ! scatter/gather operation. Ghost cells representing physical boundaries
  ! do not need such an update.
  icount = 0

  if (structured_grid%nx > 1) then
    ! west
    offset = 0
    if (structured_grid%lxs /= structured_grid%gxs) offset = -1
    i = structured_grid%istart
    do k = structured_grid%kstart, structured_grid%kend
      do j = structured_grid%jstart, structured_grid%jend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo

    ! east
    offset = 0
    if (structured_grid%lxe /= structured_grid%gxe) offset = 1
    i = structured_grid%iend
    do k = structured_grid%kstart, structured_grid%kend
      do j = structured_grid%jstart, structured_grid%jend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  endif

  if (structured_grid%ny > 1) then
    ! south
    offset = 0
    if (structured_grid%lys /= structured_grid%gys) offset = -structured_grid%ngx
    j = structured_grid%jstart
    do k = structured_grid%kstart, structured_grid%kend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  
    ! north
    offset = 0
    if (structured_grid%lye /= structured_grid%gye) offset = structured_grid%ngx
    j = structured_grid%jend
    do k = structured_grid%kstart, structured_grid%kend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  endif
  
  if (structured_grid%nz > 1) then
    ! bottom
    offset = 0
    if (structured_grid%lzs /= structured_grid%gzs) offset = -structured_grid%ngxy
    k = structured_grid%kstart
    do j = structured_grid%jstart, structured_grid%jend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  
    ! top
    offset = 0
    if (structured_grid%lze /= structured_grid%gze) offset = structured_grid%ngxy
    k = structured_grid%kend
    do j = structured_grid%jstart, structured_grid%jend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  endif
  
  deallocate(global_indices_of_local_ghosted)

  if (vector_size /= icount) then
    option%io_buffer = 'Mis-count in TVD ghosting.'
    call printErrMsgByRank(option)
  endif

  ! since global_indices_from was base-zero, global_indices_from is base-zero.
  call ISCreateBlock(option%mycomm,ndof,vector_size, &
                      global_indices_from,PETSC_COPY_VALUES,is_petsc,ierr)
  deallocate(global_indices_from)

#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_petsc_tvd.out', &
                            viewer,ierr)
  call ISView(is_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! already zero-based
  call ISCreateBlock(option%mycomm,ndof,vector_size, &
                      tvd_ghost_indices_to,PETSC_COPY_VALUES,is_ghost,ierr)
  deallocate(tvd_ghost_indices_to)

#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_ghost_tvd.out', &
                            viewer,ierr)
  call ISView(is_ghost,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecScatterCreate(global_vec,is_petsc,ghost_vec,is_ghost, &
                        scatter_ctx,ierr)

#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'tvd_ghost_scatter.out',viewer,ierr)
  call VecScatterView(scatter_ctx,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

end subroutine StructGridCreateTVDGhosts  

! ************************************************************************** !
!
! StructGetTVDGhostConnection: Returns id of tvd ghost cell for connection
! author: Glenn Hammond
! date: 02/10/11
!
! ************************************************************************** !
function StructGetTVDGhostConnection(ghosted_id,structured_grid,iface,option)

  use Option_module

  implicit none
  
  PetscInt :: ghosted_id
  type(structured_grid_type) :: structured_grid
  PetscInt :: iface
  type(option_type) :: option
  
  PetscInt :: StructGetTVDGhostConnection
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: index
  PetscInt :: offset
  PetscInt :: i, j, k
  PetscBool :: error
  
  error = PETSC_FALSE

  select case(iface)
    case(WEST_FACE,EAST_FACE)
      if (structured_grid%ngx > 1) then
        if (iface == WEST_FACE) then
          StructGetTVDGhostConnection = ghosted_id + 1
        else
          StructGetTVDGhostConnection = ghosted_id - 1
        endif
        return
      elseif (structured_grid%nx == 1) then
        option%io_buffer = 'Boundary condition cannot be assigned in X ' // &
          'dimension with nx = 1 with TVD.'
        error = PETSC_TRUE
      endif
    case(SOUTH_FACE,NORTH_FACE)
      if (structured_grid%ngy > 1) then
        if (iface == SOUTH_FACE) then
          StructGetTVDGhostConnection = ghosted_id + structured_grid%ngx
        else
          StructGetTVDGhostConnection = ghosted_id - structured_grid%ngx
        endif
        return
      elseif (structured_grid%ny == 1) then
        option%io_buffer = 'Boundary condition cannot be assigned in Y ' // &
          'dimension with ny = 1 with TVD.'
        error = PETSC_TRUE
      endif
    case(BOTTOM_FACE,TOP_FACE)
      if (structured_grid%ngz > 1) then
        if (iface == BOTTOM_FACE) then
          StructGetTVDGhostConnection = ghosted_id + structured_grid%ngxy
        else
          StructGetTVDGhostConnection = ghosted_id - structured_grid%ngxy
        endif
        return
      elseif (structured_grid%nz == 1) then
        option%io_buffer = 'Boundary condition cannot be assigned in Z ' // &
          'dimension with nz = 1 with TVD.'
        error = PETSC_TRUE
      endif
  end select

  if (error) call printErrMsg(option)
  
  call StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)
  offset = 0
  select case(iface)
    case(WEST_FACE)
      ! ensure that connection is on boundary face
      if (i /= 1) then
        error = PETSC_TRUE
        string = 'WEST'
      endif                       ! this must be in local dimension nly
      index = j + k*structured_grid%nly
    case(EAST_FACE)
      if (i /= structured_grid%ngx) then
        error = PETSC_TRUE
        string = 'EAST'
      endif
      index = j + k*structured_grid%nly + structured_grid%nlyz
    case(SOUTH_FACE)
      if (j /= 1) then
        error = PETSC_TRUE
        string = 'SOUTH'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      index = i + k*structured_grid%nlx + offset
    case(NORTH_FACE)
      if (j /= structured_grid%ngy) then
        error = PETSC_TRUE
        string = 'NORTH'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      index = i + k*structured_grid%nlx + structured_grid%nlxz + offset
    case(BOTTOM_FACE)
      if (k /= 1) then
        error = PETSC_TRUE
        string = 'BOTTOM'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      if (structured_grid%ny > 1) then
        offset = offset + 2*structured_grid%nlxz
      endif
      index = i + j*structured_grid%nlx + offset
    case(TOP_FACE)
      if (k /= structured_grid%ngz) then
        error = PETSC_TRUE
        string = 'TOP'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      if (structured_grid%ny > 1) then
        offset = offset + 2*structured_grid%nlxz
      endif
      index = i + j*structured_grid%nlx + structured_grid%nlxy + offset
  end select
  
  if (error) then
    write(option%io_buffer, '(''StructGetTVDGhostConnection not on '', a, &
    & ''face for cell:'',3i6)') trim(string), i,j,k
    call printErrMsgByRank(option)
  endif
  
  StructGetTVDGhostConnection = -index

end function StructGetTVDGhostConnection

! ************************************************************************** !
!> This routine saves indices (in ghosted order) of neighbors for all ghosted
!! cells.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 08/24/12
! ************************************************************************** !
subroutine StructGridComputeNeighbors(structured_grid,nG2L,is_bnd_vec,option)

  use Option_module

  implicit none
  
#include "finclude/petscdmda.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(structured_grid_type) :: structured_grid
  type(option_type) :: option
  PetscInt, pointer :: nG2L(:)
  Vec :: is_bnd_vec
  
  PetscInt :: ghosted_id, ncount, ghosted_neighbors(26)
  PetscInt :: local_id
  PetscScalar, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  allocate(structured_grid%cell_neighbors(0:26,structured_grid%ngmax))
  structured_grid%cell_neighbors = 0

  call VecGetArrayF90(is_bnd_vec,vec_ptr,ierr)
  do ghosted_id = 1,structured_grid%ngmax
    local_id=nG2L(ghosted_id)
    call StructGridGetGhostedNeighborsCorners(structured_grid,ghosted_id, &
                                         DMDA_STENCIL_BOX, &
                                         ONE_INTEGER_MPI, &
                                         ONE_INTEGER_MPI, &
                                         ONE_INTEGER_MPI, &
                                         ncount, &
                                         ghosted_neighbors, &
                                         option)
    structured_grid%cell_neighbors(0,ghosted_id) = ncount
    structured_grid%cell_neighbors(1:ncount,ghosted_id) = ghosted_neighbors(1:ncount)
    if(local_id>0) then
      if(ncount == 26) then
        vec_ptr(ghosted_id) = 0.d0
      else
        vec_ptr(ghosted_id) = 1.d0
      endif
    endif
  enddo
  call VecRestoreArrayF90(is_bnd_vec,vec_ptr,ierr)
  
end subroutine StructGridComputeNeighbors

end module Structured_Grid_module
