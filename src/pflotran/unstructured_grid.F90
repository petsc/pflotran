module Unstructured_Grid_module

 use pflow_gridtype_module

 implicit none

#define HASH
#define INVERT

 private 

#include "include/finclude/petsc.h"
!#include "include/petscf90.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petsclog.h"

#include "definitions.h"

  integer, pointer :: hash(:,:,:)
  integer :: num_hash = 100
  
  integer :: hdf5_err
  PetscErrorCode :: ierr
   
  public ReadUnstructuredGrid, ReadMaterials, ReadMaterials2, ReadStructuredGridHDF5
  
contains
  
! ************************************************************************** !
!
! ReadUnstructuredGrid: Reads in an unstructured grid (from a format similar
!                       to TOUGH and sets up intercelluar and boundary
!                       connectivity
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine ReadUnstructuredGrid(grid)

  use fileio_module
  use Condition_module

  implicit none

  type(pflowGrid) :: grid
  character(len=MAXSTRINGLENGTH) :: string 
 
  PetscScalar, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:), &
                          tor_p(:), volume_p(:), icap_p(:)
  PetscTruth :: option_found
  integer fid, iostatus, ierr
  integer :: connection_id, material_id, natural_id, local_id, local_ghosted_id
  integer :: natural_id_upwind, natural_id_downwind
  integer :: local_id_upwind, local_id_downwind
  integer :: local_ghosted_id_upwind, local_ghosted_id_downwind
! integer :: number_of_times
  integer :: count1, count2, irecv, i
  real*8 :: xcoord, ycoord, zcoord, xperm, yperm, zperm
  real*8 :: area, distance, distance_upwind, distance_downwind, cosB
  real*8 :: volume, porosity, tortuosity
  real*8 :: dx, dy, dz, delz
  character(len=MAXCARDLENGTH) :: card
! character(len=MAXWORDLENGTH) :: bctype
  character(len=MAXWORDLENGTH) :: time_unit
  character(len=MAXSTRINGLENGTH) :: filename
  integer :: direction, icond, icondition
  real*8 :: time, value, time_multiplier
  PetscLogDouble :: time0, time1
  type(condition_type), pointer :: new_condition

  call VecGetArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr);  CHKERRQ(ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr);  CHKERRQ(ierr)

  call PetscGetTime(time0, ierr)

  filename = "hanford.in"
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-unstructured_filename", &
    filename, option_found, ierr)

  fid = 86
  open(fid, file=filename, action="read", status="old", iostat=iostatus)
  if (iostatus /= 0) then
    print *, 'Error opening file: ', trim(filename), ' on Processor ', &
             grid%myrank
    call PetscFinalize()
    stop
  endif

  ! create hash table for fast lookup
#ifdef HASH
  call CreateNaturalToLocalGhostedHash(grid)
!  call PrintHashTable
#endif

! GRID information
  card = "GRID"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif
  
  count1 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,natural_id,ierr) 
    call fiErrorMsg('natural_id',card,ierr)

    call fiReadInt(string,material_id,ierr) 
    call fiErrorMsg('material_id',card,ierr)

    call fiReadDouble(string,xcoord,ierr)
    call fiErrorMsg('xcoord',card,ierr)
   
    call fiReadDouble(string,ycoord,ierr)
    call fiErrorMsg('ycoord',card,ierr)

    call fiReadDouble(string,zcoord,ierr)
    call fiErrorMsg('zcoord',card,ierr)
   
    call fiReadDouble(string,volume,ierr)
    call fiErrorMsg('volume',card,ierr)
   
    call fiReadDouble(string,xperm,ierr)
    call fiErrorMsg('xperm',card,ierr)

    call fiReadDouble(string,yperm,ierr)
    call fiErrorMsg('yperm',card,ierr)
   
    call fiReadDouble(string,zperm,ierr)
    call fiErrorMsg('zperm',card,ierr)

    call fiReadDouble(string,porosity,ierr)
    call fiErrorMsg('porosity',card,ierr)

    call fiReadDouble(string,tortuosity,ierr)
    call fiErrorMsg('tortuosity',card,ierr)

#ifdef HASH
    local_ghosted_id = GetLocalGhostedIdFromHash(natural_id)
#else
    local_ghosted_id = GetLocalGhostedIdFromNaturalId(natural_id,grid)
#endif
    if (local_ghosted_id > 0) then

      grid%x(local_ghosted_id)=xcoord
      grid%y(local_ghosted_id)=ycoord
!      grid%z(local_ghosted_id)=zcoord

      local_id = grid%nG2L(local_ghosted_id)
      if (local_id > 0 .and. abs(material_id) > 0) then
        perm_xx_p(local_id) = xperm
        perm_yy_p(local_id) = yperm
        perm_zz_p(local_id) = zperm
        volume_p(local_id) = volume
        por_p(local_id) = porosity
        tor_p(local_id) = tortuosity
      endif
      grid%imat(local_ghosted_id) = material_id
    endif

#if DEBUG
    print *, 'Read geom 1:',natural_id,xcoord,ycoord,zcoord,volume, &
                            xperm,yperm,zperm,porosity,tortuosity
#endif
    count1 = count1 + 1
  enddo

  call VecRestoreArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)     
  call VecRestoreArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%tor,tor_p,ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%porosity,por_p,ierr); CHKERRQ(ierr)

  if (grid%myrank == 0) print *, count1, ' cells read'
  
! CONNection information

! Count the number of local connections
  grid%nconn = GetNumberOfLocalConnections(fid,grid) ! function at bottom of file
! Allocate Arrays
  allocate(grid%nd1(grid%nconn))
  allocate(grid%nd2(grid%nconn))
  allocate(grid%dist1(grid%nconn))
  allocate(grid%dist2(grid%nconn))
  allocate(grid%area(grid%nconn))
  allocate(grid%delz(grid%nconn))
  allocate(grid%grav_ang(grid%nconn))

  allocate(grid%iperm1(grid%nconn))
  allocate(grid%iperm2(grid%nconn))

  allocate(grid%vl_loc(grid%nconn))
  allocate(grid%vvl_loc(grid%nconn))
  allocate(grid%vg_loc(grid%nconn))
  allocate(grid%vvg_loc(grid%nconn))

  grid%nconn = 0

  card = "CONN"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif
  
  count1 = 0
  count2 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

! stomp2pflotran header
! :id idup iddn distup distdn area cosB
    call fiReadInt(string,connection_id,ierr) 
    call fiErrorMsg('connection_id',card,ierr)

    call fiReadInt(string,natural_id_upwind,ierr) 
    call fiErrorMsg('natural_id_upwind',card,ierr)

    call fiReadInt(string,natural_id_downwind,ierr) 
    call fiErrorMsg('natural_id_downwind',card,ierr)

    local_id_upwind = 0
    local_id_downwind = 0
#ifdef HASH
    local_ghosted_id_upwind = GetLocalGhostedIdFromHash(natural_id_upwind)
    local_ghosted_id_downwind = GetLocalGhostedIdFromHash(natural_id_downwind)
#else
    local_ghosted_id_upwind = &
               GetLocalGhostedIdFromNaturalId(natural_id_upwind,grid)
    local_ghosted_id_downwind = &
               GetLocalGhostedIdFromNaturalId(natural_id_downwind,grid)
#endif
    if (local_ghosted_id_upwind > 0 .and. &
        local_ghosted_id_downwind > 0) then ! both cells are local ghosted
      local_id_upwind = grid%nG2L(local_ghosted_id_upwind)
      local_id_downwind = grid%nG2L(local_ghosted_id_downwind)
      ! at least 1 cell must be local (non-ghosted) since we don't want to
      ! create a connection between two ghosted cells
      if (local_id_upwind > 0 .or. local_id_downwind > 0) then
        count2 = count2 + 1

        ! Continue reading string if local
        call fiReadDouble(string,distance_upwind,ierr) 
        call fiErrorMsg('distance_upwind',card,ierr)

        call fiReadDouble(string,distance_downwind,ierr) 
        call fiErrorMsg('distance_downwind',card,ierr)

        call fiReadDouble(string,area,ierr) 
        call fiErrorMsg('area',card,ierr)

        call fiReadDouble(string,cosB,ierr) 
        call fiErrorMsg('cosB',card,ierr)

        grid%nconn = grid%nconn + 1
        grid%nd1(grid%nconn) = local_ghosted_id_upwind
        grid%nd2(grid%nconn) = local_ghosted_id_downwind
        grid% dist1(grid%nconn) = distance_upwind
        grid% dist2(grid%nconn) = distance_downwind
        if (area>0.D0) grid%area(grid%nconn)= area
        ! negate to account for left-hand rule in pflotran
#ifdef INVERT
        grid%delz(grid%nconn) = -1.d0*abs(grid%z(local_ghosted_id_downwind)-  &
                                    grid%z(local_ghosted_id_upwind))
#else
        grid%delz(grid%nconn) = +1.d0*abs(grid%z(local_ghosted_id_downwind)-  &
                                    grid%z(local_ghosted_id_upwind))
#endif
        grid%grav_ang(grid%nconn) = cosB

        ! setup direction of permeability
        dx = abs(grid%x(local_ghosted_id_upwind)-grid%x(local_ghosted_id_downwind))
        dy = abs(grid%y(local_ghosted_id_upwind)-grid%y(local_ghosted_id_downwind))
        dz = abs(grid%z(local_ghosted_id_upwind)-grid%z(local_ghosted_id_downwind))
        if (dx > dy .and. dx > dz) then
          grid%iperm1(grid%nconn) = 1
          grid%iperm2(grid%nconn) = 1
        else if (dy > dx .and. dy > dz) then
          grid%iperm1(grid%nconn) = 2
          grid%iperm2(grid%nconn) = 2
        else if (dz > dx .and. dz > dy) then
          grid%iperm1(grid%nconn) = 3
          grid%iperm2(grid%nconn) = 3
        endif
      endif
    endif
    count1 = count1 + 1
!print *, grid%myrank, local_ghosted_id_upwind, local_id_upwind, local_ghosted_id_downwind, local_id_downwind
  enddo
  call mpi_reduce(count2,irecv,1,MPI_INTEGER,MPI_SUM,0, &
                     PETSC_COMM_WORLD,ierr)
  if (grid%myrank == 0) print *, irecv, ' of ', count1, ' connections found'


! BCONection information

!geh Right now we ignore allocation since the arrays are already allocated
!geh larger than necessary.  This will need to be updated in the future.
! Count the number of boundary connections
  grid%nconnbc = GetNumberOfBoundaryConnections(fid,grid) ! function at bottom of file
! Allocate Arrays
  allocate(grid%mblkbc(grid%nconnbc)) ! id of local cell
  allocate(grid%ibconn(grid%nconnbc)) ! id of condition (boundary)
  allocate(grid%distbc(grid%nconnbc))
  allocate(grid%areabc(grid%nconnbc))
  allocate(grid%ipermbc(grid%nconnbc))
  allocate(grid%delzbc(grid%nconnbc))

  allocate(grid%vlbc(grid%nconnbc))
  allocate(grid%vvlbc(grid%nconnbc))
  allocate(grid%vgbc(grid%nconnbc))
  allocate(grid%vvgbc(grid%nconnbc))

  allocate(grid%velocitybc(grid%nphase, grid%nconnbc))
  allocate(grid%iphasebc(grid%nconnbc))
  allocate(grid%xxbc(grid%ndof,grid%nconnbc))
  allocate(grid%xphi_co2_bc(grid%nconnbc))
  allocate(grid%xxphi_co2_bc(grid%nconnbc))

  grid%nconnbc = 0

  card = "BCON"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif
  
  count1 = 0
  count2 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

! stomp2pflotran header
!:id cellid dist area delz icond

    call fiReadInt(string,connection_id,ierr) 
    call fiErrorMsg('connection_id',card,ierr)

    call fiReadInt(string,natural_id,ierr) 
    call fiErrorMsg('natural_id',card,ierr)

#ifdef HASH
    local_ghosted_id = GetLocalGhostedIdFromHash(natural_id)
#else
    local_ghosted_id = GetLocalGhostedIdFromNaturalId(natural_id,grid)
#endif
    if (local_ghosted_id > 0) then 
      local_id = grid%nG2L(local_ghosted_id)
      if (local_id > 0) then
        grid%nconnbc = grid%nconnbc + 1
        count2 = count2 + 1

        call fiReadInt(string,direction,ierr) 
        call fiErrorMsg('direction',card,ierr)

        call fiReadDouble(string,distance,ierr) 
        call fiErrorMsg('distance',card,ierr)

        call fiReadDouble(string,area,ierr) 
        call fiErrorMsg('area',card,ierr)

        call fiReadDouble(string,delz,ierr) 
        call fiErrorMsg('delz',card,ierr)

        call fiReadInt(string,icond,ierr)
        call fiErrorMsg('icond',card,ierr)
 
        grid%mblkbc(grid%nconnbc) = local_id
        grid%distbc(grid%nconnbc) = distance
        grid%areabc(grid%nconnbc) = area
        ! negate to account for left-hand rule in pflotran
#ifdef INVERT
        grid%delzbc(grid%nconnbc) = +1.d0*delz
#else
        grid%delzbc(grid%nconnbc) = +1.d0*delz
#endif
        grid%ibconn(grid%nconnbc) = icond

        ! setup direction of permeability
        grid%ipermbc(grid%nconnbc) = abs(direction) ! 1,2,3 -> x,y,z

        ! hardwired for now
        grid%iphasebc(grid%nconnbc) = 1

      endif
    endif
    count1 = count1 + 1
  enddo

  call mpi_reduce(count2,irecv,1,MPI_INTEGER,MPI_SUM,0, &
                     PETSC_COMM_WORLD,ierr)
  if (grid%myrank == 0) print *, irecv, ' of ', count1, ' connections found'

! CONDition information

  card = "COND"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif

  count1 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

! stomp2pflotran header
!:id "type" num_values (datum)

    allocate(new_condition)
    nullify(new_condition%next)
    new_condition%itype = 0
    new_condition%xdatum = 0.d0
    new_condition%ydatum = 0.d0
    new_condition%zdatum = 0.d0
    new_condition%xgrad = 0.d0
    new_condition%ygrad = 0.d0
    new_condition%zgrad = 0.d0

    call fiReadInt(string,new_condition%id,ierr) 
    call fiErrorMsg('id',card,ierr)

    call fiReadQuotedNChars(string,new_condition%ctype,MAXWORDLENGTH, &
                            .true.,ierr) 
    call fiErrorMsg('condition_type',card,ierr)
    call fiWordToLower(new_condition%ctype)

    call fiReadInt(string,new_condition%max_time_index,ierr) 
    call fiErrorMsg('max_time_index',card,ierr)

    call fiReadWord(string,time_unit,.true.,ierr) 
    call fiErrorMsg('time unit',card,ierr)

    if (fiStringCompare("hydraulic gradient",new_condition%ctype,18) .or. &
        fiStringCompare("seepage",new_condition%ctype,7)) then

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,new_condition%xdatum,ierr) 
      call fiErrorMsg('x datum',card,ierr)

      call fiReadDouble(string,new_condition%ydatum,ierr) 
      call fiErrorMsg('y datum',card,ierr)

      call fiReadDouble(string,new_condition%zdatum,ierr) 
      call fiErrorMsg('z datum',card,ierr)

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,new_condition%xgrad,ierr) 
      call fiErrorMsg('x gradient',card,ierr)

      call fiReadDouble(string,new_condition%ygrad,ierr) 
      call fiErrorMsg('y gradient',card,ierr)

      call fiReadDouble(string,new_condition%zgrad,ierr) 
      call fiErrorMsg('z gradient',card,ierr)
    endif

    allocate(new_condition%times(new_condition%max_time_index))
    allocate(new_condition%values(new_condition%max_time_index))
    new_condition%times = 0.d0
    new_condition%values = 0.d0

    call fiWordToLower(time_unit)
    if (fiStringCompare("y",time_unit,1)) then
      time_multiplier = 365.d0*24.d0*3600.d0
    else if (fiStringCompare("d",time_unit,1)) then
      time_multiplier = 24.d0*3600.d0
    else if (fiStringCompare("h",time_unit,1)) then
      time_multiplier = 3600.d0
    else if (fiStringCompare("m",time_unit,1)) then
      time_multiplier = 60.d0
    else if (fiStringCompare("s",time_unit,1)) then
      time_multiplier = 1.d0
    else
      if (grid%myrank == 0) then
        print *, 'Time unit: ', trim(time_unit), ' not recognized'
      endif
      call PetscFinalize()
      stop
    endif

    if (grid%myrank == 0) print *, new_condition%ctype
    if (grid%myrank == 0) print *, new_condition%max_time_index

    count2 = 0
    do i=1,new_condition%max_time_index

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,time,ierr)
      call fiErrorMsg('time',card,ierr)

      call fiReadDouble(string,value,ierr)
      call fiErrorMsg('value',card,ierr)
      ! times are all in seconds   yr -> seconds
      new_condition%times(i) = time*time_multiplier
      new_condition%values(i) = value

      count2 = count2 + 1

    enddo

    if (grid%myrank == 0) print *, count2, ' times read'

    call AddCondition(new_condition)
    nullify(new_condition)

    count1 = count1 + 1

  enddo

  if (grid%myrank == 0) print *, count1, ' conditions read'

! Initial Condition information

  card = "ICON"
  call fiFindStringInFile(fid,card,ierr)

  ! report warning if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'WARNING: Card (',card, ') not found in file'

  icondition = -1
  if (ierr == 0) then
    do
      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadInt(string,icondition,ierr) 
      call fiErrorMsg('icondition',card,ierr)
    enddo
  endif

! set capillary function ids based on material id
  call VecGetArrayF90(grid%icap,icap_p,ierr)
  do local_id=1,grid%nlmax
    if (abs(grid%imat(grid%nL2G(local_id))) > 0) then
      icap_p(local_id) = abs(grid%imat(grid%nL2G(local_id)))
    else
      icap_p(local_id) = 1
    endif
  enddo
  call VecRestoreArrayF90(grid%icap,icap_p,ierr)

  call UpdateGlobalToLocal(grid)
  call InitializeBoundaryConditions(grid)
  if (icondition > 0) call ComputeInitialCondition(grid,icondition)

#ifdef HASH
  deallocate(hash)
#endif
  close(fid)
  call PetscGetTime(time1, ierr)
  time1 = time1 - time0
  if (grid%myrank == 0) print *, time1, ' seconds to read unstructured grid data'
 
end subroutine ReadUnstructuredGrid
 
! ************************************************************************** !
!
! ReadMaterials: Reads just the materials from the unstructured input file
! author: Glenn Hammond
! date: 06/20/07
!
! ************************************************************************** !
subroutine ReadMaterials(grid)

  use fileio_module
  use Condition_module

  implicit none

  type(pflowGrid) :: grid
  character(len=MAXSTRINGLENGTH) :: string 
 
  PetscScalar, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:), &
                          tor_p(:), volume_p(:), icap_p(:) !, temp_ptr(:)
  PetscTruth :: option_found
  integer fid, iostatus, ierr
  integer :: material_id, natural_id, local_id, local_ghosted_id
! integer :: connection_id
! integer :: natural_id_upwind, natural_id_downwind
! integer :: local_id_upwind, local_id_downwind
! integer :: local_ghosted_id_upwind, local_ghosted_id_downwind
! integer :: number_of_times
  integer :: count1, count2, i
! integer :: irecv
  real*8 :: xcoord, ycoord, zcoord, xperm, yperm, zperm
! real*8 :: area, distance, distance_upwind, distance_downwind, cosB
  real*8 :: volume, porosity, tortuosity
! real*8 :: dx, dy, dz, delz
  character(len=MAXCARDLENGTH) :: card
! character(len=MAXWORDLENGTH) :: bctype
  character(len=MAXWORDLENGTH) :: time_unit
  character(len=MAXSTRINGLENGTH) :: filename
! integer :: direction, icond
  integer :: icondition
  real*8 :: time, value, time_multiplier
  PetscLogDouble :: time0, time1
  type(condition_type), pointer :: new_condition

  call VecGetArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr);  CHKERRQ(ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr);  CHKERRQ(ierr)

  string = " "
  filename = "hanford.in"
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-unstructured_filename", &
    filename, option_found, ierr)

  call PetscGetTime(time0, ierr)

  fid = 86
  open(fid, file=filename, action="read", status="old", iostat=iostatus)
  if (iostatus /= 0) then
    print *, 'Error opening file: ', trim(filename), ' on Processor ', &
             grid%myrank
    call PetscFinalize()
    stop
  endif

  ! create hash table for fast lookup
#ifdef HASH
  call CreateNaturalToLocalGhostedHash(grid)
!  call PrintHashTable
#endif

! GRID information
  card = "GRID"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif

  count1 = 0
  do

    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,natural_id,ierr) 
    call fiErrorMsg('natural_id',card,ierr)

    call fiReadInt(string,material_id,ierr) 
    call fiErrorMsg('material_id',card,ierr)
    call fiReadDouble(string,xcoord,ierr)
    call fiErrorMsg('xcoord',card,ierr)
   
    call fiReadDouble(string,ycoord,ierr)
    call fiErrorMsg('ycoord',card,ierr)

    call fiReadDouble(string,zcoord,ierr)
    call fiErrorMsg('zcoord',card,ierr)
   
    call fiReadDouble(string,volume,ierr)
    call fiErrorMsg('volume',card,ierr)
   
    call fiReadDouble(string,xperm,ierr)
    call fiErrorMsg('xperm',card,ierr)

    call fiReadDouble(string,yperm,ierr)
    call fiErrorMsg('yperm',card,ierr)
   
    call fiReadDouble(string,zperm,ierr)
    call fiErrorMsg('zperm',card,ierr)

    call fiReadDouble(string,porosity,ierr)
    call fiErrorMsg('porosity',card,ierr)

    call fiReadDouble(string,tortuosity,ierr)
    call fiErrorMsg('tortuosity',card,ierr)

#ifdef HASH
    local_ghosted_id = GetLocalGhostedIdFromHash(natural_id)
#else
    local_ghosted_id = GetLocalGhostedIdFromNaturalId(natural_id,grid)
#endif
    if (local_ghosted_id > 0) then

      grid%x(local_ghosted_id)=xcoord
      grid%y(local_ghosted_id)=ycoord
      grid%z(local_ghosted_id)=zcoord

      local_id = grid%nG2L(local_ghosted_id)
      if (local_id > 0 .and. abs(material_id) > 0) then
!      if (local_id > 0) then
        perm_xx_p(local_id) = xperm
        perm_yy_p(local_id) = yperm
        perm_zz_p(local_id) = zperm
        volume_p(local_id) = volume
        por_p(local_id) = porosity
        tor_p(local_id) = tortuosity
      endif
      grid%imat(local_ghosted_id) = material_id
    endif

#if DEBUG
    print *, 'Read geom 1:',natural_id,xcoord,ycoord,zcoord,volume, &
                            xperm,yperm,zperm,porosity,tortuosity
#endif
    count1 = count1 + 1
  enddo

  call VecRestoreArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%tor,tor_p,ierr);  CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%porosity,por_p,ierr);  CHKERRQ(ierr)


  if (grid%myrank == 0) print *, count1, ' cells read'
  
! CONDition information

  card = "COND"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif

  count1 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

! stomp2pflotran header
!:id "type" num_values (datum)

    allocate(new_condition)
    nullify(new_condition%next)
    new_condition%itype = 0
    new_condition%xdatum = 0.d0
    new_condition%ydatum = 0.d0
    new_condition%zdatum = 0.d0
    new_condition%xgrad = 0.d0
    new_condition%ygrad = 0.d0
    new_condition%zgrad = 0.d0

    call fiReadInt(string,new_condition%id,ierr) 
    call fiErrorMsg('id',card,ierr)

    call fiReadQuotedNChars(string,new_condition%ctype,MAXWORDLENGTH, &
                            .true.,ierr) 
    call fiErrorMsg('condition_type',card,ierr)
    call fiWordToLower(new_condition%ctype)

    call fiReadInt(string,new_condition%max_time_index,ierr) 
    call fiErrorMsg('max_time_index',card,ierr)

    call fiReadWord(string,time_unit,.true.,ierr) 
    call fiErrorMsg('time unit',card,ierr)

    if (fiStringCompare("hydraulic gradient",new_condition%ctype,18) .or. &
        fiStringCompare("seepage",new_condition%ctype,7)) then

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,new_condition%xdatum,ierr) 
      call fiErrorMsg('x datum',card,ierr)

      call fiReadDouble(string,new_condition%ydatum,ierr) 
      call fiErrorMsg('y datum',card,ierr)

      call fiReadDouble(string,new_condition%zdatum,ierr) 
      call fiErrorMsg('z datum',card,ierr)

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,new_condition%xgrad,ierr) 
      call fiErrorMsg('x gradient',card,ierr)

      call fiReadDouble(string,new_condition%ygrad,ierr) 
      call fiErrorMsg('y gradient',card,ierr)

      call fiReadDouble(string,new_condition%zgrad,ierr) 
      call fiErrorMsg('z gradient',card,ierr)
    endif

    allocate(new_condition%times(new_condition%max_time_index))
    allocate(new_condition%values(new_condition%max_time_index))
    new_condition%times = 0.d0
    new_condition%values = 0.d0

    call fiWordToLower(time_unit)
    if (fiStringCompare("y",time_unit,1)) then
      time_multiplier = 365.d0*24.d0*3600.d0
    else if (fiStringCompare("d",time_unit,1)) then
      time_multiplier = 24.d0*3600.d0
    else if (fiStringCompare("h",time_unit,1)) then
      time_multiplier = 3600.d0
    else if (fiStringCompare("m",time_unit,1)) then
      time_multiplier = 60.d0
    else if (fiStringCompare("s",time_unit,1)) then
      time_multiplier = 1.d0
    else
      if (grid%myrank == 0) then
        print *, 'Time unit: ', trim(time_unit), ' not recognized'
      endif
      call PetscFinalize()
      stop
    endif

    if (grid%myrank == 0) print *, new_condition%ctype
    if (grid%myrank == 0) print *, new_condition%max_time_index

    count2 = 0
    do i=1,new_condition%max_time_index

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,time,ierr)
      call fiErrorMsg('time',card,ierr)

      call fiReadDouble(string,value,ierr)
      call fiErrorMsg('value',card,ierr)
      ! times are all in seconds   yr -> seconds
      new_condition%times(i) = time*time_multiplier
      new_condition%values(i) = value

      count2 = count2 + 1

    enddo

    if (grid%myrank == 0) print *, count2, ' times read'

    call AddCondition(new_condition)
    nullify(new_condition)

    count1 = count1 + 1

  enddo

  if (grid%myrank == 0) print *, count1, ' conditions read'

! Initial Condition information

  card = "ICON"
  call fiFindStringInFile(fid,card,ierr)

  ! report warning if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'WARNING: Card (',card, ') not found in file'

  icondition = -1
  if (ierr == 0) then
    do
      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadInt(string,icondition,ierr) 
      call fiErrorMsg('icondition',card,ierr)
    enddo
  endif

! set capillary function ids based on material id
  call VecGetArrayF90(grid%icap,icap_p,ierr)
  do local_id=1,grid%nlmax
    if (abs(grid%imat(grid%nL2G(local_id))) > 0) then
      icap_p(local_id) = abs(grid%imat(grid%nL2G(local_id)))
    else
      icap_p(local_id) = 1
    endif
  enddo
  call VecRestoreArrayF90(grid%icap,icap_p,ierr)

  call UpdateGlobalToLocal(grid)
  call InitializeBoundaryConditions(grid)
  if (icondition > 0) call ComputeInitialCondition(grid,icondition)

#ifdef HASH
  deallocate(hash)
#endif
  close(fid)
  call PetscGetTime(time1, ierr)
  time1 = time1 - time0
  if (grid%myrank == 0) print *, time1, ' seconds to read materials data'
 
end subroutine ReadMaterials

! ************************************************************************** !
!
! ReadMaterials2: Reads just the materials from the unstructured input file
! author: Glenn Hammond
! date: 06/20/07
!
! ************************************************************************** !
subroutine ReadMaterials2(grid)

  use fileio_module
  use Condition_module

  implicit none

  type(pflowGrid) :: grid
  character(len=MAXSTRINGLENGTH) :: string 
 
  PetscScalar, pointer :: icap_p(:)

  PetscTruth :: option_found
  integer fid, iostatus, ierr
  integer :: material_id, natural_id, local_id, local_ghosted_id
! integer :: number_of_times
  integer :: count1, count2, i
! integer :: irecv
  real*8 :: xcoord, ycoord, zcoord
  character(len=MAXCARDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: time_unit
  character(len=MAXSTRINGLENGTH) :: filename
! integer :: direction, icond
  integer :: icondition
  real*8 :: time, value, time_multiplier
  PetscLogDouble :: time0, time1
  type(condition_type), pointer :: new_condition

  string = " "
  filename = "hanford.in"
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-unstructured_filename", &
    filename, option_found, ierr)

  call PetscGetTime(time0, ierr)

  fid = 86
  open(fid, file=filename, action="read", status="old", iostat=iostatus)
  if (iostatus /= 0) then
    print *, 'Error opening file: ', trim(filename), ' on Processor ', &
             grid%myrank
    call PetscFinalize()
    stop
  endif

  ! create hash table for fast lookup
#ifdef HASH
  call CreateNaturalToLocalGhostedHash(grid)
!  call PrintHashTable
#endif

! GRID information
  card = "GRID"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif

  count1 = 0
  do

    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,natural_id,ierr) 
    call fiErrorMsg('natural_id',card,ierr)

    call fiReadInt(string,material_id,ierr) 
    call fiErrorMsg('material_id',card,ierr)

    call fiReadDouble(string,xcoord,ierr)
    call fiErrorMsg('xcoord',card,ierr)
   
    call fiReadDouble(string,ycoord,ierr)
    call fiErrorMsg('ycoord',card,ierr)

    call fiReadDouble(string,zcoord,ierr)
    call fiErrorMsg('zcoord',card,ierr)
   
    grid%x(natural_id)=xcoord
    grid%y(natural_id)=ycoord
    grid%z(natural_id)=zcoord

#ifdef HASH
    local_ghosted_id = GetLocalGhostedIdFromHash(natural_id)
#else
    local_ghosted_id = GetLocalGhostedIdFromNaturalId(natural_id,grid)
#endif

    grid%x(local_ghosted_id)=xcoord
    grid%y(local_ghosted_id)=ycoord
    grid%z(local_ghosted_id)=zcoord

    if (local_ghosted_id > 0) then
      local_id = grid%nG2L(local_ghosted_id)
      grid%imat(local_ghosted_id) = material_id
    endif

#if DEBUG
    print *, 'Read geom 1:',natural_id,xcoord,ycoord,zcoord,volume, &
                            xperm,yperm,zperm,porosity,tortuosity
#endif
    count1 = count1 + 1
  enddo

  if (grid%myrank == 0) print *, count1, ' cells read'
  
! CONDition information

  card = "COND"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'ERROR: Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif

  count1 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

! stomp2pflotran header
!:id "type" num_values (datum)

    allocate(new_condition)
    nullify(new_condition%next)
    new_condition%itype = 0
    new_condition%xdatum = 0.d0
    new_condition%ydatum = 0.d0
    new_condition%zdatum = 0.d0
    new_condition%xgrad = 0.d0
    new_condition%ygrad = 0.d0
    new_condition%zgrad = 0.d0

    call fiReadInt(string,new_condition%id,ierr) 
    call fiErrorMsg('id',card,ierr)

    call fiReadQuotedNChars(string,new_condition%ctype,MAXWORDLENGTH, &
                            .true.,ierr) 
    call fiErrorMsg('condition_type',card,ierr)
    call fiWordToLower(new_condition%ctype)

    call fiReadInt(string,new_condition%max_time_index,ierr) 
    call fiErrorMsg('max_time_index',card,ierr)

    call fiReadWord(string,time_unit,.true.,ierr) 
    call fiErrorMsg('time unit',card,ierr)

    if (fiStringCompare("hydraulic gradient",new_condition%ctype,18) .or. &
        fiStringCompare("seepage",new_condition%ctype,7)) then

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,new_condition%xdatum,ierr) 
      call fiErrorMsg('x datum',card,ierr)

      call fiReadDouble(string,new_condition%ydatum,ierr) 
      call fiErrorMsg('y datum',card,ierr)

      call fiReadDouble(string,new_condition%zdatum,ierr) 
      call fiErrorMsg('z datum',card,ierr)

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,new_condition%xgrad,ierr) 
      call fiErrorMsg('x gradient',card,ierr)

      call fiReadDouble(string,new_condition%ygrad,ierr) 
      call fiErrorMsg('y gradient',card,ierr)

      call fiReadDouble(string,new_condition%zgrad,ierr) 
      call fiErrorMsg('z gradient',card,ierr)
    endif

    allocate(new_condition%times(new_condition%max_time_index))
    allocate(new_condition%values(new_condition%max_time_index))
    new_condition%times = 0.d0
    new_condition%values = 0.d0

    call fiWordToLower(time_unit)
    if (fiStringCompare("y",time_unit,1)) then
      time_multiplier = 365.d0*24.d0*3600.d0
    else if (fiStringCompare("d",time_unit,1)) then
      time_multiplier = 24.d0*3600.d0
    else if (fiStringCompare("h",time_unit,1)) then
      time_multiplier = 3600.d0
    else if (fiStringCompare("m",time_unit,1)) then
      time_multiplier = 60.d0
    else if (fiStringCompare("s",time_unit,1)) then
      time_multiplier = 1.d0
    else
      if (grid%myrank == 0) then
        print *, 'Time unit: ', trim(time_unit), ' not recognized'
      endif
      call PetscFinalize()
      stop
    endif

    if (grid%myrank == 0) print *, new_condition%ctype
    if (grid%myrank == 0) print *, new_condition%max_time_index

    count2 = 0
    do i=1,new_condition%max_time_index

      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      call fiReadDouble(string,time,ierr)
      call fiErrorMsg('time',card,ierr)

      call fiReadDouble(string,value,ierr)
      call fiErrorMsg('value',card,ierr)
      ! times are all in seconds   yr -> seconds
      new_condition%times(i) = time*time_multiplier
      new_condition%values(i) = value

      count2 = count2 + 1

    enddo

    if (grid%myrank == 0) print *, count2, ' times read'

    call AddCondition(new_condition)
    nullify(new_condition)

    count1 = count1 + 1

  enddo

  if (grid%myrank == 0) print *, count1, ' conditions read'

! Initial Condition information

  card = "ICON"
  call fiFindStringInFile(fid,card,ierr)

  ! report warning if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'WARNING: Card (',card, ') not found in file'

  icondition = -1
  if (ierr == 0) then
    do
      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(card,ierr)

      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadInt(string,icondition,ierr) 
      call fiErrorMsg('icondition',card,ierr)
    enddo
  endif

! set capillary function ids based on material id
  call VecGetArrayF90(grid%icap,icap_p,ierr)
  do local_id=1,grid%nlmax
    if (abs(grid%imat(grid%nL2G(local_id))) > 0) then
      icap_p(local_id) = abs(grid%imat(grid%nL2G(local_id)))
    else
      icap_p(local_id) = 1
    endif
  enddo
  call VecRestoreArrayF90(grid%icap,icap_p,ierr)

  call InitializeBoundaryConditions(grid)
#if 1
  if (icondition > 0) call ComputeInitialCondition(grid,icondition)
#endif

#ifdef HASH
  deallocate(hash)
#endif
  close(fid)
  call PetscGetTime(time1, ierr)
  time1 = time1 - time0
  if (grid%myrank == 0) print *, time1, ' seconds to read materials data'
 
end subroutine ReadMaterials2

#ifndef USE_HDF5
subroutine ReadStructuredGridHDF5(grid)

  implicit none
  
  type(pflowGrid) :: grid

  if (grid%myrank == 0) then
    print *
    print *, 'PFLOTRAN must be compiled with -DUSE_HDF5 to ', &
             'read HDF5 formatted structured grids.'
    print *
  endif
  stop
  
end subroutine ReadStructuredGridHDF5

#else
! ************************************************************************** !
!
! ReadHDF5StructuredGrid: Reads in a structured grid in HDF5 format
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine ReadStructuredGridHDF5(grid)

  use hdf5
  
  implicit none

  type(pflowGrid) :: grid

  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXSTRINGLENGTH) :: filename
  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id

  integer :: i, local_ghosted_id
  integer, allocatable :: indices(:)
  integer, allocatable :: integer_array(:)
  real*8, allocatable :: real_array(:)
  
  Vec :: global
  Vec :: local
  PetscScalar, pointer :: vec_ptr(:)
      
  PetscTruth :: option_found

  PetscLogDouble :: time0, time1

  call PetscGetTime(time0, ierr)

  filename = "543.h5"
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-hdf5_grid', &
                             filename, option_found, ierr)

  ! create hash table for fast lookup
#ifdef HASH
  call CreateNaturalToLocalGhostedHash(grid)
!  call PrintHashTable
#endif

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  ! open the grid cells group
  string = 'Grid Cells'
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  
  allocate(indices(grid%nlmax))
  
  call SetupCellIndices(grid,grp_id,indices)

  call DACreateGlobalVector(grid%da_1_dof,global,ierr)
  call DACreateLocalVector(grid%da_1_dof,local,ierr)
  
  allocate(integer_array(grid%nlmax))
  string = "Material Id"
  call ReadIntegerArray(grid,grp_id,grid%nlmax,indices,string,integer_array)
  call VecGetArrayF90(global,vec_ptr,ierr)
  do i=1,grid%nlmax
    vec_ptr(i) = integer_array(i)
  enddo
  deallocate(integer_array)
  call VecRestoreArrayF90(global,vec_ptr,ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call DAGlobalToLocalEnd(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call VecGetArrayF90(local,vec_ptr,ierr)
  do i=1,grid%ngmax
    grid%imat(i) = int(vec_ptr(i)+0.0001d0)
  enddo
  call VecRestoreArrayF90(local,vec_ptr,ierr)
  
  
  allocate(real_array(grid%nlmax))
  string = "X-Coordinate"
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,real_array)
  call VecGetArrayF90(global,vec_ptr,ierr)
  do i=1,grid%nlmax
    vec_ptr(i) = real_array(i)
  enddo
  call VecRestoreArrayF90(global,vec_ptr,ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call DAGlobalToLocalEnd(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call VecGetArrayF90(local,vec_ptr,ierr)
  do i=1,grid%ngmax
    grid%x(i) = vec_ptr(i)
  enddo
  call VecRestoreArrayF90(local,vec_ptr,ierr)

  string = "Y-Coordinate"
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,real_array)
  call VecGetArrayF90(global,vec_ptr,ierr)
  do i=1,grid%nlmax
    vec_ptr(i) = real_array(i)
  enddo
  call VecRestoreArrayF90(global,vec_ptr,ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call DAGlobalToLocalEnd(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call VecGetArrayF90(local,vec_ptr,ierr)
  do i=1,grid%ngmax
    grid%y(i) = vec_ptr(i)
  enddo
  call VecRestoreArrayF90(local,vec_ptr,ierr)

  string = "Z-Coordinate"
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,real_array)
  call VecGetArrayF90(global,vec_ptr,ierr)
  do i=1,grid%nlmax
    vec_ptr(i) = real_array(i)
  enddo
  call VecRestoreArrayF90(global,vec_ptr,ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call DAGlobalToLocalEnd(grid%da_1_dof,global,INSERT_VALUES,local,ierr)  
  call VecGetArrayF90(local,vec_ptr,ierr)
  do i=1,grid%ngmax
    grid%z(i) = vec_ptr(i)
  enddo
  call VecRestoreArrayF90(local,vec_ptr,ierr)
  deallocate(real_array)
  
  call VecDestroy(global,ierr)
  call VecDestroy(local,ierr)
  
  
  string = "Volume"
  call VecGetArrayF90(grid%volume,vec_ptr,ierr); CHKERRQ(ierr)
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(grid%volume,vec_ptr,ierr); CHKERRQ(ierr)
  
  string = "X-Permeability"
  call VecGetArrayF90(grid%perm_xx,vec_ptr,ierr); CHKERRQ(ierr)
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(grid%perm_xx,vec_ptr,ierr); CHKERRQ(ierr)
  
  string = "Y-Permeability"
  call VecGetArrayF90(grid%perm_yy,vec_ptr,ierr); CHKERRQ(ierr)
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(grid%perm_yy,vec_ptr,ierr); CHKERRQ(ierr)

  string = "Z-Permeability"
  call VecGetArrayF90(grid%perm_zz,vec_ptr,ierr); CHKERRQ(ierr)
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(grid%perm_zz,vec_ptr,ierr); CHKERRQ(ierr)

  string = "Porosity"
  call VecGetArrayF90(grid%porosity,vec_ptr,ierr); CHKERRQ(ierr)
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(grid%porosity,vec_ptr,ierr); CHKERRQ(ierr)

  string = "Tortuosity"
  call VecGetArrayF90(grid%tor,vec_ptr,ierr); CHKERRQ(ierr)
  call ReadRealArray(grid,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(grid%tor,vec_ptr,ierr); CHKERRQ(ierr)

  ! update local vectors
  
  ! perm_xx
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_xx, INSERT_VALUES, &
                            grid%perm_xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_xx, INSERT_VALUES, &
                          grid%perm_xx_loc, ierr)
  ! perm_yy
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_yy, INSERT_VALUES, &
                            grid%perm_yy_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_yy, INSERT_VALUES, &
                          grid%perm_yy_loc, ierr)
  ! perm_zz
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_zz, INSERT_VALUES, &
                            grid%perm_zz_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_zz, INSERT_VALUES, &
                        grid%perm_zz_loc, ierr)
  ! porosity
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                            grid%porosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                          grid%porosity_loc, ierr)
  ! tortuosity
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                            grid%tor_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                          grid%tor_loc, ierr)

  deallocate(indices)
  call h5gclose_f(grp_id,hdf5_err)
  
  string = 'Connections'
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

  allocate(indices(grid%nconn))

  call SetupConnectionIndices(grid,grp_id,indices)

  allocate(integer_array(grid%nconn))
  string = "Id Upwind"
  call ReadIntegerArray(grid,grp_id,grid%nconn,indices,string,integer_array)
  do i=1,grid%nconn
    local_ghosted_id = GetLocalGhostedIdFromHash(integer_array(i))
    grid%nd1(i) = local_ghosted_id
  enddo

  string = "Id Downwind"
  call ReadIntegerArray(grid,grp_id,grid%nconn,indices,string,integer_array)
  do i=1,grid%nconn
    local_ghosted_id = GetLocalGhostedIdFromHash(integer_array(i))
    grid%nd2(i) = local_ghosted_id
  enddo
  deallocate(integer_array)
  
  string = "Distance Upwind"
  call ReadRealArray(grid,grp_id,grid%nconn,indices,string,grid%dist1) 

  string = "Distance Downwind"
  call ReadRealArray(grid,grp_id,grid%nconn,indices,string,grid%dist2) 

  string = "Area"
  call ReadRealArray(grid,grp_id,grid%nconn,indices,string,grid%area) 

  string = "CosB"
  call ReadRealArray(grid,grp_id,grid%nconn,indices,string,grid%grav_ang) 

  deallocate(indices)

  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
   
  call h5close_f(hdf5_err)
  
  ! set up delz array
  do i=1,grid%nconn
    grid%delz(i) = grid%z(grid%nd2(i))-grid%z(grid%nd1(i))
  enddo
  
end subroutine ReadStructuredGridHDF5

! ************************************************************************** !
!
! SetupCellIndices: Set up indices array that map local cells to entries in
!                   HDF5 grid cell vectors
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine SetupCellIndices(grid,file_id,indices)

  use hdf5
  
  implicit none
  
  type(pflowGrid) :: grid
  
  integer(HID_T) :: file_id
  integer :: indices(:)
  
  character(len=MAXSTRINGLENGTH) :: string 
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: local_ghosted_id, local_id, natural_id
  integer :: index_count
  integer :: cell_count
  integer(HSIZE_T) :: num_cells_in_file
  integer ::temp_int, i
  
  integer, allocatable :: cell_ids(:)
  
  integer :: read_block_size = 100000

  string = "Cell Id"
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_cells_in_file,hdf5_err)
  if (num_cells_in_file /= grid%nmax) then
    if (grid%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_cells_in_file, ') does not match the dimensions of the ', &
               'domain (', grid%nmax, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
  
  allocate(cell_ids(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  cell_count = 0
  index_count = 0
  memory_space_id = -1
  do
    if (cell_count >= num_cells_in_file) exit
    temp_int = min(num_cells_in_file-cell_count,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
    endif
    call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
    ! offset is zero-based
    offset = cell_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset,length, &
                               hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (grid%myrank == 0) then                           
#endif
      call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,cell_ids,dims,hdf5_err, &
                     file_space_id,memory_space_id,prop_id)                     
#ifdef HDF5_BROADCAST
    endif
    if (grid%commsize > 1) &
      call mpi_bcast(cell_ids,dims(1),MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
#endif     
    do i=1,dims(1)
      cell_count = cell_count + 1
      natural_id = cell_ids(i)
      local_ghosted_id = GetLocalGhostedIdFromHash(natural_id)
      if (local_ghosted_id > 0) then
        local_id = grid%nG2L(local_ghosted_id)
        if (local_id > 0) then
          index_count = index_count + 1
          indices(index_count) = cell_count
        endif
      endif
    enddo
  enddo
  
  deallocate(cell_ids)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  if (index_count /= grid%nlmax) then
    if (grid%myrank == 0) &
      print *, 'ERROR: Number of indices read (', index_count, ') does not ', &
               'match the number of local grid cells (', grid%nlmax, ').'
      call PetscFinalize(ierr)
      stop
  endif

end subroutine SetupCellIndices

! ************************************************************************** !
!
! SetupConnectionIndices: Set up indices array that map local connection to  
!                         entriesin HDF5 grid connection vectors
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine SetupConnectionIndices(grid,file_id,indices)

  use hdf5
  
  implicit none
  
  type(pflowGrid) :: grid
  
  integer(HID_T) :: file_id
  integer :: indices(:)
  
  character(len=MAXSTRINGLENGTH) :: string 
  integer(HID_T) :: file_space_id_up, file_space_id_down
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id_up, data_set_id_down
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: local_ghosted_id_up, local_id_up, natural_id_up
  integer :: local_ghosted_id_down, local_id_down, natural_id_down
  integer :: index_count
  integer :: connection_count
  integer(HSIZE_T) :: num_connections_in_file
  integer ::temp_int, i
  
  integer, allocatable :: upwind_ids(:), downwind_ids(:)
  
  integer :: read_block_size = 100000

  string = "Id Upwind"
  call h5dopen_f(file_id,string,data_set_id_up,hdf5_err)
  call h5dget_space_f(data_set_id_up,file_space_id_up,hdf5_err)
  string = "Id Downwind"
  call h5dopen_f(file_id,string,data_set_id_down,hdf5_err)
  call h5dget_space_f(data_set_id_down,file_space_id_down,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id_up, &
                                      num_connections_in_file,hdf5_err)
#if 0
  if (num_connections_in_file /= grid%nconn) then
    if (grid%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_connections_in_file, ') does not match the dimensions ', &
               'of the domain (', grid%nconn, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
#endif
  
  allocate(upwind_ids(read_block_size),downwind_ids(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  connection_count = 0
  index_count = 0
  memory_space_id = -1
  do
    if (connection_count >= num_connections_in_file) exit
    temp_int = min(num_connections_in_file-connection_count,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
    endif
    call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
    ! offset is zero-based
    offset = connection_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id_up, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (grid%myrank == 0) then                           
#endif
      call h5dread_f(data_set_id_up,H5T_NATIVE_INTEGER,upwind_ids,dims, &
                     hdf5_err,file_space_id_up,memory_space_id,prop_id)                     
#ifdef HDF5_BROADCAST
    endif
    if (grid%commsize > 1) &
      call mpi_bcast(upwind_ids,dims(1),MPI_INTEGER,0, &
                     PETSC_COMM_WORLD,ierr)
#endif    
    call h5sselect_hyperslab_f(file_space_id_down, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (grid%myrank == 0) then                           
#endif
      call h5dread_f(data_set_id_down,H5T_NATIVE_INTEGER,downwind_ids,dims, &
                     hdf5_err,file_space_id_down,memory_space_id,prop_id)                     
#ifdef HDF5_BROADCAST
    endif
    if (grid%commsize > 1) &
      call mpi_bcast(downwind_ids,dims(1),MPI_INTEGER,0, &
                     PETSC_COMM_WORLD,ierr)
#endif    
    do i=1,dims(1)
      connection_count = connection_count + 1
      natural_id_up = upwind_ids(i)
      natural_id_down = downwind_ids(i)
      local_ghosted_id_up = GetLocalGhostedIdFromHash(natural_id_up)
      local_ghosted_id_down = GetLocalGhostedIdFromHash(natural_id_down)
      if (local_ghosted_id_up > 0 .and. local_ghosted_id_down > 0) then
        local_id_up = grid%nG2L(local_ghosted_id_up)
        local_id_down = grid%nG2L(local_ghosted_id_down)
        if (local_id_up > 0 .or. local_id_down > 0) then
          index_count = index_count + 1
          indices(index_count) = connection_count
        endif
      endif
    enddo
  enddo
  
  deallocate(upwind_ids,downwind_ids)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id_up,hdf5_err)
  call h5sclose_f(file_space_id_down,hdf5_err)
  call h5dclose_f(data_set_id_up,hdf5_err)
  call h5dclose_f(data_set_id_down,hdf5_err)

  if (index_count /= grid%nconn) then
    if (grid%myrank == 0) &
      print *, 'ERROR: Number of indices read (', index_count, ') does not ', &
               'match the number of local grid connections (', grid%nlmax, ').'
      call PetscFinalize(ierr)
      stop
  endif

end subroutine SetupConnectionIndices

! ************************************************************************** !
!
! ReadRealArray: Read in local real values from hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine ReadRealArray(grid,file_id,num_indices,indices,string,real_array)

  use hdf5
  
  implicit none
  
  type(pflowGrid) :: grid
  
  integer(HID_T) :: file_id
  integer :: num_indices
  integer :: indices(:)
  character(len=MAXSTRINGLENGTH) :: string 
  real*8 :: real_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: index_count
  integer :: real_count
  integer(HSIZE_T) :: num_reals_in_file
  integer :: temp_int, i, index
  
  real*8, allocatable :: real_buffer(:)
  
  integer :: read_block_size = 100000

  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_reals_in_file,hdf5_err)
#if 0
  if (num_reals_in_file /= grid%nmax) then
    if (grid%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_reals_in_file, ') does not match the dimensions of the ', &
               'domain (', grid%nmax, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
#endif
  allocate(real_buffer(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  real_count = 0
  index_count = 0
  memory_space_id = -1
  do i=1,num_indices
    index = indices(i)
    if (index > real_count) then
      temp_int = min(num_reals_in_file-real_count,read_block_size)
      if (dims(1) /= temp_int) then
        if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
        dims(1) = temp_int
      endif
      call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
      ! offset is zero-based
      offset = real_count
      length(1) = dims(1)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                 length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
      if (grid%myrank == 0) then                           
#endif
        call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                       hdf5_err,file_space_id,memory_space_id,prop_id)
#ifdef HDF5_BROADCAST
      endif
      if (grid%commsize > 1) &
        call mpi_bcast(real_buffer,dims(1),MPI_DOUBLE_PRECISION,0, &
                       PETSC_COMM_WORLD,ierr)
#endif
      real_count = real_count + length(1)                  
    endif
    real_array(i) = real_buffer(index)
  enddo

  deallocate(real_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

end subroutine ReadRealArray

! ************************************************************************** !
!
! ReadIntegerArray: Read in local integer values from hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine ReadIntegerArray(grid,file_id,num_indices,indices,string, &
                            integer_array)

  use hdf5
  
  implicit none
  
  type(pflowGrid) :: grid
  
  integer(HID_T) :: file_id
  integer :: num_indices
  integer :: indices(:)
  character(len=MAXSTRINGLENGTH) :: string 
  integer :: integer_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: index_count
  integer :: integer_count
  integer(HSIZE_T) :: num_integers_in_file
  integer :: temp_int, i, index
  
  integer, allocatable :: integer_buffer(:)
  
  integer :: read_block_size = 100000

  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_integers_in_file, &
                                      hdf5_err)
#if 0
  if (num_integers_in_file /= grid%nmax) then
    if (grid%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_integers_in_file, ') does not match the dimensions of ', &
               'the domain (', grid%nmax, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
#endif
  
  allocate(integer_buffer(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  integer_count = 0
  index_count = 0
  memory_space_id = -1

  do i=1,num_indices
    index = indices(i)
    if (index > integer_count) then
      temp_int = min(num_integers_in_file-integer_count,read_block_size)
      if (dims(1) /= temp_int) then
        if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
        dims(1) = temp_int
      endif
      call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
      ! offset is zero-based
      offset = integer_count
      length(1) = dims(1)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                 length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
      if (grid%myrank == 0) then                           
#endif
        call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,integer_buffer,dims, &
                       hdf5_err,file_space_id,memory_space_id,prop_id)   
#ifdef HDF5_BROADCAST
      endif
      if (grid%commsize > 1) &
        call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,0, &
                       PETSC_COMM_WORLD,ierr)
#endif
      integer_count = integer_count + length(1)                  
    endif
    integer_array(i) = integer_buffer(index)
  enddo

  deallocate(integer_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

end subroutine ReadIntegerArray

#endif
! ************************************************************************** !
!
! UpdateGlobalToLocal: Updated global vec values to local
! author: Glenn Hammond
! date: 06/20/07
!
! ************************************************************************** !
subroutine UpdateGlobalToLocal(grid)

  integer :: ierr
  type(pflowGrid) :: grid

  ! icap
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%icap, INSERT_VALUES, &
                            grid%icap_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%icap, INSERT_VALUES, &
                          grid%icap_loc, ierr)
  ! perm_xx
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_xx, INSERT_VALUES, &
                            grid%perm_xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_xx, INSERT_VALUES, &
                          grid%perm_xx_loc, ierr)
  ! perm_yy
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_yy, INSERT_VALUES, &
                            grid%perm_yy_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_yy, INSERT_VALUES, &
                          grid%perm_yy_loc, ierr)
  ! perm_zz
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_zz, INSERT_VALUES, &
                            grid%perm_zz_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_zz, INSERT_VALUES, &
                          grid%perm_zz_loc, ierr)
  ! tor
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                            grid%tor_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                          grid%tor_loc, ierr)
  ! por
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                            grid%porosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                          grid%porosity_loc, ierr)

end subroutine UpdateGlobalToLocal

! ************************************************************************** !
!
! GetNumberOfLocalConnection: Returns the number of connections involving
!                             local (ghosted) cells on the local proc
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
integer function GetNumberOfLocalConnections(fid,grid)

  use fileio_module

  implicit none

  type(pflowGrid) :: grid

  integer :: fid, ierr, num_local_connections
  integer :: connection_id
  integer :: natural_id_upwind, natural_id_downwind
  integer :: local_ghosted_id_upwind, local_ghosted_id_downwind
  integer :: local_id_upwind, local_id_downwind
  character(len=MAXCARDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: string

! CONNection information
  card = "CONN"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then 
    if (grid%myrank == 0) print *, 'Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif

  num_local_connections = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,connection_id,ierr) 
    call fiErrorMsg('connection_id',card,ierr)

    call fiReadInt(string,natural_id_upwind,ierr) 
    call fiErrorMsg('natural_id_upwind',card,ierr)

    call fiReadInt(string,natural_id_downwind,ierr) 
    call fiErrorMsg('natural_id_downwind',card,ierr)

#ifdef HASH
    local_ghosted_id_upwind = GetLocalGhostedIdFromHash(natural_id_upwind)
    local_ghosted_id_downwind = GetLocalGhostedIdFromHash(natural_id_downwind)
#else
    local_ghosted_id_upwind = &
               GetLocalGhostedIdFromNaturalId(natural_id_upwind,grid)
    local_ghosted_id_downwind = &
               GetLocalGhostedIdFromNaturalId(natural_id_downwind,grid)
#endif
    if (local_ghosted_id_upwind > 0 .and. &
        local_ghosted_id_downwind > 0) then ! both cells are local ghosted
      local_id_upwind = grid%nG2L(local_ghosted_id_upwind)
      local_id_downwind = grid%nG2L(local_ghosted_id_downwind)
      ! at least 1 cell must be local (non-ghosted) since we don't want to
      ! create a connection between two ghosted cells
      if (local_id_upwind > 0 .or. local_id_downwind > 0) then
        num_local_connections = num_local_connections + 1
      endif
    endif

  enddo

  GetNumberOfLocalConnections =  num_local_connections

end function GetNumberOfLocalConnections

! ************************************************************************** !
!
! GetNumberOfBoundaryConnections: Returns the number of connections involving
!                                 boundary cells on the local proc
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
integer function GetNumberOfBoundaryConnections(fid,grid)

  use fileio_module

  implicit none

  type(pflowGrid) :: grid

  integer :: fid, ierr, num_boundary_connections
  integer :: connection_id
  integer :: natural_id, local_id, local_ghosted_id
  character(len=MAXCARDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: string

! BCONnection information
  card = "BCON"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0) then
    if (grid%myrank == 0) print *, 'Card (',card, ') not found in file'
    call PetscFinalize()
    stop
  endif

  num_boundary_connections = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,connection_id,ierr) 
    call fiErrorMsg('connection_id',card,ierr)

    call fiReadInt(string,natural_id,ierr) 
    call fiErrorMsg('natural_id',card,ierr)

#ifdef HASH
    local_ghosted_id = GetLocalGhostedIdFromHash(natural_id)
#else
    local_ghosted_id = GetLocalGhostedIdFromNaturalId(natural_id,grid)
#endif
    if (local_ghosted_id > 0) then 
      local_id = grid%nG2L(local_ghosted_id)
      if (local_id > 0) then
        num_boundary_connections = num_boundary_connections + 1
      endif
    endif
  enddo

  GetNumberOfBoundaryConnections = num_boundary_connections

end function GetNumberOfBoundaryConnections
! ************************************************************************** !
!
! GetLocalIdFromNaturalId: Returns the local id corresponding to a natural
!                          id or 0, if the natural id is off-processor
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
integer function GetLocalIdFromNaturalId(natural_id,grid)

  implicit none

  type(pflowGrid) :: grid

  integer :: natural_id, local_id
  
  do local_id = 1, grid%nlmax
    if (natural_id == grid%nL2A(local_id)+1) then
      GetLocalIdFromNaturalId = local_id
      return
    endif
  enddo
  GetLocalIdFromNaturalId = 0

end function GetLocalIdFromNaturalId

! ************************************************************************** !
!
! GetLocalGhostedIdFromNaturalId: Returns the local ghosted id corresponding 
!                                 to a natural id or 0, if the natural id 
!                                 is off-processor
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
integer function GetLocalGhostedIdFromNaturalId(natural_id,grid)

  implicit none

  type(pflowGrid) :: grid

  integer :: natural_id, local_ghosted_id
  
  do local_ghosted_id = 1, grid%ngmax
    if (natural_id == grid%nG2A(local_ghosted_id)+1) then
      GetLocalGhostedIdFromNaturalId = local_ghosted_id
      return 
    endif
  enddo
  GetLocalGhostedIdFromNaturalId = 0

end function GetLocalGhostedIdFromNaturalId

! ************************************************************************** !
!
! CreateNaturalToLocalHash: Creates a hash table for looking up the local 
!                           ghosted id of a natural id, if it exists
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine CreateNaturalToLocalGhostedHash(grid)

  implicit none

  type(pflowGrid) :: grid

  integer :: local_ghosted_id, natural_id 
  integer :: num_in_hash, num_ids_per_hash, hash_id, id
  integer, pointer :: temp_hash(:,:,:)

  ! initial guess of 10% of ids per hash
  num_ids_per_hash = max(grid%nlmax/(num_hash/10),grid%nlmax)

  allocate(hash(2,0:num_ids_per_hash,num_hash))
  hash(:,:,:) = 0

  ! recall that natural ids are zero-based
  do local_ghosted_id = 1, grid%ngmax
    natural_id = grid%nG2A(local_ghosted_id)+1
    hash_id = mod(natural_id,num_hash)+1 
    num_in_hash = hash(1,0,hash_id)
    num_in_hash = num_in_hash+1
    ! if a hash runs out of space reallocate
    if (num_in_hash > num_ids_per_hash) then 
      allocate(temp_hash(2,0:num_ids_per_hash,0:num_hash))
      ! copy old hash
      temp_hash(1:2,0:num_ids_per_hash,num_hash) = &
                             hash(1:2,0:num_ids_per_hash,num_hash)
      deallocate(hash)
      ! recompute hash 20% larger
      num_ids_per_hash = int(dble(num_ids_per_hash)*1.2)
      allocate(hash(1:2,0:num_ids_per_hash,num_hash))
      ! copy old to new
      do hash_id = 1, num_hash
        do id = 1, temp_hash(1,0,hash_id)
          hash(1:2,id,hash_id) = temp_hash(1:2,id,hash_id)
        enddo
        hash(1,0,hash_id) = temp_hash(1,0,hash_id)
      enddo
      deallocate(temp_hash)
    endif
    hash(1,0,hash_id) = num_in_hash
    hash(1,num_in_hash,hash_id) = natural_id
    hash(2,num_in_hash,hash_id) = local_ghosted_id
  enddo

  if (grid%myrank == 0) print *, 'num_ids_per_hash:', num_ids_per_hash

end subroutine CreateNaturalToLocalGhostedHash

! ************************************************************************** !
!
! GetLocalIdFromHash: Returns the local ghosted id of a natural id, if it 
!                     exists.  Otherwise 0 is returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
integer function GetLocalGhostedIdFromHash(natural_id)

  implicit none

  integer :: natural_id
  integer :: hash_id, id

  GetLocalGhostedIdFromHash = 0
  hash_id = mod(natural_id,num_hash)+1 
  do id = 1, hash(1,0,hash_id)
    if (hash(1,id,hash_id) == natural_id) then
      GetLocalGhostedIdFromHash = hash(2,id,hash_id)
      return
    endif
  enddo

end function GetLocalGhostedIdFromHash

! ************************************************************************** !
!
! PrintHashTable: Prints the hashtable for viewing
! author: Glenn Hammond
! date: 03/09/07
!
! ************************************************************************** !
subroutine PrintHashTable()

  implicit none

  integer :: ihash, id, fid

  fid = 87 
  open(fid,file='hashtable.dat',action='write')
  do ihash=1,num_hash
    write(fid,'(a4,i3,a,i5,a2,x,200(i6,x))') 'Hash',ihash,'(', &
                         hash(1,0,ihash), &
                         '):', &
                         (hash(1,id,ihash),id=1,hash(1,0,ihash))
  enddo
  close(fid)

end subroutine PrintHashTable

end module Unstructured_Grid_module
