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
   
  public ReadUnstructuredGrid
  
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
  integer fid, ierr
  integer :: connection_id, material_id, natural_id, local_id, local_ghosted_id
  integer :: natural_id_upwind, natural_id_downwind
  integer :: local_id_upwind, local_id_downwind
  integer :: local_ghosted_id_upwind, local_ghosted_id_downwind
  integer :: number_of_times
  integer :: count1, count2, i
  real*8 :: xcoord, ycoord, zcoord, xperm, yperm, zperm
  real*8 :: area, distance, distance_upwind, distance_downwind, cosB
  real*8 :: volume, porosity, tortuosity
  real*8 :: dx, dy, dz, delz
  character(len=MAXCARDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: bctype
  character(len=MAXWORDLENGTH) :: time_unit
  integer :: direction, icond, icondition
  real*8 :: time, value
  real*8 :: time0, time1, time_multiplier
  type(condition_type), pointer :: new_condition

  call VecGetArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr)

  time0 = secnds(0.)
  fid = 86
  open(fid, file="hanford.in", action="read", status="old")

  ! create hash table for fast lookup
#ifdef HASH
  call CreateNaturalToLocalGhostedHash(grid)
  call PrintHashTable
#endif

! GRID information
  card = "GRID"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) then
    print *, 'ERROR: Card (',card, ') not found in file'
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

!geh In the future, 'natural_id' will be replaced by 'local_id'
    grid%x(natural_id)=xcoord
    grid%y(natural_id)=ycoord
    grid%z(natural_id)=zcoord

#ifdef HASH
    local_ghosted_id = GetLocalGhostedIdFromHash(natural_id)
#else
    local_ghosted_id = GetLocalGhostedIdFromNaturalId(natural_id,grid)
#endif
    if (local_ghosted_id > 0) then
      local_id = grid%nG2L(local_ghosted_id)
      if (local_id > 0) then
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
  call VecRestoreArrayF90(grid%tor,tor_p,ierr)
  call VecRestoreArrayF90(grid%porosity,por_p,ierr)

  print *, count1, ' cells read'
  
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
  if (ierr /= 0 .and. grid%myrank == 0) then
    print *, 'ERROR: Card (',card, ') not found in file'
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
        grid%delz(grid%nconn) = -1.d0*abs(grid%z(natural_id_downwind)-  &
                                    grid%z(natural_id_upwind))
#else
        grid%delz(grid%nconn) = +1.d0*abs(grid%z(natural_id_downwind)-  &
                                    grid%z(natural_id_upwind))
#endif
        grid%grav_ang(grid%nconn) = cosB

        ! setup direction of permeability
        dx = abs(grid%x(natural_id_upwind)-grid%x(natural_id_downwind))
        dy = abs(grid%y(natural_id_upwind)-grid%y(natural_id_downwind))
        dz = abs(grid%z(natural_id_upwind)-grid%z(natural_id_downwind))
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
  enddo
  print *, count2, ' of ', count1, ' connections found'


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
  if (ierr /= 0 .and. grid%myrank == 0) then
    print *, 'ERROR: Card (',card, ') not found in file'
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
    local_ghosted)id = GetLocalGhostedIdFromNaturalId(natural_id,grid)
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

  print *, count2, ' of ', count1, ' boundary connections found'

! CONDition information

  card = "COND"
  call fiFindStringInFile(fid,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) then
    print *, 'ERROR: Card (',card, ') not found in file'
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

    print *, count2, ' times read'

    call AddCondition(new_condition)
    nullify(new_condition)

    count1 = count1 + 1

  enddo

  print *, count1, ' conditions read'

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
    icap_p(local_id) = grid%imat(grid%nL2G(local_id))
  enddo
  call VecRestoreArrayF90(grid%icap,icap_p,ierr)

  call InitializeBoundaryConditions(grid)
  if (icondition > 0) call ComputeInitialCondition(grid,icondition)

#ifdef HASH
  deallocate(hash)
#endif
  close(fid)
  time1 = secnds(0.) - time0
  print *, time1, ' seconds to read unstructured grid data'
 
end subroutine ReadUnstructuredGrid
 
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
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'Card (',card, ') not found in file'

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
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'Card (',card, ') not found in file'

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
    if (natural_id == grid%nL2A(grid%nG2L(local_ghosted_id))+1) then
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
    natural_id = grid%nL2A(local_ghosted_id)+1
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

  print *, 'num_ids_per_hash:', num_ids_per_hash

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
