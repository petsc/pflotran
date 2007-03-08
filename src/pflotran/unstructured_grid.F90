module unstructured_grid_mod

 use pflow_gridtype_module

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

  implicit none

  type(pflowGrid) :: grid
  character(len=MAXSTRINGLENGTH) :: string 
 
  PetscScalar, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:), &
                          tor_p(:), volume_p(:)
  integer INUNIT1,ierr
  integer :: connection_id, material_id, natural_id, local_id, local_ghosted_id
  integer :: natural_id_upwind, natural_id_downwind
  integer :: local_id_upwind, local_id_downwind
  integer :: local_ghosted_id_upwind, local_ghosted_id_downwind
  real*8 :: xcoord, ycoord, zcoord, xperm, yperm, zperm
  real*8 :: area, distance, distance_upwind, distance_downwind, cosB
  real*8 :: volume, porosity, tortuosity, direction
  real*8 :: dx, dy, dz
  character*4 card

  call VecGetArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr)

  open(IUNIT1, file="geom_field.in", action="read", status="old")
 
  ! create hash table for fast lookup
  call CreateNaturalToLocalHash(grid)

! GRID information
  card = "GRID"
  call fiFindStringInFile(INUNIT1,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'Card (',card, ') not found in file'

  do
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,natural_id,ierr) 
    call fiDefaultMsg('natural_id',ierr)

    call fiReadInt(string,material_id,ierr) 
    call fiDefaultMsg('material_id',ierr)

    call fiReadDouble(string,xcoord,ierr)
    call fiDefaultMsg('xcoord',ierr)
   
    call fiReadDouble(string,ycoord,ierr)
    call fiDefaultMsg('ycoord',ierr)

    call fiReadDouble(string,zcoord,ierr)
    call fiDefaultMsg('zcoord',ierr)
   
    call fiReadDouble(string,volume,ierr)
    call fiDefaultMsg('volume',ierr)
   
    call fiReadDouble(string,xperm,ierr)
    call fiDefaultMsg('xperm',ierr)

    call fiReadDouble(string,yperm,ierr)
    call fiDefaultMsg('yperm',ierr)
   
    call fiReadDouble(string,zperm,ierr)
    call fiDefaultMsg('zperm',ierr)

    call fiReadDouble(string,porosity,ierr)
    call fiDefaultMsg('porosity',ierr)

    call fiReadDouble(string,tortuosity,ierr)
    call fiDefaultMsg('tortuosity',ierr)

!geh In the future, 'natural_id' will be replaced by 'local_id'
    grid%x(natural_id)=xcoord
    grid%y(natural_id)=ycoord
    grid%z(natural_id)=zcoord

    local_ghosted_id = GetLocalIdFromHash(natural_id)
    if (local_ghosted_id > 0) then
      local_id = grid%nG2L(local_ghosted_id)
      if (local_id > 0) then
        perm_xx_p(local_id) = xperm
        perm_yy_p(local_id) = yperm
        perm_zz_p(local_id) = zperm
        volume_p(local_id) = volume
        por_p(local_id) = porosity
        tor_p(local_id) = tortuosity
        ! need something in which to store material id
      endif
    endif

#if DEBUG
    print *, 'Read geom 1:',natural_id,xcoord,ycoord,zcoord,volume, &
                            xperm,yperm,zperm,porosity,tortuosity
#endif
  enddo
  
! CONNection information

!geh Right now we ignore allocation since the arrays are already allocated
!geh larger than necessary.  This will need to be updated in the future.
! Count the number of local connections
!geh  grid%nconn = GetNumberOfLocalConnections(fid) ! function at bottom of file
! Allocate Arrays
!geh  allocate(grid%nd1(grid%nconn))
!geh  allocate(grid%nd2(grid%nconn))
!geh  allocate(grid%dist1(grid%nconn))
!geh  allocate(grid%dist2(grid%nconn))
!geh  allocate(grid%area(grid%nconn))
!geh  allocate(grid%delz(grid%nconn))
!geh  allocate(grid%grav_ang(grid%nconn))
!geh
!geh  allocate(grid%iperm1(grid%nconn))
!geh  allocate(grid%iperm2(grid%nconn))
!geh
!geh  allocate(grid%vl_loc(grid%nconn))
!geh  allocate(grid%vvl_loc(grid%nconn))
!geh  allocate(grid%vg_loc(grid%nconn))
!geh  allocate(grid%vvg_loc(grid%nconn))

  grid%nconn = 0

  card = "CONN"
  call fiFindStringInFile(INUNIT1,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'Card (',card, ') not found in file'

  do
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

! :id idup iddn distup distdn area cosB
    call fiReadInt(string,connection_id,ierr) 
    call fiDefaultMsg('connection_id',ierr)

    call fiReadInt(string,natural_id_upwind,ierr) 
    call fiDefaultMsg('natural_id_upwind',ierr)

    call fiReadInt(string,natural_id_downwind,ierr) 
    call fiDefaultMsg('natural_id_downwind',ierr)

    local_ghosted_id_upwind = GetLocalIdFromHash(natural_id_upwind)
    local_ghosted_id_downwind = GetLocalIdFromHash(natural_id_downwind)
    if (local_ghosted_id_upwind > 0 .and. &
        local_ghosted_id_downwind > 0) then ! both cells are local ghosted
      local_id_upwind = grid%nG2L(local_ghosted_id_upwind)
      local_id_downwind = grid%nG2L(local_ghosted_id_downwind)
      ! at least 1 cell must be local (non-ghosted) since we don't want to
      ! create a connection between two ghosted cells
      if (local_id_upwind > 0 .or. local_id_downwind > 0) then

        ! Continue reading string if local
        call fiReadDouble(string,distance_upwind,ierr) 
        call fiDefaultMsg('distance_upwind',ierr)

        call fiReadDouble(string,distance_downwind,ierr) 
        call fiDefaultMsg('distance_downwind',ierr)

        call fiReadDouble(string,area,ierr) 
        call fiDefaultMsg('area',ierr)

        call fiReadDouble(string,cosB,ierr) 
        call fiDefaultMsg('cosB',ierr)

        grid%nconn = grid%nconn + 1
        grid%nd1(grid%nconn) = local_ghosted_id_upwind
        grid%nd2(grid%nconn) = local_ghosted_id_downwind
        grid% dist1(grid%nconn) = distance_upwind
        grid% dist2(grid%nconn) = distance_downwind
        if (area>0.D0) grid%area(grid%nconn)= area
        grid%delz(grid%nconn) = abs(grid%z(natural_id_downwind)-  &
                                    grid%z(natural_id_upwind))
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
  enddo


! BCONection information

!geh Right now we ignore allocation since the arrays are already allocated
!geh larger than necessary.  This will need to be updated in the future.
! Count the number of boundary connections
!geh  grid%nconn = GetNumberOfBoundary(fid) ! function at bottom of file
! Allocate Arrays
!geh    allocate(grid%mblkbc(grid%nconnbc))
!geh    allocate(grid%ibconn(grid%nconnbc))
!geh    allocate(grid%distbc(grid%nconnbc))
!geh    allocate(grid%areabc(grid%nconnbc))
!geh    allocate(grid%ipermbc(grid%nconnbc))
!geh    allocate(grid%delzbc(grid%nconnbc))

!    allocate(grid%velocitybc(grid%nphase,grid%nconnbc))

!geh    allocate(grid%vlbc(grid%nconnbc))
!geh    allocate(grid%vvlbc(grid%nconnbc))
!geh    allocate(grid%vgbc(grid%nconnbc))
!geh    allocate(grid%vvgbc(grid%nconnbc))

  grid%nconnbc = 0

  card = "BCON"
  call fiFindStringInFile(INUNIT1,card,ierr)

  ! report error if card does not exist
  if (ierr /= 0 .and. grid%myrank == 0) &
    print *, 'Card (',card, ') not found in file'

  do
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

!:id cellid dist area cosB "type" time value
    call fiReadInt(string,connection_id,ierr) 
    call fiDefaultMsg('connection_id',ierr)

    call fiReadInt(string,natural_id,ierr) 
    call fiDefaultMsg('natural_id',ierr)

    local_ghosted_id = GetLocalIdFromHash(natural_id)
    if (local_ghosted_id > 0) then 
      local_id = grid%nG2L(local_ghosted_id)
      if (local_id > 0) then

        ! Continue reading string if local
        call fiReadDouble(string,direction,ierr) 
        call fiDefaultMsg('direction',ierr)

        call fiReadDouble(string,distance,ierr) 
        call fiDefaultMsg('distance',ierr)

        call fiReadDouble(string,area,ierr) 
        call fiDefaultMsg('area',ierr)

        call fiReadDouble(string,cosB,ierr) 
        call fiDefaultMsg('cosB',ierr)

        grid%distbc(grid%nconnbc) = distance
        grid%areabc(grid%nconnbc) = area

        ! setup direction of permeability
        grid%ipermbc(grid%nconnbc) = abs(direction) ! 1,2,3 -> x,y,z

        if (direction == -3) then ! bottom face
          grid%delzbc(grid%nconnbc) = distance 
        else if (direction == 3) then ! top face
          grid%delzbc(grid%nconnbc) = -1.d0*distance
        else
          grid%delzbc(grid%nconnbc) = 0.d0
        endif
      endif
    endif
  enddo
  deallocate(hash)
  close(IUNIT1)
 
  call VecRestoreArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)     
  call VecRestoreArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)              
  call VecRestoreArrayF90(grid%tor,tor_p,ierr)
  call VecRestoreArrayF90(grid%porosity,por_p,ierr)

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
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg('GRID',ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,connection_id,ierr) 
    call fiDefaultMsg('connection_id',ierr)

    call fiReadInt(string,natural_id_upwind,ierr) 
    call fiDefaultMsg('natural_id_upwind',ierr)

    call fiReadInt(string,natural_id_downwind,ierr) 
    call fiDefaultMsg('natural_id_downwind',ierr)

    local_ghosted_id_upwind = GetLocalIdFromHash(natural_id_upwind)
    local_ghosted_id_downwind = GetLocalIdFromHash(natural_id_downwind)
    if (local_ghosted_id_upwind > -1 .and. &
        local_ghosted_id_downwind > -1) then ! both cells are local ghosted
      local_id_upwind = grid%nG2L(local_ghosted_id_upwind)
      local_id_downwind = grid%nG2L(local_ghosted_id_downwind)
      ! at least 1 cell must be local (non-ghosted) since we don't want to
      ! create a connection between two ghosted cells
      if (local_id_upwind > -1 .or. local_id_downwind > -1) then
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
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg('GRID',ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    call fiReadInt(string,connection_id,ierr) 
    call fiDefaultMsg('connection_id',ierr)

    call fiReadInt(string,natural_id,ierr) 
    call fiDefaultMsg('natural_id',ierr)

    local_ghosted_id = GetLocalIdFromHash(natural_id)
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
    if (natural_id == grid%nL2A(local_id)) then
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
    if (natural_id == grid%nL2A(grid%nG2L(local_ghosted_id))) then
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
subroutine CreateNaturalToLocalHash(grid)

  implicit none

  type(pflowGrid) :: grid

  integer :: local_ghosted_id, natural_id 
  integer :: num_in_hash, num_ids_per_hash, hash_id, id
  integer, pointer :: temp_hash(:,:,:)

  ! initial guess of 10% of ids per hash
  num_ids_per_hash = grid%nlmax/(num_hash/10)

  allocate(hash(2,0:num_ids_per_hash,num_hash))
  hash(:,:,:) = 0

  do local_ghosted_id = 1, grid%ngmax
    natural_id = grid%nL2A(local_ghosted_id)
    hash_id = mod(natural_id,num_hash)
    num_in_hash = hash(1,0,hash_id)+1
    ! if a hash runs out of space reallocate
    if (num_in_hash > num_ids_per_hash) then 
      allocate(temp_hash(2,0:num_ids_per_hash,num_hash))
      ! copy old hash
      temp_hash(1:2,0:num_ids_per_hash,1:num_hash) = &
                             hash(1:2,0:num_ids_per_hash,1:num_hash)
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

end subroutine CreateNaturalToLocalHash

! ************************************************************************** !
!
! GetLocalIdFromHash: Returns the local ghosted id of a natural id, if it 
!                     exists.  Otherwise 0 is returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
integer function GetLocalIdFromHash(natural_id)

  implicit none

  integer :: natural_id
  integer :: hash_id, id

  hash_id = mod(natural_id,num_hash)
  do id = 1, hash(1,0,hash_id)
    if (hash(1,id,hash_id) == natural_id) then
      GetLocalIdFromHash = hash(2,id,hash_id)
      return
    endif
  enddo
  GetLocalIdFromHash = 0

end function GetLocalIdFromHash

end module unstructured_grid_mod
