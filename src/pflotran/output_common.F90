module Output_Common_module

  use Logging_module 
  use Output_Aux_module
  
  !note: only realization_base_type can be used throughout this module.
  use Realization_Base_class, only : realization_base_type

  implicit none

  private

#include "definitions.h"

  PetscInt, save, public :: max_local_size_saved = -1

  !geh: would prefer that this be local to Output_Tecplot_module, but needed
  !     in Output_Surface_module
  PetscInt, parameter, public :: TECPLOT_INTEGER = 0
  PetscInt, parameter, public :: TECPLOT_REAL = 1  
  
  public :: OutputCommonInit, &
            OutputGetVarFromArray, &
            OutputGetVarFromArrayAtCoord, &
            OutputGetCellCenteredVelocities, &
            ConvertArrayToNatural, &
            OutputFormatInt, &
            OutputFormatDouble, &
            GetCellCoordinates, &
            GetVertexCoordinates, &
            OutputFilenameID, &
            OutputFilename, &
            GetCellConnections, &
            OutputXMFHeader, &
            OutputXMFAttribute, &
            OutputXMFFooter, &
            OutputGetFlowrates, &
            ExplicitGetCellCoordinates, &
            OutputGetExplicitFlowrates
              
contains

! ************************************************************************** !
!
! OutputCommonInit: Initializes module variables for common formats
! author: Glenn Hammond
! date: 01/16/13
!
! ************************************************************************** !
subroutine OutputCommonInit()

  use Option_module

  implicit none
  
  ! set size to -1 in order to re-initialize parallel communication blocks
  max_local_size_saved = -1

end subroutine OutputCommonInit

! ************************************************************************** !
!
! OutputFilenameID: Creates an ID for filename
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFilenameID(output_option,option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option

  character(len=MAXWORDLENGTH) :: OutputFilenameID
  
  if (output_option%plot_number < 10) then
    write(OutputFilenameID,'("00",i1)') output_option%plot_number
  else if (output_option%plot_number < 100) then
    write(OutputFilenameID,'("0",i2)') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(OutputFilenameID,'(i3)') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(OutputFilenameID,'(i4)') output_option%plot_number  
  endif 
  
  OutputFilenameID = adjustl(OutputFilenameID)

end function OutputFilenameID

! ************************************************************************** !
!
! OutputFilename: Creates a filename for a Tecplot file
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFilename(output_option,option,suffix,optional_string)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option
  character(len=*) :: suffix
  character(len=*) :: optional_string
  
  character(len=MAXSTRINGLENGTH) :: OutputFilename

  character(len=MAXWORDLENGTH) :: final_suffix
  character(len=MAXSTRINGLENGTH) :: final_optional_string


  if (len_trim(optional_string) > 0) then
    final_optional_string = '-' // optional_string
  else
    final_optional_string = ''
  endif
  final_suffix = '.' // suffix
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    OutputFilename = trim(output_option%plot_name) // &
            trim(final_optional_string) // &
            final_suffix
  else  
    OutputFilename = trim(option%global_prefix) // &
            trim(option%group_prefix) // &
            trim(final_optional_string) // &
            '-' // &
            trim(OutputFilenameID(output_option,option)) // &
            final_suffix
  endif
  
end function OutputFilename

! ************************************************************************** !
!
! OutputGetVarFromArray: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputGetVarFromArray(realization_base,vec,ivar,isubvar,isubvar1)

  use Realization_Base_class, only : RealizationGetDataset
  use Grid_module
  use Option_module
  use Field_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_output_get_var_from_array,ierr) 
                        
  call RealizationGetDataset(realization_base,vec,ivar,isubvar,isubvar1)

  call PetscLogEventEnd(logging%event_output_get_var_from_array,ierr) 
  
end subroutine OutputGetVarFromArray

! ************************************************************************** !
!
! ConvertArrayToNatural: Converts an array  to natural ordering
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine ConvertArrayToNatural(indices,array,local_size,global_size,option)

  use Option_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  PetscInt :: local_size, global_size
  PetscInt :: indices(:)
  PetscReal, pointer :: array(:)
  type(option_type) :: option
  
  Vec :: natural_vec
  PetscInt, allocatable :: indices_zero_based(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecCreate(option%mycomm,natural_vec,ierr)
  call VecSetSizes(natural_vec,PETSC_DECIDE,global_size,ierr)
  call VecSetType(natural_vec,VECMPI,ierr)

  allocate(indices_zero_based(local_size))
  indices_zero_based(1:local_size) = indices(1:local_size)-1

  call VecSetValues(natural_vec,local_size,indices_zero_based, &
                    array,INSERT_VALUES,ierr)

  call VecAssemblyBegin(natural_vec,ierr)
  call VecAssemblyEnd(natural_vec,ierr)

  call VecGetLocalSize(natural_vec,local_size,ierr)
  deallocate(array)
  allocate(array(local_size))
  
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)
  array(1:local_size) = vec_ptr(1:local_size)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)

  call VecDestroy(natural_vec,ierr)
  
end subroutine ConvertArrayToNatural

! ************************************************************************** !
!
! OutputGetVarFromArrayAtCoord: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function OutputGetVarFromArrayAtCoord(realization_base,ivar,isubvar,x,y,z, &
                                      num_cells,ghosted_ids,isubvar1)

  use Realization_Base_class, only : RealizGetDatasetValueAtCell
  use Grid_module
  use Option_module

  implicit none
  
  PetscReal :: OutputGetVarFromArrayAtCoord
  class(realization_base_type) :: realization_base
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscReal :: x,y,z
  PetscInt :: num_cells
  PetscInt :: ghosted_ids(num_cells)

  type(grid_type), pointer :: grid
  PetscInt :: icell, ghosted_id
  PetscReal :: dx, dy, dz
  PetscReal :: value, sum_value
  PetscReal :: weight, sum_weight, sum_root
  
  sum_value = 0.d0
  sum_weight = 0.d0
  
  grid => realization_base%patch%grid

  do icell=1, num_cells
    ghosted_id = ghosted_ids(icell)
    dx = x-grid%x(ghosted_id)
    dy = y-grid%y(ghosted_id)
    dz = z-grid%z(ghosted_id)
    sum_root = sqrt(dx*dx+dy*dy+dz*dz)
    value = 0.d0
    value = RealizGetDatasetValueAtCell(realization_base,ivar,isubvar,ghosted_id, &
      isubvar1)
    if (sum_root < 1.d-40) then ! bail because it is right on this coordinate
      sum_weight = 1.d0
      sum_value = value
      exit
    endif
    weight = 1.d0/sum_root
    sum_weight = sum_weight + weight
    sum_value = sum_value + weight * value
  enddo
  
  OutputGetVarFromArrayAtCoord = sum_value/sum_weight

end function OutputGetVarFromArrayAtCoord

! ************************************************************************** !
!
! OutputGetCellCenteredVelocities: Computes the cell-centered velocity component 
!                            as an averages of cell face velocities
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputGetCellCenteredVelocities(realization_base,vec,iphase,direction)

  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Patch_module
  use Logging_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: direction
  PetscInt :: iphase
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscReal :: velocity, area
  PetscReal :: average, sum, max, min, std_dev
  PetscInt :: max_loc, min_loc
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, allocatable :: sum_area(:)
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  call PetscLogEventBegin(logging%event_output_get_cell_vel,ierr) 
                            
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
    
  allocate(sum_area(grid%nlmax))
  sum_area(1:grid%nlmax) = 0.d0

  call VecSet(vec,0.d0,ierr)
  call VecGetArrayF90(vec,vec_ptr,ierr)

  ! interior velocities  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(direction,iconn)
      velocity = patch%internal_velocities(iphase,sum_connection)* &
                 area
      if (local_id_up > 0) then
        vec_ptr(local_id_up) = vec_ptr(local_id_up) + velocity
        sum_area(local_id_up) = sum_area(local_id_up) + dabs(area)
      endif
      if (local_id_dn > 0) then
        vec_ptr(local_id_dn) = vec_ptr(local_id_dn) + velocity
        sum_area(local_id_dn) = sum_area(local_id_dn) + dabs(area)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! boundary velocities
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(direction,iconn)
      vec_ptr(local_id) = vec_ptr(local_id)+ &
                          patch%boundary_velocities(iphase,sum_connection)* &
                          area
      sum_area(local_id) = sum_area(local_id) + dabs(area)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! divide by total area
  do local_id=1,grid%nlmax
    if (sum_area(local_id) > 0.d0) &
      vec_ptr(local_id) = vec_ptr(local_id)/sum_area(local_id)*output_option%tconv
  enddo
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  deallocate(sum_area)

  call PetscLogEventEnd(logging%event_output_get_cell_vel,ierr) 

end subroutine OutputGetCellCenteredVelocities

! ************************************************************************** !
!
! OutputFormatInt: Writes a integer to a string
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFormatInt(int_value)

  implicit none
  
  PetscInt :: int_value
  
  character(len=MAXWORDLENGTH) :: OutputFormatInt

  write(OutputFormatInt,'(1i12)') int_value
  
  OutputFormatInt = adjustl(OutputFormatInt)
  
end function OutputFormatInt

! ************************************************************************** !
!
! OutputFormatDouble: Writes a double or real to a string
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFormatDouble(real_value)

  implicit none
  
  PetscReal :: real_value
  
  character(len=MAXWORDLENGTH) :: OutputFormatDouble

  write(OutputFormatDouble,'(1es13.5)') real_value
  
  OutputFormatDouble = adjustl(OutputFormatDouble)
  
end function OutputFormatDouble

! ************************************************************************** !
!
! GetCellCoordinates: Extracts coordinates of cells into a PetscVec
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GetCellCoordinates(grid,vec,direction)

  use Grid_module
  use Variables_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  PetscErrorCode :: ierr
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  
  if (direction == X_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%x(grid%nL2G(i))
    enddo
  else if (direction == Y_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%y(grid%nL2G(i))
    enddo
  else if (direction == Z_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%z(grid%nL2G(i))
    enddo
  endif
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine GetCellCoordinates

! ************************************************************************** !
!
! GetVertexCoordinates: Extracts vertex coordinates of cells into a PetscVec
! author: Gautam Bisht
! date: 11/01/2011
!
! ************************************************************************** !
subroutine GetVertexCoordinates(grid,vec,direction,option)

  use Grid_module
  use Option_module
  use Variables_module, only : X_COORDINATE, Y_COORDINATE, Z_COORDINATE
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  type(option_type) :: option
  
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: indices(:)
  PetscReal, allocatable :: values(:)
  PetscErrorCode :: ierr
  
  if (option%mycommsize == 1) then
    call VecGetArrayF90(vec,vec_ptr,ierr)
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%z
        enddo
    end select
    call VecRestoreArrayF90(vec,vec_ptr,ierr)
  else
    ! initialize to -999 to catch bugs
    call VecSet(vec,-999.d0,ierr)
    allocate(values(grid%unstructured_grid%num_vertices_local))
    allocate(indices(grid%unstructured_grid%num_vertices_local))
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%z
        enddo
    end select
    indices(:) = grid%unstructured_grid%vertex_ids_natural(:)-1
    call VecSetValues(vec,grid%unstructured_grid%num_vertices_local, &
                      indices,values,INSERT_VALUES,ierr)
    call VecAssemblyBegin(vec,ierr)
    deallocate(values)
    deallocate(indices)
    call VecAssemblyEnd(vec,ierr)
  endif
  
end subroutine GetVertexCoordinates

! ************************************************************************** !
!
! ExplicitGetCellCoordinates: Extracts cell coordinates for explicit grid
! into a PetscVec
! author: Satish Karra, LANL
! date: 12/11/12
!
! ************************************************************************** !
subroutine ExplicitGetCellCoordinates(grid,vec,direction,option)

  use Grid_module
  use Option_module
  use Variables_module, only : X_COORDINATE, Y_COORDINATE, Z_COORDINATE
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  type(option_type) :: option
  
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: indices(:)
  PetscReal, allocatable :: values(:)
  PetscErrorCode :: ierr  
  
  if (option%mycommsize == 1) then
    call VecGetArrayF90(vec,vec_ptr,ierr)
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%ngmax
          vec_ptr(ivertex) = grid%unstructured_grid%explicit_grid%vertex_coordinates(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%ngmax
          vec_ptr(ivertex) = grid%unstructured_grid%explicit_grid%vertex_coordinates(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%ngmax
          vec_ptr(ivertex) = grid%unstructured_grid%explicit_grid%vertex_coordinates(ivertex)%z
        enddo
    end select
    call VecRestoreArrayF90(vec,vec_ptr,ierr)
  else
    ! initialize to -999 to catch bugs
    call VecSet(vec,-999.d0,ierr)
    allocate(values(grid%nlmax))
    allocate(indices(grid%nlmax))
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%nlmax
          values(ivertex) = grid%unstructured_grid%explicit_grid%vertex_coordinates(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%nlmax
          values(ivertex) = grid%unstructured_grid%explicit_grid%vertex_coordinates(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%nlmax
          values(ivertex) = grid%unstructured_grid%explicit_grid%vertex_coordinates(ivertex)%z
        enddo
    end select
    indices(:) = grid%unstructured_grid%cell_ids_natural(:)-1
    call VecSetValues(vec,grid%nlmax,indices,values,INSERT_VALUES,ierr)
    call VecAssemblyBegin(vec,ierr)
    deallocate(values)
    deallocate(indices)
    call VecAssemblyEnd(vec,ierr)
  endif
  
  
end subroutine ExplicitGetCellCoordinates

! ************************************************************************** !
!> This routine returns a vector containing vertex ids in natural order of
!! local cells for unstructured grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/31/12
! ************************************************************************** !
subroutine GetCellConnections(grid, vec)

  use Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(grid_type) :: grid
  type(unstructured_grid_type),pointer :: ugrid
  Vec :: vec
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  ugrid => grid%unstructured_grid
  
  call GridVecGetArrayF90(grid, vec, vec_ptr, ierr)

  ! initialize
  vec_ptr = -999.d0
  do local_id=1, ugrid%nlmax
    ghosted_id = local_id
    select case(ugrid%cell_type(ghosted_id))
      case(HEX_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 8
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case(WEDGE_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 6
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        vec_ptr(offset + 7) = 0
        vec_ptr(offset + 8) = 0
      case (PYR_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 5
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        do ivertex = 6, 8
          vec_ptr(offset + ivertex) = 0
        enddo
      case (TET_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        do ivertex = 5, 8
          vec_ptr(offset + ivertex) = 0
        enddo
      case (QUAD_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case (TRI_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        ivertex = 4
        vec_ptr(offset + 4) = 0
    end select
  enddo

  call GridVecRestoreArrayF90(grid, vec, vec_ptr, ierr)

end subroutine GetCellConnections

! ************************************************************************** !
!> This subroutine writes header to a .xmf file
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/29/12
! ************************************************************************** !
subroutine OutputXMFHeader(fid,time,nmax,xmf_vert_len,ngvert,filename)

  implicit none

  PetscInt :: fid, vert_count
  PetscReal :: time
  PetscInt :: nmax,xmf_vert_len,ngvert
  character(len=MAXSTRINGLENGTH) :: filename

  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: comma_count, quote_count, variable_count
  PetscInt :: i
  
  string="<?xml version=""1.0"" ?>"
  write(fid,'(a)') trim(string)
  
  string="<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
  write(fid,'(a)') trim(string)

  string="<Xdmf>"
  write(fid,'(a)') trim(string)

  string="  <Domain>"
  write(fid,'(a)') trim(string)

  string="    <Grid Name=""Mesh"">"
  write(fid,'(a)') trim(string)

  write(string2,'(es13.5)') time
  string="      <Time Value = """ // trim(adjustl(string2)) // """ />"
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string="      <Topology Type=""Mixed"" NumberOfElements=""" // &
    trim(adjustl(string2)) // """ >"
  write(fid,'(a)') trim(string)

  write(string2,*) xmf_vert_len
  string="        <DataItem Format=""HDF"" DataType=""Int"" Dimensions=""" // &
    trim(adjustl(string2)) // """>"
  write(fid,'(a)') trim(string)

  string="          "//trim(filename) //":/Domain/Cells"
  write(fid,'(a)') trim(string)

  string="        </DataItem>"
  write(fid,'(a)') trim(string)

  string="      </Topology>"
  write(fid,'(a)') trim(string)

  string="      <Geometry GeometryType=""XYZ"">"
  write(fid,'(a)') trim(string)

  write(string2,*) ngvert
  string="        <DataItem Format=""HDF"" Dimensions=""" // trim(adjustl(string2)) // " 3"">"
  write(fid,'(a)') trim(string)

  string="          "//trim(filename) //":/Domain/Vertices"
  write(fid,'(a)') trim(string)

  string="        </DataItem>"
  write(fid,'(a)') trim(string)

  string="      </Geometry>"
  write(fid,'(a)') trim(string)

  string="      <Attribute Name=""XC"" AttributeType=""Scalar""  Center=""Cell"">"
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string="        <DataItem Dimensions=""" // trim(adjustl(string2)) // " 1"" Format=""HDF""> "
  write(fid,'(a)') trim(string)

  string="          "//trim(filename) //":/Domain/XC"
  write(fid,'(a)') trim(string)

  string="        </DataItem> " 
  write(fid,'(a)') trim(string)

  string="      </Attribute>"
  write(fid,'(a)') trim(string)

  string="      <Attribute Name=""YC"" AttributeType=""Scalar""  Center=""Cell"">"
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string="        <DataItem Dimensions=""" // trim(adjustl(string2)) // " 1"" Format=""HDF""> "
  write(fid,'(a)') trim(string)

  string="          "//trim(filename) //":/Domain/YC"
  write(fid,'(a)') trim(string)

  string="        </DataItem> " 
  write(fid,'(a)') trim(string)

  string="      </Attribute>"
  write(fid,'(a)') trim(string)

  string="      <Attribute Name=""ZC"" AttributeType=""Scalar""  Center=""Cell"">"
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string="        <DataItem Dimensions=""" // trim(adjustl(string2)) // " 1"" Format=""HDF""> "
  write(fid,'(a)') trim(string)

  string="          "//trim(filename) //":/Domain/ZC"
  write(fid,'(a)') trim(string)

  string="        </DataItem> " 
  write(fid,'(a)') trim(string)

  string="      </Attribute>"
  write(fid,'(a)') trim(string)

end subroutine OutputXMFHeader

! ************************************************************************** !
!> This subroutine writes footer to a .xmf file
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/29/12
! ************************************************************************** !
subroutine OutputXMFFooter(fid)

  implicit none

  PetscInt :: fid

  character(len=MAXSTRINGLENGTH) :: string

  string="    </Grid>"
  write(fid,'(a)') trim(string)

  string="  </Domain>"
  write(fid,'(a)') trim(string)

  string="</Xdmf>"
  write(fid,'(a)') trim(string)

end subroutine OutputXMFFooter

! ************************************************************************** !
!> This subroutine writes an attribute to a .xmf file
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/29/12
! ************************************************************************** !
subroutine OutputXMFAttribute(fid,nmax,attname,att_datasetname)

  implicit none

  PetscInt :: fid,nmax
  
  character(len=MAXSTRINGLENGTH) :: attname, att_datasetname
  character(len=MAXSTRINGLENGTH) :: string,string2
  string="      <Attribute Name=""" // trim(attname) // &
    """ AttributeType=""Scalar""  Center=""Cell"">"
  write(fid,'(a)') trim(string)

!  write(string2,*) grid%nmax
  write(string2,*) nmax
  string="        <DataItem Dimensions=""" // trim(adjustl(string2)) // " 1"" Format=""HDF""> "
  write(fid,'(a)') trim(string)

  string="        " // trim(att_datasetname)
  write(fid,'(a)') trim(string)

  string="        </DataItem> " 
  write(fid,'(a)') trim(string)

  string="      </Attribute>"
  write(fid,'(a)') trim(string)

end subroutine OutputXMFAttribute

! ************************************************************************** !
!> This returns mass/energy flowrate at all faces of a control volume
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/21/2013
! ************************************************************************** !
subroutine OutputGetFlowrates(realization_base)

  use HDF5_module
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Variables_module
  use Connection_module
  use Coupler_module
  use HDF5_Aux_module
  use Output_Aux_module
  use Field_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"
#include "definitions.h"

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(unstructured_grid_type),pointer :: ugrid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition
  type(ugdm_type),pointer :: ugdm
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: idual
  PetscInt :: iconn
  PetscInt :: face_id
  PetscInt :: local_id_up,local_id_dn
  PetscInt :: ghosted_id_up,ghosted_id_dn
  PetscInt :: iface_up,iface_dn
  PetscInt :: dof
  PetscInt :: sum_connection
  PetscInt :: offset
  PetscInt :: cell_type
  PetscInt :: local_size
  PetscInt :: i
  PetscInt :: iface
  PetscInt :: ndof

  PetscReal, pointer :: flowrates(:,:,:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal, pointer :: vec_ptr3(:)
  PetscReal, pointer :: double_array(:)
  PetscReal :: dtime

  Vec :: natural_flowrates_vec

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: unit_string

  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  output_option =>realization_base%output_option
  option => realization_base%option
  field => realization_base%field

  ! Create UGDM for
  call UGridCreateUGDM(grid%unstructured_grid,ugdm, &
                       (option%nflowdof*MAX_FACE_PER_CELL + 1),option)

  ! Create a flowrate vector in natural order
  call UGridDMCreateVector(grid%unstructured_grid,ugdm,natural_flowrates_vec, &
                           NATURAL,option)


  allocate(flowrates(option%nflowdof,MAX_FACE_PER_CELL,ugrid%nlmax))
  flowrates = 0.d0
  call VecGetArrayF90(field%flowrate_inst,vec_ptr,ierr)
  vec_ptr = 0.d0
  
  offset = 1+option%nflowdof*MAX_FACE_PER_CELL
  ! Save the number of faces of all cell
  do local_id = 1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    cell_type = ugrid%cell_type(ghosted_id)
    vec_ptr((local_id-1)*offset+1) = UCellGetNFaces(cell_type,option)
  enddo

  ! Interior Flowrates Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      face_id = cur_connection_set%face_id(iconn)
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      do iface_up = 1,MAX_FACE_PER_CELL
        if(face_id==ugrid%cell_to_face_ghosted(iface_up,local_id_up)) exit
      enddo
      iface_dn=-1
      if(local_id_dn>0) then
        do iface_dn = 1,MAX_FACE_PER_CELL
          if(face_id==ugrid%cell_to_face_ghosted(iface_dn,local_id_dn)) exit
        enddo
      endif
      
      do dof=1,option%nflowdof
        ! Save flowrate for iface_up of local_id_up cell using flowrate up-->dn
        flowrates(dof,iface_up,local_id_up) = patch%internal_fluxes(dof,1,sum_connection)
        vec_ptr((local_id_up-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_up + 1) = &
          patch%internal_fluxes(dof,1,sum_connection)

        if(iface_dn>0) then
          ! Save flowrate for iface_dn of local_id_dn cell using -ve flowrate up-->dn
          flowrates(dof,iface_dn,local_id_dn) = -patch%internal_fluxes(dof,1,sum_connection)
          vec_ptr((local_id_dn-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_dn + 1) = &
            -patch%internal_fluxes(dof,1,sum_connection)
        endif
      enddo
      
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flowrates Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0
  do 
    if (.not.associated(boundary_condition)) exit

    cur_connection_set => boundary_condition%connection_set
    sum_connection = 0

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      face_id = cur_connection_set%face_id(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      do iface_dn = 1,MAX_FACE_PER_CELL
        if(face_id==ugrid%cell_to_face_ghosted(iface_dn,local_id_dn)) exit
      enddo

      do dof=1,option%nflowdof
        ! Save flowrate for iface_dn of local_id_dn cell using -ve flowrate up-->dn
        flowrates(dof,iface_dn,local_id_dn) = -patch%boundary_fluxes(dof,1,sum_connection)
        vec_ptr((local_id_dn-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_dn + 1) = &
          -patch%boundary_fluxes(dof,1,sum_connection)
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

  deallocate(flowrates)

  call VecRestoreArrayF90(field%flowrate_inst,vec_ptr,ierr)

  ! Scatter flowrate from Global --> Natural order
  call VecScatterBegin(ugdm%scatter_gton,field%flowrate_inst,natural_flowrates_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm%scatter_gton,field%flowrate_inst,natural_flowrates_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(natural_flowrates_vec,vec_ptr,ierr)
  call VecGetArrayF90(field%flowrate_inst,vec_ptr2,ierr)

  ! Copy the vectors
  vec_ptr2 = vec_ptr
  
  if(output_option%print_hdf5_aveg_mass_flowrate.or. &
    output_option%print_hdf5_aveg_energy_flowrate) then

    dtime = option%time-output_option%aveg_var_time
    call VecGetArrayF90(field%flowrate_aveg,vec_ptr3,ierr)
    vec_ptr3 = vec_ptr3 + vec_ptr2/dtime
    call VecRestoreArrayF90(field%flowrate_aveg,vec_ptr3,ierr)
      
  endif

  call VecRestoreArrayF90(natural_flowrates_vec,vec_ptr,ierr)
  call VecRestoreArrayF90(field%flowrate_inst,vec_ptr2,ierr)

  call VecDestroy(natural_flowrates_vec,ierr)
  call UGridDMDestroy(ugdm)
  
end subroutine OutputGetFlowrates

! ************************************************************************** !
!
! OutputGetExplicitFlowrates: Forms a vector of magnitude of flowrates
! which will be printed out to file for particle tracking.
! author: Satish Karra, LANL
! date: 04/24/13
!
! ************************************************************************** !
subroutine OutputGetExplicitFlowrates(realization_base)

  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Unstructured_Grid_Aux_module
  use Field_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"
#include "definitions.h"

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(unstructured_grid_type),pointer :: ugrid
  type(field_type), pointer :: field
  
  PetscInt :: dof
  PetscInt :: offset
  PetscInt :: istart, iend
  PetscInt :: iconn
  PetscErrorCode :: ierr
  PetscReal :: val
  
  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  option => realization_base%option
  field => realization_base%field


  call VecGetOwnershipRange(field%flowrate_inst,istart,iend,ierr)
  
  offset = option%nflowdof

  do iconn = istart,iend-1
    do dof = 1,option%nflowdof
      val = abs(patch%internal_fluxes(dof,1,iconn))
      call VecSetValues(field%flowrate_inst,ONE_INTEGER,(iconn-1)*offset + dof, &
                        val,INSERT_VALUES,ierr) 
    enddo
  enddo
   
 call VecAssemblyBegin(field%flowrate_inst,ierr)
 call VecAssemblyEnd(field%flowrate_inst,ierr)

end subroutine OutputGetExplicitFlowrates


end module Output_Common_module
