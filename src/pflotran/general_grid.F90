module General_Grid_module

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
  
  integer :: hdf5_err
  PetscErrorCode :: ierr
  public :: ReadStructuredGridHDF5
  
contains
  
#ifndef USE_HDF5
subroutine ReadStructuredGridHDF5(realization)

  use Realization_module
  use Option_module
  
  implicit none

  type(realization_type) :: realization

  if (realization%option%myrank == 0) then
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
subroutine ReadStructuredGridHDF5(realization)

  use hdf5
  
  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Connection_module
  use HDF5_module
  
  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field

  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXSTRINGLENGTH) :: filename
  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id

  integer :: i, local_ghosted_id, iconn
  integer, allocatable :: indices(:)
  integer, allocatable :: integer_array(:)
  real*8, allocatable :: real_array(:)
  
  Vec :: global
  Vec :: local
  PetscScalar, pointer :: vec_ptr(:)

  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
      
  PetscTruth :: option_found

  PetscLogDouble :: time0, time1, time3, time4

  field => realization%field
  option => realization%option
  grid => realization%grid

  call PetscGetTime(time0, ierr)

  filename = option%generalized_grid
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-hdf5_grid', &
                             filename, option_found, ierr)
 
  ! grab connection object
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  
  ! create hash table for fast lookup
#ifdef HASH
  call GridCreateNaturalToGhostedHash(grid,option)
!  call GridPrintHashTable
#endif

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)

  if (option%myrank == 0) print *, 'Opening hdf5 file: ', trim(filename)
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  ! open the grid cells group
  string = 'Grid Cells'
  if (option%myrank == 0) print *, 'Opening group: ', trim(string)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  
  allocate(indices(grid%nlmax))
  
  if (option%myrank == 0) print *, 'Setting up grid cell indices'
  call PetscGetTime(time3, ierr)
  call HDF5MapLocalToNaturalIndices(grid,option,grp_id,indices)
  call PetscGetTime(time4, ierr)
  time4 = time4 - time3
  if (option%myrank == 0) print *, time4, &
       ' seconds to set up grid cell indices for hdf5 file'

  call GridCreateVector(grid,ONEDOF,global,GLOBAL)
  call GridCreateVector(grid,ONEDOF,local,LOCAL)
  
  allocate(integer_array(grid%nlmax))
  string = "Material Id"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadIntegerArray(grid,option,grp_id,grid%nlmax,indices,string,integer_array)
  call GridCopyIntegerArrayToPetscVec(integer_array,global,grid%nlmax)
  deallocate(integer_array)
  call GridGlobalToLocal(grid,global,local,ONEDOF)
  call GridCopyPetscVecToIntegerArray(field%imat,local,grid%ngmax)
  
  allocate(real_array(grid%nlmax))
  string = "X-Coordinate"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,real_array)
  call GridCopyRealArrayToPetscVec(real_array,global,grid%nlmax)
  call GridGlobalToLocal(grid,global,local,ONEDOF)  
  call GridCopyPetscVecToRealArray(grid%x,local,grid%ngmax)

  string = "Y-Coordinate"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,real_array)
  call GridCopyRealArrayToPetscVec(real_array,global,grid%nlmax)
  call GridGlobalToLocal(grid,global,local,ONEDOF)  
  call GridCopyPetscVecToRealArray(grid%y,local,grid%ngmax)

  string = "Z-Coordinate"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,real_array)
  call GridCopyRealArrayToPetscVec(real_array,global,grid%nlmax)
  call GridGlobalToLocal(grid,global,local,ONEDOF)  
  call GridCopyPetscVecToRealArray(grid%z,local,grid%ngmax)

  deallocate(real_array)
  
  string = "Volume"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call VecGetArrayF90(grid%volume,vec_ptr,ierr); CHKERRQ(ierr)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(grid%volume,vec_ptr,ierr); CHKERRQ(ierr)
  
  string = "X-Permeability"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call VecGetArrayF90(field%perm0_xx,vec_ptr,ierr); CHKERRQ(ierr)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(field%perm0_xx,vec_ptr,ierr); CHKERRQ(ierr)
  
  string = "Y-Permeability"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call VecGetArrayF90(field%perm0_yy,vec_ptr,ierr); CHKERRQ(ierr)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(field%perm0_yy,vec_ptr,ierr); CHKERRQ(ierr)

  string = "Z-Permeability"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call VecGetArrayF90(field%perm0_zz,vec_ptr,ierr); CHKERRQ(ierr)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(field%perm0_zz,vec_ptr,ierr); CHKERRQ(ierr)

  string = "Porosity"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call VecGetArrayF90(field%porosity0,vec_ptr,ierr); CHKERRQ(ierr)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(field%porosity0,vec_ptr,ierr); CHKERRQ(ierr)

  string = "Tortuosity"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call VecGetArrayF90(global,vec_ptr,ierr); CHKERRQ(ierr)
  call HDF5ReadRealArray(grid,option,grp_id,grid%nlmax,indices,string,vec_ptr)
  call VecRestoreArrayF90(global,vec_ptr,ierr); CHKERRQ(ierr)
  call GridGlobalToLocal(grid,global,field%tor_loc,ONEDOF)

  call VecDestroy(global,ierr)
  call VecDestroy(local,ierr)

  ! update local vectors
  call UpdateGlobalToLocal(grid,field)

  deallocate(indices)
  if (option%myrank == 0) print *, 'Closing group: Grid Cells'
  call h5gclose_f(grp_id,hdf5_err)
  
  string = 'Connections'
  if (option%myrank == 0) print *, 'Opening group: ', trim(string)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

  allocate(indices(cur_connection_set%num_connections))

  if (option%myrank == 0) print *, 'Setting up connection indices'
  call PetscGetTime(time3, ierr)
  call SetupConnectionIndices(grid,option,grp_id,indices)
  call PetscGetTime(time4, ierr)
  time4 = time4 - time3
  if (option%myrank == 0) print *, time4, &
       ' seconds to set up connection indices for hdf5 file'

  allocate(integer_array(cur_connection_set%num_connections))
  string = "Id Upwind"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadIntegerArray(grid,option,grp_id,cur_connection_set%num_connections, &
                        indices,string,integer_array)
  do i=1,cur_connection_set%num_connections
    local_ghosted_id = GridGetLocalGhostedIdFromHash(grid,integer_array(i))
    cur_connection_set%id_up(i) = local_ghosted_id
  enddo

  string = "Id Downwind"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadIntegerArray(grid,option,grp_id,cur_connection_set%num_connections, &
                        indices,string,integer_array)
  do i=1,cur_connection_set%num_connections
    local_ghosted_id = GridGetLocalGhostedIdFromHash(grid,integer_array(i))
    cur_connection_set%id_dn(i) = local_ghosted_id
  enddo
  deallocate(integer_array)
  
  allocate(real_array(cur_connection_set%num_connections))
  string = "Distance Upwind"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadRealArray(grid,option,grp_id,cur_connection_set%num_connections, &
                     indices,string,real_array)
                     
  cur_connection_set%dist(-1,1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      
  cur_connection_set%dist(0,1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      

  string = "Distance Downwind"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadRealArray(grid,option,grp_id,cur_connection_set%num_connections, &
                     indices,string,real_array) 

  cur_connection_set%dist(0,1:cur_connection_set%num_connections) = &
    cur_connection_set%dist(0,1:cur_connection_set%num_connections) + &
    real_array(1:cur_connection_set%num_connections)                      
  ! compute upwind fraction
  cur_connection_set%dist(-1,1:cur_connection_set%num_connections) = &
    cur_connection_set%dist(-1,1:cur_connection_set%num_connections) / &
    cur_connection_set%dist(0,1:cur_connection_set%num_connections)

  string = "Area"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadRealArray(grid,option,grp_id,cur_connection_set%num_connections, &
                     indices,string,real_array) 

  cur_connection_set%area(1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      

  string = "CosB"
  if (option%myrank == 0) print *, 'Reading dataset: ', trim(string)
  call HDF5ReadRealArray(grid,option,grp_id,cur_connection_set%num_connections, &
                     indices,string,real_array) 
  
  ! for now, store cosB in z component of distance array
  cur_connection_set%dist(3,1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      
  
  deallocate(real_array)
  deallocate(indices)

  if (option%myrank == 0) print *, 'Closing group: Connections'
  call h5gclose_f(grp_id,hdf5_err)
  if (option%myrank == 0) print *, 'Closing hdf5 file: ', filename
  call h5fclose_f(file_id,hdf5_err)
   
  call h5close_f(hdf5_err)
  
  ! surely a kludge here!
!  do iconn = 1, cur_connection_set%num_connections
!    if (cur_connection_set%dist(3,iconn) > 0.99d0) then
!      cur_connection_set%dist(3,iconn) = &
!        grid%z(cur_connection_set%id_dn(iconn))- &
!        grid%z(cur_connection_set%id_up(iconn))
!    endif
!  enddo
  
  ! set up delz array
!  do i=1,num_internal_connections
!    grid%delz(i) = grid%z(grid%nd2(i))-grid%z(grid%nd1(i))
!  enddo

  call GridDestroyHashTable(grid)

  call PetscGetTime(time1, ierr)
  time1 = time1 - time0
  if (option%myrank == 0) print *, time1, &
       ' seconds to read unstructured grid data from hdf5 file'
     
end subroutine ReadStructuredGridHDF5

! ************************************************************************** !
!
! SetupConnectionIndices: Set up indices array that map local connection to  
!                         entries in HDF5 grid connection vectors
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine SetupConnectionIndices(grid,option,file_id,indices)

  use hdf5
  use Connection_module
  use Grid_module
  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
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
  integer :: temp_int, i, num_internal_connections
  
  integer, allocatable :: upwind_ids(:), downwind_ids(:)
  
  integer :: read_block_size = HDF5_READ_BUFFER_SIZE

  num_internal_connections = ConnectionGetNumberInList(grid%internal_connection_list)
  
  string = "Id Upwind"
  call h5dopen_f(file_id,string,data_set_id_up,hdf5_err)
  call h5dget_space_f(data_set_id_up,file_space_id_up,hdf5_err)
  string = "Id Downwind"
  call h5dopen_f(file_id,string,data_set_id_down,hdf5_err)
  call h5dget_space_f(data_set_id_down,file_space_id_down,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id_up, &
                                      num_connections_in_file,hdf5_err)
#if 1
  if (num_connections_in_file /= num_internal_connections) then
    if (option%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_connections_in_file, ') does not match the dimensions ', &
               'of the domain (', num_internal_connections, ').'
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
    temp_int = num_connections_in_file-connection_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
    endif
    ! offset is zero-based
    offset(1) = connection_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id_up, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (option%myrank == 0) then                           
#endif
      call h5dread_f(data_set_id_up,H5T_NATIVE_INTEGER,upwind_ids,dims, &
                     hdf5_err,memory_space_id,file_space_id_up,prop_id)                     
#ifdef HDF5_BROADCAST
    endif
    if (option%commsize > 1) &
      call mpi_bcast(upwind_ids,dims(1),MPI_INTEGER,0, &
                     PETSC_COMM_WORLD,ierr)
#endif    
    call h5sselect_hyperslab_f(file_space_id_down, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (option%myrank == 0) then                           
#endif
      call h5dread_f(data_set_id_down,H5T_NATIVE_INTEGER,downwind_ids,dims, &
                     hdf5_err,memory_space_id,file_space_id_down,prop_id)                     
#ifdef HDF5_BROADCAST
    endif
    if (option%commsize > 1) &
      call mpi_bcast(downwind_ids,dims(1),MPI_INTEGER,0, &
                     PETSC_COMM_WORLD,ierr)
#endif    
    do i=1,dims(1)
      connection_count = connection_count + 1
      natural_id_up = upwind_ids(i)
      natural_id_down = downwind_ids(i)
      local_ghosted_id_up = GridGetLocalGhostedIdFromHash(grid,natural_id_up)
      local_ghosted_id_down = GridGetLocalGhostedIdFromHash(grid,natural_id_down)
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

  if (index_count /= num_internal_connections) then
    if (option%myrank == 0) &
      print *, 'ERROR: Number of indices read (', index_count, ') does not ', &
               'match the number of local grid connections (', num_internal_connections, ').'
      call PetscFinalize(ierr)
      stop
  endif

end subroutine SetupConnectionIndices
#endif

! ************************************************************************** !
!
! UpdateGlobalToLocal: Updated global vec values to local
! author: Glenn Hammond
! date: 06/20/07
!
! ************************************************************************** !
subroutine UpdateGlobalToLocal(grid,field)

  use Grid_module
  use Field_module
  
  implicit none
  
  integer :: ierr
  type(grid_type) :: grid
  type(field_type) :: field

  ! icap
  call GridLocalToLocal(grid,field%icap_loc,field%icap_loc,ONEDOF)

  ! ithrm
  call GridLocalToLocal(grid,field%ithrm_loc,field%ithrm_loc,ONEDOF)

  ! perm_xx, perm_yy, perm_zz
  call GridGlobalToLocal(grid,field%perm0_xx,field%perm_xx_loc,ONEDOF)
  call GridGlobalToLocal(grid,field%perm0_yy,field%perm_yy_loc,ONEDOF)
  call GridGlobalToLocal(grid,field%perm0_zz,field%perm_zz_loc,ONEDOF)

  ! tor
  call GridLocalToLocal(grid,field%tor_loc,field%tor_loc,ONEDOF)

  ! por
  call GridGlobalToLocal(grid,field%porosity0,field%porosity_loc,ONEDOF)

end subroutine UpdateGlobalToLocal


end module General_Grid_module
