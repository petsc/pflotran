module General_Grid_module

#include "finclude/petscvec.h"
  use petscvec
  use PFLOTRAN_Constants_module

 implicit none

#define INVERT

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE  
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

 private 

  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.

  PetscMPIInt :: hdf5_err
  public :: ReadStructuredGridHDF5
  
contains
  
#if !defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine ReadStructuredGridHDF5(realization)

  use Realization_class
  use Option_module
  
  implicit none

  type(realization_type) :: realization

  call printMsg(realization%option,'')
  write(realization%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 &
          &to read HDF5 formatted structured grids.",/)')
  call printErrMsg(realization%option)
  
end subroutine ReadStructuredGridHDF5

#else

! ************************************************************************** !

subroutine ReadStructuredGridHDF5(realization)
  ! 
  ! ReadHDF5StructuredGrid: Reads in a structured grid in HDF5 format
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/21/07
  ! 

  use hdf5
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Connection_module
  use HDF5_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXSTRINGLENGTH) :: filename
  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id

  PetscInt :: i, local_ghosted_id, iconn
  PetscInt, pointer :: indices(:)
  PetscInt, allocatable :: integer_array(:)
  PetscReal, pointer :: real_array(:)
  
  Vec :: global_vec
  Vec :: local_vec
  PetscReal, pointer :: vec_ptr(:)

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
      
  PetscBool :: option_found
  PetscErrorCode :: ierr

  PetscLogDouble :: time0, time1, time3, time4

  field => realization%field
  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  call PetscTime(time0, ierr);CHKERRQ(ierr)

  filename = option%generalized_grid
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-hdf5_grid', &
                             filename, option_found, ierr);CHKERRQ(ierr)
 
  ! grab connection object
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  
  ! create hash table for fast lookup
  call GridCreateNaturalToGhostedHash(grid,option)
!  call GridPrintHashTable

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)

  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option)
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  ! open the grid cells group
  string = 'Grid Cells'
  option%io_buffer = 'Opening group: Grid Cells'
  call printMsg(option)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  
  allocate(indices(grid%nlmax))
  
  option%io_buffer = 'Setting up grid cell indices'
  call printMsg(option)
  call PetscTime(time3, ierr);CHKERRQ(ierr)
  string = 'Cell Id'
  call HDF5MapLocalToNaturalIndices(grid,option,grp_id,string,grid%nmax, &
                                    indices,grid%nlmax)
  call PetscTime(time4, ierr);CHKERRQ(ierr)
  time4 = time4 - time3
  write(option%io_buffer, &
        '(f6.2," seconds to set up grid cell indices for hdf5 file")') time4
  call printMsg(option)

  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,local_vec,LOCAL, &
                                  option)
  
  allocate(integer_array(grid%nlmax))
  string = "Material Id"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call HDF5ReadIntegerArray(option,grp_id,string,grid%nmax,indices, &
                            grid%nlmax,integer_array)
  call GridCopyIntegerArrayToVec(grid,integer_array,global_vec,grid%nlmax)
  deallocate(integer_array)
  call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)

  call GridCopyVecToIntegerArray(grid,patch%imat,local_vec,grid%ngmax)

  
  allocate(real_array(grid%nlmax))
  string = "X-Coordinate"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         real_array)

  call GridCopyRealArrayToVec(grid,real_array,global_vec,grid%nlmax)

  call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)  

  call GridCopyVecToRealArray(grid,grid%x,local_vec,grid%ngmax)

  string = "Y-Coordinate"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         real_array)

  call GridCopyRealArrayToVec(grid,real_array,global_vec,grid%nlmax)

  call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)  

  call GridCopyVecToRealArray(grid,grid%y,local_vec,grid%ngmax)

  string = "Z-Coordinate"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         real_array)

  call GridCopyRealArrayToVec(grid,real_array,global_vec,grid%nlmax)
  call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)  

  call GridCopyVecToRealArray(grid,grid%z,local_vec,grid%ngmax)

  deallocate(real_array)
  
  string = "Volume"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call VecGetArrayF90(field%volume,vec_ptr,ierr);CHKERRQ(ierr)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         vec_ptr)
  call VecRestoreArrayF90(field%volume,vec_ptr,ierr);CHKERRQ(ierr)
  
  string = "X-Permeability"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call VecGetArrayF90(field%perm0_xx,vec_ptr,ierr);CHKERRQ(ierr)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         vec_ptr)
  call VecRestoreArrayF90(field%perm0_xx,vec_ptr,ierr);CHKERRQ(ierr)
  
  string = "Y-Permeability"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call VecGetArrayF90(field%perm0_yy,vec_ptr,ierr);CHKERRQ(ierr)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         vec_ptr)
  call VecRestoreArrayF90(field%perm0_yy,vec_ptr,ierr);CHKERRQ(ierr)

  string = "Z-Permeability"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call VecGetArrayF90(field%perm0_zz,vec_ptr,ierr);CHKERRQ(ierr)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         vec_ptr)
  call VecRestoreArrayF90(field%perm0_zz,vec_ptr,ierr);CHKERRQ(ierr)

  string = "Porosity"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call VecGetArrayF90(field%porosity0,vec_ptr,ierr);CHKERRQ(ierr)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         vec_ptr)
  call VecRestoreArrayF90(field%porosity0,vec_ptr,ierr);CHKERRQ(ierr)

  string = "Tortuosity"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)
  call VecGetArrayF90(field%tortuosity0,vec_ptr,ierr);CHKERRQ(ierr)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         vec_ptr)
  call VecRestoreArrayF90(field%tortuosity0,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(local_vec,ierr);CHKERRQ(ierr)

  ! update local vectors
  call UpdateGlobalToLocal(discretization,field)

  deallocate(indices)
  option%io_buffer = 'Closing group: Grid Cells'
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)
  
  option%io_buffer = 'Opening group: Connections' 
  call printMsg(option)   
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

  allocate(indices(cur_connection_set%num_connections))

  option%io_buffer = 'Setting up connection indices'
  call printMsg(option)   
  call PetscTime(time3, ierr);CHKERRQ(ierr)
  call SetupConnectionIndices(grid,option,grp_id,indices)
  call PetscTime(time4, ierr);CHKERRQ(ierr)
  time4 = time4 - time3
  write(option%io_buffer, &
        '(f6.2," seconds to set up connection indices for hdf5 file")') &
        time4
  call printMsg(option)     

  allocate(integer_array(cur_connection_set%num_connections))
  string = "Id Upwind"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIntegerArray(option,grp_id,string, &
                            cur_connection_set%num_connections, &
                            indices,cur_connection_set%num_connections, &
                            integer_array)
  do i=1,cur_connection_set%num_connections
    local_ghosted_id = GridGetLocalGhostedIdFromHash(grid,integer_array(i))
    cur_connection_set%id_up(i) = local_ghosted_id
  enddo

  string = "Id Downwind"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIntegerArray(option,grp_id,string, &
                            cur_connection_set%num_connections, &
                            indices,cur_connection_set%num_connections, &
                            integer_array)
  do i=1,cur_connection_set%num_connections
    local_ghosted_id = GridGetLocalGhostedIdFromHash(grid,integer_array(i))
    cur_connection_set%id_dn(i) = local_ghosted_id
  enddo
  deallocate(integer_array)
  
  allocate(real_array(cur_connection_set%num_connections))
  string = "Distance Upwind"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         real_array)
                     
  cur_connection_set%dist(-1,1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      
  cur_connection_set%dist(0,1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      

  string = "Distance Downwind"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         real_array) 

  cur_connection_set%dist(0,1:cur_connection_set%num_connections) = &
    cur_connection_set%dist(0,1:cur_connection_set%num_connections) + &
    real_array(1:cur_connection_set%num_connections)                      
  ! compute upwind fraction
  cur_connection_set%dist(-1,1:cur_connection_set%num_connections) = &
    cur_connection_set%dist(-1,1:cur_connection_set%num_connections) / &
    cur_connection_set%dist(0,1:cur_connection_set%num_connections)

  string = "Area"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         real_array) 

  cur_connection_set%area(1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      

  string = "CosB"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices,grid%nlmax, &
                         real_array) 
  
  ! for now, store cosB in z component of distance array
  cur_connection_set%dist(3,1:cur_connection_set%num_connections) = &
    real_array(1:cur_connection_set%num_connections)                      
  
  deallocate(real_array)
  deallocate(indices)

  option%io_buffer = 'Closing group: Connections'
  call printMsg(option)
  call h5gclose_f(grp_id,hdf5_err)
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)  
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

  call PetscTime(time1, ierr);CHKERRQ(ierr)
  time1 = time1 - time0
  write(option%io_buffer, &
        '(f6.2," seconds to read unstructured grid data from hdf5 file")') &
        time1
  call printMsg(option) 
       
end subroutine ReadStructuredGridHDF5

! ************************************************************************** !

subroutine SetupConnectionIndices(grid,option,file_id,indices)
  ! 
  ! Set up indices array that map local connection to
  ! entries in HDF5 grid connection vectors
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/21/07
  ! 

  use hdf5
  use Connection_module
  use Grid_module
  use Option_module
  use Logging_module
  use HDF5_Aux_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  PetscInt :: indices(:)
  
  character(len=MAXSTRINGLENGTH) :: string 
  integer(HID_T) :: file_space_id_up, file_space_id_down
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id_up, data_set_id_down
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank
  PetscInt :: local_ghosted_id_up, local_id_up, natural_id_up
  PetscInt :: local_ghosted_id_down, local_id_down, natural_id_down
  PetscInt :: index_count
  PetscInt :: connection_count
  integer(HSIZE_T) :: num_connections_in_file
  PetscInt :: temp_int, i, num_internal_connections
  PetscErrorCode :: ierr
  PetscMPIInt :: int_mpi
  
!  PetscInt, allocatable :: upwind_ids(:), downwind_ids(:)
  integer, allocatable :: upwind_ids_i4(:), downwind_ids_i4(:)
  
  PetscInt :: read_block_size

  read_block_size = HDF5_READ_BUFFER_SIZE
  num_internal_connections = ConnectionGetNumberInList(grid%internal_connection_set_list)
  
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
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimension",&
           &" of the domain (",i9,").")') trim(string), &
           num_connections_in_file, num_internal_connections
    call printErrMsg(option)
  endif
#endif
  
  allocate(upwind_ids_i4(read_block_size),downwind_ids_i4(read_block_size))
  
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
    temp_int = int(num_connections_in_file)-connection_count
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
    if (option%myrank == option%io_rank) then                           
#endif
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      call h5dread_f(data_set_id_up,HDF_NATIVE_INTEGER,upwind_ids_i4,dims, &
                     hdf5_err,memory_space_id,file_space_id_up,prop_id)                     
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
#ifdef HDF5_BROADCAST
    endif
    if (option%mycommsize > 1) then
      int_mpi = dims(1)
!geh      call MPI_Bcast(upwind_ids,int_mpi,MPIU_INTEGER,option%io_rank, &
      call MPI_Bcast(upwind_ids_i4,int_mpi,MPI_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    endif
#endif    
    call h5sselect_hyperslab_f(file_space_id_down, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (option%myrank == option%io_rank) then                           
#endif
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      call h5dread_f(data_set_id_down,HDF_NATIVE_INTEGER,downwind_ids_i4,dims, &
                     hdf5_err,memory_space_id,file_space_id_down,prop_id)                     
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
#ifdef HDF5_BROADCAST
    endif
    if (option%mycommsize > 1) then
      int_mpi = dims(1)
!geh      call MPI_Bcast(downwind_ids,int_mpi,MPIU_INTEGER,option%io_rank, &
      call MPI_Bcast(downwind_ids_i4,int_mpi,MPI_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    endif
#endif    
    do i=1,int(dims(1))
      connection_count = connection_count + 1
      natural_id_up = upwind_ids_i4(i)
      natural_id_down = downwind_ids_i4(i)
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
  
  deallocate(upwind_ids_i4,downwind_ids_i4)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id_up,hdf5_err)
  call h5sclose_f(file_space_id_down,hdf5_err)
  call h5dclose_f(data_set_id_up,hdf5_err)
  call h5dclose_f(data_set_id_down,hdf5_err)

  if (index_count /= num_internal_connections) then
    write(option%io_buffer, &
          '("Number of indices read (",i9,") does not match the number of",&
           &" local grid connections (",i9,").")') index_count, &
           num_internal_connections
    call printErrMsg(option)      
  endif

end subroutine SetupConnectionIndices
#endif

! ************************************************************************** !

subroutine UpdateGlobalToLocal(discretization,field)
  ! 
  ! Updated global vec values to local
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/20/07
  ! 

  use Discretization_module
  use Field_module
  
  implicit none
  
  type(discretization_type) :: discretization
  type(field_type) :: field

  ! icap
  call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                  field%icap_loc,ONEDOF)

  ! ithrm
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                  field%ithrm_loc,ONEDOF)

  ! perm_xx, perm_yy, perm_zz
  call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                   field%perm_xx_loc,ONEDOF)
  call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                   field%perm_yy_loc,ONEDOF)
  call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                   field%perm_zz_loc,ONEDOF)

  ! tor
  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
                                   field%tortuosity_loc,ONEDOF)

  ! por
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                                   field%porosity_loc,ONEDOF)

end subroutine UpdateGlobalToLocal




end module General_Grid_module
