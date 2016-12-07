module HDF5_Vamsi_module

#include "finclude/petscsys.h"
  use petscsys
  use Logging_module

  use PFLOTRAN_Constants_module

  implicit none

  private
  
  PetscErrorCode :: ierr

  PetscBool, public :: trick_hdf5 = PETSC_FALSE

#if defined(PETSC_HAVE_HDF5)
  PetscMPIInt :: hdf5_err
  PetscMPIInt :: io_rank
      
! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#endif
      
#if defined(PETSC_HAVE_HDF5)
  public :: HDF5MapLocToNatIndicesVamsi, &
            HDF5ReadIntegerArrayVamsi, &
            HDF5WriteStructuredDataSetVamsi, &
            HDF5ReadRegionFromFileVamsi, &       
            HDF5ReadCellIndxIntArrayVamsi, & 
            HDF5ReadCellIndxRealArrayVamsi
#else
  public :: HDF5ReadRegionFromFileVamsi, &
            HDF5ReadUnstructuredGridRegionFromFile, &
            HDF5ReadCellIndxIntArrayVamsi, &
            HDF5ReadCellIndxRealArrayVamsi
#endif            

contains

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine HDF5MapLocToNatIndicesVamsi(grid,option,file_id, &
                                        dataset_name,dataset_size, &
                                        indices,num_indices)
  ! 
  ! Set up indices array that maps local cells to
  ! entries in HDF5 grid cell vectors
  ! 
  ! Author: Vamsi Sripath
  ! Date: ??/??/09
  ! 

  use hdf5
  
  use Option_module
  use Grid_module
  
  implicit none

! Vamsi's HDF5 Mechanism 
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices

  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscInt :: local_ghosted_id, local_id, natural_id
  PetscInt :: index_count
  PetscInt :: cell_count
  integer(HSIZE_T) :: num_cells
  PetscInt :: num_cells_in_file
  PetscInt :: temp_int, i, counter
  PetscMPIInt :: int_mpi
  
  integer, allocatable :: cell_ids_i4(:)
  PetscInt, allocatable :: temp(:)
  
  PetscInt :: read_block_size
  PetscInt :: indices_array_size

  call PetscLogEventBegin(logging%event_map_indices_hdf5,ierr);CHKERRQ(ierr)

  dims = 0
  cell_count = 0
  index_count = 0
  memory_space_id = -1
  length = 0
  read_block_size = HDF5_READ_BUFFER_SIZE

  allocate(cell_ids_i4(read_block_size))
  if (num_indices > 0) then
    allocate(indices(num_indices))
  else
    indices_array_size = read_block_size
    allocate(indices(indices_array_size))
  endif
  
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
     call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
     call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
     ! should be a rank=1 data space
     call h5sget_simple_extent_npoints_f(file_space_id,num_cells,hdf5_err)
     if (dataset_size > 0 .and. num_cells_in_file /= dataset_size) then
        write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimension",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_cells_in_file, dataset_size
        call printErrMsg(option)    
     endif
     call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
     call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
     rank_mpi = 1
     offset = 0
     length = 0
     stride = 1
     num_cells_in_file = int(num_cells)
  endif
    
  call MPI_Bcast(num_cells_in_file,ONE_INTEGER_MPI,MPIU_INTEGER, &
                 ZERO_INTEGER_MPI,option%read_group,ierr)
                     
  do
    if (cell_count >= num_cells_in_file) exit
    temp_int = num_cells_in_file-cell_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
    endif
  
    if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
       call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
       ! offset is zero-based
       offset(1) = cell_count
       length(1) = dims(1)
       call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset,length, &
                                  hdf5_err,stride,stride) 
       call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
       call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,cell_ids_i4,dims,hdf5_err, &
                      memory_space_id,file_space_id,prop_id)                     
       call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    endif
    if (option%mycommsize > 1) then
      int_mpi = dims(1)
!geh      call MPI_Bcast(cell_ids,int_mpi,MPIU_INTEGER,ZERO_INTEGER_MPI, &
      call MPI_Bcast(cell_ids_i4,int_mpi,MPI_INTEGER,ZERO_INTEGER_MPI, &
                     option%read_group,ierr)
    endif
  
    call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)

    do i=1,dims(1)
      cell_count = cell_count + 1
      natural_id = cell_ids_i4(i)
      local_ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (local_ghosted_id > 0) then
        local_id = grid%nG2L(local_ghosted_id)
        if (local_id > 0) then
          index_count = index_count + 1
          if (index_count > indices_array_size .and. num_indices <= 0) then
            !reallocate if array grows too large and num_indices <= 0
            allocate(temp(index_count))
            temp(1:index_count) = indices(1:index_count)
            deallocate(indices)
            indices_array_size = 2*indices_array_size
            allocate(indices(indices_array_size))
            indices = 0
            indices(1:index_count) = temp(1:index_count)
            deallocate(temp)
          endif
          indices(index_count) = cell_count
        endif
      endif
    enddo
  enddo

  call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)

  deallocate(cell_ids_i4)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
     call h5pclose_f(prop_id,hdf5_err)
     call h5sclose_f(memory_space_id,hdf5_err)
     call h5sclose_f(file_space_id,hdf5_err)
     call h5dclose_f(data_set_id,hdf5_err)
  endif

  if (num_indices > 0 .and. index_count /= num_indices) then
     write(option%io_buffer, &
          '("Number of indices read (",i9,") does not match the number",&
           &" of indices requested (",i9,").")') &
           index_count, num_indices
     call printErrMsg(option)        
  endif
  
  if (index_count < indices_array_size .and. num_indices <= 0) then
    ! resize to index count
    allocate(temp(index_count))
    temp(1:index_count) = indices(1:index_count)
    deallocate(indices)
    allocate(indices(index_count))
    indices(1:index_count) = temp(1:index_count)
    deallocate(temp)
  endif
  
  if (num_indices <= 0) num_indices = index_count

  call PetscLogEventEnd(logging%event_map_indices_hdf5,ierr);CHKERRQ(ierr)
! End of Vamsi's HDF5 Mechanism  

end subroutine HDF5MapLocToNatIndicesVamsi

! ************************************************************************** !

subroutine HDF5ReadIntegerArrayVamsi(option,file_id,dataset_name, &
                                     dataset_size, &
                                     indices,num_indices,integer_array)
  ! 
  ! Read in local integer values from hdf5 global file
  ! 
  ! Author: Vamsi Sripathi
  ! Date: ??/??/09
  ! 

  use hdf5
  
  use Grid_module
  use Option_module
  
  implicit none

! Vamsi's HDF5 Mechanism 
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt :: indices(:)
  PetscInt :: num_indices
  PetscInt :: integer_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscInt :: index_count
  PetscInt :: integer_count, prev_integer_count
  integer(HSIZE_T) :: num_integers
  PetscInt :: num_integers_in_file
  PetscInt :: temp_int, i, index
  PetscMPIInt :: int_mpi
  
  integer, allocatable :: integer_buffer_i4(:)
  PetscInt :: read_block_size

  call PetscLogEventBegin(logging%event_read_int_array_hdf5, &
                          ierr);CHKERRQ(ierr)

  read_block_size = HDF5_READ_BUFFER_SIZE
  allocate(integer_buffer_i4(read_block_size))
  dims = 0
  integer_count = 0
  prev_integer_count = 0
  index_count = 0
  memory_space_id = -1
  length = 0
  num_integers_in_file = 0

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
     call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
     call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
     ! should be a rank=1 data space
     call h5sget_simple_extent_npoints_f(file_space_id,num_integers, &
                                         hdf5_err)
     call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
     call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
     rank_mpi = 1
     offset = 0
     stride = 1
     num_integers_in_file = int(num_integers) 
  endif  

  call MPI_Bcast(num_integers_in_file,ONE_INTEGER_MPI,MPIU_INTEGER, &
                 ZERO_INTEGER_MPI,option%read_group,ierr) 
                  
  do i=1,num_indices
    index = indices(i)
    if (index > integer_count) then
      do
        if (index <= integer_count) exit
        temp_int = num_integers_in_file-integer_count
        temp_int = min(temp_int,read_block_size)
        if (dims(1) /= temp_int) then
          if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
          dims(1) = temp_int
          length(1) = dims(1)
        endif
        if (mod(option%myrank,option%hdf5_read_group_size) == 0) then     
           call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
           ! offset is zero-based
           offset(1) = integer_count
           length(1) = dims(1)
           call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                      length,hdf5_err,stride,stride) 
           call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
           call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                          hdf5_err,memory_space_id,file_space_id,prop_id)   
           call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
        endif  
        if (option%mycommsize > 1) then
          int_mpi = dims(1)
          call MPI_Bcast(integer_buffer_i4,int_mpi,MPI_INTEGER, &
                         ZERO_INTEGER_MPI,option%read_group,ierr) 
        endif
        prev_integer_count = integer_count
        integer_count = integer_count + length(1)                  
      enddo
    endif
    integer_array(i) = integer_buffer_i4(index-prev_integer_count)
  enddo

  do
    if (integer_count >= num_integers_in_file) exit
    temp_int = num_integers_in_file-integer_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      length(1) = dims(1)
    endif
    if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
       call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
       ! offset is zero-based
       offset(1) = integer_count
       length(1) = dims(1)
       call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                  length,hdf5_err,stride,stride) 
       call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
       call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)   
       call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    endif 
    if (option%mycommsize > 1) then 
      int_mpi = dims(1)
      call MPI_Bcast(integer_buffer_i4,int_mpi,MPI_INTEGER,ZERO_INTEGER_MPI, &
                     option%read_group,ierr) 
    endif                    
    integer_count = integer_count + length(1)                  
  enddo
  deallocate(integer_buffer_i4)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
     call h5pclose_f(prop_id,hdf5_err)
     if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
     call h5sclose_f(file_space_id,hdf5_err)
     call h5dclose_f(data_set_id,hdf5_err)
  endif

  call PetscLogEventEnd(logging%event_read_int_array_hdf5,ierr);CHKERRQ(ierr)
! End of Vamsi's HDF5 Mechanism  
                          
end subroutine HDF5ReadIntegerArrayVamsi

! ************************************************************************** !

subroutine HDF5WriteStructuredDataSetVamsi(name,array,file_id,data_type,option, &
                                      nx_global,ny_global,nz_global, &
                                      nx_local,ny_local,nz_local, &
                                      istart_local,jstart_local,kstart_local)
  ! 
  ! Writes data from an array into HDF5 file
  ! 
  ! Author: Vamsi Sripathi
  ! Date: ??/??/09
  ! 

  use hdf5
  use Option_module
  
  implicit none
  
  type(option_type) :: option

  character(len=32) :: name
  PetscReal :: array(:)

  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3),mem_dims(3)
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscInt :: istart_local, jstart_local, kstart_local
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscInt :: i, j, k, count, id
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscInt :: ny_local_X_nz_local
  PetscMPIInt :: num_to_write_mpi
  
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscMPIInt :: hdf5_flag

  integer, pointer :: int_array_i4(:)
  PetscReal, pointer :: double_array(:)
  
! Vamsi's HDF5 Write 

  PetscMPIInt, allocatable :: group_count_mpi(:),disp_mpi(:)
  PetscInt, allocatable :: group_xyz(:)
  PetscInt :: xyz(0:6),group_num_to_write
  PetscInt :: group_size
  integer, pointer :: group_int_array_i4(:)
  PetscReal, pointer :: group_double_array(:)

  call PetscLogEventBegin(logging%event_write_struct_dataset_hdf5, &
                          ierr);CHKERRQ(ierr)
#define INVERT  

  group_size = option%write_grp_size                      
  group_num_to_write = 0

  ny_local_X_nz_local = ny_local*nz_local
  num_to_write_mpi = nx_local*ny_local_X_nz_local

  if (num_to_write_mpi > 0) then
    if (data_type == HDF_NATIVE_INTEGER) then
      allocate(int_array_i4(nx_local*ny_local*nz_local))
#ifdef INVERT
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            int_array_i4(id) = int(array(count))
          enddo
        enddo
      enddo
#else
      do i=1,grid%nlmax
        int_array_i4(i) = int(array(i))
      enddo
#endif
    else
#ifdef INVERT
      allocate(double_array(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            double_array(id) = array(count)
          enddo
        enddo
      enddo
#else
      do i=1,grid%nlmax
        double_array(i) = array(i)
      enddo
#endif
    endif
  endif

  if (mod(option%myrank,option%hdf5_write_group_size) == 0) then 
     ! I/O - Masters allocate memory
     ! write(*,'(" My Global Rank = ",i5," My Reader Rank = ",i5," Group Size = ",i5,"dataset is ",A)') option%myrank,option%writers_rank,group_size,name
     allocate(group_xyz(0:(group_size*7)-1))
     allocate(group_count_mpi(0:(group_size-1)))
     allocate(disp_mpi(0:(group_size-1)))
  endif

  xyz(0) = nx_local
  xyz(1) = ny_local
  xyz(2) = nz_local
  xyz(3) = istart_local
  xyz(4) = jstart_local
  xyz(5) = kstart_local
  xyz(6) = num_to_write_mpi

  call MPI_Gather(xyz(0),SEVEN_INTEGER_MPI,MPIU_INTEGER, &
                  group_xyz(0),SEVEN_INTEGER_MPI,MPIU_INTEGER, &
                  ZERO_INTEGER_MPI,option%write_group,ierr)

  if (mod(option%myrank,option%hdf5_write_group_size) == 0) then

     do i=0,group_size-1,1
        group_num_to_write = group_num_to_write + group_xyz(i*7+6)
        group_count_mpi(i) = group_xyz(i*7+6)
     enddo

     disp_mpi(0) = 0
     do i=1,group_size-1,1
        disp_mpi(i) = disp_mpi(i-1) + group_count_mpi(i-1)
     enddo

     if (data_type == HDF_NATIVE_INTEGER) then
        allocate(group_int_array_i4(group_num_to_write))
     else
        allocate(group_double_array(group_num_to_write))
     endif

  endif


  if (data_type == HDF_NATIVE_INTEGER) then
!geh      call MPI_Gatherv(int_array,num_to_write_mpi,MPIU_INTEGER,group_int_array,&
!geh                       group_count_mpi(0),disp_mpi(0),MPIU_INTEGER,ZERO_INTEGER_MPI, &
      call MPI_Gatherv(int_array_i4,num_to_write_mpi,MPI_INTEGER,group_int_array_i4,&
                       group_count_mpi(0),disp_mpi(0),MPI_INTEGER,ZERO_INTEGER_MPI, &
                       option%write_group,ierr)
  else
      call MPI_Gatherv(double_array,num_to_write_mpi,MPI_DOUBLE_PRECISION, &
                       group_double_array,&
                       group_count_mpi(0),disp_mpi(0),MPI_DOUBLE_PRECISION, &
                       ZERO_INTEGER_MPI,option%write_group,ierr)
  endif

  ! write (*,'(" Starting HDF5 stuff in ",A," dataset,my rank = ",i6)') name,option%myrank
  if (mod(option%myrank,option%hdf5_write_group_size) == 0) then
     ! write (*,'(" HDF5 Masters at P-1 in ",A," dataset,my rank = ",i6)') name,option%myrank
  
     ! file space which is a 3D block
     file_space_rank_mpi = 3
#ifndef INVERT
     dims(1) = nx_global
     dims(2) = ny_global
     dims(3) = nz_global
#else
     ! have to trick hdf5 for now with inverted ordering
     dims(3) = nx_global
     dims(2) = ny_global
     dims(1) = nz_global
#endif
     call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

     call h5eset_auto_f(OFF,hdf5_err)
     call h5dopen_f(file_id,name,data_set_id,hdf5_err)
     hdf5_flag = hdf5_err
     call h5eset_auto_f(ON,hdf5_err)
     if (hdf5_flag < 0) then 
        call h5screate_simple_f(file_space_rank_mpi,dims,file_space_id,hdf5_err,dims)
        call h5dcreate_f(file_id,name,data_type,file_space_id, &
                         data_set_id,hdf5_err,prop_id)
     else
        call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
     endif

     call h5pclose_f(prop_id,hdf5_err)
  
     call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
     if (trick_hdf5) then
        call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                                hdf5_err) 
     else
        call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &   ! H5FD_MPIO_COLLECTIVE_F
                                hdf5_err)
     endif
#endif
     ! write (*,'(" HDF5 Masters at P-2 in ",A," dataset,my rank = ",i6)') name,option%myrank

     do i=0,group_size-1,1 
        rank_mpi = 1
        mem_dims = 0
        mem_dims(1) = group_count_mpi(i)
        ! if (num_to_write == 0) dims(1) = 1  ----> Need to handle this exception -- Vamsi.
        call h5screate_simple_f(rank_mpi,mem_dims,memory_space_id,hdf5_err,mem_dims)
        ! if (option%myrank == 0) write (*,'(" Created memory space - ",i4)') i   

     if (hdf5_flag < 0) then 
        call h5screate_simple_f(file_space_rank_mpi,dims,file_space_id,hdf5_err,dims)
     else
        call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
     endif
        ! if (option%myrank == 0) write (*,'(" Created file space id - ",i4)') i  

        ! create the hyperslab
#ifndef INVERT
        start(1) = group_xyz(i*7+3) ! istart_local
        start(2) = group_xyz(i*7+4) ! jstart_local
        start(3) = group_xyz(i*7+5) ! kstart_local
        length(1) = group_xyz(i*7)  ! nx_local
        length(2) = group_xyz(i*7+1) ! ny_local
        length(3) = group_xyz(i*7+2) ! nz_local
#else
        start(3) = group_xyz(i*7+3) ! istart_local
        start(2) = group_xyz(i*7+4) ! jstart_local
        start(1) = group_xyz(i*7+5) ! kstart_local
        length(3) = group_xyz(i*7)  ! nx_local
        length(2) = group_xyz(i*7+1) ! ny_local
        length(1) = group_xyz(i*7+2) ! nz_local
#endif
        !  if (num_to_write == 0) length(1) = 1
        if (group_count_mpi(i) .NE. length(1)*length(2)*length(3)) &
          write (*,'("My Rank is ",i8," Memory space and Hyperslab space do not match!!")') option%myrank
        stride = 1
        call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                                   hdf5_err,stride,stride)
        ! if (option%myrank == 0) write (*,'(" Created  hyperslab - ",i4)') i  

        ! write the data
        if (num_to_write_mpi > 0) then
           if (data_type == HDF_NATIVE_INTEGER) then
              call PetscLogEventBegin(logging%event_h5dwrite_f, &
                                      ierr);CHKERRQ(ierr)
              call h5dwrite_f(data_set_id,data_type,group_int_array_i4(disp_mpi(i)+1),mem_dims, &
                              hdf5_err,memory_space_id,file_space_id,prop_id)
              call PetscLogEventEnd(logging%event_h5dwrite_f, &
                                    ierr);CHKERRQ(ierr)
           else
              call PetscLogEventBegin(logging%event_h5dwrite_f, &
                                      ierr);CHKERRQ(ierr)
              call h5dwrite_f(data_set_id,data_type,group_double_array(disp_mpi(i)+1),mem_dims, &
                              hdf5_err,memory_space_id,file_space_id,prop_id)  
              call PetscLogEventEnd(logging%event_h5dwrite_f, &
                                    ierr);CHKERRQ(ierr)
           endif
        endif
        ! if (option%myrank == 0) write (*,'("  Finished h5dwrite_f - ",i4)') i   
         call h5sclose_f(file_space_id,hdf5_err)
         call h5sclose_f(memory_space_id,hdf5_err)

     enddo
     ! write (*,'(" Finished HDF5 write in ",A," dataset,my rank = ",i6)') name,option%myrank
     ! call MPI_Barrier(option%readers,ierr)
     if (data_type == HDF_NATIVE_INTEGER) then
         deallocate(group_int_array_i4)
     else
         deallocate(group_double_array)
     endif

     deallocate(group_xyz)
     deallocate(group_count_mpi)
     deallocate(disp_mpi)
        !call h5sclose_f(file_space_id,hdf5_err)
        !call h5sclose_f(memory_space_id,hdf5_err)
     call h5pclose_f(prop_id,hdf5_err)
     ! call h5sclose_f(file_space_id,hdf5_err)
     call h5dclose_f(data_set_id,hdf5_err)
  endif

  call PetscLogEventEnd(logging%event_write_struct_dataset_hdf5, &
                        ierr);CHKERRQ(ierr)
                         
end subroutine HDF5WriteStructuredDataSetVamsi

! ************************************************************************** !

subroutine HDF5ReadIndicesVamsi(grid,option,file_id,dataset_name,dataset_size, &
                           indices)
  ! 
  ! Reads cell indices from an hdf5 dataset
  ! 
  ! Author: Vamsi Sripathi
  ! Date: ??/??/09
  ! 

  use hdf5
  
  use Option_module
  use Grid_module
  
  implicit none

! Vamsi's HDF5 Mechanism  
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), group_length(3), stride(3) 
  PetscMPIInt :: rank_mpi
  PetscMPIInt, allocatable :: displacement_mpi(:)  
  PetscInt, allocatable :: glength(:)  
  PetscMPIInt, allocatable :: glength_mpi(:)  
  integer, allocatable :: group_indices_i4(:)
  integer(HSIZE_T) :: num_data
  PetscInt :: num_data_in_file
  PetscInt :: istart, iend, id, i
  PetscInt :: length
  PetscMPIInt :: length_mpi
  PetscInt :: group_size 
  
  call PetscLogEventBegin(logging%event_read_indices_hdf5,ierr);CHKERRQ(ierr)

  istart = 0  ! this will be zero-based
  iend = 0
  group_length = 0
  length = 0
  num_data_in_file = 0
  i = 0
  id = 0                      
  group_size = option%read_grp_size

  ! first determine upper and lower bound on PETSc global array
  call MPI_Exscan(grid%nlmax,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                  MPI_SUM,option%mycomm,ierr)
  call MPI_Scan(grid%nlmax,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                MPI_SUM,option%mycomm,ierr)
  if (iend /= istart + grid%nlmax) then
    call printErrMsg(option,'iend /= istart+grid%nlmax')
  endif

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
     call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
     call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
     ! should be a rank=1 data space
     call h5sget_simple_extent_npoints_f(file_space_id,num_data,hdf5_err)
     if (dataset_size > 0 .and. num_data /= dataset_size) then
         write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_data_in_file,dataset_size
         call printErrMsg(option)   
     else
        dataset_size = num_data
        num_data_in_file = int(num_data)
     endif  
     ! Allocate arrays that would hold the length and offset values for the group
     allocate(glength(0:group_size-1))
     allocate(glength_mpi(0:group_size-1))
     allocate(displacement_mpi(0:group_size-1))
  endif


  call MPI_Bcast(num_data_in_file,ONE_INTEGER_MPI,MPIU_INTEGER, &
                 ZERO_INTEGER_MPI,option%read_group,ierr) 

  length = iend - istart 
  call MPI_Gather(length,ONE_INTEGER_MPI,MPIU_INTEGER,glength(0),ONE_INTEGER_MPI, &
                  MPIU_INTEGER,ZERO_INTEGER_MPI,option%read_group,ierr)
 
  if (istart < num_data_in_file) then
     if (mod(option%myrank,option%hdf5_read_group_size) .NE. 0) then 
         allocate(indices(-1:iend-istart))
         indices(-1) = istart
         indices(0) = iend
     else    
         do i = 0,group_size-1,1
            group_length(1) = group_length(1) + glength(i)
         enddo  
      
         displacement_mpi(0) = 0
         id = glength(0) 

         do i = 1,group_size-1,1
            displacement_mpi(i) = id
            id = id + glength(i)
         enddo

         allocate(group_indices_i4(1:group_length(1)))
         allocate(indices(-1:length))  
         indices(-1) = istart
         indices(0) = iend
         dims = 0
         offset = 0
         dims(1) = group_length(1)
         rank_mpi = 1
         stride = 1
         memory_space_id = -1
         offset(1) = istart
     endif   

    if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
        call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
        call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
        call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,offset, &
                                    group_length(1),hdf5_err,stride,stride)
        call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
        call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,group_indices_i4(1:group_length(1)), &
                       dims,hdf5_err,memory_space_id,file_space_id,prop_id)
        call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    endif
  endif

  length_mpi = length
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) glength_mpi = glength

  call MPI_Scatterv(group_indices_i4(1:group_length(1)),glength_mpi,displacement_mpi, &
                    MPIU_INTEGER,indices(1:length),length_mpi,MPIU_INTEGER, &
                    ZERO_INTEGER_MPI,option%read_group,ierr)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
      deallocate(group_indices_i4)
      deallocate(glength)
      deallocate(glength_mpi) 
      deallocate(displacement_mpi)

      call h5pclose_f(prop_id,hdf5_err)
      call h5sclose_f(memory_space_id,hdf5_err)
      call h5sclose_f(file_space_id,hdf5_err)
      call h5dclose_f(data_set_id,hdf5_err)
  endif

  call PetscLogEventEnd(logging%event_read_indices_hdf5,ierr);CHKERRQ(ierr)
! End of Vamsi's HDF5 Mechansim  
  
end subroutine HDF5ReadIndicesVamsi

! ************************************************************************** !

subroutine HDF5ReadArrayVamsi(discretization,grid,option,file_id,dataset_name, &
                         dataset_size, &
                         indices,global_vec,data_type)
  ! 
  ! Read an hdf5 array into a Petsc Vec
  ! 
  ! Author: Vamsi Sripathi
  ! Date: ??/??/09
  ! 
#include "finclude/petscvec.h"
 use petscvec
  use hdf5
  
  use Option_module
  use Grid_module
  use Discretization_module
  
  implicit none

! Vamsi's HDF5 Mechanism 
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  Vec :: global_vec
  integer(HID_T) :: data_type
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), group_length(3), stride(3)
  PetscMPIInt :: rank_mpi
  integer(HSIZE_T) :: num_data_in_file
  Vec :: natural_vec
  PetscInt :: i, istart, iend
  PetscReal, allocatable :: real_buffer(:)
  PetscReal, allocatable :: real_group_buffer(:)  
  integer, allocatable :: integer_buffer_i4(:)
  integer, allocatable :: integer_group_buffer_i4(:)
  PetscInt, allocatable :: indices0(:)
  PetscMPIInt, allocatable :: glength_mpi(:) 
  PetscMPIInt, allocatable :: displacement_mpi(:) 
  PetscInt :: id
  PetscInt :: length
  PetscMPIInt :: length_mpi
  PetscInt :: group_size

  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)
                          
  istart = 0
  iend = 0
  group_length = 0
  length_mpi = 0
  id = 0
  i = 0
  group_size = option%read_grp_size
  
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then 
     call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
     call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
     ! should be a rank=1 data space
     call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)
     if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
         write(option%io_buffer, &
              '(a," data space dimension (",i9,") does not match the dimensions",&
               &" of the domain (",i9,").")') trim(dataset_name), &
               num_data_in_file,grid%nmax
         call printErrMsg(option)   
     endif
     call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
     call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
     ! must initialize here to avoid error below when closing memory space
     memory_space_id = -1
     rank_mpi = 1
     offset = 0
     stride = 1
     ! Allocate arrays that would hold the length and offset values for the group
     allocate(glength_mpi(0:group_size-1))
     allocate(displacement_mpi(0:group_size-1))
  endif

  if (associated(indices)) then
     istart = indices(-1)
     iend = indices(0)
     length = iend - istart
     ! Gather the length values from the group to rank 0 in each group
     call MPI_Gather(length,ONE_INTEGER_MPI,MPIU_INTEGER,glength_mpi(0), &
                     ONE_INTEGER_MPI,MPIU_INTEGER,ZERO_INTEGER_MPI, &
                     option%read_group,ierr)

     if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
        do i = 0,group_size-1,1
           group_length(1) = group_length(1) + glength_mpi(i)
        enddo
      
        displacement_mpi(0) = 0 
        id = glength_mpi(0)
        do i = 1,group_size-1,1
           displacement_mpi(i) = id
           id = id + glength_mpi(i)
        enddo
      
        dims = 0
        dims(1) = group_length(1) 
        ! offset is zero-based
        offset(1) = istart
      
        call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
        call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                   group_length(1),hdf5_err,stride,stride) 
    
        if (data_type == H5T_NATIVE_DOUBLE) then
           allocate(real_group_buffer(group_length(1)))

           call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
           call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_group_buffer,dims, &
                           hdf5_err,memory_space_id,file_space_id,prop_id)
           call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
        else if (data_type == HDF_NATIVE_INTEGER) then
           allocate(integer_group_buffer_i4(group_length(1)))

           call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
           call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_group_buffer_i4,dims, &
                          hdf5_err,memory_space_id,file_space_id,prop_id)
           call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
        endif
     endif

     length_mpi = length
     if (data_type == H5T_NATIVE_DOUBLE) then
        allocate(real_buffer(length_mpi))
        call MPI_Scatterv(real_group_buffer,glength_mpi,displacement_mpi, &
                          MPI_DOUBLE_PRECISION,real_buffer,length_mpi, &
                          MPI_DOUBLE_PRECISION,ZERO_INTEGER_MPI, &
                          option%read_group,ierr)          
        if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
           deallocate(real_group_buffer)
           deallocate(glength_mpi)
           deallocate(displacement_mpi)
        endif  
     else if (data_type == HDF_NATIVE_INTEGER) then
        allocate(integer_buffer_i4(length_mpi))
        allocate(real_buffer(length_mpi))
!geh        call MPI_Scatterv(integer_group_buffer_i4,glength_mpi,displacement_mpi, &
!geh                          MPIU_INTEGER,integer_buffer,length_mpi,MPIU_INTEGER, &
        call MPI_Scatterv(integer_group_buffer_i4,glength_mpi,displacement_mpi, &
                          MPI_INTEGER,integer_buffer_i4,length_mpi,MPI_INTEGER, &
                          ZERO_INTEGER_MPI,option%read_group,ierr)          
        do i = 1,length_mpi,1
           real_buffer(i) = real(integer_buffer_i4(i))
        enddo
        deallocate(integer_buffer_i4)

        if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
           deallocate(integer_group_buffer_i4)
           deallocate(glength_mpi)
           deallocate(displacement_mpi)
        endif  
     endif  

   ! must convert indices to zero based for VecSetValues
   allocate(indices0(iend-istart))
   indices0 = indices(1:iend-istart)-1
   
   call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)
   call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
   call VecSetValues(natural_vec,iend-istart,indices0, &
                     real_buffer,INSERT_VALUES,ierr);CHKERRQ(ierr)
   deallocate(indices0)
   deallocate(real_buffer)

 endif

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then 
     call h5pclose_f(prop_id,hdf5_err)
     if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
     call h5sclose_f(file_space_id,hdf5_err)
     call h5dclose_f(data_set_id,hdf5_err)
  endif

  call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
  call DiscretizationNaturalToGlobal(discretization,natural_vec,global_vec, &
                                     ONEDOF)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  
  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)
! End of Vamsi's HDF5 Mechanism  
                            
end subroutine HDF5ReadArrayVamsi

#endif !PETSC_HAVE_HDF5

! ************************************************************************** !

subroutine HDF5ReadRegionFromFileVamsi(realization,region,filename)
  ! 
  ! Reads a region from an hdf5 file
  ! 
  ! Author: Vamsi Sripathi
  ! Date: ??/??/09
  ! 
#include "finclude/petscvec.h"
  use petscvec
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_module
  use Option_module
  use Grid_module
  use Region_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  type(region_type) :: region
  character(len=MAXSTRINGLENGTH) :: filename

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
#endif

  PetscInt :: num_indices, i, local_id
  PetscInt, pointer :: indices(:)
  PetscInt, pointer :: integer_array(:)
  
  option => realization%option

#if !defined(PETSC_HAVE_HDF5)
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  patch => realization%patch
  grid => patch%grid

  call PetscLogEventBegin(logging%event_region_read_hdf5,ierr);CHKERRQ(ierr)
                          
  ! create hash table for fast lookup
  call GridCreateNaturalToGhostedHash(grid,option)

#if defined(SCORPIO)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
      option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
      call printMsg(option)
  endif

  filename = trim(filename) // CHAR(0)
  call scorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, &
          file_id, ierr)
  string = '/Regions/' // trim(region%name) // '/Cell Ids' //CHAR(0)
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)

 allocate(indices(grid%nlmax))
  ! Read Cell Ids  
  string = 'Regions' // '/' // trim(region%name) // "Cell Ids" // CHAR(0)
  ! num_indices <= 0 indicates that the array size is uncertain and
  ! the size will be returned in num_indices
  num_indices = -1
  call HDF5MapLocalToNaturalIndices(grid,option,file_id,string,ZERO_INTEGER,indices, &
                                    num_indices)
  allocate(integer_array(num_indices))
  integer_array = 0
  string = '/Regions/' // trim(region%name) // '/Cell Ids' //CHAR(0)
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIntegerArray(option,file_id,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)

  ! convert cell ids from natural to local
  call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)
  do i=1,num_indices
    integer_array(i) = grid%nG2L(GridGetLocalGhostedIdFromHash(grid,integer_array(i))) 
  enddo
  call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)
  region%cell_ids => integer_array
                            
  allocate(integer_array(num_indices))
  integer_array = 0
  string = '/Regions/' // trim(region%name) // '/Face Ids' //CHAR(0)
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)  
  call HDF5ReadIntegerArray(option,file_id,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)
                            
  region%faces => integer_array
  region%num_cells = num_indices
  deallocate(indices)
  nullify(indices)

   if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
       option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
       call printMsg(option)   
   endif
   call scorpio_close_file(file_id, option%ioread_group_id, ierr)

  call GridDestroyHashTable(grid)

! SCORPIO
#else   

  ! initialize fortran hdf5 interface 
  call h5open_f(hdf5_err)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
    option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
    call printMsg(option)
    call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(prop_id,option%readers,MPI_INFO_NULL,hdf5_err) 
#endif
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
    call h5pclose_f(prop_id,hdf5_err)

    ! Open the Regions group
    string = 'Regions' 
    option%io_buffer = 'Opening group: ' // trim(string)
    call printMsg(option)  
    call h5gopen_f(file_id,string,grp_id,hdf5_err)

    ! Open the Regions group
    string = trim(region%name)
    option%io_buffer = 'Opening group: ' // trim(string)
    call printMsg(option)  
    call h5gopen_f(grp_id,string,grp_id2,hdf5_err)
  endif

  allocate(indices(grid%nlmax))
 ! Read Cell Ids
  string = "Cell Ids"
  ! num_indices <= 0 indicates that the array size is uncertain and
  ! the size will be returned in num_indices
  num_indices = -1
  call HDF5MapLocalToNaturalIndices(grid,option,grp_id2,string,ZERO_INTEGER,indices, &
                                    num_indices)
  allocate(integer_array(num_indices))
  integer_array = 0
  string = "Cell Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIntegerArray(option,grp_id2,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)

  ! convert cell ids from natural to local
  call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)
  do i=1,num_indices
    integer_array(i) = grid%nG2L(GridGetLocalGhostedIdFromHash(grid,integer_array(i))) 
  enddo
  call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)
  region%cell_ids => integer_array
                            
  allocate(integer_array(num_indices))
  integer_array = 0
  string = "Face Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)  
  call HDF5ReadIntegerArray(option,grp_id2,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)
                            
  region%faces => integer_array
  region%num_cells = num_indices
  deallocate(indices)
  nullify(indices)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
    option%io_buffer = 'Closing group: ' // trim(region%name)
    call printMsg(option)  
    call h5gclose_f(grp_id2,hdf5_err)
    option%io_buffer = 'Closing group: Regions'
    call printMsg(option)   
    call h5gclose_f(grp_id,hdf5_err)
    option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
    call printMsg(option)   
    call h5fclose_f(file_id,hdf5_err)
  endif
  call h5close_f(hdf5_err)

  call GridDestroyHashTable(grid)

#endif !PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_region_read_hdf5,ierr);CHKERRQ(ierr)

end subroutine HDF5ReadRegionFromFileVamsi

! ************************************************************************** !

subroutine HDF5ReadCellIndxIntArrayVamsi(realization,global_vec,filename, &
                                           group_name, &
                                           dataset_name,append_realization_id)
  ! 
  ! Reads an array of integer values from an
  ! hdf5 file
  ! 
  ! Author: Vamsi Sripathi
  ! Date: ??/??/09
  ! 

#include "finclude/petscvec.h"
  use petscvec
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  Vec :: global_vec
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscInt, allocatable :: integer_array(:)
  
#if !defined(PETSC_HAVE_HDF5)
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_cell_indx_int_read_hdf5, &
                          ierr);CHKERRQ(ierr)
  
 ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
    option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
    call printMsg(option) 
    call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(prop_id,option%readers,MPI_INFO_NULL,hdf5_err) 
#endif
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
    call h5pclose_f(prop_id,hdf5_err)

    option%io_buffer = 'Setting up grid cell indices'
    call printMsg(option) 

    ! Open group if necessary
    if (len_trim(group_name) > 1) then
       option%io_buffer = 'Opening group: ' // trim(group_name)
       call printMsg(option)   
       call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
    else
       grp_id = file_id
    endif
  endif

  ! Read Cell Ids
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = "Cell Ids"
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)


  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif
  string = trim(dataset_name) // adjustl(trim(string))
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,grp_id,string,grid%nmax, &
                     indices,global_vec,HDF_NATIVE_INTEGER)
  
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to read integer array.")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
    if (file_id /= grp_id) then
      option%io_buffer = 'Closing group: ' // trim(group_name)
      call printMsg(option)   
      call h5gclose_f(grp_id,hdf5_err)
    endif
    option%io_buffer = 'Closing hdf5 file: ' // filename
    call printMsg(option)   
    call h5fclose_f(file_id,hdf5_err)
  endif
  call h5close_f(hdf5_err)

#endif  ! PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_cell_indx_int_read_hdf5, &
                        ierr);CHKERRQ(ierr)
                          
end subroutine HDF5ReadCellIndxIntArrayVamsi

! ************************************************************************** !

subroutine HDF5ReadCellIndxRealArrayVamsi(realization,global_vec,filename, &
                                        group_name, &
                                        dataset_name,append_realization_id)
  ! 
  ! Reads an array of real values from an hdf5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/09, 02/18/09
  ! 

#include "finclude/petscvec.h"
  use petscvec
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  Vec :: global_vec
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscReal, allocatable :: real_array(:)
  
#if !defined(PETSC_HAVE_HDF5)
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_cell_indx_real_read_hdf5, &
                          ierr);CHKERRQ(ierr)

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
    option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
    call printMsg(option) 
    call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
     call h5pset_fapl_mpio_f(prop_id,option%readers,MPI_INFO_NULL,hdf5_err) 
#endif
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
    call h5pclose_f(prop_id,hdf5_err)

    option%io_buffer = 'Setting up grid cell indices'
    call printMsg(option) 

    ! Open group if necessary
    if (len_trim(group_name) > 1) then
      option%io_buffer = 'Opening group: ' // trim(group_name)
      call printMsg(option)   
      call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
    else
      grp_id = file_id
    endif
  endif

  ! Read Cell Ids
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = "Cell Ids"
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif
  string = trim(dataset_name) // adjustl(trim(string))
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,file_id,string,grid%nmax, &
                     indices,global_vec,H5T_NATIVE_DOUBLE)
  
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to read real array")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
    if (file_id /= grp_id) then
      option%io_buffer = 'Closing group: ' // trim(group_name)
      call printMsg(option)   
      call h5gclose_f(grp_id,hdf5_err)
    endif
    option%io_buffer = 'Closing hdf5 file: ' // filename
    call printMsg(option)   
    call h5fclose_f(file_id,hdf5_err)
 endif  
  call h5close_f(hdf5_err)

#endif ! PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_cell_indx_real_read_hdf5, &
                        ierr);CHKERRQ(ierr)
                          
end subroutine HDF5ReadCellIndxRealArrayVamsi

#endif !if defined(PETSC_HAVE_HDF5)
      
end module HDF5_Vamsi_module



