module hdf5_output_module

! if USE_HDF5 is not defined, this is merely a dummy module
#ifdef USE_HDF5

  use pflow_gridtype_module
  use hdf5

  implicit none

  private

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petsclog.h"

#define TEMPERATURE 1
#define PRESSURE 2
#define LIQUID_SATURATION 3
#define GAS_SATURATION 4
#define LIQUID_ENERGY 5
#define GAS_ENERGY 6
#define LIQUID_MOLE_FRACTION 7
#define GAS_MOLE_FRACTION 8
#define VOLUME_FRACTION 9
#define PHASE 10

  integer :: hdferr
  logical, save :: first = .true.
  PetscErrorCode :: ierr

  public :: OutputHDF5
  
contains

subroutine OutputHDF5(grid)

  implicit none

  type(pflowGrid) :: grid

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: rank
  integer(HSIZE_T) :: dims(3)
  
  Vec :: global
  Vec :: natural
  PetscScalar, pointer :: v_ptr
  
  character(len=32) :: filename = "pflow.h5"
  character(len=128) :: string
  real*8, pointer :: array(:)
  integer :: i

  ! initialize fortran interface
  call h5open_f(hdferr)

!  string = "pflow.h5"
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdferr)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL,hdferr);
#endif
  if (.not.first) call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdferr,prop_id)
  if (hdferr < 0 .or. first) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdferr,H5P_DEFAULT_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdferr)
  
  if (first) then

    ! create a group for the coordinates data set
    string = "Coordinates"
    call h5gcreate_f(file_id,string,grp_id,hdferr,OBJECT_NAMELEN_DEFAULT_F)
  
    ! write out coordinates in x, y, and z directions
    string = "X-coordinates"
    allocate(array(grid%nx))
    array = grid%x(1:grid%nx)
    call WriteCoordinate(string,grid,grid%nx,array,grp_id)
    deallocate(array)
  
    string = "Y-coordinates"
    allocate(array(grid%ny))
    do i=1,grid%ny
      array(i) = grid%y(i*grid%nx)
    enddo
    call WriteCoordinate(string,grid,grid%ny,array,grp_id)
    deallocate(array)
  
    string = "Z-coordinates"
    allocate(array(grid%nz))
    do i=1,grid%nz
      array(i) = grid%z(i*grid%nx*grid%ny)
    enddo
    call WriteCoordinate(string,grid,grid%nz,array,grp_id)
    deallocate(array)
    
    call h5gclose_f(grp_id,hdferr)
    
  endif

  ! create a group for the data set
  write(string,'('' Time('',i4,''):'',pg12.4,x,a1)') &
        grid%flowsteps,grid%t/grid%tconv,grid%tunit
  call h5gcreate_f(file_id,string,grp_id,hdferr,OBJECT_NAMELEN_DEFAULT_F)
  
  ! write out data sets  
  call DACreateGlobalVector(grid%da_1_dof,global,ierr)
  !call DACreateNaturalVector(grid%da_1_dof,natural,ierr)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE .or. &
      grid%use_richards == PETSC_TRUE) then
  
    ! temperature
    call GetVarFromArray(grid,global,TEMPERATURE,0)
    !call DAGlobalToNaturalBegin(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)
    !call DAGlobalToNaturalEnd(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)
    string = "Temperature"
    call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)

    ! pressure
    call GetVarFromArray(grid,global,PRESSURE,0)
    string = "Pressure"
    call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)

    ! liquid saturation
    call GetVarFromArray(grid,global,LIQUID_SATURATION,0)
    string = "Liquid Saturation"
    call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)  

    ! gas saturation
    call GetVarFromArray(grid,global,GAS_SATURATION,0)
    string = "Gas Saturation"
    call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! liquid energy
    call GetVarFromArray(grid,global,LIQUID_ENERGY,0)
    string = "Liquid Energy"
    call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! gas energy
    call GetVarFromArray(grid,global,GAS_ENERGY,0)
    string = "Gas Energy"
    call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! liquid mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,LIQUID_MOLE_FRACTION,i-1)
      write(string,'(''Liquid Mole Fraction('',i4,'')'')') i
      call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    enddo
    
    ! gas mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,GAS_MOLE_FRACTION,i-1)
      write(string,'(''Gas Mole Fraction('',i4,'')'')') i
      call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    enddo
    
    ! Volume Fraction
    if (grid%rk > 0.d0) then
      call GetVarFromArray(grid,global,VOLUME_FRACTION,0)
      string = "Volume Fraction"
      call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    endif
    
    ! phase
    call GetVarFromArray(grid,global,PHASE,0)
    string = "Phase"
    call WriteDataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_INTEGER) 
  
  else
  
    ! temperature
    string = "Temperature"
    call WriteDataSetFromVec(string,grid,grid%temp,grp_id,H5T_NATIVE_DOUBLE)

    ! pressure
    string = "Pressure"
    call WriteDataSetFromVec(string,grid,grid%pressure,grp_id,H5T_NATIVE_DOUBLE)

    ! saturation
    string = "Saturation"
    call WriteDataSetFromVec(string,grid,grid%sat,grp_id,H5T_NATIVE_DOUBLE)

    ! concentration
    string = "Concentration"
    call WriteDataSetFromVec(string,grid,grid%conc,grp_id,H5T_NATIVE_DOUBLE)
  
  endif
    
  ! call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  call h5gclose_f(grp_id,hdferr)
  call h5fclose_f(file_id,hdferr)
  call h5close_f(hdferr)
  first = .false.

end subroutine OutputHDF5

subroutine WriteCoordinate(name,grid,length,array,file_id)

  implicit none
  
  character(len=32) :: name
  type(pflowGrid) :: grid
  integer :: length
  real*8 :: array(:)
  integer(HID_T) :: file_id
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer :: rank
  
  ! write out grid structure
  rank = 1
  dims = 0
  ! x-direction
  dims(1) = length
  call h5screate_simple_f(rank,dims,file_space_id,hdferr,dims);
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdferr)
  call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,file_space_id, &
                   data_set_id,hdferr,prop_id)
  call h5pclose_f(prop_id,hdferr)

  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdferr)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdferr); ! must be independent and only from p0
#endif
  if (grid%myrank == 0) then
    call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                    hdferr,H5S_ALL_F,H5S_ALL_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdferr)
  call h5dclose_f(data_set_id,hdferr)
  call h5sclose_f(file_space_id,hdferr)

end subroutine WriteCoordinate
  
subroutine WriteDataSetFromVec(name,grid,vector,file_id,data_type)

  implicit none

  character(len=32) :: name
  type(pflowGrid) :: grid
  Vec :: vector
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  PetscScalar, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  call WriteDataSet(name,grid,vec_ptr,file_id,data_type)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine WriteDataSetFromVec
  
subroutine WriteDataSet(name,grid,array,file_id,data_type)

  implicit none
  
  character(len=32) :: name
  type(pflowGrid) :: grid
  real*8 :: array(:)
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer :: rank
  
  integer, pointer :: int_array(:)
  real*8, pointer :: double_array(:)
  integer :: i, j, k, count, id
  integer(HSIZE_T) :: start(3), length(3), stride(3)

  ! memory space which is a 1D vector  
  rank = 1
  dims = 0
  dims(1) = grid%nlmax
  call h5screate_simple_f(rank,dims,memory_space_id,hdferr,dims);

  ! file space which is a 3D block
  rank = 3
#define INVERT
#ifndef INVERT
  dims(1) = grid%nx
  dims(2) = grid%ny
  dims(3) = grid%nz
#else
! have to trick hdf5 for now with inverted ordering
  dims(3) = grid%nx
  dims(2) = grid%ny
  dims(1) = grid%nz
#endif
  call h5screate_simple_f(rank,dims,file_space_id,hdferr,dims);


  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdferr)
  call h5dcreate_f(file_id,name,data_type,file_space_id, &
                   data_set_id,hdferr,prop_id)
  call h5pclose_f(prop_id,hdferr)
  
  ! create the hyperslab
#ifndef INVERT
  start(1) = grid%nxs
  start(2) = grid%nys
  start(3) = grid%nzs
  length(1) =  grid%nlx
  length(2) =  grid%nly
  length(3) =  grid%nlz
#else
  start(3) = grid%nxs
  start(2) = grid%nys
  start(1) = grid%nzs
  length(3) =  grid%nlx
  length(2) =  grid%nly
  length(1) =  grid%nlz
#endif
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdferr,stride,stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdferr)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F,hdferr); ! must be independent and only from p0
#endif
  if (data_type == H5T_NATIVE_INTEGER) then
    allocate(int_array(grid%nlmax))
#ifdef INVERT
    count = 0
    do k=1,grid%nlz
      do j=1,grid%nly
        do i=1,grid%nlx
          id = k+(j-1)*grid%nlz+(i-1)*grid%nlyz
          count = count+1
          int_array(id) = int(array(count))
        enddo
      enddo
    enddo
#else
    do i=1,grid%nlmax
      int_array(i) = int(array(i))
    enddo
#endif
    call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                    hdferr,memory_space_id,file_space_id,prop_id)
    deallocate(int_array)
  else
#ifdef INVERT
    allocate(double_array(grid%nlmax))
    count = 0
    do k=1,grid%nlz
      do j=1,grid%nly
        do i=1,grid%nlx
          id = k+(j-1)*grid%nlz+(i-1)*grid%nlyz
          count = count+1
          double_array(id) = array(count)
        enddo
      enddo
    enddo
    call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                    hdferr,memory_space_id,file_space_id,prop_id)  
    deallocate(double_array)
#else
    call h5dwrite_f(data_set_id,data_type,array,dims, &
                    hdferr,memory_space_id,file_space_id,prop_id)  
#endif
  endif
  call h5pclose_f(prop_id,hdferr)
  call h5dclose_f(data_set_id,hdferr)
  call h5sclose_f(memory_space_id,hdferr)
  call h5sclose_f(file_space_id,hdferr)

end subroutine WriteDataSet

subroutine GetVarFromArray(grid,vector,ivar,isubvar)

  use pflow_gridtype_module

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: vector
  integer :: ivar
  integer :: isubvar

  integer :: i
  integer :: offset
  integer :: size_var_use
  integer :: size_var_node
  PetscScalar, pointer :: var_ptr(:)
  PetscScalar, pointer :: vec_ptr(:)

  call VecGetArrayF90(vector,vec_ptr,ierr)
      
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_ENERGY,GAS_ENERGY,LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION)
      select case(ivar)
        case(TEMPERATURE)
          offset = 1
        case(PRESSURE)
          offset = 2
        case(LIQUID_SATURATION)
          offset = 3
        case(GAS_SATURATION)
          offset = 4
        case(LIQUID_ENERGY)
          offset = 11
        case(GAS_ENERGY)
          offset = 12    
        case(LIQUID_MOLE_FRACTION)
          offset = 17+isubvar
        case(GAS_MOLE_FRACTION)
          offset = 17+grid%nspec+isubvar
      end select
    
      size_var_use = 2 + 7*grid%nphase + 2* grid%nphase*grid%nspec
      size_var_node = (grid%ndof + 1) * size_var_use
        
      call VecGetArrayF90(grid%var,var_ptr,ierr)
      do i=1,grid%nlmax
        vec_ptr(i) = var_ptr((i-1)*size_var_node+offset)
      enddo
      call VecRestoreArrayF90(grid%var,var_ptr,ierr)

    case(VOLUME_FRACTION)
    
      ! need to set minimum to 0.
      call VecGetArrayF90(grid%phis,var_ptr,ierr)
      vec_ptr(1:grid%nlmax) = var_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(grid%phis,var_ptr,ierr)
     
    case(PHASE)
    
      call VecGetArrayF90(grid%iphas,var_ptr,ierr)
      vec_ptr(1:grid%nlmax) = var_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(grid%iphas,var_ptr,ierr)
     
  end select
  
  call VecRestoreArrayF90(vector,vec_ptr,ierr)

end subroutine GetVarFromArray
 
#endif

end module hdf5_output_module
