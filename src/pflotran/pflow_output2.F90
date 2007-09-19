module pflow_output2_module

  use pflow_gridtype_module

  implicit none
  
  private

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petsclog.h"

#define X_COORDINATE 1
#define Y_COORDINATE 2
#define Z_COORDINATE 3
#define TEMPERATURE 4
#define PRESSURE 5
#define LIQUID_SATURATION 6
#define GAS_SATURATION 7
#define LIQUID_ENERGY 8
#define GAS_ENERGY 9
#define LIQUID_MOLE_FRACTION 10
#define GAS_MOLE_FRACTION 11
#define VOLUME_FRACTION 12
#define PHASE 13

#define TECPLOT_INTEGER 0
#define TECPLOT_REAL 1

#define TECPLOT_FILE 0
#define HDF5_FILE 1

  integer :: hdf5_err
  logical, save :: hdf5_first = .true.
  PetscErrorCode :: ierr
  
  public :: OutputTecplot, OutputHDF5
  
contains
  
subroutine OutputTecplot(grid,kplot)
 
  implicit none

#include "definitions.h"

  type(pflowGrid) :: grid
  integer :: kplot
  
  integer :: i
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  Vec :: global
  Vec :: natural
  
  ! open file
  if (kplot < 10) then
    write(filename,'("pflow00",i1,".tec")') kplot  
  else if (kplot < 100) then
    write(filename,'("pflow0",i2,".tec")') kplot  
  else if (kplot < 1000) then
    write(filename,'("pflow",i3,".tec")') kplot  
  else if (kplot < 10000) then
    write(filename,'("pflow",i4,".tec")') kplot  
  endif
  
  if (grid%myrank == 0) then
    print *, '--> write tecplot output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 grid%t/grid%tconv, grid%tunit
    ! write variables
    if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
        grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE .or. &
        grid%use_richards == PETSC_TRUE) then
      string = 'VARIABLES=' // &
               '"X-coordinates",' // &
               '"Y-coordinates",' // &
               '"Z-coordinates",' // &
               '"Temperature",' // &
               '"Pressure",' // &
               '"Liquid Saturation",' // &
               '"Gas Saturation",' // &
               '"Liquid Energy",' // &
               '"Gas Energy",'
      do i=1,grid%nspec
        write(string2,'(''"Liquid Mole Fraction('',i2,'')",'')') i
        string = trim(string) // trim(string2)
      enddo
      do i=1,grid%nspec
        write(string2,'(''"Gas Mole Fraction('',i2,'')",'')') i
        string = trim(string) // trim(string2)
      enddo
      if (grid%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction",'
      endif
      string = trim(string) // '"Phase"'
    else
      string = '"X-coordinates",' // &
               '"Y-coordinates",' // &
               '"Z-coordinates",' // &
               '"Temperature",' // &
               '"Pressure",' // &
               '"Saturation",' // &
               '"Concentration",'
      if (grid%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction",'
      endif
      string = trim(string) // '"Phase"'
    endif
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') &
                 grid%t/grid%tconv,grid%nx,grid%ny,grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call DACreateGlobalVector(grid%da_1_dof,global,ierr)
  call DACreateNaturalVector(grid%da_1_dof,natural,ierr)

  ! write out coorindates
  call GetCoordinates(grid,global,X_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Y_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Z_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE .or. &
      grid%use_richards == PETSC_TRUE) then

    ! temperature
    call GetVarFromArray(grid,global,TEMPERATURE,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! pressure
    call GetVarFromArray(grid,global,PRESSURE,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! liquid saturation
    call GetVarFromArray(grid,global,LIQUID_SATURATION,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! gas saturation
    call GetVarFromArray(grid,global,GAS_SATURATION,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    
    ! liquid energy
    call GetVarFromArray(grid,global,LIQUID_ENERGY,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    
    ! gas energy
    call GetVarFromArray(grid,global,GAS_ENERGY,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    
    ! liquid mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,LIQUID_MOLE_FRACTION,i-1)
      call ConvertGlobalToNatural(grid,global,natural)
      call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    enddo
    
    ! gas mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,GAS_MOLE_FRACTION,i-1)
      call ConvertGlobalToNatural(grid,global,natural)
      call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    enddo
    
    ! Volume Fraction
    if (grid%rk > 0.d0) then
      call GetVarFromArray(grid,global,VOLUME_FRACTION,0)
      call ConvertGlobalToNatural(grid,global,natural)
      call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    endif
    
    ! phase
    call GetVarFromArray(grid,global,PHASE,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_INTEGER)
  
  else
  
    ! temperature
    call ConvertGlobalToNatural(grid,grid%temp,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! pressure
    call ConvertGlobalToNatural(grid,grid%pressure,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! saturation
    call ConvertGlobalToNatural(grid,grid%sat,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! concentration
    call ConvertGlobalToNatural(grid,grid%conc,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  endif
  
  call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  close(IUNIT3)
      
end subroutine OutputTecplot

subroutine WriteTecplotDataSetFromVec(fid,grid,vec,datatype)

  implicit none

  integer :: fid
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: datatype
  
  PetscScalar, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,grid,vec_ptr,datatype)
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteTecplotDataSetFromVec

subroutine WriteTecplotDataSet(fid,grid,vec_ptr,datatype)

  implicit none
  
  integer :: fid
  type(pflowGrid) :: grid
  PetscReal, pointer :: vec_ptr(:)
  integer, save :: max_local_size = -1
  integer :: datatype
  
  integer :: i, iproc, recv_size
  integer :: istart, iend, num_in_array
  integer :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)
  
  ! if first time, determine the maximum size of any local array across all procs
  if (max_local_size < 0) then
    call MPI_Allreduce(grid%nlmax,max_local_size,1,MPI_INTEGER,MPI_MAX, &
                       PETSC_COMM_WORLD,ierr)
    if (grid%myrank == 0) print *, 'max_local_size: ', max_local_size
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == TECPLOT_INTEGER) then
    allocate(integer_data(max_local_size+10))
    allocate(integer_data_recv(max_local_size))
    do i=1,grid%nlmax
      integer_data(i) = int(vec_ptr(i))
    enddo
  else
    allocate(real_data(max_local_size+10))
    allocate(real_data_recv(max_local_size))
    do i=1,grid%nlmax
      real_data(i) = vec_ptr(i)
    enddo
  endif
  
  ! communicate data to processor 0, round robin style
  if (grid%myrank == 0) then
    if (datatype == TECPLOT_INTEGER) then
      ! This approach makes output files identical, regardless of processor
      ! distribution.  It is necessary when diffing files.
      iend = 0
      do
        istart = iend+1
        if (iend+10 > grid%nlmax) exit
        iend = istart+9
        write(IUNIT3,'(10(i3,x))') integer_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      integer_data(1:grid%nlmax-iend) = integer_data(iend+1:grid%nlmax)
      num_in_array = grid%nlmax-iend
    else
      iend = 0
      do
        istart = iend+1
        if (iend+10 > grid%nlmax) exit
        iend = istart+9
        write(IUNIT3,'(10(es11.4,x))') real_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      real_data(1:grid%nlmax-iend) = real_data(iend+1:grid%nlmax)
      num_in_array = grid%nlmax-iend
    endif
    do iproc=1,grid%commsize-1
      call MPI_Probe(iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
      recv_size = status(MPI_TAG)
      if (datatype == 0) then
        call MPI_Recv(integer_data_recv,recv_size,MPI_INTEGER,iproc, &
                      MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
        if (recv_size > 0) then
          integer_data(num_in_array+1:num_in_array+recv_size) = &
                                             integer_data_recv(1:recv_size)
          num_in_array = num_in_array+recv_size
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(IUNIT3,'(10(i3,x))') integer_data(istart:iend)
        enddo
        if (iend > 0) then
          integer_data(1:num_in_array-iend) = integer_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      else
        call MPI_Recv(real_data_recv,recv_size,MPI_DOUBLE_PRECISION,iproc, &
                      MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
        if (recv_size > 0) then
          real_data(num_in_array+1:num_in_array+recv_size) = &
                                             real_data_recv(1:recv_size)
          num_in_array = num_in_array+recv_size
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(IUNIT3,'(10(es11.4,x))') real_data(istart:iend)
        enddo
        if (iend > 0) then
          real_data(1:num_in_array-iend) = real_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      endif
    enddo
    ! Print the remaining values, if they exist
    if (datatype == 0) then
      if (num_in_array > 0) &
        write(IUNIT3,'(10(i3,x))') integer_data(1:num_in_array)
    else
      if (num_in_array > 0) &
        write(IUNIT3,'(10(es11.4,x))') real_data(1:num_in_array)
    endif
  else
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,grid%nlmax,MPI_INTEGER,0,grid%nlmax, &
                    PETSC_COMM_WORLD,ierr)
    else
      call MPI_Send(real_data,grid%nlmax,MPI_DOUBLE_PRECISION,0,grid%nlmax, &
                    PETSC_COMM_WORLD,ierr)
    endif
  endif
      
  if (datatype == TECPLOT_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

end subroutine WriteTecplotDataSet

#ifndef USE_HDF5
subroutine OutputHDF5(grid)

  implicit none
  
  type(pflowGrid) :: grid

  if (grid%myrank == 0) then
    print *
    print *, 'PFLOTRAN must be compiled with -DUSE_HDF5 to ', &
             'write to an HDF5 format.'
    print *
  endif
  stop

end subroutine OutputHDF5

#else
subroutine OutputHDF5(grid)

  use hdf5
  
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
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL,hdf5_err);
#endif
  if (.not.hdf5_first) call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id, &
                                      hdf5_err,prop_id)
  if (hdf5_err < 0 .or. hdf5_first) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                     prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (hdf5_first) then
    if (grid%myrank == 0) print *, '--> creating hdf5 output file: ', filename
  else
    if (grid%myrank == 0) print *, '--> appending to hdf5 output file: ', &
                                   filename
  endif
  
  if (hdf5_first) then

    ! create a group for the coordinates data set
    string = "Coordinates"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  
    ! write out coordinates in x, y, and z directions
    string = "X-coordinates"
    allocate(array(grid%nx))
    do i=1,grid%nx
      if (i == 1) then
        array(i) = 0.5d0*grid%dx0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%dx0(i-1)+grid%dx0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,grid,grid%nx,array,grp_id)
    deallocate(array)
  
    string = "Y-coordinates"
    allocate(array(grid%ny))
    do i=1,grid%ny
      if (i == 1) then
        array(i) = 0.5d0*grid%dy0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%dy0(i-1)+grid%dy0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,grid,grid%ny,array,grp_id)
    deallocate(array)
  
    string = "Z-coordinates"
    allocate(array(grid%nz))
    do i=1,grid%nz
      if (i == 1) then
        array(i) = 0.5d0*grid%dz0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%dz0(i-1)+grid%dz0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,grid,grid%nz,array,grp_id)
    deallocate(array)
    
    call h5gclose_f(grp_id,hdf5_err)
    
  endif

  ! create a group for the data set
  write(string,'('' Time('',i4,''):'',es12.4,x,a1)') &
        grid%flowsteps,grid%t/grid%tconv,grid%tunit
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  
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
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)

    ! pressure
    call GetVarFromArray(grid,global,PRESSURE,0)
    string = "Pressure"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)

    ! liquid saturation
    call GetVarFromArray(grid,global,LIQUID_SATURATION,0)
    string = "Liquid Saturation"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)  

    ! gas saturation
    call GetVarFromArray(grid,global,GAS_SATURATION,0)
    string = "Gas Saturation"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! liquid energy
    call GetVarFromArray(grid,global,LIQUID_ENERGY,0)
    string = "Liquid Energy"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! gas energy
    call GetVarFromArray(grid,global,GAS_ENERGY,0)
    string = "Gas Energy"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! liquid mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,LIQUID_MOLE_FRACTION,i-1)
      write(string,'(''Liquid Mole Fraction('',i4,'')'')') i
      call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)
    enddo
    
    ! gas mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,GAS_MOLE_FRACTION,i-1)
      write(string,'(''Gas Mole Fraction('',i4,'')'')') i
      call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)
    enddo
    
    ! Volume Fraction
    if (grid%rk > 0.d0) then
      call GetVarFromArray(grid,global,VOLUME_FRACTION,0)
      string = "Volume Fraction"
      call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)
    endif
    
    ! phase
    call GetVarFromArray(grid,global,PHASE,0)
    string = "Phase"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_INTEGER) 
  
  else
  
    ! temperature
    string = "Temperature"
    call WriteHDF5DataSetFromVec(string,grid,grid%temp,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    ! pressure
    string = "Pressure"
    call WriteHDF5DataSetFromVec(string,grid,grid%pressure,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    ! saturation
    string = "Saturation"
    call WriteHDF5DataSetFromVec(string,grid,grid%sat,grp_id,H5T_NATIVE_DOUBLE)

    ! concentration
    string = "Concentration"
    call WriteHDF5DataSetFromVec(string,grid,grid%conc,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  
  endif
    
  ! call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  hdf5_first = .false.

end subroutine OutputHDF5

subroutine WriteHDF5Coordinates(name,grid,length,array,file_id)

  use hdf5
  
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
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims);
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err); ! must be independent and only from p0
#endif
  if (grid%myrank == 0) then
    call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                    hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

end subroutine WriteHDF5Coordinates

subroutine WriteHDF5DataSetFromVec(name,grid,vec,file_id,data_type)

  use hdf5
  
  implicit none

  character(len=32) :: name
  type(pflowGrid) :: grid
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  PetscScalar, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteHDF5DataSet(name,grid,vec_ptr,file_id,data_type)
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteHDF5DataSetFromVec

subroutine WriteHDF5DataSet(name,grid,array,file_id,data_type)

  use hdf5
  
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
  call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims);

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
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims);


  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,data_type,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)
  
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
                             hdf5_err,stride,stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F,hdf5_err); ! must be independent and only from p0
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
                    hdf5_err,memory_space_id,file_space_id,prop_id)
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
                    hdf5_err,memory_space_id,file_space_id,prop_id)  
    deallocate(double_array)
#else
    call h5dwrite_f(data_set_id,data_type,array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)  
#endif
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

end subroutine WriteHDF5DataSet
#endif

subroutine GetCoordinates(grid,vec,direction)

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: direction
  
  integer :: i
  PetscScalar, pointer :: vec_ptr(:)
  
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
  
end subroutine GetCoordinates

subroutine ConvertGlobalToNatural(grid,global,natural)

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: global
  Vec :: natural
  
  call DAGlobalToNaturalBegin(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)
  call DAGlobalToNaturalEnd(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)

end subroutine 

subroutine GetVarFromArray(grid,vec,ivar,isubvar)

  use pflow_gridtype_module

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: ivar
  integer :: isubvar

  integer :: i
  integer :: offset, saturation_offset
  integer :: size_var_use
  integer :: size_var_node
  PetscScalar, pointer :: var_ptr(:)
  PetscScalar, pointer :: vec_ptr(:)

  call VecGetArrayF90(vec,vec_ptr,ierr)
      
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION)
      select case(ivar)
        case(TEMPERATURE)
          offset = 1
        case(PRESSURE)
          offset = 2
        case(LIQUID_SATURATION)
          offset = 3
        case(GAS_SATURATION)
          offset = 4
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

    case(LIQUID_ENERGY,GAS_ENERGY)

      select case (ivar)
        case(LIQUID_ENERGY)
          offset = 11
          saturation_offset = 3
        case(GAS_ENERGY)
          offset = 12
          saturation_offset = 4  
      end select

      size_var_use = 2 + 7*grid%nphase + 2* grid%nphase*grid%nspec
      size_var_node = (grid%ndof + 1) * size_var_use
        
      call VecGetArrayF90(grid%var,var_ptr,ierr)
      do i=1,grid%nlmax
        if (var_ptr((i-1)*size_var_node+saturation_offset) > 1.d-30) then
          vec_ptr(i) = var_ptr((i-1)*size_var_node+offset)
        else
          vec_ptr(i) = 0.d0
        endif
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
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

end subroutine GetVarFromArray

end module pflow_output2_module

    
       
