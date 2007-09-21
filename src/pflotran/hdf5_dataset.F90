module dataset_write_module

public :: WriteHDF5DataSet

contains

subroutine WriteHDF5DataSet(string,int_array,real_array,length,file_id, &
                            data_type)

  use hdf5

  implicit none
  
  character(len=512) :: string
  integer :: int_array(:)
  real*8 :: real_array(:)
  integer :: length
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3), compress_dim(3)
  integer :: rank
  integer :: hdf5_err

  rank = 1
  dims = 0
  dims(1) = length
  ! create the file space (layout in hdf5 file)
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims)

  ! create the data set
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

#define COMPRESS  
#ifdef COMPRESS
  compress_dim(1) = 32768
  if (length > 2*compress_dim(1)) then
    call h5pset_chunk_f(prop_id,1,compress_dim,hdf5_err)
    call h5pset_shuffle_f(prop_id,hdf5_err)
    call h5pset_szip_f(prop_id,H5_SZIP_NN_OM_F,16,hdf5_err)
    call h5pset_deflate_f(prop_id,9,hdf5_err)
  endif
#endif  
  
  call h5dcreate_f(file_id,trim(string),data_type,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)
  
  ! write the data 
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
  if (data_type == H5T_NATIVE_INTEGER) then
    call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                     hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
  else ! H5T_NATIVE_DOUBLE
    call h5dwrite_f(data_set_id,data_type,real_array,dims, &
                    hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
  endif
  
  ! free everything up
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

end subroutine WriteHDF5DataSet

end module dataset_write_module
