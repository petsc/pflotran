#ifdef USE_HDF5

#include "HDF.h"

static int neg_one = -1;

HDF::HDF(char *filename, int overwrite) {

  file_id = -1;
  for (int i=0; i<20; i++)
    grp_id[i] = -1;
  file_space_id = -1;
  memory_space_id = -1;
  data_set_id = -1;
  ngrp = 0;

  for (int i=0; i<3; i++) {
    hyperslab_start[i] = 0;
    hyperslab_stride[i] = 0;
    hyperslab_count[i] = 0;
    hyperslab_block[i] = 0;
  }

  hid_t prop_id = H5Pcreate(H5P_FILE_ACCESS);
#ifndef SERIAL
  H5Pset_fapl_mpio(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL);
#endif

// Create a new file collectively and release property list identifier
  /* Access flags:
     H5F_ACC_TRUNC  - overwrites existing file, deleting original contents
     H5F_ACC_EXCL   - open should fail if file exist
     H5F_ACC_RDONLY - read access only
     H5F_ACC_RDWR   - read and write access
  */

  if (!overwrite) { // solely edit file if it exists
    herr_t (*old_func)(void*);
    void  *old_client_data;
    H5Eget_auto(&old_func, &old_client_data);
    H5Eset_auto(NULL,NULL);
    file_id = H5Fopen(filename,H5F_ACC_RDWR,prop_id);
    H5Eset_auto(old_func,old_client_data);
  }
  if (file_id < 0 || overwrite) {
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, prop_id);
  }
  H5Pclose(prop_id);

  if (file_id < 0) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: HDF file %s failed\n",filename);
  }
}

void HDF::createGroup(char *group_name) {

  herr_t (*old_func)(void*);
  void  *old_client_data;
  H5Eget_auto(&old_func, &old_client_data);
  H5Eset_auto(NULL,NULL);
  if (ngrp > 0)
    grp_id[ngrp] = H5Gopen(grp_id[ngrp-1],group_name);
  else 
    grp_id[ngrp] = H5Gopen(file_id,group_name);
  H5Eset_auto(old_func,old_client_data);
  if (grp_id[ngrp] < 0) {              // < 0 sets to default hint
    if (ngrp > 0)
      grp_id[ngrp] = H5Gcreate(grp_id[ngrp-1],group_name,neg_one);
    else
      grp_id[ngrp] = H5Gcreate(file_id,group_name,neg_one);
  }
  ngrp++;
}

void HDF::closeGroup() {
  if (grp_id[ngrp-1] > -1) {
    H5Gclose(grp_id[ngrp-1]);
    ngrp--;
    grp_id[ngrp] = -1;
  }
}

void HDF::createFileSpace(int rank, int dim0, int dim1, int dim2) {
  HDF::createDataSpace(&file_space_id,rank,dim0,dim1,dim2,dim0,dim1,dim2);
}

void HDF::createMemorySpace(int rank, int dim0, int dim1, int dim2) {
  int product = dim0;
  if (rank > 1) product *= dim1;
  if (rank > 2) product *= dim2;
  if (!product) dim0 = 1;
  HDF::createDataSpace(&memory_space_id,rank,dim0,dim1,dim2,dim0,dim1,dim2);
}

void HDF::createDataSpace(int rank, int dim0, int dim1, int dim2) {
  HDF::createDataSpace(&file_space_id,rank,dim0,dim1,dim2,dim0,dim1,dim2);
}

void HDF::createDataSpace(hid_t *space_id, int rank, int dim0, int dim1, 
                           int dim2) {
  HDF::createDataSpace(space_id,rank,dim0,dim1,dim2,dim0,dim1,dim2);
}

void HDF::createDataSpace(hid_t *space_id, int rank, int dim0, int dim1,
                           int dim2, int max_dim0, int max_dim1, int max_dim2) {

  if (*space_id > -1) H5Sclose(*space_id);
  hsize_t *dims = new hsize_t[rank];
  hsize_t *max_dims = new hsize_t[rank];
  dims[0] = dim0;
  max_dims[0] = max_dim0;
  if (rank > 1) {
    dims[1] = dim1;
    max_dims[1] = max_dim1;
  }
  if (rank > 2) {
    dims[2] = dim2;
    max_dims[2] = max_dim2;
  }
//  if (dims[0] == 0) dims[0] = 1;
//  if (max_dims[0] == 0) max_dims[0] = 1;
  *space_id = H5Screate_simple(rank,dims,max_dims);
  delete [] dims;
  delete [] max_dims;

//printf("DataSpace %d rank: %d  dim[0]: %d",*space_id,rank,(int)dims[0]);
//if (rank > 1) printf(" %d",(int)dims[1]);
//if (rank > 2) printf(" %d",(int)dims[2]);
//printf("\n");
}

void HDF::closeDataSpaces() {
  HDF::closeDataSpace(&file_space_id);
// just do this automatically
  HDF::closeDataSpace(&memory_space_id);
}

void HDF::closeDataSpace(hid_t *space_id) {
  if (*space_id > -1) H5Sclose(*space_id);
  *space_id = -1;
}

void HDF::setHyperSlab(int n) {
  setHyperSlab(n,1);
}

void HDF::setHyperSlab(int n, int stride0) {
  int offset = 0;
  MPI_Exscan(&n,&offset,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
  int start[3] = {0,0,0};
  int stride[3] = {1,1,1};
  int count[3] = {1,1,1};
  start[0] = offset;
  stride[0] = stride0;
  count[0] = n > 0 ? n : 1;
  setHyperSlab(start,stride,count,NULL);
}

void HDF::setHyperSlab(int *start, int *stride, int *count, int *block) {

  if (start) {
    hyperslab_start[0] = (hsize_t)start[0];
    hyperslab_start[1] = (hsize_t)start[1];
    hyperslab_start[2] = (hsize_t)start[2];
  }
  else {
    for (int i=0; i<3; i++)
      hyperslab_start[i] = (hsize_t)0;
  }
  if (stride) {
    hyperslab_stride[0] = stride[0];
    hyperslab_stride[1] = stride[1];
    hyperslab_stride[2] = stride[2];
  }
  else {
    for (int i=0; i<3; i++)
      hyperslab_stride[i] = 1;
  }
  if (count) {
    hyperslab_count[0] = count[0];
    hyperslab_count[1] = count[1];
    hyperslab_count[2] = count[2];
  }
  else {
    for (int i=0; i<3; i++)
      hyperslab_count[i] = 1;
  }
  if (block) {
    hyperslab_block[0] = 1;
    hyperslab_block[1] = 1;
    hyperslab_block[2] = 1;
  }
  else {
    for (int i=0; i<3; i++)
      hyperslab_block[i] = 1;
  }

  status = H5Sselect_hyperslab(file_space_id,H5S_SELECT_SET,
                               hyperslab_start,
                               hyperslab_stride,hyperslab_count,NULL);

//  int myrank;
//  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
//  printf("%d %d %d\n",myrank,(int)hyperslab_start[0],(int)hyperslab_count[0]);
//  printf("%d H5Sselect_hyperslab %d\n", myrank,status);

}

void HDF::createDataSet(char *data_set_name, hid_t type, int compress) {
  // Create data set
  hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);

  // here is where you set chunking, shuffle, compression, deflate, fill_value, etc.
  // chunking
  hsize_t dims[3];
  int ndim = H5Sget_simple_extent_dims(file_space_id,dims,NULL);
//  int size = 65536;
  int size = 32768;
//  int size = 16384;
//  int size = 8192;
//  int size = 4092;
  int min_size = 2048;
  while (size > dims[0])
    size /= 2;
  if (compress && size > min_size) {
    hsize_t *dim = new hsize_t[ndim]; 
    dim[0] = size/ndim;
    if (ndim > 1) dim[1] = ndim;
    if (ndim > 2) dim[2] = 1;
    status = H5Pset_chunk(prop_id,ndim,dim);
    delete [] dim;
    // shuffle
    status = H5Pset_shuffle(prop_id);
  // compression
    status = H5Pset_szip(prop_id,H5_SZIP_NN_OPTION_MASK,16);
    status = H5Pset_deflate(prop_id,9); // 0 - 9
  }

  if (type == H5T_NATIVE_INT) {
    int i = -999;
    H5Pset_fill_value(prop_id,type,&i);
  }
  else {
    double d = -999.;
    H5Pset_fill_value(prop_id,type,&d);
  }

  if (ngrp > 0 && grp_id[ngrp-1] > -1)
    data_set_id = H5Dcreate(grp_id[ngrp-1],data_set_name,type,file_space_id,
                            prop_id);
  else
    data_set_id = H5Dcreate(file_id,data_set_name,type,file_space_id,prop_id);
  H5Pclose(prop_id);
}


void HDF::closeDataSet() {
  if (data_set_id > -1) H5Dclose(data_set_id);
  data_set_id = -1;
}

void HDF::writeInt(int *values) {
  writeInt(values,1);
}

void HDF::writeInt(int *values, int collective) {

  hid_t prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifndef SERIAL
  if (collective)
    H5Pset_dxpl_mpio(prop_id,H5FD_MPIO_COLLECTIVE);
  else
    H5Pset_dxpl_mpio(prop_id,H5FD_MPIO_INDEPENDENT);
#endif
  if (memory_space_id > -1) {

/*
    if (collective) {
      PetscErrorCode ierr;
    }
    int myrank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
    hsize_t dims[3];
    int rank = H5Sget_simple_extent_dims(memory_space_id,dims,NULL);
    ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
    printf("proc: %d - writeInt - ",myrank);
    int count = 1;
    for (int i=0; i< rank; i++)
      count *= (int)dims[i];
    for (int i=0; i<count; i++)
      printf(" %d",values[i]);
    printf("\n");
    printDataSpaceInfo(); 
    if (collective) {
      ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
    }
//*/

    H5Dwrite(data_set_id,H5T_NATIVE_INT,memory_space_id,file_space_id,
             prop_id,values);
  }
  else
    H5Dwrite(data_set_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,prop_id,values);
  H5Pclose(prop_id);

}

void HDF::writeDouble(double *values) {
  writeDouble(values,1);
}

void HDF::writeDouble(double *values, int collective) {

  hid_t prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifndef SERIAL
  if (collective)
    H5Pset_dxpl_mpio(prop_id,H5FD_MPIO_COLLECTIVE);
  else
    H5Pset_dxpl_mpio(prop_id,H5FD_MPIO_INDEPENDENT);
#endif
  if (memory_space_id > -1) {
/*
    printDataSpaceInfo(); 
    PetscErrorCode ierr;
    int myrank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
    hsize_t dims[3];
    int rank = H5Sget_simple_extent_dims(memory_space_id,dims,NULL);
    if (collective) {
      ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
    }
    printf("proc: %d - writeInt - ",myrank);
    int count = 1;
    for (int i=0; i< rank; i++)
      count *= dims[i];
    for (int i=0; i<count; i++)
      printf(" %f",values[i]);
    printf("\n");
    if (collective) {
      ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
    }
//*/
    H5Dwrite(data_set_id,H5T_NATIVE_DOUBLE,memory_space_id,file_space_id,
             prop_id,values);
  }
  else
    H5Dwrite(data_set_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,prop_id,values);
  H5Pclose(prop_id);

}

void HDF::printDataSpaceInfo() {
  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  hsize_t dims[3];
  char string[512];
  char string2[512];
  char string3[512];
  int rank = H5Sget_simple_extent_dims(file_space_id,dims,NULL);
  sprintf(string,"filespace (%d)  - proc: %d  rank: %d  dim0: %d",
          (int)file_space_id,myrank,rank,(int)dims[0]);
  if (rank > 1) {
    sprintf(string2," dim1: %d",(int)dims[1]);
    strcat(string,string2);
  }
  if (rank > 2) {
    sprintf(string2," dim2: %d",(int)dims[2]);
    strcat(string,string2);
  }
  rank = H5Sget_simple_extent_dims(memory_space_id,dims,NULL);
  sprintf(string2,"memoryspace (%d)  - proc: %d  rank: %d  dim0: %d",
          (int)memory_space_id,myrank,rank,(int)dims[0]);
  if (rank > 1) {
    sprintf(string3," dim1: %d",(int)dims[1]);
    strcat(string2,string3);
  }
  if (rank > 2) {
    sprintf(string3," dim2: %d",(int)dims[2]);
    strcat(string2,string3);
  }
  printf("%s\n%s\n",string,string2);
}

void HDF::writeString(char *title, char *string) {
  writeString(title,string,1);
}

void HDF::writeString(char *title, char *string, int collective) {
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_strpad(string_type,H5T_STR_NULLTERM);
  H5Tset_size(string_type,strlen(string)+1);
  createDataSpace(1,1,0,0);
  createDataSet(title,string_type,0);
  hid_t prop_id = H5Pcreate(H5P_DATASET_XFER);
#ifndef SERIAL
  if (collective)
    H5Pset_dxpl_mpio(prop_id,H5FD_MPIO_COLLECTIVE);
  else
    H5Pset_dxpl_mpio(prop_id,H5FD_MPIO_INDEPENDENT);
#endif
  H5Dwrite(data_set_id,string_type,H5S_ALL,H5S_ALL,prop_id,string);
  H5Pclose(prop_id);
  closeDataSet();
  closeDataSpaces();
}

HDF::~HDF() {
  if (file_space_id > -1) 
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: file_space not freed.\n");
  if (memory_space_id > -1) 
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: memory_space not freed.\n");
  if (data_set_id > -1) 
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: data_set not freed.\n");
  if (ngrp > 0)
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: %d groups not freed.\n",ngrp);
  H5Fclose(file_id);
}

#endif
