#ifndef HDF_H_
#define HDF_H_

#include "hdf5.h"
//#include <stdio.h>
//#include <string.h>

//#ifdef PETSC_USE_64BIT_INDICES
//#define HDF_NATIVE_INT H5T_NATIVE_LLONG 
//#else
#define HDF_NATIVE_INT H5T_NATIVE_LONG
//#endif

class HDF {

public:

  HDF(char *filename, int overwrite);
  virtual ~HDF();

  void createGroup(char *group_name);
  void createFileSpace(int rank, int dim0, int dim1, int idim2);
  void createMemorySpace(int rank, int dim0, int dim1, int idim2);
  void createDataSpace(int rank, int dim0, int dim1, int idim2);
  void createDataSpace(hid_t *space_id, int rank, int dim0, int dim1, 
                       int idim2);
  static void createDataSpace(hid_t *space_id, int rank, int dim0, int dim1, 
                              int idim2, int max_dim0, int max_dim1, 
                              int max_dim2);
  void setHyperSlab(int n);
  void setHyperSlab(int n, int stride0);
  void setHyperSlab(int *start, int *stride, int *count, int *block);
  void createDataSet(char *data_set_name, hid_t type, int compress);
  void closeGroup();
  void closeDataSpaces();
  void printDataSpaceInfo(); 
  static void closeDataSpace(hid_t *space_id);
  void closeDataSet();
  void writeInt(int *values);
  void writeInt(int *values, int collective);
  void writeDouble(double *values);
  void writeDouble(double *values, int collective);
  void writeString(char *title, char *string);
  void writeString(char *title, char *string, int collective);
  void writeAttribute(char *title, char *string);
  void writeAttribute(char *title, int value);
  void writeAttribute(char *title, double value);

private:

  hid_t file_id; 
  hid_t grp_id[20]; 
  hid_t file_space_id;
  hid_t memory_space_id;
  hid_t data_set_id; 
  herr_t status;

  int ngrp;

  hsize_t hyperslab_start[3];
  hsize_t hyperslab_stride[3];
  hsize_t hyperslab_count[3];
  hsize_t hyperslab_block[3];
  
};

#endif /*HDF_H_*/
