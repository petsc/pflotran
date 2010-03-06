#ifndef HDF_H_
#define HDF_H_

#include "petscsys.h"

#if defined(PETSC_HAVE_HDF5)

#include "hdf5.h"

//#ifdef PETSC_USE_64BIT_INDICES
//#define HDF_NATIVE_INT H5T_NATIVE_LLONG 
//#else
#define HDF_NATIVE_INT H5T_NATIVE_INT
//#endif

class HDF {

public:

  HDF(char *filename, PetscInt overwrite);
  virtual ~HDF();

  void createGroup(char *group_name);
  void createFileSpace(PetscInt rank, PetscInt dim0, PetscInt dim1, PetscInt idim2);
  void createMemorySpace(PetscInt rank, PetscInt dim0, PetscInt dim1, PetscInt idim2);
  void createDataSpace(PetscInt rank, PetscInt dim0, PetscInt dim1, PetscInt idim2);
  void createDataSpace(hid_t *space_id, PetscInt rank, PetscInt dim0, PetscInt dim1, 
                       PetscInt idim2);
  static void createDataSpace(hid_t *space_id, PetscInt rank, PetscInt dim0, PetscInt dim1, 
                              PetscInt idim2, PetscInt max_dim0, PetscInt max_dim1, 
                              PetscInt max_dim2);
  void setHyperSlab(PetscInt n);
  void setHyperSlab(PetscInt n, PetscInt stride0);
  void setHyperSlab(PetscInt *start, PetscInt *stride, PetscInt *count, PetscInt *block);
  void createDataSet(char *data_set_name, hid_t type, PetscInt compress);
  void closeGroup();
  void closeDataSpaces();
  void printDataSpaceInfo(); 
  static void closeDataSpace(hid_t *space_id);
  void closeDataSet();
  void writeInt(PetscInt *values);
  void writeInt(PetscInt *values, PetscInt collective);
  void writeDouble(double *values);
  void writeDouble(double *values, PetscInt collective);
  void writeString(char *title, char *string);
  void writeString(char *title, char *string, PetscInt collective);
  void writeAttribute(char *title, char *string);
  void writeAttribute(char *title, PetscInt value);
  void writeAttribute(char *title, double value);

private:

  hid_t file_id; 
  hid_t grp_id[20]; 
  hid_t file_space_id;
  hid_t memory_space_id;
  hid_t data_set_id; 
  herr_t status;

  PetscInt ngrp;

  hsize_t hyperslab_start[3];
  hsize_t hyperslab_stride[3];
  hsize_t hyperslab_count[3];
  hsize_t hyperslab_block[3];
  
};

#endif /*PETSC_HAVE_HDF5*/
#endif /*HDF_H_*/
