#include "FileIO.h"
#include "HDF.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  
  int nx, ny, nz;

  cout << "Enter nx ny nz: ";
  cin >> nx;
  cin >> ny;
  cin >> nz;
//  nx = 1;
//  ny = 2;
//  nz = 3;

  int n = nx*ny*nz;

  int sample_nx, sample_ny, sample_nz;
  sample_nx = 32;
  sample_ny = 32;
  sample_nz = 32;
  double x_scale = ((double)nx)/((double)sample_nx);
  double y_scale = ((double)ny)/((double)sample_ny);
  double z_scale = ((double)nz)/((double)sample_nz);
  int sample_n = sample_nx*sample_ny*sample_nz;

  char filename[1024];
  sprintf(filename,"parameters-%d.h5",nx);

  HDF *file = new HDF(filename,1);
  int compress = 0;

  file->createFileSpace(1,n,NULL,NULL);
  printf("Cell Ids\n");
  file->createDataSet("Cell Ids",H5T_NATIVE_INT,compress);

  int *cell_ids = new int[n];
  for (int i=0; i<n; i++)
    cell_ids[i] = i+1;

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeInt(cell_ids);

  delete [] cell_ids;
  cell_ids = NULL;

  file->closeDataSet();
 // file->closeDataSpaces();

  printf("Permeability\n");
  file->createDataSet("Permeability",H5T_NATIVE_DOUBLE,compress);

  double *values = new double[n];
  for (int i=0; i<n; i++) values[i] = 0.;

  double *sample_values = new double[sample_n];
  for (int i=0; i<sample_n; i++) sample_values[i] = 0.;

//  strcpy(filename,"permeability.final");
  strcpy(filename,"permeability.sample");
  FileIO *datafile = new FileIO(filename);
//  for (int i=0; i<n; i++) {
  for (int i=0; i<sample_n; i++) {
    datafile->getLine();
//    datafile->readDouble(&values[i]);
    datafile->readDouble(&sample_values[i]);
  }
  delete datafile;

  int count = 0;
  for (int k=0; k<nz; k++) {
    int kk = (int)(((double)k)/z_scale);
    for (int j=0; j<ny; j++) {
      int jj = (int)(((double)j)/y_scale);
      for (int i=0; i<nx; i++) {
        int ii = (int)(((double)i)/x_scale);
        int index = ii + jj*sample_nx + kk*sample_nx*sample_ny;
        values[count++] = sample_values[index];
      }
    }
  }

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();
 // file->closeDataSpaces();

  // use same data space
  printf("Porosity\n");
  file->createDataSet("Porosity",H5T_NATIVE_DOUBLE,compress);

  for (int i=0; i<sample_n; i++) sample_values[i] = 0.;

//  strcpy(filename,"porosity.final");
  strcpy(filename,"porosity.sample");
  datafile = new FileIO(filename);
  for (int i=0; i<n; i++) values[i] = 0.;
//  for (int i=0; i<n; i++) {
  for (int i=0; i<sample_n; i++) {
    datafile->getLine();
//    datafile->readDouble(&values[i]);
    datafile->readDouble(&sample_values[i]);
  }
  delete datafile;

  count = 0;
  for (int k=0; k<nz; k++) {
    int kk = (int)(((double)k)/z_scale);
    for (int j=0; j<ny; j++) {
      int jj = (int)(((double)j)/y_scale);
      for (int i=0; i<nx; i++) {
        int ii = (int)(((double)i)/x_scale);
        int index = ii + jj*sample_nx + kk*sample_nx*sample_ny;
        values[count++] = sample_values[index];
      }
    }
  }

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();
  file->closeDataSpaces();

  delete file;

  delete [] values;
  delete [] sample_values;
 
  printf("Done!\n");
}  
