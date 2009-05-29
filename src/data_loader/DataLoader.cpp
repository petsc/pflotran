#include "FileIO.h"
#include "HDF.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  
  int nx, ny, nz;

  cout << "Enter nx ny nz: ";
//  cin >> nx;
//  cin >> ny;
//  cin >> nz;
  nx = 1;
  ny = 2;
  nz = 3;
  int n = nx*ny*nz;

  char filename[1024];
  sprintf(filename,"parameters-%d.h5",nx);

  HDF *file = new HDF(filename,1);
  int compress = 0;

  file->createFileSpace(1,n,NULL,NULL);
  printf("Cell Ids\n");
  file->createDataSet("Cell Ids",HDF_NATIVE_INT,compress);

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

  strcpy(filename,"permeability.final");
  FileIO *datafile = new FileIO(filename);
  for (int i=0; i<n; i++) {
    datafile->getLine();
    datafile->readDouble(&values[i]);
  }
  delete datafile;

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();
 // file->closeDataSpaces();

  // use same data space
  printf("Porosity\n");
  file->createDataSet("Porosity",H5T_NATIVE_DOUBLE,compress);

  strcpy(filename,"Porosity.final");
  datafile = new FileIO(filename);
  for (int i=0; i<n; i++) values[i] = 0.;
  for (int i=0; i<n; i++) {
    datafile->getLine();
    datafile->readDouble(&values[i]);
  }
  delete datafile;

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();
  file->closeDataSpaces();

  delete file;

  printf("Done!\n");
}  
