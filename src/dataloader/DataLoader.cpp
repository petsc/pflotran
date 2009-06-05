#include "FileIO.h"
#include "HDF.h"

#include <math.h>

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

  strcpy(filename,"permeability.final");
  FileIO *datafile = new FileIO(filename);
  // remove top 3 lines
  datafile->getLine();
  datafile->getLine();
  datafile->getLine();
  for (int i=0; i<n; i++) {
    datafile->getLine();
    double value;
    datafile->readDouble(&value);
    values[i] = exp(value);
  }
  delete datafile;

  double min = 1.e20;
  double max = -1.e20;
  double average = 0.;
  for (int i=0; i<n; i++) {
    if (values[i] < min) min = values[i];
    if (values[i] > max) max = values[i];
    average += values[i];
  }

  average /= double(n);
  double stdev = 0.;
  for (int i=0; i<n; i++) {
    stdev += (values[i]-average)*(values[i]-average);
  }
  stdev = sqrt(stdev/double(n));

  printf("Average Permeability: %e\n", average);
  printf("Stdev Permeability: %e\n", stdev);
  printf("Minimum Permeability: %e\n", min);
  printf("Maximum Permeability: %e\n", max);

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();
 // file->closeDataSpaces();

  // use same data space
  printf("Porosity\n");
  file->createDataSet("Porosity",H5T_NATIVE_DOUBLE,compress);

  strcpy(filename,"porosity.final");
  datafile = new FileIO(filename);
  // remove top 3 lines
  datafile->getLine();
  datafile->getLine();
  datafile->getLine();
  for (int i=0; i<n; i++) values[i] = 0.;
  for (int i=0; i<n; i++) {
    datafile->getLine();
    datafile->readDouble(&values[i]);
  }
  delete datafile;

  min = 1.e20;
  max = -1.e20;
  average = 0.;
  int num_below_0 = 0;
  int num_below_05 = 0;
  double average_trunc = 0.;
  for (int i=0; i<n; i++) {
    if (values[i] < min) min = values[i];
    if (values[i] > max) max = values[i];
    average += values[i];
    if (values[i] < 0.05) {
      if (values[i] < 0.) num_below_0++;
      num_below_05++;
      values[i] = 0.05;
    }
    average_trunc += values[i];
  }

  average /= double(n);
  stdev = 0.;
  for (int i=0; i<n; i++) {
    stdev += (values[i]-average)*(values[i]-average);
  }
  stdev = sqrt(stdev/double(n));

  printf("Average Permeability: %e\n", average);
  printf("Stdev Permeability: %e\n", stdev);
  printf("Number of porosities below 0.: %d\n", num_below_0);
  printf("Number of porosities below 0.05: %d\n", num_below_05);
  printf("Minimum Porosity: %e\n", min);
  printf("Maximum Porosity: %e\n", max);

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();
  file->closeDataSpaces();

  delete file;

  printf("Done!\n");
}  
