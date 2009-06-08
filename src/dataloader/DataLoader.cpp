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

  int sample_n = 0;
  int sample_nx, sample_ny, sample_nz;
  sample_nx = nx/2;
  sample_ny = ny/2;
  sample_nz = nz/2;
  double x_scale = ((double)nx)/((double)sample_nx);
  double y_scale = ((double)ny)/((double)sample_ny);
  double z_scale = ((double)nz)/((double)sample_nz);
  sample_n = sample_nx*sample_ny*sample_nz;

  char filename[1024];
  sprintf(filename,"parameters-%d.h5",nx);

  HDF *file = new HDF(filename,1);
  int compress = 0;

  // file #1
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

  // file #1
  printf("Permeability #1\n");
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

  double *sample_perm_values = NULL;
  if (sample_n > 0) {
    sample_perm_values  = new double[sample_n];
    int count = 0;
    for (int k=0; k<sample_nz; k++) {
      for (int j=0; j<sample_ny; j++) {
        for (int i=0; i<sample_nx; i++) {
          int index = i + j*nx + k*nx*ny;
          sample_perm_values[count++] = values[index];
        }
      }
    }
  }

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();

  // use same data space
  // file #1
  printf("Porosity #1\n");
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
  printf("Minimum Porosity: %e (prior to truncation)\n", min);
  printf("Maximum Porosity: %e\n", max);

  double *sample_poro_values = NULL;
  if (sample_n > 0) {
    sample_poro_values = new double[sample_n];
    int count = 0;
    for (int k=0; k<sample_nz; k++) {
      for (int j=0; j<sample_ny; j++) {
        for (int i=0; i<sample_nx; i++) {
          int index = i + j*nx + k*nx*ny;
          sample_poro_values[count++] = values[index];
        }
      }
    }
  }

  file->setHyperSlab(n);
  file->createMemorySpace(1,n,NULL,NULL);
  file->writeDouble(values);

  file->closeDataSet();
  file->closeDataSpaces();
  delete file;

  // file #2
  if (sample_n > 0) {
    char filename2[1024];
    sprintf(filename2,"parameters-%d.h5",sample_nx);
    HDF *file2 = new HDF(filename2,1);
    file2->createFileSpace(1,sample_n,NULL,NULL);
    printf("Cell Ids\n");
    file2->createDataSet("Cell Ids",H5T_NATIVE_INT,compress);

    int *cell_ids = new int[sample_n];
    for (int i=0; i<sample_n; i++)
      cell_ids[i] = i+1;

    file2->setHyperSlab(sample_n);
    file2->createMemorySpace(1,sample_n,NULL,NULL);
    file2->writeInt(cell_ids);

    delete [] cell_ids;
    cell_ids = NULL;
    file2->closeDataSet();

    printf("Permeability #2\n");
    file2->createDataSet("Permeability",H5T_NATIVE_DOUBLE,compress);
    file2->setHyperSlab(sample_n);
    file2->createMemorySpace(1,sample_n,NULL,NULL);
    file2->writeDouble(sample_perm_values);
    file2->closeDataSet();

    printf("Porosity #2\n");
    file2->createDataSet("Porosity",H5T_NATIVE_DOUBLE,compress);
    file2->setHyperSlab(sample_n);
    file2->createMemorySpace(1,sample_n,NULL,NULL);
    file2->writeDouble(sample_poro_values);
    file2->closeDataSet();
    file2->closeDataSpaces();
    delete [] sample_perm_values;
    delete [] sample_poro_values;
    delete file2;
  }

  printf("Done!\n");
}  
