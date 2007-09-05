#include "AsciiGrid.h"

int AsciiGrid::nasciigrids = 0;

AsciiGrid::AsciiGrid() {
}

AsciiGrid::AsciiGrid(char *filename) {
  nullify();
  strcpy(name,filename);
  readAsciiGridFile(filename);
}

AsciiGrid::AsciiGrid(char *name_, int ncols_, int nrows_, double xllcorner_, 
                     double yllcorner_, double cellsize_, double nodata_, 
                     double default_elev_, int material_id_) {
  strcpy(name,name_);
  ncols = ncols_;
  nrows = nrows_;
  ndata = ncols*nrows;
  xllcorner = xllcorner_;
  yllcorner = yllcorner_;
  cellsize = cellsize_;
  nodata = nodata_;
  values = new double[ndata];
  for (int i=0; i<ndata; i++) 
    values[i] = default_elev_;
  material_id = material_id_;
}


void AsciiGrid::nullify() {
  name[0] = '\0';
  ncols=0; nrows=0;
  xllcorner=0.;yllcorner=0.;
  nodata = -9999.;
  cellsize = 0.;
  values = NULL;
  xcoord = NULL;
  ycoord = NULL;
  material_id = 0;
}

void AsciiGrid::setName(char *name_) { strcpy(name,name_); }
void AsciiGrid::setMaterialId(int id) { material_id = id; }
void AsciiGrid::readAsciiGridFile(char *filename) {

  char word[32];

  FileIO *file = new FileIO(filename);
  for (int iline=0; iline < 6; iline++) {
    file->getLine();
    file->readWord(word);
  //  cout << word << "\n";
    if (!strcmp(word,"ncols"))
      file->readInt(&ncols);
    else if (!strcmp(word,"nrows"))
      file->readInt(&nrows);
    else if (!strcmp(word,"xllcorner"))
      file->readDouble(&xllcorner);
    else if (!strcmp(word,"yllcorner"))
      file->readDouble(&yllcorner);
    else if (!strcmp(word,"cellsize"))
      file->readDouble(&cellsize);
    else if (!strcmp(word,"NODATA_value"))
      file->readDouble(&nodata);
  }

  ndata = ncols*nrows;
  values = new double[ndata];
  for (int i=0; i<ndata; i++)
    values[i] = nodata;
// file is set up to start reading at upper left hand corner.
  for (int irow=nrows-1; irow>-1; irow--) {
    file->getLine();
    for (int icol=0; icol<ncols; icol++) {
      int id = icol+irow*ncols;
      if (id > ndata) {
        PetscPrintf(PETSC_COMM_WORLD,"ERROR: Number of ASCII Grid data read higher than number in file.\n");
        PetscFinalize();
        exit(0);
      }
      file->readDouble(&(values[id]));
    }
  }
}

void AsciiGrid::getName(char *name_) { strcpy(name_,name); }
int AsciiGrid::getMaterialId() { return material_id; }

void AsciiGrid::computeCoordinates() {
  double sum = xllcorner-0.5*cellsize;
  xcoord = new double[ncols];
  for (int icol=0; icol<ncols; icol++) {
    sum += cellsize;
    xcoord[icol] = sum;
  }
  
  sum = yllcorner-0.5*cellsize;
  xcoord = new double[nrows];
  for (int irow=0; irow<nrows; irow++) {
    sum += cellsize;
    ycoord[irow] = sum;
  }
}

double AsciiGrid::computeElevationFromCoordinate(double x, double y) {
  double half_cellsize = 0.5*cellsize;
  int icol = (int)((x-xllcorner+half_cellsize)/cellsize)-1;
  int irow = (int)((y-yllcorner+half_cellsize)/cellsize)-1;
  if (irow < 0 || icol < 0 || irow >= nrows-1 || icol >= ncols+1) {
    PetscPrintf(PETSC_COMM_WORLD,
                "ERROR:  row or column index outsite ASCII Grid bounds %d %d\n",
                irow,icol);
  }
  else {
    double z1 = values[icol+irow*ncols];
    double z2 = values[icol+1+irow*ncols];
    double z3 = values[icol+(irow+1)*ncols];
    double z4 = values[icol+1+(irow+1)*ncols];

    double x1 = icol*cellsize+half_cellsize+xllcorner;
    double x2 = x1+cellsize;

    double y1 = irow*cellsize+half_cellsize+yllcorner;
    double y2 = y1+cellsize;

    if (x < x1 || x > x2) {
      PetscPrintf(PETSC_COMM_WORLD,
                  "ERROR:  x out of bounds for interpolation %f %f %f\n",
                  x,x1,x2);
    }

    if (y < y1 || y > y2) {
      PetscPrintf(PETSC_COMM_WORLD,
                  "ERROR:  y out of bounds for interpolation %f %f %f\n",
                  y,y1,y2);
    }

    // bilinear interpolation
    return (z1*(x2-x)*(y2-y)+z2*(x-x1)*(y2-y)+
            z3*(x-x1)*(y-y1)+z4*(x2-x)*(y-y1))/
            ((x2-x1)*(y2-y1));
  }
}

void AsciiGrid::printInfo() {
  PetscPrintf(PETSC_COMM_WORLD,"name:         %s\n",name);
  PetscPrintf(PETSC_COMM_WORLD,"nrows:        %d\n",nrows);
  PetscPrintf(PETSC_COMM_WORLD,"ncols:        %d\n",ncols);
  PetscPrintf(PETSC_COMM_WORLD,"xllcorner:    %f\n",xllcorner);
  PetscPrintf(PETSC_COMM_WORLD,"yllcorner:    %f\n",yllcorner);
  PetscPrintf(PETSC_COMM_WORLD,"cellsize:     %f\n",cellsize);
  PetscPrintf(PETSC_COMM_WORLD,"NODATA_value: %f\n",nodata);
  if (values) {
    PetscPrintf(PETSC_COMM_WORLD,"value0:       %f\n",values[0]);
    PetscPrintf(PETSC_COMM_WORLD,"valueN:       %f\n",values[ndata-1]);
  }
  PetscPrintf(PETSC_COMM_WORLD,"material id:  %d\n",material_id);
}

AsciiGrid::~AsciiGrid() {
}
