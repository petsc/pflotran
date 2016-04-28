#include "AsciiGrid.h"

PetscInt AsciiGrid::nasciigrids = 0;

AsciiGrid::AsciiGrid() {
}

AsciiGrid::AsciiGrid(char *filename) {
  nullify();
  strcpy(name,filename);
  readAsciiGridFile(filename);
}

AsciiGrid::AsciiGrid(char *name_, PetscInt ncols_, PetscInt nrows_, PetscReal xllcorner_, 
                     PetscReal yllcorner_, PetscReal cellsize_, PetscReal nodata_, 
                     PetscReal default_elev_, PetscInt material_id_) {
  strcpy(name,name_);
  ncols = ncols_;
  nrows = nrows_;
  ndata = ncols*nrows;
  xllcorner = xllcorner_;
  yllcorner = yllcorner_;
  cellsize = cellsize_;
  nodata = nodata_;
  values = new PetscReal[ndata];
  for (PetscInt i=0; i<ndata; i++) 
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
void AsciiGrid::setMaterialId(PetscInt id) { material_id = id; }
void AsciiGrid::readAsciiGridFile(char *filename) {

  char word[32];

  FileIO *file = new FileIO(filename);
  for (PetscInt iline=0; iline < 6; iline++) {
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
  values = new PetscReal[ndata];
  for (PetscInt i=0; i<ndata; i++)
    values[i] = nodata;
// file is set up to start reading at upper left hand corner.
  for (PetscInt irow=nrows-1; irow>-1; irow--) {
    file->getLine();
    for (PetscInt icol=0; icol<ncols; icol++) {
      PetscInt id = icol+irow*ncols;
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
PetscInt AsciiGrid::getMaterialId() { return material_id; }

void AsciiGrid::computeCoordinates() {
  PetscReal sum = xllcorner-0.5*cellsize;
  xcoord = new PetscReal[ncols];
  for (PetscInt icol=0; icol<ncols; icol++) {
    sum += cellsize;
    xcoord[icol] = sum;
  }
  
  sum = yllcorner-0.5*cellsize;
  xcoord = new PetscReal[nrows];
  for (PetscInt irow=0; irow<nrows; irow++) {
    sum += cellsize;
    ycoord[irow] = sum;
  }
}

void AsciiGrid::printRegion(PetscReal origx, PetscReal origy, 
                            PetscReal lenx, PetscReal leny) {
  PetscReal half_cellsize = 0.5*cellsize;
  for (PetscInt j=0; j<nrows; j++) {
    PetscReal y = yllcorner + half_cellsize + j*cellsize;
    for (PetscInt i=0; i<ncols; i++) {
      PetscReal x = xllcorner + half_cellsize + i*cellsize;
      if (x - origx >= 0 && x - origx <= lenx &&
          y - origy >= 0 && y - origy <= leny) {
        printf("%.8f %.8f %.8f\n",x,y,values[i+j*ncols]);
      }
    }
  }
}

PetscReal AsciiGrid::computeElevationFromCoordinate(PetscReal x, PetscReal y) {
  PetscReal half_cellsize = 0.5*cellsize;
  PetscInt icol = (PetscInt)((x-xllcorner+half_cellsize)/cellsize)-1;
  PetscInt irow = (PetscInt)((y-yllcorner+half_cellsize)/cellsize)-1;
  if (irow < 0 || icol < 0 || irow >= nrows-1 || icol >= ncols+1) {
    PetscPrintf(PETSC_COMM_WORLD,
                "ERROR:  row or column index outsite ASCII Grid bounds %d %d\n",
                irow,icol);
    return -999.;
  }
  else {
    PetscReal z1 = values[icol+irow*ncols];
    PetscReal z2 = values[icol+1+irow*ncols];
    PetscReal z3 = values[icol+(irow+1)*ncols];
    PetscReal z4 = values[icol+1+(irow+1)*ncols];

    PetscReal x1 = icol*cellsize+half_cellsize+xllcorner;
    PetscReal x2 = x1+cellsize;

    PetscReal y1 = irow*cellsize+half_cellsize+yllcorner;
    PetscReal y2 = y1+cellsize;

#if 0
    // for debugging
    printf("--------------------\n");
    printf("%d %d %.1f %.1f %.8f\n",icol,irow,x1,y1,z1);
    printf("%d %d %.1f %.1f %.8f\n",icol+1,irow,x2,y1,z2);
    printf("%d %d %.1f %.1f %.8f\n",icol,irow+1,x1,y2,z3);
    printf("%d %d %.1f %.1f %.8f\n",icol+1,irow+1,x2,y2,z4);
#endif

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

#if 1
    // bilinear interpolation
    return (z1*(x2-x)*(y2-y)+z2*(x-x1)*(y2-y)+
            z3*(x2-x)*(y-y1)+z4*(x-x1)*(y-y1))/
            ((x2-x1)*(y2-y1));
#else
    // inverse distance weighted interpolation
    PetscReal power = 2.;
    PetscReal temp = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
    if (temp < 1.e-20) return z1;
    PetscReal w1 = 1./pow(temp,power);
    temp = sqrt((x-x2)*(x-x2)+(y-y1)*(y-y1));
    if (temp < 1.e-20) return z2;
    PetscReal w2 = 1./pow(temp,power);
    temp = sqrt((x-x1)*(x-x1)+(y-y2)*(y-y2));
    if (temp < 1.e-20) return z3;
    PetscReal w3 = 1./pow(temp,power);
    temp = sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2));
    if (temp < 1.e-20) return z4;
    PetscReal w4 = 1./pow(temp,power);
    return (w1*z1+w2*z2+w3*z3+w4*z4)/(w1+w2+w3+w4);
#endif
    
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
