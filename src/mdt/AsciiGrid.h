#ifndef ASCIIGRID_H_
#define ASCIIGRID_H_

#include <string.h>

#include "petscsys.h"

#include "FileIO.h"

class AsciiGrid {
  
public:
  AsciiGrid();
  AsciiGrid(char *filename);
  AsciiGrid(char *name_, PetscInt ncols_, PetscInt nrows_, PetscReal xllcorner_, 
            PetscReal yllcorner_, PetscReal cellsize_, PetscReal nodata_, 
            PetscReal default_elev_, PetscInt material_id_);
  virtual ~AsciiGrid();

  void nullify();
  void setName(char *name_);
  void setMaterialId(PetscInt id);
  void getName(char *name_);
  PetscInt getMaterialId();
  void readAsciiGridFile(char *filename);
  void computeCoordinates();
  PetscReal computeElevationFromCoordinate(PetscReal x, PetscReal y);
  void printRegion(PetscReal origx, PetscReal origy,
                   PetscReal lenx, PetscReal leny);

  void printInfo();

  static PetscInt nasciigrids;
  
//private:

  char name[32];
  PetscInt ncols;
  PetscInt nrows;
  PetscInt ndata;
  PetscReal xllcorner;
  PetscReal yllcorner;
  PetscReal nodata;
  PetscReal cellsize;
  PetscReal *values;
  PetscReal *xcoord;
  PetscReal *ycoord;
  PetscInt material_id;

};

#endif /*ASCIIGRID_H_*/
