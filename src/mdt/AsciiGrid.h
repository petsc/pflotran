#ifndef ASCIIGRID_H_
#define ASCIIGRID_H_

#include <string.h>

#include "include/petsc.h"

#include "FileIO.h"

class AsciiGrid {
  
public:
  AsciiGrid();
  AsciiGrid(char *filename);
  AsciiGrid(char *name_, int ncols_, int nrows_, double xllcorner_, 
            double yllcorner_, double cellsize_, double nodata_, 
            double default_elev_, int material_id_);
  virtual ~AsciiGrid();

  void nullify();
  void setName(char *name_);
  void setMaterialId(int id);
  void getName(char *name_);
  int getMaterialId();
  void readAsciiGridFile(char *filename);
  void computeCoordinates();
  double computeElevationFromCoordinate(double x, double y);
  void printInfo();

  static int nasciigrids;
  
//private:

  char name[32];
  int ncols;
  int nrows;
  int ndata;
  double xllcorner;
  double yllcorner;
  double nodata;
  double cellsize;
  double *values;
  double *xcoord;
  double *ycoord;
  int material_id;

};

#endif /*ASCIIGRID_H_*/
