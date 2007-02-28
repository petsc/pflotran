#ifndef GRIDCELL_H_
#define GRIDCELL_H_

#include "include/petsc.h"

class GridCell {
public:
  int num_cells;

  GridCell();
  virtual ~GridCell();
  void setIdLocal(int i);
  void setIdGhosted(int i);
  void setVolume(double d);
  void setPermX(double d);
  void setPermY(double d);
  void setPermZ(double d);
  void setViscosity(double d);
  void setDensity(double d);
  
  int getIdLocal();
  int getIdGhosted();
  double getVolume();
  double getPermX();
  double getPermY();
  double getPermZ();
  double *getPermPtr();
  double getViscosity();
  double getDensity();
  
  void printInfo();
  
private:
  int id_local, id_ghosted;
  double volume;
  double centoid[3];
  GridCell *neighbor_cells;
  
// matrix variables
  double permeability[3];

// fluid variables
  double viscosity;
  double density;

// solute variables

};

#endif /*GRIDCELL_H_*/
