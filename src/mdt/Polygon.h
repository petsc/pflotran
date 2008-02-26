#ifndef POLYGON_H_
#define POLYGON_H_

#include <string.h>

#include "include/petsc.h"

class Polygon {
  
public:
  Polygon();
  virtual ~Polygon();

  void createRiverEdgePolygon();
  void createNorthPondWestTrPolygon();
  void createNorthPondEastTrPolygon();
  void createPlumePolygon();
  PetscInt pointInPolygon(PetscReal x_, PetscReal y_);
  
private:

  PetscInt num_points;
  PetscReal *x;
  PetscReal *y;

};

#endif /*POLYGON_H_*/
