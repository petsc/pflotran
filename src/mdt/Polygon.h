#ifndef POLYGON_H_
#define POLYGON_H_

#include <string.h>

#include "petscsys.h"

class Polygon {
  
public:
  Polygon();
  virtual ~Polygon();

  void createRiverEdgePolygon();
  void createNorthPondWestTrPolygon();
  void createNorthPondEastTrPolygon();
  void createPlumePolygon();
  void createIFCPolygon();
  void createSPPPolygon();
  void createSPPPolygonCorrect();
  void createMidIFCPolygon();
  PetscInt pointInPolygon(PetscReal x_, PetscReal y_);
  
private:

  PetscInt num_points;
  PetscReal *x;
  PetscReal *y;

};

#endif /*POLYGON_H_*/
