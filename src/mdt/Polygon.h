#ifndef POLYGON_H_
#define POLYGON_H_

#include <string.h>

#include "include/petsc.h"

class Polygon {
  
public:
  Polygon();
  virtual ~Polygon();

  void createRiverEdgePolygon();
  int pointInPolygon(double x_, double y_);
  
private:

  int num_points;
  double *x;
  double *y;

};

#endif /*POLYGON_H_*/
