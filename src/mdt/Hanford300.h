#ifndef HANFORD300_H_
#define HANFORD300_H_

#include <string.h>

#include "AsciiGrid.h"
#include "BoundarySet.h"
#include "Connection.h"
#include "Grid.h"
#include "Polygon.h"

#include "include/petsc.h"

class Hanford300 {
  
public:
  Hanford300(Grid **grid);
  virtual ~Hanford300();

  void computeRiverBoundary(Grid *grid);
  void computeWestBoundary(Grid *grid);
  void computeSouthBoundary(Grid *grid);
  void computeRechargeBoundary(Grid *grid);

private:

  Polygon *river_polygon;
  AsciiGrid **ascii_grids;

};

#endif /*HANFORD300_H_*/
