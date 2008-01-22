#ifndef HANFORD300_H_
#define HANFORD300_H_

#include <string.h>

#include "include/petsc.h"

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

  void computeTopBoundary(Grid *grid, int complete);
  void computeBottomBoundary(Grid *grid, int complete);
  void computeNorthBoundary(Grid *grid, int complete);
  void computeSouthBoundary(Grid *grid, int complete);
  void computeEastBoundary(Grid *grid, int complete);
  void computeWestBoundary(Grid *grid, int complete);
  void flagGridCells(Grid *grid);

private:

  Polygon *river_polygon;
  AsciiGrid **ascii_grids;

  void setMaterialIdBasedOnNaturalId(int natural_id, int material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(int natural_id, int active,
                                 Grid *grid);

};

#endif /*HANFORD300_H_*/
