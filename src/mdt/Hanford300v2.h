#ifndef HANFORD300V2_H_
#define HANFORD300V2_H_

#include <string.h>

#include "petscsys.h"

#include "Globals.h"
#include "AsciiGrid.h"
#include "BoundarySet.h"
#include "Connection.h"
#include "Grid.h"
#include "Polygon.h"

class Hanford300v2 {
  
public:
  Hanford300v2(Grid **grid);
  virtual ~Hanford300v2();

  void computeTopBoundary(Grid *grid, PetscInt complete);
  void computeBottomBoundary(Grid *grid, PetscInt complete);
  void computeNorthBoundary(Grid *grid, PetscInt complete);
  void computeSouthBoundary(Grid *grid, PetscInt complete);
  void computeEastBoundary(Grid *grid, PetscInt complete);
  void computeWestBoundary(Grid *grid, PetscInt complete);
  void computeNorthPondWestTrBoundary(Grid *grid, Polygon *p);
  void computeNorthPondEastTrBoundary(Grid *grid, Polygon *p);
  void computePlumeBoundary(Grid *grid, Polygon *p);
  void computePlumeSource(Grid *grid, Polygon *p);
  void computePlumeCells(Grid *grid, Polygon *p);
  void flagGridCells(Grid *grid);

private:

  Polygon *river_polygon;
  Polygon *north_pond_west_trench;
  Polygon *north_pond_east_trench;
  Polygon *plume;
  AsciiGrid **ascii_grids;

  void setMaterialIdBasedOnNaturalId(PetscInt natural_id, PetscInt material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(PetscInt natural_id, PetscInt active,
                                 Grid *grid);

};

#endif /*HANFORD300V2_H_*/
