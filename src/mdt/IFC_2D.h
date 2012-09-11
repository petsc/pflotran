#ifndef IFC_2D_H_
#define IFC_2D_H_

#include <string.h>

#include "petscsys.h"

#include "Globals.h"
#include "AsciiGrid.h"
#include "BoundarySet.h"
#include "Connection.h"
#include "Grid.h"
#include "Polygon.h"

class IFC_2D {
  
public:
  IFC_2D(Grid **grid);
  virtual ~IFC_2D();

  void computeTopBoundary(Grid *grid, PetscInt complete);
  void computeBottomBoundary(Grid *grid, PetscInt complete);
  void computeEastBoundary(Grid *grid, PetscInt complete);
  void computeWestBoundary(Grid *grid, PetscInt complete);
  void computeIFCRegion(Grid *grid, Polygon *p);
  void flagGridCells(Grid *grid);

private:

  Polygon *ifc_polygon;
  Polygon *river_polygon;
  AsciiGrid **ascii_grids;

  void setMaterialIdBasedOnNaturalId(PetscInt natural_id, PetscInt material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(PetscInt natural_id, PetscInt active,
                                 Grid *grid);

};

#endif /*IFC_2D_H_*/
