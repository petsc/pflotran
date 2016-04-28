#ifndef IFC_Plot_120x120x15_H_
#define IFC_Plot_120x120x15_H_

#include <string.h>

#include "petscsys.h"

#include "Globals.h"
#include "AsciiGrid.h"
#include "BoundarySet.h"
#include "Connection.h"
#include "Grid.h"
#include "Polygon.h"

class IFC_Plot_120x120x15 {
  
public:
  IFC_Plot_120x120x15(Grid **grid);
  virtual ~IFC_Plot_120x120x15();

  void computeTopBoundary(Grid *grid, PetscInt complete);
  void computeBottomBoundary(Grid *grid, PetscInt complete);
  void computeNorthBoundary(Grid *grid, PetscInt complete);
  void computeSouthBoundary(Grid *grid, PetscInt complete);
  void computeEastBoundary(Grid *grid, PetscInt complete);
  void computeWestBoundary(Grid *grid, PetscInt complete);
  void computeIFCBoundary(Grid *grid, Polygon *p);
  void flagGridCells(Grid *grid);

private:

  Polygon *ifc_polygon;
  Polygon *boundary_polygon;
  AsciiGrid **ascii_grids;

  void setMaterialIdBasedOnNaturalId(PetscInt natural_id, PetscInt material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(PetscInt natural_id, PetscInt active,
                                 Grid *grid);

};

#endif /*IFC_Plot_120x120x15_H_*/
