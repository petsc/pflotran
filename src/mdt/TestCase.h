#ifndef TESTCASE_H_
#define TESTCASE_H_

#include <string.h>

#include "BoundarySet.h"
#include "Connection.h"
#include "Grid.h"

#include "petscsys.h"

class TestCase {
  
public:
  TestCase(Grid **grid);
  virtual ~TestCase();

  void computeTopBoundary(Grid *grid, PetscInt complete);
  void computeBottomBoundary(Grid *grid, PetscInt complete);
  void computeNorthBoundary(Grid *grid, PetscInt complete);
  void computeSouthBoundary(Grid *grid, PetscInt complete);
  void computeEastBoundary(Grid *grid, PetscInt complete);
  void computeWestBoundary(Grid *grid, PetscInt complete);
  void flagGridCells(Grid *grid);

private:

  void setMaterialIdBasedOnNaturalId(PetscInt natural_id, PetscInt material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(PetscInt natural_id, PetscInt active,
                                 Grid *grid);

};

#endif /*TESTCASE_H_*/
