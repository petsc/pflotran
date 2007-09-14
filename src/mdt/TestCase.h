#ifndef TESTCASE_H_
#define TESTCASE_H_

#include <string.h>

#include "BoundarySet.h"
#include "Connection.h"
#include "Grid.h"

#include "include/petsc.h"

class TestCase {
  
public:
  TestCase(Grid **grid);
  virtual ~TestCase();

  void computeTopBoundary(Grid *grid, int complete);
  void computeBottomBoundary(Grid *grid, int complete);
  void computeNorthBoundary(Grid *grid, int complete);
  void computeSouthBoundary(Grid *grid, int complete);
  void computeEastBoundary(Grid *grid, int complete);
  void computeWestBoundary(Grid *grid, int complete);
  void flagGridCells(Grid *grid);

private:

  void setMaterialIdBasedOnNaturalId(int natural_id, int material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(int natural_id, int active,
                                 Grid *grid);

};

#endif /*TESTCASE_H_*/
