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

  void computeRiverBoundary(Grid *grid);
  void computeWestBoundary(Grid *grid);
  void computeSouthBoundary(Grid *grid);
  void computeRechargeBoundary(Grid *grid);

private:

  void setMaterialIdBasedOnNaturalId(int natural_id, int material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(int natural_id, int active,
                                 Grid *grid);

};

#endif /*TESTCASE_H_*/
