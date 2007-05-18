#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "include/petsc.h"
#include "include/petscvec.h"
#include "include/petscviewer.h"

#include "Grid.h"

class Output {
public:
  Output();
  void printGMSGrid(Grid *g);
  void printGMSOutput(Grid *g, double time);
  void printGMSPressure(Grid *g, double time);
  void printGMSFlux(Grid *g, double time);
  
  virtual ~Output();

  PetscViewer pressure_handle;
  PetscViewer flux_handle;


};

#endif /*OUTPUT_H_*/
