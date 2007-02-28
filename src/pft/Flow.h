#ifndef FLOW_H_
#define FLOW_H_

#include "Globals.h"
#include "Grid.h"

#include "include/petsc.h"
#include "include/petscvec.h"
#include "include/petscmat.h"

class Flow {
public:

  Flow(int fdof_, Grid *g);
  virtual ~Flow();
  
private:

  void computeFlow(Grid *g, double dt);
  void computeAccumulation(Grid *g, double dt);
  void computeFlux(Grid *g);
  void computeBoundaryFlux(Grid *g);
  void printMatrix();
  
  PetscErrorCode ierr;
  int fdof; // number of flow degrees of freedom
  Vec residual_vec;
  Mat Jac;
  PetscScalar *pressure_ptr;
  PetscScalar *ppressure_ptr;
  PetscScalar *residual_ptr;
  
};

#endif /*FLOW_H_*/
