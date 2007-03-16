#ifndef FLOW_H_
#define FLOW_H_

#include "Globals.h"
#include "Grid.h"

#include "include/petsc.h"
#include "include/petscvec.h"
#include "include/petscmat.h"
#include "include/petscsnes.h"
#include "include/petscksp.h"

class Flow {
public:

  Flow(int fdof_, Grid *g);
  virtual ~Flow();

  void solve(double dt);
  PetscErrorCode computeFlowR(SNES snes, Vec pressure_vec, Vec f_vec, void *ctx);
  PetscErrorCode computeFlowJ(SNES snes, Vec pressure_vec, Mat *Jac, Mat *B, 
                              MatStructure *m, void *ctx);
  
private:

  void computeAccumulationR(Grid *g, double dt);
  void computeAccumulationJ(Grid *g, double dt);
  void computeFluxR(Grid *g);
  void computeFluxJ(Grid *g);
  void computeBoundaryFluxR(Grid *g);
  void computeBoundaryFluxJ(Grid *g);
  void computeSourceFluxR(Grid *g);
  void printMatrix();
  
  PetscErrorCode ierr;
  int fdof; // number of flow degrees of freedom
  Vec residual_vec;
  Vec ppressure_vec;
  Vec pressure_vec;
  Vec work_vec;
  Mat Jac;
  PetscScalar *pressure_ptr;
  PetscScalar *ppressure_ptr;
  PetscScalar *residual_ptr;
  
  SNES snes;
  KSP ksp;
  PC pc;
  Grid *grid;
  double flow_dt;

};

#endif /*FLOW_H_*/
