#ifndef FLOW_H_
#define FLOW_H_

#include "Globals.h"
#include "Grid.h"

#include "include/petsc.h"
#include "include/petscvec.h"
#include "include/petscmat.h"
#include "include/petscsnes.h"
#include "include/petscksp.h"
#include "include/petscviewer.h"

class Flow {
public:

  Flow(int fdof_, Grid *g);
  virtual ~Flow();

  void solve(double dt);
  PetscErrorCode computeFlowR(SNES snes, Vec pressure_vec, Vec f_vec, void *ctx);
  PetscErrorCode computeFlowJ(SNES snes, Vec pressure_vec, Mat *Jac, Mat *B, 
                              MatStructure *m, void *ctx);
  
private:

  void init(Grid *g);
  void computeAccumulationR(Grid *g, double dt);
  void computeAccumulationJ(Grid *g);
  void computeFluxR(Grid *g);
  void computeFluxJ(Grid *g);
  void computeBoundaryFluxR(Grid *g);
  void computeBoundaryFluxJ(Grid *g);
  void computeSourceFluxR(Grid *g);
  void computeNumericalJacobian(Grid *g);

  void computeAccumulationRVSat(Grid *g, double dt);
  void computeAccumulationJVSat(Grid *g);
  void computeFluxRVSat(Grid *g);
  void computeFluxJVSat(Grid *g);
  void computeBoundaryFluxRVSat(Grid *g);
  void computeBoundaryFluxJVSat(Grid *g);

  double convertP(double p);

  void printMatrix();
  
  PetscErrorCode ierr;
  int fdof; // number of flow degrees of freedom
  Vec residual_vec;
  Vec work_vec;
  Mat Jac;
  PetscScalar *pressure_ptr;
  PetscScalar *ppressure_ptr;
  PetscScalar *residual_ptr;
  
  SNES snes;
  KSP ksp;
  PC pc;
  PetscViewer viewer;
  Grid *grid;
  double flow_dt;

};

#endif /*FLOW_H_*/
