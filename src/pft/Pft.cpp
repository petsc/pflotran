//#include <iostream>
#include "include/petsc.h"
#include "include/petscts.h"

int myrank = 0;

#include "Globals.h"
#include "Grid.h"
#include "Flow.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  
  PetscErrorCode ierr;
  
  int nx,ny,nz,n,f_dof,t_dof;
  double dx,dy,dz;
  int commsize;

// test case 1
  nx = 3;
  ny = 2;
  nz = 2;

/*
// test case 2
  nx = 60;
  ny = 40;
  nz = 40;
*/

  dx = 1.;
  dy = 1.;
  dz = 0.25;

  n = nx*ny*nz;
  f_dof = 1;
  t_dof = 1;
  
  // initialize petsc
  ierr = PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
//  printf("here 0 - proc: %d\n",myrank);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);
  if (myrank == 0) printf("commsize: %d\n",commsize);

  // initialize problem domain
//  printf("here 1 - proc: %d\n",myrank);
  Grid *grid = new Grid(nx,ny,nz,dx,dy,dz,f_dof,t_dof);
//  grid->printConnectivity();
//  printf("here 2 - proc: %d\n",myrank);
//  grid->printCells();

  // initialize boundary conditions and sources
  grid->addBoundaryCondition(1,1,1,ny,1,nz,"west","hydrostatic head",100.);
  grid->addBoundaryCondition(nx,nx,1,ny,1,nz,"east","hydrostatic head",90.);
  grid->addSource(nx/2,nx/2+1,ny/2,ny/2+1,nz/2,nz/2+1,"mass",90.);
//  BoundaryCondition::printBCs();
//  Source::printSrcs();
  
  // initialize flow
  Flow *flow = new Flow(f_dof,grid);
  // initialize transport
  //Transport *tran = new Transport(t_dof,grid); // needs to be implemented
  // initialize reaction
  //Reaction *react = new Reaction(); // needs to be implemented
  
  // initialized time step
//  TS ts_flow;
//  ierr = TSCreate(PETSC_COMM_WORLD,&ts_flow);
//  ierr = TSSetProblemType(ts,TS_NONLINEAR);
//  ierr = TSSetRHSFunction(ts,flow->computeRHS(),&grid)
  
  // Manual time stepping for now
  double final_time = 1.;
  double time = 0.;
  double dt = 0.01;
  if (myrank == 0) printf("Time: %f\n",time);
  while (time < final_time) {
    flow->solve(dt);
    if (time+dt > final_time) dt = final_time-time;
    else time += dt;
    if (myrank == 0) printf("Time: %f\n",time);
  }
  
  // clean up
  delete grid;
  delete flow;
  //delete tran;
  //delete react;

  if (myrank == 0) printf("Done!\n");
  ierr = PetscFinalize();CHKERRQ(ierr);

}  
