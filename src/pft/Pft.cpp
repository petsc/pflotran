//#include <iostream>
#include "include/petsc.h"

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
  nx = 6;
  ny = 4;
  nz = 4;

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
  
  PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);
  if (myrank == 0) printf("commsize: %d\n",commsize);

  Grid *grid = new Grid(nx,ny,nz,dx,dy,dz,f_dof,t_dof);
  grid->printConnectivity();
  grid->printCells();
  grid->addBoundaryCondition(1,1,1,ny,1,nz,"west","hydrostatic head",100.);
  grid->addBoundaryCondition(nx,nx,1,ny,1,nz,"east","hydrostatic head",90.);
  grid->addSource(nx/2,nx/2+1,ny/2,ny/2+1,nz/2,nz/2+1,"mass",90.);
  BoundaryCondition::printBCs();
  Source::printSrcs();
  Flow *flow = new Flow(f_dof,grid);

  delete grid;
  delete flow;
  
  delete flow;
  delete grid;

  ierr = PetscFinalize();CHKERRQ(ierr);

}  
