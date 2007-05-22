//#include <iostream>
#include "include/petsc.h"
#include "include/petscts.h"

int myrank = 0;

#include "Globals.h"
#include "Grid.h"
#include "Flow.h"
#include "Output.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  
  PetscErrorCode ierr;
  
  int nx,ny,nz,n,f_dof,t_dof;
  double dx,dy,dz;
  int commsize;

/*
// test case 1
  nx = 3;
  ny = 2;
  nz = 2;
*/
/*
// test case 2
  nx = 60;
  ny = 40;
  nz = 40;
*/
// test case 3
  nx = 1;
  ny = 1;
  nz = 10;

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
//  grid->addBoundaryCondition(1,1,1,ny,1,nz,"west","hydrostatic head",100.);
//  grid->addBoundaryCondition(nx,nx,1,ny,1,nz,"east","hydrostatic head",90.);

  // partially saturated (bottom half of column)
  //grid->addBoundaryCondition(1,nx,1,ny,1,1,"bottom","dirichlet",113566.899); 
  //grid->addBoundaryCondition(1,nx,1,ny,nz,nz,"top","dirichlet",89083.101);

//#define SAT
#ifdef SAT
  // for saturated
  grid->addBoundaryCondition(1,nx,1,ny,1,1,"bottom","dirichlet",150000.); 
  grid->addBoundaryCondition(1,nx,1,ny,nz,nz,"top","dirichlet",125516.202);
  grid->setInitialHydrostaticPressure(150000.,0.);

#else
  // for half-saturated
  grid->addBoundaryCondition(1,nx,1,ny,1,1,"bottom","dirichlet",113566.899); 
  grid->addBoundaryCondition(1,nx,1,ny,nz,nz,"top","dirichlet",89083.101);
  grid->setInitialHydrostaticPressure(113566.899,0.);
#endif

  //  grid->addSource(nx/2,nx/2+1,ny/2,ny/2+1,nz/2,nz/2+1,"mass",90.);
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
  
  Output *output = new Output();

  // Manual time stepping for now
  double final_time = 1000.;
  double time = 0.;
  double dt = 10;
  if (myrank == 0) printf("Time: %f\n",time);
  output->printGMSGrid(grid);
  output->printGMSOutput(grid,time);
  while (time < final_time) {
    if (time+dt > final_time) dt = final_time-time;
    flow->solve(dt);
    time += dt;
    if (int(time/dt+1.e-5)%(int(final_time/dt+1.e-5)/10) == 0) output->printGMSOutput(grid,time);
    if (myrank == 0) printf("Time: %f\n",time);
  }
  
  // clean up
  delete grid;
  delete flow;
  delete output;
  //delete tran;
  //delete react;

  if (myrank == 0) printf("Done!\n");
  ierr = PetscFinalize();CHKERRQ(ierr);

}  
