#include "include/petsc.h"
#include "include/petscda.h"

int myrank = 0;
int commsize = 0;

#include "Globals.h"
#include "Grid.h"
#include "Hanford300.h"
#include "TestCase.h"
#include "Output.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  
  PetscErrorCode ierr;
  
  // initialize petsc
  ierr = PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
//  printf("here 0 - proc: %d\n",myrank);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);
  if (myrank == 0) printf("commsize: %d\n",commsize);

  Grid *grid = NULL;
  Hanford300 *hanford300 = NULL;
  TestCase *testcase = NULL;

//  hanford300 = new Hanford300(&grid);
  testcase = new TestCase(&grid);

  Output *out = new Output(grid);
  out->printGMSGrid();
//  out->printBoundarySets();

  out->printHDFMesh();

  delete out;
  delete hanford300;
  delete testcase;
  delete grid;

  if (myrank == 0) printf("Done!\n");
  ierr = PetscFinalize();CHKERRQ(ierr);

}  
