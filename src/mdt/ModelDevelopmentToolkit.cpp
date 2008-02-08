#include "include/petsc.h"
#include "include/petscda.h"

PetscMPIInt myrank = 0;
PetscMPIInt commsize = 0;

#include "Globals.h"
#include "Grid.h"
#include "Hanford300.h"
#include "TestCase.h"
#include "Output.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  
  PetscErrorCode ierr;
  PetscLogDouble start, end;
  
  // initialize petsc
  ierr = PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
//  printf("here 0 - proc: %d\n",myrank);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);
  if (myrank == 0) printf("commsize: %d\n",commsize);

  Grid *grid = NULL;
  Hanford300 *hanford300 = NULL;
  TestCase *testcase = NULL;

  hanford300 = new Hanford300(&grid);
//  testcase = new TestCase(&grid);

  Output *out = new Output(grid);

#if 0
  PetscGetTime(&start);
  out->printGMSGrid();
  PetscGetTime(&end);
  if (myrank == 0) printf("  %f seconds to print to GMS\n",end-start); 
#endif

//  out->printBoundarySets();

#if 0
  PetscGetTime(&start);
  out->printHDFMesh();
  PetscGetTime(&end);
  if (myrank == 0) printf("  %f seconds to print to HDF5 Mesh\n",end-start); 
#endif

#if 1
  PetscGetTime(&start);
  out->printHDFMaterialsAndRegions();
  PetscGetTime(&end);
  if (myrank == 0) printf("  %f seconds to print HDF5 Materials and Regions\n",end-start); 
#endif

#if 0
  PetscGetTime(&start);
  out->printHDFSieveMesh();
  PetscGetTime(&end);
  if (myrank == 0) printf("  %f seconds to print to HDF5 Sieve Mesh\n",end-start); 
#endif

  delete out;
  delete hanford300;
  delete testcase;
  delete grid;

  if (myrank == 0) printf("Done!\n");
  ierr = PetscFinalize();CHKERRQ(ierr);

}  
