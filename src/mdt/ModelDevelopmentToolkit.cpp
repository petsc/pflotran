#include "petscsys.h"
#include "petscdm.h"

PetscMPIInt myrank = 0;
PetscMPIInt commsize = 0;

#include "Globals.h"
#include "Grid.h"
#include "Hanford300.h"
#include "Hanford300v2.h"
#include "TestCase.h"
#include "IFC.h"
#include "IFC_Scope3.h"
#include "MAD.h"
#include "MidIFC.h"
#include "IFC_2D.h"
#include "IFC_Plot_120x120x20.h"
#include "IFC_Plot_120x120x15.h"
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
//  Hanford300 *hanford300 = NULL;
  Hanford300v2 *hanford300 = NULL;
  IFC *ifc = NULL;
  IFC_Scope3 *ifc3 = NULL;
  IFC_Plot_120x120x20 *ifc_plot_120x120x20 = NULL;
  IFC_Plot_120x120x15 *ifc_plot_120x120x15 = NULL;
  MAD *mad = NULL;
  MidIFC *midifc = NULL;
  TestCase *testcase = NULL;

//  hanford300 = new Hanford300(&grid);
//  hanford300 = new Hanford300v2(&grid);
//  ifc = new IFC(&grid);
//  ifc3 = new IFC_Scope3(&grid);
//  ifc_plot_120x120x20 = new IFC_Plot_120x120x20(&grid);
  ifc_plot_120x120x15 = new IFC_Plot_120x120x15(&grid);
//  mad = new MAD(&grid);
//  midifc = new MidIFC(&grid);
//    ifc_2d = new IFC_2D(&grid);
  //  testcase = new TestCase(&grid);

  if (commsize == 1) {
    if (myrank == 0) printf("\n %d inactive cells out of %d cells total\n\n",
                            grid->getNumInactiveCells(),
                            grid->getNumberOfCellsGlobal());
  }

  Output *out = new Output(grid);

#if 1
  if (grid->getN() <= 50000000) {
    PetscGetTime(&start);
    out->printGMSGrid();
    PetscGetTime(&end);
    if (myrank == 0) printf("  %f seconds to print to GMS\n",end-start); 
  }
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
