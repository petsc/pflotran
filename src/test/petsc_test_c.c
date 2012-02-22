#include <stdio.h>

#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscis.h"
#include "petscviewer.h"

int main(int argc,char* argv[]){

  PetscMPIInt rank;
  PetscMPIInt size;

  PetscInt ndof;
  PetscInt my_ndof;

  Vec vec;
  PetscViewer viewer;

  PetscInitialize(&argc,&argv,(char *)0,(char *)0);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) printf("Beginning of C test program\n");

  ndof = 100000;
  my_ndof = ((int)(ndof / size)) + (ndof % size) > rank ? 1 : 0;

  VecCreateMPI(PETSC_COMM_WORLD, my_ndof, ndof, &vec);
  VecSet(vec, -999.);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, "vec.bin", FILE_MODE_WRITE, 
                        &viewer);
  PetscViewerBinarySetFlowControl(viewer,2);
  if (!rank) printf("Before VecView\n");
  VecView(vec, viewer);
  if (!rank) printf("After VecView\n");
  PetscViewerDestroy(viewer);
  VecDestroy(vec);

  if (!rank) printf("End of C test program\n");

  MPI_Finalize();

}
