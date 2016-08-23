#include <stdio.h>

#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscis.h"
#include "petscviewer.h"

int main(int argc,char* argv[]){

  PetscMPIInt rank;
  PetscMPIInt size;

  Vec x;
  Vec b;
  Mat A;
  KSP ksp;
  PC pc;
  PetscReal value;
  PetscReal tolerance;
  PetscInt i;
  PetscInt n;
  PetscInt offset;
  PetscScalar *x_ptr;

  PetscInitialize(&argc,&argv,(char *)0,(char *)0);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) printf("Beginning of C test program\n");

  n = 2/size;
  MatCreateAIJ(PETSC_COMM_WORLD,n,n,
               PETSC_DETERMINE,PETSC_DETERMINE,
               1,NULL,
               0,NULL,&A);
  VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DETERMINE,&x);
  VecDuplicate(x,&b);
  offset = rank * n;
  value = 1.e-16;
  for(i=0;i<n;i++) {
    MatSetValue(A,i+offset,i+offset,value,INSERT_VALUES);
    VecSetValue(b,i+offset,value,INSERT_VALUES);
  }

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetErrorIfNotConverged(ksp,PETSC_TRUE);
  KSPSetOperators(ksp,A,A);
  KSPGetPC(ksp,&pc);
  KSPSetFromOptions(ksp);
  PCSetFromOptions(pc);
//  KSPSetUp(ksp);

  tolerance = 1.e-20;
  PCFactorSetZeroPivot(pc,tolerance);

//  KSPSetUp(ksp);
  KSPSolve(ksp,b,x);

  VecGetArray(x,&x_ptr);
  if (!rank) printf("These values should be ~1:");
  for(i=0;i<n;i++)
    printf(" %f",x_ptr[i]);
  if (!rank) printf("\n");
  VecRestoreArray(x,&x_ptr);

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&b);

  if (!rank) printf("End of C test program\n");

  MPI_Finalize();

}
