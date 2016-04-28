#include <stdio.h>
#include <mpi.h>

int main(int argc,char* argv[]){

  int ierr;
  int rank;
  int commsize;

  ierr = MPI_Init(&argc,&argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&commsize);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  printf("%d of %d processors\n",rank,commsize);
  
  ierr = MPI_Finalize();

}
