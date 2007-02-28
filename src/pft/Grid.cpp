#include "Grid.h"

// these static members of the BoundaryCondition class must be initialized HERE!
BoundaryCondition *BoundaryCondition::list = NULL;
BoundaryCondition *BoundaryCondition::end_of_list = NULL;
int BoundaryCondition::num_bcs = 0;

// these static members of the Source.h class must be initialized HERE!
Source *Source::list = NULL;
Source *Source::end_of_list = NULL;
int Source::num_srcs = 0;


Grid::Grid() {
  nullifyArrays();
  structuredGrid = NULL;
}

Grid::Grid(int nx, int ny, int nz, int fdof, int tdof) {
  nullifyArrays();
  setUp(nx,ny,nz,fdof,tdof);
} 

Grid::Grid(int nx, int ny, int nz, double dx, double dy, double dz, int fdof, int tdof) {
  nullifyArrays();
  setUp(nx,ny,nz,dx,dy,dz,fdof,tdof);
} 

void Grid::nullifyArrays() {
  mapping_local_to_ghosted = NULL;
  mapping_ghosted_to_local = NULL;
  connections = NULL;
  cells = NULL;
  pressure_vec = PETSC_NULL;
  ppressure_vec = PETSC_NULL;
  density_vec = PETSC_NULL;

}

void Grid::setUp(int nx, int ny, int nz, double dx, double dy, double dz, 
                 int fdof, int tdof) {
  structuredGrid = new StructuredGrid(nx,ny,nz,dx,dy,dz,fdof,tdof);
  setUpConnectivity();
  setUpMapping();
  setUpCells();
  getFdofVectorLocal(&pressure_vec);
  ierr = VecDuplicate(pressure_vec,&ppressure_vec);
  ierr = VecDuplicate(pressure_vec,&density_vec);
}

void Grid::setUp(int nx, int ny, int nz, int fdof, int tdof) {
  setUp(nx,ny,nz,1.,1.,1.,fdof,tdof);
}

void Grid::setUpCells() {
  num_cells = num_nodes_ghosted;
  cells = new GridCell[num_cells];
  for (int i=0; i<num_cells; i++) {
    cells[i].setIdGhosted(i);
    cells[i].setIdLocal(mapping_ghosted_to_local[i]);
  }
  if (structuredGrid)
    structuredGrid->setUpCells(num_cells,cells);
}

void Grid::setUpConnectivity() {
  if (structuredGrid)
    structuredGrid->setUpConnectivity(&num_connections,&connections);
}

void Grid::setUpMapping() {
  if (structuredGrid)
    structuredGrid->setUpMapping(&num_nodes_local,&num_nodes_ghosted,
                                 &mapping_local_to_ghosted,
                                 &mapping_ghosted_to_local);
}

void Grid::addBoundaryCondition(int is, int ie, int js, int je, int ks, 
                                int ke, char *face, char *type, 
                                double scalar) {

  // set a pointer to the last bc in the list, it exists                                  
  BoundaryCondition *lastbc = BoundaryCondition::end_of_list ?
                              BoundaryCondition::end_of_list : NULL;

  if (structuredGrid)
  // set up boundary connection
    structuredGrid->mapBoundaryCondition(is,ie,js,je,ks,ke,face);

  // set the type and scalar values for the newly created bcs, which
  // are between the lastbc (from above) and the end of the 
  // list.  if lastbc is null, BoundaryCondition::list points to the first.
  BoundaryCondition *cur_bc = lastbc ? lastbc : BoundaryCondition::list;

  // add the bc type and scalar; these are independent of the grid type
  while (cur_bc) {
  // set boundary condition type
    cur_bc->setType(type);
  // set boundary condition parameters (e.g. Neumann flux, Dirchlet scalar) 
    cur_bc->setScalar(scalar);
    cur_bc = cur_bc->getNext();
  }
}

BoundaryCondition *Grid::getBoundaryConditions() {
  return BoundaryCondition::list;
}

void Grid::addSource(int is, int ie, int js, int je, int ks, 
                     int ke, char *type, double scalar) {

  // set a pointer to the last bc in the list, it exists                                  
  Source *lastbc = Source::end_of_list ? Source::end_of_list : NULL;

  if (structuredGrid)
  // set up boundary connection
    structuredGrid->mapSource(is,ie,js,je,ks,ke);

  // set the type and scalar values for the newly created bcs, which
  // are between the lastbc (from above) and the end of the 
  // list.  if lastbc is null, Source::list points to the first.
  Source *cur_src = lastbc ? lastbc : Source::list;

  // add the bc type and scalar; these are independent of the grid type
  while (cur_src) {
  // set boundary condition type
    cur_src->setType(type);
  // set boundary condition parameters  
    cur_src->setScalar(scalar);
    cur_src = cur_src->getNext();
  }
}

Source *Grid::getSources() {
  return Source::list;
}

void Grid::get1dofVectorLocal(Vec *v) {
  if (structuredGrid)
    structuredGrid->get1dofVectorLocal(v);
}

void Grid::get1dofVectorGlobal(Vec *v) {
  if (structuredGrid)
    structuredGrid->get1dofVectorGlobal(v);
}

void Grid::getFdofMatrix(Mat *m, MatType mtype) {
  if (structuredGrid)
    structuredGrid->getFdofMatrix(m,mtype);
}

void Grid::getFdofVectorLocal(Vec *v) {
  if (structuredGrid)
    structuredGrid->getFdofVectorLocal(v);
}

void Grid::getFdofVectorGlobal(Vec *v) {
  if (structuredGrid)
    structuredGrid->getFdofVectorGlobal(v);
}

void Grid::getTdofVectorLocal(Vec *v) {
  if (structuredGrid)
    return structuredGrid->getTdofVectorLocal(v);
}

void Grid::getTdofVectorGlobal(Vec *v) {
  if (structuredGrid)
    return structuredGrid->getTdofVectorGlobal(v);
}

#if 0
void Grid::printISLocalToGlobalMapping() {
  
//#include "include/petscis.h"

  PetscErrorCode ierr;
  
  PetscInt num_neighbor_proc, *neighbor_procs, *num_indices, **indices;
  
  ISLocalToGlobalMapping *map = getISLocalToGlobalMapping();
  ierr = ISLocalToGlobalMappingGetInfo(map,
                                       &num_neighbor_proc,
                                       &neighbor_procs,
                                       &num_indices,
                                       &indices);

  if (rank == 0) printf("/**************************************************/\n");
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,rank);
  printf("proc%d\n",rank);
  for (int i=0; i<num_neighbor_proc; i++) {
    printf("neighbor_proc: %d - %d nodes\n",neighbor_procs[i],num_indices[i]);
    for (int j=0; j<num_indices[i]; j++) {
      printf("node %d\n",indices[i][j]);
    }
  }
  printf("\n");
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,rank);
                                         
  ierr = ISLocalToGlobalMappingRestoreInfo(map,&num_neighbor_proc,&neighbor_procs,
                                           &num_indices,&indices);
}

void Grid::printAO() {

#include "include/petscao.h"

  PetscErrorCode ierr;
  AO ao;
  ierr = DAGetAO(da,&ao);
  AOView(ao,PETSC_VIEWER_STDOUT_WORLD);
  AODestroy(ao);

}
#endif

void Grid::printConnectivity() {
  PetscErrorCode ierr;
  if (!connections) return;
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD);
  if (myrank == 0) printf("\nConnectivity:\n");
  printf("Processor[%d]\n",myrank);
  for (int i=0; i<num_connections; i++)
    connections[i].printInfo(); 
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD);
}

void Grid::printCells() {
  PetscErrorCode ierr;
  if (!cells) return;
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD);
  if (myrank == 0) printf("\nCells:\n");
  printf("Processor[%d]\n",myrank);
  for (int i=0; i<num_cells; i++)
    cells[i].printInfo();
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD);
}


void Grid::test(int **array) {
  printf("test1\n");
  array = new int*[1];
  array[0] = new int[1];
  array[0][0] = 5;
}
void Grid::test2(int ***array) {
  printf("test2\n");
  *array = new int*[1];
  *array[0] = new int[1];
  *array[0][0] = 5;
}

Grid::~Grid() {
  delete [] mapping_local_to_ghosted;
  delete [] mapping_ghosted_to_local;
  delete [] connections;
  delete [] cells;
  delete structuredGrid;
  VecDestroy(ppressure_vec);
  VecDestroy(pressure_vec);
  VecDestroy(density_vec);

  BoundaryCondition *cur_bc = BoundaryCondition::list;
  while (cur_bc) {
    BoundaryCondition *next_bc = cur_bc->getNext();
    delete cur_bc;
    cur_bc = next_bc;
  }

  Source *cur_src = Source::list;
  while (cur_src) {
    Source *next_src = cur_src->getNext();
    delete cur_src;
    cur_src = next_src;
  }

}