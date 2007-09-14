#include "Grid.h"


// these static members of the BoundaryConnection class must be initialized HERE!
BoundaryConnection *BoundaryConnection::list = NULL;
BoundaryConnection *BoundaryConnection::end_of_list = NULL;
int BoundaryConnection::num_bcs = 0;
BoundaryConnection **BoundaryConnection::_array = NULL;

Condition *Condition::list = NULL;
Condition *Condition::end_of_list = NULL;
int Condition::num_conditions = 0;
Condition **Condition::_array = NULL;
int Condition::initial_condition_id = -1; 

// these static members of the Source.h class must be initialized HERE!
Source *Source::list = NULL;
Source *Source::end_of_list = NULL;
int Source::num_srcs = 0;

Grid::Grid(int nx, int ny, int nz) {
  nullifyArrays();
  createStructured(nx,ny,nz);
}

Grid::Grid() {
  nullifyArrays();
}

void Grid::zeroGridCellFlags() {
  for (int i=0; i<num_cells_ghosted; i++)
    cells[i].flag = ZERO_FLAG;
}

void Grid::createStructured(int nx, int ny, int nz) {
  structuredGrid = new StructuredGrid(nx,ny,nz);
  structuredGrid->createDA();
}

int *Grid::getCellIds() {

  int *ids = new int[num_cells_local];
  for (int icell=0; icell<num_cells_local; icell++)
    ids[icell] = -999;
  for (int icell=0; icell<num_cells_ghosted; icell++) {
    int cell_local_id = cells[icell].getIdLocal();
    if (cell_local_id > -1)
      ids[cell_local_id] = cell_local_id;
  }
  return ids;

}

int *Grid::getCellIdsNatural() {

  int *natural_ids = new int[num_cells_local];
  for (int icell=0; icell<num_cells_local; icell++)
    natural_ids[icell] = -999;
  for (int icell=0; icell<num_cells_ghosted; icell++) {
    int cell_local_id = cells[icell].getIdLocal();
    if (cell_local_id > -1)
      natural_ids[cell_local_id] = cells[icell].getIdNatural();
  }
  return natural_ids;

}

int *Grid::getCellMaterialIds() {

  int *material_ids = new int[num_cells_local];
  for (int icell=0; icell<num_cells_local; icell++)
    material_ids[icell] = -999;
  for (int icell=0; icell<num_cells_ghosted; icell++) {
    int cell_local_id = cells[icell].getIdLocal();
    if (cell_local_id > -1)
      material_ids[cell_local_id] = cells[icell].getMaterialId();
  }
  return material_ids;

}

int *Grid::getCellVertexIds(int ivert) {

  // the grid vertices array is of size num_vertices_ghosted
  // the cell vertices array is of size 9, where value 0 is the # of vertices in the cell
  // ivert is the index of the vertex in the cell N vertices

  int *vertex_ids = new int[num_cells_local];
  for (int icell=0; icell<num_cells_local; icell++)
    vertex_ids[icell] = -999;
  for (int icell=0; icell<num_cells_ghosted; icell++) {
    int cell_local_id = cells[icell].getIdLocal();
    if (cell_local_id > -1) {
      if (ivert <= cells[icell].vertices[0]) {
        vertex_ids[cell_local_id] = 
              vertices[cells[icell].vertices[ivert]].getIdNatural();
      }
    }
  }
  return vertex_ids;

}

int Grid::getVertexIdsNaturalLocal(int **natural_ids) {

/* need to decide whether we want to print hte vertices truly in natural ordering.  I guess
  natural ordering will always start at zero, thus we can create a vector of size num_vertices_global
  and set the values in the vector.  This vector will automatically be numbered naturally and
  sequentially from proc 0 to proc size-1. */

  Vec vec;
  PetscScalar *vec_ptr = NULL;
  VecCreate(PETSC_COMM_WORLD,&vec);
  VecSetSizes(vec,PETSC_DECIDE,getNumberOfVerticesGlobal());
  VecSetType(vec,VECMPI);

  *natural_ids = new int[num_vertices_local];
  for (int ivert=0; ivert<num_vertices_local; ivert++)
    (*natural_ids)[ivert] = -999;
  for (int ivert=0; ivert<num_vertices_ghosted; ivert++) {
    int vert_local_id = vertices[ivert].getIdLocal();
    if (vert_local_id > -1) {
      (*natural_ids)[vert_local_id] = vertices[ivert].getIdNatural();
    }
  }

  PetscScalar *double_ids = new double[num_vertices_local];
  for (int ivert=0; ivert<num_vertices_local; ivert++)
    double_ids[ivert] = (double)(*natural_ids)[ivert];
  VecSetValues(vec,num_vertices_local,*natural_ids,double_ids,INSERT_VALUES);
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);
  delete [] double_ids;
  double_ids = NULL;
  delete [] *natural_ids;
  *natural_ids = NULL;

  int new_local_size = VecGetLocalSize(vec);
  *natural_ids = new int[new_local_size];
  VecGetArray(vec,&vec_ptr);
  for (int ivert=0; ivert<new_local_size; ivert++)
    (*natural_ids)[ivert] = (int)(vec_ptr[ivert]+0.0001);
  VecRestoreArray(vec,&vec_ptr);
  VecDestroy(vec);

  return new_local_size;

}

int Grid::getVertexCoordinatesNaturalLocal(double **coordinates, int direction) {

  Vec vec;
  PetscScalar *vec_ptr = NULL;
  VecCreate(PETSC_COMM_WORLD,&vec);
  VecSetSizes(vec,PETSC_DECIDE,getNumberOfVerticesGlobal());
  VecSetType(vec,VECMPI);

  int *natural_ids = new int[num_vertices_local];
  *coordinates = new double[num_vertices_local];
  for (int ivert=0; ivert<num_vertices_local; ivert++) {
    (*coordinates)[ivert] = -999.;
    natural_ids[ivert] = -999;
  }
  for (int ivert=0; ivert<num_vertices_ghosted; ivert++) {
    int vert_local_id = vertices[ivert].getIdLocal();
    if (vert_local_id > -1) {
      natural_ids[vert_local_id] = vertices[ivert].getIdNatural();
      if (direction == 0) // x-direction
        (*coordinates)[vert_local_id] = vertices[ivert].getX();
      else if (direction == 1) // y-direction
        (*coordinates)[vert_local_id] = vertices[ivert].getY();
      else // z-direction
        (*coordinates)[vert_local_id] = vertices[ivert].getZ();
    }
  }

  VecSetValues(vec,num_vertices_local,natural_ids,*coordinates,INSERT_VALUES);
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);
  delete [] natural_ids;
  natural_ids = NULL;
  delete [] *coordinates;
  *coordinates = NULL;

  int new_local_size = VecGetLocalSize(vec);
  *coordinates = new double[new_local_size];
  VecGetArray(vec,&vec_ptr);
  for (int ivert=0; ivert<new_local_size; ivert++)
    (*coordinates)[ivert] = vec_ptr[ivert];
  VecRestoreArray(vec,&vec_ptr);
  VecDestroy(vec);
  return new_local_size;

}


void Grid::convertLocalCellDataGtoN(int *data) {

  double *d_data = new double[num_cells_local];
  for (int i=0; i<num_cells_local; i++)
    d_data[i] = (double)data[i];
  convertLocalCellDataGtoN(d_data);
  for (int i=0; i<num_cells_local; i++)
    data[i] = (int)(d_data[i]+0.0001); // avoid roundoff & truncation

}

void Grid::convertLocalCellDataGtoN(double *data) {

  if (structuredGrid) structuredGrid->convertLocalCellDataGtoN(data);

}

void Grid::nullifyArrays() {
  cell_mapping_local_to_ghosted = NULL;
  cell_mapping_ghosted_to_local = NULL;
  cells = NULL;
  vertex_mapping_local_to_ghosted = NULL;
  vertex_mapping_ghosted_to_local = NULL;
  vertices = NULL;
  boundary_sets = NULL;
}

void Grid::setGridSpacing(double *dx, double *dy, double *dz) {
  structuredGrid->setGridSpacing(dx,dy,dz);
}

void Grid::setGridSpacing(double dx, double dy, double dz) {
  structuredGrid->setGridSpacing(dx,dy,dz);
}

void Grid::setLocalGridSpacing() {
  structuredGrid->setLocalGridSpacing();
}

void Grid::setOrigin(double x, double y, double z) {
  structuredGrid->setOrigin(x,y,z);
}

void Grid::setRotation(double r) {
  structuredGrid->setRotation(r);
}


void Grid::setUpCells() {
  cells = new GridCell[num_cells_ghosted];
  for (int i=0; i<num_cells_ghosted; i++) {
    cells[i].setIdGhosted(i);
    cells[i].setIdLocal(cell_mapping_ghosted_to_local[i]);
    cells[i].setIdNatural(cell_mapping_ghosted_to_natural[i]);
  }
  if (structuredGrid)
    structuredGrid->setUpCells(num_cells_ghosted,cells);
}

void Grid::setUpVertices() {
  vertices = new GridVertex[num_vertices_ghosted];
//  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
  for (int i=0; i<num_vertices_ghosted; i++) {
    vertices[i].setIdGhosted(i);
    vertices[i].setIdLocal(vertex_mapping_ghosted_to_local[i]);
    vertices[i].setIdNatural(vertex_mapping_ghosted_to_natural[i]);
//    printf("Proc[%d]: ghst id: %d  loc id: %d  nat id: %d\n", myrank,i,vertex_mapping_ghosted_to_local[i], vertex_mapping_ghosted_to_natural[i]);
  }
//  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
  if (structuredGrid)
    structuredGrid->setUpVertices(num_vertices_ghosted,vertices);
}

void Grid::mapVerticesToCells() {
  if (structuredGrid)
    structuredGrid->mapVerticesToCells(cells,vertices);
}

void Grid::computeCoordinates() {
  if (structuredGrid) structuredGrid->computeCoordinates();
}
#if 0
void Grid::computeConnectivity() {
  if (structuredGrid)
    structuredGrid->computeConnectivity(&num_connections,&connections);
}
#endif
void Grid::computeCellMapping() {
  if (structuredGrid)
    structuredGrid->computeCellMapping(&num_cells_local,&num_cells_ghosted,
                                       &cell_mapping_local_to_ghosted,
                                       &cell_mapping_ghosted_to_local,
                                       &cell_mapping_ghosted_to_natural);
}

void Grid::computeVertexMapping() {
  if (structuredGrid)
    structuredGrid->computeVertexMapping(&num_vertices_local,&num_vertices_ghosted,
                                       &vertex_mapping_local_to_ghosted,
                                       &vertex_mapping_ghosted_to_local,
                                       &vertex_mapping_ghosted_to_natural);
}
#if 0
void Grid::addBoundaryConnection(int is, int ie, int js, int je, int ks, 
                                int ke, char *face, char *type, 
                                double scalar) {

  // set a pointer to the last bc in the list, it exists                                  
  BoundaryConnection *lastbc = BoundaryConnection::end_of_list ?
                              BoundaryConnection::end_of_list : NULL;

  if (structuredGrid)
  // set up boundary connection
    structuredGrid->mapBoundaryConnection(is,ie,js,je,ks,ke,face);

  // set the type and scalar values for the newly created bcs, which
  // are between the lastbc (from above) and the end of the 
  // list.  if lastbc is null, BoundaryConnection::list points to the first.
  BoundaryConnection *cur_bc = lastbc ? lastbc->getNext() : BoundaryConnection::list;

  // add the bc type and scalar; these are independent of the grid type
  while (cur_bc) {
  // set boundary condition type
    cur_bc->setType(type);
  // set boundary condition parameters (e.g. Neumann flux, Dirchlet scalar) 
    cur_bc->setScalar(scalar);
    cur_bc = cur_bc->getNext();
  }
}

BoundaryConnection *Grid::getBoundaryConnections() {
  return BoundaryConnection::list;
}
#endif

void Grid::addBoundarySet(BoundarySet *new_set) {
  if (!boundary_sets) boundary_sets = new_set;
  else {
    BoundarySet *cur_set = boundary_sets;
    while (cur_set->next)
      cur_set = cur_set->next;
    cur_set->next = new_set;
  }
}

BoundarySet *Grid::getBoundarySet(char *name) {
  BoundarySet *cur_set = boundary_sets;
  while (cur_set->next) {
    if (!strcmp(name,cur_set->name)) break;
    cur_set = cur_set->next;
  }
  return cur_set;
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

void Grid::getVectorLocal(Vec *v) {
  if (structuredGrid)
    structuredGrid->getVectorLocal(v);
}

void Grid::getVectorGlobal(Vec *v) {
  if (structuredGrid)
    structuredGrid->getVectorGlobal(v);
}

void Grid::getVectorNatural(Vec *v) {
  if (structuredGrid)
    structuredGrid->getVectorNatural(v);
}

void Grid::globalToNatural(Vec global, Vec natural) {
  if (structuredGrid)
    structuredGrid->globalToNatural(global,natural);
}
#if 0
void Grid::printConnectivity() {
  PetscErrorCode ierr;
  if (!connections) return;
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
  if (myrank == 0) printf("\nConnectivity:\n");
  printf("Processor[%d]\n",myrank);
  for (int i=0; i<num_connections; i++)
    connections[i].printInfo(); 
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}
#endif
void Grid::printCells() {
  PetscErrorCode ierr;
  if (!cells) return;
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
  if (myrank == 0) printf("\nCells:\n");
  printf("Processor[%d]\n",myrank);

//  for (int i=0; i<num_cells_ghosted; i++)
//    printf("%d ",cell_mapping_ghosted_to_natural[i]);
//  printf("\n");

  for (int i=0; i<num_cells_ghosted; i++)
    cells[i].printInfo();
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

void Grid::printVertices() {
  PetscErrorCode ierr;
  if (!vertices) return;
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
  if (myrank == 0) printf("\nVertice:\n");
  printf("Processor[%d]\n",myrank);

//  for (int i=0; i<num_vertices_ghosted; i++)
//    printf("%d ",vertex_mapping_ghosted_to_natural[i]);
//  printf("\n");

  for (int i=0; i<num_vertices_ghosted; i++)
    vertices[i].printInfo();
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

int Grid::getNumberOfCellsGlobal() {
  if (structuredGrid) return structuredGrid->getN();
  else return -1;
}

int Grid::getNumberOfCellsGhosted() {
  return num_cells_ghosted;
}

int Grid::getNumberOfCellsLocal() {
  return num_cells_local;
}

int Grid::getNumberOfVerticesGlobal() {
  if (structuredGrid) return (structuredGrid->getNx()+1)*
                             (structuredGrid->getNy()+1)*
                             (structuredGrid->getNz()+1);
  else return -1;
}

int Grid::getNumberOfVerticesGhosted() {
  return num_vertices_ghosted;
}


int Grid::getNumberOfVerticesLocal() {
  return num_vertices_local;
}

int Grid::getNx() {
  if (structuredGrid) return structuredGrid->getNx();
  else return -999;
}

int Grid::getNy() {
  if (structuredGrid) return structuredGrid->getNy();
  else return -999;
}

int Grid::getNz() {
  if (structuredGrid) return structuredGrid->getNz();
  else return -999;
}

int Grid::getN() {
  if (structuredGrid) return structuredGrid->getN();
  else return -999;
}


void Grid::getCorners(int *xs, int *ys, int *zs,
                      int *nx, int *ny, int *nz) {
  if (structuredGrid) structuredGrid->getCorners(xs,ys,zs,nx,ny,nz);
  else { *xs=-999; *ys=-999; *zs=-999; *nx=-999; *ny=-999; *nz=-999; }
}

void Grid::getGhostCorners(int *xs, int *ys, int *zs, 
                           int *nx, int *ny, int *nz) {
  if (structuredGrid) structuredGrid->getGhostCorners(xs,ys,zs,nx,ny,nz);
  else { *xs=-999; *ys=-999; *zs=-999; *nx=-999; *ny=-999; *nz=-999; }
}

double Grid::getDx(int i) {
  if (structuredGrid) return structuredGrid->getDx(i);
  else return -999.;
}

double Grid::getDy(int j) {
  if (structuredGrid) return structuredGrid->getDy(j);
  else return -999.;
}

double Grid::getDz(int k) {
  if (structuredGrid) return structuredGrid->getDz(k);
  else return -999.;
}

double *Grid::getOriginPtr() {
  if (structuredGrid) return structuredGrid->getOriginPtr();
  else return NULL;
}

double Grid::getRotationDegrees() {
  if (structuredGrid) return structuredGrid->getRotationDegrees();
  else return 0.;
}

Vec Grid::getGridCellMaterialIDs() {
  Vec v;
  PetscScalar *ptr = NULL;
  structuredGrid->getVectorGlobal(&v);
  VecGetArray(v,&ptr);
  for (int i=0; i<num_cells_local; i++)
    ptr[i] = -999;
  for (int icell=0; icell<num_cells_ghosted; icell++) {
    int local_id = cells[icell].getIdLocal();
    if (local_id > -1) ptr[local_id] = cells[icell].getMaterialId();
  }
  for (int i=0; i<num_cells_local; i++) {
    if (ptr[i] < -998) 
      printf("ERROR: Grid material ids (%d,%d) not set correctly on processor %d\n",i,(int)ptr[i],myrank);
  }
  VecRestoreArray(v,&ptr);
  return v;
}

Vec Grid::getGridCellActivities() {
  Vec v;
  PetscScalar *ptr = NULL;
  structuredGrid->getVectorGlobal(&v);
  VecGetArray(v,&ptr);
  for (int i=0; i<num_cells_local; i++)
    ptr[i] = -999;
  for (int icell=0; icell<num_cells_ghosted; icell++) {
    int local_id = cells[icell].getIdLocal();
    if (local_id > -1) ptr[local_id] = cells[icell].getActive();
  }
  for (int i=0; i<num_cells_local; i++) {
    if (ptr[i] < -998) 
      printf("ERROR: Grid activity (%d,%d) not set correctly on processor %d\n",i,(int)ptr[i],myrank);
  }
  VecRestoreArray(v,&ptr);
  return v;
}

void Grid::receiveFlag(int *flag, int direction) {
  structuredGrid->receiveFlag(flag,direction);
}

void Grid::sendFlag(int *flag, int direction) {
  structuredGrid->sendFlag(flag,direction);
}

Grid::~Grid() {
  delete [] cell_mapping_local_to_ghosted;
  delete [] cell_mapping_ghosted_to_local;
//  delete [] connections;
  delete [] cells;

  delete [] vertex_mapping_local_to_ghosted;
  delete [] vertex_mapping_ghosted_to_local;
  delete [] vertices;

  BoundarySet *cur_set = boundary_sets;
  BoundarySet *prev_set = NULL;
  while (cur_set) {
    prev_set = cur_set;
    cur_set = cur_set->next;
    delete prev_set;
  }

#if 0
  BoundaryConnection *cur_bc = BoundaryConnection::list;
  while (cur_bc) {
    BoundaryConnection *next_bc = cur_bc->getNext();
    delete cur_bc;
    cur_bc = next_bc;
  }
#endif
  Source *cur_src = Source::list;
  while (cur_src) {
    Source *next_src = cur_src->getNext();
    delete cur_src;
    cur_src = next_src;
  }

}
