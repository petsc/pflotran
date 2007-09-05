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

void Grid::createStructured(int nx, int ny, int nz) {
  structuredGrid = new StructuredGrid(nx,ny,nz);
  structuredGrid->createDA();
}

int *Grid::getLocalCellNaturalIDs() {

  int *natural_ids = new int[num_cells_local];
  for (int icell=0; icell<num_cells_local; icell++)
    natural_ids[icell] = -999;
  for (int icell=0; icell<num_cells; icell++) {
    int cell_local_id = cells[icell].getIdLocal();
    if (cell_local_id > -1)
      natural_ids[cell_local_id] = cells[icell].getIdNatural();
  }
  convertLocalCellDataGtoN(natural_ids);
  return natural_ids;

}

int *Grid::getLocalCellMaterialNaturalIDs() {

  int *material_ids = new int[num_cells_local];
  for (int icell=0; icell<num_cells_local; icell++)
    material_ids[icell] = -999;
  for (int icell=0; icell<num_cells; icell++) {
    int cell_local_id = cells[icell].getIdLocal();
    if (cell_local_id > -1)
      material_ids[cell_local_id] = cells[icell].getMaterialId();
  }
  convertLocalCellDataGtoN(material_ids);
  return material_ids;

}


int *Grid::getLocalCellVertexNaturalIDs() {
// return a list of vertex ids for each grid cells 
// currently assuming only hex cells
  if (structuredGrid) return structuredGrid->getLocalCellVertexNaturalIDs(cells,vertices);
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

void Grid::setOrigin(double x, double y, double z) {
  structuredGrid->setOrigin(x,y,z);
}

void Grid::setRotation(double r) {
  structuredGrid->setRotation(r);
}


void Grid::setUpCells() {
  num_cells = num_cells_ghosted;
  cells = new GridCell[num_cells];
  for (int i=0; i<num_cells; i++) {
    cells[i].setIdGhosted(i);
    cells[i].setIdLocal(cell_mapping_ghosted_to_local[i]);
    cells[i].setIdNatural(cell_mapping_ghosted_to_natural[i]);
  }
  if (structuredGrid)
    structuredGrid->setUpCells(num_cells,cells);
}

void Grid::setUpVertices() {
  num_vertices = num_vertices_ghosted;
  vertices = new GridVertex[num_vertices];
  for (int i=0; i<num_vertices; i++) {
    vertices[i].setIdGhosted(i);
    vertices[i].setIdLocal(vertex_mapping_ghosted_to_local[i]);
    vertices[i].setIdNatural(vertex_mapping_ghosted_to_natural[i]);
  }
  if (structuredGrid)
    structuredGrid->setUpVertices(num_vertices,vertices);
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
  for (int i=0; i<num_cells; i++)
    cells[i].printInfo();
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

void Grid::printVertices() {
  PetscErrorCode ierr;
  if (!vertices) return;
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
  if (myrank == 0) printf("\nVertice:\n");
  printf("Processor[%d]\n",myrank);
  for (int i=0; i<num_vertices; i++)
    vertices[i].printInfo();
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

int Grid::getGlobalNumberOfCells() {
  if (structuredGrid) return structuredGrid->getN();
  else return -1;
}

int Grid::getNumberOfCells() {
  return num_cells;
}

int Grid::getNumberOfLocalCells() {
  return num_cells_local;
}

int Grid::getGlobalNumberOfVertices() {
  if (structuredGrid) return (structuredGrid->getNx()+1)*
                             (structuredGrid->getNy()+1)*
                             (structuredGrid->getNz()+1);
  else return -1;
}

int Grid::getNumberOfVertices() {
  return num_vertices;
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
  for (int icell=0; icell<num_cells_local; icell++) {
    int local_id = cells[icell].getIdLocal();
    if (local_id > -1) ptr[local_id] = cells[icell].getMaterialId();
  }
  for (int i=0; i<num_cells_local; i++) {
    if (ptr[i] < -998) 
      printf("ERROR: Grid material ids not set correctly on processor %d\n",myrank);
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
  for (int icell=0; icell<num_cells_local; icell++) {
    int local_id = cells[icell].getIdLocal();
    if (local_id > -1) ptr[local_id] = cells[icell].getActive();
  }
  for (int i=0; i<num_cells_local; i++) {
    if (ptr[i] < -998) 
      printf("ERROR: Grid activity not set correctly on processor %d\n",myrank);
  }
  VecRestoreArray(v,&ptr);
  return v;
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