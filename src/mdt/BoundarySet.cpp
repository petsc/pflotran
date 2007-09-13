#include "BoundarySet.h"

BoundarySet::BoundarySet(char *name_) {

  PetscPrintf(PETSC_COMM_WORLD,"Creating new Boundary Set (%s).\n",name_);
  strcpy(name,name_);
  list = NULL;
  end_of_list = NULL;
  num_connections_local = 0;
  next = NULL;
  _array = NULL;
  condition = NULL;

}


void BoundarySet::addConnection(Connection *new_connection) {
  if (!list) list = new_connection;
  if (end_of_list) end_of_list->next = new_connection;
  end_of_list = new_connection;
  num_connections_local++;
}

void BoundarySet::convertListToArray() {
  int count = 0;
  Connection **_array = new Connection *[num_connections_local];

  Connection *cur_connection = list;
  while (cur_connection) {
    _array[count++] = cur_connection;
    cur_connection = cur_connection->next;
  }
}

int BoundarySet::getNumberOfConnectionsLocal() { 
  return num_connections_local; 
}

int BoundarySet::getNumberOfConnectionsGlobal() { 
  int num_connections_global = 0;
  MPI_Allreduce(&num_connections_local,&num_connections_global,1,MPI_INTEGER,
                MPI_SUM,PETSC_COMM_WORLD);
  return num_connections_global;
}

int *BoundarySet::getCellIdsLocal() {
  int count = 0;
  int *cell_ids = new int[num_connections_local];
  Connection *cur_connection = list;
  while (cur_connection) {
    cell_ids[count++] = cur_connection->cell;
    cur_connection = cur_connection->next;
  }
  if (count != num_connections_local) {
    printf("Error: Number of boundary connections incorrect\n");
  }
  return cell_ids;
}

int *BoundarySet::getFaceVertexIds(int ivert) {
  int count = 0;
  int *vertex_ids = new int[num_connections_local];
  Connection *cur_connection = list;
  while (cur_connection) {
    vertex_ids[count++] = cur_connection->getFaceVertex(ivert);
    cur_connection = cur_connection->next;
  }
  if (count != num_connections_local) {
    printf("Error: Number of boundary connections incorrect\n");
  }
  return vertex_ids;
}

void BoundarySet::printInfo() {
  PetscPrintf(PETSC_COMM_WORLD,"BoundarySet::%s\n",name);
  Connection *cur_connection = list;
  while (cur_connection) {
    cur_connection->printInfo();
    cur_connection = cur_connection->next;
  }
}

void BoundarySet::setCondition(Condition *condition_) {
  condition = condition_;
}

BoundarySet::~BoundarySet() {

  if (_array) delete [] _array;

  Connection *cur_connection = list;
  Connection *prev_connection = NULL;
  while (cur_connection) {
    prev_connection = cur_connection;
    cur_connection = cur_connection->next;
    delete prev_connection;
    prev_connection = NULL;
  }
}
