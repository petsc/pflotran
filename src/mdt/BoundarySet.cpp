#include "BoundarySet.h"

BoundarySet::BoundarySet(char *name_) {

  PetscPrintf(PETSC_COMM_WORLD,"Creating new Boundary Set (%s).\n",name_);
  strcpy(name,name_);
  list = NULL;
  end_of_list = NULL;
  num_connections = 0;
  next = NULL;
  _array = NULL;
  condition = NULL;

}


void BoundarySet::addConnection(Connection *new_connection) {
  if (!list) list = new_connection;
  if (end_of_list) end_of_list->next = new_connection;
  end_of_list = new_connection;
  num_connections++;
}

void BoundarySet::convertListToArray() {
  int count = 0;
  Connection **_array = new Connection *[num_connections];

  Connection *cur_connection = list;
  while (cur_connection) {
    _array[count++] = cur_connection;
    cur_connection = cur_connection->next;
  }
}

int BoundarySet::getNumberOfConnections() { return num_connections; }

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
