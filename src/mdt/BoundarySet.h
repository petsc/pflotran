#ifndef BOUNDARYSET_H_
#define BOUNDARYSET_H_

#include <stddef.h>
#include <string.h>

//#include "include/petscsys.h"
#include "petscsys.h"

#include "Connection.h"
#include "Condition.h"

class BoundarySet {
  
public:
  BoundarySet(char *name_);
  virtual ~BoundarySet();

  void addConnection(Connection *new_connection);
  void convertListToArray();
  PetscInt getNumberOfConnectionsLocal();
  PetscInt getNumberOfConnectionsGlobal();
  void setCondition(Condition *condition_);
  PetscInt *getCellIdsLocal();
  PetscInt *getCellIdsLocal1Based();
  PetscInt *getFaceIdsLocal();
  PetscInt *getFaceVertexIds(PetscInt ivert);
  void printInfo();

  Connection *list;
  Condition *condition;
  BoundarySet *next;
  char name[32];
  
private:

  Connection *end_of_list;
  Connection **_array;
  Connection *connections;
  PetscInt num_connections_local;

};

#endif /*BOUNDARYSET_H_*/
