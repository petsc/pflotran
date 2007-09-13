#ifndef BOUNDARYSET_H_
#define BOUNDARYSET_H_

#include <stddef.h>
#include <string.h>

#include "include/petsc.h"

#include "Connection.h"
#include "Condition.h"

class BoundarySet {
  
public:
  BoundarySet(char *name_);
  virtual ~BoundarySet();

  void addConnection(Connection *new_connection);
  void convertListToArray();
  int getNumberOfConnectionsLocal();
  int getNumberOfConnectionsGlobal();
  void setCondition(Condition *condition_);
  int *getCellIdsLocal();
  int *getFaceVertexIds(int ivert);
  void printInfo();

  Connection *list;
  Condition *condition;
  BoundarySet *next;
  char name[32];
  
private:

  Connection *end_of_list;
  Connection **_array;
  Connection *connections;
  int num_connections_local;

};

#endif /*BOUNDARYSET_H_*/
