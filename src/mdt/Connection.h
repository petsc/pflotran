#ifndef CONNECTION_H_
#define CONNECTION_H_

#include <stddef.h>

#include "include/petsc.h"

class Connection {
  
public:

  Connection(int icell, int iface);
  virtual ~Connection();

  void printInfo();

  int cell;
  int face;
  Connection *next;

};

#endif /*CONNECTION_H_*/