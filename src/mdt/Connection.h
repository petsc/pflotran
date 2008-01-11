#ifndef CONNECTION_H_
#define CONNECTION_H_

#include <stddef.h>

#include "include/petsc.h"

class Connection {
  
public:

//  Connection(int icell, int iface);
  Connection(int icell, int *vertices);
  Connection(int icell, int *vertices, int iface);
  virtual ~Connection();

  int getFaceVertex(int ivert);
  int getFace();

  void printInfo();
// these are local ids
  int cell;
  int face;
  int vertices[5];
  Connection *next;

};

#endif /*CONNECTION_H_*/