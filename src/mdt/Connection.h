#ifndef CONNECTION_H_
#define CONNECTION_H_

#include <stddef.h>

#include "petscsys.h"

class Connection {
  
public:

//  Connection(PetscInt icell, PetscInt iface);
  Connection(PetscInt icell, PetscInt *vertices);
  Connection(PetscInt icell, PetscInt *vertices, PetscInt iface);
  virtual ~Connection();

  PetscInt getFaceVertex(PetscInt ivert);
  PetscInt getFace();

  void printInfo();
// these are local ids
  PetscInt cell;
  PetscInt face;
  PetscInt vertices[5];
  Connection *next;

};

#endif /*CONNECTION_H_*/
