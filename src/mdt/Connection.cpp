#include "Connection.h"

//Connection::Connection(int icell, int iface) {
//
//  cell = icell;
//  face = iface;
//  next = NULL;
//
//}

Connection::Connection(int icell, int *vertices) {

  cell = icell;
  for (int i=0; i<=vertices[0]; i++)
    vertices[i] = vertices[i];
  next = NULL;

}

int Connection::getFaceVertex(int ivert) {
  if (ivert <= vertices[0]) return vertices[ivert+1];
  else return -1;
}

void Connection::printInfo() {
//  PetscPrintf(PETSC_COMM_WORLD,"cell: %d   face %d%",cell,face);
  printf("cell: %d   vertices",cell);
  for (int i=0; i<vertices[0]; i++)
    printf(" %d",vertices[i+1]);
  printf("\n");
}

Connection::~Connection() {
}
