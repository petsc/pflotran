#include "Connection.h"

//Connection::Connection(PetscInt icell, PetscInt iface) {
//
//  cell = icell;
//  face = iface;
//  next = NULL;
//
//}

Connection::Connection(PetscInt icell, PetscInt *vertices_) {

  cell = icell;
  for (PetscInt i=0; i<=vertices_[0]; i++)
    vertices[i] = vertices_[i];
  next = NULL;

}

Connection::Connection(PetscInt icell, PetscInt *vertices_, PetscInt iface) {

  cell = icell;
  face = iface;
  for (PetscInt i=0; i<=vertices_[0]; i++)
    vertices[i] = vertices_[i];
  next = NULL;

}

PetscInt Connection::getFaceVertex(PetscInt ivert) {
  if (ivert <= vertices[0]) return vertices[ivert+1];
  else return -1;
}

PetscInt Connection::getFace() {
  return face;
}

void Connection::printInfo() {
//  PetscPrintf(PETSC_COMM_WORLD,"cell: %d   face %d%",cell,face);
  printf("cell: %d   ",cell);
  if (face > 0) printf("face: %d   ",face);
  printf("vertices:");
  for (PetscInt i=0; i<vertices[0]; i++)
    printf(" %d",vertices[i+1]);
  printf("\n");
}

Connection::~Connection() {
}
