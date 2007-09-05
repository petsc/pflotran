#include "Connection.h"

Connection::Connection(int icell, int iface) {

  cell = icell;
  face = iface;
  next = NULL;

}
 
void Connection::printInfo() {
  PetscPrintf(PETSC_COMM_WORLD,"cell: %d   face %d%",cell,face);
}

Connection::~Connection() {
}
