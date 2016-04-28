#include "GridVertex.h"

GridVertex::GridVertex() {
  
  id_local = -999;
  id_ghosted = -999;
  id_natural = -999;

  x = 0.;
  y = 0.;
  z = 0.;

  cells[0] = 0;
  for (PetscInt i=1; i<9; i++)
    cells[i] = -999;

}

void GridVertex::setIdLocal(PetscInt i) { id_local = i; }
void GridVertex::setIdGhosted(PetscInt i) { id_ghosted = i; }
void GridVertex::setIdNatural(PetscInt i) { id_natural = i; }
void GridVertex::setX(PetscReal x_) { x = x_; }
void GridVertex::setY(PetscReal y_) { y = y_; }
void GridVertex::setZ(PetscReal z_) { z = z_; }

PetscInt GridVertex::getIdLocal() { return id_local; }
PetscInt GridVertex::getIdGhosted() { return id_ghosted; }
PetscInt GridVertex::getIdNatural() { return id_natural; }
PetscReal GridVertex::getX() { return x; }
PetscReal GridVertex::getY() { return y; }
PetscReal GridVertex::getZ() { return z; }

void GridVertex::printInfo() {
  printf("lid: %d gid: %d nid: %d - cells[%1d]:",id_local,id_ghosted,
         id_natural,cells[0]);
  for (PetscInt i=1; i<=cells[0]; i++)
    printf(" %d",cells[i]);
  printf("\n");
  printf("x: %f  y: %f  z: %f\n",x,y,z);
}

GridVertex::~GridVertex() {
}
