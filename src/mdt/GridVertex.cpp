#include "GridVertex.h"

GridVertex::GridVertex() {
  
  id_local = -999;
  id_ghosted = -999;
  id_natural = -999;

  x = 0.;
  y = 0.;
  z = 0.;

  cells[0] = 0;
  for (int i=1; i<9; i++)
    cells[i] = -999;

}

void GridVertex::setIdLocal(int i) { id_local = i; }
void GridVertex::setIdGhosted(int i) { id_ghosted = i; }
void GridVertex::setIdNatural(int i) { id_natural = i; }
void GridVertex::setX(double x_) { x = x_; }
void GridVertex::setY(double y_) { y = y_; }
void GridVertex::setZ(double z_) { z = z_; }

int GridVertex::getIdLocal() { return id_local; }
int GridVertex::getIdGhosted() { return id_ghosted; }
int GridVertex::getIdNatural() { return id_natural; }
double GridVertex::getX() { return x; }
double GridVertex::getY() { return y; }
double GridVertex::getZ() { return z; }

void GridVertex::printInfo() {
  printf("lid: %d gid: %d nid: %d - cells[%1d]:",id_local,id_ghosted,
         id_natural,cells[0]);
  for (int i=1; i<=cells[0]; i++)
    printf(" %d",cells[i]);
  printf("\n");
  printf("x: %f  y: %f  z: %f\n",x,y,z);
}

GridVertex::~GridVertex() {
}
