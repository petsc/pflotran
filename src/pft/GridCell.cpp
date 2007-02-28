#include "GridCell.h"

GridCell::GridCell() {
  
  id_local = -1;
  id_ghosted = -1;
  volume = -1;
  viscosity = 1.;
  density = 1.;
  for (int i=0; i<3; i++) {
    centoid[i] = -1;
    permeability[i] = 1.;
  }
  neighbor_cells = NULL;
  
}

void GridCell::setIdLocal(int i) { id_local = i; }
void GridCell::setIdGhosted(int i) { id_ghosted = i; }
void GridCell::setVolume(double d) { volume = d; }
void GridCell::setPermX(double d) { permeability[0] = d; }
void GridCell::setPermY(double d) { permeability[1] = d; }
void GridCell::setPermZ(double d) { permeability[2] = d; }
void GridCell::setViscosity(double d) { viscosity = d; }
void GridCell::setDensity(double d) { density = d; }

int GridCell::getIdLocal() { return id_local; }
int GridCell::getIdGhosted() { return id_ghosted; }
double GridCell::getVolume() { return volume; }
double GridCell::getPermX() { return permeability[0]; }
double GridCell::getPermY() { return permeability[1]; }
double GridCell::getPermZ() { return permeability[2]; }
double GridCell::getViscosity() { return viscosity; }
double GridCell::getDensity() { return density; }

double *GridCell::getPermPtr() { return &permeability[0]; }

void GridCell::printInfo() {
  printf("lid: %d gid: %d vol: %f\n",
         id_local,id_ghosted,volume);
}

void computeAccumulation(double dt) {
}

GridCell::~GridCell() {
  if (neighbor_cells) delete [] neighbor_cells;
}
