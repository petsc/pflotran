#include "GridCell.h"

GridCell::GridCell() {
  
  id_local = -999;
  id_ghosted = -999;
  id_natural = -999;
  material_id = 0;
  active = 1;
  volume = -999.;
  porosity = -999.;

  for (int i=0; i<3; i++) {
    centroid[i] = -999.;
    permeability[i] = -999.;
  }
  vertices[0] = 0;
  for (int i=1; i<9; i++)
    vertices[i] = -999;

}

void GridCell::setIdLocal(int i) { id_local = i; }
void GridCell::setIdGhosted(int i) { id_ghosted = i; }
void GridCell::setIdNatural(int i) { id_natural = i; }
void GridCell::setCentroid(double x, double y, double z) {
  centroid[0] = x;
  centroid[1] = y;
  centroid[2] = z;
}
void GridCell::setX(double x) { centroid[0] = x; }
void GridCell::setY(double y) { centroid[1] = y; }
void GridCell::setZ(double z) { centroid[2] = z; }
void GridCell::setVolume(double d) { volume = d; }
void GridCell::setPermX(double d) { permeability[0] = d; }
void GridCell::setPermY(double d) { permeability[1] = d; }
void GridCell::setPermZ(double d) { permeability[2] = d; }
void GridCell::setPorosity(double d) { porosity = d; }
void GridCell::setMaterialId(int i) { material_id = i; }
void GridCell::setActive(int i) { active = i; }
void GridCell::setActive() { active = 1; }
void GridCell::setInactive() { active = 0; }
void GridCell::negateMaterialId() { material_id = -abs(material_id); }

int GridCell::getIdLocal() { return id_local; }
int GridCell::getIdGhosted() { return id_ghosted; }
int GridCell::getIdNatural() { return id_natural; }
double *GridCell::getCentroidPtr() { return &centroid[0]; }
double GridCell::getX() { return centroid[0]; }
double GridCell::getY() { return centroid[1]; }
double GridCell::getZ() { return centroid[2]; }
double GridCell::getVolume() { return volume; }
double GridCell::getPermX() { return permeability[0]; }
double GridCell::getPermY() { return permeability[1]; }
double GridCell::getPermZ() { return permeability[2]; }
int GridCell::getMaterialId() { return material_id; }
int GridCell::getActive() { return active; }

void GridCell::getHexFaceVertices(int face, int *vertex_list) {
  vertex_list[0] = 4;
  if (face == WEST) {
    vertex_list[1] = vertices[1];
    vertex_list[2] = vertices[5];
    vertex_list[3] = vertices[8];
    vertex_list[4] = vertices[4];
  }
  else if (face == EAST) {
    vertex_list[1] = vertices[2];
    vertex_list[2] = vertices[3];
    vertex_list[3] = vertices[7];
    vertex_list[4] = vertices[6];
  }
  else if (face == SOUTH) {
    vertex_list[1] = vertices[1];
    vertex_list[2] = vertices[2];
    vertex_list[3] = vertices[6];
    vertex_list[4] = vertices[5];
  }
  else if (face == NORTH) {
    vertex_list[1] = vertices[3];
    vertex_list[2] = vertices[4];
    vertex_list[3] = vertices[8];
    vertex_list[4] = vertices[7];
  }
  else if (face == BOTTOM) {
    vertex_list[1] = vertices[1];
    vertex_list[2] = vertices[4];
    vertex_list[3] = vertices[3];
    vertex_list[4] = vertices[2];
  }
  else if (face == TOP) {
    vertex_list[1] = vertices[5];
    vertex_list[2] = vertices[6];
    vertex_list[3] = vertices[7];
    vertex_list[4] = vertices[8];
  }
}

void GridCell::printInfo() {
  printf("lid: %d gid: %d nid: %d - vertices[%1d]:",id_local,id_ghosted,
         id_natural,vertices[0]);
  for (int i=1; i<=vertices[0]; i++)
    printf(" %d",vertices[i]);
  printf("\n");
  printf("x: %f  y: %f  z: %f\n",centroid[0],centroid[1],centroid[2]);
  printf("mat: %d  active: %d\n",material_id,active);
  printf("---\n");
}

GridCell::~GridCell() {
}
