#include "GridCell.h"

GridCell::GridCell() {
  
  id_local = -999;
  id_ghosted = -999;
  id_natural = -999;
  material_id = 0;
  active = 1;
  volume = -999.;
  porosity = -999.;

  for (PetscInt i=0; i<3; i++) {
    centroid[i] = -999.;
    centroidlocal[i] = -999.;
    permeability[i] = -999.;
  }
  vertices[0] = 0;
  for (PetscInt i=1; i<9; i++)
    vertices[i] = -999;

}

void GridCell::setIdLocal(PetscInt i) { id_local = i; }
void GridCell::setIdGhosted(PetscInt i) { id_ghosted = i; }
void GridCell::setIdNatural(PetscInt i) { id_natural = i; }
void GridCell::setCentroid(PetscReal x, PetscReal y, PetscReal z) {
  centroid[0] = x;
  centroid[1] = y;
  centroid[2] = z;
}
void GridCell::setCentroidLocal(PetscReal x, PetscReal y, PetscReal z) {
  centroidlocal[0] = x;
  centroidlocal[1] = y;
  centroidlocal[2] = z;
}
void GridCell::setX(PetscReal x) { centroid[0] = x; }
void GridCell::setY(PetscReal y) { centroid[1] = y; }
void GridCell::setZ(PetscReal z) { centroid[2] = z; }
void GridCell::setXLocal(PetscReal x) { centroidlocal[0] = x; }
void GridCell::setYLocal(PetscReal y) { centroidlocal[1] = y; }
void GridCell::setZLocal(PetscReal z) { centroidlocal[2] = z; }
void GridCell::setVolume(PetscReal d) { volume = d; }
void GridCell::setPermX(PetscReal d) { permeability[0] = d; }
void GridCell::setPermY(PetscReal d) { permeability[1] = d; }
void GridCell::setPermZ(PetscReal d) { permeability[2] = d; }
void GridCell::setPorosity(PetscReal d) { porosity = d; }
void GridCell::setMaterialId(PetscInt i) { material_id = i; }
void GridCell::setActive(PetscInt i) { active = i; }
void GridCell::setActive() { active = 1; }
void GridCell::setInactive() { active = 0; }
void GridCell::negateMaterialId() { material_id = -PetscAbsInt(material_id); }

PetscInt GridCell::getIdLocal() { return id_local; }
PetscInt GridCell::getIdGhosted() { return id_ghosted; }
PetscInt GridCell::getIdNatural() { return id_natural; }
PetscReal *GridCell::getCentroidPtr() { return &centroid[0]; }
PetscReal *GridCell::getCentroidLocalPtr() { return &centroidlocal[0]; }
PetscReal GridCell::getX() { return centroid[0]; }
PetscReal GridCell::getY() { return centroid[1]; }
PetscReal GridCell::getZ() { return centroid[2]; }
PetscReal GridCell::getXLocal() { return centroidlocal[0]; }
PetscReal GridCell::getYLocal() { return centroidlocal[1]; }
PetscReal GridCell::getZLocal() { return centroidlocal[2]; }
PetscReal GridCell::getVolume() { return volume; }
PetscReal GridCell::getPermX() { return permeability[0]; }
PetscReal GridCell::getPermY() { return permeability[1]; }
PetscReal GridCell::getPermZ() { return permeability[2]; }
PetscInt GridCell::getMaterialId() { 
  if (active) return material_id; 
  else return 0;
}
PetscInt GridCell::getActive() { return active; }

void GridCell::getHexFaceVertices(PetscInt face, PetscInt *vertex_list) {
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
  for (PetscInt i=1; i<=vertices[0]; i++)
    printf(" %d",vertices[i]);
  printf("\n");
  printf("x: %f  y: %f  z: %f\n",centroid[0],centroid[1],centroid[2]);
  printf("mat: %d  active: %d\n",material_id,active);
  printf("---\n");
}

GridCell::~GridCell() {
}
