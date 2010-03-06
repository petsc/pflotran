#ifndef GRIDVERTEX_H_
#define GRIDVERTEX_H_

#include "petscsys.h"

class GridVertex {
public:
  PetscInt num_vertices;

  GridVertex();
  GridVertex(PetscReal xcoord, PetscReal ycoord, PetscReal zcoord);
  virtual ~GridVertex();

  void setIdLocal(PetscInt i);
  void setIdGhosted(PetscInt i);
  void setIdNatural(PetscInt i);
  void setX(PetscReal x);
  void setY(PetscReal y);
  void setZ(PetscReal z);
  
  PetscInt getIdLocal();
  PetscInt getIdGhosted();
  PetscInt getIdNatural();
  PetscReal getX();
  PetscReal getY();
  PetscReal getZ();

  void printInfo();

  PetscInt cells[9];

private:
  PetscInt id_local, id_ghosted, id_natural;
  PetscReal x, y, z;
  
};

#endif /*GRIDVERTEX_H_*/
