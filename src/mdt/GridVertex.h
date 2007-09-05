#ifndef GRIDVERTEX_H_
#define GRIDVERTEX_H_

#include "include/petsc.h"

class GridVertex {
public:
  int num_vertices;

  GridVertex();
  GridVertex(double xcoord, double ycoord, double zcoord);
  virtual ~GridVertex();

  void setIdLocal(int i);
  void setIdGhosted(int i);
  void setIdNatural(int i);
  void setX(double x);
  void setY(double y);
  void setZ(double z);
  
  int getIdLocal();
  int getIdGhosted();
  int getIdNatural();
  double getX();
  double getY();
  double getZ();

  void printInfo();

  int cells[9];

private:
  int id_local, id_ghosted, id_natural;
  double x, y, z;
  
};

#endif /*GRIDVERTEX_H_*/
