#ifndef GRIDCELL_H_ 
#define GRIDCELL_H_ 
 
#include <stdlib.h>

#include "include/petsc.h" 
 
#include "Globals.h" 
 
class GridCell { 
public: 
  int num_cells; 
 
  GridCell(); 
  GridCell(double xcoord, double ycoord, double zcoord, double volume_,
                   double porosity_, int material_id, double xperm,
                   double yperm, double zperm); 
  virtual ~GridCell(); 
  void setIdLocal(int i); 
  void setIdGhosted(int i); 
  void setIdNatural(int i); 
  void setCentroid(double x, double y, double z); 
  void setX(double x); 
  void setY(double y); 
  void setZ(double z); 
  void setVolume(double d); 
  void setPermX(double d); 
  void setPermY(double d); 
  void setPermZ(double d); 
  void setPorosity(double d); 
  void setMaterialId(int i); 
  void setActive(int i); 
  void setActive(); 
  void setInactive(); 
  void negateMaterialId();
   
  int getIdLocal(); 
  int getIdGhosted(); 
  int getIdNatural(); 
  double *getCentroidPtr(); 
  double getX(); 
  double getY(); 
  double getZ(); 
  double getVolume(); 
  double getPermX(); 
  double getPermY(); 
  double getPermZ(); 
  double *getPermPtr(); 
  int getMaterialId(); 
  int getActive(); 
  void getHexFaceVertices(int face, int *vertex_list); 
 
  void printInfo(); 
   
  int vertices[9]; 
  unsigned long int flag;
 
private: 
  int id_local, id_ghosted, id_natural; 
  // materials of inactive cells outside the boundary are set to -material_id
  int material_id; 
  int active; 

  double volume; 
  double centroid[3]; 
 
// matrix variables 
  double permeability[3]; 
  double porosity; 
 
}; 
 
#endif /*GRIDCELL_H_*/ 
