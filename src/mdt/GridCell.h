#ifndef GRIDCELL_H_ 
#define GRIDCELL_H_ 
 
#include <stdlib.h>

#include "petscsys.h" 
 
#include "Globals.h" 
 
class GridCell { 
public: 
  PetscInt num_cells; 
 
  GridCell(); 
  GridCell(PetscReal xcoord, PetscReal ycoord, PetscReal zcoord, 
           PetscReal volume_, PetscReal porosity_, PetscInt material_id, 
           PetscReal xperm, PetscReal yperm, PetscReal zperm); 
  virtual ~GridCell(); 
  void setIdLocal(PetscInt i); 
  void setIdGhosted(PetscInt i); 
  void setIdNatural(PetscInt i); 
  void setCentroid(PetscReal x, PetscReal y, PetscReal z); 
  void setCentroidLocal(PetscReal x, PetscReal y, PetscReal z); 
  void setX(PetscReal x); 
  void setY(PetscReal y); 
  void setZ(PetscReal z); 
  void setXLocal(PetscReal x); 
  void setYLocal(PetscReal y); 
  void setZLocal(PetscReal z); 
  void setVolume(PetscReal d); 
  void setPermX(PetscReal d); 
  void setPermY(PetscReal d); 
  void setPermZ(PetscReal d); 
  void setPorosity(PetscReal d); 
  void setMaterialId(PetscInt i); 
  void setActive(PetscInt i); 
  void setActive(); 
  void setInactive(); 
  void negateMaterialId();
   
  PetscInt getIdLocal(); 
  PetscInt getIdGhosted(); 
  PetscInt getIdNatural(); 
  PetscReal *getCentroidPtr(); 
  PetscReal *getCentroidLocalPtr(); 
  PetscReal getX(); 
  PetscReal getY(); 
  PetscReal getZ(); 
  PetscReal getXLocal(); 
  PetscReal getYLocal(); 
  PetscReal getZLocal(); 
  PetscReal getVolume(); 
  PetscReal getPermX(); 
  PetscReal getPermY(); 
  PetscReal getPermZ(); 
  PetscReal *getPermPtr(); 
  PetscInt getMaterialId(); 
  PetscInt getActive(); 
  void getHexFaceVertices(PetscInt face, PetscInt *vertex_list); 
 
  void printInfo(); 
   
  PetscInt vertices[9]; 
  unsigned long int flag;
 
private: 
  PetscInt id_local, id_ghosted, id_natural; 
  // materials of inactive cells outside the boundary are set to -material_id
  PetscInt material_id; 
  PetscInt active; 

  PetscReal volume; 
  PetscReal centroid[3]; 
  PetscReal centroidlocal[3]; 
 
// matrix variables 
  PetscReal permeability[3]; 
  PetscReal porosity; 
 
}; 
 
#endif /*GRIDCELL_H_*/ 
