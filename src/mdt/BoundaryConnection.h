#ifndef BOUNDARYCONNECTION_H_
#define BOUNDARYCONNECTION_H_

#include <string.h>

#include "petscsys.h"

#include "Condition.h"

class BoundaryConnection {
  
public:
  BoundaryConnection();
  virtual ~BoundaryConnection();

  void setId(PetscInt i);
  void setConditionId(PetscInt i);
  PetscInt getConditionType();
  void setDistance(PetscReal *d);
  void setCenter(PetscReal *d);
  void setNormal(PetscReal *d);
  void setArea(PetscReal d);
  void setScalar(PetscReal d);
  void setType(char *str);
  void setNext(BoundaryConnection *bc);
    
  void printInfo();

  static void convertListToArray();
  static void printBCs();

  PetscInt getId();
  PetscInt getConditionId();
  PetscReal *getDistancePtr();
  PetscReal *getCenterPtr();
  PetscReal *getNormalPtr();
  char *getTypePtr();
  PetscReal getArea();
  PetscReal getScalar();
  BoundaryConnection *getNext();
  
  static BoundaryConnection *list;
  static BoundaryConnection *end_of_list;
  static PetscInt num_bcs;
  static BoundaryConnection **_array;
  
private:

  // all boundary ids are local nonghosted numbering 
  PetscInt idlocal;
  PetscInt condition_id;
  PetscReal area;
  PetscReal center[3],dist[3];
  PetscReal normal[3];
  char *type;
  PetscReal scalar;
  BoundaryConnection *next;

};

#endif /*BOUNDARYCONNECTION_H_*/
