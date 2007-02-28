#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include <string.h>

#include "include/petsc.h"

class BoundaryCondition {
  
public:
  BoundaryCondition();
  virtual ~BoundaryCondition();

  void setId(int i);
  void setDistance(double *d);
  void setCenter(double *d);
  void setNormal(double *d);
  void setArea(double d);
  void setScalar(double d);
  void setType(char *str);
  void setNext(BoundaryCondition *bc);
    
  void printInfo();
  static void printBCs();

  int getId();
  double *getDistancePtr();
  double *getCenterPtr();
  double *getNormalPtr();
  char *getTypePtr();
  double getArea();
  double getScalar();
  BoundaryCondition *getNext();

  static BoundaryCondition *list;
  static BoundaryCondition *end_of_list;
  static int num_bcs;
  
private:

  // all boundary ids are local nonghosted numbering 
  int idlocal;
  double area;
  double center[3],dist[3];
  double normal[3];
  char *type;
  double scalar;
  BoundaryCondition *next;

};

#endif /*BOUNDARYCONDITION_H_*/
