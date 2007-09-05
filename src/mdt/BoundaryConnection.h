#ifndef BOUNDARYCONNECTION_H_
#define BOUNDARYCONNECTION_H_

#include <string.h>

#include "include/petsc.h"

#include "Condition.h"
#include "GridCell.h"

class BoundaryConnection {
  
public:
  BoundaryConnection();
  virtual ~BoundaryConnection();

  void setId(int i);
  void setConditionId(int i);
  int getConditionType();
  void setDistance(double *d);
  void setCenter(double *d);
  void setNormal(double *d);
  void setArea(double d);
  void setScalar(double d);
  void setType(char *str);
  void setNext(BoundaryConnection *bc);
    
  void printInfo();

  static void BoundaryConnection::convertListToArray();
  static void printBCs();

  int getId();
  int getConditionId();
  double *getDistancePtr();
  double *getCenterPtr();
  double *getNormalPtr();
  char *getTypePtr();
  double getArea();
  double getScalar();
  BoundaryConnection *getNext();
  
  static BoundaryConnection *list;
  static BoundaryConnection *end_of_list;
  static int num_bcs;
  static BoundaryConnection **_array;
  
private:

  // all boundary ids are local nonghosted numbering 
  int idlocal;
  int condition_id;
  double area;
  double center[3],dist[3];
  double normal[3];
  char *type;
  double scalar;
  BoundaryConnection *next;

};

#endif /*BOUNDARYCONNECTION_H_*/
