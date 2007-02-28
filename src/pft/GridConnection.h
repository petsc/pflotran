#ifndef GridConnection_H_
#define GridConnection_H_

#include "include/petsc.h"

class GridConnection {
  
public:
  GridConnection();
  virtual ~GridConnection();
  
  void setIdUpwind(int i);
  void setIdDownwind(int i);
  void setDistanceUpwind(double *d);
  void setDistanceDownwind(double *d);
  void setArea(double d);
  void printInfo();

  int getIdUpwind();
  int getIdDownwind();
  double *getDistancePtrUpwind();
  double *getDistancePtrDownwind();
  double getArea();

private:
  // all boundary ids are ghosted
  int id1, id2;
  double area;
  double center[3],dist1[3], dist2[3];
  double normal_vector[3];

};

#endif /*GridConnection_H_*/
