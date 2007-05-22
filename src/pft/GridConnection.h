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
  void setNormal(double *d);
  void setArea(double d);
  void setFlowFluxCoef(double d);
  void setFlowFluxCoefVSat(double dup, double ddn);
  void printInfo();

  int getIdUpwind();
  int getIdDownwind();
  double *getDistancePtrUpwind();
  double *getDistancePtrDownwind();
  double *getNormalPtr();
  double getArea();
  double getFlowFluxCoef();
  void getFlowFluxCoefVSat(double *dup, double *ddn);

private:
  // all boundary ids are ghosted
  int idup, iddown;
  double area;
  double center[3],distup[3], distdown[3];
  double normal[3];
  
  double flow_flux_coef;
  double flow_flux_coef_vsat_up;
  double flow_flux_coef_vsat_dn;

};

#endif /*GridConnection_H_*/
