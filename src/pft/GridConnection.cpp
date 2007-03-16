#include "GridConnection.h"

GridConnection::GridConnection() {
  area = 0.;
  idup = -1;
  iddown = -1;
  for (int i=0; i<3; i++) {
    center[i] = 0.;
    distup[i] = 0.;
    distdown[i] = 0.;
    normal[i] = 0.;
  }
}

void GridConnection::setIdUpwind(int i) { idup = i; }
void GridConnection::setIdDownwind(int i) { iddown = i; }
void GridConnection::setArea(double d) { area = d; }
void GridConnection::setDistanceUpwind(double *d) { 
  distup[0] = d[0];
  distup[1] = d[1];
  distup[2] = d[2];
}
void GridConnection::setDistanceDownwind(double *d) { 
  distdown[0] = d[0];
  distdown[1] = d[1];
  distdown[2] = d[2];
}
void GridConnection::setNormal(double *d) { 
  normal[0] = d[0];
  normal[1] = d[1];
  normal[2] = d[2];
}
void GridConnection::setFlowFluxCoef(double d) { flow_flux_coef = d; }

int GridConnection::getIdUpwind() { return idup; }
int GridConnection::getIdDownwind() { return iddown; }
double GridConnection::getArea() { return area; }
double *GridConnection::getDistancePtrUpwind() { return &distup[0]; }
double *GridConnection::getDistancePtrDownwind() { return &distdown[0]; }
double *GridConnection::getNormalPtr() { return &normal[0]; }
double GridConnection::getFlowFluxCoef() { return flow_flux_coef; }

void GridConnection::printInfo() {
  printf("id1: %d id2: %d area: %f \n",idup,iddown,area);
  printf("dx1: %f dy1: %f dz1: %f \n",distup[0],distup[1],distup[2]);
  printf("dx2: %f dy2: %f dz2: %f\n",distdown[0],distdown[1],distdown[2]);
}

GridConnection::~GridConnection() {

}
