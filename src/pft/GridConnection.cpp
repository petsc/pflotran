#include "GridConnection.h"

GridConnection::GridConnection() {
  area = 0.;
  id1 = -1;
  id2 = -1;
  for (int i=0; i<3; i++) {
    center[i] = 0.;
    dist1[i] = 0.;
    dist2[i] = 0.;
    normal_vector[i] = 0.;
  }
}

void GridConnection::setIdUpwind(int i) { id1 = i; }
void GridConnection::setIdDownwind(int i) { id2 = i; }
void GridConnection::setArea(double d) { area = d; }
void GridConnection::setDistanceUpwind(double *d) { 
  dist1[0] = d[0];
  dist1[1] = d[1];
  dist1[2] = d[2];
}
void GridConnection::setDistanceDownwind(double *d) { 
  dist2[0] = d[0];
  dist2[1] = d[1];
  dist2[2] = d[2];
}

int GridConnection::getIdUpwind() { return id1; }
int GridConnection::getIdDownwind() { return id2; }
double GridConnection::getArea() { return area; }
double *GridConnection::getDistancePtrUpwind() { return &dist1[0]; }
double *GridConnection::getDistancePtrDownwind() { return &dist2[0]; }

void GridConnection::printInfo() {
  printf("id1: %d id2: %d area: %f \n",id1,id2,area);
  printf("dx1: %f dy1: %f dz1: %f \n",dist1[0],dist1[1],dist1[2]);
  printf("dx2: %f dy2: %f dz2: %f\n",dist2[0],dist2[1],dist2[2]);
}

GridConnection::~GridConnection() {

}
