#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition() {
  area = 0.;
  idlocal = -1;
  for (int i=0; i<3; i++) {
    dist[i] = 0.;
    center[i] = 0.;
    normal[i] = 0.;
  }
  num_bcs++;
  if (!list) list = this;
  if (end_of_list) end_of_list->next = this;
  end_of_list = this;
//  printf("%d %d %d\n",list,end_of_list,this);
  next = NULL;
}

void BoundaryCondition::setId(int i) { idlocal = i; }
void BoundaryCondition::setArea(double d) { area = d; }
void BoundaryCondition::setScalar(double d) { scalar = d; }
void BoundaryCondition::setType(char *str) { 
  type = new char[32]; 
  strcpy(type,str);
}
void BoundaryCondition::setDistance(double *d) { 
  dist[0] = d[0];
  dist[1] = d[1];
  dist[2] = d[2];
}
void BoundaryCondition::setCenter(double *d) { 
  center[0] = d[0];
  center[1] = d[1];
  center[2] = d[2];
}
void BoundaryCondition::setNormal(double *d) { 
  normal[0] = d[0];
  normal[1] = d[1];
  normal[2] = d[2];
}
void BoundaryCondition::setFlowFluxCoef(double d) { flow_flux_coef = d; }

int BoundaryCondition::getId() { return idlocal; }
double BoundaryCondition::getArea() { return area; }
double BoundaryCondition::getScalar() { return scalar; }
char *BoundaryCondition::getTypePtr() { return &type[0]; }
double *BoundaryCondition::getDistancePtr() { return &dist[0]; }
double *BoundaryCondition::getCenterPtr() { return &center[0]; }
double *BoundaryCondition::getNormalPtr() { return &normal[0]; }
BoundaryCondition *BoundaryCondition::getNext() { return next; };
double BoundaryCondition::getFlowFluxCoef() { return flow_flux_coef; }
 
void BoundaryCondition::printInfo() {
/*
    int idlocal;
  double area;
  double center[3],dist[3];
  double normal[3];
  char *type;
  double scalar;
  BoundaryCondition *next;
*/  
  printf("\n  id: %5d\n",idlocal);
  printf("area: %9.3f\n",area);
  printf("type: %s scalar: %9.3e\n",type,scalar);
  printf("  dx: %9.3f dy: %9.3f dz: %9.3f\n",dist[0],dist[1],dist[2]);
  printf("  dx: %9.3f dy: %9.3f dz: %9.3f\n",normal[0],normal[1],normal[2]);
}

void BoundaryCondition::printBCs() {
  PetscErrorCode ierr;
 
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
#include "Globals.h" 
  if (myrank == 0) printf("\nBoundary Connections\n");
  printf("Processor[%d]: %d connections\n",myrank,num_bcs);
  BoundaryCondition *cur_bc = list;
  while (cur_bc) {
    cur_bc->printInfo();
    cur_bc = cur_bc->next;
  }

  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

BoundaryCondition::~BoundaryCondition() {
  if (type) delete [] type;
}
