#include "BoundaryConnection.h"

BoundaryConnection::BoundaryConnection() {
  area = 0.;
  idlocal = -1;
  condition_id = -1;
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

void BoundaryConnection::convertListToArray() {

  int count = 0;
  BoundaryConnection::_array = new BoundaryConnection *[num_bcs];

  BoundaryConnection *cur_bc = list;
  while (cur_bc) {
    _array[count++] = cur_bc;
    cur_bc = cur_bc->next;
  }
}

void BoundaryConnection::setId(int i) { idlocal = i; }
void BoundaryConnection::setConditionId(int i) { condition_id = i; }
void BoundaryConnection::setArea(double d) { area = d; }
void BoundaryConnection::setScalar(double d) { scalar = d; }
void BoundaryConnection::setType(char *str) { 
  type = new char[32]; 
  strcpy(type,str);
}
void BoundaryConnection::setDistance(double *d) { 
  dist[0] = d[0];
  dist[1] = d[1];
  dist[2] = d[2];
}
void BoundaryConnection::setCenter(double *d) { 
  center[0] = d[0];
  center[1] = d[1];
  center[2] = d[2];
}
void BoundaryConnection::setNormal(double *d) { 
  normal[0] = d[0];
  normal[1] = d[1];
  normal[2] = d[2];
}

int BoundaryConnection::getId() { return idlocal; }
int BoundaryConnection::getConditionId() { return condition_id; }
int BoundaryConnection::getConditionType() { 
  return Condition::_array[condition_id]->getType();
}
double BoundaryConnection::getArea() { return area; }
double BoundaryConnection::getScalar() { return scalar; }
char *BoundaryConnection::getTypePtr() { return &type[0]; }
double *BoundaryConnection::getDistancePtr() { return &dist[0]; }
double *BoundaryConnection::getCenterPtr() { return &center[0]; }
double *BoundaryConnection::getNormalPtr() { return &normal[0]; }
BoundaryConnection *BoundaryConnection::getNext() { return next; };

void BoundaryConnection::printInfo() {
/*
    int idlocal;
  double area;
  double center[3],dist[3];
  double normal[3];
  char *type;
  double scalar;
  BoundaryConnection *next;
*/  
  printf("\n  id: %5d\n",idlocal);
  printf("area: %9.3f\n",area);
  printf("type: %s scalar: %9.3e\n",type,scalar);
  printf("  dx: %9.3f dy: %9.3f dz: %9.3f\n",dist[0],dist[1],dist[2]);
  printf("  dx: %9.3f dy: %9.3f dz: %9.3f\n",normal[0],normal[1],normal[2]);
}

void BoundaryConnection::printBCs() {
  PetscErrorCode ierr;
 
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
#include "Globals.h" 
  if (myrank == 0) printf("\nBoundary Connections\n");
  printf("Processor[%d]: %d connections\n",myrank,num_bcs);
  BoundaryConnection *cur_bc = list;
  while (cur_bc) {
    cur_bc->printInfo();
    cur_bc = cur_bc->next;
  }

  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

BoundaryConnection::~BoundaryConnection() {
  if (type) delete [] type;
}
