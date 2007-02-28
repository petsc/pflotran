#include "Source.h"

Source::Source() {
  idlocal = -1;
  num_srcs++;
  if (!list) list = this;
  if (end_of_list) end_of_list->next = this;
  end_of_list = this;
//  printf("%d %d %d\n",list,end_of_list,this);
  next = NULL;
}

void Source::setId(int i) { idlocal = i; }
void Source::setScalar(double d) { scalar = d; }
void Source::setType(char *str) { 
  type = new char[32]; 
  strcpy(type,str);
}

int Source::getId() { return idlocal; }
double Source::getScalar() { return scalar; }
char *Source::getTypePtr() { return &type[0]; }
Source *Source::getNext() { return next; };
 
void Source::printInfo() {
  printf("\n  id: %5d\n",idlocal);
  printf("type: %s scalar: %9.3e\n",type,scalar);
}

void Source::printSrcs() {
  PetscErrorCode ierr;
 
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD);
#include "Globals.h" 
  if (myrank == 0) printf("\nSource Connections\n");
  printf("Processor[%d]: %d connections\n",myrank,num_srcs);
  Source *cur_src = list;
  while (cur_src) {
    cur_src->printInfo();
    cur_src = cur_src->next;
  }

  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD);
}

Source::~Source() {
  if (type) delete [] type;
}
