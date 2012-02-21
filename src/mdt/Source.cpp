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

void Source::setId(PetscInt i) { idlocal = i; }
void Source::setScalar(PetscReal d) { scalar = d; }
void Source::setType(char *str) { 
  type = new char[32]; 
  strcpy(type,str);
}

PetscInt Source::getId() { return idlocal; }
PetscReal Source::getScalar() { return scalar; }
char *Source::getTypePtr() { return &type[0]; }
Source *Source::getNext() { return next; };
 
void Source::printInfo() {
  printf("\n  id: %5d\n",idlocal);
  printf("type: %s scalar: %9.3e\n",type,scalar);
}

void Source::printSrcs() {
  PetscErrorCode ierr;
 
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
#include "Globals.h" 
  if (myrank == 0) printf("\nSource Connections\n");
  printf("Processor[%d]: %d connections\n",myrank,num_srcs);
  Source *cur_src = list;
  while (cur_src) {
    cur_src->printInfo();
    cur_src = cur_src->next;
  }

  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

Source::~Source() {
  if (type) delete [] type;
}
