#ifndef SOURCE_H_
#define SOURCE_H_

#include <string.h>

#include "petscsys.h"

class Source {
  
public:
  Source();
  virtual ~Source();

  void setId(PetscInt i);
  void setScalar(PetscReal d);
  void setType(char *str);
  void setNext(Source *src);
  
  void printInfo();
  static void printSrcs();

  PetscInt getId();
  char *getTypePtr();
  PetscReal getScalar();
  Source *getNext();

  static Source *list;
  static Source *end_of_list;
  static PetscInt num_srcs;
  
private:

  // all boundary ids are local nonghosted numbering 
  PetscInt idlocal;
  char *type;
  PetscReal scalar;
  Source *next;

};

#endif /*SOURCE_H_*/
