#ifndef SOURCE_H_
#define SOURCE_H_

#include <string.h>

#include "include/petsc.h"

class Source {
  
public:
  Source();
  virtual ~Source();

  void setId(int i);
  void setScalar(double d);
  void setType(char *str);
  void setNext(Source *src);
  
  void printInfo();
  static void printSrcs();

  int getId();
  char *getTypePtr();
  double getScalar();
  Source *getNext();

  static Source *list;
  static Source *end_of_list;
  static int num_srcs;
  
private:

  // all boundary ids are local nonghosted numbering 
  int idlocal;
  char *type;
  double scalar;
  Source *next;

};

#endif /*SOURCE_H_*/
