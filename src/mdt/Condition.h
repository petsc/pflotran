#ifndef CONDITION_H_
#define CONDITION_H_

#include <string.h>

#include "petscsys.h"

#include "FileIO.h"
#include "Globals.h"

class Condition {
  
public:
  Condition();
  Condition(Condition *old_condition);
  Condition(char *filename);
  virtual ~Condition();

  void nullify();
  void addToList();
  void setId(PetscInt i);
  void setType(PetscInt i);
  void setType(char *str);
  void setDatum(PetscReal *d);
  void setGradient(PetscReal *d);
  void setScalar(PetscReal d);
  void setTimeUnit(char *c);
  void setNext(Condition *cond);

  PetscReal computeHydrostaticPressure(PetscReal *coord);

  static void convertListToArray();
  static void updateConditions(PetscReal time);
  static void initializeConditions();

  void printInfo();

  PetscInt getId();
  PetscInt getType();
  char *getTypePtr();
  PetscReal *getDatumPtr();
  PetscReal *getGradientPtr();
  PetscReal getScalar();
  char *getTimeUnitPtr();
  Condition *getNext();
  char *getName();
  void printConditions();

  PetscReal *scalars;
  PetscReal *times;
  PetscInt max_time_index;
  PetscReal cur_value;

  static Condition *list;
  static Condition *end_of_list;
  static PetscInt num_conditions;
  static Condition **_array;
  static PetscInt initial_condition_id;
  Condition *next;  

private:

  // all boundary ids are local nonghosted numbering 
  PetscInt id;
  PetscInt itype;
  char *ctype;
  char name[32];
  PetscReal datum[3];
  PetscReal gradient[3];
  char time_unit[MAXWORDLENGTH];


  PetscInt cur_time_index;

};

#endif /*CONDITION_H_*/
