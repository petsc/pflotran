#ifndef CONDITION_H_
#define CONDITION_H_

#include <string.h>

#include "include/petsc.h"

#include "FileIO.h"
#include "Globals.h"

class Condition {
  
public:
  Condition();
  Condition::Condition(Condition *old_condition);
  Condition::Condition(char *filename);
  virtual ~Condition();

  void nullify();
  void addToList();
  void setId(int i);
  void setType(int i);
  void setType(char *str);
  void setDatum(double *d);
  void setGradient(double *d);
  void setScalar(double d);
  void setTimeUnit(char *c);
  void setNext(Condition *cond);

  double computeHydrostaticPressure(double *coord);

  static void Condition::convertListToArray();
  static void Condition::updateConditions(double time);
  static void Condition::initializeConditions();

  void printInfo();

  int getId();
  int getType();
  char *getTypePtr();
  double *getDatumPtr();
  double *getGradientPtr();
  double getScalar();
  char *getTimeUnitPtr();
  Condition *getNext();
  char *getName();
  void printConditions();

  double *scalars;
  double *times;
  int max_time_index;
  double cur_value;

  static Condition *list;
  static Condition *end_of_list;
  static int num_conditions;
  static Condition **_array;
  static int initial_condition_id;
  Condition *next;  

private:

  // all boundary ids are local nonghosted numbering 
  int id;
  int itype;
  char *ctype;
  char name[32];
  double datum[3];
  double gradient[3];
  char time_unit[MAXWORDLENGTH];


  int cur_time_index;

};

#endif /*CONDITION_H_*/
