#ifndef SPECIATION_H_
#define SPECIATION_H_

#include "Block.h"
#include "LU.h"

#include <iostream>
#include<string>
#include <iomanip>
using namespace std;

class Speciation {
  
public:
  Speciation(int n);
  virtual ~Speciation();

  // class function prototypes
  void calculateAQComplexes();
  void createCarbonateSystem();
  void calculateTotal();
  void calculateActivityCoefficients(int flag);
  int speciate(double *total);

  // static class variables
  static int ncomp;
  static string *primary_species_names;
  static double *primary_species_Z;
  static double *primary_species_a0;
  static int neqcplx;
  static string *secondary_species_names;
  static double *secondary_species_Z;
  static double *secondary_species_a0;
  static int **eqcplxspecid;
  static double **eqcplxstoich;
  static int *eqcplxh2oid;
  static double *h2ostoich;
  static double *eqcplx_logK;

  static double debyeA;
  static double debyeB;
  static double debyeBdot;

  static double speciation_tolerance;

private:

  // instance variables
  double *total;
  Block *dtotal;
  double *pri_molal;
  double *sec_molal;
  double *pri_act;
  double *sec_act;

};

#endif /*SPECIATION_H_*/