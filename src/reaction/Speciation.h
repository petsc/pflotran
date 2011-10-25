#ifndef SPECIATION_H_
#define SPECIATION_H_

#include "Block.h"
#include "LU.h"
#include "FileIO.h"

#include <iostream>
#include<string>
#include <iomanip>

#include <fstream>
#include <sstream>
using namespace std;

//#define DEBUG

class Speciation {
  
public:
  Speciation(int n);
  Speciation();
  virtual ~Speciation();

  // class function prototypes
  void calculateAQComplexes();
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

public:

static void createCarbonateSystem() {

  ncomp = 2;
  int count = 0;

  primary_species_names = new string[ncomp];
  primary_species_names[count++] = "H+";
  primary_species_names[count++] = "HCO3-";
  count = 0;
  primary_species_Z = new double[ncomp];
  primary_species_Z[0] = 1.;
  primary_species_Z[count++] = -1.;
  count = 0;
  primary_species_a0 = new double[ncomp];
  primary_species_a0[count++] = 9.;
  primary_species_a0[count++] = 4.;

  neqcplx = 3;
  count = 0;
  secondary_species_names = new string[neqcplx];
  secondary_species_names[count++] = "OH-";
  secondary_species_names[count++] = "CO3--";
  secondary_species_names[count++] = "CO2(aq)";
  count = 0;
  secondary_species_Z = new double[neqcplx];
  secondary_species_Z[count++] = -1.;
  secondary_species_Z[count++] = -2.;
  secondary_species_Z[count++] = 0.;
  count = 0;
  secondary_species_a0 = new double[neqcplx];
  secondary_species_a0[count++] = 3.5;
  secondary_species_a0[count++] = 4.5;
  secondary_species_a0[count++] = 3.;

  eqcplxspecid = new int*[neqcplx];
  for (int i=0; i<neqcplx; i++) {
    eqcplxspecid[i] = new int[ncomp+1];
    for (int j=0; j<ncomp+1; j++)
      eqcplxspecid[i][j] = 0;
  }
  eqcplxstoich = new double*[neqcplx];
  for (int i=0; i<neqcplx; i++) {
    eqcplxstoich[i] = new double[ncomp];
    for (int j=0; j<ncomp; j++)
      eqcplxstoich[i][j] = 0.;
  }
  eqcplxh2oid = new int[neqcplx];
  h2ostoich = new double[neqcplx];
  eqcplx_logK = new double[neqcplx];
  for (int i=0; i<neqcplx; i++) {
    eqcplxh2oid[i] = -1;
    h2ostoich[i] = 0.;
    eqcplx_logK[i] = 0.;
  }

  // OH-
  eqcplxspecid[0][0] = 1; // # of components in rxn
  eqcplxspecid[0][1] = 0; // H+ id
  eqcplxstoich[0][0] = -1.; // H+ stoich
  eqcplxh2oid[0] = 1; // id of h2o in rxn
  h2ostoich[0] = -1.; // stoich of h2o in rxn
  eqcplx_logK[0] = 13.9951; // equilibrium constant
  // CO3--
  eqcplxspecid[1][0] = 2; // # of components in rxn
  eqcplxspecid[1][1] = 0; // H+ id
  eqcplxspecid[1][2] = 1; // HCO3- id
  eqcplxstoich[1][0] = -1.; // H+ stoich
  eqcplxstoich[1][1] = 1.; // HCO3- stoich
  eqcplxh2oid[1] = -1; // id of h2o in rxn
  h2ostoich[1] = 0.; // stoich of h2o in rxn
  eqcplx_logK[1] = 10.3288; // equilibrium constant
  // CO2(aq)
  eqcplxspecid[2][0] = 2; // # of components in rxn
  eqcplxspecid[2][1] = 0; // H+ id
  eqcplxspecid[2][2] = 1; // HCO3- id
  eqcplxstoich[2][0] = 1.; // H+ stoich
  eqcplxstoich[2][1] = 1.; // HCO3- stoich
  eqcplxh2oid[2] = 1; // id of h2o in rxn
  h2ostoich[2] = -1.; // stoich of h2o in rxn
  eqcplx_logK[2] = -6.3447; // equilibrium constant

};

static void readReactionFromFile(string filename) {

  char word[32];

  // open file with FileIO buffer
  FileIO *file = new FileIO(filename);
  // first line indicates number of primary and secondary components
  file->getLine();
  file->readInt(&ncomp);
  file->readInt(&neqcplx);

  // allocate arrays
  primary_species_names = new string[ncomp];
  primary_species_Z = new double[ncomp];
  primary_species_a0 = new double[ncomp];
  // initialize to zero
  for (int i=0; i<ncomp; i++) {
    primary_species_names[i] = '\0';
    primary_species_Z[i] = 0.;
    primary_species_a0[i] = 0.;
  }

  secondary_species_names = new string[neqcplx];
  secondary_species_Z = new double[neqcplx];
  secondary_species_a0 = new double[neqcplx];
  // initialize to zero
  for (int i=0; i<neqcplx; i++) {
    secondary_species_names[i] = '\0';
    secondary_species_Z[i] = 0.;
    secondary_species_a0[i] = 0.;
  }

  eqcplxspecid = new int*[neqcplx];
  for (int i=0; i<neqcplx; i++) {
    eqcplxspecid[i] = new int[ncomp+1];
    for (int j=0; j<ncomp+1; j++)
      eqcplxspecid[i][j] = 0;
  }
  eqcplxstoich = new double*[neqcplx];
  for (int i=0; i<neqcplx; i++) {
    eqcplxstoich[i] = new double[ncomp];
    for (int j=0; j<ncomp; j++)
      eqcplxstoich[i][j] = 0.;
  }
  eqcplxh2oid = new int[neqcplx];
  h2ostoich = new double[neqcplx];
  eqcplx_logK = new double[neqcplx];
  for (int i=0; i<neqcplx; i++) {
    eqcplxh2oid[i] = -1;
    h2ostoich[i] = 0.;
    eqcplx_logK[i] = 0.;
  }

  // for ncomp primary speces, each line lists name, Z, a0
  for (int i=0; i<ncomp; i++) {
    file->getLine();
    file->readWord(word);
    primary_species_names[i] = word;
    file->readDouble(&(primary_species_Z[i]));
    file->readDouble(&(primary_species_a0[i]));
#ifdef DEBUG
    cout << primary_species_names[i] << " " << primary_species_Z[i] << " " <<
            primary_species_a0[i] << endl;
#endif
  }

  // for neqcplx secondary speces, each line lists name, Z, a0, etc.
  for (int i=0; i<neqcplx; i++) {
    file->getLine();
    file->readWord(word);
    secondary_species_names[i] = word;
    file->readDouble(&(secondary_species_Z[i]));
    file->readDouble(&(secondary_species_a0[i]));
#ifdef DEBUG
    cout << secondary_species_names[i] << " " << secondary_species_Z[i] << " " <<
            secondary_species_a0[i] << endl;
#endif
    // species ids
    file->getLine();
    int ncomp_in_rxn;
    file->readInt(&ncomp_in_rxn);
    eqcplxspecid[i][0] = ncomp_in_rxn;
    for (int j=1; j<ncomp_in_rxn+1; j++) {
      file->readInt(&(eqcplxspecid[i][j]));
      // decrement for zero-based indexing
      if (j > 0) eqcplxspecid[i][j]--;
    }
#ifdef DEBUG
    for (int j=0; j<ncomp_in_rxn+1; j++)
      cout << eqcplxspecid[i][j] << " " ;
    cout << endl;
#endif
    // species stoichiometries
    double temp;
    file->getLine();
    // skip first value (PFLOTRAN eqcplxstroich mirrors eqcplxspecid)
    file->readDouble(&temp);
    for (int j=0; j<ncomp_in_rxn; j++)
      file->readDouble(&(eqcplxstoich[i][j]));
#ifdef DEBUG
    for (int j=0; j<ncomp_in_rxn; j++)
      cout << eqcplxstoich[i][j] << " " ;
    cout << endl;
#endif
    // h2o id
    file->getLine();
    file->readInt(&(eqcplxh2oid[i]));
    // decrement for zero-based indexing
    eqcplxh2oid[i]--;
#ifdef DEBUG
    cout << eqcplxh2oid[i] << endl;
#endif
    // h2o stoichiometry
    file->getLine();
    file->readDouble(&(h2ostoich[i]));
#ifdef DEBUG
    cout << h2ostoich[i] << endl;
#endif
    // logK
    file->getLine();
    file->readDouble(&(eqcplx_logK[i]));
#ifdef DEBUG
    cout << eqcplx_logK[i] << endl;
    cout << "-------------------------\n";
#endif
  }
  
  delete file;
};

};

#endif /*SPECIATION_H_*/

