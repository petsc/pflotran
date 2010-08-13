#include "Speciation.h"

  int Speciation::ncomp = 0;

  int Speciation::neqcplx = 0;
  string *Speciation::primary_species_names = NULL;
  double *Speciation::primary_species_Z = NULL;
  double *Speciation::primary_species_a0 = NULL;

  string *Speciation::secondary_species_names = NULL;
  double *Speciation::secondary_species_Z = NULL;
  double *Speciation::secondary_species_a0 = NULL;
  int **Speciation::eqcplxspecid = NULL;
  double **Speciation::eqcplxstoich = NULL;
  int *Speciation::eqcplxh2oid = NULL;
  double *Speciation::h2ostoich = NULL;
  double *Speciation::eqcplx_logK = NULL;

  const double LOG_TO_LN = 2.30258509299;

  // based on 25C
  double Speciation::debyeA = 0.4939;
  double Speciation::debyeB = 0.3253;
  double Speciation::debyeBdot = 0.0374;

  // tolerance for speciation convergence
  double Speciation::speciation_tolerance = 1.e-12;

Speciation::Speciation(int n) {

  ncomp = n;
  total = new double[ncomp];
  dtotal = new Block(ncomp);
  pri_molal = new double[ncomp];
  sec_molal = new double[ncomp];
  pri_act = new double[ncomp];
  sec_act = new double[ncomp];
}

void Speciation::calculateAQComplexes() {
  cout << "Calculated!\n";
}

void Speciation::createCarbonateSystem() {

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
    eqcplxspecid[i] = new int[3];
    for (int j=0; j<3; j++)
      eqcplxspecid[i][j] = 0;
  }
  eqcplxstoich = new double*[neqcplx];
  for (int i=0; i<neqcplx; i++) {
    eqcplxstoich[i] = new double[3];
    for (int j=0; j<3; j++)
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

}

void Speciation::calculateTotal() {

  const double den_kg_per_L = 1.;
  double *ln_conc = new double[ncomp];
  double *ln_act = new double[ncomp];

  for (int i=0; i<ncomp; i++) {
    total[i] = pri_molal[i];
    ln_conc[i] = log(pri_molal[i]);
    ln_act[i] = ln_conc[i]+log(pri_act[i]);
  }
  // zero dtotal and set diagonal = 1.
  dtotal->zero();
  dtotal->setDiagonal(1.);

  // calculate aqueous complex concentrations
  for (int icplx=0; icplx<neqcplx; icplx++) {
    double lnQK = -eqcplx_logK[icplx]*LOG_TO_LN;
    //lnQK += eqcplxh2ostoich(icplx)*ln_act_h2o;
    int ncomp = eqcplxspecid[icplx][0];
    for (int i=0; i<ncomp; i++) {
      int icomp = eqcplxspecid[icplx][i+1]; // add 1 due to 1-based indexing
      lnQK += eqcplxstoich[icplx][i]*ln_act[icomp];
    }
    sec_molal[icplx] = exp(lnQK)/sec_act[icplx];
 
     // add contribution to primary totals
    for (int i=0; i<ncomp; i++) {
      total[eqcplxspecid[icplx][i+1]] += eqcplxstoich[icplx][i]*sec_molal[icplx];
    }

    // add derivative with repect to free
    for (int j=0; j<ncomp; j++) {
      int jcomp = eqcplxspecid[icplx][j+1]; // add 1 due to 1-based indexing
      double tempd = eqcplxstoich[icplx][j]*exp(lnQK-ln_conc[jcomp])/ 
                     sec_act[icplx];
      for (int i=0; i<ncomp; i++) {
        int icomp = eqcplxspecid[icplx][i+1]; // add 1 due to 1-based indexing
        dtotal->addValue(icomp,jcomp,eqcplxstoich[icplx][i]*tempd);
      }
    }

  }

  for (int i=0; i<ncomp; i++)
    total[i] *= den_kg_per_L;
  dtotal->scale(den_kg_per_L);
 
  delete [] ln_conc;
  delete [] ln_act;
}

void Speciation::calculateActivityCoefficients(int flag) {

  double I; // ionic strength

  if (flag < 0) {
    for (int i=0; i<ncomp; i++)
      pri_act[i] = 1.;
    for (int i=0; i<neqcplx; i++)
      sec_act[i] = 1.;
  }
  else if (flag == 0) {
    cout << "Add Newton's methods for calculating activity coeffiencts\n";
    exit(1);
  }
  else {
    // calculate ionic strength
    // primary species
    I = 0.;
    for (int i=0; i<ncomp; i++)
      I += pri_molal[i]*primary_species_Z[i]*primary_species_Z[i];
    // secondary species
    for (int i=0; i<neqcplx; i++)
      I += sec_molal[i]*secondary_species_Z[i]*secondary_species_Z[i];
    I *= 0.5;
    double sqrt_I = sqrt(I);
    // calculate activity coefficients
    // primary species
    for (int i=0; i<ncomp; i++) {
      double Z = primary_species_Z[i];
      if (fabs(Z) > 1.-10) {
        pri_act[i] = exp(-debyeA*Z*Z*sqrt_I/
                         (1.+primary_species_a0[i]*debyeB*sqrt_I)*
                         LOG_TO_LN);
      }
      else {
        pri_act[i] = 1.;
      }
    }
    // secondary species
    for (int i=0; i<neqcplx; i++) {
      double Z = secondary_species_Z[i];
      if (fabs(Z) > 1.-10) {
        sec_act[i] = exp(-debyeA*Z*Z*sqrt_I/
                         (1.+secondary_species_a0[i]*debyeB*sqrt_I)*
                         LOG_TO_LN);
      }
      else {
        sec_act[i] = 1.;
      }
    }

    // need to add calculation for activity of water

  }
  
}

int Speciation::speciate(double *target_total) {

    // initialize free-ion concentration s
  for (int i=0; i<ncomp; i++) 
    pri_molal[i] = 1.e-9;

  double *residual = new double[ncomp];
  double *rhs = new double[ncomp];
  double *prev_molal = new double[ncomp];
  double *update = new double[ncomp];
  Block *J = new Block(ncomp);

  int *indices = new int[ncomp];

  double max_rel_change;
  int num_iterations = 0;
  do {
    calculateActivityCoefficients(-1);
    calculateTotal();
    J->zero();
    J->addValues(0,0,dtotal);
    for (int i=0; i<ncomp; i++)
      residual[i] = total[i]-target_total[i];

//#define DEBUG
#ifdef DEBUG
    cout << "before scale\n";
    J->print();
#endif
    // scale the Jacobian
    for (int i=0; i<ncomp; i++) {
      double max = J->getRowAbsMax(i);
      if (max > 1.) {
        double scale = 1./max;
        rhs[i] = residual[i]*scale;
        J->scaleRow(i,scale);
      }
      else {
        rhs[i] = residual[i];
      }
    }

#ifdef DEBUG
    cout << "after scale\n";
    J->print();
#endif
    // for derivatives with respect to ln concentration
    for (int i=0; i<ncomp; i++)
      J->scaleColumn(i,pri_molal[i]);

#ifdef DEBUG
    cout << "before solve\n";
    J->print();
#endif
    double D;
    ludcmp(J->getValues(),ncomp,indices,&D);
    lubksb(J->getValues(),ncomp,indices,rhs);

    for (int i=0; i<ncomp; i++) {
      update[i] = rhs[i] > 0. ? 
        (rhs[i] > 5. ? 5. : rhs[i]) : (rhs[i] < -5. ? -5. : rhs[i]);
      prev_molal[i] = pri_molal[i];
      pri_molal[i] *= exp(-update[i]);
    }

    max_rel_change = 0.;
    for (int i=0; i<ncomp; i++) {
      double delta = fabs(pri_molal[i]-prev_molal[i])/prev_molal[i];
      max_rel_change = delta > max_rel_change ? delta : max_rel_change;
    }

    for (int i=0; i<ncomp; i++)
      cout << primary_species_names[i] << " " << pri_molal[i] << " " << total[i] << "\n";

    num_iterations++;

  } while (max_rel_change > speciation_tolerance);

  delete J;
  delete [] residual;
  delete [] rhs;
  delete [] update;
  delete [] prev_molal;
  delete [] indices;

  return num_iterations;
}

Speciation::~Speciation() {
  if (dtotal) delete dtotal;
  if (total) delete [] total;
  if (pri_molal) delete [] pri_molal;
  if (sec_molal) delete [] sec_molal;
  if (pri_act) delete [] pri_act;
  if (sec_act) delete [] sec_act;
}
