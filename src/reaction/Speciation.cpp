#include "Speciation.h"

  // initialization of static class variables
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

// constructor
Speciation::Speciation(int n) {
  ncomp = n;
  total = new double[ncomp];
  dtotal = new Block(ncomp);
  pri_molal = new double[ncomp];
  sec_molal = new double[neqcplx];
  pri_act = new double[ncomp];
  sec_act = new double[neqcplx];
}

Speciation::Speciation() {
  total = new double[ncomp];
  dtotal = new Block(ncomp);
  pri_molal = new double[ncomp];
  sec_molal = new double[neqcplx];
  pri_act = new double[ncomp];
  sec_act = new double[neqcplx];
}

void Speciation::calculateAQComplexes() {
  cout << "Calculated!\n";
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

  // scale by density of water to convert total from molality to molarity
  for (int i=0; i<ncomp; i++)
    total[i] *= den_kg_per_L;
  dtotal->scale(den_kg_per_L);
 
  delete [] ln_conc;
  delete [] ln_act;
}

void Speciation::calculateActivityCoefficients(int flag) {

  double I; // ionic strength

  if (flag < 0) {  // activity coefficients = 1.
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

    // still need to add calculation for activity of water

  }
  
}

int Speciation::speciate(double *target_total) {

  cout << endl;
  cout << "Target Total Component Concentrations\n";
  for (int i=0; i<ncomp; i++) {
    cout << setw(10) << primary_species_names[i];
    cout << ": " << target_total[i] << endl;
  }
  cout << endl;

  // initialize free-ion concentration s
  for (int i=0; i<ncomp; i++) 
    pri_molal[i] = 1.e-9;

  // allocate arrays for Newton-Raphson
  double *residual = new double[ncomp];
  double *rhs = new double[ncomp];
  double *prev_molal = new double[ncomp];
  double *update = new double[ncomp];
  Block *J = new Block(ncomp);

  // allocate pivoting array for LU
  int *indices = new int[ncomp];

  double max_rel_change;
  int num_iterations = 0;

  do {
    
    calculateActivityCoefficients(-1);
    calculateTotal();

    // add derivatives of total with respect to free to Jacobian
    J->zero();
    J->addValues(0,0,dtotal);

    // calculate residual
    for (int i=0; i<ncomp; i++)
      residual[i] = total[i]-target_total[i];

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
    // for derivatives with respect to ln concentration for log formulation
    for (int i=0; i<ncomp; i++)
      J->scaleColumn(i,pri_molal[i]);

#ifdef DEBUG
    cout << "before solve\n";
    J->print();
#endif
    // LU direct solve
    double D;
    ludcmp(J->getValues(),ncomp,indices,&D);
    lubksb(J->getValues(),ncomp,indices,rhs);

    // calculate update truncating at a maximum of 5 in log space
    for (int i=0; i<ncomp; i++) {
      update[i] = rhs[i] > 0. ? 
        (rhs[i] > 5. ? 5. : rhs[i]) : (rhs[i] < -5. ? -5. : rhs[i]);
      prev_molal[i] = pri_molal[i];
      pri_molal[i] *= exp(-update[i]);
    }

    // calculate maximum relative change in concentration over all species
    max_rel_change = 0.;
    for (int i=0; i<ncomp; i++) {
      double delta = fabs(pri_molal[i]-prev_molal[i])/prev_molal[i];
      max_rel_change = delta > max_rel_change ? delta : max_rel_change;
    }

#ifdef DEBUG
    for (int i=0; i<ncomp; i++)
      cout << primary_species_names[i] << " " << pri_molal[i] << " " << total[i] << "\n";
#endif

    num_iterations++;

    // exist if maximum relative change is below tolerance
  } while (max_rel_change > speciation_tolerance);

  // output for testing purposes
  cout << endl;
  cout << "Primary Species ---------------------\n";
  for (int i=0; i<ncomp; i++) {
    cout << "  " << primary_species_names[i] << endl;
    cout << "       Total: " << total[i] << endl;
    cout << "    Free-Ion: " << pri_molal[i] << endl;
  }
  cout << endl;
  cout << "Secondary Species -------------------\n";
  for (int i=0; i<neqcplx; i++) {
    cout << "  " << secondary_species_names[i] << endl;
    cout << "    Free-Ion: " << sec_molal[i] << endl;
  }
  cout << "-------------------------------------\n";
  cout << endl;

  // free up memory
  delete J;
  delete [] residual;
  delete [] rhs;
  delete [] update;
  delete [] prev_molal;
  delete [] indices;

  return num_iterations;
}

// Destructor for freeing instance variables
Speciation::~Speciation() {
  if (dtotal) delete dtotal;
  if (total) delete [] total;
  if (pri_molal) delete [] pri_molal;
  if (sec_molal) delete [] sec_molal;
  if (pri_act) delete [] pri_act;
  if (sec_act) delete [] sec_act;
}
