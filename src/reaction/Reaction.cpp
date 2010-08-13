#include "Speciation.h"


int main (int argc, char **args) {

  Speciation *speciation = new Speciation(2);
  speciation->createCarbonateSystem();

  // currently hardwired to 2 dof carbonate system (H+, HCO3-)
  double *total = new double[speciation->ncomp];
  total[0] = 1.e-6;
  total[1] = 1.e-3;
  speciation->speciate(total);

  cout << "Done!\n";

}