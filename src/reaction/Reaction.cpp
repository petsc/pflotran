#include "Speciation.h"


int main (int argc, char **args) {

#if 1
  // currently hardwired to 2 dof carbonate system (H+, HCO3-)
  Speciation::createCarbonateSystem();
  Speciation *speciation = new Speciation();
  double *total = new double[speciation->ncomp];
  total[0] = 1.e-6;
  total[1] = 1.e-3;
#else
  string filename("reaction.dat");
  Speciation::readReactionFromFile(filename);
  Speciation *speciation = new Speciation();

  // open file with FileIO buffer
  double *total = new double[speciation->ncomp];
  string filename2("target_total.dat");
  FileIO *file = new FileIO(filename2);
  // first line indicates number of primary and secondary components
  char word[32];
  double temp;
  for (int i=0; i<speciation->ncomp; i++) {
    file->getLine();
    file->readWord(word);
    file->readDouble(&temp);
    file->readDouble(&total[i]);
  }
  delete file;
#endif

  speciation->speciate(total);

  delete speciation;
  delete [] total;

  cout << "Done!\n";

}
