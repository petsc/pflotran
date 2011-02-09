#ifndef __Reaction_hpp__
#define __Reaction_hpp__

#include <vector>

#include "Species.hpp"

class Reaction : public Species {

 public: 
  Reaction();
  Reaction(const int);
  void Zero(void);
  void CalculatelnQK1(std::vector<Species>);
  void CalculatelnQK2(std::vector<double>);
  void CalculatelnQK3(double *);
  
  double GetResult1(void) const { return this->get_concentration(); }
  double GetResult2or3(void) const { return this->concentration; }

  virtual ~Reaction();

 private:
  
  int ncomp_;
  double lnQK_;
  std::vector<int> species_ids_;
  std::vector<double> stoichiometry_;

  double concentration;

};

#endif // __Reaction_hpp__

