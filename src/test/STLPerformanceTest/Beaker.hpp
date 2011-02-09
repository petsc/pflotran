#ifndef __Beaker_hpp__
#define __Beaker_hpp__

#include <cmath>

#include "Species.hpp"
#include "Reaction.hpp"

class Beaker : Species {

 public: 
  Beaker();
  void Setup(const int, const int);
  void Test1(void);
  void Test2(void);
  void Test3(void);
  double GetSummationOfResults1(void);
  double GetSummationOfResults2or3(void);
  virtual ~Beaker();

 private:

   int num_calls_;

   int ncomp_;
   std::vector<Species> species_;
   std::vector<Reaction> reactions_;

   std::vector<double> concentrations_;
   std::vector<double> ln_concentrations_;

   double *concentration_array_;
   double *ln_concentration_array_;
  
};

#endif // __Beaker_hpp__

