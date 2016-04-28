#ifndef __Species_hpp__
#define __Species_hpp__

#include <cmath>

class Species {

 public: 
  Species();
  Species(const double);
  virtual ~Species();

  double get_concentration(void) const { return this->concentration_; }
  double get_ln_concentration(void) const { return this->ln_concentration_; }
  void set_concentration(const double d) { this->concentration_ = d; }
  void set_ln_concentration(const double d) { this->ln_concentration_ = d; }

 private:
  double concentration_;
  double ln_concentration_;
};

#endif // __Species_hpp__

