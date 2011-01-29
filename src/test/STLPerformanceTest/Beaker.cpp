
#include "Beaker.hpp"

Beaker::Beaker() : ncomp_(0)
{
  species_.clear();
  reactions_.clear();

  concentrations_.clear();
  ln_concentrations_.clear();

  concentration_array_ = NULL;
  ln_concentration_array_ = NULL;
}

void Beaker::Setup(const int num_calls, const int ncomp) {
  
  num_calls_ = num_calls;
  ncomp_ = ncomp;

  for (int i = 0; i < ncomp_; i++) {
    species_.push_back(Species((float)(i+1)*1.e-5));
    reactions_.push_back(Reaction(i+1));
  }

  concentrations_.resize(ncomp);
  ln_concentrations_.resize(ncomp);

  for (int i = 0; i < ncomp_; i++) {
    concentrations_[i] = (float)(i+1)*1.e-5;
    ln_concentrations_[i] = std::log(concentrations_[i]);
  }

  concentration_array_ = new double[ncomp_];
  ln_concentration_array_ = new double[ncomp_];

  for (int i = 0; i < ncomp_; i++) {
    concentration_array_[i] = (float)(i+1)*1.e-5;
    ln_concentration_array_[i] = std::log(concentration_array_[i]);
  }

}

void Beaker::Test1(void) {
  for (int i = 0; i < num_calls_; i++) {
    for (std::vector<Reaction>:: iterator r = reactions_.begin();
         r != reactions_.end(); r++) {
      r->CalculatelnQK1(species_);
    }
  }
}

void Beaker::Test2(void) {
  for (int i = 0; i < num_calls_; i++) {
    for (std::vector<Reaction>:: iterator r = reactions_.begin();
         r != reactions_.end(); r++) {
      r->CalculatelnQK2(ln_concentrations_);
    }
  }
}

void Beaker::Test3(void) {
  for (int i = 0; i < num_calls_; i++) {
    for (std::vector<Reaction>:: iterator r = reactions_.begin();
         r != reactions_.end(); r++) {
      r->CalculatelnQK3(ln_concentration_array_);
    }
  }
}

double Beaker::GetSummationOfResults1(void) {
  double sum = 0.;
  for (std::vector<Reaction>:: iterator r = reactions_.begin();
       r != reactions_.end(); r++) {
    sum += std::log(r->GetResult1());
  }
  return sum;
}

double Beaker::GetSummationOfResults2or3(void) {
  double sum = 0.;
  for (std::vector<Reaction>:: iterator r = reactions_.begin();
       r != reactions_.end(); r++) {
    sum += std::log(r->GetResult2or3());
  }
  return sum;
}

Beaker::~Beaker()
{
  if (concentration_array_) delete concentration_array_;
  concentration_array_ = NULL;
  if (ln_concentration_array_) delete ln_concentration_array_;
  ln_concentration_array_ = NULL;
}

