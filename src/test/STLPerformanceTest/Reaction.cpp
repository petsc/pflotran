
#include "Reaction.hpp"

Reaction::Reaction() : Species(),
                       lnQK_(0.),
                       ncomp_(0),
                       concentration(0.)
{
  species_ids_.clear();
  stoichiometry_.clear();
}

Reaction::Reaction(const int ncomp) : Species(),
                       lnQK_(0.),
                       ncomp_(ncomp),
                       concentration(0.)
{
  species_ids_.clear();
  stoichiometry_.clear();

  for (int i = 0; i < ncomp; i++) {
    species_ids_.push_back(i);
  }
  for (int i = 0; i < ncomp; i++) {
    stoichiometry_.push_back((i % 2 == 0 ? 1.*(i+1) : -1.*(i+1)));
  }
}

void Reaction::Zero(void)
{
  lnQK_ = 0.;
  set_concentration(0.);
  concentration = 0.;
}

void Reaction::CalculatelnQK1(std::vector<Species> species)
{
  double lnQK = 0.;
  for (int i = 0; i < ncomp_; i++) {
    lnQK += stoichiometry_[i] * species[species_ids_[i]].get_ln_concentration();
  }
  lnQK_ = lnQK;
  set_concentration(std::exp(lnQK));
}

void Reaction::CalculatelnQK2(std::vector<double> ln_concentrations)
{
  double lnQK = 0.;
  for (int i = 0; i < ncomp_; i++) {
    lnQK += stoichiometry_[i] * ln_concentrations[species_ids_[i]];
  }
  lnQK_ = lnQK;
  concentration = std::exp(lnQK);
}


void Reaction::CalculatelnQK3(double *ln_concentrations)
{
  double lnQK = 0.;
  for (int i = 0; i < ncomp_; i++) {
    lnQK += stoichiometry_[i] * ln_concentrations[species_ids_[i]];
  }
  lnQK_ = lnQK;
  concentration = std::exp(lnQK);
}

Reaction::~Reaction()
{
}

