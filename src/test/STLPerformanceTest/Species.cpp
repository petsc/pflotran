#include "Species.hpp"

Species::Species() {
}

Species::Species(const double concentration) : concentration_(concentration),
                                               ln_concentration_(std::log(concentration))
{
}

Species::~Species()
{
}

