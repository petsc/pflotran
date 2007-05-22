#ifndef GLOBALS_H_
#define GLOBALS_H_

extern int myrank;

const double density0 = 998.32;
const double viscosity0 = 8.9e-4;
//const double gravity = 9.8068;
const double gravity = 9.81;
//const double gravity = 0.;
const double patm = 101325.;
const double p_threshhold = 0.03*gravity*density0;
const double permeability0 = 1.e-6;
const double porosity0 = 0.33;
const double betap = -4./(density0*gravity);

#endif /*GLOBALS_H_*/
