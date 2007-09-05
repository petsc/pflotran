#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <string.h>

#include "include/petsc.h"

#include "globals.h"

class Material {
  
public:

  Material();
  Material(char *name_);
  virtual ~Material();
  void nullify();
  void setName(char *name_);
  void setId(int id_);
  void setHydraulicConductivity_H(double hc_h);
  void setHydraulicConductivity_V(double hc_v);
  void setPermeability_H(double perm_h);
  void setPermeability_V(double perm_v);
  void setPorosity(double por);
  void setDensity(double den);
  void setAirEntryHead(double air_entry_head_);
  void setAirEntryPressure(double air_entry_pres_);
  void setLambda(double l);
  void setResidualSaturation(double Sr);

  void getName(char *name_);
  int getId();
  double getHydraulicConductivity_H();
  double getHydraulicConductivity_V();
  double getPermeability_H();
  double getPermeability_V();
  double getPorosity();
  double getDensity();
  double getAirEntryHead();
  double getAirEntryPressure();
  double getLambda();
  double getResidualSaturation();

private:

  char name[32];
  int id;
  double hc_horz;    // m/d
  double hc_vert;    // m/d
  double perm_horz;  // m^2
  double perm_vert;  // m^2
  double porosity;
  double density;    // kg/m^3
  double air_entry_head; // m
  double air_entry_pressure; // Pa
  double lambda;
  double residual_saturation;

};

#endif /*MATERIAL_H_*/
