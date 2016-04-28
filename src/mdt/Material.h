#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <string.h>

#include "petscsys.h"

#include "Globals.h"

class Material {
  
public:

  Material();
  Material(char *name_);
  virtual ~Material();
  void nullify();
  void setName(char *name_);
  void setId(PetscInt id_);
  void setHydraulicConductivity_H(PetscReal hc_h);
  void setHydraulicConductivity_V(PetscReal hc_v);
  void setPermeability_H(PetscReal perm_h);
  void setPermeability_V(PetscReal perm_v);
  void setPorosity(PetscReal por);
  void setDensity(PetscReal den);
  void setAirEntryHead(PetscReal air_entry_head_);
  void setAirEntryPressure(PetscReal air_entry_pres_);
  void setLambda(PetscReal l);
  void setResidualSaturation(PetscReal Sr);

  void getName(char *name_);
  PetscInt getId();
  PetscReal getHydraulicConductivity_H();
  PetscReal getHydraulicConductivity_V();
  PetscReal getPermeability_H();
  PetscReal getPermeability_V();
  PetscReal getPorosity();
  PetscReal getDensity();
  PetscReal getAirEntryHead();
  PetscReal getAirEntryPressure();
  PetscReal getLambda();
  PetscReal getResidualSaturation();

private:

  char name[32];
  PetscInt id;
  PetscReal hc_horz;    // m/d
  PetscReal hc_vert;    // m/d
  PetscReal perm_horz;  // m^2
  PetscReal perm_vert;  // m^2
  PetscReal porosity;
  PetscReal density;    // kg/m^3
  PetscReal air_entry_head; // m
  PetscReal air_entry_pressure; // Pa
  PetscReal lambda;
  PetscReal residual_saturation;

};

#endif /*MATERIAL_H_*/
