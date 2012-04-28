#include "Material.h"

Material::Material() {
  nullify();
}

Material::Material(char *name_) {
  nullify();
  strcpy(name,name_);
}

void Material::nullify() {
  name[0] = '\0';
  hc_horz = 0.;
  hc_vert = 0.;
  perm_horz = 0.;
  perm_vert = 0.;
  porosity = 0.;
  density = 0.;
  air_entry_head = 0.;
  air_entry_pressure = 0.;
  lambda = 0.;
  residual_saturation = 0.;
}

void Material::setName(char *name_) { strcpy(name,name_); }
void Material::setId(PetscInt id_) { id = id_; }
void Material::setHydraulicConductivity_H(PetscReal hc_h) {
  hc_horz = hc_h;
  perm_horz = hc_h*hc_to_perm;
}
void Material::setHydraulicConductivity_V(PetscReal hc_v) {
  hc_vert = hc_v;
  perm_vert = hc_v*hc_to_perm;
}
void Material::setPermeability_H(PetscReal perm_h) {
  perm_horz = perm_h;
  hc_horz = perm_h*perm_to_hc;
}
void Material::setPermeability_V(PetscReal perm_v) {
  perm_vert = perm_v;
  hc_vert = perm_v*perm_to_hc;
}
void Material::setPorosity(PetscReal por) { porosity = por; }
void Material::setDensity(PetscReal den) { density = den; }
void Material::setAirEntryHead(PetscReal air_entry_head_) { 
  air_entry_head = air_entry_head_;
  air_entry_pressure = air_entry_head*density0*gravity;
}
void Material::setAirEntryPressure(PetscReal air_entry_pres_) { 
  air_entry_pressure = air_entry_pres_;
  air_entry_head = air_entry_pres_/(density0*gravity);
}
void Material::setLambda(PetscReal l) { lambda = l; }
void Material::setResidualSaturation(PetscReal Sr) { residual_saturation = Sr; }

void Material::getName(char *name_) { strcpy(name_,name); }
PetscInt Material::getId() { return id; }
PetscReal Material::getHydraulicConductivity_H() { return hc_horz; }
PetscReal Material::getHydraulicConductivity_V() { return hc_vert; }
PetscReal Material::getPermeability_H() { return perm_horz; }
PetscReal Material::getPermeability_V() { return perm_vert; }
PetscReal Material::getPorosity() { return porosity; }
PetscReal Material::getDensity() { return density; }
PetscReal Material::getAirEntryHead() { return air_entry_head; }
PetscReal Material::getAirEntryPressure() { return air_entry_pressure; }
PetscReal Material::getLambda() { return lambda; }
PetscReal Material::getResidualSaturation() { return residual_saturation; }

Material::~Material() {
}
