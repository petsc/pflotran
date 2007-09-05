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
void Material::setId(int id_) { id = id_; }
void Material::setHydraulicConductivity_H(double hc_h) {
  hc_horz = hc_h;
  perm_horz = hc_h*hc_to_perm;
}
void Material::setHydraulicConductivity_V(double hc_v) {
  hc_vert = hc_v;
  perm_vert = hc_v*hc_to_perm;
}
void Material::setPermeability_H(double perm_h) {
  perm_horz = perm_h;
  hc_horz = perm_h*perm_to_hc;
}
void Material::setPermeability_V(double perm_v) {
  perm_vert = perm_v;
  hc_vert = perm_v*perm_to_hc;
}
void Material::setPorosity(double por) { porosity = por; }
void Material::setDensity(double den) { density = den; }
void Material::setAirEntryHead(double air_entry_head_) { 
  air_entry_head = air_entry_head_;
  air_entry_pressure = air_entry_head*density0*gravity;
}
void Material::setAirEntryPressure(double air_entry_pres_) { 
  air_entry_pressure = air_entry_pres_;
  air_entry_head = air_entry_pres_/(density0*gravity);
}
void Material::setLambda(double l) { lambda = l; }
void Material::setResidualSaturation(double Sr) { residual_saturation = Sr; }

void Material::getName(char *name_) { strcpy(name_,name); }
int Material::getId() { return id; }
double Material::getHydraulicConductivity_H() { return hc_horz; }
double Material::getHydraulicConductivity_V() { return hc_vert; }
double Material::getPermeability_H() { return perm_horz; }
double Material::getPermeability_V() { return perm_vert; }
double Material::getPorosity() { return porosity; }
double Material::getDensity() { return density; }
double Material::getAirEntryHead() { return air_entry_head; }
double Material::getAirEntryPressure() { return air_entry_pressure; }
double Material::getLambda() { return lambda; }
double Material::getResidualSaturation() { return residual_saturation; }

Material::~Material() {
}
