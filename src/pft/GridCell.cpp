#include "GridCell.h"

const double GridCell::beta = -4.;  // m^-1
Saturation GridCell::saturation;

GridCell::GridCell() {
  
  id_local = -1;
  id_ghosted = -1;
  volume = -1;
  viscosity = 1.;
  density = 1.;
  moisture_content = 1.;
  relative_permeability = 1.;
  pressure_adj = 1.;
  specific_moisture_capacity = 1.;
  specific_moisture_capacity_adj = 1.;
  saturated = 1;
  for (int i=0; i<3; i++) {
    centoid[i] = -1;
    permeability[i] = 1.;
  }
  neighbor_cells = NULL;
  
}

void GridCell::setIdLocal(int i) { id_local = i; }
void GridCell::setIdGhosted(int i) { id_ghosted = i; }
void GridCell::setVolume(double d) { volume = d; }
void GridCell::setPermX(double d) { permeability[0] = d; }
void GridCell::setPermY(double d) { permeability[1] = d; }
void GridCell::setPermZ(double d) { permeability[2] = d; }
void GridCell::setViscosity(double d) { viscosity = d; }
void GridCell::setDensity(double d) { density = d; }
void GridCell::setMoistureContent(double d) { moisture_content = d; }
void GridCell::updateMoistureContent() { 
  moisture_content_0 = moisture_content; 
}
void GridCell::updateRelativePermeability(double pressure) {
  double capillary_pressure = pressure-patm;
  GridCell::saturation.computeRelativePermeability(capillary_pressure,
                                         &moisture_content,
                                         &relative_permeability,
                                         &specific_moisture_capacity);
  double one_plus_beta_p = 1.+beta*capillary_pressure;
  double one_plus_beta_p_sq = one_plus_beta_p * one_plus_beta_p;
  if (capillary_pressure < 0)  {
    saturated = 0;
    // _adj = transformed
    pressure_adj = pressure/one_plus_beta_p;
  }
  else {
    saturated = 1;
    pressure_adj = pressure;
  }
  relative_permeability_adj = relative_permeability * one_plus_beta_p_sq;
  specific_moisture_capacity_adj = specific_moisture_capacity * one_plus_beta_p_sq;
}
void GridCell::setFlowAccumulationCoef(double d) { 
  flow_accumulation_coef = d; 
}

int GridCell::getIdLocal() { return id_local; }
int GridCell::getIdGhosted() { return id_ghosted; }
double GridCell::getVolume() { return volume; }
double GridCell::getPermX() { return permeability[0]; }
double GridCell::getPermY() { return permeability[1]; }
double GridCell::getPermZ() { return permeability[2]; }
double GridCell::getViscosity() { return viscosity; }
double GridCell::getDensity() { return density; }
double GridCell::getMoistureContent() { return moisture_content; }
double GridCell::getFlowAccumulationCoef() { 
  return flow_accumulation_coef; 
}

double GridCell::computeFlowAccumulationCoef(double one_over_dt) {
  double vol_over_dt = volume*one_over_dt;
  double temp = vol_over_dt;
  if (!saturated) temp += (moisture_content_0-moisture_content)*vol_over_dt;
  return temp;
}

double *GridCell::getPermPtr() { return &permeability[0]; }

void GridCell::printInfo() {
  printf("lid: %d gid: %d vol: %f\n",
         id_local,id_ghosted,volume);
}

void computeAccumulation(double dt) {
}

GridCell::~GridCell() {
  if (neighbor_cells) delete [] neighbor_cells;
}
