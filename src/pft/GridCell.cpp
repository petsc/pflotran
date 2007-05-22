#include "GridCell.h"

GridCell::GridCell() {
  
  id_local = -1;
  id_ghosted = -1;
  volume = -1;
  viscosity = 1.;
  density = 1.;
  moisture_content = 1.;
  relative_permeability = 1.;
  deriv_rel_perm = 1.;
  pressure_n = 0.;
  specific_moisture_capacity = 1.;
  viscosity = viscosity0;
  density = density0;
  porosity = porosity0;
  moisture_content_0 = 1.;
  moisture_content_n = 1.;
  moisture_content = 1.;

  for (int i=0; i<3; i++) {
    centroid[i] = -1;
    permeability[i] = permeability0;
  }
  neighbor_cells = NULL;
  sat_func = NULL;
  
}

void GridCell::setIdLocal(int i) { id_local = i; }
void GridCell::setIdGhosted(int i) { id_ghosted = i; }
void GridCell::setCentroid(double x, double y, double z) {
  centroid[0] = x;
  centroid[1] = y;
  centroid[2] = z;
}
void GridCell::setX(double x) { centroid[0] = x; }
void GridCell::setY(double y) { centroid[1] = y; }
void GridCell::setZ(double z) { centroid[2] = z; }
void GridCell::setVolume(double d) { volume = d; }
void GridCell::setPermX(double d) { permeability[0] = d; }
void GridCell::setPermY(double d) { permeability[1] = d; }
void GridCell::setPermZ(double d) { permeability[2] = d; }
void GridCell::setViscosity(double d) { viscosity = d; }
void GridCell::setDensity(double d) { density = d; }
void GridCell::setMoistureContent(double d) { moisture_content = d; }
//void GridCell::setPressure(double d) { pressure_0_t = d; }

void GridCell::updateMoistureContent0() { 
  moisture_content_0 = moisture_content; 
}
void GridCell::updateMoistureContentN() { 
  moisture_content_n = moisture_content; 
}

void GridCell::setSaturationFunction(SaturationFunction *sat_func_) {
  sat_func = sat_func_;
}

void GridCell::updateRelativePermeability(double pressure) {

  moisture_content_n = moisture_content;

  double capillary_pressure = pressure-patm;
  sat_func->computeRelativePermeability(capillary_pressure,
                                         &moisture_content,
                                         &relative_permeability,
                                         &specific_moisture_capacity,
                                         &deriv_rel_perm,
                                         &deriv_spec_moist_cap);
 // relative_permeability_t = relative_permeability;
//  specific_moisture_capacity_t = specific_moisture_capacity;
//  pressure_n_t = capillary_pressure;
//  sat_func->transformPressure(&pressure_n_t,
//                              &relative_permeability_t,
//                              &specific_moisture_capacity_t);
//  for (int i=0; i<3; i++)
//    permeability_t[i] = permeability[i]*relative_permeability_t;
}

void GridCell::setFlowAccumulationCoef(double d) { 
  flow_accumulation_coef = d; 
}

void GridCell::zeroFlux() {
  q[0] = 0.;
  q[1] = 0.;
  q[2] = 0.;
}

void GridCell::addFlux(double qDarcy, double *norm) {
  for (int i=0; i<3; i++)
    q[i] += qDarcy*norm[i];
}

double *GridCell::getFlux() { return &q[0]; }

int GridCell::getIdLocal() { return id_local; }
int GridCell::getIdGhosted() { return id_ghosted; }
double *GridCell::getCentroidPtr() { return &centroid[0]; }
double GridCell::getX() { return centroid[0]; }
double GridCell::getY() { return centroid[1]; }
double GridCell::getZ() { return centroid[2]; }
double GridCell::getVolume() { return volume; }
double GridCell::getPermX() { return permeability[0]; }
double GridCell::getPermY() { return permeability[1]; }
double GridCell::getPermZ() { return permeability[2]; }
double GridCell::getViscosity() { return viscosity; }
double GridCell::getDensity() { return density; }
double GridCell::getMoistureContent() { return moisture_content; }
double GridCell::getPressure0() { return pressure_0; }
double GridCell::getPressureN() { return pressure_n; }
double GridCell::getMoistureContent0() { return moisture_content_0; }
double GridCell::getMoistureContentN() { return moisture_content_n; }
double GridCell::getSpecificMoistureCapacity() { return specific_moisture_capacity; }
double GridCell::getDerivSpecMoistCap() { return deriv_spec_moist_cap; }
double GridCell::getFlowAccumulationCoef() { 
  return flow_accumulation_coef; 
}

double GridCell::computeFlowAccumulationCoef(double one_over_dt) {
  return volume*one_over_dt;
}

double *GridCell::getPermPtr() { return &permeability[0]; }
double GridCell::getRPerm() { return relative_permeability; }
double GridCell::getDRelPerm() { return deriv_rel_perm; }

void GridCell::printInfo() {
  printf("lid: %d gid: %d vol: %f\n",
         id_local,id_ghosted,volume);
}

void computeAccumulation(double dt) {
}

GridCell::~GridCell() {
  if (neighbor_cells) delete [] neighbor_cells;
}
