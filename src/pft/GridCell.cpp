#include "GridCell.h"

GridCell::GridCell() {
  
  id_local = -1;
  id_ghosted = -1;
  volume = -1;
  viscosity = 1.;
  density = 1.;
  moisture_content = 1.;
  relative_permeability = 1.;
  relative_permeability_t = 1.;
  pressure_0_t = 0.;
  pressure_n_t = 0.;
  specific_moisture_capacity = 1.;
  specific_moisture_capacity_t = 1.;
  viscosity = viscosity0;
  density = density0;
  moisture_content_0 = 1.;
  moisture_content_n = 1.;
  moisture_content = 1.;

  for (int i=0; i<3; i++) {
    centroid[i] = -1;
    permeability[i] = permeability0;
  }
  neighbor_cells = NULL;
  saturationfunction = new SaturationFunction();
  
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
void GridCell::setPressure(double d) { pressure_0_t = d; }
void GridCell::updateMoistureContent0() { 
  moisture_content_0 = moisture_content; 
}
void GridCell::updateMoistureContentN() { 
  moisture_content_n = moisture_content; 
}
void GridCell::updateRelativePermeability(double pressure) {

  double capillary_pressure = pressure-patm;
  saturationfunction->computeRelativePermeability(capillary_pressure,
                                         &moisture_content,
                                         &relative_permeability,
                                         &specific_moisture_capacity);
  relative_permeability_t = relative_permeability;
  specific_moisture_capacity_t = specific_moisture_capacity;
  saturationfunction->transformPressure(&capillary_pressure,
                                        &relative_permeability_t,
                                        &specific_moisture_capacity_t);
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
double GridCell::getPressure0_t() { return pressure_0_t; }
double GridCell::getPressureN_t() { return pressure_n_t; }
double GridCell::getMoistureContent0() { return moisture_content_0; }
double GridCell::getMoistureContentN() { return moisture_content_n; }
double GridCell::getSpecificMoistureCapacity_t() { return specific_moisture_capacity_t; }
double GridCell::getFlowAccumulationCoef() { 
  return flow_accumulation_coef; 
}

double GridCell::computeFlowAccumulationCoef(double one_over_dt) {

  return volume*one_over_dt;
  /*
  double temp = vol_over_dt;
  if (capillary_pressure < 
  if (capillary_pressure < 0.) 
    temp += (moisture_content_0-moisture_content)*vol_over_dt;

  return temp;*/
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
