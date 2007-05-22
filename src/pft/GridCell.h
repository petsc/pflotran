#ifndef GRIDCELL_H_
#define GRIDCELL_H_

#include "include/petsc.h"

#include "SaturationFunction.h"

class GridCell {
public:
  int num_cells;

  GridCell();
  virtual ~GridCell();
  void setIdLocal(int i);
  void setIdGhosted(int i);
  void setCentroid(double x, double y, double z);
  void setX(double x);
  void setY(double y);
  void setZ(double z);
  void setVolume(double d);
  void setPermX(double d);
  void setPermY(double d);
  void setPermZ(double d);
  void setViscosity(double d);
  void setDensity(double d);
  void setMoistureContent(double d);
  void setPressure(double d);
  void updateMoistureContent0();
  void updateMoistureContentN(); 

  void setSaturationFunction(SaturationFunction *sat_func_);
  void updateRelativePermeability(double pressure);

  void setFlowAccumulationCoef(double d);
  void addFlux(double qDarcy, double *norm);
  void zeroFlux();
  double *getFlux();
  
  int getIdLocal();
  int getIdGhosted();
  double *getCentroidPtr();
  double getX();
  double getY();
  double getZ();
  double getVolume();
  double getPermX();
  double getPermY();
  double getPermZ();
  double *getPermPtr();
  double getRPerm();
  double getDRelPerm();
  double getViscosity();
  double getDensity();
  double getMoistureContent();
  double getFlowAccumulationCoef();
  double getPressure0();
  double getPressureN();
  double getMoistureContent0();
  double getMoistureContentN();
  double getSpecificMoistureCapacity(); 
  double getDerivSpecMoistCap();
  
  double computeFlowAccumulationCoef(double one_over_dt);
  
  void printInfo();
  
private:
  int id_local, id_ghosted;
  double volume;
  double centroid[3];
  GridCell *neighbor_cells;
  
// matrix variables
  double permeability[3];
  double relative_permeability;
  double deriv_rel_perm;
  double deriv_spec_moist_cap;

// fluid variables
  SaturationFunction *sat_func;
  double flow_accumulation_coef;
  double pressure_0;
  double pressure_n;
  double specific_moisture_capacity;
  double viscosity;
  double porosity;
  double density;
  double moisture_content_0;
  double moisture_content_n;
  double moisture_content;

  double q[3];

  // solute variables

};

#endif /*GRIDCELL_H_*/
