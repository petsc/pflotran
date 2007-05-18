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
  double getViscosity();
  double getDensity();
  double getMoistureContent();
  double getFlowAccumulationCoef();
  double getPressure0_t();
  double getPressureN_t();
  double getMoistureContent0();
  double getMoistureContentN();
  double getSpecificMoistureCapacity_t(); 
  
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
  double relative_permeability_t;

// fluid variables
  SaturationFunction *saturationfunction;
  double flow_accumulation_coef;
  double pressure_0_t;
  double pressure_n_t;
  double specific_moisture_capacity;
  double specific_moisture_capacity_t;
  double viscosity;
  double density;
  double moisture_content_0;
  double moisture_content_n;
  double moisture_content;

  double q[3];

  // solute variables

};

#endif /*GRIDCELL_H_*/
