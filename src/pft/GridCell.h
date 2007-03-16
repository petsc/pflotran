#ifndef GRIDCELL_H_
#define GRIDCELL_H_

#include "include/petsc.h"

#include "Saturation.h"

class GridCell {
public:
  int num_cells;

  GridCell();
  virtual ~GridCell();
  void setIdLocal(int i);
  void setIdGhosted(int i);
  void setVolume(double d);
  void setPermX(double d);
  void setPermY(double d);
  void setPermZ(double d);
  void setViscosity(double d);
  void setDensity(double d);
  void setMoistureContent(double d);
  void updateMoistureContent();
  void updateRelativePermeability(double pressure);
  void setFlowAccumulationCoef(double d);
  
  int getIdLocal();
  int getIdGhosted();
  double getVolume();
  double getPermX();
  double getPermY();
  double getPermZ();
  double *getPermPtr();
  double getViscosity();
  double getDensity();
  double getMoistureContent();
  double getFlowAccumulationCoef();
  
  double computeFlowAccumulationCoef(double one_over_dt);
  
  void printInfo();
  
private:
  int id_local, id_ghosted;
  double volume;
  double centoid[3];
  GridCell *neighbor_cells;
  
// matrix variables
  double permeability[3];
  double relative_permeability;
  double relative_permeability_adj;

// fluid variables
  static Saturation saturation;
  double flow_accumulation_coef;
  double pressure_adj;
  double specific_moisture_capacity;
  double specific_moisture_capacity_adj;
  static const double beta;
  double one_plus_beta_h;
  double viscosity;
  double density;
  double moisture_content;
  double moisture_content_0;
  int saturated;

  // solute variables

};

#endif /*GRIDCELL_H_*/
