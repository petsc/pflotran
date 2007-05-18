#ifndef STRUCTUREDGRID_H_
#define STRUCTUREDGRID_H_

#include <stdlib.h>

#include "include/petsc.h"
#include "include/petscda.h"
#include "include/petscvec.h"
#include "include/petscmat.h"

#include "BoundaryCondition.h"
#include "Globals.h"
#include "GridCell.h"
#include "GridConnection.h"
#include "Source.h"

class StructuredGrid {

public:

  StructuredGrid();
  StructuredGrid(int nx, int ny, int nz, int fdof, int tdof);
  StructuredGrid(int nx, int ny, int nz, double dx, double dy, double dz, 
                 int fdof, int tdof);
  virtual ~StructuredGrid();
  
  void setUpDA(double dx, double dy, double dz);
  void setUpConnectivity(int *num_connections, GridConnection **connections);
  void setUpMapping(int *num_nodes_local, int *num_nodes_ghosted,
                    int **mapping_local_to_ghosted,
                    int **mapping_ghosted_to_local);
  void mapBoundaryCondition(int is, int ie, int js, int je, int ks, int ke, 
                            char *face);
  void mapSource(int is, int ie, int js, int je, int ks, int ke);
  void setUpCells(int num_cells, GridCell *cells);
  void printDACorners();
  void printAO();
  
  
  void get1dofVectorLocal(Vec *v);
  void get1dofVectorGlobal(Vec *v);
  void get1dofVectorNatural(Vec *v);
  void getFdofMatrix(Mat *m, MatType mtype);
  void getFdofVectorGlobal(Vec *v);
  void getFdofVectorLocal(Vec *v);
  void getTdofVectorGlobal(Vec *v);
  void getTdofVectorLocal(Vec *v);
  void globalToNatural(Vec global, Vec natural);
  int getNx();
  int getNy();
  int getNz();
  int getN();
  double getDx(int i);
  double getDy(int j);
  double getDz(int k);
  
private:

  int nx, ny, nz, fdof, tdof;
  double *dx, *dy, *dz;
  DA da_1dof, da_fdof, da_tdof;
  int gnx, gny, gnz, lnx, lny, lnz;
  int gxs, gys, gzs, lxs, lys, lzs;
  int gxe, gye, gze, lxe, lye, lze;
  int lnxXny, gnxXny;
};

#endif /*STRUCTUREDGRID_H_*/
