#ifndef STRUCTUREDGRID_H_
#define STRUCTUREDGRID_H_

#include <stdlib.h>

#include "include/petsc.h"
#include "include/petscda.h"
#include "include/petscvec.h"
#include "include/petscmat.h"

#include "BoundaryConnection.h"
#include "Globals.h"
#include "GridCell.h"
#include "GridVertex.h"
#include "Source.h"

#define PI 3.14159265 

class StructuredGrid {

public:

  StructuredGrid();
  StructuredGrid(int nx, int ny, int nz);
  virtual ~StructuredGrid();
  
  void createDA();
  void createDA(int ndof);
  void setGridSpacing(double *dx_, double *dy_, double *dz_);
  void setGridSpacing(double dx_, double dy_, double dz_);
  void setLocalGridSpacing();
  void setOrigin(double origx, double origy, double origz);
  void setRotation(double r);
//  void computeConnectivity(int *num_connections, GridConnection **connections);
  void computeCoordinates();
  void computeCellMapping(int *num_cells_local, int *num_cells_ghosted,
                          int **cell_mapping_local_to_ghosted,
                          int **cell_mapping_ghosted_to_local,
                          int **cell_mapping_ghosted_to_natural);
  void computeVertexMapping(int *num_vertices_local, int *num_vertices_ghosted,
                          int **vertex_mapping_local_to_ghosted,
                          int **vertex_mapping_ghosted_to_local,
                          int **vertex_mapping_ghosted_to_natural);
  void mapVerticesToCells(GridCell *cells, GridVertex *vertices);
//  void mapBoundaryConnection(int is, int ie, int js, int je, int ks, int ke, 
//                            char *face);
  void mapSource(int is, int ie, int js, int je, int ks, int ke);
  void setUpCells(int num_cells, GridCell *cells);
  void setUpVertices(int num_vertices, GridVertex *vertices);
  void printDACorners();
  void printAO();
  
  
  void getVectorLocal(Vec *v);
  void getVectorGlobal(Vec *v);
  void getVectorNatural(Vec *v);
  void globalToNatural(Vec global, Vec natural);
  int getNx();
  int getNy();
  int getNz();
  int getN();

  int getNeighboringProcessor(int direction);
  void sendFlag(int *flag, int direction);
  void receiveFlag(int *flag, int direction);

  void getCorners(int *xs, int *ys, int *zs, int *nx, int *ny, int *nz);
  void getGhostCorners(int *xs, int *ys, int *zs, int *nx, int *ny, int *nz);
  
  int getGnx();
  int getGny();
  int getGnz();
  int getGxs();
  int getGys();
  int getGzs();
  int getGxe();
  int getGye();
  int getGze();
 
  double getDx(int i);
  double getDy(int j);
  double getDz(int k);
  double *getOriginPtr();
  double getRotationDegrees();
  void computeGridSpacingFromCellCentroids(GridCell *cells);
  
  void convertLocalCellDataGtoN(double *data);
  int *getLocalCellVertexNaturalIDs(GridCell *cells, GridVertex *vertices);

private:

  int nx, ny, nz, fdof, tdof;
  double *dx, *dy, *dz;
  double *ldx, *ldy, *ldz;
  double *gdx, *gdy, *gdz;
  double local_origin[3];
  double global_origin[3];
  double rotationZdegrees, rotationZradians;
  DA da;
  int gnx, gny, gnz, lnx, lny, lnz;
  int gxs, gys, gzs, lxs, lys, lzs;
  int gxe, gye, gze, lxe, lye, lze;
  int lnxXny, gnxXny;
};

#endif /*STRUCTUREDGRID_H_*/
