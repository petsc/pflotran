#ifndef GRID_H_
#define GRID_H_

#include "include/petsc.h"
#include "include/petscvec.h"
#include "include/petscmat.h"
#include "include/petscis.h"

#include "BoundarySet.h"
#include "Globals.h"
#include "GridCell.h"
#include "GridVertex.h"
#include "Source.h"
#include "StructuredGrid.h"

class Grid {
  
public:

  Grid();
  Grid(int nx, int ny, int nz);
  virtual ~Grid();

  void nullifyArrays();
  void createStructured(int nx, int ny, int nz);
  void createDA(int ndof);
//  void computeConnectivity();
  void computeCellMapping();
  void computeVertexMapping();
  void setUpCells();
  void setUpVertices();
  void mapVerticesToCells();
  void printISLocalToGlobalMapping();
  void printAO();
  void addBoundarySet(BoundarySet *new_set);
//  void printConnectivity();
  void printCells();
  void printVertices();
//  void addBoundaryConnection(int is, int ie, int js, int je, int ks, int ke, 
//                            char *face, char *type, double scalar);
  void addSource(int is, int ie, int js, int je, int ks, int ke, 
                 char *type, double scalar);

  void getVectorNatural(Vec *v);
  void getVectorLocal(Vec *v);
  void getVectorGlobal(Vec *v);

  void globalToNatural(Vec global, Vec natural);
//  BoundaryConnection *getBoundaryConnections();
  Source *getSources();

  int getNumberOfCellsGlobal();
  int getNumberOfCellsGhosted();
  int getNumberOfCellsLocal();

  int getNumberOfVerticesGlobal();
  int getNumberOfVerticesGhosted();
  int getNumberOfVerticesLocal();

  int getNx();
  int getNy();
  int getNz();
  int getN();
  
  void getCorners(int *xs, int *ys, int *zs, int *nx, int *ny, int *nz);
  void getGhostCorners(int *xs, int *ys, int *zs, int *nx, int *ny, int *nz);

  double getDx(int i);
  double getDy(int j);
  double getDz(int k);
  double *getOriginPtr();
  double getRotationDegrees();
  Vec getGridCellMaterialIDs();
  Vec getGridCellActivities();
  void setGridSpacing(double *dx, double *dy, double *dz);
  void setGridSpacing(double dx, double dy, double dz);
  void setOrigin(double x, double y, double z);
  void setRotation(double r);
  void computeCoordinates();
  void addBoundarySet();
  BoundarySet *getBoundarySet(char *name);

  void convertLocalCellDataGtoN(int *data);
  void convertLocalCellDataGtoN(double *data);
  int getVertexIdsNaturalLocal(int **natural_ids);
  int getVertexCoordinatesNaturalLocal(double **coordinates, int direction);
  int *getCellMaterialIds();
  int *getCellIds();
  int *getCellIdsNatural();
  int *getCellVertexIds(int ivert);

  int num_cells_global;
  int num_cells_local;
  int num_cells_ghosted;
//  int num_connections;
  int *cell_mapping_local_to_ghosted;
  int *cell_mapping_ghosted_to_local;
  int *cell_mapping_ghosted_to_natural;
  GridCell *cells;

  int num_vertices_global;
  int num_vertices_local;
  int num_vertices_ghosted;
  int *vertex_mapping_local_to_ghosted;
  int *vertex_mapping_ghosted_to_local;
  int *vertex_mapping_ghosted_to_natural;
  GridVertex *vertices;
//  GridConnection *connections;

  BoundarySet *boundary_sets;
  
  IS inactive_is;

private:
  PetscErrorCode ierr;
  StructuredGrid *structuredGrid;

};

#endif /*GRID_H_*/
