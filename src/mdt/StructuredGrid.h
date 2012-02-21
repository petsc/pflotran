#ifndef STRUCTUREDGRID_H_
#define STRUCTUREDGRID_H_

#include <stdlib.h>

#include "petscsys.h"
#include "petscdmda.h"
#include "petscvec.h"
#include "petscmat.h"

#include "BoundaryConnection.h"
#include "Globals.h"
#include "GridCell.h"
#include "GridVertex.h"
#include "Source.h"

#define PI 3.14159265 

class StructuredGrid {

public:

  StructuredGrid();
  StructuredGrid(PetscInt nx, PetscInt ny, PetscInt nz);
  virtual ~StructuredGrid();
  
  void createDA();
  void createDA(PetscInt ndof);
  void setGridSpacing(PetscReal *dx_, PetscReal *dy_, PetscReal *dz_);
  void setGridSpacing(PetscReal dx_, PetscReal dy_, PetscReal dz_);
  void setLocalGridSpacing();
  void setOrigin(PetscReal origx, PetscReal origy, PetscReal origz);
  void setRotation(PetscReal r);
//  void computeConnectivity(PetscInt *num_connections, GridConnection **connections);
  void computeCoordinates();
  void computeCoordinate(PetscReal *x, PetscReal*y);
  void computeCellMapping(PetscInt *num_cells_local, PetscInt *num_cells_ghosted,
                          PetscInt **cell_mapping_local_to_ghosted,
                          PetscInt **cell_mapping_ghosted_to_local,
                          PetscInt **cell_mapping_ghosted_to_natural);
  void computeVertexMapping(PetscInt *num_vertices_local, PetscInt *num_vertices_ghosted,
                          PetscInt **vertex_mapping_local_to_ghosted,
                          PetscInt **vertex_mapping_ghosted_to_local,
                          PetscInt **vertex_mapping_ghosted_to_natural);
  void mapVerticesToCells(GridCell *cells, GridVertex *vertices);
//  void mapBoundaryConnection(PetscInt is, PetscInt ie, PetscInt js, 
//                             PetscInt je, PetscInt ks, PetscInt ke, 
//                            char *face);
  void mapSource(PetscInt is, PetscInt ie, PetscInt js, PetscInt je, 
                 PetscInt ks, PetscInt ke);
  void setUpCells(PetscInt num_cells, GridCell *cells);
  void setUpVertices(PetscInt num_vertices, GridVertex *vertices);
  void printDACorners();
  void printAO();
  
  
  void getVectorLocal(Vec *v);
  void getVectorGlobal(Vec *v);
  void getVectorNatural(Vec *v);
  void globalToNatural(Vec global, Vec natural);
  PetscInt getNx();
  PetscInt getNy();
  PetscInt getNz();
  PetscInt getN();

  PetscInt getNeighboringProcessor(PetscInt direction);
  void sendFlag(PetscInt *flag, PetscInt direction);
  void receiveFlag(PetscInt *flag, PetscInt direction);

  void getCorners(PetscInt *xs, PetscInt *ys, PetscInt *zs, 
                  PetscInt *nx, PetscInt *ny, PetscInt *nz);
  void getGhostCorners(PetscInt *xs, PetscInt *ys, PetscInt *zs, 
                       PetscInt *nx, PetscInt *ny, PetscInt *nz);
  
  PetscInt getGnx();
  PetscInt getGny();
  PetscInt getGnz();
  PetscInt getGxs();
  PetscInt getGys();
  PetscInt getGzs();
  PetscInt getGxe();
  PetscInt getGye();
  PetscInt getGze();
 
  PetscReal getDx(PetscInt i);
  PetscReal getDy(PetscInt j);
  PetscReal getDz(PetscInt k);
  PetscReal *getOriginPtr();
  PetscReal getRotationDegrees();
  void computeGridSpacingFromCellCentroids(GridCell *cells);
  
  void convertLocalCellDataGtoN(PetscReal *data);
  PetscInt *getLocalCellVertexNaturalIDs(GridCell *cells, GridVertex *vertices);

private:

  PetscInt nx, ny, nz, fdof, tdof;
  PetscReal *dx, *dy, *dz;
  PetscReal *ldx, *ldy, *ldz;
  PetscReal *gdx, *gdy, *gdz;
  PetscReal local_origin[3];
  PetscReal global_origin[3];
  PetscReal rotationZdegrees, rotationZradians;
  DM da;
  PetscInt gnx, gny, gnz, lnx, lny, lnz;
  PetscInt gxs, gys, gzs, lxs, lys, lzs;
  PetscInt gxe, gye, gze, lxe, lye, lze;
  PetscInt lnxXny, gnxXny;
};

#endif /*STRUCTUREDGRID_H_*/
