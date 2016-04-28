#include "TestCase.h"

#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif
#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif

TestCase::TestCase(Grid **grid_) {

  PetscReal mx = 1.;
  PetscReal my = 1.;
  PetscReal mz = 1.;

//#define GRID_3X3X3
#ifdef GRID_3X3X3

  PetscInt nx = 3;
  PetscInt ny = 3;
  PetscInt nz = 3;

  PetscReal dx = 1.;
  PetscReal dy = 1.;
  PetscReal dz = 1.;

#else

  PetscInt nx = 5;
  PetscInt ny = 4;
  PetscInt nz = 3;

  PetscReal dx[5] = {10.,11.,12.,13.,14.};
  PetscReal dy[4] = {13.,12.,11.,10.};
  PetscReal dz[3] = {15.,20.,25.};

#endif

  PetscInt n = nx*ny*nz;

  *grid_ = new Grid(nx,ny,nz);
  Grid *grid = *grid_;
  grid->setGridSpacing(dx,dy,dz);
  grid->setLocalGridSpacing();

  grid->setOrigin(0.,0.,0.);
  grid->computeCoordinates();
//  grid->computeConnectivity();
  grid->computeCellMapping();
  grid->setUpCells();
  grid->computeVertexMapping();
  grid->setUpVertices();
  grid->mapVerticesToCells();

#ifdef GRID_3X3X3
  for (PetscInt ia=0; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,1,grid);
  for (PetscInt ia=5; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,2,grid);
  for (PetscInt ia=11; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,3,grid);

  setMaterialIdBasedOnNaturalId(9,1,grid);
  setMaterialIdBasedOnNaturalId(12,2,grid);
  setMaterialIdBasedOnNaturalId(13,2,grid);
  setMaterialIdBasedOnNaturalId(15,2,grid);
  setMaterialIdBasedOnNaturalId(26,4,grid);
  
  setActiveBasedOnNaturalId(10,0,grid);
  setActiveBasedOnNaturalId(17,0,grid);
  setActiveBasedOnNaturalId(24,0,grid);
#else
  for (PetscInt ia=0; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,1,grid);
  for (PetscInt ia=9; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,2,grid);
  for (PetscInt ia=23; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,3,grid);

  setMaterialIdBasedOnNaturalId(11,1,grid);
  setMaterialIdBasedOnNaturalId(12,1,grid);
  setMaterialIdBasedOnNaturalId(21,1,grid);
  setMaterialIdBasedOnNaturalId(22,1,grid);
  setMaterialIdBasedOnNaturalId(27,1,grid);

  setMaterialIdBasedOnNaturalId(25,2,grid);
  setMaterialIdBasedOnNaturalId(26,2,grid);
  setMaterialIdBasedOnNaturalId(31,2,grid);
  setMaterialIdBasedOnNaturalId(32,2,grid);
  setMaterialIdBasedOnNaturalId(37,2,grid);
  setMaterialIdBasedOnNaturalId(40,2,grid);

  setMaterialIdBasedOnNaturalId(24,4,grid);
  setMaterialIdBasedOnNaturalId(43,4,grid);
  setMaterialIdBasedOnNaturalId(44,4,grid);
  setMaterialIdBasedOnNaturalId(48,4,grid);
  setMaterialIdBasedOnNaturalId(49,4,grid);

  setActiveBasedOnNaturalId(21,0,grid);
  setActiveBasedOnNaturalId(23,0,grid);
  setActiveBasedOnNaturalId(26,0,grid);
  setActiveBasedOnNaturalId(33,0,grid);
  setActiveBasedOnNaturalId(34,0,grid);
  setActiveBasedOnNaturalId(38,0,grid);
  setActiveBasedOnNaturalId(39,0,grid);
  setActiveBasedOnNaturalId(46,0,grid);
  setActiveBasedOnNaturalId(47,0,grid);
  setActiveBasedOnNaturalId(51,0,grid);
  setActiveBasedOnNaturalId(52,0,grid);
  setActiveBasedOnNaturalId(55,0,grid);

#endif    

//  grid->printCells();
//  grid->printVertices();

  flagGridCells(grid);

  computeWestBoundary(grid,1);
  computeEastBoundary(grid,1);
  computeSouthBoundary(grid,0);
  computeNorthBoundary(grid,0);
  computeBottomBoundary(grid,0);
  computeTopBoundary(grid,0);

  BoundarySet *west = grid->getBoundarySet("West");
  BoundarySet *east = grid->getBoundarySet("East");
  BoundarySet *south = grid->getBoundarySet("South");
  BoundarySet *north = grid->getBoundarySet("North");
  BoundarySet *bottom = grid->getBoundarySet("Bottom");
  BoundarySet *top = grid->getBoundarySet("Top");

  Condition *new_condition = NULL;
//  new_condition = new Condition("river.bc");
  east->condition = new_condition;

//  new_condition = new Condition("west.bc");
  west->condition = new_condition;
  south->condition = new_condition;
  north->condition = new_condition;

//  new_condition = new Condition("recharge.bc");
  top->condition = new_condition;
  
  new_condition = NULL;

}


void TestCase::computeWestBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *west = new BoundarySet("West");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & WEST_DIR_WEST_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(WEST,vertex_list);
        west->addConnection(new Connection(local_id,vertex_list,WEST));
      }
      if (complete) {
        if (grid->cells[i].flag & WEST_DIR_SOUTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,SOUTH));
        }
        if (grid->cells[i].flag & WEST_DIR_NORTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,NORTH));
        }
        if (grid->cells[i].flag & WEST_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & WEST_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(west);
  west = NULL;

}

void TestCase::computeEastBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *east = new BoundarySet("East");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & EAST_DIR_EAST_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(EAST,vertex_list);
        east->addConnection(new Connection(local_id,vertex_list,EAST));
      }
      if (complete) {
        if (grid->cells[i].flag & EAST_DIR_SOUTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,SOUTH));
        }
        if (grid->cells[i].flag & EAST_DIR_NORTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,NORTH));
        }
        if (grid->cells[i].flag & EAST_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & EAST_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(east);
  east = NULL;

}

void TestCase::computeSouthBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *south = new BoundarySet("South");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & SOUTH_DIR_SOUTH_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
        south->addConnection(new Connection(local_id,vertex_list,SOUTH));
      }
      if (complete) {
        if (grid->cells[i].flag & SOUTH_DIR_WEST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,WEST));
        }
        if (grid->cells[i].flag & SOUTH_DIR_EAST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,EAST));
        }
        if (grid->cells[i].flag & SOUTH_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & SOUTH_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(south);
  south = NULL;
}

void TestCase::computeNorthBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *north = new BoundarySet("North");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & NORTH_DIR_NORTH_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
        north->addConnection(new Connection(local_id,vertex_list,NORTH));
      }
      if (complete) {
        if (grid->cells[i].flag & NORTH_DIR_WEST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,WEST));
        }
        if (grid->cells[i].flag & NORTH_DIR_EAST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,EAST));
        }
        if (grid->cells[i].flag & NORTH_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & NORTH_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(north);
  north = NULL;

}

void TestCase::computeBottomBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *bottom = new BoundarySet("Bottom");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & BOTTOM_DIR_BOTTOM_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
        bottom->addConnection(new Connection(local_id,vertex_list,BOTTOM));
      }
      if (complete) {
        if (grid->cells[i].flag & BOTTOM_DIR_WEST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,WEST));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_EAST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,EAST));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_SOUTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,SOUTH));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_NORTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,NORTH));
        }
      }
    }
  }

  grid->addBoundarySet(bottom);
  bottom = NULL;

}

void TestCase::computeTopBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *top = new BoundarySet("Top");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & TOP_DIR_TOP_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(TOP,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
#if 0
      if (grid->cells[i].flag & TOP_DIR_WEST_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(WEST,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
      if (grid->cells[i].flag & TOP_DIR_EAST_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(EAST,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
      if (grid->cells[i].flag & TOP_DIR_SOUTH_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
      if (grid->cells[i].flag & TOP_DIR_NORTH_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
#endif
    }
  }

  grid->addBoundarySet(top);
  top = NULL;

}

void TestCase::flagGridCells(Grid *grid) {

  PetscReal top_stage = 1.e20;
  PetscReal bottom_stage = 1.e20;
  PetscReal north_stage = 115.;
  PetscReal south_stage = 115.;
  PetscReal west_stage = 115.;
  PetscReal east_stage = 106.;

  PetscInt nx = grid->getNx();
  PetscInt ny = grid->getNy();
  PetscInt nz = grid->getNz();

  PetscInt lnx,lny,lnz,lxs,lys,lzs,lxe,lye,lze;
  PetscInt gnx,gny,gnz,gxs,gys,gzs,gxe,gye,gze;

  grid->getCorners(&lxs,&lys,&lzs,&lnx,&lny,&lnz);
  grid->getGhostCorners(&gxs,&gys,&gzs,&gnx,&gny,&gnz);

  lxe = lxs+lnx;
  lye = lys+lny;
  lze = lzs+lnz;
  gxe = gxs+gnx;
  gye = gys+gny;
  gze = gzs+gnz;

  PetscInt gnxXny = gnx*gny;

  PetscInt istart = lxs-gxs;
  PetscInt jstart = lys-gys;
  PetscInt kstart = lzs-gzs;
  PetscInt iend = istart+lnx;
  PetscInt jend = jstart+lny;
  PetscInt kend = kstart+lnz;


  grid->zeroGridCellFlags();

  PetscInt *matrix = NULL;

  // west
  matrix = new PetscInt[ny*nz];
  for (PetscInt i=0; i<ny*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
        if (lys <= jglobal && jglobal < lye) {
          PetscInt j = jglobal-gys;
          PetscInt k = kglobal-gzs;

          PetscInt flag[3];
          flag[0] = 0;
          flag[1] = 0;
          flag[2] = 0;
          if (lxs > 0) grid->receiveFlag(flag,WEST);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt i=0; i<gnx; i++) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= west_stage) {
                grid->cells[ghosted_id].flag |= WEST_DIR_WEST_FACE;
                flag[0] = 1;
                matrix[jglobal+kglobal*ny] = i+gxs;
                break;
              }
            }
          }
          if (lxe < nx) grid->sendFlag(flag,EAST);
        }
      }
    }
  }
  PetscMPIInt nyXnz = ny*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nyXnz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt j=0; j<gny; j++) {
      PetscInt jglobal = j+gys;
      // y-direction
      if (j < gny-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+1+kglobal*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=idown; i<iup; i++)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= WEST_DIR_SOUTH_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=iup; i<idown; i++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= WEST_DIR_NORTH_FACE;     
          }
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+(kglobal+1)*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=idown; i<iup; i++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= WEST_DIR_BOTTOM_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=iup; i<idown; i++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= WEST_DIR_TOP_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // east
  matrix = new PetscInt[ny*nz];
  for (PetscInt i=0; i<ny*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
        if (lys <= jglobal && jglobal < lye) {
          PetscInt j = jglobal-gys;
          PetscInt k = kglobal-gzs;

          PetscInt flag[3] = {0,0,0};

          if (lxe < nx) grid->receiveFlag(flag,EAST);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt i=gnx-1; i>=0; i--) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= east_stage) {
                grid->cells[ghosted_id].flag |= EAST_DIR_EAST_FACE;
                flag[0] = 1;
                matrix[jglobal+kglobal*ny] = i+gxs;
                break;
              }
            }
          }
          if (lxs > 0) grid->sendFlag(flag,WEST);
        }
      }
    }
  }
  nyXnz = ny*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nyXnz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt j=0; j<gny; j++) {
      PetscInt jglobal = j+gys;
      // y-direction
      if (j < gny-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+1+kglobal*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=iup; i>idown; i--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= EAST_DIR_NORTH_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=idown; i>iup; i--)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= EAST_DIR_SOUTH_FACE;     
          }
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+(kglobal+1)*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=iup; i>idown; i--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= EAST_DIR_TOP_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=idown; i>iup; i--)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= EAST_DIR_BOTTOM_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // south
  matrix = new PetscInt[nx*nz];
  for (PetscInt i=0; i<nx*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt k = kglobal-gzs;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lys > 0) grid->receiveFlag(flag,SOUTH);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt j=0; j<gny; j++) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= south_stage) {
                grid->cells[ghosted_id].flag |= SOUTH_DIR_SOUTH_FACE;
                flag[0] = 1;
                matrix[iglobal+kglobal*nx] = j+gys;
                break;
              }
            }
          }
          if (lye < ny) grid->sendFlag(flag,NORTH);
        }
      }
    }
  }
  PetscMPIInt nxXnz = nx*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nxXnz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+1+kglobal*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jdown; j<jup; j++)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= SOUTH_DIR_WEST_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jup; j<jdown; j++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= SOUTH_DIR_EAST_FACE;     
          }                                         
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+(kglobal+1)*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jdown; j<jup; j++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= SOUTH_DIR_BOTTOM_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jup; j<jdown; j++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= SOUTH_DIR_TOP_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;
      
  // north
  matrix = new PetscInt[nx*nz];
  for (PetscInt i=0; i<nx*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt k = kglobal-gzs;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lye < ny) grid->receiveFlag(flag,NORTH);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt j=gny-1; j>=0; j--) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= north_stage) {
                grid->cells[ghosted_id].flag |= NORTH_DIR_NORTH_FACE;
                flag[0] = 1;
                matrix[iglobal+kglobal*nx] = j+gys;
                break;
              }
            }
          }
          if (lys > 0) grid->sendFlag(flag,SOUTH);
        }
      }
    }
  }
  nxXnz = nx*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nx*nz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+1+kglobal*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jup; j>jdown; j--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= NORTH_DIR_EAST_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jdown; j>jup; j--)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= NORTH_DIR_WEST_FACE;     
          }                                         
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+(kglobal+1)*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jup; j>jdown; j--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= NORTH_DIR_TOP_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jdown-1; j>=jup; j--)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= NORTH_DIR_BOTTOM_FACE;
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // bottom
  matrix = new PetscInt[nx*ny];
  for (PetscInt i=0; i<nx*ny; i++)
    matrix[i] = -1;

  for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
    if (lys <= jglobal && jglobal < lye) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt j = jglobal-gys;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lzs > 0) grid->receiveFlag(flag,BOTTOM);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt k=0; k<gnz; k++) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= bottom_stage) {
                grid->cells[ghosted_id].flag |= BOTTOM_DIR_BOTTOM_FACE;
                flag[0] = 1;
                matrix[iglobal+jglobal*nx] = k+gzs;
                break;
              }
            }
          }
          if (lze < nz) grid->sendFlag(flag,TOP);
        }
      }
    }
  }
  PetscMPIInt nxXny = nx*ny;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nx*ny,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);
 
  for (PetscInt j=0; j<gny; j++) {
    PetscInt jglobal = j+gys;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+1+jglobal*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kdown; k<kup; k++)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_WEST_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kup; k<kdown; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_EAST_FACE;     
          }                                         
        }
      }
      // y-direction
      if (j < gny-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+(jglobal+1)*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kdown; k<kup; k++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= BOTTOM_DIR_SOUTH_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kup; k<kdown; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_NORTH_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // top
  matrix = new PetscInt[nx*ny];
  for (PetscInt i=0; i<nx*ny; i++)
    matrix[i] = -1;

  for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
    if (lys <= jglobal && jglobal < lye) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt j = jglobal-gys;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lze < nz) grid->receiveFlag(flag,TOP);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt k=gnz-1; k>=0; k--) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= top_stage) {
                grid->cells[ghosted_id].flag |= TOP_DIR_TOP_FACE;
                flag[0] = 1;
                matrix[iglobal+jglobal*nx] = k+gzs;
                break;
              }
            }
          }
          if (lzs > 0) grid->sendFlag(flag,BOTTOM);
        }
      }
    }
  }
#if 0
  nxXny = nx*ny;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nx*ny,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);
 
  for (PetscInt j=0; j<gny; j++) {
    PetscInt jglobal = j+gys;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+1+jglobal*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kup; k>kdown; k--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= TOP_DIR_EAST_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kdown; k>kup; k--)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= TOP_DIR_WEST_FACE;     
          }                                         
        }
      }
      // y-direction
      if (j < gny-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+(jglobal+1)*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kdown; k<kup; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= TOP_DIR_NORTH_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kup; k<kdown; k++)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= TOP_DIR_SOUTH_FACE;     
          }                                         
        }
      }
    }
  }
#endif
  delete [] matrix;
  matrix = NULL;

}

void TestCase::setMaterialIdBasedOnNaturalId(PetscInt natural_id, PetscInt material_id,
                                             Grid *grid) {
  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) 
    if (grid->cells[i].getIdNatural() == natural_id) {
      grid->cells[i].setMaterialId(material_id);
      break;
//printf("%d %d %d %d\n",myrank,i,grid->cells[i].getIdNatural(),grid->cells[i].getMaterialId());
  }
}

void TestCase::setActiveBasedOnNaturalId(PetscInt natural_id, PetscInt active,
                                         Grid *grid) {
  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) 
    if (grid->cells[i].getIdNatural() == natural_id) {
      grid->cells[i].setActive(active);
      break;
    }
}

TestCase::~TestCase() {
}
