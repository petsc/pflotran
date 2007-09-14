#include "TestCase.h"

#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif
#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif

TestCase::TestCase(Grid **grid_) {

  double mx = 1.;
  double my = 1.;
  double mz = 1.;

//#define GRID_3X3X3
#ifdef GRID_3X3X3

  int nx = 3;
  int ny = 3;
  int nz = 3;

  double dx = 1.;
  double dy = 1.;
  double dz = 1.;

#else

  int nx = 5;
  int ny = 4;
  int nz = 3;

  double dx[5] = {1.,2.,3.,2.,1.};
  double dy[4] = {1.,2.,3.,4.};
  double dz[3] = {1.,2.,1.};

#endif

  int n = nx*ny*nz;

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
  for (int ia=0; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,1,grid);
  for (int ia=5; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,2,grid);
  for (int ia=11; ia<grid->getN(); ia++)
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
  for (int ia=0; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,1,grid);
  for (int ia=9; ia<grid->getN(); ia++)
    setMaterialIdBasedOnNaturalId(ia,2,grid);
  for (int ia=23; ia<grid->getN(); ia++)
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


void TestCase::computeWestBoundary(Grid *grid, int complete) {

  BoundarySet *west = new BoundarySet("West");

  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    int local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & WEST_DIR_WEST_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(WEST,vertex_list);
        west->addConnection(new Connection(local_id,vertex_list));
      }
      if (complete) {
        if (grid->cells[i].flag & WEST_DIR_SOUTH_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & WEST_DIR_NORTH_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & WEST_DIR_BOTTOM_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & WEST_DIR_TOP_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list));
        }
      }
    }
  }

  grid->addBoundarySet(west);
  west = NULL;

}

void TestCase::computeEastBoundary(Grid *grid, int complete) {

  BoundarySet *east = new BoundarySet("East");

  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    int local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & EAST_DIR_EAST_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(EAST,vertex_list);
        east->addConnection(new Connection(local_id,vertex_list));
      }
      if (complete) {
        if (grid->cells[i].flag & EAST_DIR_SOUTH_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & EAST_DIR_NORTH_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & EAST_DIR_BOTTOM_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & EAST_DIR_TOP_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list));
        }
      }
    }
  }

  grid->addBoundarySet(east);
  east = NULL;

}

void TestCase::computeSouthBoundary(Grid *grid, int complete) {

  BoundarySet *south = new BoundarySet("South");

  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    int local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & SOUTH_DIR_SOUTH_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
        south->addConnection(new Connection(local_id,vertex_list));
      }
      if (complete) {
        if (grid->cells[i].flag & SOUTH_DIR_WEST_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & SOUTH_DIR_EAST_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & SOUTH_DIR_BOTTOM_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & SOUTH_DIR_TOP_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list));
        }
      }
    }
  }

  grid->addBoundarySet(south);
  south = NULL;
}

void TestCase::computeNorthBoundary(Grid *grid, int complete) {

  BoundarySet *north = new BoundarySet("North");

  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    int local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & NORTH_DIR_NORTH_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
        north->addConnection(new Connection(local_id,vertex_list));
      }
      if (complete) {
        if (grid->cells[i].flag & NORTH_DIR_WEST_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & NORTH_DIR_EAST_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & NORTH_DIR_BOTTOM_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & NORTH_DIR_TOP_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list));
        }
      }
    }
  }

  grid->addBoundarySet(north);
  north = NULL;

}

void TestCase::computeBottomBoundary(Grid *grid, int complete) {

  BoundarySet *bottom = new BoundarySet("Bottom");

  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    int local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & BOTTOM_DIR_BOTTOM_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
        bottom->addConnection(new Connection(local_id,vertex_list));
      }
      if (complete) {
        if (grid->cells[i].flag & BOTTOM_DIR_WEST_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_EAST_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_SOUTH_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_NORTH_FACE) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list));
        }
      }
    }
  }

  grid->addBoundarySet(bottom);
  bottom = NULL;

}

void TestCase::computeTopBoundary(Grid *grid, int complete) {

  BoundarySet *top = new BoundarySet("Top");

  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    int local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & TOP_DIR_TOP_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(TOP,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
#if 0
      if (grid->cells[i].flag & TOP_DIR_WEST_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(WEST,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
      if (grid->cells[i].flag & TOP_DIR_EAST_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(EAST,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
      if (grid->cells[i].flag & TOP_DIR_SOUTH_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list));
      }
      if (grid->cells[i].flag & TOP_DIR_NORTH_FACE) {
        int vertex_list[5] = {4,0,0,0,0};
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

  double top_stage = 1.e20;
  double bottom_stage = 1.e20;
  double north_stage = 115.;
  double south_stage = 115.;
  double west_stage = 115.;
  double east_stage = 106.;

  int nx = grid->getNx();
  int ny = grid->getNy();
  int nz = grid->getNz();

  int lnx,lny,lnz,lxs,lys,lzs,lxe,lye,lze;
  int gnx,gny,gnz,gxs,gys,gzs,gxe,gye,gze;

  grid->getCorners(&lxs,&lys,&lzs,&lnx,&lny,&lnz);
  grid->getGhostCorners(&gxs,&gys,&gzs,&gnx,&gny,&gnz);

  lxe = lxs+lnx;
  lye = lys+lny;
  lze = lzs+lnz;
  gxe = gxs+gnx;
  gye = gys+gny;
  gze = gzs+gnz;

  int gnxXny = gnx*gny;

  int istart = lxs-gxs;
  int jstart = lys-gys;
  int kstart = lzs-gzs;
  int iend = istart+lnx;
  int jend = jstart+lny;
  int kend = kstart+lnz;


  grid->zeroGridCellFlags();

  int *matrix = NULL;

  // west
  matrix = new int[ny*nz];
  for (int i=0; i<ny*nz; i++)
    matrix[i] = -1;

  for (int kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (int jglobal=0; jglobal<ny; jglobal++) {
        if (lys <= jglobal && jglobal < lye) {
          int j = jglobal-gys;
          int k = kglobal-gzs;

          int flag[3];
          flag[0] = 0;
          flag[1] = 0;
          flag[2] = 0;
          if (lxs > 0) grid->receiveFlag(flag,WEST);
          if (flag[0] == 0) {

          // loop over local cells
            for (int i=0; i<gnx; i++) {
              int ghosted_id = i+j*gnx+k*gnxXny;
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
  MPI_Allreduce(MPI_IN_PLACE,matrix,ny*nz,MPI_INTEGER,MPI_MAX,
                PETSC_COMM_WORLD);

  for (int k=0; k<gnz; k++) {
    int kglobal = k+gzs;
    for (int j=0; j<gny; j++) {
      int jglobal = j+gys;
      // y-direction
      if (j < gny-1) {
        int iup = matrix[jglobal+kglobal*ny];
        int idown = matrix[jglobal+1+kglobal*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (int i=idown; i<iup; i++)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= WEST_DIR_SOUTH_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (int i=iup; i<idown; i++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= WEST_DIR_NORTH_FACE;     
          }
        }
      }
      // z-direction
      if (k < gnz-1) {
        int iup = matrix[jglobal+kglobal*ny];
        int idown = matrix[jglobal+(kglobal+1)*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (int i=idown; i<iup; i++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= WEST_DIR_BOTTOM_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (int i=iup; i<idown; i++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= WEST_DIR_TOP_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // east
  matrix = new int[ny*nz];
  for (int i=0; i<ny*nz; i++)
    matrix[i] = -1;

  for (int kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (int jglobal=0; jglobal<ny; jglobal++) {
        if (lys <= jglobal && jglobal < lye) {
          int j = jglobal-gys;
          int k = kglobal-gzs;

          int flag[3] = {0,0,0};

          if (lxe < nx) grid->receiveFlag(flag,EAST);
          if (flag[0] == 0) {

          // loop over local cells
            for (int i=gnx-1; i>=0; i--) {
              int ghosted_id = i+j*gnx+k*gnxXny;
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
  MPI_Allreduce(MPI_IN_PLACE,matrix,ny*nz,MPI_INTEGER,MPI_MAX,
                PETSC_COMM_WORLD);

  for (int k=0; k<gnz; k++) {
    int kglobal = k+gzs;
    for (int j=0; j<gny; j++) {
      int jglobal = j+gys;
      // y-direction
      if (j < gny-1) {
        int iup = matrix[jglobal+kglobal*ny];
        int idown = matrix[jglobal+1+kglobal*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (int i=iup; i>idown; i--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= EAST_DIR_NORTH_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (int i=idown; i>iup; i--)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= EAST_DIR_SOUTH_FACE;     
          }
        }
      }
      // z-direction
      if (k < gnz-1) {
        int iup = matrix[jglobal+kglobal*ny];
        int idown = matrix[jglobal+(kglobal+1)*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (int i=iup; i>idown; i--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= EAST_DIR_TOP_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (int i=idown; i>iup; i--)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= EAST_DIR_BOTTOM_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // south
  matrix = new int[nx*nz];
  for (int i=0; i<nx*nz; i++)
    matrix[i] = -1;

  for (int kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (int iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          int k = kglobal-gzs;
          int i = iglobal-gxs;

          int flag[3] = {0,0,0};

          if (lys > 0) grid->receiveFlag(flag,SOUTH);
          if (flag[0] == 0) {

          // loop over local cells
            for (int j=0; j<gny; j++) {
              int ghosted_id = i+j*gnx+k*gnxXny;
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
  MPI_Allreduce(MPI_IN_PLACE,matrix,nx*nz,MPI_INTEGER,MPI_MAX,
                PETSC_COMM_WORLD);

  for (int k=0; k<gnz; k++) {
    int kglobal = k+gzs;
    for (int i=0; i<gnx; i++) {
      int iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        int jup = matrix[iglobal+kglobal*nx];
        int jdown = matrix[iglobal+1+kglobal*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (int j=jdown; j<jup; j++)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= SOUTH_DIR_WEST_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (int j=jup; j<jdown; j++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= SOUTH_DIR_EAST_FACE;     
          }                                         
        }
      }
      // z-direction
      if (k < gnz-1) {
        int jup = matrix[iglobal+kglobal*nx];
        int jdown = matrix[iglobal+(kglobal+1)*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (int j=jdown; j<jup; j++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= SOUTH_DIR_BOTTOM_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (int j=jup; j<jdown; j++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= SOUTH_DIR_TOP_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;
      
  // north
  matrix = new int[nx*nz];
  for (int i=0; i<nx*nz; i++)
    matrix[i] = -1;

  for (int kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (int iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          int k = kglobal-gzs;
          int i = iglobal-gxs;

          int flag[3] = {0,0,0};

          if (lye < ny) grid->receiveFlag(flag,NORTH);
          if (flag[0] == 0) {

          // loop over local cells
            for (int j=gny-1; j>=0; j--) {
              int ghosted_id = i+j*gnx+k*gnxXny;
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
  MPI_Allreduce(MPI_IN_PLACE,matrix,nx*nz,MPI_INTEGER,MPI_MAX,
                PETSC_COMM_WORLD);

  for (int k=0; k<gnz; k++) {
    int kglobal = k+gzs;
    for (int i=0; i<gnx; i++) {
      int iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        int jup = matrix[iglobal+kglobal*nx];
        int jdown = matrix[iglobal+1+kglobal*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (int j=jup; j>jdown; j--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= NORTH_DIR_EAST_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (int j=jdown; j>jup; j--)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= NORTH_DIR_WEST_FACE;     
          }                                         
        }
      }
      // z-direction
      if (k < gnz-1) {
        int jup = matrix[iglobal+kglobal*nx];
        int jdown = matrix[iglobal+(kglobal+1)*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (int j=jup; j>jdown; j--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= NORTH_DIR_TOP_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (int j=jdown-1; j>=jup; j--)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= NORTH_DIR_BOTTOM_FACE;
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // bottom
  matrix = new int[nx*ny];
  for (int i=0; i<nx*ny; i++)
    matrix[i] = -1;

  for (int jglobal=0; jglobal<ny; jglobal++) {
    if (lys <= jglobal && jglobal < lye) {
      for (int iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          int j = jglobal-gys;
          int i = iglobal-gxs;

          int flag[3] = {0,0,0};

          if (lzs > 0) grid->receiveFlag(flag,BOTTOM);
          if (flag[0] == 0) {

          // loop over local cells
            for (int k=0; k<gnz; k++) {
              int ghosted_id = i+j*gnx+k*gnxXny;
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
  MPI_Allreduce(MPI_IN_PLACE,matrix,nx*ny,MPI_INTEGER,MPI_MAX,
                PETSC_COMM_WORLD);
 
  for (int j=0; j<gny; j++) {
    int jglobal = j+gys;
    for (int i=0; i<gnx; i++) {
      int iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        int kup = matrix[iglobal+jglobal*nx];
        int kdown = matrix[iglobal+1+jglobal*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (int k=kdown; k<kup; k++)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_WEST_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (int k=kup; k<kdown; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_EAST_FACE;     
          }                                         
        }
      }
      // y-direction
      if (j < gny-1) {
        int kup = matrix[iglobal+jglobal*nx];
        int kdown = matrix[iglobal+(jglobal+1)*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (int k=kdown; k<kup; k++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= BOTTOM_DIR_SOUTH_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (int k=kup; k<kdown; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_NORTH_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // top
  matrix = new int[nx*ny];
  for (int i=0; i<nx*ny; i++)
    matrix[i] = -1;

  for (int jglobal=0; jglobal<ny; jglobal++) {
    if (lys <= jglobal && jglobal < lye) {
      for (int iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          int j = jglobal-gys;
          int i = iglobal-gxs;

          int flag[3] = {0,0,0};

          if (lze < nz) grid->receiveFlag(flag,TOP);
          if (flag[0] == 0) {

          // loop over local cells
            for (int k=gnz-1; k>=0; k--) {
              int ghosted_id = i+j*gnx+k*gnxXny;
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
  MPI_Allreduce(MPI_IN_PLACE,matrix,nx*ny,MPI_INTEGER,MPI_MAX,
                PETSC_COMM_WORLD);
 
  for (int j=0; j<gny; j++) {
    int jglobal = j+gys;
    for (int i=0; i<gnx; i++) {
      int iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        int kup = matrix[iglobal+jglobal*nx];
        int kdown = matrix[iglobal+1+jglobal*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (int k=kup; k>kdown; k--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= TOP_DIR_EAST_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (int k=kdown; k>kup; k--)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= TOP_DIR_WEST_FACE;     
          }                                         
        }
      }
      // y-direction
      if (j < gny-1) {
        int kup = matrix[iglobal+jglobal*nx];
        int kdown = matrix[iglobal+(jglobal+1)*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (int k=kdown; k<kup; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= TOP_DIR_NORTH_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (int k=kup; k<kdown; k++)
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

void TestCase::setMaterialIdBasedOnNaturalId(int natural_id, int material_id,
                                             Grid *grid) {
  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) 
    if (grid->cells[i].getIdNatural() == natural_id)//  {
      grid->cells[i].setMaterialId(material_id);
//printf("%d %d %d %d\n",myrank,i,grid->cells[i].getIdNatural(),grid->cells[i].getMaterialId());
//}
}

void TestCase::setActiveBasedOnNaturalId(int natural_id, int active,
                                         Grid *grid) {
  for (int i=0; i<grid->getNumberOfCellsGhosted(); i++) 
    if (grid->cells[i].getIdNatural() == natural_id) 
      grid->cells[i].setActive(active);
}

TestCase::~TestCase() {
}
