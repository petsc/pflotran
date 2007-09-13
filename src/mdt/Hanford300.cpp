#include "Hanford300.h"

Hanford300::Hanford300(Grid **grid_) {

  river_polygon = NULL;
  ascii_grids = NULL;

  double mx = 1.;
  double my = 1.;
  double mz = 1.;

#if 1
  int nx = 135;
  int ny = 250;
  int nz = 60;
#if 0
  nx *= 2;
  ny *= 2;
  nz *= 2;
  mx = .5;
  my = mx;
  mz = .5;
#endif
#else
  int nx = 17;
  int ny = 32;
  int nz = 15;
  mx = 8.;
  my = mx;
  mz = 4.;
#endif
//  int nx = 1; //debug1
//  int ny = 1; //debug1
//  int nz = 15; //debug1
  int n = nx*ny*nz;

  double dx = 10.*mx;
  double dy = 10.*my;
  double dz = 1.*mz;// */

  *grid_ = new Grid(nx,ny,nz);
  Grid *grid = *grid_;
  grid->setGridSpacing(dx,dy,dz);

  grid->setOrigin(593618.9,114565.1,70.);
//  grid->setOrigin(594860.,114875.,100.); //debug1
  grid->setRotation(14.);

  grid->computeCoordinates();
//  grid->computeConnectivity();
  grid->computeCellMapping();
  grid->setUpCells();
  grid->computeVertexMapping();
  grid->setUpVertices();
  grid->mapVerticesToCells();

  river_polygon = new Polygon();
  river_polygon->createRiverEdgePolygon();

#if 0
  char ascii_filename[1024];
  strcpy(ascii_filename,"test.asc");
  AsciiGrid **ascii_grid = new AsciiGrid*[2];
  ascii_grid[0] = new AsciiGrid(ascii_filename);
  ascii_grid[0]->setMaterialId(1);
  ascii_grid[1] = new AsciiGrid("default",2,2,-1.e10,-1.e10,1.e10,-9999.,
                                1.e20,0); 
#else

  AsciiGrid::nasciigrids = 7;
  string *grid_filenames = new string[AsciiGrid::nasciigrids];
  grid_filenames[0].append("basalt_PNNL_grid_20m.asc");
  grid_filenames[1].append("u9PNNL_grid_20m.asc");
  grid_filenames[2].append("u8PNNL_grid_20m.asc");
  grid_filenames[3].append("u7PNNL_grid_20m.asc");
  grid_filenames[4].append("u6PNNL_grid_20m.asc");
  grid_filenames[5].append("u5PNNL_grid_20m.asc");
  grid_filenames[6].append("u1PNNL_grid_20m.asc");

  ascii_grids = new AsciiGrid*[AsciiGrid::nasciigrids];
  for (int i=0; i<AsciiGrid::nasciigrids; i++) {
    char filename[32];
    strcpy(filename,grid_filenames[i].c_str());
    ascii_grids[i] = new AsciiGrid(filename);
  }
  ascii_grids[0]->setMaterialId(10);
  ascii_grids[1]->setMaterialId(9);
  ascii_grids[2]->setMaterialId(8);
  ascii_grids[3]->setMaterialId(7);
  ascii_grids[4]->setMaterialId(6);
  ascii_grids[5]->setMaterialId(5);
  ascii_grids[6]->setMaterialId(1);

  int mod = grid->num_cells_ghosted/10;
  for (int i=0; i<grid->num_cells_ghosted; i++) {
    int material_id = 0;
    double x = grid->cells[i].getX();
    double y = grid->cells[i].getY();
    double z = grid->cells[i].getZ();
    for (int ilayer=0; ilayer<AsciiGrid::nasciigrids; ilayer++) {
      double zlayer = ascii_grids[ilayer]->computeElevationFromCoordinate(x,y);
      if (zlayer > ascii_grids[ilayer]->nodata && zlayer >= z) {
        material_id = ascii_grids[ilayer]->getMaterialId();
        break;
      }
    }
    if (material_id == 0) grid->cells[i].setActive(0);
    grid->cells[i].setMaterialId(material_id);
    if (!river_polygon->pointInPolygon(x,y)) grid->cells[i].setActive(0);
    if (i%mod == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%d of %d cells mapped with materials and activity.\n",
                  i,grid->num_cells_ghosted);
    }
  }

#endif


  computeRiverBoundary(grid);
  computeWestBoundary(grid);
  computeRechargeBoundary(grid);
//  computeSouthBoundary(grid);

  BoundarySet *river = grid->getBoundarySet("River");
  BoundarySet *west = grid->getBoundarySet("West");
  BoundarySet *recharge = grid->getBoundarySet("Recharge");

  Condition *new_condition = new Condition("river.bc");
  river->condition = new_condition;
  new_condition = new Condition("west.bc");
  west->condition = new_condition;
  new_condition = new Condition("recharge.bc");
  recharge->condition = new_condition;
  
  new_condition = NULL;

}

void Hanford300::computeRiverBoundary(Grid *grid) {

  double river_stage = 106.;

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

  BoundarySet *river = new BoundarySet("River");

  int istart = lxs-gxs;
  int jstart = lys-gys;
  int kstart = lzs-gzs;
  int iend = istart+lnx;
  int jend = jstart+lny;
  int kend = kstart+lnz;

  // loop over local cells
  for (int k=kstart; k<kend; k++) {
    for (int j=jstart; j<jend; j++) {
      // start from east side
      for (int i=iend-1; i>=istart; i--) {
        int ghosted_id = i+j*gnx+k*gnxXny;
        if (((i == iend-1 && gxe-lxe == 0) || 
             (i<iend-1 && !grid->cells[ghosted_id+1].getActive())) &&
            grid->cells[ghosted_id].getActive() && 
            grid->cells[ghosted_id].getZ() <= river_stage) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[ghosted_id].getHexFaceVertices(EAST,vertex_list);
          river->addConnection(new Connection(grid->cells[ghosted_id].getIdLocal(),vertex_list));
          // Look downwind in y and z to see if the bc cell is located at
          // the same distance x.  If not, we need to set the same bc along
          // the y and z faces to the next cell in x
          // y-direction
          if (j < jend-1+(gye-lye)) {
            for (int ijp1=iend-1; ijp1>=istart; ijp1--) {
              int ghosted_id_jp1 = ijp1+(j+1)*gnx+k*gnxXny;
              if (((ijp1 == iend-1 && gxe-lxe == 0) || 
                   (ijp1<iend-1 && !grid->cells[ghosted_id_jp1+1].getActive())) &&
                  grid->cells[ghosted_id_jp1].getActive() && 
                  grid->cells[ghosted_id_jp1].getZ() <= river_stage) {
                if (ijp1 != i) {
                  int iistart = i > ijp1 ? i : ijp1;
                  int iiend = i < ijp1 ? i : ijp1;
                  for (int ii=iistart; ii>iiend; ii--) {
                    if (i > ijp1) {
                      int id = ii+j*gnx+k*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(NORTH,vertex_list);
                      river->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                    else {
                      int id = ii+(j+1)*gnx+k*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(SOUTH,vertex_list);
                      river->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                  }
                }
                break;
              }
            }
          }
          if (k < kend-1+(gze-lze)) {
            for (int ikp1=iend-1; ikp1>=istart; ikp1--) {
              int ghosted_id_kp1 = ikp1+j*gnx+(k+1)*gnxXny;
              if (((ikp1 == iend-1 && gxe-lxe == 0) || 
                   (ikp1<iend-1 && !grid->cells[ghosted_id_kp1+1].getActive())) &&
                  grid->cells[ghosted_id_kp1].getActive() && 
                  grid->cells[ghosted_id_kp1].getZ() <= river_stage) {
                if (ikp1 != i) {
                  int iistart = i > ikp1 ? i : ikp1;
                  int iiend = i < ikp1 ? i : ikp1;
                  for (int ii=iistart; ii>iiend; ii--) {
                    if (i > ikp1) {
                      int id = ii+j*gnx+k*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(TOP,vertex_list);
                      river->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                    else {
                      int id = ii+j*gnx+(k+1)*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(BOTTOM,vertex_list);
                      river->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                  }
                }
                break;
              }
            }
          }
          break;
        }
      }
    }
  }
  grid->addBoundarySet(river);
  river = NULL;
}

void Hanford300::computeWestBoundary(Grid *grid) {

  double west_stage = 115.;

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

  BoundarySet *west = new BoundarySet("West");

  int istart = lxs-gxs;
  int jstart = lys-gys;
  int kstart = lzs-gzs;
  int iend = istart+lnx;
  int jend = jstart+lny;
  int kend = kstart+lnz;

  // loop over local cells
  for (int k=kstart; k<kend; k++) {
    for (int j=jstart; j<jend; j++) {
      for (int i=istart; i<iend; i++) {
        int ghosted_id = i+j*gnx+k*gnxXny;
        if (((i == 0 && lxs-gxs == 0) || 
             (i>istart && !grid->cells[ghosted_id-1].getActive())) &&
            grid->cells[ghosted_id].getActive() && 
            grid->cells[ghosted_id].getZ() <= west_stage) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[ghosted_id].getHexFaceVertices(WEST,vertex_list);
          west->addConnection(new Connection(grid->cells[ghosted_id].getIdLocal(),
                                             vertex_list));
          // Look downwind in y and z to see if the bc cell is located at
          // the same distance x.  If not, we need to set the same bc along
          // the y and z faces to the next cell in x
          // y-direction
          if (j < jend-1+(gye-lye)) {
            for (int ijp1=istart; ijp1<iend; ijp1++) {
              int ghosted_id_jp1 = ijp1+(j+1)*gnx+k*gnxXny;
              if (((ijp1 == 0 && lxs-gxs == 0) || 
                   (ijp1>istart && !grid->cells[ghosted_id_jp1-1].getActive())) &&
                  grid->cells[ghosted_id_jp1].getActive() && 
                  grid->cells[ghosted_id_jp1].getZ() <= west_stage) {
                if (ijp1 != i) {
                  int iistart = i < ijp1 ? i : ijp1;
                  int iiend = i > ijp1 ? i : ijp1;
                  for (int ii=iistart; ii<iiend; ii++) {
                    if (i < ijp1) {
                      int id = ii+j*gnx+k*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(NORTH,vertex_list);
                      west->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                    else {
                      int id = ii+(j+1)*gnx+k*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(SOUTH,vertex_list);
                      west->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                  }
                }
                break;
              }
            }
          }
          if (k < kend-1+(gze-lze)) {
            for (int ikp1=istart; ikp1<iend; ikp1++) {
              int ghosted_id_kp1 = ikp1+j*gnx+(k+1)*gnxXny;
              if (((ikp1 == 0 && lxs-gxs == 0) || 
                   (ikp1>istart && !grid->cells[ghosted_id_kp1-1].getActive())) &&
                  grid->cells[ghosted_id_kp1].getActive() && 
                  grid->cells[ghosted_id_kp1].getZ() <= west_stage) {
                if (ikp1 != i) {
                  int iistart = i < ikp1 ? i : ikp1;
                  int iiend = i > ikp1 ? i : ikp1;
                  for (int ii=iistart; ii<iiend; ii++) {
                    if (i < ikp1) {
                      int id = ii+j*gnx+k*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(TOP,vertex_list);
                      west->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                    else {
                      int id = ii+j*gnx+(k+1)*gnxXny;
                      grid->cells[ghosted_id].getHexFaceVertices(BOTTOM,vertex_list);
                      west->addConnection(new Connection(grid->cells[id].getIdLocal(),vertex_list));
                    }
                  }
                }
                break;
              }
            }
          }
          break;
        }
      }
    }
  }
  grid->addBoundarySet(west);
  west = NULL;
}

void Hanford300::computeRechargeBoundary(Grid *grid) {

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

  BoundarySet *recharge = new BoundarySet("Recharge");

  int istart = lxs-gxs;
  int jstart = lys-gys;
  int kstart = lzs-gzs;
  int iend = istart+lnx;
  int jend = jstart+lny;
  int kend = kstart+lnz;

  // loop over local cells
  for (int j=jstart; j<jend; j++) {
    for (int i=istart; i<iend; i++) {
      // start from top side
      for (int k=kend-1; k>=kstart; k--) {
        int ghosted_id = i+j*gnx+k*gnxXny;
        if (((k == kend-1 && gze-lze == 0) || 
             (k<kend-1 && !grid->cells[ghosted_id+gnxXny].getActive())) &&
            grid->cells[ghosted_id].getActive()) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[ghosted_id].getHexFaceVertices(TOP,vertex_list);
          recharge->addConnection(new Connection(grid->cells[ghosted_id].getIdLocal(),
                                                 vertex_list));
          break;
        }
      }
    }
  }
  grid->addBoundarySet(recharge);
  recharge = NULL;
}

void Hanford300::computeSouthBoundary(Grid *grid) {

  double south_stage = 115.;

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

  BoundarySet *south = new BoundarySet("South");

  int istart = lxs-gxs;
  int jstart = lys-gys;
  int kstart = lzs-gzs;
  int iend = istart+lnx;
  int jend = jstart+lny;
  int kend = kstart+lnz;

  // loop over local cells
  for (int k=kstart; k<kend; k++) {
    for (int i=istart; i<iend; i++) {
      for (int j=jstart; j<jend; j++) {
        int ghosted_id = i+j*gnx+k*gnxXny;
        if (((j == 0 && lys-gys == 0) || 
             (j>jstart && !grid->cells[ghosted_id-gnx].getActive())) &&
            grid->cells[ghosted_id].getActive() && 
            grid->cells[ghosted_id].getZ() <= south_stage) {
          int vertex_list[5] = {4,0,0,0,0};
          grid->cells[ghosted_id].getHexFaceVertices(SOUTH,vertex_list);
          south->addConnection(new Connection(grid->cells[ghosted_id].getIdLocal(),
                                              vertex_list));
          // Look downwind in y and z to see if the bc cell is located at
          // the same distance x.  If not, we need to set the same bc along
          // the y and z faces to the next cell in x
          // y-direction
          if (i < iend-1+(gxe-lxe)) {
            for (int jip1=jstart; jip1<jend; jip1++) {
              int ghosted_id_ip1 = (i+1)+jip1*gnx+k*gnxXny;
              if (((jip1 == 0 && lys-gys == 0) || 
                   (jip1>jstart && !grid->cells[ghosted_id_ip1-gnx].getActive())) &&
                  grid->cells[ghosted_id_ip1].getActive() && 
                  grid->cells[ghosted_id_ip1].getZ() <= south_stage) {
                if (jip1 != j) {
                  int jjstart = j < jip1 ? j : jip1;
                  int jjend = j > jip1 ? j : jip1;
                  for (int jj=jjstart; jj<jjend; jj++) {
                    if (j < jip1) {
                      int id = i+jj*gnx+k*gnxXny;
                      grid->cells[id].getHexFaceVertices(EAST,vertex_list);
                      south->addConnection(new Connection(grid->cells[id].getIdLocal(),
                                                          vertex_list));
                    }
                    else {
                      int id = (i+1)+jj*gnx+k*gnxXny;
                      grid->cells[id].getHexFaceVertices(WEST,vertex_list);
                      south->addConnection(new Connection(grid->cells[id].getIdLocal(),
                                                          vertex_list));
                    }
                  }
                }
                break;
              }
            }
          }
          if (k < kend-1+(gze-lze)) {
            for (int jkp1=istart; jkp1<iend; jkp1++) {
              int ghosted_id_kp1 = i+jkp1*gnx+(k+1)*gnxXny;
              if (((jkp1 == 0 && lys-gys == 0) || 
                   (jkp1>jstart && !grid->cells[ghosted_id_kp1-gnx].getActive())) &&
                  grid->cells[ghosted_id_kp1].getActive() && 
                  grid->cells[ghosted_id_kp1].getZ() <= south_stage) {
                if (jkp1 != j) {
                  int jjstart = j < jkp1 ? j : jkp1;
                  int jjend = j > jkp1 ? j : jkp1;
                  for (int jj=jjstart; jj<jjend; jj++) {
                    if (j < jkp1) {
                      int id = i+jj*gnx+k*gnxXny;
                      grid->cells[id].getHexFaceVertices(TOP,vertex_list);
                      south->addConnection(new Connection(grid->cells[id].getIdLocal(),
                                                          vertex_list));
                    }
                    else {
                      int id = i+jj*gnx+(k+1)*gnxXny;
                      grid->cells[id].getHexFaceVertices(BOTTOM,vertex_list);
                      south->addConnection(new Connection(grid->cells[id].getIdLocal(),
                                                          vertex_list));
                    }
                  }
                }
                break;
              }
            }
          }
          break;
        }
      }
    }
  }
  grid->addBoundarySet(south);
  south = NULL;
}

Hanford300::~Hanford300() {
  if (ascii_grids) {
    for (int i=0; i<AsciiGrid::nasciigrids; i++)
      delete ascii_grids[i];
    delete [] ascii_grids;
  }
  if (river_polygon) delete river_polygon;
}
