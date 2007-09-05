#include "StructuredGrid.h"

StructuredGrid::StructuredGrid() {
}

StructuredGrid::StructuredGrid(int nx_, int ny_, int nz_) {
  nx = nx_; ny = ny_; nz = nz_;
  global_origin[0] = 0.; global_origin[1] = 0.; global_origin[2] = 0.;
  local_origin[0] = -999.; local_origin[1] = -999.; local_origin[2] = -999.;
  rotationZdegrees = 0.; rotationZradians = 0.;
}

void StructuredGrid::createDA() {

  PetscErrorCode ierr;

  int commsize;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);
  int pnx = PETSC_DECIDE;
  int pnz = PETSC_DECIDE;
  if (commsize > 1) pnx = 2; 
  if (commsize > 3) pnz = 2; 

  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
                    nx,ny,nz,pnx,PETSC_DECIDE,pnz,
                    1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&da);

  ierr = DAGetCorners(da,&lxs,&lys,&lzs,&lnx,&lny,&lnz);
  lxe = lxs+lnx;
  lye = lys+lny;
  lze = lzs+lnz;
  lnxXny = lnx*lny;
  
  ierr = DAGetGhostCorners(da,&gxs,&gys,&gzs,&gnx,&gny,&gnz);
  gxe = gxs+gnx;
  gye = gys+gny;
  gze = gzs+gnz;
  gnxXny = gnx*gny;
  
  printf(" local: %d %d %d %d %d %d\n",lxs,lxe,lys,lye,lzs,lze);
  printf("global: %d %d %d %d %d %d\n",gxs,gxe,gys,gye,gzs,gze);

}

void StructuredGrid::setGridSpacing(double *dx_, double *dy_, double *dz_) {

  dx = new double[nx];
  for (int i=0; i<nx; i++)
    dx[i] = dx_[i];
  dy = new double[ny];
  for (int j=0; j<ny; j++)
    dy[j] = dy_[j];
  dz = new double[nz];
  for (int k=0; k<nz; k++)
    dz[k] = dz_[k];

}

void StructuredGrid::setGridSpacing(double dx_, double dy_, double dz_) {

  dx = new double[nx];
  for (int i=0; i<nx; i++)
    dx[i] = dx_;
  dy = new double[ny];
  for (int j=0; j<ny; j++)
    dy[j] = dy_;
  dz = new double[nz];
  for (int k=0; k<nz; k++)
    dz[k] = dz_;

  setLocalGridSpacing();

}

void StructuredGrid::setLocalGridSpacing() {

  if (!dx || !dy || !dz) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Global grid spacing must be set before local.\n");
    PetscFinalize();
    exit(0);
  }

  gdx = new double[gnx];
  for (int i=gxs; i<gxe; i++)
    gdx[i-gxs] = dx[i];
  gdy = new double[gny];
  for (int j=gys; j<gye; j++)
    gdy[j-gys] = dy[j];
  gdz = new double[gnz];
  for (int k=gzs; k<gze; k++)
    gdz[k-gzs] = dz[k];

  ldx = new double[lnx];
  for (int i=lxs; i<lxe; i++)
    ldx[i-lxs] = dx[i];
  ldy = new double[lny];
  for (int j=lys; j<lye; j++)
    ldy[j-lys] = dy[j];
  ldz = new double[lnz];
  for (int k=lzs; k<lze; k++)
    ldz[k-lzs] = dz[k];

}

void StructuredGrid::setOrigin(double origx, double origy, double origz) {
  
  if (!dx || !dy || !dz) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Global grid spacing must be set before origin.\n");
    PetscFinalize();
    exit(0);
  }

  global_origin[0] = origx;
  global_origin[1] = origy;
  global_origin[2] = origz;

  double sumx = 0.;

  for (int i=0; i<lxs; i++)
    sumx += dx[i];

  double sumy = 0.;
  for (int j=0; j<lys; j++)
    sumy += dy[j];

  double sumz = 0.;
  for (int k=0; k<lzs; k++)
    sumz += dz[k];

  local_origin[0] = global_origin[0] + sumx*cos(rotationZradians) - 
                                       sumy*sin(rotationZradians);
  local_origin[1] = global_origin[1] + sumx*sin(rotationZradians) + 
                                       sumy*cos(rotationZradians);
  local_origin[2] = global_origin[2] + sumz;

  /*
  delete(dx);
  dx = NULL;
  delete(dy);
  dy = NULL;
  delete(dz);
  dz = NULL;
  */
}

void StructuredGrid::setRotation(double r) {
  rotationZdegrees = r;
  rotationZradians = r/180.*PI; // convert from degrees to radians
}

#if 0
void StructuredGrid::computeConnectivity(int *num_connections, 
                                         GridConnection **connections) {


// all connectivity is in the local ghosted context
  *num_connections = (gnx-1)*gny*gnz+gnx*(gny-1)*gnz+gnx*gny*(gnz-1);

  *connections = new GridConnection[*num_connections];
  GridConnection *conn = *connections;
  // conn[count].xx is the same as (*connections)[count].xx

  double vec[3];
  
// x-direction
  vec[1] = 0.; vec[2] = 0.;
  int count = 0;
  for (int k=0; k<gnz; k++) {
    for (int j=0; j<gny; j++) {
      for (int i=0; i<gnx-1; i++) {
        int id = i+j*gnx+k*gnxXny;
        conn[count].setIdUpwind(id);
        conn[count].setIdDownwind(id+1);
        vec[0] = 0.5*gdx[i]; 
        conn[count].setDistanceUpwind(vec);
        vec[0] = 0.5*gdx[i+1]; 
        conn[count].setDistanceDownwind(vec);
        vec[0] = 1.;
        conn[count].setNormal(vec);
        conn[count].setArea(gdy[j]*gdz[k]);
        count++;
      }
    }
  }  
// y-direction
  vec[0] = 0.; vec[2] = 0.;
  for (int k=0; k<gnz; k++) {
    for (int i=0; i<gnx; i++) {
      for (int j=0; j<gny-1; j++) {
        int id = i+j*gnx+k*gnxXny;
        conn[count].setIdUpwind(id);
        conn[count].setIdDownwind(id+gny);
        vec[1] = 0.5*gdy[j];
        conn[count].setDistanceUpwind(vec);
        vec[1] = 0.5*gdy[j+1];
        conn[count].setDistanceDownwind(vec);
        vec[1] = 1.;
        conn[count].setNormal(vec);
        conn[count].setArea(gdx[i]*gdz[k]);
        count++;
      }
    }
  } 
// z-direction
  vec[0] = 0.; vec[1] = 0.; 
  for (int j=0; j<gny; j++) {
    for (int i=0; i<gnx; i++) {
      for (int k=0; k<gnz-1; k++) {
        int id = i+j*gnx+k*gnxXny;
        conn[count].setIdUpwind(id);
        conn[count].setIdDownwind(id+gnxXny);
        vec[2] = 0.5*gdz[k];
        conn[count].setDistanceUpwind(vec);
        vec[2] = 0.5*gdz[k+1];
        conn[count].setDistanceDownwind(vec);
        vec[2] = 1.;
        conn[count].setNormal(vec);
        conn[count].setArea(gdx[i]*gdy[j]);
        count++;
      }
    }
  } 
}
#endif

void StructuredGrid::computeCoordinates() {

  if (fabs(local_origin[0] - -999.) < 1.e-40 &&
      fabs(local_origin[1] - -999.) < 1.e-40 &&
      fabs(local_origin[2] - -999.) < 1.e-40) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Origin must be set before computing coordinates.\n");
    PetscFinalize();
    exit(0);
  }

  Vec v;
  PetscScalar *v_ptr = NULL;
  VecCreate(PETSC_COMM_WORLD,&v);
  VecSetSizes(v,lnx*lny*lnz*3,PETSC_DECIDE);
  VecSetType(v,VECMPI);
  VecGetArray(v,&v_ptr);

  double sumz = 0.5*ldz[0];     
  for (int k=1; k<lnz; k++) {   
    double sumy = 0.5*ldy[0];
    for (int j=1; j<lny; j++) {
      double sumx = 0.5*ldx[0];
      for (int i=1; i<lnx; i++) {
        int local_id = i+j*lnx+k*lnxXny;
        v_ptr[local_id*3] =   local_origin[0] + sumx*cos(rotationZradians) - 
                                                sumy*sin(rotationZradians);
        v_ptr[local_id*3+1] = local_origin[1] + sumx*sin(rotationZradians) + 
                                                sumy*cos(rotationZradians);
        v_ptr[local_id*3+2] = local_origin[2] + sumz;
        if (i < lnx) sumx += 0.5*(ldx[i-1]+ldx[i]);
      }
      if (j < lny) sumy += 0.5*(ldy[j-1]+ldy[j]);
    }
    if (k < lnz) sumz += 0.5*(ldz[k-1]+ldz[k]);
  }
  
  VecRestoreArray(v,&v_ptr);

  DASetCoordinates(da,v);

}

#if 0
void StructuredGrid::mapBoundaryConnection(int istart, int iend, int jstart,
                                          int jend, int kstart, int kend,
                                          char *face) {

  // ensure monotonic increasing delineation of region
  if (istart > iend) {
    int temp = istart;
    istart =  iend;
    iend = temp;
  }
  if (jstart > jend) {
    int temp = jstart;
    jstart =  jend;
    jend = temp;
  }
  if (kstart > kend) {
    int temp = kstart;
    kstart =  kend;
    kend = temp;
  }

// decrement 'Xstarts' to zero-based numbering
  istart--;
  jstart--;
  kstart--;

//  printf("1 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);
// all connectivity is in the local non-ghosted context
// istart, iend, ..., kend are in global DA coordinates
  // check wether region is local to processor; if not, return;
  if (!((istart >= lxs && istart < lxe) &&
        (jstart >= lys && jstart < lye) &&
        (kstart >= lzs && kstart < lze))) return; 

//  printf("2 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);
  // clip region to portion on local processor
  istart = istart > lxs ? istart : lxs;
  iend = iend < lxe ? iend : lxe;
  jstart = jstart > lys ? jstart : lys;
  jend = jend < lye ? jend : lye;
  kstart = kstart > lzs ? kstart : lzs;
  kend = kend < lze ? kend : lze;

  int num_connections = abs(iend-istart)*abs(jend-jstart)*abs(kend-kstart);

//  printf("3 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);

  int count = 0;
  double dist[3];
  double center[3];
  double normal[3];
//  printf("4 %d %d %d %d %d %d\n",lxs,lxe,lys,lye,lzs,lze);
  for (int k=kstart-lzs; k<kend-lzs; k++) {
    for (int j=jstart-lys; j<jend-lys; j++) {
      for (int i=istart-lxs; i<iend-lxs; i++) {
        int ighosted = i+lxs-gxs;
        int jghosted = j+lys-gys;
        int kghosted = k+lzs-gzs;
        BoundaryConnection *newbc = new BoundaryConnection();
        newbc->setId(i+j*lnx+k*lnxXny);
        dist[0] = 0.; dist[1] = 0.; dist[2] = 0.; 
        normal[0] = 0.; normal[1] = 0.; normal[2] = 0.; 
        // init at center of cell and adjust single coordinate below
        center[0] = 0.5*gdx[ighosted]; center[1] = 0.5*gdy[jghosted]; 
        center[2] = 0.5*gdz[kghosted];
        if (!strcmp(face,"west") || !strcmp(face,"east")) {
          newbc->setArea(gdy[jghosted]*gdz[kghosted]);
          if (!strcmp(face,"west")) {
            normal[0] = -1.;
            dist[0] = -0.5*gdx[ighosted];
          }
          else if (!strcmp(face,"east")) {
            normal[0] = 1.;
            dist[0] = 0.5*gdx[ighosted];
          }
        }
        else if (!strcmp(face,"north") || !strcmp(face,"south")) {
          newbc->setArea(gdx[ighosted]*gdz[kghosted]);
          if (!strcmp(face,"north")) {
            normal[1] = -1.;
            dist[1] = -0.5*gdy[jghosted];
          }
          else if (!strcmp(face,"south")) {
            normal[1] = 1.;
            dist[1] = 0.5*gdy[jghosted];
          }
        }
        else if (!strcmp(face,"top") || !strcmp(face,"bottom")) {
          newbc->setArea(gdx[ighosted]*gdy[jghosted]);
          if (!strcmp(face,"top")) {
            normal[2] = -1.;
            dist[2] = -0.5*gdz[kghosted];
          }
          else if (!strcmp(face,"bottom")) {
            normal[2] = 1.;
            dist[2] = 0.5*gdz[kghosted];
          }
        }
        newbc->setDistance(dist);
        newbc->setNormal(normal);
        newbc = NULL;
        count++;
      }
    }
  }

#include "Globals.h"
  if (count != num_connections) {
    printf("Error setting up boundary connections on Processors[%d] %d %d\n",
           myrank,count,num_connections);
  }
}
#endif

void StructuredGrid::computeGridSpacingFromCellCentroids(GridCell *cells) {

  if (gnx > 1) {
    for (int i=0; i<gnx; i++) {
      int icell = i;
      double dx_up;
      double dx_dn;
      if (i > 0) dx_up = cells[icell].getX()-cells[icell-1].getX();
      else dx_up = cells[1].getX()-cells[0].getX();
      if (i < gnx-1) dx_dn = cells[icell+1].getX()-cells[icell].getX();
      else dx_dn = cells[gnx-1].getX()-cells[gnx-2].getX();
      dx[i] = 0.5*(dx_up+dx_dn);
    }
  }
  else dx[0] = 1.;
  if (gny > 1) {
    for (int j=0; j<gny; j++) {
      int icell = j*gnx;
      double dy_up;
      double dy_dn;
      if (j > 0) dy_up = cells[icell].getY()-cells[icell-gnx].getY();
      else dy_up = cells[gnx].getY()-cells[0].getY();
      if (j < gny-1) dy_dn = cells[icell+gnx].getY()-cells[icell].getY();
      else dy_dn = cells[(gny-1)*gnx].getY()-cells[(gny-2)*gnx].getY();
      dy[j] = 0.5*(dy_up+dy_dn);
    }
  }
  else dy[0] = 1.;
  if (gnz > 1) {
    for (int k=0; k<gnz; k++) {
      int icell = k*gnxXny;
      double dz_up;
      double dz_dn;
      if (k > 0) dz_up = cells[icell].getZ()-cells[icell-gnxXny].getZ();
      else dz_up = cells[gnxXny].getZ()-cells[0].getZ();
      if (k < gnz-1) dz_dn = cells[icell+gnxXny].getZ()-cells[icell].getZ();
      else dz_dn = cells[(gnz-1)*gnxXny].getZ()-cells[(gnz-2)*gnxXny].getZ();
      dz[k] = 0.5*(dz_up+dz_dn);
    }
  }
  else dz[0] = 1.;

}

void StructuredGrid::mapSource(int istart, int iend, int jstart,
                               int jend, int kstart, int kend) {

  // ensure monotonic increasing delineation of region
  if (istart > iend) {
    int temp = istart;
    istart =  iend;
    iend = temp;
  }
  if (jstart > jend) {
    int temp = jstart;
    jstart =  jend;
    jend = temp;
  }
  if (kstart > kend) {
    int temp = kstart;
    kstart =  kend;
    kend = temp;
  }

// decrement 'Xstarts' to zero-based numbering
  istart--;
  jstart--;
  kstart--;

//  printf("1 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);
// all connectivity is in the local non-ghosted context
// istart, iend, ..., kend are in global DA coordinates
  // check wether region is local to processor; if not, return;
  if (!((istart >= lxs && istart < lxe) &&
        (jstart >= lys && jstart < lye) &&
        (kstart >= lzs && kstart < lze))) return; 

//  printf("2 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);
  // clip region to portion on local processor
  istart = istart > lxs ? istart : lxs;
  iend = iend < lxe ? iend : lxe;
  jstart = jstart > lys ? jstart : lys;
  jend = jend < lye ? jend : lye;
  kstart = kstart > lzs ? kstart : lzs;
  kend = kend < lze ? kend : lze;

  int num_connections = abs(iend-istart)*abs(jend-jstart)*abs(kend-kstart);

//  printf("3 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);

  int count = 0;
//  printf("4 %d %d %d %d %d %d\n",lxs,lxe,lys,lye,lzs,lze);
  for (int k=kstart-lzs; k<kend-lzs; k++) {
    for (int j=jstart-lys; j<jend-lys; j++) {
      for (int i=istart-lxs; i<iend-lxs; i++) {
        Source *newsrc = new Source();
        newsrc->setId(i+j*lnx+k*lnxXny);
        newsrc = NULL;
        count++;
      }
    }
  }

#include "Globals.h"
  if (count != num_connections) {
    printf("Error setting up boundary connections on Processors[%d] %d %d\n",
           myrank,count,num_connections);
  }
}

void StructuredGrid::setUpCells(int num_cells, GridCell *cells) {
  double z = 0.5*gdz[0];
  for (int k=0; k<gnz; k++) {
    double y = 0.5*gdy[0];
    for (int j=0; j<gny; j++) {
      double x = 0.5*gdx[0];
      for (int i=0; i<gnx; i++) {
        int id = i+j*gnx+k*gnxXny;
        cells[id].setVolume(gdx[i]*gdy[j]*gdz[k]);
        double x_rot = local_origin[0] + (x-(lxs-gxs>0?gdx[0]:0.))*cos(rotationZradians) - 
                                         (y-(lys-gys>0?gdy[0]:0.))*sin(rotationZradians);
        double y_rot = local_origin[1] + (x-(lxs-gxs>0?gdx[0]:0.))*sin(rotationZradians) + 
                                         (y-(lys-gys>0?gdy[0]:0.))*cos(rotationZradians);
        double z_adj = local_origin[2] + (z-(lzs-gzs>0?gdz[0]:0.));
        cells[id].setCentroid(x_rot,y_rot,z_adj);
        if (i < gnx-1) x += 0.5*(gdx[i]+gdx[i+1]);
      }
      if (j < gny-1) y += 0.5*(gdy[j]+gdy[j+1]);
    }
    if (k < gnz-1) z += 0.5*(gdz[k]+gdz[k+1]);
  }
}

void StructuredGrid::setUpVertices(int num_vertices, GridVertex *vertices) {
  int gnxp1 = gnx+1;
  int gnyp1 = gny+1;
  int gnzp1 = gnz+1;

  double z = 0.;
  for (int k=0; k<gnzp1; k++) {
    double y = 0.;
    for (int j=0; j<gnyp1; j++) {
      double x = 0.;
      for (int i=0; i<gnxp1; i++) {
        int id = i+j*gnxp1+k*gnxp1*gnyp1;
        double x_rot = local_origin[0] + (x-(lxs-gxs>0?gdx[0]:0.))*cos(rotationZradians) - 
                                         (y-(lys-gys>0?gdy[0]:0.))*sin(rotationZradians);
        double y_rot = local_origin[1] + (x-(lxs-gxs>0?gdx[0]:0.))*sin(rotationZradians) + 
                                         (y-(lys-gys>0?gdy[0]:0.))*cos(rotationZradians);
        double z_adj = local_origin[2] + (z-(lzs-gzs>0?gdz[0]:0.));
        vertices[id].setX(x_rot);
        vertices[id].setY(y_rot);
        vertices[id].setZ(z_adj);
        if (i < gnx) x += gdx[i];
      }
      if (j < gny) y += gdy[j];
    }
    if (k < gnz) z += gdz[k];
  }
}

void StructuredGrid::mapVerticesToCells(GridCell *cells, GridVertex *vertices) {
  int gnxp1 = gnx+1;
  int gnyp1 = gny+1;
  int gnzp1 = gnz+1;
  for (int k=0; k<gnz; k++) {
    for (int j=0; j<gny; j++) {
      for (int i=0; i<gnx; i++) {
        int cell_id = i+j*gnx+k*gnxXny;
        int id[8];
        id[0] = i+j*gnxp1+k*gnxp1*gnyp1;
        id[1] = id[0]+1;
        id[2] = id[0]+gnxp1+1;
        id[3] = id[0]+gnxp1;
        id[4] = id[0]+gnxp1*gnyp1;
        id[5] = id[4]+1;
        id[6] = id[4]+gnxp1+1;
        id[7] = id[4]+gnxp1;
        cells[cell_id].vertices[0] = 8;
        for (int ii=0; ii<8; ii++) {
          cells[cell_id].vertices[ii+1] = id[ii];
          vertices[id[ii]].cells[0]++;
          vertices[id[ii]].cells[vertices[id[ii]].cells[0]] = cell_id;
        }
      }
    }
  }
}

void StructuredGrid::computeCellMapping(int *num_cells_local, int *num_cells_ghosted,
                                      int **cell_mapping_local_to_ghosted,
                                      int **cell_mapping_ghosted_to_local,
                                      int **cell_mapping_ghosted_to_natural) {

  *num_cells_local = lnx*lny*lnz;
  *num_cells_ghosted = gnx*gny*gnz;

  *cell_mapping_local_to_ghosted = new int[*num_cells_local];
  *cell_mapping_ghosted_to_local = new int[*num_cells_ghosted];
  *cell_mapping_ghosted_to_natural = new int[*num_cells_ghosted];

  for (int i=0; i<*num_cells_local; i++) (*cell_mapping_local_to_ghosted)[i] = -1;
  for (int i=0; i<*num_cells_ghosted; i++) (*cell_mapping_ghosted_to_local)[i] = -1;
  for (int i=0; i<*num_cells_ghosted; i++) (*cell_mapping_ghosted_to_natural)[i] = -1;
  
  int istart = lxs-gxs;
  int jstart = lys-gys;
  int kstart = lzs-gzs;
  int iend = istart+lnx;
  int jend = jstart+lny;
  int kend = kstart+lnz;
  int local_id = 0;
  for (int k=kstart; k<kend; k++) {
    for (int j=jstart; j<jend; j++) {
      for (int i=istart; i<iend; i++) {
        int ghosted_id = i+j*gnx+k*gnxXny;
        int natural_id = (i+gxs)+(j+gys)*nx+(k+gzs)*nx*ny;
        (*cell_mapping_local_to_ghosted)[local_id] = ghosted_id;
        (*cell_mapping_ghosted_to_local)[ghosted_id] = local_id;
        (*cell_mapping_ghosted_to_natural)[ghosted_id] = natural_id;
        local_id = local_id + 1;
      }
    }
  }

}

void StructuredGrid::computeVertexMapping(int *num_vertices_local, 
                                          int *num_vertices_ghosted,
                                          int **vertex_mapping_local_to_ghosted,
                                          int **vertex_mapping_ghosted_to_local,
                                          int **vertex_mapping_ghosted_to_natural) {

  *num_vertices_local = (lnx+1)*(lny+1)*(lnz+1);
  *num_vertices_ghosted = (gnx+1)*(gny+1)*(gnz+1);

  *vertex_mapping_local_to_ghosted = new int[*num_vertices_local];
  *vertex_mapping_ghosted_to_local = new int[*num_vertices_ghosted];
  *vertex_mapping_ghosted_to_natural = new int[*num_vertices_ghosted];

  for (int i=0; i<*num_vertices_local; i++) (*vertex_mapping_local_to_ghosted)[i] = -1;
  for (int i=0; i<*num_vertices_ghosted; i++) (*vertex_mapping_ghosted_to_local)[i] = -1;
  for (int i=0; i<*num_vertices_ghosted; i++) (*vertex_mapping_ghosted_to_natural)[i] = -1;
  
  int nxp1 = nx+1;
  int nyp1 = ny+1;
  int nzp1 = nz+1;

  int gnxp1 = gnx+1;
  int gnyp1 = gny+1;
  int gnzp1 = gnz+1;
  int gnxp1Xgnyp1 = gnxp1*gnyp1;

  int istart = lxs-gxs;
  int jstart = lys-gys;
  int kstart = lzs-gzs;
  int iend = istart+lnx+1; // +1 for vertex 
  int jend = jstart+lny+1;
  int kend = kstart+lnz+1;
  int local_id = 0;
  for (int k=kstart; k<kend; k++) {
    for (int j=jstart; j<jend; j++) {
      for (int i=istart; i<iend; i++) {
        int ghosted_id = i+j*gnxp1+k*gnxp1Xgnyp1;
        int natural_id = (i+gxs)+(j+gys)*nxp1+(k+gzs)*nxp1*nyp1;
        (*vertex_mapping_local_to_ghosted)[local_id] = ghosted_id;
        (*vertex_mapping_ghosted_to_local)[ghosted_id] = local_id;
        (*vertex_mapping_ghosted_to_natural)[ghosted_id] = natural_id;
        local_id = local_id + 1;
      }
    }
  }

}

int *StructuredGrid::getLocalCellVertexNaturalIDs(GridCell *cells, GridVertex *vertices) {

  Vec global;
  Vec natural;
  PetscScalar *v_ptr;

  int num_cells_local = lnx*lny*lnz;
  int num_cells_ghosted = gnx*gny*gnz;

  DACreateGlobalVector(da,&global);
  DACreateNaturalVector(da,&natural);

  int *vertex_ids = new int[num_cells_local*8];
  for (int ivert=0; ivert<8; ivert++) {

    VecGetArray(global,&v_ptr);
    for (int icell=0; icell<num_cells_ghosted; icell++) {
      int cell_id_local = cells[icell].getIdLocal();
      if (cell_id_local > -1) {
          v_ptr[cell_id_local] = (double)vertices[cells[icell].vertices[ivert+1]].getIdNatural();
      }
    }
    VecRestoreArray(global,&v_ptr);

    globalToNatural(global,natural);

    VecGetArray(natural,&v_ptr);
    for (int i=0; i<num_cells_local; i++)
      vertex_ids[ivert+i*8] = (int)(v_ptr[i]+0.0001);
    VecRestoreArray(natural,&v_ptr);
  }

  VecDestroy(global);
  VecDestroy(natural);

  return vertex_ids;

}

void StructuredGrid::getVectorLocal(Vec *v) {
  DACreateLocalVector(da,v);
}

void StructuredGrid::getVectorGlobal(Vec *v) {
  DACreateGlobalVector(da,v);
}

void StructuredGrid::getVectorNatural(Vec *v) {
  DACreateNaturalVector(da,v);
}

void StructuredGrid::globalToNatural(Vec global, Vec natural) {
  DAGlobalToNaturalBegin(da,global,INSERT_VALUES,natural);
  DAGlobalToNaturalEnd(da,global,INSERT_VALUES,natural);
}

void StructuredGrid::convertLocalCellDataGtoN(double *data) {

  Vec global;
  Vec natural;
  PetscScalar *v_ptr = NULL;

  int num_cells_local = lnx*lny*lnz;

  DACreateGlobalVector(da,&global);
  DACreateNaturalVector(da,&natural);

  VecGetArray(global,&v_ptr);
  for (int i=0; i<num_cells_local; i++)
    v_ptr[i] = data[i];
  VecRestoreArray(global,&v_ptr);

  globalToNatural(global,natural);

  VecGetArray(global,&v_ptr);
  for (int i=0; i<num_cells_local; i++)
    data[i] = v_ptr[i];
  VecRestoreArray(global,&v_ptr);

  VecDestroy(global);
  VecDestroy(natural);

}

void StructuredGrid::printDACorners() {
  
  PetscErrorCode ierr;
  
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
  printf("local         %d lxs:%d lxe:%d lnx:%d\n",myrank,lxs,lxe,lnx);
  printf("local         %d lys:%d lye:%d lny:%d\n",myrank,lys,lye,lny);
  printf("local         %d lzs:%d lze:%d lnz:%d\n",myrank,lzs,lze,lnz);

  printf("local ghosted %d gxs:%d gxe:%d gnx:%d\n",myrank,gxs,gxe,gnx);
  printf("local ghosted %d gys:%d gye:%d gny:%d\n",myrank,gys,gye,gny);
  printf("local ghosted %d gzs:%d gze:%d gnz:%d\n\n",myrank,gzs,gze,gnz);
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
  
}

int StructuredGrid::getNx() { return nx; }
int StructuredGrid::getNy() { return ny; }
int StructuredGrid::getNz() { return nz; }
int StructuredGrid::getN() { return nx*ny*nz; }


void StructuredGrid::getCorners(int *xs, int *ys, int *zs, 
                                int *nx, int *ny, int *nz) {
  *xs = lxs; *ys = lys; *zs = lzs; 
  *nx = lnx; *ny = lny; *nz = lnz; 
}

void StructuredGrid::getGhostCorners(int *xs, int *ys, int *zs, 
                                     int *nx, int *ny, int *nz) {
  *xs = gxs; *ys = gys; *zs = gzs; 
  *nx = gnx; *ny = gny; *nz = gnz; 
}

double StructuredGrid::getDx(int i) { return dx[i]; }
double StructuredGrid::getDy(int j) { return dy[j]; }
double StructuredGrid::getDz(int k) { return dz[k]; }
double *StructuredGrid::getOriginPtr() { return &(global_origin[0]); }
double StructuredGrid::getRotationDegrees() { return rotationZdegrees; }

void StructuredGrid::printAO() {

#include "include/petscao.h"

  PetscErrorCode ierr;
  AO ao;
  ierr = DAGetAO(da,&ao);
  AOView(ao,PETSC_VIEWER_STDOUT_WORLD);
  AODestroy(ao);

}



StructuredGrid::~StructuredGrid() {
  DADestroy(da);
}