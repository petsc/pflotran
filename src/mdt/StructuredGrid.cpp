#include "StructuredGrid.h"

StructuredGrid::StructuredGrid() {
}

StructuredGrid::StructuredGrid(PetscInt nx_, PetscInt ny_, PetscInt nz_) {
  nx = nx_; ny = ny_; nz = nz_;
  global_origin[0] = 0.; global_origin[1] = 0.; global_origin[2] = 0.;
  local_origin[0] = -999.; local_origin[1] = -999.; local_origin[2] = -999.;
  rotationZdegrees = 0.; rotationZradians = 0.;
}

void StructuredGrid::createDA() {

  PetscErrorCode ierr;

  PetscMPIInt commsize;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);
  PetscInt pnx = PETSC_DECIDE;
  PetscInt pnz = PETSC_DECIDE;
  if (commsize > 1) pnx = 2; 
  if (commsize > 3) pnz = 2; 

  ierr = DMDACreate3d(PETSC_COMM_WORLD,
                    DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
                    DMDA_STENCIL_STAR,
                    nx,ny,nz,pnx,PETSC_DECIDE,pnz,
                    1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&da);

  ierr = DMDAGetCorners(da,&lxs,&lys,&lzs,&lnx,&lny,&lnz);
  lxe = lxs+lnx;
  lye = lys+lny;
  lze = lzs+lnz;
  lnxXny = lnx*lny;
  
  ierr = DMDAGetGhostCorners(da,&gxs,&gys,&gzs,&gnx,&gny,&gnz);
  gxe = gxs+gnx;
  gye = gys+gny;
  gze = gzs+gnz;
  gnxXny = gnx*gny;
  
  printf(" local(%d): %d %d %d %d %d %d\n",myrank,lxs,lxe,lys,lye,lzs,lze);
  printf("global(%d): %d %d %d %d %d %d\n",myrank,gxs,gxe,gys,gye,gzs,gze);

}

void StructuredGrid::setGridSpacing(PetscReal *dx_, PetscReal *dy_, PetscReal *dz_) {

  dx = new PetscReal[nx];
  for (PetscInt i=0; i<nx; i++)
    dx[i] = dx_[i];
  dy = new PetscReal[ny];
  for (PetscInt j=0; j<ny; j++)
    dy[j] = dy_[j];
  dz = new PetscReal[nz];
  for (PetscInt k=0; k<nz; k++)
    dz[k] = dz_[k];

  setLocalGridSpacing();

}

void StructuredGrid::setGridSpacing(PetscReal dx_, PetscReal dy_, PetscReal dz_) {

  dx = new PetscReal[nx];
  for (PetscInt i=0; i<nx; i++)
    dx[i] = dx_;
  dy = new PetscReal[ny];
  for (PetscInt j=0; j<ny; j++)
    dy[j] = dy_;
  dz = new PetscReal[nz];
  for (PetscInt k=0; k<nz; k++)
    dz[k] = dz_;

  setLocalGridSpacing();

}

void StructuredGrid::setLocalGridSpacing() {

  if (!dx || !dy || !dz) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Global grid spacing must be set before local.\n");
    PetscFinalize();
    exit(0);
  }

  gdx = new PetscReal[gnx];
  for (PetscInt i=gxs; i<gxe; i++)
    gdx[i-gxs] = dx[i];
  gdy = new PetscReal[gny];
  for (PetscInt j=gys; j<gye; j++)
    gdy[j-gys] = dy[j];
  gdz = new PetscReal[gnz];
  for (PetscInt k=gzs; k<gze; k++)
    gdz[k-gzs] = dz[k];

  ldx = new PetscReal[lnx];
  for (PetscInt i=lxs; i<lxe; i++)
    ldx[i-lxs] = dx[i];
  ldy = new PetscReal[lny];
  for (PetscInt j=lys; j<lye; j++)
    ldy[j-lys] = dy[j];
  ldz = new PetscReal[lnz];
  for (PetscInt k=lzs; k<lze; k++)
    ldz[k-lzs] = dz[k];

}

void StructuredGrid::setOrigin(PetscReal origx, PetscReal origy, PetscReal origz) {
  
  if (!dx || !dy || !dz) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Global grid spacing must be set before origin.\n");
    PetscFinalize();
    exit(0);
  }

  global_origin[0] = origx;
  global_origin[1] = origy;
  global_origin[2] = origz;

  PetscReal sumx = 0.;

  for (PetscInt i=0; i<lxs; i++)
    sumx += dx[i];

  PetscReal sumy = 0.;
  for (PetscInt j=0; j<lys; j++)
    sumy += dy[j];

  PetscReal sumz = 0.;
  for (PetscInt k=0; k<lzs; k++)
    sumz += dz[k];

  local_origin[0] = global_origin[0] + sumx*cos(rotationZradians) - 
                                       sumy*sin(rotationZradians);
  local_origin[1] = global_origin[1] + sumx*sin(rotationZradians) + 
                                       sumy*cos(rotationZradians);
  local_origin[2] = global_origin[2] + sumz;

  /*
  delete [] dx;
  dx = NULL;
  delete [] dy;
  dy = NULL;
  delete [] dz;
  dz = NULL;
  */
}

void StructuredGrid::setRotation(PetscReal r) {
  
  if (PetscAbsReal(local_origin[0] - -999.) > 1.e-40 &&
      PetscAbsReal(local_origin[1] - -999.) > 1.e-40 &&
      PetscAbsReal(local_origin[2] - -999.) > 1.e-40) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Rotation must be set before setting origin.\n");
    PetscFinalize();
    exit(0);
  }

  rotationZdegrees = r;
  rotationZradians = r/180.*PI; // convert from degrees to radians
}

#if 0
void StructuredGrid::computeConnectivity(PetscInt *num_connections, 
                                         GridConnection **connections) {


// all connectivity is in the local ghosted context
  *num_connections = (gnx-1)*gny*gnz+gnx*(gny-1)*gnz+gnx*gny*(gnz-1);

  *connections = new GridConnection[*num_connections];
  GridConnection *conn = *connections;
  // conn[count].xx is the same as (*connections)[count].xx

  PetscReal vec[3];
  
// x-direction
  vec[1] = 0.; vec[2] = 0.;
  PetscInt count = 0;
  for (PetscInt k=0; k<gnz; k++) {
    for (PetscInt j=0; j<gny; j++) {
      for (PetscInt i=0; i<gnx-1; i++) {
        PetscInt id = i+j*gnx+k*gnxXny;
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
  for (PetscInt k=0; k<gnz; k++) {
    for (PetscInt i=0; i<gnx; i++) {
      for (PetscInt j=0; j<gny-1; j++) {
        PetscInt id = i+j*gnx+k*gnxXny;
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
  for (PetscInt j=0; j<gny; j++) {
    for (PetscInt i=0; i<gnx; i++) {
      for (PetscInt k=0; k<gnz-1; k++) {
        PetscInt id = i+j*gnx+k*gnxXny;
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

  if (PetscAbsReal(local_origin[0] - -999.) < 1.e-40 &&
      PetscAbsReal(local_origin[1] - -999.) < 1.e-40 &&
      PetscAbsReal(local_origin[2] - -999.) < 1.e-40) {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Origin must be set before computing coordinates.\n");
    PetscFinalize();
    exit(0);
  }

  Vec v;
  PetscReal *v_ptr = NULL;
  VecCreate(PETSC_COMM_WORLD,&v);
  VecSetSizes(v,lnx*lny*lnz*3,PETSC_DECIDE);
  VecSetType(v,VECMPI);
  VecGetArray(v,&v_ptr);

  PetscReal sumz = 0.5*ldz[0];     
  for (PetscInt k=1; k<lnz; k++) {   
    PetscReal sumy = 0.5*ldy[0];
    for (PetscInt j=1; j<lny; j++) {
      PetscReal sumx = 0.5*ldx[0];
      for (PetscInt i=1; i<lnx; i++) {
        PetscInt local_id = i+j*lnx+k*lnxXny;
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

  DMDASetCoordinates(da,v);

}

void StructuredGrid::computeCoordinate(PetscReal *x, PetscReal*y) {

  PetscReal distx = *x - global_origin[0];
  PetscReal disty = *y - global_origin[1];
  *x = global_origin[0] + distx*cos(rotationZradians) - 
                          disty*sin(rotationZradians);
  *y = global_origin[1] + distx*sin(rotationZradians) + 
                          disty*cos(rotationZradians);
}

#if 0
void StructuredGrid::mapBoundaryConnection(PetscInt istart, PetscInt iend, PetscInt jstart,
                                          PetscInt jend, PetscInt kstart, PetscInt kend,
                                          char *face) {

  // ensure monotonic increasing delineation of region
  if (istart > iend) {
    PetscInt temp = istart;
    istart =  iend;
    iend = temp;
  }
  if (jstart > jend) {
    PetscInt temp = jstart;
    jstart =  jend;
    jend = temp;
  }
  if (kstart > kend) {
    PetscInt temp = kstart;
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

  PetscInt num_connections = PetscAbsInt(iend-istart)*PetscAbsInt(jend-jstart)*
                             PetscAbsInt(kend-kstart);

//  printf("3 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);

  PetscInt count = 0;
  PetscReal dist[3];
  PetscReal center[3];
  PetscReal normal[3];
//  printf("4 %d %d %d %d %d %d\n",lxs,lxe,lys,lye,lzs,lze);
  for (PetscInt k=kstart-lzs; k<kend-lzs; k++) {
    for (PetscInt j=jstart-lys; j<jend-lys; j++) {
      for (PetscInt i=istart-lxs; i<iend-lxs; i++) {
        PetscInt ighosted = i+lxs-gxs;
        PetscInt jghosted = j+lys-gys;
        PetscInt kghosted = k+lzs-gzs;
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
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt icell = i;
      PetscReal dx_up;
      PetscReal dx_dn;
      if (i > 0) dx_up = cells[icell].getX()-cells[icell-1].getX();
      else dx_up = cells[1].getX()-cells[0].getX();
      if (i < gnx-1) dx_dn = cells[icell+1].getX()-cells[icell].getX();
      else dx_dn = cells[gnx-1].getX()-cells[gnx-2].getX();
      dx[i] = 0.5*(dx_up+dx_dn);
    }
  }
  else dx[0] = 1.;
  if (gny > 1) {
    for (PetscInt j=0; j<gny; j++) {
      PetscInt icell = j*gnx;
      PetscReal dy_up;
      PetscReal dy_dn;
      if (j > 0) dy_up = cells[icell].getY()-cells[icell-gnx].getY();
      else dy_up = cells[gnx].getY()-cells[0].getY();
      if (j < gny-1) dy_dn = cells[icell+gnx].getY()-cells[icell].getY();
      else dy_dn = cells[(gny-1)*gnx].getY()-cells[(gny-2)*gnx].getY();
      dy[j] = 0.5*(dy_up+dy_dn);
    }
  }
  else dy[0] = 1.;
  if (gnz > 1) {
    for (PetscInt k=0; k<gnz; k++) {
      PetscInt icell = k*gnxXny;
      PetscReal dz_up;
      PetscReal dz_dn;
      if (k > 0) dz_up = cells[icell].getZ()-cells[icell-gnxXny].getZ();
      else dz_up = cells[gnxXny].getZ()-cells[0].getZ();
      if (k < gnz-1) dz_dn = cells[icell+gnxXny].getZ()-cells[icell].getZ();
      else dz_dn = cells[(gnz-1)*gnxXny].getZ()-cells[(gnz-2)*gnxXny].getZ();
      dz[k] = 0.5*(dz_up+dz_dn);
    }
  }
  else dz[0] = 1.;

}

void StructuredGrid::mapSource(PetscInt istart, PetscInt iend, PetscInt jstart,
                               PetscInt jend, PetscInt kstart, PetscInt kend) {

  // ensure monotonic increasing delineation of region
  if (istart > iend) {
    PetscInt temp = istart;
    istart =  iend;
    iend = temp;
  }
  if (jstart > jend) {
    PetscInt temp = jstart;
    jstart =  jend;
    jend = temp;
  }
  if (kstart > kend) {
    PetscInt temp = kstart;
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

  PetscInt num_connections = PetscAbsInt(iend-istart)*PetscAbsInt(jend-jstart)*
                             PetscAbsInt(kend-kstart);

//  printf("3 %d %d %d %d %d %d\n",istart,iend,jstart,jend,kstart,kend);

  PetscInt count = 0;
//  printf("4 %d %d %d %d %d %d\n",lxs,lxe,lys,lye,lzs,lze);
  for (PetscInt k=kstart-lzs; k<kend-lzs; k++) {
    for (PetscInt j=jstart-lys; j<jend-lys; j++) {
      for (PetscInt i=istart-lxs; i<iend-lxs; i++) {
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

void StructuredGrid::setUpCells(PetscInt num_cells, GridCell *cells) {
  PetscReal z;
  if (lzs-gzs > 0) z = -0.5*gdz[0];
  else z = 0.5*gdz[0];
  for (PetscInt k=0; k<gnz; k++) {
    PetscReal y;
    if (lys-gys > 0) y = -0.5*gdy[0];
    else y = 0.5*gdy[0];
    for (PetscInt j=0; j<gny; j++) {
      PetscReal x;
      if (lxs-gxs > 0) x = -0.5*gdx[0];
      else x = 0.5*gdx[0];
      for (PetscInt i=0; i<gnx; i++) {
        PetscInt id = i+j*gnx+k*gnxXny;
        cells[id].setCentroidLocal(x,y,z);
        cells[id].setVolume(gdx[i]*gdy[j]*gdz[k]);
        PetscReal x_rot = local_origin[0] + x*cos(rotationZradians) -
                                            y*sin(rotationZradians);
        PetscReal y_rot = local_origin[1] + x*sin(rotationZradians) +
                                            y*cos(rotationZradians);
        PetscReal z_adj = local_origin[2] + z;
        cells[id].setCentroid(x_rot,y_rot,z_adj);
        if (i < gnx-1) x += 0.5*(gdx[i]+gdx[i+1]);
      }
      if (j < gny-1) y += 0.5*(gdy[j]+gdy[j+1]);
    }
    if (k < gnz-1) z += 0.5*(gdz[k]+gdz[k+1]);
  }
}

void StructuredGrid::setUpVertices(PetscInt num_vertices, GridVertex *vertices) {
  PetscInt gnxp1 = gnx+1;
  PetscInt gnyp1 = gny+1;
  PetscInt gnzp1 = gnz+1;

  PetscReal z = 0.;
  for (PetscInt k=0; k<gnzp1; k++) {
    PetscReal y = 0.;
    for (PetscInt j=0; j<gnyp1; j++) {
      PetscReal x = 0.;
      for (PetscInt i=0; i<gnxp1; i++) {
        PetscInt id = i+j*gnxp1+k*gnxp1*gnyp1;
        PetscReal x_rot = local_origin[0] + (x-(lxs-gxs>0?gdx[0]:0.))*cos(rotationZradians) - 
                                         (y-(lys-gys>0?gdy[0]:0.))*sin(rotationZradians);
        PetscReal y_rot = local_origin[1] + (x-(lxs-gxs>0?gdx[0]:0.))*sin(rotationZradians) + 
                                         (y-(lys-gys>0?gdy[0]:0.))*cos(rotationZradians);
        PetscReal z_adj = local_origin[2] + (z-(lzs-gzs>0?gdz[0]:0.));
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
  PetscInt gnxp1 = gnx+1;
  PetscInt gnyp1 = gny+1;
  PetscInt gnzp1 = gnz+1;
  for (PetscInt k=0; k<gnz; k++) {
    for (PetscInt j=0; j<gny; j++) {
      for (PetscInt i=0; i<gnx; i++) {
        PetscInt cell_id = i+j*gnx+k*gnxXny;
        PetscInt id[8];
        id[0] = i+j*gnxp1+k*gnxp1*gnyp1;
        id[1] = id[0]+1;
        id[2] = id[0]+gnxp1+1;
        id[3] = id[0]+gnxp1;
        id[4] = id[0]+gnxp1*gnyp1;
        id[5] = id[4]+1;
        id[6] = id[4]+gnxp1+1;
        id[7] = id[4]+gnxp1;
        cells[cell_id].vertices[0] = 8;
        for (PetscInt ii=0; ii<8; ii++) {
          cells[cell_id].vertices[ii+1] = id[ii];
          vertices[id[ii]].cells[0]++;
          vertices[id[ii]].cells[vertices[id[ii]].cells[0]] = cell_id;
        }
      }
    }
  }
}

void StructuredGrid::computeCellMapping(PetscInt *num_cells_local, PetscInt *num_cells_ghosted,
                                      PetscInt **cell_mapping_local_to_ghosted,
                                      PetscInt **cell_mapping_ghosted_to_local,
                                      PetscInt **cell_mapping_ghosted_to_natural) {

  *num_cells_local = lnx*lny*lnz;
  *num_cells_ghosted = gnx*gny*gnz;

  *cell_mapping_local_to_ghosted = new PetscInt[*num_cells_local];
  *cell_mapping_ghosted_to_local = new PetscInt[*num_cells_ghosted];
  *cell_mapping_ghosted_to_natural = new PetscInt[*num_cells_ghosted];

  for (PetscInt i=0; i<*num_cells_local; i++) (*cell_mapping_local_to_ghosted)[i] = -1;
  for (PetscInt i=0; i<*num_cells_ghosted; i++) (*cell_mapping_ghosted_to_local)[i] = -1;
  for (PetscInt i=0; i<*num_cells_ghosted; i++) (*cell_mapping_ghosted_to_natural)[i] = -1;
  
  PetscInt istart = lxs-gxs;
  PetscInt jstart = lys-gys;
  PetscInt kstart = lzs-gzs;
  PetscInt iend = istart+lnx;
  PetscInt jend = jstart+lny;
  PetscInt kend = kstart+lnz;
  PetscInt local_id = 0;
  for (PetscInt k=kstart; k<kend; k++) {
    for (PetscInt j=jstart; j<jend; j++) {
      for (PetscInt i=istart; i<iend; i++) {
        PetscInt ghosted_id = i+j*gnx+k*gnxXny;
        (*cell_mapping_local_to_ghosted)[local_id] = ghosted_id;
        (*cell_mapping_ghosted_to_local)[ghosted_id] = local_id;
        local_id = local_id + 1;
      }
    }
  }
  PetscInt ghosted_id = 0;
  for (PetscInt k=0; k<gnz; k++) {
    for (PetscInt j=0; j<gny; j++) {
      for (PetscInt i=0; i<gnx; i++) {
        PetscInt natural_id = (i+gxs)+(j+gys)*nx+(k+gzs)*nx*ny;
        (*cell_mapping_ghosted_to_natural)[ghosted_id] = natural_id;
        ghosted_id = ghosted_id + 1;
      }
    }
  }

}

void StructuredGrid::computeVertexMapping(PetscInt *num_vertices_local, 
                                          PetscInt *num_vertices_ghosted,
                                          PetscInt **vertex_mapping_local_to_ghosted,
                                          PetscInt **vertex_mapping_ghosted_to_local,
                                          PetscInt **vertex_mapping_ghosted_to_natural) {

  *num_vertices_local = (lnx+1)*(lny+1)*(lnz+1);
  *num_vertices_ghosted = (gnx+1)*(gny+1)*(gnz+1);

  *vertex_mapping_local_to_ghosted = new PetscInt[*num_vertices_local];
  *vertex_mapping_ghosted_to_local = new PetscInt[*num_vertices_ghosted];
  *vertex_mapping_ghosted_to_natural = new PetscInt[*num_vertices_ghosted];

  for (PetscInt i=0; i<*num_vertices_local; i++) (*vertex_mapping_local_to_ghosted)[i] = -1;
  for (PetscInt i=0; i<*num_vertices_ghosted; i++) (*vertex_mapping_ghosted_to_local)[i] = -1;
  for (PetscInt i=0; i<*num_vertices_ghosted; i++) (*vertex_mapping_ghosted_to_natural)[i] = -1;
  
  PetscInt nxp1 = nx+1;
  PetscInt nyp1 = ny+1;
  PetscInt nzp1 = nz+1;

  PetscInt gnxp1 = gnx+1;
  PetscInt gnyp1 = gny+1;
  PetscInt gnzp1 = gnz+1;
  PetscInt gnxp1Xgnyp1 = gnxp1*gnyp1;

  PetscInt istart = lxs-gxs;
  PetscInt jstart = lys-gys;
  PetscInt kstart = lzs-gzs;
  PetscInt iend = istart+lnx+1; // +1 for vertex 
  PetscInt jend = jstart+lny+1;
  PetscInt kend = kstart+lnz+1;
  PetscInt local_id = 0;
  for (PetscInt k=kstart; k<kend; k++) {
    for (PetscInt j=jstart; j<jend; j++) {
      for (PetscInt i=istart; i<iend; i++) {
        PetscInt ghosted_id = i+j*gnxp1+k*gnxp1Xgnyp1;
        PetscInt natural_id = (i+gxs)+(j+gys)*nxp1+(k+gzs)*nxp1*nyp1;
        (*vertex_mapping_local_to_ghosted)[local_id] = ghosted_id;
        (*vertex_mapping_ghosted_to_local)[ghosted_id] = local_id;
        local_id = local_id + 1;
      }
    }
  }
  PetscInt ghosted_id = 0;
  for (PetscInt k=0; k<gnzp1; k++) {
    for (PetscInt j=0; j<gnyp1; j++) {
      for (PetscInt i=0; i<gnxp1; i++) {
        PetscInt natural_id = (i+gxs)+(j+gys)*nxp1+(k+gzs)*nxp1*nyp1;
        (*vertex_mapping_ghosted_to_natural)[ghosted_id] = natural_id;
        ghosted_id = ghosted_id + 1;
      }
    }
  }

}
#if 0
PetscInt *StructuredGrid::getLocalCellVertexNaturalIDs(GridCell *cells, GridVertex *vertices) {

  Vec vec;
  PetscReal *v_ptr;

  PetscInt num_vertices_local = (lnx+1)*(lny+1)*(lnz+1);

  PetscInt *vertex_ids = new PetscInt[num_verts_local];
  for (PetscInt ivert=0; ivert<8; ivert++) {

    VecGetArray(global,&v_ptr);
    for (PetscInt icell=0; icell<num_cells_ghosted; icell++) {
      PetscInt cell_id_local = cells[icell].getIdLocal();
      if (cell_id_local > -1) {
          vertex_ids[cell_id_local] = (PetscReal)vertices[cells[icell].vertices[ivert+1]].getIdNatural();
      }
    }
    VecRestoreArray(global,&v_ptr);

    globalToNatural(global,natural);

    VecGetArray(natural,&v_ptr);
    for (PetscInt i=0; i<num_cells_local; i++)
      vertex_ids[ivert+i*8] = (int)(v_ptr[i]+0.0001);
    VecRestoreArray(natural,&v_ptr);
  }

  VecDestroy(&global);
  VecDestroy(&natural);

  return vertex_ids;

}
#endif
void StructuredGrid::getVectorLocal(Vec *v) {
  DMCreateLocalVector(da,v);
}

void StructuredGrid::getVectorGlobal(Vec *v) {
  DMCreateGlobalVector(da,v);
}

void StructuredGrid::getVectorNatural(Vec *v) {
  DMDACreateNaturalVector(da,v);
}

void StructuredGrid::globalToNatural(Vec global, Vec natural) {
  DMDAGlobalToNaturalBegin(da,global,INSERT_VALUES,natural);
  DMDAGlobalToNaturalEnd(da,global,INSERT_VALUES,natural);
}

void StructuredGrid::convertLocalCellDataGtoN(PetscReal *data) {

  Vec global;
  Vec natural;
  PetscReal *v_ptr = NULL;

  PetscInt num_cells_local = lnx*lny*lnz;

  DMCreateGlobalVector(da,&global);
  DMDACreateNaturalVector(da,&natural);

  VecGetArray(global,&v_ptr);
  for (PetscInt i=0; i<num_cells_local; i++)
    v_ptr[i] = data[i];
  VecRestoreArray(global,&v_ptr);

  globalToNatural(global,natural);

  VecGetArray(natural,&v_ptr);
  for (PetscInt i=0; i<num_cells_local; i++)
    data[i] = v_ptr[i];
  VecRestoreArray(natural,&v_ptr);

/*
  PetscErrorCode ierr;
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
  VecGetArray(global,&v_ptr);
  printf("proc: %d - global - ",myrank);
  for (PetscInt i=0; i<num_cells_local; i++)
    printf(" %.1f",v_ptr[i]);
  printf("\n");
  VecRestoreArray(global,&v_ptr);
  VecGetArray(natural,&v_ptr);
  printf("proc: %d - natural - ",myrank);
  for (PetscInt i=0; i<num_cells_local; i++)
    printf(" %.1f",v_ptr[i]);
  printf("\n");
  VecRestoreArray(natural,&v_ptr);
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
*/

  VecDestroy(&global);
  VecDestroy(&natural);

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

PetscInt StructuredGrid::getNx() { return nx; }
PetscInt StructuredGrid::getNy() { return ny; }
PetscInt StructuredGrid::getNz() { return nz; }
PetscInt StructuredGrid::getN() { return nx*ny*nz; }


void StructuredGrid::getCorners(PetscInt *xs, PetscInt *ys, PetscInt *zs, 
                                PetscInt *nx, PetscInt *ny, PetscInt *nz) {
  *xs = lxs; *ys = lys; *zs = lzs; 
  *nx = lnx; *ny = lny; *nz = lnz; 
}

void StructuredGrid::getGhostCorners(PetscInt *xs, PetscInt *ys, PetscInt *zs, 
                                     PetscInt *nx, PetscInt *ny, PetscInt *nz) {
  *xs = gxs; *ys = gys; *zs = gzs; 
  *nx = gnx; *ny = gny; *nz = gnz; 
}

PetscReal StructuredGrid::getDx(PetscInt i) { return dx[i]; }
PetscReal StructuredGrid::getDy(PetscInt j) { return dy[j]; }
PetscReal StructuredGrid::getDz(PetscInt k) { return dz[k]; }
PetscReal *StructuredGrid::getOriginPtr() { return &(global_origin[0]); }
PetscReal StructuredGrid::getRotationDegrees() { return rotationZdegrees; }

PetscInt StructuredGrid::getNeighboringProcessor(PetscInt direction) {

  PetscInt pnx, pny, pnz;
  DMDAGetInfo(da,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
            &pnx,&pny,&pnz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
            PETSC_NULL,PETSC_NULL);
  PetscInt iproc = myrank%pnx;
  PetscInt jproc = myrank%(pnx*pny)/pnx;
  PetscInt kproc = myrank/(pnx*pny);
  PetscInt neighbor;
  switch(direction) {
    case(NORTH):
      neighbor = jproc < pny-1 ? iproc+(jproc+1)*pnx+kproc*(pnx*pny) : 
                                 iproc+kproc*(pnx*pny); // cycle
      break;
    case(SOUTH):
      neighbor = jproc > 0 ? iproc+(jproc-1)*pnx+kproc*(pnx*pny) : 
                             iproc+(pny-1)*pnx+kproc*(pnx*pny);
      break;
    case(EAST):
      neighbor = iproc < pnx-1 ? iproc+1+jproc*pnx+kproc*(pnx*pny) :
                                 jproc*pnx+kproc*(pnx*pny);
      break;
    case(WEST):
      neighbor = iproc > 0 ? iproc-1+jproc*pnx+kproc*(pnx*pny) : 
                             pnx-1+jproc*pnx+kproc*(pnx*pny);
      break;
    case(TOP):
      neighbor = kproc < pnz-1 ? iproc+jproc*pnx+(kproc+1)*(pnx*pny) :
                                 iproc+jproc*pnx;
      break;
    case(BOTTOM):
      neighbor = kproc > 0 ? iproc+jproc*pnx+(kproc-1)*(pnx*pny) :
                             iproc+jproc*pnx+(pnz-1)*(pnx*pny);
      break;
  }
  return neighbor;
}

void StructuredGrid::sendFlag(PetscInt *flag, PetscInt direction) {
  MPI_Send(&flag,3,MPIU_INT,getNeighboringProcessor(direction),0,
           PETSC_COMM_WORLD);
}

void StructuredGrid::receiveFlag(PetscInt *flag, PetscInt direction) {
  MPI_Status status;
  MPI_Recv(flag,3,MPIU_INT,getNeighboringProcessor(direction),0,
           PETSC_COMM_WORLD,&status);
}

void StructuredGrid::printAO() {

#include "petscao.h"

  PetscErrorCode ierr;
  AO ao;
  ierr = DMDAGetAO(da,&ao);
  AOView(ao,PETSC_VIEWER_STDOUT_WORLD);
  AODestroy(&ao);

}



StructuredGrid::~StructuredGrid() {
  DMDestroy(&da);
}
