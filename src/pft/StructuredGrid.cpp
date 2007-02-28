#include "StructuredGrid.h"

StructuredGrid::StructuredGrid() {
}

StructuredGrid::StructuredGrid(int nx_, int ny_, int nz_, int fdof_, int tdof_) {
  nx = nx_; ny = ny_; nz = nz_; fdof = fdof_; tdof = tdof_;
  setUpDA(1.,1.,1.);
}

StructuredGrid::StructuredGrid(int nx_, int ny_, int nz_, double dx_,
                               double dy_, double dz_, int fdof_, int tdof_) {
  nx = nx_; ny = ny_; nz = nz_; fdof = fdof_; tdof = tdof_;
  setUpDA(dx_,dy_,dz_);
}

void StructuredGrid::setUpDA(double dx_, double dy_, double dz_) {

  PetscErrorCode ierr;

  int commsize;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);
  int pnx = PETSC_DECIDE;
  int pnz = PETSC_DECIDE;
  if (commsize > 1) pnx = 2; 
  if (commsize > 3) pnz = 2; 

  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
                    nx,ny,nz,pnx,PETSC_DECIDE,pnz,
                    1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&da_1dof);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
                    nx,ny,nz,pnx,PETSC_DECIDE,pnz,
                    fdof,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&da_fdof);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
                    nx,ny,nz,pnx,PETSC_DECIDE,pnz,
                    tdof,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&da_tdof);

  ierr = DAGetCorners(da_1dof,&lxs,&lys,&lzs,&lnx,&lny,&lnz);
  lxe = lxs+lnx;
  lye = lys+lny;
  lze = lzs+lnz;
  lnxXny = lnx*lny;
  
  ierr = DAGetGhostCorners(da_1dof,&gxs,&gys,&gzs,&gnx,&gny,&gnz);
  gxe = gxs+gnx;
  gye = gys+gny;
  gze = gzs+gnz;
  gnxXny = gnx*gny;
  
  printf(" local: %d %d %d %d %d %d\n",lxs,lxe,lnx,lys,lye,lny,lzs,lze,lnz);
  printf("global: %d %d %d %d %d %d\n",gxs,gxe,gnx,gys,gye,gny,gzs,gze,gnz);

  dx = new double[gnx];
  for (int i=0; i<gnx; i++)
    dx[i] = dx_;
  dy = new double[gny];
  for (int j=0; j<gny; j++)
    dy[j] = dy_;
  dz = new double[gnz];
  for (int k=0; k<gnz; k++)
    dz[k] = dz_;

}

void StructuredGrid::setUpConnectivity(int *num_connections, 
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
        vec[0] = 0.5*dx[i]; 
        conn[count].setDistanceUpwind(vec);
        vec[0] = 0.5*dx[i+1]; 
        conn[count].setDistanceDownwind(vec);
        conn[count].setArea(dy[j]*dz[k]);
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
        vec[1] = 0.5*dy[j];
        conn[count].setDistanceUpwind(vec);
        vec[1] = 0.5*dy[j+1];
        conn[count].setDistanceDownwind(vec);
        conn[count].setArea(dx[i]*dz[k]);
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
        vec[2] = 0.5*dz[k];
        conn[count].setDistanceUpwind(vec);
        vec[2] = 0.5*dz[k+1];
        conn[count].setDistanceDownwind(vec);
        conn[count].setArea(dx[i]*dy[j]);
        count++;
      }
    }
  } 

}

void StructuredGrid::mapBoundaryCondition(int istart, int iend, int jstart,
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
        int ighosted = istart+lxs-gxs;
        int jghosted = jstart+lys-gys;
        int kghosted = kstart+lzs-gzs;
        BoundaryCondition *newbc = new BoundaryCondition();
        newbc->setId(i+j*lnx+k*lnxXny);
        dist[0] = 0.; dist[1] = 0.; dist[2] = 0.; 
        normal[0] = 0.; normal[1] = 0.; normal[2] = 0.; 
        // init at center of cell and adjust single coordinate below
//        center[0] = 0.5*dx[ighosted]; center[1] = 0.5*dy[jghosted]; 
//        center[2] = 0.5*dz[kghosted];
        if (!strcmp(face,"west") || !strcmp(face,"east")) {
          newbc->setArea(dy[jghosted]*dz[kghosted]);
          dist[0] = 0.5*dx[ighosted];
          if (!strcmp(face,"west"))
            normal[0] = 1.;
          else if (!strcmp(face,"east"))
            normal[0] = -1.;
        }
        else if (!strcmp(face,"north") || !strcmp(face,"south")) {
          newbc->setArea(dx[ighosted]*dz[kghosted]);
          dist[1] = 0.5*dy[jghosted];
          if (!strcmp(face,"north"))
            normal[1] = -1.;
          else if (!strcmp(face,"south"))
            normal[1] = 1.;
        }
        else if (!strcmp(face,"top") || !strcmp(face,"bottom")) {
          newbc->setArea(dx[ighosted]*dy[jghosted]);
          dist[2] = 0.5*dz[kghosted];
          if (!strcmp(face,"top"))
            normal[2] = -1.;
          else if (!strcmp(face,"bottom"))
            normal[2] = 1.;
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
  for (int k=0; k<gnz; k++) 
    for (int j=0; j<gny; j++) 
      for (int i=0; i<gnx; i++) 
        cells[i].setVolume(dx[i]*dy[j]*dz[k]);
}

void StructuredGrid::setUpMapping(int *num_nodes_local, int *num_nodes_ghosted,
                                  int **mapping_local_to_ghosted,
                                  int **mapping_ghosted_to_local) {

  *num_nodes_local = lnx*lny*lnz;
  *num_nodes_ghosted = gnx*gny*gnz;

  *mapping_local_to_ghosted = new int[*num_nodes_local];
  *mapping_ghosted_to_local = new int[*num_nodes_ghosted];

  for (int i=0; i<*num_nodes_local; i++) (*mapping_local_to_ghosted)[i] = -1;
  for (int i=0; i<*num_nodes_ghosted; i++) (*mapping_ghosted_to_local)[i] = -1;
  
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
        (*mapping_local_to_ghosted)[local_id] = ghosted_id;
        (*mapping_ghosted_to_local)[ghosted_id] = local_id;
        local_id = local_id + 1;
      }
    }
  }

}

void StructuredGrid::get1dofVectorLocal(Vec *v) {
  DACreateLocalVector(da_1dof,v);
}

void StructuredGrid::get1dofVectorGlobal(Vec *v) {
  DACreateGlobalVector(da_1dof,v);
}

void StructuredGrid::getFdofMatrix(Mat *m, MatType mtype) {
  DAGetMatrix(da_fdof,mtype,m);
}

void StructuredGrid::getFdofVectorLocal(Vec *v) {
  DACreateLocalVector(da_fdof,v);
}

void StructuredGrid::getFdofVectorGlobal(Vec *v) {
  DACreateGlobalVector(da_fdof,v);
}

void StructuredGrid::getTdofVectorLocal(Vec *v) {
  DACreateLocalVector(da_fdof,v);
}

void StructuredGrid::getTdofVectorGlobal(Vec *v) {
  DACreateGlobalVector(da_fdof,v);
}

void StructuredGrid::printDACorners() {
  
  PetscErrorCode ierr;
  
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD);
  printf("local         %d lxs:%d lxe:%d lnx:%d\n",myrank,lxs,lxe,lnx);
  printf("local         %d lys:%d lye:%d lny:%d\n",myrank,lys,lye,lny);
  printf("local         %d lzs:%d lze:%d lnz:%d\n",myrank,lzs,lze,lnz);

  printf("local ghosted %d gxs:%d gxe:%d gnx:%d\n",myrank,gxs,gxe,gnx);
  printf("local ghosted %d gys:%d gye:%d gny:%d\n",myrank,gys,gye,gny);
  printf("local ghosted %d gzs:%d gze:%d gnz:%d\n\n",myrank,gzs,gze,gnz);
  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD);
  
}

void StructuredGrid::printAO() {

#include "include/petscao.h"

  PetscErrorCode ierr;
  AO ao;
  ierr = DAGetAO(da_1dof,&ao);
  AOView(ao,PETSC_VIEWER_STDOUT_WORLD);
  AODestroy(ao);

}

StructuredGrid::~StructuredGrid() {
  DADestroy(da_1dof);
  DADestroy(da_fdof);
  DADestroy(da_tdof);
}