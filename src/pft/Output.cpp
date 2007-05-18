#include "Output.h"

Output::Output() {

  PetscErrorCode ierr;

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"pressure.dat",&pressure_handle);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"flux.dat",&flux_handle);

}

void Output::printGMSGrid(Grid *g) {

  PetscErrorCode ierr;
  PetscViewer viewer = NULL;

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"grid.3dg",&viewer);
  ierr = PetscViewerASCIIPrintf(viewer,"GRID3D\n");
  ierr = PetscViewerASCIIPrintf(viewer,"TYPE 1\n");
  ierr = PetscViewerASCIIPrintf(viewer,"IJK +y +x +z\n");
  ierr = PetscViewerASCIIPrintf(viewer,"NUMBERING 0\n");
  ierr = PetscViewerASCIIPrintf(viewer,"Origin 0. 0. 0.\n");
  ierr = PetscViewerASCIIPrintf(viewer,"ROTZ 0.\n");
  ierr = PetscViewerASCIIPrintf(viewer,"DIM %d %d %d\n",g->getNx()+1,g->getNy()+1,g->getNz()+1);
  int nx = g->getNx();
  int ny = g->getNy();
  int nz = g->getNz();
  double sum = 0.;
  for (int i=0; i<=nx; i++) {
    ierr = PetscViewerASCIIPrintf(viewer,"%f\n",sum);
    if (i < nx) sum += g->getDx(i);
  }
  sum = 0.;
  for (int j=0; j<=ny; j++) {
    ierr = PetscViewerASCIIPrintf(viewer,"%f\n",sum);
    if (j < ny) sum += g->getDy(j);
  }
  sum = 0.;
  for (int k=0; k<=nz; k++) {
    ierr = PetscViewerASCIIPrintf(viewer,"%f\n",sum);
    if (k < nz) sum += g->getDz(k);
  }
//  ierr = PetscViewerASCIIPrintf(viewer,"MAT\n");
//  ierr = PetscViewerASCIIPrintf(viewer,"ACTIVE\n");
  ierr = PetscViewerDestroy(viewer);


  ierr = PetscViewerASCIIPrintf(pressure_handle,"DATASET\n");
  ierr = PetscViewerASCIIPrintf(pressure_handle,"OBJTYPE \"grid3d\"\n");
  ierr = PetscViewerASCIIPrintf(pressure_handle,"BEGSCL\n");
  ierr = PetscViewerASCIIPrintf(pressure_handle,"ND %d\n",g->getN());
  ierr = PetscViewerASCIIPrintf(pressure_handle,"NC %d\n",g->getN());
  ierr = PetscViewerASCIIPrintf(pressure_handle,"Name \"pressure\"\n");

  ierr = PetscViewerASCIIPrintf(flux_handle,"DATASET\n");
  ierr = PetscViewerASCIIPrintf(flux_handle,"OBJTYPE \"grid3d\"\n");
  ierr = PetscViewerASCIIPrintf(flux_handle,"BEGVEC\n");
  ierr = PetscViewerASCIIPrintf(flux_handle,"VECTYPE 1\n");
  ierr = PetscViewerASCIIPrintf(flux_handle,"ND %d\n",g->getN());
  ierr = PetscViewerASCIIPrintf(flux_handle,"NC %d\n",g->getN());
  ierr = PetscViewerASCIIPrintf(flux_handle,"Name \"Darcy flux\"\n");

}

void Output::printGMSOutput(Grid *g, double time) {
  printGMSPressure(g,time);
  printGMSFlux(g,time);
}

void Output::printGMSPressure(Grid *g, double time) {

  PetscErrorCode ierr;
  Vec natural_vec = NULL;
  PetscScalar *vec_ptr = NULL;

  g->get1dofVectorNatural(&natural_vec);
  g->globalToNatural(g->pressure_vec,natural_vec);
  ierr = PetscViewerASCIIPrintf(pressure_handle,"TS 0 %f\n", time);
  ierr = VecGetArray(natural_vec,&vec_ptr);
  for (int i=0; i<g->getN(); i++) {
    ierr = PetscViewerASCIIPrintf(pressure_handle,"%f\n",vec_ptr[i]);
  }
  ierr = VecRestoreArray(natural_vec,&vec_ptr);
  VecDestroy(natural_vec);

}

void Output::printGMSFlux(Grid *g, double time) {

  PetscErrorCode ierr;
  double *double_ptr = NULL;

  ierr = PetscViewerASCIIPrintf(flux_handle,"TS 0 %f\n", time);
  for (int i=0; i<g->getN(); i++) {
    double_ptr = g->cells[i].getFlux();
    ierr = PetscViewerASCIIPrintf(flux_handle,"%e %e %e\n",
                                  double_ptr[0],double_ptr[1],double_ptr[2]);
  }

}

Output::~Output() {

  PetscErrorCode ierr;

  ierr = PetscViewerASCIIPrintf(pressure_handle,"ENDDS\n");
  ierr = PetscViewerDestroy(pressure_handle);

  ierr = PetscViewerASCIIPrintf(flux_handle,"ENDDS\n");
  ierr = PetscViewerDestroy(flux_handle);

}
