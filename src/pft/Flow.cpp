#include "Flow.h"

PetscErrorCode dummyFunction(SNES snes, Vec p, Vec f, void *ctx) {
  return ((Flow *)ctx)->computeFlowR(snes,p,f,PETSC_NULL);
}
PetscErrorCode dummyJacobian(SNES snes, Vec x, Mat *J, Mat *B,
                             MatStructure *m, void *ctx) {
  return ((Flow *)ctx)->computeFlowJ(snes,x,J,B,m,PETSC_NULL);
}

int solve_count = 0;

Flow::Flow(int fdof_, Grid *g) {
  
  fdof = fdof_;
  grid = g;
  pressure_ptr = NULL;
  ppressure_ptr = NULL;
  residual_ptr = NULL;
  grid->getFdofVectorGlobal(&work_vec);
  ierr = VecDuplicate(work_vec,&residual_vec);
  grid->getFdofMatrix(&Jac,MATMPIAIJ); 

  PetscOptionsSetValue("-snes_ls","basic");

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);
  ierr = SNESSetFromOptions(snes);
  ierr = SNESSetFunction(snes,residual_vec,dummyFunction,this);
  ierr = SNESSetJacobian(snes,Jac,Jac,dummyJacobian,this);
  ierr = SNESGetKSP(snes,&ksp);
  ierr = KSPGetPC(ksp,&pc);
  ierr = PCSetType(pc,PCBJACOBI);

  init(g);

}

void Flow::init(Grid *g) {
  
  ierr = 0;

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"initial_pressure.out",&viewer);
  ierr = VecView(g->pressure_vec,viewer);
  ierr = PetscViewerDestroy(viewer);

}

void Flow::solve(double dt) {

  solve_count = 0;
  flow_dt = dt;
  grid->zeroFlux();
//  ierr = SNESView(snes,PETSC_VIEWER_STDOUT_WORLD);
  ierr = VecCopy(grid->pressure_vec,work_vec);
  ierr = SNESSolve(snes,PETSC_NULL,work_vec);
//  ierr = SNESView(snes,PETSC_VIEWER_STDOUT_WORLD);
  ierr = VecCopy(work_vec,grid->pressure_vec);

}

PetscErrorCode Flow::computeFlowR(SNES snes, Vec p_vec, Vec f_vec, 
                                  void *ctx) {

  ierr = VecZeroEntries(f_vec);
  ierr = VecGetArray(p_vec,&ppressure_ptr);
  ierr = VecGetArray(grid->pressure_vec,&pressure_ptr);
  ierr = VecGetArray(f_vec,&residual_ptr);

  computeAccumulationR(grid,flow_dt);
  computeFluxR(grid);
  computeBoundaryFluxR(grid);
  CHKMEMQ;

//  for (int i=0; i<grid->num_nodes_local; i++)
//    printf("%d %f\n",i,residual_ptr[i]);

  ierr = VecRestoreArray(p_vec,&ppressure_ptr);
  ierr = VecRestoreArray(grid->pressure_vec,&pressure_ptr);CHKERRQ(ierr);
  ierr = VecRestoreArray(f_vec,&residual_ptr);CHKERRQ(ierr);
  //ierr = VecScale(f_vec,-1.);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"residual.out",&viewer);
  ierr = VecView(f_vec,viewer);
  ierr = PetscViewerDestroy(viewer);

  double norm;
  VecNorm(f_vec,NORM_INFINITY,&norm);
  solve_count++;
  printf("Solve: %d %f\n",solve_count,norm);
 
  return ierr;
}

PetscErrorCode Flow::computeFlowJ(SNES snes, Vec p_vec, Mat *J, Mat *B,
                                  MatStructure *m, void *ctx) {

  MatZeroEntries(Jac);
//  printMatrix();

  computeAccumulationJ(grid,flow_dt);
  computeFluxJ(grid);
  computeBoundaryFluxJ(grid);

  ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);

//  printMatrix();
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"jacobian.out",&viewer);
  ierr = MatView(Jac,viewer);
  ierr = PetscViewerDestroy(viewer);

 // computeNumericalJacobian(grid);

  return ierr;

}

void Flow::computeAccumulationR(Grid *g, double dt) {

  double one_over_dt = 1./dt;
  for (int inodelocal=0; inodelocal < g->num_nodes_local; inodelocal++) {
    int inodeghosted = g->mapping_local_to_ghosted[inodelocal];
    int offset = inodelocal*fdof;
    double accumulation_coef = 
      g->cells[inodeghosted].computeFlowAccumulationCoef(one_over_dt);
    double accumulation_coef_t = accumulation_coef*
                          g->cells[inodeghosted].getSpecificMoistureCapacity_t();
    g->cells[inodeghosted].setFlowAccumulationCoef(accumulation_coef_t);
    for (int idof=0; idof<fdof; idof++) {
      double pp = ppressure_ptr[inodeghosted];
      double p = pressure_ptr[inodeghosted];
      
      if (abs(pp-p) < p_threshhold) {
//        residual_ptr[offset+idof] = (g->cells[inodeghosted].getPressureN_t()-
//                                     g->cells[inodeghosted].getPressure0_t())*
        residual_ptr[offset+idof] = (pp-p)*
                                    accumulation_coef_t;
      }
      else {
        // the below is correct and will always be zero, just using as a place holder
//        residual_ptr[offset+idof] = (g->cells[inodeghosted].getPressureN_t()-
//                                     g->cells[inodeghosted].getPressureN_t())*
//                                    accumulation_coef_t;
//        residual_ptr[offset+idof] -= (g->cells[inodeghosted].getMoistureContentN()-
//                                      g->cells[inodeghosted].getMoistureContent0())*
        residual_ptr[offset+idof] = (pp-p)*
                                     accumulation_coef;
      }
    }
  }
}

void Flow::computeAccumulationJ(Grid *g, double dt) {

  for (int inodelocal=0; inodelocal < g->num_nodes_local; inodelocal++) {
    int inodeghosted = g->mapping_local_to_ghosted[inodelocal];
    double accumulation_coef = 
                       g->cells[inodeghosted].getFlowAccumulationCoef();
    int offset = inodeghosted*fdof;
    for (int idof=0; idof<fdof; idof++) {
      int id = offset+idof;
      ierr = MatSetValuesLocal(Jac,1,&id,1,&id,&accumulation_coef,ADD_VALUES);
    }
  }

}

void Flow::computeFluxR(Grid *g) {

  for (int i=0; i<g->num_connections; i++) {

    int inodeghostedup = g->connections[i].getIdUpwind();
    int inodeghosteddn = g->connections[i].getIdDownwind();
    double area = g->connections[i].getArea();
    double *distup = g->connections[i].getDistancePtrUpwind();
    double *distdn = g->connections[i].getDistancePtrDownwind();
    double *norm = g->connections[i].getNormalPtr();

    double *permup = g->cells[inodeghostedup].getPermPtr();
    double *permdn = g->cells[inodeghosteddn].getPermPtr();
    double muup = g->cells[inodeghostedup].getViscosity();
    double mudn = g->cells[inodeghostedup].getViscosity();
    double densityup = g->cells[inodeghostedup].getDensity();
    double densitydn = g->cells[inodeghostedup].getDensity();

    double dx = abs(distup[0]+distdn[0]); // terribly inefficient for now
    double dy = abs(distup[1]+distdn[1]);
    double dz = abs(distup[2]+distdn[2]);
    double dup = sqrt(distup[0]*distup[0]+distup[1]*distup[1]+
                      distup[2]*distup[2]);
    double ddn = sqrt(distdn[0]*distdn[0]+distdn[1]*distdn[1]+
                      distdn[2]*distdn[2]);
    double dist = dup+ddn;
    
//    printf("%f %f %f %f %f %f\n",dx,dy,dz,dup,ddn,dist);
    
    double density_ave = dist*densityup*densitydn/
                         (dup*densitydn+ddn*densityup);
    double permx_ave = dist*permup[0]*permdn[0]/(dup*permdn[0]+ddn*permup[0]);
    double permy_ave = dist*permup[1]*permdn[1]/(dup*permdn[1]+ddn*permup[1]);
    double permz_ave = dist*permup[2]*permdn[2]/(dup*permdn[2]+ddn*permup[2]);
    double mu_ave = dist*muup*mudn/(dup*mudn+ddn*muup);
    
    double coef = 0.;
    if (dx > 0.) coef += abs(norm[0]*permx_ave/dx);
    if (dy > 0.) coef += abs(norm[1]*permy_ave/dy);
    if (dz > 0.) coef += abs(norm[2]*permz_ave/dz);
    coef /= mu_ave;
    // the above is essentially:
    //    coef = -(permx_ave/dx+permy_ave/dy+permz_ave/dz)/mu_ave;


    double vDarcy = -coef*(ppressure_ptr[inodeghosteddn]-
                           ppressure_ptr[inodeghostedup]+
                           gravity*density_ave*dz);
    double qDarcy = vDarcy*area;
    double flux = qDarcy; //*density_ave;

    g->connections[i].setFlowFluxCoef(coef);

    g->cells[inodeghostedup].addFlux(-1.*qDarcy,norm);
    g->cells[inodeghosteddn].addFlux(qDarcy,norm);

    // upwind row (downwind flux for upwind cell)
    int inodelocalup = g->cells[inodeghostedup].getIdLocal();
    if (inodelocalup > -1) { // not ghosted
      int offset = inodelocalup*fdof;
      for (int idof=0; idof<fdof; idof++) {
        residual_ptr[offset+idof] += flux;
      }
    }
    // downwind row (upwind flux for downwind cell)
    int inodelocaldn = g->cells[inodeghosteddn].getIdLocal();
    if (inodelocaldn > -1) { // not ghosted
      int offset = inodelocaldn*fdof;
      for (int idof=0; idof<fdof; idof++) {
        residual_ptr[offset+idof] -= flux;
      }
    }
  }
}

void Flow::computeFluxJ(Grid *g) {

  for (int i=0; i<g->num_connections; i++) {

    int inodeghostedup = g->connections[i].getIdUpwind();
    int inodeghosteddn = g->connections[i].getIdDownwind();

    double coef = g->connections[i].getFlowFluxCoef();
    double neg_coef = -coef;

    int offsetup = inodeghostedup*fdof;
    int offsetdn = inodeghosteddn*fdof;
    // upwind row (downwind flux for upwind cell)
    int inodelocalup = g->mapping_ghosted_to_local[inodeghostedup];
    if (inodelocalup > -1) { // not ghosted
      for (int idof=0; idof<fdof; idof++) {
        int idup = offsetup+idof;
        int iddn = offsetdn+idof;
        // upwind col
        ierr = MatSetValuesLocal(Jac,1,&idup,1,&idup,&coef,ADD_VALUES);
        // downwind col
        ierr = MatSetValuesLocal(Jac,1,&idup,1,&iddn,&neg_coef,ADD_VALUES);
      }
    }
    // downwind row (upwind flux for downwind cell)
    int inodelocaldn = g->mapping_ghosted_to_local[inodeghosteddn];
    if (inodelocaldn > -1) { // not ghosted
      for (int idof=0; idof<fdof; idof++) {
        int idup = offsetup+idof;
        int iddn = offsetdn+idof;
        // upwind col
        ierr = MatSetValuesLocal(Jac,1,&iddn,1,&idup,&neg_coef,ADD_VALUES);
        // downwind col
        ierr = MatSetValuesLocal(Jac,1,&iddn,1,&iddn,&coef,ADD_VALUES);
      }
    }
  }
}

void Flow::computeBoundaryFluxR(Grid *g) {

  BoundaryCondition *cur_bc = g->getBoundaryConditions();

  while (cur_bc) {
    int inodelocal = cur_bc->getId();
    int inodeghosted = g->mapping_local_to_ghosted[inodelocal];
    double area = cur_bc->getArea();
    double *dist = cur_bc->getDistancePtr();
    double *norm = cur_bc->getNormalPtr();
    double *perm = g->cells[inodeghosted].getPermPtr();
    double mu = g->cells[inodeghosted].getViscosity();
    double density = g->cells[inodeghosted].getDensity();

    double dx = dist[0]; // terribly inefficient for now
    double dy = dist[1];
    double dz = dist[2];
    double distance = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
    
//    printf("%f %f %f %f %f %f\n",dx,dy,dz,dup,ddn,dist);
    
    double permx = perm[0];
    double permy = perm[1];
    double permz = perm[2];
    
    double coef = 0.;
    if (abs(dx) > 0.) coef += abs(norm[0]*permx/dx);
    if (abs(dy) > 0.) coef += abs(norm[1]*permy/dy);
    if (abs(dz) > 0.) coef += abs(norm[2]*permz/dz);
    coef /= mu;
    // the above is essentially:
    //    coef = -(permx_ave/dx+permy_ave/dy+permz_ave/dz)/mu_ave;
 
    cur_bc->setFlowFluxCoef(coef);

    double pressure_bc = cur_bc->getScalar();

    double vDarcy = -coef*(ppressure_ptr[inodeghosted]-pressure_bc+
                           gravity*density*dz);
    double qDarcy = vDarcy*area;
    double flux = qDarcy; //*density;

    g->cells[inodeghosted].addFlux(qDarcy,norm);
 
    int offset = inodelocal*fdof;
    for (int idof=0; idof<fdof; idof++) {
      residual_ptr[offset+idof] -= flux;
    }
    
    cur_bc = cur_bc->getNext();
  }
}

void Flow::computeBoundaryFluxJ(Grid *g) {

  BoundaryCondition *cur_bc = g->getBoundaryConditions();

  while (cur_bc) {

    double coef = cur_bc->getFlowFluxCoef();

    int inodeghosted = g->mapping_local_to_ghosted[cur_bc->getId()];
    int offset = inodeghosted*fdof;
    for (int idof=0; idof<fdof; idof++) {
      int id = offset+idof;
      ierr = MatSetValuesLocal(Jac,1,&id,1,&id,&coef,ADD_VALUES);
    }
    cur_bc = cur_bc->getNext();
  }

}

void Flow::computeSourceFluxR(Grid *g) {

  Source *cur_src = g->getSources();

  while (cur_src) {
    cur_src = cur_src->getNext();
  }
}

void Flow::computeNumericalJacobian(Grid *g) {

  double *temp_ptr = NULL;

  Vec temp_f_baseline_vec;
  Vec temp_p_vec;
  Vec temp_f_vec;
  Mat nJac;

  g->getFdofVectorGlobal(&temp_p_vec);
  g->getFdofMatrix(&nJac,MATMPIAIJ); 

  ierr = VecDuplicate(temp_p_vec,&temp_f_baseline_vec);
  ierr = VecDuplicate(temp_p_vec,&temp_f_vec);

  computeFlowR(PETSC_NULL,grid->pressure_vec,temp_f_baseline_vec,NULL);

  double tol = 1.e-6;
  for (int inodelocal=0; inodelocal < g->num_nodes_local; inodelocal++) {
    ierr = VecCopy(grid->pressure_vec,temp_p_vec);
    ierr = VecGetArray(temp_p_vec,&temp_ptr);

    double perturbation;
    if (temp_ptr[inodelocal] > 0.)
      perturbation = temp_ptr[inodelocal]*tol;
    else
      perturbation = 1.;
    temp_ptr[inodelocal] += perturbation;

    ierr = VecRestoreArray(temp_p_vec,&temp_ptr);

    computeFlowR(PETSC_NULL,temp_p_vec,temp_f_vec,NULL);

    ierr = VecAXPY(temp_f_vec,-1.,temp_f_baseline_vec);
    ierr = VecScale(temp_f_vec,1./perturbation);
    ierr = VecGetArray(temp_f_vec,&temp_ptr);

    for (int j=0; j<g->num_nodes_local; j++) {
      if (abs(temp_ptr[j]) > 0.)
        ierr = MatSetValuesLocal(nJac,1,&j,1,&inodelocal,&(temp_ptr[j]),INSERT_VALUES);
    }

    ierr = VecRestoreArray(temp_f_vec,&temp_ptr);
  }

  ierr = MatAssemblyBegin(nJac,MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(nJac,MAT_FINAL_ASSEMBLY);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"nJac.out",&viewer);
  ierr = MatView(nJac,viewer);
  ierr = PetscViewerDestroy(viewer);

  ierr = MatDestroy(nJac);
  ierr = VecDestroy(temp_f_baseline_vec);
  ierr = VecDestroy(temp_p_vec);
  ierr = VecDestroy(temp_f_vec);
}

void Flow::printMatrix() {
  MatView(Jac,PETSC_VIEWER_STDOUT_SELF);
}

Flow::~Flow() {
  VecDestroy(work_vec);
  MatDestroy(Jac);
  SNESDestroy(snes);
}


