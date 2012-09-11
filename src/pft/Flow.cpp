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

  init(grid);

}

void Flow::init(Grid *g) {
  
  ierr = 0;

  // parameters from unit 3 in Socorro EMSP site
  double theta_res = porosity0*0.1825;
  double theta_sat = porosity0;
  double alpha = 0.0405/(density0*gravity); // 1/Pa
  double n = 5.1330;
  g->initSaturationFunction(theta_res,theta_sat,alpha,n);

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

  grid->updateRelativePermeability(ppressure_ptr);

#if 0
  computeAccumulationR(grid,flow_dt);
  computeFluxR(grid);
  computeBoundaryFluxR(grid);
#else
//  computeAccumulationRVSat(grid,flow_dt);
  computeFluxRVSat(grid);
  computeBoundaryFluxRVSat(grid);
#endif
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
#if 0
  computeAccumulationJ(grid);
  computeFluxJ(grid);
  computeBoundaryFluxJ(grid);
#else
//  computeAccumulationJVSat(grid);
  computeFluxJVSat(grid);
  computeBoundaryFluxJVSat(grid);
#endif
  ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);

//  printMatrix();
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"jacobian.out",&viewer);
  ierr = MatView(Jac,viewer);
  ierr = PetscViewerDestroy(viewer);

  computeNumericalJacobian(grid);

  return ierr;

}

void Flow::computeAccumulationR(Grid *g, double dt) {

  double one_over_dt = 1./dt;
  for (int inodelocal=0; inodelocal < g->num_nodes_local; inodelocal++) {
    int inodeghosted = g->mapping_local_to_ghosted[inodelocal];
    int offset = inodelocal*fdof;
    double accumulation_coef = 
      g->cells[inodeghosted].computeFlowAccumulationCoef(one_over_dt);
    g->cells[inodeghosted].setFlowAccumulationCoef(accumulation_coef);
    for (int idof=0; idof<fdof; idof++) {
      double pp = ppressure_ptr[inodeghosted];
      double p = pressure_ptr[inodeghosted];
      residual_ptr[offset+idof] = (pp-p)*accumulation_coef;
    }
  }
}

void Flow::computeAccumulationJ(Grid *g) {

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

  //  g->connections[i].setFlowFluxCoef(coef);
    g->connections[i].setFlowFluxCoefVSat(coef,-coef);

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

  //  double coef = g->connections[i].getFlowFluxCoef();
  //  double neg_coef = -coef;
    double coef;
    double neg_coef;
    g->connections[i].getFlowFluxCoefVSat(&coef,&neg_coef);

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

  computeFlowR(PETSC_NULL,g->pressure_vec,temp_f_baseline_vec,NULL);

  double tol = 1.e-6;
  for (int inodelocal=0; inodelocal < g->num_nodes_local; inodelocal++) {
    ierr = VecCopy(g->pressure_vec,temp_p_vec);
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

double Flow::convertP(double p) {
  return p/(1.+betap*p);
}

void Flow::computeAccumulationRVSat(Grid *g, double dt) {

  double one_over_dt = 1./dt;
  for (int inodelocal=0; inodelocal < g->num_nodes_local; inodelocal++) {
    int inodeghosted = g->mapping_local_to_ghosted[inodelocal];
    int offset = inodelocal*fdof;
    double accumulation_coef = 
      g->cells[inodeghosted].computeFlowAccumulationCoef(one_over_dt);
    double spec_moist_cap = g->cells[inodeghosted].getSpecificMoistureCapacity();
    double pp = ppressure_ptr[inodeghosted];
    double p = pressure_ptr[inodeghosted];
    double ppt = convertP(pp);
    double pt = convertP(p);
    double one_plus_beta_p = 1.+betap*p;
    double one_plus_beta_p_sq = one_plus_beta_p*one_plus_beta_p;
    double one_over_one_plus_beta_p_sq = 1./one_plus_beta_p_sq;
    double deriv_one_plus_beta_p_sq = 2.*betap*(1+betap*p);
    double deriv = 0.;
    for (int idof=0; idof<fdof; idof++) {
  //    if (abs(pp-p) < p_threshhold) {
      if (1) {
        residual_ptr[offset+idof] = (ppt-pt)*spec_moist_cap*one_plus_beta_p_sq*
                                    accumulation_coef;
        deriv += (ppt-pt)*g->cells[inodeghosted].getDerivSpecMoistCap()*
                 one_plus_beta_p_sq*accumulation_coef;
        deriv += (ppt-pt)*spec_moist_cap*
                 deriv_one_plus_beta_p_sq*accumulation_coef;
        deriv += one_over_one_plus_beta_p_sq*spec_moist_cap*
                 one_plus_beta_p_sq*accumulation_coef;
      }
      else {
        // the below is correct and will always be zero, just using as a place holder
        double pnt = convertP(g->cells[inodeghosted].getPressureN());
        residual_ptr[offset+idof] = (ppt-pnt)*spec_moist_cap*one_plus_beta_p_sq*
                                    accumulation_coef;
        residual_ptr[offset+idof] += (g->cells[inodeghosted].getMoistureContentN()-
                                     g->cells[inodeghosted].getMoistureContent0())*
                                     accumulation_coef;
        deriv += (ppt-pnt)*g->cells[inodeghosted].getDerivSpecMoistCap()*
                 one_plus_beta_p_sq*accumulation_coef;
        deriv += (ppt-pnt)*spec_moist_cap*
                 deriv_one_plus_beta_p_sq*accumulation_coef;
        deriv += one_over_one_plus_beta_p_sq*spec_moist_cap*
                 one_plus_beta_p_sq*accumulation_coef;
      }
    }
    g->cells[inodeghosted].setFlowAccumulationCoef(deriv);
  }
}

void Flow::computeAccumulationJVSat(Grid *g) {

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

void Flow::computeFluxRVSat(Grid *g) {

  for (int i=0; i<g->num_connections; i++) {

    int inodeghostedup = g->connections[i].getIdUpwind();
    int inodeghosteddn = g->connections[i].getIdDownwind();
    double area = g->connections[i].getArea();
    double *distup = g->connections[i].getDistancePtrUpwind();
    double *distdn = g->connections[i].getDistancePtrDownwind();
    double *norm = g->connections[i].getNormalPtr();

    double *permup = g->cells[inodeghostedup].getPermPtr();
    double *permdn = g->cells[inodeghosteddn].getPermPtr();
    double relpermup = g->cells[inodeghostedup].getRPerm();
    double relpermdn = g->cells[inodeghosteddn].getRPerm();
    double drelpermup = g->cells[inodeghostedup].getDRelPerm();
    double drelpermdn = g->cells[inodeghosteddn].getDRelPerm();
    double muup = g->cells[inodeghostedup].getViscosity();
    double mudn = g->cells[inodeghosteddn].getViscosity();
    double densityup = g->cells[inodeghostedup].getDensity();
    double densitydn = g->cells[inodeghosteddn].getDensity();

    double dx = abs(distup[0]+distdn[0]); // terribly inefficient for now
    double dy = abs(distup[1]+distdn[1]);
    double dz = abs(distup[2]+distdn[2]);
    double dup = sqrt(distup[0]*distup[0]+distup[1]*distup[1]+
                      distup[2]*distup[2]);
    double ddn = sqrt(distdn[0]*distdn[0]+distdn[1]*distdn[1]+
                      distdn[2]*distdn[2]);
    double distance = dup+ddn;
    
//    printf("%f %f %f %f %f %f\n",dx,dy,dz,dup,ddn,dist);
    
    double density_ave = distance*densityup*densitydn/
                         (dup*densitydn+ddn*densityup);
    double mu_ave = distance*muup*mudn/(dup*mudn+ddn*muup);
    
    double permup_ = abs(norm[0]*permup[0])+ abs(norm[1]*permup[1])
                   + abs(norm[2]*permup[2]);
    double permdn_ = abs(norm[0]*permdn[0])+ abs(norm[1]*permdn[1])
                   + abs(norm[2]*permdn[2]);

    permup_ /= mu_ave;
    permdn_ /= mu_ave;

    double pup = ppressure_ptr[inodeghostedup];
    double pdn = ppressure_ptr[inodeghosteddn];

    double betaup = 0.;//pup < patm ? betap : 0.;
    double betadn = 0.;//pdn < patm ? betap : 0.;

    double A1 = 0.5*permup_*relpermup*(1.+betaup*pup)/(1.+betadn*pdn);
    double A2 = -0.5*permdn_*relpermdn*(1.+betadn*pdn)/(1.+betaup*pup);
    double A = A1+A2;
    double B = (pdn-pup)/distance;
    double C = (1.+pup*pdn*(betaup-betadn))/(pdn-pup);

    // not sure about the mult by norm[2]....
    double vDarcy = A*B*C-0.5*(permup_*relpermup+permdn_*relpermdn)*
                    density_ave*gravity*norm[2];

    double qDarcy = vDarcy*area;
    double flux = qDarcy; //*density_ave;

    double dA1_dpup = 0.5*permup_*relpermup*betaup/(1+betadn*pdn);
    double dA2_dpup = 0.5*permdn_*relpermdn*betaup*(1.+betadn*pdn)/
                          ((1.+betaup*pup)*(1.+betaup*pup));
    double dA1_dkup = 0.5*permup_*drelpermup*(1.+betaup*pup)/(1.+betadn*pdn);
    double dA1_dpdn = -0.5*permup_*relpermup*betadn*(1.+betaup*pup)/
                          ((1.+betadn*pdn)*(1.+betadn*pdn));
    double dA2_dpdn = -0.5*permdn_*relpermdn*betadn/(1.+betaup*pup);
    double dA2_dkdn = -0.5*permdn_*drelpermdn*(1.+betadn*pdn)/(1.+betaup*pup);
    double dB_dpup = -1./distance;
    double dB_dpdn = 1./distance;
    double dC_dpup = pdn*(betaup-betadn)/(pdn-pup)+
                     (1.+pup*pdn*(betaup-betadn))/((pdn-pup)*(pdn-pup));
    double dC_dpdn = pup*(betaup-betadn)/(pdn-pup)-
                     (1.+pup*pdn*(betaup-betadn))/((pdn-pup)*(pdn-pup));

    double dR_dup = (dA1_dpup+dA2_dpup+dA1_dkup)*B*C+
                    A*dB_dpup*C+
                    A*B*dC_dpup-
                    0.5*permup_*drelpermup*density_ave*gravity*norm[2];
    double dR_ddn = (dA1_dpdn+dA2_dpdn+dA2_dkdn)*B*C+
                    A*dB_dpdn*C+
                    A*B*dC_dpdn-
                    0.5*permdn_*drelpermdn*density_ave*gravity*norm[2];

    g->connections[i].setFlowFluxCoefVSat(-dR_dup,-dR_ddn);

    g->cells[inodeghostedup].addFlux(-1.*qDarcy,norm);
    g->cells[inodeghosteddn].addFlux(qDarcy,norm);

    // upwind row (downwind flux for upwind cell)
    int inodelocalup = g->cells[inodeghostedup].getIdLocal();
    if (inodelocalup > -1) { // not ghosted
      int offset = inodelocalup*fdof;
      for (int idof=0; idof<fdof; idof++) {
        residual_ptr[offset+idof] -= flux;
      }
    }
    // downwind row (upwind flux for downwind cell)
    int inodelocaldn = g->cells[inodeghosteddn].getIdLocal();
    if (inodelocaldn > -1) { // not ghosted
      int offset = inodelocaldn*fdof;
      for (int idof=0; idof<fdof; idof++) {
        residual_ptr[offset+idof] += flux;
      }
    }
  }
}

// technically, the first implementation above should be adequate since
// there is essentially no difference between the two
void Flow::computeFluxJVSat(Grid *g) {

  for (int i=0; i<g->num_connections; i++) {

    int inodeghostedup = g->connections[i].getIdUpwind();
    int inodeghosteddn = g->connections[i].getIdDownwind();

    double dR_dup;
    double dR_ddn;
    g->connections[i].getFlowFluxCoefVSat(&dR_dup,&dR_ddn);

    int offsetup = inodeghostedup*fdof;
    int offsetdn = inodeghosteddn*fdof;
    // upwind row (downwind flux for upwind cell)
    int inodelocalup = g->mapping_ghosted_to_local[inodeghostedup];
    if (inodelocalup > -1) { // not ghosted
      for (int idof=0; idof<fdof; idof++) {
        int idup = offsetup+idof;
        int iddn = offsetdn+idof;
        // upwind col
        ierr = MatSetValuesLocal(Jac,1,&idup,1,&idup,&dR_dup,ADD_VALUES);
        // downwind col
        ierr = MatSetValuesLocal(Jac,1,&idup,1,&iddn,&dR_ddn,ADD_VALUES);
      }
    }
    // downwind row (upwind flux for downwind cell)
    int inodelocaldn = g->mapping_ghosted_to_local[inodeghosteddn];
    if (inodelocaldn > -1) { // not ghosted
      for (int idof=0; idof<fdof; idof++) {
        int idup = offsetup+idof;
        int iddn = offsetdn+idof;
        // upwind col
        ierr = MatSetValuesLocal(Jac,1,&iddn,1,&idup,&dR_ddn,ADD_VALUES);
        // downwind col
        ierr = MatSetValuesLocal(Jac,1,&iddn,1,&iddn,&dR_dup,ADD_VALUES);
      }
    }
  }
}

void Flow::computeBoundaryFluxRVSat(Grid *g) {

  BoundaryCondition *cur_bc = g->getBoundaryConditions();

  while (cur_bc) {
    int inodelocal = cur_bc->getId();
    int inodeghosted = g->mapping_local_to_ghosted[inodelocal];
    double area = cur_bc->getArea();
    double *dist = cur_bc->getDistancePtr();
    double *norm = cur_bc->getNormalPtr();
    double *perm = g->cells[inodeghosted].getPermPtr();
    double relperm = g->cells[inodeghosted].getRPerm();
    double drelperm = g->cells[inodeghosted].getDRelPerm();
 
    double mu = g->cells[inodeghosted].getViscosity();
    double density = g->cells[inodeghosted].getDensity();

    double dx = dist[0]; // terribly inefficient for now
    double dy = dist[1];
    double dz = dist[2];
    double distance = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
    
//    printf("%f %f %f %f %f %f\n",dx,dy,dz,dup,ddn,dist);
    
    double perm_ = abs(norm[0]*perm[0])+ abs(norm[1]*perm[1])
                 + abs(norm[2]*perm[2]);

    perm_ /= mu;
 
    double pup = cur_bc->getScalar(); // bc pressure
    double pdn = ppressure_ptr[inodeghosted];

    double betaup = pup < patm ? betap : 0.;
    double betadn = pdn < patm ? betap : 0.;

    // since no averaging of perm takes place, perm_ is overkill,
    // but left alone to match the internodal flux calcs.
  //  double A = 0.5*(perm_*relperm*(1.+betaup*pup)/(1.+betadn*pdn)-
  //                  perm_*relperm*(1.+betadn*pdn)/(1.+betaup*pup));
    double A = perm_*relperm;
    double B = (pdn-pup)/distance;
    double C = (1.+pup*pdn*(betaup-betadn))/(pdn-pup);

    // not sure about the mult by norm[2]....
    double vDarcy = A*B*C-(perm_*relperm)*density*gravity*norm[2];

    double qDarcy = vDarcy*area;
    double flux = qDarcy; 

    // dup stuff not needed since it is the bc pressure
//    double dA_dpup = 0.5*(permup_*betaup/(1+betadn*pdn)+
//                          permdn_*betadn*(1.+betadn*pdn)/
//                          ((1.+betaup*pup)*(1.+betaup*pup)));
/*    double dA_dpdn = 0.5*(-perm_*relperm*betadn*(1.+betaup*pup)/
                          ((1+betadn*pdn)*(1+betadn*pdn))-
                           perm_*relperm*betadn/(1.+betaup*pup));
    double dA_dkdn = 0.5*(perm_*drelperm*(1.+betaup*pup)/(1.+betadn*pdn)-
                          perm_*drelperm*(1.+betadn*pdn)/(1.+betaup*pup));*/
    // since the perm is within the cell, no weighting.
    double dA_dkdn = perm_*drelperm;
//    double dB_dpup = -1./dist;
    double dB_dpdn = 1./distance;
//    double dC_dpup = pdn*(betaup-betadn)/(pdn-pup)+
//                     (1.+pup*pdn*(betaup-betadn))/((pdn-pup)*(pdn-pup));
    double dC_dpdn = pup*(betaup-betadn)/(pdn-pup)-
                     (1.+pup*pdn*(betaup-betadn))/((pdn-pup)*(pdn-pup));

//    double dR_dup = dA_dpup*B*C+A*dB_dpup*C+A*B*dC_dpup;
//    double dR_ddn = (dA_dpdn+dA_dkdn)*B*C+
    double dR_ddn = dA_dkdn*B*C+
                    A*dB_dpdn*C+
                    A*B*dC_dpdn-
                    perm_*drelperm*density*gravity*norm[2];

    cur_bc->setFlowFluxCoef(dR_ddn);

    g->cells[inodeghosted].addFlux(qDarcy,norm);
 
    int offset = inodelocal*fdof;
    for (int idof=0; idof<fdof; idof++) {
      residual_ptr[offset+idof] -= flux;
    }
    
    cur_bc = cur_bc->getNext();
  }
}

void Flow::computeBoundaryFluxJVSat(Grid *g) {

  BoundaryCondition *cur_bc = g->getBoundaryConditions();

  while (cur_bc) {

    double dR_ddn = cur_bc->getFlowFluxCoef();

    int inodeghosted = g->mapping_local_to_ghosted[cur_bc->getId()];
    int offset = inodeghosted*fdof;
    for (int idof=0; idof<fdof; idof++) {
      int id = offset+idof;
      ierr = MatSetValuesLocal(Jac,1,&id,1,&id,&dR_ddn,ADD_VALUES);
    }
    cur_bc = cur_bc->getNext();
  }

}

Flow::~Flow() {
  VecDestroy(work_vec);
  MatDestroy(Jac);
  SNESDestroy(snes);
}


