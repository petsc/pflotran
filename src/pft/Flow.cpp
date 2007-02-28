#include "Flow.h"

Flow::Flow(int fdof_, Grid *g) {
  
  fdof = fdof_;
  g->getFdofVectorGlobal(&residual_vec);
  g->getFdofMatrix(&Jac,MATMPIAIJ); 
  computeFlow(g,1.0);
  
}

void Flow::computeFlow(Grid *g, double dt) {

  ierr = VecGetArray(g->pressure_vec,&pressure_ptr);
  ierr = VecGetArray(g->ppressure_vec,&ppressure_ptr);
  ierr = VecGetArray(residual_vec,&residual_ptr);

  computeAccumulation(g,dt);
  computeFlux(g);

  ierr = VecRestoreArray(g->pressure_vec,&pressure_ptr);
  ierr = VecRestoreArray(g->ppressure_vec,&ppressure_ptr);
  ierr = VecRestoreArray(residual_vec,&residual_ptr);

  ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);

//  printMatrix();

}

void Flow::computeAccumulation(Grid *g, double dt) {

  double one_over_dt = 1./dt;
  for (int i=0; i<g->num_nodes_local; i++) {
    int offset = i*fdof;
    int inodeghosted = g->cells[i].getIdGhosted();
    double vol_over_dt = g->cells[inodeghosted].getVolume()*one_over_dt;
    double accum = (ppressure_ptr[inodeghosted] - 
                    pressure_ptr[inodeghosted])*one_over_dt;
    residual_ptr[offset] += accum;
    ierr = MatSetValuesLocal(Jac,1,&inodeghosted,1,&inodeghosted,
                             &vol_over_dt,ADD_VALUES);
  }

}

void Flow::computeFlux(Grid *g) {

  for (int i=0; i<g->num_connections; i++) {

    int inodeghostedup = g->connections[i].getIdUpwind();
    int inodeghosteddn = g->connections[i].getIdDownwind();
    double area = g->connections[i].getArea();
    double *distup = g->connections[i].getDistancePtrUpwind();
    double *distdn = g->connections[i].getDistancePtrDownwind();

    double *permup = g->cells[inodeghostedup].getPermPtr();
    double *permdn = g->cells[inodeghosteddn].getPermPtr();
    double muup = g->cells[inodeghostedup].getViscosity();
    double mudn = g->cells[inodeghostedup].getViscosity();
    double densityup = g->cells[inodeghostedup].getDensity();
    double densitydn = g->cells[inodeghostedup].getDensity();

    double dx = distup[0]+distdn[0]; // terribly inefficient for now
    double dy = distup[1]+distdn[1];
    double dz = distup[2]+distdn[2];
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
    if (dx > 0.) coef += permx_ave/dx;
    if (dy > 0.) coef += permy_ave/dy;
    if (dz > 0.) coef += permz_ave/dz;
    coef /= mu_ave;
    // the above is essentially:
    //    coef = -(permx_ave/dx+permy_ave/dy+permz_ave/dz)/mu_ave;


    double vDarcy = -coef*(ppressure_ptr[inodeghosteddn]-
                           ppressure_ptr[inodeghostedup]-
                           gravity*density_ave*dz);
    double qDarcy = vDarcy*area;
    double flux = qDarcy*density_ave;

    double neg_coef = -coef;
    
    // upwind row (downwind flux for upwind cell)
    int inodelocalup = g->cells[inodeghostedup].getIdLocal();
    if (inodelocalup > -1) { // not ghosted
      int offset = inodelocalup*fdof;
      residual_ptr[offset] -= flux;
      // upwind col
      ierr = MatSetValuesLocal(Jac,1,&inodeghostedup,1,&inodeghostedup,
                               &coef,ADD_VALUES);
      // downwind col
      ierr = MatSetValuesLocal(Jac,1,&inodeghostedup,1,&inodeghosteddn,
                               &neg_coef,ADD_VALUES);
    }
    // downwind row (upwind flux for downwind cell)
    int inodelocaldn = g->cells[inodeghosteddn].getIdLocal();
    if (inodelocaldn > -1) { // not ghosted
      int offset = inodelocaldn*fdof;
      residual_ptr[offset] += flux;
      // upwind col
      ierr = MatSetValuesLocal(Jac,1,&inodeghosteddn,1,&inodeghostedup,
                               &neg_coef,ADD_VALUES);
      // downwind col
      ierr = MatSetValuesLocal(Jac,1,&inodeghosteddn,1,&inodeghosteddn,
                               &coef,ADD_VALUES);
    }
  }
}

void Flow::computeBoundaryFlux(Grid *g) {

  BoundaryCondition *cur_bc = g->getBoundaryConditions();

  while (cur_bc) {
    int inodelocal = cur_bc->getId();
    int inodeghosted = g->mapping_local_to_ghosted[inodelocal];
    double area = cur_bc->getArea();
    double *dist = cur_bc->getDistancePtr();
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
    if (dx > 0.) coef += permx/dx;
    if (dy > 0.) coef += permy/dy;
    if (dz > 0.) coef += permz/dz;
    coef /= mu;
    // the above is essentially:
    //    coef = -(permx_ave/dx+permy_ave/dy+permz_ave/dz)/mu_ave;
 
    double pressure_bc = cur_bc->getScalar();

    double vDarcy = -coef*(ppressure_ptr[inodeghosted]-pressure_bc-
                          gravity*density*dz);
    double qDarcy = vDarcy*area;
    double flux = qDarcy*density;

    int offset = inodelocal*fdof;
    residual_ptr[offset] += flux;
    
    ierr = MatSetValuesLocal(Jac,1,&inodeghosted,1,&inodeghosted,
                             &coef,ADD_VALUES);

    cur_bc = cur_bc->getNext();
  }
}

void Flow::printMatrix() {
  MatView(Jac,PETSC_VIEWER_STDOUT_SELF);
}

Flow::~Flow() {
  VecDestroy(residual_vec);
  MatDestroy(Jac);
}


