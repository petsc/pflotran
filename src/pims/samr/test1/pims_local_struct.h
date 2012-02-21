#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscda.h"

extern "C"{
struct pflow_local_patch_info{

  
   /* Local quantities*/
  
   
   int nlx, nly, nlz; /* Local grid dimension w/o ghost nodes.*/
   int ngx, ngy, ngz;  /* Local grid dimension with ghost nodes.*/
   int nxs, nys, nzs; 
   /* Global indices of non-ghosted corner (starting) of local domain.*/
   int ngxs, ngys, ngzs; 
   /* Global indices of ghosted starting corner of local domain.*/
   int nxe, nye, nze, ngxe, ngye, ngze; 
   /* Global indices of non-ghosted/ghosted ending corner of local domain.*/
   int nlxy, nlxz, nlyz; 
   int ngxy, ngxz, ngyz; 
   int nlmax;   /* Total number of non-ghosted nodes in local domain.*/
   int ngmax;   /* Number of ghosted & non-ghosted nodes in local domain.*/
   int nldof;   /* nlmax times the number of phases.*/
   int ngdof;   /* ngmax times the number of phases.*/
   int istart, jstart, kstart, iend, jend, kend; 
   /* istart gives the local x-index of the non-ghosted starting (lower left)*/
   /* corner. iend gives the local x-index of the non-ghosted ending */
   /* corner. jstart, jend correspond to y-index, kstart, kend to z-index.*/
   

   /* Grid connections*/
   int nconn, nconnx, nconny; 
   int *nd1;
   int *nd2;
   /* Nodes upstream and downstream of a connection (assuming flow in */
   /* positive direction.  These are local, ghosted indices.*/
   
   int *iperm1;
   int *iperm2;
   int *ipermbc;
      
   double *dist1;
   double *dist2;
   double *distbc;
   double *area;
   double *areabc;
   double *delzbc;
   double *vlbc;
   double *vvlbc;
   double *vgbc;
   double *vvgbc;
   double* delz; 


   int *nL2G;
   int *nG2L;
   int *nL2A;
   int *nG2N;
   
   /* Boundary conditions (BC's)*/
   int nblkbc;
   /* The number of "blocks" of boundary conditions that are defined.*/
   /* Such a block is a specification of a set of boundary conditions.*/
   /* This set of boundary conditions can apply to any number of regions,*/
   /* so nblkbc does NOT equal the number of boundary condition regions.*/
   int nconnbc;  /* The number of interfaces along boundaries.*/
   int *ibconn;
   /* ibconn(nc) specifies the index of the boundary condition block that*/
   /* applies at boundary interface nc.  */
   int *ibndtyp;
   /* ibndtyp(ibc) specifies the type of boundary condition that applies*/
   /* for boundary condition block ibc.*/
   int *iface;
   /* iface(ibc) specifies the face (left, right, top, bottom, etc.) on*/
   /* which BC block ibc lies.*/
   int *mblkbc;
   /* mblkbc(nc) gives the local, non-ghosted index of the cell that has*/
   /* boundary connection nc.*/
   double *pressurebc;
   /* For a Dirichlet BC, pressurebc(j,ibc) gives the partial pressure */
   /* for phase j along the BC block ibc.*/
   double *velocitybc;
   /* For a Neumann BC, velocitybc(j,ibc) gives the velocity q for phase*/
   /* j along BC block ibc.*/
   double *tempbc;
   double *xxbc;
   double *varbc;
   
   /*   double*vl_loc(:), vvl_loc(:), vg_loc(:), vvg_loc(:)*/
   double *rtot;
   double *rate;
   double *area_var;
   double *delx;


   double *var;
   PetscScalar *accum_p;

   PetscScalar *r_p;
   PetscScalar *xx_loc_p;
   PetscScalar *xx_p;
   PetscScalar *yy_p;
   PetscScalar *porosity_loc_p;
   PetscScalar *volume_p;
   PetscScalar *phis_p;
   PetscScalar *tor_loc_p;
   PetscScalar *perm_xx_loc_p;
   PetscScalar *perm_yy_loc_p;
   PetscScalar *perm_zz_loc_p;
   PetscScalar *vl_p;
                          
               
   PetscScalar *pc_p;
   PetscScalar *pc_loc_p;
   PetscScalar *kvr_p;
   PetscScalar *kvr_loc_p;
   
   PetscScalar *icap_p;
   PetscScalar *icap_loc_p;
   PetscScalar *ithrm_loc_p;
   PetscScalar *ithrm_p;


};

struct pflowGrid{

/* Note that preprocessor directives MUST start in the first column!*/
/*#ifndef DEBUG*/
/*   private*/
/*#endif*/

   int myrank, commsize;  /* Rank in PETSC_COMM_WORLD.*/
   int npx, npy, npz; /* Processor partition in each direction.*/
   int nxy, nmax;     /* nx * ny, nx * ny * nz*/
   int nphase, nvar, ndof;  /* Number of phases we are dealing with.*/
   int size_var_use, size_var_node;
   int jh2o, jgas, joil; /* specific phase indices*/


    /* Program options*/
   PetscTruth  use_analytical;  /* If true, use analytical Jacobian.*/
   PetscTruth  use_matrix_free;  /* If true, do not form the Jacobian.*/
   /* Note that if 'use_analytical' and 'use_matrix_free' are both false,*/
   /* the Jacobian will be computed numerically and stored.*/
   PetscTruth  print_hhistory;
   /* If true, and if use_matrix_free is true, then store the differencing*/
   /* values h and print them out at the end of the simulation.*/
   
   PetscScalar *hhistory;
   PetscTruth  monitor_h;
   /* If true, print the value of h at the end of each SNES iteration.*/
   PetscTruth  use_ksp;
   PetscTruth  use_isoth, use_debug;	
   /* If using_pflowGrid == PETSC_TRUE, then some parts of ptran_init */
   /* will not be executed, since they are made redundant by */
   /* pflowGrid_new() and pflowGrid_setup().*/
   
   double t;  /* The time elapsed in the simulation.*/
   double dt; /* The size of the time step.*/
   double dt_min;  /* Maximum size of the time step.*/
   double dt_max;  /* Maximum size of the time step.*/
   double tconv; /* Input time conversion factor*/
   char tunit; /* Input time units*/
   double *tplot;
   double *tstep;
   double *dtstep;
   double *tfac;
   /* An array of multiplicative factors that specify how to increase time step.*/
   int flowsteps;  /* The number of time-steps taken by the flow code.*/
   int stepmax;    /* The maximum number of time-steps taken by the flow code.*/
   int nstpmax;    /* The maximum number of time-step increments.*/
   int kplot;      /* Printout steps.*/
   int write_init; /* Flag to printout initial conditions.*/
   int iprint; /* Print level (-1-none, 0-fields, >=1-vel, 2-perm/por, 3-pflow.bc)*/
   int imod;   /* screen printout  modulus*/
   int itecplot; /* tecplot print format (1-interchange x and z)*/
   int iblkfmt; /* blocked format*/
   int isync;  /* Synchronize pflow and ptran time steps (1)*/
   int ndtcmx; /* Steps needed after cutting to increase time step*/
   int newtcum;    /* Total number of Newton steps taken.*/
   int icutcum;    /* Total number of cuts in the timestep taken.*/
   int newton_max; /* Max number of Newton steps for one time step.*/
   int icut_max;   /* Max number of dt cuts for one time step.*/
   int iaccel,iphch;
   int iread_init; /* flag for reading initial conditions.*/
   /* Basically our target number of newton iterations per time step.*/
   double dpmxe,dsmxe; /*maximum allowed changes in field vars.*/
   double dpmax,dsmax;
   
   /* Grid topology*/
   int igeom;
   int nx, ny, nz  ;  /* Global domain dimensions of the grid.*/
   
   /* Arrays for indexing between local ghosted and non-ghosted, local to natural arrays.*/
   DA da_1_dof, da_3np_dof, da_ndof;
   /* DA's for 1, 3, and multiple (number of phases) degrees of freedom.*/
   /* da_ndof = total degrees of freedom per node*/
   
   int *i1bc;
   int *i2bc;
   int *j1bc;
   int *j2bc;
   int *k1bc;  
   int *k2bc;

   int *iregbc1;
   int *iregbc2;
   /* iregbc1(ibc) and iregbc2(ibc) give the id of the first region and */
   /* last region, respectively, that utilizes the boundary conditions in */
   /* boundary condition block ibc.*/
   
   
   /*block BC values read from input*/
   double *velocitybc0;
   double *xxbc0;
   double radius_0;
   /*   phik*/
   int iregperm, iran_por, iread_perm;
   double ran_fac;
   int *i1reg;
   int *i2reg;
   int *j1reg;
   int *j2reg;
   int *k1reg;
   int *k2reg;
   double *por_reg;
   double *tor_reg;
   double *perm_reg;

/*   initial conditions*/
    int iregini;
   int *i1ini;
   int *i2ini;
   int *j1ini;
   int *j2ini;
   int *k1ini;
   int *k2ini;
   double *xx_ini;
   
   
   /*   source term*/
   int nblksrc, ntimsrc, isrc1;
   int * i1src;
   int *i2src;
   int *j1src;
   int *j2src;
   int *k1src;
   int *k2src;
   double *timesrc;
   double *tempsrc;
   double *qsrc;

/*   solid reaction rate*/
   int ityprxn;
   double rk, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts;
   double qu_kin, yh2o_in_co2;

   /*   breakthrough curves*/
   int ibrkcrv;
   int *i1brk;
   int *i2brk;
   int *j1brk;
   int *j2brk;
   int *k1brk;
   int *k2brk;
   int *ibrktyp;
   int *ibrkface;

/*   dual continuum*/
   int idcdm, idcmblk;
   int *i1dcm;
   int *i2dcm;
   int *j1dcm;
   int *j2dcm;
   int *k1dcm;
   int *k2dcm;
   double *fracture_aperture;
   double *matrix_block;

   int * icap_reg;
   int *ithrm_reg;
   double scale; 
   double *rock_density;
   double *cpr;
   double *dencpr;
   double *ckdry;
   double *ckwet; 
   double *tau;
   double *cdiff;
   double *cexp;

   double *swir;
   double *lambda;
   double *alpha;
   double *pckrm;
   double *pcwmax;
   double *pcbetac;
   double *pwrprm;
   double *sir;
   int *icaptype;
   

   int ihydrostatic,ideriv;
   double dTdz,beta,tref,pref, gravity ;
   
   /*   table lookup*/
   int itable;

   double *dx0;
   double *dy0;
   double *dz0;
   double *rd;
   double *x;
   double *y;
   double *z;

   /*-------------------------------------------------------------------*/
   /* Quantities defined at each grid point.*/
   /* NOTE: I adopt the convention that _loc indicates the local portion*/
   /* of any global vector.*/
   /*-------------------------------------------------------------------*/
   
   /* One degree of freedom: Physical coordinates.*/
   Vec conc;
   Vec porosity, porosity0, porosity_loc, tor, tor_loc;
   Vec dx, dy, dz, dx_loc, dy_loc, dz_loc;  /* Grid spacings*/
   Vec volume;  /* Volume of a cell in the grid*/
   Vec ithrm, ithrm_loc, icap, icap_loc;
   Vec ttemp, ttemp_loc, temp; /* 1 dof*/

   /* Three degrees of freedom:*/
   /*   Vec perm, perm_loc*/
   Vec perm_xx, perm_xx_loc, perm_yy, perm_yy_loc, perm_zz, perm_zz_loc;
   Vec perm0_xx, perm0_yy, perm0_zz, perm_pow;
   /* Multiple degrees of freedom (equal to number of phases present):*/
   Vec r;            /* The residual.  (NOT the negative of the residual.)*/
   
   Vec vl;/* , vvl, vg, vvg ! phase (liquid and gas) velocities stored at interfaces*/
   
   /* Solution vectors*/
   Vec xx, xx_loc, dxx, yy, accum;
   /* Jacobian matrix*/
   Mat J;
   MatFDColoring matfdcoloring;
   /* Coloring used for computing the Jacobian via finite differences.*/
   
   
   double atol, rtol, stol, dtol;
   /* Absolute, relative, and "change in norm of solution" tolerances.*/
   int maxit, maxf;
   /* The maximum number of iterations and function evaluations, respectively*/
   
   /*  Vec p_nat, t_nat, c_nat, phis_nat, por_nat, vl_nat, s_nat, x_nat !, perm_nat*/
   /* Holds contents in natural ordering.*/
   /*   Vec p_all, t_all, c_all, phis_all, por_all, vl_all, s_all !, perm_all */
   /* Used to hold all values on processor 0.*/

   void *locpat;			

};
}
