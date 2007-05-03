#ifndef fc_interface_h
#define fc_interface_h
extern "C"{
   void f_create_local_patch_data_(void **p_data);
   void f_create_hierarchy_data_(void **p_data);
   void f_initialize_hierarchy_data_(void **p_data);

struct gridparameters{
   void *p_grid;
   void *p_timestep;
   int igeom;
   int nx, ny, nz;
   int npx, npy, npz;
   int nphase;
   bool usesamrai;
   void *p_samr_hierarchy;
};

struct time_integrator{
   double *tfac;
   double dt_min;
   double dt_max;
   char tunit[2];
   int iaccel;
   int icut_max;
   int nstpmax;
   int kplot; 
   double *tplot;
   double *tstep;
   double *dtstep;
   double dpmxe;
   double dsmxe;
};

}


#endif
