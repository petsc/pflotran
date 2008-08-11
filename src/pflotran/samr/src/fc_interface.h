

#ifndef fc_interface_h
#define fc_interface_h

extern "C"{

#include "petsc.h"

   void f_create_simulation_(void **p_data, void **p_application);
   void f_create_local_patch_data_(void **p_data);
   void f_create_hierarchy_data_(void **p_data);
   void f_create_integrator_(void **p_data);
   void f_initialize_hierarchy_data_(void **p_data);
   void f_setup_hierarchy_data_(void **p_data, void **p_integrator);

   void f_initialize_simulation_(void **p_simulation);

   void f_stepper_run_(void **p_simulation);

}


#endif
