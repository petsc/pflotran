
subroutine f_create_local_patch_data(patch_obj)

 use pflow_gridtype_module

 integer :: patch_obj
 type(pflow_localpatch_info), pointer :: grid_data
 allocate(grid_data)
 call pflow_loc(patch_obj, grid_data)

end subroutine f_create_local_patch_data

subroutine f_create_hierarchy_data(hierarchy_obj)

 use pflow_gridtype_module
 use pflow_grid_module
 
#include "include/finclude/petsc.h"


 external Pflow_allocate_Vec
 PetscFortranAddr :: hierarchy_obj
 type(pflowGrid), pointer :: grid_data

 allocate(grid_data)
 call pflow_loc(hierarchy_obj, grid_data)

end subroutine f_create_hierarchy_data

subroutine f_initialize_hierarchy_data(params)

 use pflow_gridtype_module
 use pflow_grid_module
 
#include "include/finclude/petsc.h"

 external pflowGrid_newpm
 type(pflowGridParameters), pointer :: params

 call pflowGrid_newpm(params)

end subroutine f_initialize_hierarchy_data
