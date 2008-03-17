subroutine f_create_simulation(simulation_obj)
 use Simulation_module
 use Discretization_module
 use AMR_Grid_module

 integer :: grid_obj
 type(simulation_type), pointer :: sim
 type(discretization_type), pointer :: discretization

 sim=>SimulationCreate()
 discretization => sim%realization%discretization
 discretization%amrgrid => AMRGridCreate()

 call assign_c_ptr(simulation_obj, sim)

end subroutine f_create_simulation

subroutine f_create_grid_data(grid_obj)

 use Grid_module

 integer :: grid_obj
 type(grid_type), pointer :: grid
 allocate(grid)
 call assign_c_ptr(grid_obj, grid)

end subroutine f_create_grid_data

subroutine f_create_hierarchy_data(amrgrid_obj)

 use AMR_Grid_module
 
#include "include/finclude/petsc.h"

 PetscFortranAddr :: amrgrid_obj
 type(amrgrid_type), pointer :: amrgrid_data

 allocate(amrgrid_data)
 call assign_c_ptr(amrgrid_obj, amrgrid_data)

end subroutine f_create_hierarchy_data
