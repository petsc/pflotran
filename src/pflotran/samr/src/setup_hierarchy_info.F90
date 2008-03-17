subroutine f_set_hierarchy_ptr(simulation, hierarchy_ptr)
  use Simulation_module
  use Discretization_module
  use AMR_Grid_module

  implicit none

#include "include/finclude/petsc.h"

 type(simulation_type), pointer :: simulation
 type(discretization_type), pointer :: discretization
 PetscFortranAddr :: hierarchy_ptr

 discretization => simulation%realization%discretization
 discretization%amrgrid%p_samr_hierarchy = hierarchy_ptr

end subroutine f_set_hierarchy_ptr
