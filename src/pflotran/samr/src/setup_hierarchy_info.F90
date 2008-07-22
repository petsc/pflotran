subroutine f_set_application_ptr(simulation, application_ptr)
  use Simulation_module
  use Discretization_module
  use AMR_Grid_module

  implicit none

#include "include/finclude/petsc.h"

 type(simulation_type), pointer :: simulation
 type(discretization_type), pointer :: discretization
 PetscFortranAddr :: application_ptr

 discretization => simulation%realization%discretization
 discretization%amrgrid%p_application = application_ptr

end subroutine f_set_application_ptr
