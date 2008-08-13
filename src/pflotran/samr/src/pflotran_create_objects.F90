subroutine f_create_simulation(simulation_obj, application_ptr)
 
 use Simulation_module
 use Discretization_module
 use AMR_Grid_module

 implicit none
#include "include/finclude/petsc.h"

 PetscFortranAddr :: simulation_obj
 PetscFortranAddr :: application_ptr
 type(simulation_type), pointer :: sim
 type(discretization_type), pointer :: discretization

 sim=>SimulationCreate()
 discretization => sim%realization%discretization

 if(.not.(application_ptr.eq.0)) then
    discretization%amrgrid => AMRGridCreate()
    nullify(discretization%grid)
    discretization%amrgrid%p_application = application_ptr
 else
    nullify(discretization%amrgrid)
 endif

 call assign_c_ptr(simulation_obj, sim)

end subroutine f_create_simulation

