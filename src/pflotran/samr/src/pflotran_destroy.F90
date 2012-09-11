subroutine f_simulation_destroy(simulation)
  use Simulation_module
  implicit none

#include "definitions.h"

  type(simulation_type), pointer :: simulation

  call SimulationDestroy(simulation)

end subroutine f_simulation_destroy
