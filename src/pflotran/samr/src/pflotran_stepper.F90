subroutine f_stepper_run(simulation)
  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Init_module
  use Logging_module

  implicit none

#include "definitions.h"
#include "include/finclude/petsclog.h"

  type(simulation_type), pointer :: simulation

  call StepperRun(simulation%realization,simulation%flow_stepper, &
                  simulation%tran_stepper)

end subroutine f_stepper_run
