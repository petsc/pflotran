subroutine f_initialize_simulation(simulation)

  use Simulation_module
  use Init_module

  implicit none

  type(simulation_type), pointer :: simulation

  call Init(simulation)

end subroutine f_initialize_simulation
