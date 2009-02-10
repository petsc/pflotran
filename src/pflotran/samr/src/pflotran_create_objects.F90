subroutine f_create_simulation(simulation_obj, application_ptr)
 
 use Simulation_module
 use Realization_module
 use Discretization_module
 use Option_module
 use Input_module
 use Init_module
 use Logging_module
 use Stochastic_module
 use Stochastic_Aux_module  
 use AMR_Grid_module

 implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

 PetscFortranAddr :: simulation_obj
 PetscFortranAddr :: application_ptr
 PetscLogDouble :: timex(4), timex_wall(4)
 
 PetscTruth :: truth
 PetscTruth :: option_found  
 PetscInt :: test_int
 PetscErrorCode :: ierr
 character(len=MAXSTRINGLENGTH) :: string
 PetscMPIInt :: myrank, commsize
 
 type(stochastic_type), pointer :: stochastic
 type(simulation_type), pointer :: simulation
 type(realization_type), pointer :: realization
 type(option_type), pointer :: option
 type(discretization_type), pointer :: discretization
 
 option => OptionCreate()
 option%fid_out = IUNIT2

 call MPI_Comm_rank(PETSC_COMM_WORLD,myrank, ierr)
 call MPI_Comm_size(PETSC_COMM_WORLD,commsize,ierr)
 
 option%mycomm    = PETSC_COMM_WORLD
 option%myrank = myrank
 option%mycommsize = commsize

 option%input_filename = "pflotran_well.in"

 string = '-pflotranin'
 call InputGetCommandLineString(string,option%input_filename,option_found,option)
 
 string = '-screen_output'
 call InputGetCommandLineTruth(string,option%print_to_screen,option_found,option)
 
 string = '-file_output'
 call InputGetCommandLineTruth(string,option%print_to_file,option_found,option)

 if (associated(stochastic)) then
    call StochasticInit(stochastic,option)
    call StochasticRun(stochastic,option)
 else
    call LoggingCreate()
    
    simulation => SimulationCreate(option)
    realization => simulation%realization
    
    call OptionCheckCommandLine(option)
    
    discretization => simulation%realization%discretization
    
    if(.not.(application_ptr.eq.0)) then
       discretization%amrgrid => AMRGridCreate()
       nullify(discretization%grid)
       discretization%amrgrid%p_application = application_ptr
    else
       nullify(discretization%amrgrid)
    endif
    
    call assign_c_ptr(simulation_obj, simulation)
 endif

end subroutine f_create_simulation

