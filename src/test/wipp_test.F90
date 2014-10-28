
program wipp_test
  

  use PFLOTRAN_Constants_module
  use Creep_Closure_module
  use Input_Aux_module
  use Option_module

  implicit none

#include "finclude/petscsys.h"

  class(creep_closure_type), pointer :: creep_closure
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  
  character(len=MAXSTRINGLENGTH) :: string
  
  option => OptionCreate()
  call OptionInitMPI(option)
  creep_closure => CreepClosureCreate()
  string = 'input_block.txt'
  input => InputCreate(IUNIT_TEMP,string,option)
  call creep_closure%Read(input,option)
  call InputDestroy(input)
  call CreepClosureDestroy(creep_closure)
  
end program wipp_test
