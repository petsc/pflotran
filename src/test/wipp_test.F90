
program wipp_test
  

  use PFLOTRAN_Constants_module
  use Creep_Closure_module
  use Input_Aux_module
  use Option_module

  implicit none

#include "finclude/petscsys.h"

  class(creep_closure_type), pointer :: cc
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  
  PetscReal :: time
  PetscReal :: pressure
  character(len=MAXSTRINGLENGTH) :: string
  
  option => OptionCreate()
  call OptionInitMPI(option)
  call CreepClosureInit()
  creep_closure => CreepClosureCreate()
  string = 'input_block.txt'
  input => InputCreate(IUNIT_TEMP,string,option)
  call creep_closure%Read(input,option)
  call InputDestroy(input)

  time = 500.d0*365.d0*24.d0*3600.d0
  pressure = 12.d6
  call creep_closure%Test(time,pressure)
  time = 1000.d0*365.d0*24.d0*3600.d0
  pressure = 12.d7
  call creep_closure%Test(time,pressure)
  time = 1500.d0*365.d0*24.d0*3600.d0
  pressure = 1.d5
  call creep_closure%Test(time,pressure)
  time = 2000.d0*365.d0*24.d0*3600.d0
  pressure = 1.d8
  call creep_closure%Test(time,pressure)
  time = 2000.d0*365.d0*24.d0*3600.d0
  pressure = 18.5d6
  call creep_closure%Test(time,pressure)
  
  call CreepClosureDestroy(creep_closure)
  
end program wipp_test
