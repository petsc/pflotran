program pflotran_interface_main
#include "finclude/petscsys.h"

use pflotran_model_module

implicit none

  
  type(pflotran_model_type),pointer :: pflotran_m

  allocate(pflotran_m)

  pflotran_m => pflotranModelCreate()	


  call pflotranModelStepperRunInit(pflotran_m)
  call pflotranModelStepperRunTillPauseTime(pflotran_m,1.0d0)  
  
  call pflotranModelUpdateTopBCHomogeneous(pflotran_m, 2.0d0 * 3.171d-10 )
  
  call pflotranModelStepperRunTillPauseTime(pflotran_m,2.0d0)  
  call pflotranModelStepperRunTillPauseTime(pflotran_m,3.0d0)  

  call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)
  

end program pflotran_interface_main
