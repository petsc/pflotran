module Well_module

  use PFLOTRAN_Constants_module
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class
  use WellSpec_Base_class
  use Well_TOilIms_class
  !add here other well classes, e.g. Wells_XXXX_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: CreateWell

contains

! ************************************************************************** !

function CreateWell(well_spec,option)
  ! 
  ! Create a toil ims well object based on type specified in the well_spec
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module

  implicit none

  class(well_spec_base_type), pointer :: well_spec
  type(option_type) :: option
  class(well_base_type), pointer :: CreateWell
  
  select case(option%iflowmode)
    case(TOIL_IMS_MODE)
      CreateWell => CreateTOilImsWell(well_spec,option)
      !add here create well for other flow modes
    case default
      option%io_buffer = 'Well model supported supported for TOIL_IMS only'
      call printErrMsg(option)
  end select

  CreateWell%spec => well_spec

end function CreateWell


! ************************************************************************** !

end module Well_module


