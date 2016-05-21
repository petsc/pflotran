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

  use Connection_module
  use Option_module

  implicit none

  class(well_spec_base_type), pointer :: well_spec
  character(len=MAXWORDLENGTH) :: coupler_name
  type(option_type) :: option
  class(well_base_type), pointer :: CreateWell

  character(len=MAXWORDLENGTH) :: wfile_name
  
  select case(option%iflowmode)
    case(TOIL_IMS_MODE)
      CreateWell => CreateTOilImsWell(well_spec,option)
      !add here create well for other flow modes
    case default
      option%io_buffer = 'Well model supported supported for TOIL_IMS only'
      call printErrMsg(option)
  end select

  !Debug printing 
  write(*,*) "well_factor type = ", CreateWell%spec%well_fact_itype 
  write(*,*) "well type = ", CreateWell%spec%ctype
  write(*,*) "radius = ", CreateWell%spec%radius

  select type(CreateWell)
    class is(well_toil_ims_type)
      write(*,*) "temp", CreateWell%tw_ref
      write(*,*) "well_press", CreateWell%pw_ref
  end select 

  call CreateWell%PrintMsg(); 

  !Create well outfile and write its header 

end function CreateWell

! ************************************************************************** !

end module Well_module


