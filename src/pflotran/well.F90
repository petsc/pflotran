module Well_module

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class
  use Well_WaterInjector_class
  use Well_TOilIms_class
  !add here other well classes, e.g. Wells_XXXX_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: CreateWell, WellAuxVarSetUp

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

  !Debug printing 
  write(*,*) "well_factor type = ", CreateWell%spec%well_fact_itype 
  write(*,*) "well type = ", CreateWell%spec%ctype
  write(*,*) "radius = ", CreateWell%spec%radius

  select type(CreateWell)
    class is(well_toil_ims_wat_inj_type)
      write(*,*) "temp", CreateWell%tw_ref
      write(*,*) "well_press", CreateWell%pw_ref
  end select 

  call CreateWell%PrintMsg(); 

  !Create well outfile and write its header 
  !not here - otherwise will attempt to create a file for each process
  !at this stage a well is created in each process - even if empty   


end function CreateWell

! ************************************************************************** !
subroutine WellAuxVarSetUp(well,connection_set,flow_condition,aux,option)
  ! 
  ! Create a toil ims well object based on type specified in the well_spec
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/03/16
  ! 

  use Option_module
  use Auxiliary_module
  use Condition_module
  use Connection_module

  implicit none

  class(well_base_type), pointer :: well
  type(connection_set_type), pointer :: connection_set 
  type(flow_condition_type), pointer :: flow_condition
  type(auxiliary_type) :: aux  
  type(option_type) :: option

  write(*,"('WS d11 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%den(1)
  write(*,"('WS d12 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%den(2) 
  write(*,"('WS p11 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%pres(1) 
  write(*,"('WS t1 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%temp 

  select type(well)
    !if only auxvar_flow_energy needed can use class is(well_flow_energy_type)
    class is(well_toil_ims_wat_inj_type)
      well%flow_energy_auxvars => aux%TOil_ims%auxvars   
      well%flow_auxvars => aux%TOil_ims%auxvars
    !when well implmented for other flow modes - add below
  end select 

  select type(well)
    class is(well_flow_type)
      well%flow_condition => flow_condition

      nullify(well%well_conn_den_kg)
      allocate(well%well_conn_den_kg(well%well_num_conns))
      well%well_conn_den_kg = 0.0d0
      nullify(well%well_conn_h_sorted)
      allocate(well%well_conn_h_sorted(well%well_num_conns))
      well%well_conn_h_sorted = 0.0d0

  end select

  !for well base
  well%connection_set => connection_set


end subroutine WellAuxVarSetUp

! ************************************************************************** !

end module Well_module


