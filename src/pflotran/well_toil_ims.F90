module Well_TOilIms_class

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_energy_type) :: well_toil_ims_type
    !class(auxvar_toil_ims_type), pointer :: toil_ims_auxvars(:,:)
    ! .................
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintTOilIms
    procedure, public :: ConnInit => WellTOilImsConnInit
    procedure, public  :: PrintOutputHeader => PrintOutputHeaderWellTOilIms
    procedure, public :: VarsExplUpdate => TOilImsVarsExplUpdate
  end type  well_toil_ims_type

  type, public, extends(well_toil_ims_type) :: well_toil_ims_wat_inj_type
    ! ......................
    contains
      procedure, public :: PrintMsg => PrintTOilImsWatInj
      procedure, public :: VarsExplUpdate => TOilImsWatInjVarsExplUpdate
  end type

  type, public, extends(well_toil_ims_type) :: well_toil_ims_oil_inj_type
    ! ......................
  end type

  type, public, extends(well_toil_ims_type) :: well_toil_ims_oil_prod_type
    ! ......................
  end type

  type, public, extends(well_toil_ims_type) :: well_toil_ims_wat_prod_type
    ! ......................
  end type

  public :: CreateTOilImsWell

contains

! ************************************************************************** !

subroutine PrintTOilIms(this)

  implicit none

  class(well_toil_ims_type) :: this

  write(*,*) "Well TOilIms Printing message"

end subroutine PrintTOilIms

! ************************************************************************** !

subroutine PrintOutputHeaderWellTOilIms(this,output_option,file_unit)
  ! 
  ! Write header for well_TOilIms output file
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Output_Aux_module

  implicit none

  class(well_toil_ims_type) :: this
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: file_unit

  character(len=MAXWORDLENGTH) :: tunit

  tunit = trim(output_option%tunit)

  !write(IUNIT_TEMP,*) " VARIABLES = " // &
  !    '"Time [' // trim(output_option%tunit) // ']", ' // &
  !              """Pw[Pa]"", " // &
  !      """Tw[C]"", ""dh2o[kg/m3]"", ""doil[kg/[t]]"", "// &
  !      """Qwat[m3/[t]]"", ""Qoil[m3/[t]]"", " // &
  !      """Mwat[kg/[t]]"", ""Moil[kg/[t]]""" 

  !TODO: can do something more clever than this: 
  !      e.g. small loop to add well vars
  write(IUNIT_TEMP,*) " VARIABLES = " // &
      '"Time [' // trim(tunit) // ']", ' // &
                '""Pw[Pa]"", ' // &
        '"Tw[C]", "dh2o[kg/m3]", ' // &
        '"doil[kg/' // trim(tunit) // ']", ' // &
        '"Qwat[m3/' // trim(tunit) // ']", ' // &
        '"Qoil[m3/' // trim(tunit) //  ']", ' // &
        '"Mwat[kg/' // trim(tunit) // ']", ' // &
        '"Moil[kg/' // trim(tunit) // ']"' 


end subroutine PrintOutputHeaderWellTOilIms


! ************************************************************************** !

subroutine PrintTOilImsWatInj(this)

  implicit none

  class(well_toil_ims_wat_inj_type) :: this

  write(*,*) "Well TOilImsWatInj Printing message"

end subroutine PrintTOilImsWatInj

! ************************************************************************** !

!subroutine CreateTOilImsWell(well_spec,option)
function CreateTOilImsWell(well_spec,option)
  ! 
  ! Create a toil ims well object based on type specified in the well_spec
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module
  use WellSpec_Base_class

  implicit none

  class(well_spec_base_type), pointer :: well_spec
  type(option_type) :: option

  class(well_toil_ims_type), pointer :: CreateTOilImsWell

  class(well_toil_ims_wat_inj_type), pointer :: well_toil_ims_wat_inj
  class(well_toil_ims_oil_prod_type), pointer :: well_toil_ims_oil_prod
  !class(well_toil_ims_type), pointer :: well_toil_ims


  select case(well_spec%itype)
    case( WATER_INJ_WELL_TYPE )
      allocate(well_toil_ims_wat_inj);
      !well_toil_ims => well_toil_ims_wat_inj;
      CreateTOilImsWell => well_toil_ims_wat_inj;
    case( OIL_PROD_WELL_TYPE )
      allocate(well_toil_ims_oil_prod);
      !well_toil_ims => well_toil_ims_oil_prod;
      CreateTOilImsWell => well_toil_ims_oil_prod

    ! need to add water producer and oil injector
    case default
      option%io_buffer = 'Well type not recognize in CreateTOilImsWell'
      call printErrMsg(option)
  end select

  !initialise different well components 
  call WellBaseInit(CreateTOilImsWell,well_spec,option);
  call WellFlowInit(CreateTOilImsWell,option);
  call WellFlowEnergyInit(CreateTOilImsWell,option);

!end subroutine CreateTOilImsWell
end function CreateTOilImsWell

! ************************************************************************** !

subroutine WellTOilImsConnInit(this,num_connections,option)
  ! 
  ! Allocate and initilize well_base connections arrays
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/20/2016

  use Option_module

  implicit none

  class(well_toil_ims_type) :: this
  PetscInt, intent(in) :: num_connections 
  type(option_type) :: option  

  call WellBaseConnInit(this,num_connections,option);
  call WellFlowConnInit(this,num_connections,option);

end subroutine WellTOilImsConnInit

! ************************************************************************** !

subroutine TOilImsVarsExplUpdate(this,grid,option)

  use Grid_module
  use Option_module

  class(well_toil_ims_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  print *, "TOilImsVarsExplUpdate must be extended"
  stop

end subroutine TOilImsVarsExplUpdate

! ************************************************************************** !

subroutine TOilImsWatInjVarsExplUpdate(this,grid,option)

  use Grid_module
  use Option_module

  class(well_toil_ims_wat_inj_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  write(*,"('TOilImsWatInj d11 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%den(1)
  write(*,"('TOilImsWatInj d12 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%den(2) 
  write(*,"('TOilImsWatInj p11 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%pres(1) 
  write(*,"('TOilImsWatInj t1 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%temp 

  write(*,"('TOilImsWatInj rate = ',e10.4)"), &
          !this%flow_condition%toil_ims%rate%dataset%rarray(1)
          this%flow_condition%flow_well%rate%dataset%rarray(1)
  write(*,"('TOilImsWatInj temp = ',e10.4)"), &
          !this%flow_condition%toil_ims%temperature%dataset%rarray(1)
          this%flow_condition%flow_well%temperature%dataset%rarray(1) 

  write(*,"('TOilImsWatInj press = ',e10.4)"), &
          !this%flow_condition%toil_ims%pressure%dataset%rarray(1)
          this%flow_condition%flow_well%pressure%dataset%rarray(1) 

end subroutine TOilImsWatInjVarsExplUpdate

! ************************************************************************** !


end module Well_TOilIms_class




