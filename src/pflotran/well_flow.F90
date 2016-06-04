module Well_Flow_class

  use PFLOTRAN_Constants_module
  use Well_Base_class
  use Condition_module
  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_base_type) :: well_flow_type
    PetscReal :: pw_ref                        ! [Pa] well pressure at reference elevation
    PetscReal, pointer :: dw_ref(:)            ! dw_ref(iphase) [kg/m3] well fluid density of iphase at reference elevation
    PetscReal, pointer :: q_fld(:)             ! q_fld(iphase)  [m3/s] well fluid flow rates of iphase
    PetscReal, pointer :: mr_fld(:)            ! mr_fld(iphase) [kg/s] well fluid mass rates of iphase
    PetscReal, pointer :: conn_h(:)            ! connection hydrostatic pressure corrections
    type(flow_condition_type), pointer :: flow_condition ! pointer to flow_condition associated with the well
    !PetscReal, pointer :: conn_mobs(:,:)       ! well connection mobilities ! TO REMOVE - computed when needed flight
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintFlow
    procedure, public :: ConnInit => WellFlowConnInit
    procedure, public :: ExplUpdate => FlowExplUpdate
    !procedure, public :: Init => WellAuxVarBaseInit
    !procedure, public :: Read => WellAuxVarBaseRead
    !procedure, public :: WellAuxVarClear => WellAuxVarBaseClear
    !procedure, public :: WellInit => WellBaseInit
    !procedure, public :: UpdateConnFactor
    !procedure, public :: Output
    !procedure  WellConnInit ! init all vars related to well connections
    !procedure  :: InitWellZRefCntrlConn
    !procedure  :: WellConnSort
  end type  well_flow_type

  public :: WellFlowInit, WellFlowConnInit

contains

! ************************************************************************** !

subroutine PrintFlow(this)

  implicit none

  class(well_flow_type) :: this

  write(*,*) "Well Flow Printing message"

end subroutine PrintFlow

! ************************************************************************** !

subroutine WellFlowInit(this,option)
  ! 
  ! Initializes variables/objects in flow well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module

  implicit none

  class(well_flow_type) :: this

  type(option_type) :: option

  this%pw_ref = 0.0d0;

  allocate( this%dw_ref(option%nphase) );
  this%dw_ref = 0.0d0;
  allocate( this%q_fld(option%nphase) );
  this%q_fld = 0.0d0;
  allocate( this%mr_fld(option%nphase) );
  this%mr_fld = 0.0d0;
  
end subroutine WellFlowInit

! ************************************************************************** !

subroutine WellFlowConnInit(this,num_connections,option)
  ! 
  ! Allocate and initilize well_base connections arrays
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/20/2016
  !

  use Option_module

  implicit none

  class(well_flow_type) :: this
  PetscInt, intent(in) :: num_connections 
  type(option_type) :: option  

  nullify(this%conn_h);
  allocate(this%conn_h(num_connections));
  this%conn_h = 0.0d0; 

end subroutine WellFlowConnInit

! ************************************************************************** !

subroutine FlowExplUpdate(this,grid,option)
  ! 
  ! - Update FlowEnergy well vars
  ! - Perform a limit on well checks 
  ! - Update well control variable in case of switch when a limit is reached
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 6/03/2016
  ! 

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  print *, "Well => FlowExplUpdate must be extended"
  stop  

end subroutine FlowExplUpdate

! ************************************************************************** !

end module Well_Flow_class
