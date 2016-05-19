module Well_Flow_class

  use PFLOTRAN_Constants_module
  use Well_Base_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_base_type) :: well_flow_type
    PetscReal :: pw_ref                        ! [Pa] well pressure at reference elevation
    PetscReal, pointer :: dw_ref(:)            ! dw_ref(iphase) [kg/m3] well fluid density of iphase at reference elevation
    PetscReal, pointer :: q_fld(:)             ! q_fld(iphase)  [m3/s] well fluid flow rates of iphase
    PetscReal, pointer :: mr_fld(:)            ! mr_fld(iphase) [kg/s] well fluid mass rates of iphase
    PetscReal, pointer :: conn_h(:)            ! connection hydrostatic pressure corrections
    !PetscReal, pointer :: conn_mobs(:,:)       ! well connection mobilities ! TO REMOVE - computed when needed flight
  contains  ! add here type-bound procedure 
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

  public :: WellFlowInit

contains

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

  allocate( this%dw_ref(option%nphase) );
  this%dw_ref = 0.0d0;
  allocate( this%q_fld(option%nphase) );
  this%q_fld = 0.0d0;
  allocate( this%mr_fld(option%nphase) );
  this%mr_fld = 0.0d0
  nullify(this%conn_h);
  !nullify(this%conn_mobs);

end subroutine WellFlowInit

! ************************************************************************** !

end module Well_Flow_class
