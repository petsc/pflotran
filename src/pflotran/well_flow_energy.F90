module Well_FlowEnergy_class

  use PFLOTRAN_Constants_module
  use Well_Flow_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_type) :: well_flow_energy_type
      PetscReal :: tw_ref                        ! [°C] well temperature at reference elevation
      PetscReal, pointer :: ent_ref(:)           ! ent_ref(iphase) MJ/kmol, well fluid enthalpy of iphase
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintFlowEnergy
    !procedure, public :: Init => WellAuxVarBaseInit
    !procedure, public :: Read => WellAuxVarBaseRead
    !procedure, public :: WellAuxVarClear => WellAuxVarBaseClear
    !procedure, public :: WellInit => WellBaseInit
    !procedure, public :: UpdateConnFactor
    !procedure, public :: Output
    !procedure  WellConnInit ! init all vars related to well connections
    !procedure  :: InitWellZRefCntrlConn
    !procedure  :: WellConnSort
  end type  well_flow_energy_type

  public :: WellFlowEnergyInit

contains

! ************************************************************************** !

subroutine PrintFlowEnergy(this)

  implicit none

  class(well_flow_energy_type) :: this

  write(*,*) "Well FlowEnergy Printing message"

end subroutine PrintFlowEnergy

! ************************************************************************** !

subroutine WellFlowEnergyInit(this,option)
  ! 
  ! Initializes variables/objects in flow and energy well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(option_type) :: option

  this%tw_ref = 0.0d0; 

  allocate( this%ent_ref(option%nphase) );
  this%ent_ref = 0.0d0;

end subroutine WellFlowEnergyInit

! ************************************************************************** !

end module Well_FlowEnergy_class


