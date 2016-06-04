module Well_FlowEnergy_class

  use PFLOTRAN_Constants_module
  use Well_Flow_class
  use AuxVars_FlowEnergy_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_type) :: well_flow_energy_type
      PetscReal :: tw_ref                        ! [°C] well temperature at reference elevation
      PetscReal, pointer :: ent_ref(:)           ! ent_ref(iphase) MJ/kmol, well fluid enthalpy of iphase
      class(auxvar_flow_energy_type), pointer :: auxvar_flow_energy(:,:)
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintFlowEnergy
    procedure, public :: ExplUpdate => FlowEnergyExplUpdate
    procedure, public :: VarsExplUpdate => FlowEnergyVarsExplUpdate
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
  ! Date: 5/20/2016
  ! 

  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(option_type) :: option

  this%tw_ref = 0.0d0; 

  allocate( this%ent_ref(option%nphase) );
  this%ent_ref = 0.0d0;

  nullify(this%auxvar_flow_energy);

end subroutine WellFlowEnergyInit

! ************************************************************************** !

subroutine FlowEnergyExplUpdate(this,grid,option)
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

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  write(*,"('FlowEnergyExpl d11 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%den(1)
  write(*,"('FlowEnergyExpl d12 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%den(2) 
  write(*,"('FlowEnergyExpl p11 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%pres(1) 
  write(*,"('FlowEnergyExpl t1 before = ',e10.4)"), this%auxvar_flow_energy(0,1)%temp 


  call this%VarsExplUpdate(grid,option)

end subroutine FlowEnergyExplUpdate

! ************************************************************************** !

subroutine FlowEnergyVarsExplUpdate(this,grid,option)

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  print *, "FlowEnergyVarsExplUpdate must be extended"
  stop

end subroutine FlowEnergyVarsExplUpdate

! ************************************************************************** !

end module Well_FlowEnergy_class


