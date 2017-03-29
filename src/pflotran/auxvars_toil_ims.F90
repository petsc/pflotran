module AuxVars_TOilIms_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module

  implicit none
  
  private 

  type, public, extends(auxvar_flow_energy_type) :: auxvar_toil_ims_type
    ! no data at the moment
  contains
    procedure, public :: Init => AuxVarTOilImsInit
    procedure, public :: Strip => AuxVarTOilImsStrip
  end type auxvar_toil_ims_type

  public :: AuxVarTOilImsStrip

contains

! ************************************************************************** !

subroutine AuxVarTOilImsInit(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  ! 

  use Option_module

  implicit none
  
  class(auxvar_toil_ims_type) :: this
  type(option_type) :: option

  this%effective_porosity = 0.d0
  this%pert = 0.d0

  !PO to do: 
  !allign TOIL_IMS AuxVarCompute to new pc array (see AuxVarFlowInit), 
  !then remove the following iniitlisation block, and replace with a call to 
  !AuxVarFlowInit  
  !the indended part below is the block to replace with AuxVarFlowInit
    !two phases (water,oil) and capillary pressure
    allocate(this%pres(option%nphase+ONE_INTEGER))
    this%pres = 0.d0
    allocate(this%sat(option%nphase))
    this%sat = 0.d0
    allocate(this%den(option%nphase))
    this%den = 0.d0
    allocate(this%den_kg(option%nphase))
    this%den_kg = 0.d0
    allocate(this%mobility(option%nphase))
    this%mobility = 0.d0
    !if visc needed only for toil_ims, use AuxVarFlowInit to allocate sat, 
    !den, den_kg, mobility then allocate viscosity in here
    !more likley to allocate viscosity for all flow module though
    allocate(this%viscosity(option%nphase))
    this%viscosity = 0.d0
  ! end block to replace with AuxVarFlowInit

  !PO to do:
  !replace the indended block below with AuxVarFlowEnergyInit
    this%temp = 0.d0
    allocate(this%H(option%nphase))
    this%H = 0.d0
    allocate(this%U(option%nphase))
    this%U = 0.d0
  !end block to replace with AuxVarFlowEnergyInit 

end subroutine AuxVarTOilImsInit

! ************************************************************************** !

subroutine AuxVarTOilImsStrip(this)
  ! 
  ! TOilImsAuxVarDestroy: Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/30/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_toil_ims_type) :: this

  call AuxVarFlowStrip(this)

  call AuxVarFlowEnergyStrip(this)

end subroutine AuxVarTOilImsStrip
! ************************************************************************** !

end module AuxVars_TOilIms_module

