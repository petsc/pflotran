module AuxVars_TOWG_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public, extends(auxvar_flow_energy_type) :: auxvar_towg_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    ! no data at the moment
  contains
    procedure, public :: Init => AuxVarTOWGInit
    procedure, public :: Strip => AuxVarTOWGStrip
  end type auxvar_towg_type

  public :: AuxVarTOWGStrip

contains

! ************************************************************************** !

! ************************************************************************** !

subroutine AuxVarTOWGInit(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 11/07/16
  ! 

  use Option_module

  implicit none
  
  class(auxvar_towg_type) :: this
  type(option_type) :: option

  this%effective_porosity = 0.d0
  this%pert = 0.d0
  
  call AuxVarFlowInit(this,option)

  call AuxVarFlowEnergyInit(this,option)

  this%istate_store = 0

end subroutine AuxVarTOWGInit

! ************************************************************************** !

subroutine AuxVarTOWGStrip(this)
  ! 
  ! TOilImsAuxVarDestroy: Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/30/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_towg_type) :: this

  call AuxVarFlowStrip(this)

  call AuxVarFlowEnergyStrip(this)

end subroutine AuxVarTOWGStrip
! ************************************************************************** !

end module AuxVars_TOWG_module

