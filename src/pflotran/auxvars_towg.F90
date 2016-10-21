module AuxVars_TOWG_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public, extends(auxvar_flow_energy_type) :: auxvar_towg_type
    ! no data at the moment
  contains
    !procedure, public :: Init => AuxVarTOilImsInit
    !procedure, public :: Strip => AuxVarTOilImsStrip
  end type auxvar_towg_type

  !public :: AuxVarTOilImsStrip

contains

! ************************************************************************** !


end module AuxVars_TOWG_module

