module AuxVars_FlowEnergy_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public, extends(auxvar_flow_type) :: auxvar_flow_energy_type
    PetscReal :: temp
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
  contains
   !..............
  end type auxvar_flow_energy_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!  interface TOilImsAuxVarDestroy
!    module procedure TOilImsAuxVarSingleDestroy
!    module procedure TOilImsAuxVarArray1Destroy
!    module procedure TOilImsAuxVarArray2Destroy
!  end interface TOilImsAuxVarDestroy
  
!  public :: TOilImsAuxCreate, &
!            TOilImsAuxDestroy, &
!            TOilImsAuxVarInit, &
!            TOilImsAuxVarCompute, &
!            TOilImsAuxVarPerturb, &
!            TOilImsAuxVarDestroy, &
!            TOilImsAuxVarStrip

!contains


! ************************************************************************** !

! ************************************************************************** !

end module AuxVars_FlowEnergy_module

