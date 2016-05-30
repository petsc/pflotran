module AuxVars_Base_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  !BEGINNING-Paramters to move into toil_ims_paramters 

  ! Primary DOF indices 
  ! Indices used to map aux_real for condition values 
  ! these variables, which are global to general, can be modified

  !phase mapping:

  type, public :: auxvar_base_type
    PetscReal :: effective_porosity ! factors in compressibility - common to all modes??
    PetscReal :: pert ! common to all modes to (perturbation for numerical jacobian)
  contains
    procedure, public :: Init => InitAuxVarBase
  end type auxvar_base_type


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
            

contains

! ************************************************************************** !

subroutine InitAuxVarBase(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  ! 

  use Option_module

  implicit none
  
  class(auxvar_base_type) :: this
  type(option_type) :: option

  !currently does nothing - could init the base members
  print *, 'Must extend InitBaseAuxVars '
  stop    

end subroutine InitAuxVarBase


! ************************************************************************** !

end module AuxVars_Base_module

