module Saturation_Function2_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  ! Saturation Function
  type :: sat_func_base_type
    PetscReal :: Sr
  contains
!    procedure, public :: Init => SF_Base_Init
    procedure, public :: CapillaryPressure => SF_Base_CapillaryPressure
    procedure, public :: Saturation => SF_Base_Saturation
  end type sat_func_base_type
  type, public, extends(sat_func_base_type) :: sat_func_VG_type
    PetscReal :: alpha
    PetscReal :: m
  contains
!    procedure, public :: Init => SF_VG_Init
    procedure, public :: CapillaryPressure => SF_VG_CapillaryPressure
    procedure, public :: Saturation => SF_VG_Saturation
  end type sat_func_VG_type  
  
  ! Relative Permeability Function
  type :: rel_perm_func_base_type
    PetscReal :: Sr
  contains
!    procedure, public :: Init => RPF_Base_Init
    procedure, public :: RelativePermeability => RPF_Base_RelPerm
  end type rel_perm_func_base_type
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_Mualem_type
    PetscReal :: m
  contains
!    procedure, public :: Init => RPF_Base_Init
    procedure, public :: RelativePermeability => RPF_Mualem_RelPerm
  end type rel_perm_func_Mualem_type

  type, public :: characteristic_curve_type
    class(sat_func_base_type), pointer :: saturation_function
    class(rel_perm_func_base_type), pointer :: liq_rel_perm_func
    class(rel_perm_func_base_type), pointer :: gas_rel_perm_func
  end type characteristic_curve_type
  
contains

! ************************************************************************** !

function CharacteristicCurveCreate()
  ! 
  ! Creates a characteristic curve object that holds parameters and pointers
  ! to functions for calculating saturation, capillary pressure, relative
  ! permeability, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/23/14
  ! 

  implicit none

  class(characteristic_curve_type), pointer :: CharacteristicCurveCreate
  
  class(characteristic_curve_type), pointer :: characteristic_curve
  
  allocate(characteristic_curve)
  nullify(characteristic_curve%saturation_function)
  nullify(characteristic_curve%liq_rel_perm_func)
  nullify(characteristic_curve%gas_rel_perm_func)

  CharacteristicCurveCreate => characteristic_curve

end function CharacteristicCurveCreate

!!!BASE ROUTINES

! ************************************************************************** !

subroutine SFBaseCreate(this)

  implicit none
  
  class(sat_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  this%Sr = 0.d0
  
end subroutine SFBaseCreate

! ************************************************************************** !

subroutine RPFBaseCreate(this)

  implicit none
  
  class(rel_perm_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  this%Sr = 0.d0
  
end subroutine RPFBaseCreate

! ************************************************************************** !

subroutine SF_Base_CapillaryPressure(this,saturation,capillary_pressure, &
                                     option)
  
  use Option_module
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SF_Base_CapillaryPressure must be extended.'
  call printErrMsg(option)
  
end subroutine SF_Base_CapillaryPressure

! ************************************************************************** !

subroutine SF_Base_Saturation(this,capillary_pressure,saturation,dsat_pres, &
                              option)

  use Option_module
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: saturation
  PetscReal, intent(out) :: dsat_pres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SF_Base_Saturation must be extended.'
  call printErrMsg(option)
  
end subroutine SF_Base_Saturation

! ************************************************************************** !

subroutine RPF_Base_RelPerm(this,saturation,relative_permeability,dkr_Se, &
                            option)
  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  PetscReal, intent(in) :: saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'RPF_Base_RelPerm must be extended.'
  call printErrMsg(option)
  
end subroutine RPF_Base_RelPerm

! ************************************************************************** !

function SF_VG_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_VG_type), pointer :: SF_VG_Create
  
  class(sat_func_VG_type), pointer :: sf

  allocate(sf)
  call SFBaseCreate(sf)
  sf%alpha = 0.d0
  sf%m = 0.d0
  
  SF_VG_Create => sf
  
end function SF_VG_Create

! ************************************************************************** !

subroutine SF_VG_CapillaryPressure(this,saturation,capillary_pressure, &
                                   option)
  ! 
  ! Computes the capillary_pressure as a function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  
  implicit none
  
  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: n
  PetscReal :: Se
  PetscReal :: one_plus_pc_alpha_n
  PetscReal :: pc_alpha_n
  PetscReal :: pc_alpha
  
  n = 1.d0/(1.d0-this%m)
  Se = (saturation-this%Sr)/(1.d0-this%Sr)
  one_plus_pc_alpha_n = Se**(-1.d0/this%m)
  pc_alpha_n = one_plus_pc_alpha_n - 1.d0
  pc_alpha = pc_alpha_n**(1.d0/n)
  capillary_pressure = pc_alpha/this%alpha
  
end subroutine SF_VG_CapillaryPressure

! ************************************************************************** !

subroutine SF_VG_Saturation(this,capillary_pressure,saturation,dsat_pres, &
                            option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  
  implicit none

  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: saturation
  PetscReal, intent(out) :: dsat_pres
  type(option_type), intent(inout) :: option
  
  PetscReal, parameter :: pc_alpha_n_epsilon = 1.d-15
  PetscReal :: n
  PetscReal :: pc_alpha
  PetscReal :: pc_alpha_n
  PetscReal :: one_plus_pc_alpha_n
  PetscReal :: Se
  PetscReal :: dSe_pc
  PetscReal :: dsat_pc
  
  dsat_pres = 0.d0
  
#if 0
  if (capillary_pressure < this%spline%low) then
    saturation = 1.d0
    relative_perm = 1.d0
    switch_to_saturated = PETSC_TRUE
    return
  else if (capillary_pressure < this%spline%high) then
    call CubicPolynomialEvaluate(this%spline%coefficients, &
                                 capillary_pressure,Se,dSe_pc)
    saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_pc = (1.d0-this%Sr)*dSe_pc
#else
  if (capillary_pressure <= 0.d0) then
    saturation = 1.d0
!    switch_to_saturated = PETSC_TRUE
    return
#endif        
  else
    n = 1.d0/(1.d0-this%m)
    pc_alpha = capillary_pressure*this%alpha
    pc_alpha_n = pc_alpha**n
    !geh:  This conditional does not catch potential cancelation in 
    !      the dkr_Se deriviative calculation.  Therefore, I am setting
    !      an epsilon here
!        if (1.d0 + pc_alpha_n == 1.d0) then ! check for zero perturbation
    if (pc_alpha_n < pc_alpha_n_epsilon) then 
      saturation = 1.d0
!      switch_to_saturated = PETSC_TRUE
      return
    endif
    one_plus_pc_alpha_n = 1.d0+pc_alpha_n
    Se = one_plus_pc_alpha_n**(-this%m)
    dSe_pc = -this%m*n*this%alpha*pc_alpha_n/ &
            (pc_alpha*one_plus_pc_alpha_n**(this%m+1.d0))
    saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_pc = (1.d0-this%Sr)*dSe_pc
  endif
  
end subroutine SF_VG_Saturation

! ************************************************************************** !

subroutine RPF_Mualem_RelPerm(this,saturation,relative_permeability,dkr_Se, &
                              option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rel_perm_func_Mualem_type) :: this
  PetscReal, intent(in) :: saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: one_over_m
  PetscReal :: Se_one_over_m
  
  Se = (saturation - this%Sr) / (1.d0 - this%Sr)
#ifdef MUALEM_SPLINE
  if (Se > this%spline%low) then
    call CubicPolynomialEvaluate(this%spline%coefficients, &
                                 Se,relative_permeability,dkr_Se)
  else
#endif          
  one_over_m = 1.d0/this%m
  Se_one_over_m = Se**one_over_m
  relative_permeability = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**this%m)**2.d0
  dkr_Se = 0.5d0*relative_permeability/Se+ &
            2.d0*Se**(one_over_m-0.5d0)* &
                (1.d0-Se_one_over_m)**(this%m-1.d0)* &
                (1.d0-(1.d0-Se_one_over_m)**this%m)
#ifdef MUALEM_SPLINE
  endif
#endif          
  
end subroutine RPF_Mualem_RelPerm

end module Saturation_Function2_module