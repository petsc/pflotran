module EOS_Oil_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "finclude/petscsys.h"

  ! public oil eos variables
  PetscReal, public :: fmw_oil

  ! module variables
  PetscReal :: constant_density
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity
  PetscReal :: constant_sp_heat

  ! In order to support generic EOS subroutines, we need the following:
  ! 1. An interface declaration that defines the argument list (best to have 
  !    "Dummy" appended.
  ! 2. A procedure pointer that is initially set to null.  This pointer is
  !    pointed to the appropriate subroutine later on (e.g. EOSOilInit())
  ! 3. An interface for derivative/non-derivative versions

  ! procedure pointers

  procedure(EOSOilViscosityDummy), pointer :: EOSOilViscosityPtr => null()
  procedure(EOSOilDensityDummy), pointer :: EOSOilDensityPtr => null()
  procedure(EOSOilEnthalpyDummy), pointer :: EOSOilEnthalpyPtr => null()
  procedure(EOSOilDensityEnergyDummy), pointer :: &
    EOSOilDensityEnergyPtr => null()

  ! these should be define as astract interfaces, because there are no   
  ! precedures named as xxxDummy that have such interfaces
  interface
    subroutine EOSOilViscosityDummy(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! oil pressure [Pa]
      PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]  
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Vis     ! oil viscosity 
      PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt temperature
      PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt Pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSOilViscosityDummy
    subroutine EOSOilDensityDummy(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSOilDensityDummy
    subroutine EOSOilEnthalpyDummy(T,P,deriv,H,dH_dT,dH_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSOilEnthalpyDummy
    subroutine EOSOilDensityEnergyDummy(T,P,deriv,Rho,dRho_dT,dRho_dP, &
                                        H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscReal, intent(out) :: U       ! internal energy [J/kmol]
      PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
      PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSOilDensityEnergyDummy
  end interface 

  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  interface EOSOilViscosity
    procedure EOSOilViscosityNoDerive
  ! procedure EOSOilViscosityDerive
  end interface
  interface EOSOilDensity
    procedure EOSOilDensityNoDerive
    procedure EOSOilDensityDerive
  end interface
  interface EOSOilEnthalpy
    procedure EOSOilEnthalpyNoDerive
   ! procedure EOSOilEnthalpyDerive
  end interface
  interface EOSOilDensityEnergy
    procedure EOSOilDenEnergyNoDerive
   ! procedure EOSOilDenEnergyDerive
  end interface  


  public :: EOSOilInit, &
            EOSOilViscosity, &
            EOSOilDensity, &
            EOSOilEnthalpy, & 
            EOSOilDensityEnergy

  public :: EOSOilSetViscosityConstant, &
            EOSOilSetDensityConstant, &
            EOSOilSetEnthalpyConstant, &
            EOSOilSetEnthalpyLinearTemp

contains

! ************************************************************************** !

subroutine EOSOilInit()

  implicit none
  
  constant_density = UNINITIALIZED_DOUBLE
  constant_viscosity = UNINITIALIZED_DOUBLE
  constant_enthalpy = UNINITIALIZED_DOUBLE
  fmw_oil = FMWOIL

  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms

  ! to be replaced with the tough model for which parameters 
  ! are required, or another simple model
  ! once decided the default model, add the routine that verify the input 
  EOSOilViscosityPtr => EOSOilViscosityConstant
  EOSOilDensityPtr => EOSOilDensityConstant
  EOSOilEnthalpyPtr => EOSOilEnthalpyConstant  
  
end subroutine EOSOilInit


! ************************************************************************** !

subroutine EOSOilSetViscosityConstant(viscosity)

  implicit none
  
  PetscReal :: viscosity
  
  constant_viscosity = viscosity  
  EOSOilViscosityPtr => EOSOilViscosityConstant
  
end subroutine EOSOilSetViscosityConstant

! ************************************************************************** !

subroutine EOSOilSetDensityConstant(density)

  implicit none
  
  PetscReal :: density
  
  constant_density = density  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilDensityPtr => EOSOilDensityConstant
  
end subroutine EOSOilSetDensityConstant

! ************************************************************************** !

subroutine EOSOilSetEnthalpyConstant(enthalpy)

  implicit none
  
  PetscReal :: enthalpy
  
  constant_enthalpy = enthalpy  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilEnthalpyPtr => EOSOilEnthalpyConstant
  
end subroutine EOSOilSetEnthalpyConstant

! ************************************************************************** !

subroutine EOSOilSetEnthalpyLinearTemp(specific_heat)

  implicit none
  
  PetscReal :: specific_heat 
  
  constant_sp_heat = specific_heat  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilEnthalpyPtr => EOSOilEnthalpyLinearTemp

  !write(*,*) "I am in EOS oil linear set up"  

end subroutine EOSOilSetEnthalpyLinearTemp

! ************************************************************************** !

subroutine EOSOilViscosityConstant(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]  
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Vis     ! oil viscosity 
  PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt temperature
  PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt Pressure
  PetscErrorCode, intent(out) :: ierr

  Vis = constant_viscosity

  dVis_dT = 0.0d0
  dVis_dT = 0.0d0
  
end subroutine EOSOilViscosityConstant

! ************************************************************************** !

subroutine EOSOilViscosityNoDerive(T,P,Rho,Vis,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]  
  PetscReal, intent(out) :: Vis     ! oil viscosity 
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2
  
  call EOSOilViscosityPtr(T,P,Rho,PETSC_FALSE,Vis,dum1,dum2,ierr)
  
end subroutine EOSOilViscosityNoDerive

! ************************************************************************** !

subroutine EOSOilDensityConstant(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr

  Rho = constant_density ! kmol/m^3

  dRho_dT = 0.d0
  dRho_dP = 0.d0

end subroutine EOSOilDensityConstant

! ************************************************************************** !

subroutine EOSOilDensityNoDerive(T,P,Rho,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2
  
  ! derivatives are so cheap, just compute them
  call EOSOilDensityPtr(T, P, PETSC_FALSE, Rho, dum1, dum2, ierr)
  
end subroutine EOSOilDensityNoDerive

! ************************************************************************** !

subroutine EOSOilDensityDerive(T,P,Rho,dRho_dT,dRho_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  call EOSOilDensityPtr(T,P,PETSC_TRUE,Rho,dRho_dT,dRho_dP,ierr)
                          
end subroutine EOSOilDensityDerive

! ************************************************************************** !

subroutine EOSOilEnthalpyConstant(T,P,deriv,H,dH_dT,dH_dP,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  H = constant_enthalpy ! J/kmol

  dH_dP = 0.d0
  dH_dT = 0.d0

end subroutine EOSOilEnthalpyConstant

! ************************************************************************** !

subroutine EOSOilEnthalpyLinearTemp(T,P,deriv,H,dH_dT,dH_dP,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  H = constant_sp_heat * T * fmw_oil ! J/(kg °C) °C * Kg/Kmol = J/Kmol 

  dH_dT = UNINITIALIZED_DOUBLE
  dH_dP = UNINITIALIZED_DOUBLE

  if(deriv) then
    dH_dP = 0.d0
    dH_dT = constant_sp_heat * fmw_oil
  end if

end subroutine EOSOilEnthalpyLinearTemp

! ************************************************************************** !

subroutine EOSOilEnthalpyNoDerive(T,P,H,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2

  call EOSOilEnthalpyPtr(T,P,PETSC_FALSE,H,dum1,dum2,ierr)

end subroutine EOSOilEnthalpyNoDerive

! ************************************************************************** !

subroutine EOSOilDensityEnergyTOilIms(T,P,deriv,Rho,dRho_dT,dRho_dP, &
                                    H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  call EOSOilDensityPtr(T,P,deriv,Rho,dRho_dT,dRho_dP,ierr)
  call EOSOilEnthalpyPtr(T,P,deriv,H,dH_dT,dH_dP,ierr)  

  U = H - P/Rho

  dU_dT = UNINITIALIZED_DOUBLE
  dU_dP = UNINITIALIZED_DOUBLE

  if(deriv) then
    print*, "EOSOilDensityEnergyTOilIms - U derivatives not supported"
    stop  
  end if   


end subroutine EOSOilDensityEnergyTOilIms

! ************************************************************************** !

subroutine EOSOilDenEnergyNoDerive(T,P,Rho,H,U,ierr)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2, dum3, dum4, dum5, dum6
 
  call EOSOilDensityEnergyPtr(T,P,PETSC_FALSE,Rho,dum1,dum2, &
                              H,dum3,dum4,U,dum5,dum6,ierr) 
 

end subroutine EOSOilDenEnergyNoDerive

! ************************************************************************** !

end module EOS_Oil_module
