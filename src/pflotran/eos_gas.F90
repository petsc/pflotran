module EOS_Gas_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "finclude/petscsys.h"

  ! module variables
  PetscReal :: constant_density
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity

  ! exponential
  PetscReal :: exponent_reference_density
  PetscReal :: exponent_reference_pressure
  PetscReal :: exponent_gas_compressibility

  ! In order to support generic EOS subroutines, we need the following:
  ! 1. An interface declaration that defines the argument list (best to have 
  !    "Dummy" appended.
  ! 2. A procedure pointer that is initially set to null.  This pointer is
  !    pointed to the appropriate subroutine later on (e.g. EOSGasInit())
  ! 3. An interface for derivative/non-derivative versions

  ! procedure pointer declarations
  procedure(EOSGasViscosityDummy), pointer :: EOSGasViscosityPtr => null()
  procedure(EOSGasDensityEnergyDummy), pointer :: &
    EOSGasDensityEnergyPtr => null()
  procedure(EOSGasDensityDummy), pointer :: EOSGasDensityPtr => null()
  procedure(EOSGasEnergyDummy), pointer :: EOSGasEnergyPtr => null()
  
  ! interface blocks
  interface
    subroutine EOSGasViscosityDummy(T, P_comp, P_gas, Rho_comp, V_mix, ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
      PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
      PetscReal, intent(in) :: Rho_comp ! air density [C]
      PetscReal, intent(out) :: V_mix   ! mixture viscosity
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSGasViscosityDummy
    subroutine EOSGasDensityDummy(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSGasDensityDummy
    subroutine EOSGasEnergyDummy(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscReal, intent(out) :: U       ! internal energy [J/kmol]
      PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
      PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSGasEnergyDummy
    subroutine EOSGasDensityEnergyDummy(T,P,Rho_gas,dRho_dT,dRho_dP, &
                                        H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscReal, intent(out) :: U       ! internal energy [J/kmol]
      PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
      PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSGasDensityEnergyDummy
  end interface
  
  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  interface EOSGasViscosity
    procedure EOSGasViscosityNoDerive
!    procedure EOSGasViscosityDerive
  end interface
  interface EOSGasDensity
    procedure EOSGasDensityNoDerive
    procedure EOSGasDensityDerive
  end interface
  interface EOSGasEnergy
    procedure EOSGasEnergyNoDerive
    procedure EOSGasEnergyDerive
  end interface
  interface EOSGasDensityEnergy
    procedure EOSGasDenEnthNoDerive
    procedure EOSGasDenEnthDerive
  end interface

  ! the "public" definition that makes subroutines visible outside.
  public :: EOSGasInit, &
            EOSGasVerify, &
            EOSGasViscosity, &
            EOSGasDensity, &
            EOSGasEnergy, &
            EOSGasDensityEnergy
            
  public :: EOSGasSetDensityIdeal, &
            EOSGasSetEnergyIdeal, &
!            EOSGasSetDensityRKS, &
            EOSGasSetDensityConstant, &
            EOSGasSetEnergyConstant, &
            EOSGasSetViscosityConstant, &
            EOSGasHenry_air_noderiv, &
            EOSGasHenry_air
 
  contains

! ************************************************************************** !

subroutine EOSGasInit()

  implicit none
  
  constant_density = -999.d0
  constant_viscosity = -999.d0
  constant_enthalpy = -999.d0

  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityIdeal
  EOSGasEnergyPtr => EOSGasEnergyIdeal
  EOSGasViscosityPtr => EOSGasViscosity1
  
end subroutine EOSGasInit

! ************************************************************************** !

subroutine EOSGasVerify(ierr,error_string)

  implicit none
  
  PetscErrorCode, intent(out) :: ierr
  character(len=MAXSTRINGLENGTH), intent(out) :: error_string
  
  ierr = 0
  
  error_string = ''
  if ((associated(EOSGasDensityPtr,EOSGasDensityIdeal) .and. &
        constant_density > -998.d0) .or. &
!      (associated(EOSGasDensityPtr,EOSGasDensityRKS) .and. &
!        constant_density > -998.d0) .or. &
      (associated(EOSGasEnergyPtr,EOSGasEnergyIdeal) .and. &
        constant_enthalpy > -998.d0) &
     ) then
    ierr = 1
  endif

  if (associated(EOSGasDensityPtr,EOSGasDensityConstant) .and. &
      constant_density < -998.d0) then
    error_string = trim(error_string) // &
      ' CONSTANT density not set.'
    ierr = 1
  endif
  
  if (associated(EOSGasEnergyPtr,EOSGasEnergyConstant) .and. &
      constant_enthalpy < -998.d0) then
    error_string = trim(error_string) // &
      ' CONSTANT enthalpy not set.'
    ierr = 1
  endif
  
  if ((associated(EOSGasViscosityPtr, &
                  EOSGasViscosityConstant) .and. &
       constant_viscosity < -998.d0) .or. &
      (associated(EOSGasViscosityPtr, &
                  EOSGasViscosity1) .and. &
       constant_viscosity > -998.d0)) then
    ierr = 1
  endif
  
end subroutine EOSGasVerify

! ************************************************************************** !

subroutine EOSGasSetDensityIdeal()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityIdeal
  
end subroutine EOSGasSetDensityIdeal

! ************************************************************************** !

subroutine EOSGasSetEnergyIdeal()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasEnergyPtr => EOSGasEnergyIdeal
  
end subroutine EOSGasSetEnergyIdeal

! ************************************************************************** !

subroutine EOSGasSetDensityConstant(density)

  implicit none
  
  PetscReal :: density
  
  constant_density = density  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityConstant
  
end subroutine EOSGasSetDensityConstant

! ************************************************************************** !

subroutine EOSGasSetEnergyConstant(enthalpy)

  implicit none
  
  PetscReal :: enthalpy
  
  constant_enthalpy = enthalpy  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasEnergyPtr => EOSGasEnergyConstant
  
end subroutine EOSGasSetEnergyConstant

! ************************************************************************** !

subroutine EOSGasSetViscosityConstant(viscosity)

  implicit none
  
  PetscReal :: viscosity
  
  constant_viscosity = viscosity  
  EOSGasViscosityPtr => EOSGasViscosityConstant
  
end subroutine EOSGasSetViscosityConstant


! ************************************************************************** !

subroutine EOSGasViscosityNoDerive(T, P_comp, P_gas, Rho_comp, V_mix, ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscErrorCode, intent(out) :: ierr
  
  call EOSGasViscosityPtr(T, P_comp, P_gas, Rho_comp, V_mix, ierr)
  
end subroutine EOSGasViscosityNoDerive

! ************************************************************************** !
#if 0
subroutine EOSGasViscosityDerive(T, P_comp, P_gas, Rho_comp, V_mix, ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscErrorCode, intent(out) :: ierr
 
  ! not yet supported
  print *, 'EOSGasViscosityDerive() not yet supported.'
  stop
  
end subroutine EOSGasViscosityDerive
#endif
! ************************************************************************** !

subroutine EOSGasViscosity1(T, P_comp, P_gas, Rho_comp, V_mix, ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscErrorCode, intent(out) :: ierr
  
  !geh: copied from gas_eos_mod.F90
  !
  ! REFERENCES
  ! THIS ROUTINE IS LARGELY ADAPTED FROM THE TOUGH CODE.
  ! this routine computes the viscosity of vapor-air mixtures.
  ! it uses a modified version of a formulation based on kinetic
  ! gas theory, as given by j.o. hirschfelder, c.f. curtiss, and
  ! r.b. bird, molecular theory of gases and liquids, john wiley
  ! & sons, 1954, pp. 528-530.
  ! the modification made to the hirschfelder et al. expressions is
  ! that for vapor viscosity accurate (empirical) values are used,
  ! rather than the first order expression of kinetic theory.
  ! the formulation matches experimental data on viscosities of
  ! vapor-air mixtures in the temperature range from 100 to 150
  ! deg. c, for all compositions, to better than 4%.
  ! 
  !PetscReal, intent(in) :: t     ! [C]
  !PetscReal, intent(in) :: p_air ! [Pa]
  !PetscReal, intent(in) :: p_gas ! [Pa]
  !PetscReal, intent(in) :: d_air ! [kmol/m^3]
  !PetscReal, intent(out) :: visg ! [Pa-s]

  PetscReal ::  fair,fwat,cair,cwat
  PetscReal :: p_air, d_air, visg

  data  fair,   fwat,    cair,  cwat &
        /97.d0, 363.d0, 3.617d0, 2.655d0/
 
  PetscReal fmix,cmix,d,xga,xg1,tk,trd1,trd3,ome1,ome3,ard,fmw3,vis1, &
          v1,vs,vis2,vis3,z1,g,h,e,z2,z3


!c======================================================================

  p_air = P_comp
  d_air = Rho_comp
          
  fmix = sqrt (fair*fwat)
  cmix = (cair+cwat)*0.5d0

!      do k = 1,nb
!       if (iphas(k).eq.2 .or. iphas(k).eq.0) then

      d   = d_air *FMWAIR       
      xga = p_air / p_gas ! for debug, set x constant
      xg1 = 1.D0 - xga
      tk  = t +273.15d0

      trd1 = tk/fair
      trd3 = tk/fmix
      ome1 = (1.188d0-0.051d0*trd1)/trd1
      ome3 = (1.480d0-0.412d0*log(trd3))/trd3
      ard  = 1.095d0/trd3
      fmw3 = 2.d0*FMWAIR*FMWH2O/(FMWAIR+FMWH2O)
      vis1 = 266.93d-7*sqrt(FMWAIR*trd1*fair)/(cair*cair*ome1*trd1)
 
      v1 = .407d0*t +80.4d0
      if (t .le.350.d0) then
        vs = 1.d-7*(v1-d*(1858.d0-5.9d0*t )*1.d-3)
      else
!             if (t .gt.350.d0) 
!cpcl .      vs = 1.d-7*(v1 + 0.353d0*d + 676.5d-6*d**2 + 102.1d-9*d**3)
        vs = 1.d-7*(v1 + (0.353d0 + (676.5d-6 + 102.1d-9*d)*d)*d)
      endif

      vis2 = 10.d0*vs
      vis3 = 266.93d-7*sqrt(fmw3*trd3*fmix)/(cmix*cmix*ome3*trd3)
      z1   = xga*xga/vis1+2.d0*xg1*xga/vis3+xg1*xg1/vis2
      g    = xga*xga*FMWAIR/FMWH2O
      h    = xg1*xg1*FMWH2O/FMWAIR
      e    = (2.d0*xga*xg1*FMWAIR*FMWH2O/fmw3**2)*vis3/(vis1*vis2)
      z2   = 0.6d0*ard*(g/vis1+e+h/vis2)
      z3   = 0.6d0*ard*(g+e*(vis1+vis2)-2.d0*xga*xg1+h)
      visg  = (1.d0+z3)/(z1+z2)*.1d0 
      
  V_mix = visg
  
end subroutine EOSGasViscosity1

! ************************************************************************** !

subroutine EOSGasViscosityConstant(T, P_comp, P_gas, Rho_comp, V_mix, ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P_comp   ! air pressure [Pa]
  PetscReal, intent(in) :: P_gas    ! gas pressure [Pa]
  PetscReal, intent(in) :: Rho_comp ! air density [C]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscErrorCode, intent(out) :: ierr

  V_mix = constant_viscosity
  
end subroutine EOSGasViscosityConstant

! ************************************************************************** !

subroutine EOSGasDensityNoDerive(T,P,Rho_gas,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2
  
  ! derivatives are so cheap, just compute them
  call EOSGasDensityPtr(T,P,Rho_gas,dum1,dum2,ierr)
  
end subroutine EOSGasDensityNoDerive

! ************************************************************************** !

subroutine EOSGasDensityDerive(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  call EOSGasDensityPtr(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)
                          
end subroutine EOSGasDensityDerive

! ************************************************************************** !

subroutine EOSGasEnergyNoDerive(T,P,H,U,ierr)

  implicit none
  
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2, dum3, dum4
  
  call EOSGasEnergyPtr(T,P,H,dum1,dum2,U,dum3,dum4,ierr)
  
end subroutine EOSGasEnergyNoDerive

! ************************************************************************** !

subroutine EOSGasEnergyDerive(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
    
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  call EOSGasEnergyPtr(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
  
end subroutine EOSGasEnergyDerive

! ************************************************************************** !

subroutine EOSGasDenEnthNoDerive(T,P,Rho_gas,H,U,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2, dum3, dum4, dum5, dum6
  
  call EOSGasDensityEnergyPtr(T,P,Rho_gas,dum1,dum2, &
                              H,dum3,dum4,U,dum5,dum6,ierr)
  
end subroutine EOSGasDenEnthNoDerive

! ************************************************************************** !

subroutine EOSGasDenEnthDerive(T,P,Rho_gas,dRho_dT,dRho_dP, &
                               H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  call EOSGasDensityEnergyPtr(T,P,Rho_gas,dRho_dT,dRho_dP, &
                              H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
  
end subroutine EOSGasDenEnthDerive

! ************************************************************************** !

subroutine EOSGasDensityEnergyGeneral(T,P,Rho_gas,dRho_dT,dRho_dP, &
                                      H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
                                      
  implicit none
  
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressuret
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  call EOSGasDensityPtr(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)
  call EOSGasEnergyPtr(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)  
  
end subroutine EOSGasDensityEnergyGeneral

! ************************************************************************** !

subroutine EOSGasDensityIdeal(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr

  PetscReal, parameter:: Rg = 8.31415 
  PetscReal  T_kelvin

  T_kelvin = T + 273.15d0
  Rho_gas = P / T_kelvin / Rg * 1.d-3 ! mol/m^3 -> kmol/m^3

  dRho_dP =  Rho_gas / P
  dRho_dT = -Rho_gas / T_kelvin

end subroutine EOSGasDensityIdeal

! ************************************************************************** !

subroutine EOSGasDensityRKS(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)
! Redlich-Kwong-Soave (RKS) equation of state used in BRAGFLO.  See
! Soave, Giorgio, 1972, "Equilibrium constants from a modified Redlich-Kwong
! equation of state", Chem. Eng. Sci., V27, pp 1197-1203.
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(inout) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr

  
  
  PetscReal, parameter:: Rg = 8.31415 
  PetscReal  T_kelvin

  T_kelvin = T + 273.15d0
  
  ! to be completed later
  print *, 'RKS gas density not yet implemented.'
  stop
  
end subroutine EOSGasDensityRKS

! ************************************************************************** !

subroutine EOSGasEnergyIdeal(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
    
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  PetscReal, parameter:: Rg = 8.31415 
  ! Cpg units: J/mol-K
  PetscReal, parameter:: Cv_air = 20.85 ! head capacity wiki
  PetscReal  T_kelvin

  T_kelvin = T + 273.15d0
  H = Cv_air * T_kelvin * 1.d3  ! J/mol -> J/kmol
  U = (Cv_air - Rg) * T_kelvin * 1.d3 ! J/mol -> J/kmol

  dH_dP = 0.d0
  dH_dT = Cv_air * 1.d3
  dU_dP = 0.d0
  dU_dT = (Cv_air - Rg) * 1.d3
    
end subroutine EOSGasEnergyIdeal

! ************************************************************************** !

subroutine EOSGasDensityConstant(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)

  implicit none
  
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  Rho_gas = constant_density ! kmol/m^3

  dRho_dT = 0.d0
  dRho_dP = 0.d0

end subroutine EOSGasDensityConstant

! ************************************************************************** !

subroutine EOSGasEnergyConstant(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
    
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  H = constant_enthalpy ! J/kmol
  
  dH_dP = 0.d0
  dH_dT = 0.d0
  
  print *, 'Calcuation of internal gas energy not set up in ' // &
    'EOSGasEnergyConstant.'
  stop
  
end subroutine EOSGasEnergyConstant

! ************************************************************************** !

subroutine EOSGasHenry_air_noderiv(p,tc,ps,Henry)
! Calculate Henry Coefficient for N2
! t in K
! Henry have the same unit as p and ps, then make it dimensionless by
! devide it with p

    implicit none
    PetscReal,intent(in) ::  p,tc,ps
    PetscReal,intent(out)::  Henry

    PetscReal  Tr,tao,tmp,t
    PetscReal, parameter :: a=-9.67578, b=4.72162, c=11.70585
    PetscReal, parameter :: Tcl=647.096 ! H2O critical temp(K) from IAPWS(1995b)

    t=tc+273.15D0
    Tr=t/Tcl
    tao=1.D0-Tr
    tmp= a/Tr + B * tao**0.355/Tr + c * (Tr**(-0.41)) * exp(tao)
    Henry=exp(tmp)*ps

   return 
end subroutine EOSGasHenry_air_noderiv

! ************************************************************************** !

subroutine EOSGasHenry_air(p,tc,ps,ps_p,ps_t,Henry,Henry_p,Henry_t)
   implicit none
    PetscReal,intent(in) ::  p,tc,ps,ps_p,ps_t
    PetscReal,intent(out)::  Henry,Henry_p,Henry_t
! note t/K, p/Pa, Henry/Pa 

    PetscReal  Tr,tao,tmp,t
    PetscReal, parameter :: a=-9.67578, b=4.72162, c=11.70585
    PetscReal, parameter :: Tcl=647.096 ! H2O critical temp from IAPWS(1995b)

    t=tc+273.15D0
    Tr=t/Tcl
    tao=1.D0-Tr
    tmp= a/Tr + b * tao**0.355/Tr + c * (Tr**(-0.41)) * exp(tao)
    Henry=exp(tmp)*ps

    tmp =((-a/Tr+b*(-0.355*tao**(-0.645)-tao**0.355/Tr))/Tr - &
         c*exp(tao)*(tao**(-.41))*(0.41/Tr-1.))/Tcl
    Henry_t=Henry*(tmp +ps_t/ps)
    Henry_p=ps_p*Henry/ps

  
   return 
end subroutine EOSGasHenry_air

end module EOS_Gas_module
