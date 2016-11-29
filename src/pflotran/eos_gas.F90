module EOS_Gas_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

  ! module variables
  PetscReal :: fmw_gas           !kg/Kmol
  PetscReal :: constant_density
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity
  PetscReal :: constant_henry
  PetscBool :: hydrogen
  PetscReal :: Tc
  PetscReal :: Pc
  PetscReal :: acentric
  PetscReal :: coeff_a
  PetscReal :: coeff_b
  
  ! exponential
  PetscReal :: exponent_reference_density
  PetscReal :: exponent_reference_pressure
  PetscReal :: exponent_gas_compressibility

  ! This is the offset added to temperature [C] used to calculate the energy
  ! equation of state.  Swithing between 0. and 273.15 greatly changes results.
#if defined(MATCH_TOUGH2)
  PetscReal, parameter :: T_energy_offset = 0.d0
#else
  PetscReal, parameter :: T_energy_offset = 273.15d0
#endif

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
  procedure(EOSGasHenryDummy), pointer :: EOSGasHenryPtr => null()
  
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
    subroutine EOSGasHenryDummy(T,Psat,Hc)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: Psat     ! saturation pressure
      PetscReal, intent(out) :: Hc      ! Henry's constant
    end subroutine EOSGasHenryDummy
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
  interface EOSGasHenry
    procedure EOSGasHenryNoDerive
!    procedure EOSGasHenryDerive
  end interface

  ! the "public" definition that makes subroutines visible outside.
  public :: EOSGasInit, &
            EOSGasVerify, &
            EOSGasViscosity, &
            EOSGasDensity, &
            EOSGasEnergy, &
            EOSGasDensityEnergy, &
            EOSGasHenry, &
            EOSGasInputRecord, &
            EOSGasTest
            
  public :: EOSGasSetDensityIdeal, &
            EOSGasSetEnergyIdeal, &
            EOSGasSetDensityRKS, &
            EOSGasSetDensityPRMethane, &
            EOSGasSetEnergyIdealMethane, &
            EOSGasSetFMWConstant, &
            EOSGasGetFMW, & 
            EOSGasSetDensityConstant, &
            EOSGasSetEnergyConstant, &
            EOSGasSetViscosityConstant, &
            EOSGasSetHenry, &
            EOSGasSetHenryConstant
 
  contains

! ************************************************************************** !

subroutine EOSGasInit()

  implicit none

  fmw_gas = UNINITIALIZED_DOUBLE  
  constant_density = UNINITIALIZED_DOUBLE
  constant_viscosity = UNINITIALIZED_DOUBLE
  constant_enthalpy = UNINITIALIZED_DOUBLE
  constant_henry = UNINITIALIZED_DOUBLE

  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityIdeal
  EOSGasEnergyPtr => EOSGasEnergyIdeal
  EOSGasViscosityPtr => EOSGasViscosity1
  EOSGasHenryPtr => EOSGasHenry_air_noderiv
  
end subroutine EOSGasInit

! ************************************************************************** !

subroutine EOSGasVerify(ierr,error_string)

  implicit none
  
  PetscErrorCode, intent(out) :: ierr
  character(len=MAXSTRINGLENGTH), intent(out) :: error_string
  
  ierr = 0
  
  error_string = ''
  if ((associated(EOSGasDensityPtr,EOSGasDensityIdeal) .and. &
        Initialized(constant_density)) .or. &
      (associated(EOSGasDensityPtr,EOSGasDensityRKS) .and. &
        Initialized(constant_density)) .or. &
      (associated(EOSGasEnergyPtr,EOSGasEnergyIdeal) .and. &
        Initialized(constant_enthalpy)) &
     ) then
    ierr = 1
  endif

  if (associated(EOSGasDensityPtr,EOSGasDensityConstant) .and. &
      Uninitialized(constant_density)) then
    error_string = trim(error_string) // &
      ' CONSTANT density not set.'
    ierr = 1
  endif
  
  if (associated(EOSGasEnergyPtr,EOSGasEnergyConstant) .and. &
      Uninitialized(constant_enthalpy)) then
    error_string = trim(error_string) // &
      ' CONSTANT enthalpy not set.'
    ierr = 1
  endif
  
  if ((associated(EOSGasViscosityPtr, &
                  EOSGasViscosityConstant) .and. &
       Uninitialized(constant_viscosity)) .or. &
      (associated(EOSGasViscosityPtr, &
                  EOSGasViscosity1) .and. &
       Initialized(constant_viscosity))) then
    ierr = 1
  endif
  
  if (associated(EOSGasHenryPtr, &
                 EOSGasHenryConstant) .and. &
      Uninitialized(constant_henry)) then
    error_string = trim(error_string) // " Henry's constant not set"
    ierr = 1
  endif
  
  if (associated(EOSGasDensityPtr,EOSGasDensityRKS)) then
    if (hydrogen) then
      ! Assign default hydrogen parameters if not assigned
      if (Uninitialized(Tc)) then
        Tc = 41.67d0
        ierr = 5
        error_string = trim(error_string) // " Tc"
      endif
      if (Uninitialized(Pc)) then
        Pc = 2.1029d6
        ierr = 5
        error_string = trim(error_string) // " Pc"
      endif
      if (Uninitialized(coeff_a)) then
        coeff_a = 0.42747d0
        ierr = 5
        error_string = trim(error_string) // " omega_a"
      endif
      if (Uninitialized(coeff_b)) then
        coeff_b = 0.08664d0
        ierr = 5
        error_string = trim(error_string) // " omega_b"
      endif
    else
      if (Uninitialized(Tc) .or. Uninitialized(coeff_a) .or. &
        Uninitialized(Pc) .or. Uninitialized(coeff_b) .or. &
        Uninitialized(acentric)) then
          error_string = trim(error_string) // &
          " RKS parameters not set for non-hydrogen gas"
          ierr = 1
      endif
    endif
  endif

      
end subroutine EOSGasVerify

! ************************************************************************** !

subroutine EOSGasSetDensityIdeal()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityIdeal
  
end subroutine EOSGasSetDensityIdeal

! ************************************************************************** !

subroutine EOSGasSetDensityRKS(h,rks_tc,rks_pc,rks_acen,rks_omegaa,rks_omegab)

  implicit none

  PetscBool :: h
  PetscReal :: rks_tc
  PetscReal :: rks_pc
  PetscReal :: rks_acen
  PetscReal :: rks_omegaa
  PetscReal :: rks_omegab
  
  hydrogen = h
  Tc = rks_tc
  Pc = rks_pc
  acentric = rks_acen
  coeff_a = rks_omegaa
  coeff_b = rks_omegab
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityRKS
  
end subroutine EOSGasSetDensityRKS

! ************************************************************************** !

subroutine EOSGasSetDensityPRMethane()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasDensityPtr => EOSGasDensityPRMethane
  
end subroutine EOSGasSetDensityPRMethane

! ************************************************************************** !

subroutine EOSGasSetEnergyIdeal()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasEnergyPtr => EOSGasEnergyIdeal
  
end subroutine EOSGasSetEnergyIdeal

! ************************************************************************** !

subroutine EOSGasSetEnergyIdealMethane()

  implicit none
  
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyGeneral
  EOSGasEnergyPtr => EOSGasEnergyIdealMethane
  
end subroutine EOSGasSetEnergyIdealMethane

! ************************************************************************** !

subroutine EOSGasSetFMWConstant(input_fmw_gas)

  implicit none

  PetscReal :: input_fmw_gas
  
  fmw_gas = input_fmw_gas

end subroutine EOSGasSetFMWConstant

! ************************************************************************** !

function EOSGasGetFMW()

  implicit none
  
  PetscReal :: EOSGasGetFMW
  
  EOSGasGetFMW = fmw_gas
  
end function EOSGasGetFMW

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

subroutine EOSGasSetHenry()

  implicit none
  
  EOSGasHenryPtr => EOSGasHenry_air_noderiv
  
end subroutine EOSGasSetHenry

! ************************************************************************** !

subroutine EOSGasSetHenryConstant(henrys_constant)

  implicit none
  
  PetscReal :: henrys_constant
  
  constant_henry = henrys_constant
  EOSGasHenryPtr => EOSGasHenryConstant
  
end subroutine EOSGasSetHenryConstant

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

  PetscReal  T_kelvin

  T_kelvin = T + 273.15d0
  Rho_gas = P / T_kelvin / IDEAL_GAS_CONSTANT * 1.d-3 ! mol/m^3 -> kmol/m^3

  dRho_dP =  Rho_gas / P
  dRho_dT = -Rho_gas / T_kelvin

end subroutine EOSGasDensityIdeal

! ************************************************************************** !

subroutine EOSGasDensityRKS(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)
! Redlich-Kwong-Soave (RKS) equation of state used in BRAGFLO.  See
! Soave, Giorgio, 1972, "Equilibrium constants from a modified Redlich-Kwong
! equation of state", Chem. Eng. Sci., V27, pp 1197-1203.
! current version is for hydrogen only.
!
! Author: Heeho Park
! Date: 05/08/14
!
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: T_kelvin, RT, alpha, a, B , a_RT, p_RT
  PetscReal :: b2, V, f, dfdV, dVd
  PetscInt :: i
  PetscReal :: coeff_alpha
  
  !for hydrogen
  ! American Petroleum Institute,
  ! "Technical Data Book-Petroleum Refining" (1977)
  
  !solver
  PetscReal :: coef(4)
  PetscInt, parameter :: maxit = 50
  dRho_dT = UNINITIALIZED_DOUBLE
  dRho_dP = UNINITIALIZED_DOUBLE
  
  T_kelvin = T + 273.15d0
  RT = IDEAL_GAS_CONSTANT * T_kelvin
  
  if (hydrogen) then
    ! suggested by Peter Lichtner
    ! this function honors alpha(Tc)=1 while closely fitting Graboski curve from
    ! "A Modified Soave Equation of State for Phase Equilibrium
    !  Calculations. 3. Systems Containing Hydrogen"
    alpha = EXP(0.340d0*(1-T_kelvin/Tc))
  else
    coeff_alpha = 0.48508d0 + acentric*(1.55171d0 - 0.15613*acentric)
    alpha = (1.d0 + coeff_alpha*(1.d0 - SQRT(T_Kelvin/Tc)))**2
  end if
  
  a = coeff_a * alpha * (IDEAL_GAS_CONSTANT * Tc)**2 / Pc
  b = coeff_b * IDEAL_GAS_CONSTANT * Tc / Pc
  
  a_RT = a / RT
  P_RT = P / RT
  coef(4) = P / RT
  coef(3) = 1.0d0
  coef(2) = a_RT - b - P_RT*b**2
  coef(1) = a_RT*b
  V = RT/P  ! initial guess
  
  !Newton iteration to find a Volume of gas
  do i = 1, maxit
    f = V*(V*(coef(4)*V - coef(3)) + coef(2)) - coef(1)
    dfdV = V*(3.0d0*coef(4)*V - 2.0d0*coef(3)) + coef(2)
    if (dfdV .ne. 0.d0) then
      dVd = F/dfdV
      V = V - dVd
      if (V .ne. 0.d0) then
        if (abs(dVd/V) .lt. 1.d-10) then
          Rho_gas = 1/V * 1d-3 ! mol/m^3 -> kmol/m^3
          exit
        end if
      else
        print *, 'Error: RKS gas density cannot be solved'
      end if 
    else
      print *, 'Error: Zero Slope in Equation of State'
    end if
  end do
  
end subroutine EOSGasDensityRKS

! ************************************************************************** !

subroutine EOSGasDensityPRMethane(T,P,Rho_gas,dRho_dT,dRho_dP,ierr)
! Peng-Robinson (PR) equation of state. This is the uncorrected version 
! sufficient for methane. For other higher hydrocarbons one might need
! the corrected version published in 1978.
! Reference: Peng, D. Y., and Robinson, D. B. (1976).
! A New Two-Constant Equation of State Industrial and Engineering Chemistry: 
! Fundamentals 15: 59–64 DOI:10.1021/i160057a011
! current version is for methane only.
! if omega is greater than 0.49 then alpha needs to be modified below
! Satish Karra, LANL
! 09/30/2014

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_gas ! gas density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative gas density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative gas density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: T_kelvin, RT, alpha,  T_r, Z
  PetscReal :: f, dfdZ, dZd
  PetscReal :: a, b
  PetscInt :: i
  PetscReal, parameter :: coeff_a = 0.45724
  PetscReal, parameter :: coeff_b = 0.0778
  
  ! for methane critical temp, pressure 
  PetscReal, parameter :: Tc = 191.15    ! in K
  PetscReal, parameter :: Pc = 4.641e6   ! in Pa
  PetscReal, parameter :: omega = 0.0115 ! acentric factor 
  
  !solver
  PetscReal :: coef(4)
  PetscInt, parameter :: maxit = 50
  
  T_kelvin = T + 273.15d0
  RT = IDEAL_GAS_CONSTANT*T_kelvin
  T_r = T_kelvin/Tc
  alpha = (1 + (0.3744 + 1.5422*omega - 0.26992*omega**2)*(1-T_r**(0.5)))**2
  
  a = coeff_a * alpha * (IDEAL_GAS_CONSTANT * Tc)**2 / Pc
  b = coeff_b * IDEAL_GAS_CONSTANT * Tc / Pc

  A = alpha * a / (RT)**2
  B = b * P / RT 

  coef(1) = 1.0d0
  coef(2) = B - 1.0d0
  coef(3) = A - 2 * B - 3 * B**2
  coef(4) = B**3 +  B**2 - A * B 
  Z = 1.0d0  ! initial guess for compressibility
  
  !Newton iteration to find a Volume of gas
  do i = 1, maxit
    f = coef(1) * Z**3 + coef(2) * Z**2 + coef(3) * Z + coef(4)
    dfdZ = coef(1) * 3 * Z**2  + coef(2) * 2 * Z + coef(3)
    if (dfdZ .ne. 0.d0) then
      dZd = f/dfdZ
      Z = Z - dZd
      if (Z .ne. 0.d0) then
        if (abs(dZd/Z) .lt. 1.d-10) then
          Rho_gas = 1.d0/Z * P/RT * 1d-3 ! mol/m^3 -> kmol/m^3
          exit
        end if
      else
        print *, 'Error: PR76 gas density cannot be solved'
      end if 
    else
      print *, 'Error: Zero Slope in PR76 Equation of State'
    end if
  end do
  
end subroutine EOSGasDensityPRMethane

! ************************************************************************** !

subroutine EOSGasFugacity(T,P,Rho_gas,phi,ierr)
! current version is for hydrogen only. See
! Soave, Giorgio, 1972, "Equilibrium constants from a modified Redlich-Kwong
! equation of state", Chem. Eng. Sci., V27, pp 1197-1203.

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(in) :: Rho_gas  ! gas density [kmol/m^3]
  PetscReal, intent(out) :: phi     ! fugacity coefficient
  PetscErrorCode, intent(out) :: ierr

  !for hydrogen
  ! American Petroleum Institute,
  ! "Technical Data Book-Petroleum Refining" (1977)
  PetscReal, parameter :: Tc = 41.67d0   ! critical temperature
  PetscReal, parameter :: Pc = 2.1029d6  ! critical pressure

  PetscReal :: T_kelvin, Z, V, Tr, LHS
  PetscReal :: A, B, alpha

  V = 1.d0 / Rho_gas * 1.d-3
  T_kelvin = T + 273.15d0
  Tr = T_kelvin/Tc  ! reduced temperature
  alpha = EXP(0.340d0*(1.d0-Tr))  ! alpha(T)

  A = 0.42747 * alpha * (P/Pc) / (T/Tc)**2
  B = 0.08664 * (P/Pc) / (T/Tc)
  Z = P*V / (IDEAL_GAS_CONSTANT*T)
  
  ! Fugacity Coefficient
  LHS = Z - 1.d0 - LOG(Z-B) - A/B * LOG((Z+B)/Z)
  phi = EXP(LHS)  ! dimensionless (Pa/Pa)
  
end subroutine EOSGasFugacity

! ************************************************************************** !

subroutine EOSGasEnergyIdealMethane(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
  
! Based on EOSGasEnergyIdeal. Modified for Methane
! Satish Karra, LANL
! 09/30/2014
  
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

  ! Cpg units: J/mol-K
  PetscReal, parameter :: Cp_methane = 35.69 ! at 298.15 K and 1 bar http://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Units=SI&Mask=1
  PetscReal :: T_energy
  PetscReal :: T_k

  ierr = 0
  ! T_energy is either T or T + 273.15
  ! do not change below
  T_energy = T + T_energy_offset
  T_k = T + 273.15d0
  H = Cp_methane * T_energy * 1.d3  ! J/mol -> J/kmol
  U = H - IDEAL_GAS_CONSTANT * T_k * 1.d3 ! J/mol -> J/kmol

  dH_dP = 0.d0
  dH_dT = Cp_methane * 1.d3
  dU_dP = 0.d0
  dU_dT = (Cp_methane - IDEAL_GAS_CONSTANT) * 1.d3
    
end subroutine EOSGasEnergyIdealMethane

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

  ! Cv_air units: J/mol-K
  PetscReal, parameter :: Cv_air = 20.85 ! heat capacity wiki [J/mol-K]
  PetscReal :: T_energy
  PetscReal :: T_k

  ! T_energy is either T or T + 273.15
  ! do not change below
  T_energy = T + T_energy_offset
  T_k = T + 273.15d0
  U = Cv_air * T_energy * 1.d3  ! J/mol -> J/kmol
  H = U + IDEAL_GAS_CONSTANT * T_k * 1.d3 ! J/mol -> J/kmol

  dU_dP = 0.d0
  dU_dT = Cv_air * 1.d3
  dH_dP = 0.d0
  dH_dT = dU_dT + IDEAL_GAS_CONSTANT * 1.d3
    
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
  
  ierr = 0
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
  
  PetscReal :: T_energy

  ierr = 0
  H = constant_enthalpy ! J/kmol
  ! T_energy is either T or T + 273.15
  ! do not change below
  T_energy = T + T_energy_offset
  U = H - IDEAL_GAS_CONSTANT * T_energy * 1.d3 ! J/mol -> J/kmol
  
  dH_dP = 0.d0
  dH_dT = 0.d0
  dU_dP = 0.d0
  dU_dT = -1.d0 * IDEAL_GAS_CONSTANT * 1.d3
  
  print *, 'Calculation of internal gas energy not set up properly in ' // &
    'EOSGasEnergyConstant.'
  stop  
  
end subroutine EOSGasEnergyConstant


! ************************************************************************** !

subroutine EOSGasHenryNoDerive(T,Psat,Hc)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: Psat     ! saturation pressure
  PetscReal, intent(out) :: Hc      ! Henry's constant
  
  call EOSGasHenryPtr(T,Psat,Hc)
                          
end subroutine EOSGasHenryNoDerive

! ************************************************************************** !

subroutine EOSGasHenry_air_noderiv(tc,ps,Henry)
! 
!   Calculates Henry's constant as a function of temperature [C], 
!   and saturation pressure [Pa].  
!
!   Fernandez-Prini, F., J. Alvarez and A. Harvey (2003) Henry's Constants and
!   Vapor-Liquid Distribution Constants for Gaseous Solutes in H2O and D2O at
!   High Temperatures, J. Phys. Chem. Ref. Data, Vol. 32, No. 2, Equation 15
!   with coefficients A,B,C from Table 3 for N2(g)

    implicit none
    
    PetscReal, intent(in) :: tc
    PetscReal, intent(in) :: ps
    PetscReal, intent(out) :: Henry

    PetscReal  Tr,tao,tmp,t
    PetscReal, parameter :: a=-9.67578d0, b=4.72162d0, c=11.70585d0
    PetscReal, parameter :: Tcl=647.096d0 ! H2O critical temp(K) from IAPWS(1995b)

    t = tc+273.15D0
    Tr = t/Tcl
    tao = 1.D0-Tr
    tmp = a/Tr + b * tao**0.355d0/Tr + c * (Tr**(-0.41d0)) * exp(tao)
    Henry = exp(tmp)*ps

end subroutine EOSGasHenry_air_noderiv

! ************************************************************************** !

subroutine EOSGasHenry_air(tc,ps,ps_p,ps_t,Henry,Henry_p,Henry_t)
! 
!   Calculates Henry's constant as a function of temperature [C], 
!   and saturation pressure [Pa].  
!
!   Fernandez-Prini, F., J. Alvarez and A. Harvey (2003) Henry's Constants and
!   Vapor-Liquid Distribution Constants for Gaseous Solutes in H2O and D2O at
!   High Temperatures, J. Phys. Chem. Ref. Data, Vol. 32, No. 2, Equation 15
!   with coefficients A,B,C from Table 3 for N2(g)

    implicit none
    
    PetscReal, intent(in) :: tc       ! [C]
    PetscReal, intent(in) :: ps       ! [Pa]
    PetscReal, intent(in) :: ps_p,ps_t
    PetscReal,intent(out) :: Henry,Henry_p,Henry_t

    PetscReal  Tr,tao,tmp,t
    PetscReal, parameter :: a=-9.67578d0, b=4.72162d0, c=11.70585d0
    PetscReal, parameter :: Tcl=647.096d0 ! H2O critical temp from IAPWS(1995b)

    t = tc+273.15D0
    Tr = t/Tcl
    tao = 1.D0-Tr
    tmp = a/Tr + b * tao**0.355d0/Tr + c * (Tr**(-0.41d0)) * exp(tao)
    Henry = exp(tmp)*ps

    tmp = ((-a/Tr+b*(-0.355d0*tao**(-0.645d0)-tao**0.355d0/Tr))/Tr - &
         c*exp(tao)*(tao**(-0.41d0))*(0.41d0/Tr-1.d0))/Tcl
    Henry_t = Henry*(tmp+ps_t/ps)
    Henry_p = ps_p*Henry/ps

end subroutine EOSGasHenry_air

! ************************************************************************** !

subroutine EOSGasHenryConstant(T,Psat,Hc)
   implicit none
    PetscReal, intent(in) :: T    ! temperature [C]
    PetscReal, intent(in) :: Psat ! saturation pressure [Pa]
    PetscReal, intent(out) :: Hc   ! Henry's constant [Pa]

    Hc = constant_henry
    
end subroutine EOSGasHenryConstant

! **************************************************************************** !

subroutine EOSGasInputRecord()
  ! 
  ! Prints ingested equation of state information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/04/2016
  ! 
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'GAS'
  
  ! gas density [kg/m^3]
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityConstant)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(word1,*) constant_density
    write(id,'(a)') 'constant, ' // adjustl(trim(word1)) // ' kg/m^3'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityRKS)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(id,'(a)') 'rsk'
    write(id,'(a29)',advance='no') 'critical temperature: '
    write(word1,*) Tc
    write(id,'(a)') adjustl(trim(word1)) // ' C'
    write(id,'(a29)',advance='no') 'critical pressure: '
    write(word1,*) Pc
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    if (.not.hydrogen) then
      write(id,'(a29)',advance='no') 'acentric factor: '
      write(word1,*) acentric
      write(id,'(a)') adjustl(trim(word1)) 
    endif
    write(id,'(a29)',advance='no') 'omega A: '
    write(word1,*) coeff_a
    write(id,'(a)') adjustl(trim(word1))
    write(id,'(a29)',advance='no') 'omega B: '
    write(word1,*) coeff_b
    write(id,'(a)') adjustl(trim(word1))
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityPRMethane)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(id,'(a)') 'pr methane'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasDensityPtr,EOSGasDensityIdeal)) then
    write(id,'(a29)',advance='no') 'gas density: '
    write(id,'(a)') 'default, ideal'
  endif
  
  ! gas viscosity [Pa-s]
  if (associated(EOSGasViscosityPtr,EOSGasViscosityConstant)) then
    write(id,'(a29)',advance='no') 'gas viscosity: '
    write(word1,*) constant_viscosity
    write(id,'(a)') 'constant, ' // trim(word1) // ' Pa-sec'
  endif
  if (associated(EOSGasViscosityPtr,EOSGasViscosity1)) then
    write(id,'(a29)',advance='no') 'gas viscosity: '
    write(id,'(a)') 'default'
  endif
  
  ! gas enthalpy [J/kmol]
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasEnergyPtr,EOSGasEnergyConstant)) then
    write(id,'(a29)',advance='no') 'gas enthalpy: '
    write(word1,*) constant_enthalpy
    write(id,'(a)') 'constant, ' // trim(word1) // 'J/kmol'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasEnergyPtr,EOSGasEnergyIdeal)) then
    write(id,'(a29)',advance='no') 'gas enthalpy: '
    write(id,'(a)') 'default, ideal'
  endif
  if (associated(EOSGasDensityEnergyPtr,EOSGasDensityEnergyGeneral) .and. &
      associated(EOSGasEnergyPtr,EOSGasEnergyIdealMethane)) then
    write(id,'(a29)',advance='no') 'gas enthalpy: '
    write(id,'(a)') 'ideal methane'
  endif
  
  ! henry's constant
  if (associated(EOSGasHenryPtr,EOSGasHenryConstant)) then
    write(id,'(a29)',advance='no') "henry's constant: "
    write(word1,*) constant_henry
    write(id,'(a)') 'constant, ' // trim(word1) 
  endif
  if (associated(EOSGasHenryPtr,EOSGasHenry_air_noderiv)) then
    write(id,'(a29)',advance='no') "henry's constant: "
    write(id,'(a)') 'default, air'
  endif
  
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  
end subroutine EOSGasInputRecord

! ************************************************************************** !

subroutine EOSGasTest(temp_low,temp_high,pres_low,pres_high, &
                        ntemp,npres,uniform_temp,uniform_pres,filename)

  use EOS_Water_module, only : EOSWaterSaturationPressure

  implicit none

  PetscReal :: temp_low
  PetscReal :: temp_high
  PetscReal :: pres_low
  PetscReal :: pres_high
  PetscInt :: npres
  PetscInt :: ntemp
  PetscBool :: uniform_temp
  PetscBool :: uniform_pres
  character(len=MAXWORDLENGTH) :: filename

  PetscReal, allocatable :: temp(:)
  PetscReal, allocatable :: pres(:)
  PetscReal, allocatable :: density_kg(:,:)
  PetscReal, allocatable :: enthalpy(:,:)
  PetscReal, allocatable :: internal_energy(:,:)
  PetscReal, allocatable :: viscosity(:,:)
  PetscReal, allocatable :: saturation_pressure_array(:)
  PetscReal, allocatable :: henry(:)
  PetscReal :: dum1, dum2, dum3, dum4
  PetscInt :: itemp, ipres
  PetscReal :: ln_low, ln_high
  PetscReal :: saturation_pressure
  PetscReal :: air_pressure
  PetscReal :: NaN
  character(len=MAXWORDLENGTH) :: eos_density_name
  character(len=MAXWORDLENGTH) :: eos_energy_name
  character(len=MAXWORDLENGTH) :: eos_viscosity_name
  character(len=MAXWORDLENGTH) :: eos_saturation_pressure_name
  character(len=MAXWORDLENGTH) :: eos_henry_name
  character(len=MAXSTRINGLENGTH) :: header

  PetscErrorCode :: ierr

  NaN = 0.d0
  NaN = 1.d0/NaN
  NaN = 0.d0*NaN

  allocate(temp(ntemp))
  temp = UNINITIALIZED_DOUBLE
  allocate(pres(ntemp))
  pres = UNINITIALIZED_DOUBLE
  allocate(density_kg(npres,ntemp))
  density_kg = UNINITIALIZED_DOUBLE
  allocate(viscosity(npres,ntemp))
  viscosity = UNINITIALIZED_DOUBLE
  allocate(enthalpy(npres,ntemp))
  enthalpy = UNINITIALIZED_DOUBLE
  allocate(internal_energy(npres,ntemp))
  internal_energy = UNINITIALIZED_DOUBLE
  allocate(saturation_pressure_array(ntemp))
  saturation_pressure_array = UNINITIALIZED_DOUBLE
  allocate(henry(ntemp))
  henry = UNINITIALIZED_DOUBLE

  if (uniform_pres) then
    do ipres = 1, npres
      pres(ipres) = (pres_high-pres_low)/dble(npres-1) * (ipres-1) + pres_low
    enddo
  else
    ln_high = log(pres_high)
    ln_low = log(pres_low)
    do ipres = 1, npres
      pres(ipres) = exp((ln_high-ln_low)/dble(npres-1) * (ipres-1) + ln_low)
    enddo
  endif

  if (uniform_temp) then
    do itemp = 1, ntemp
      temp(itemp) = (temp_high-temp_low)/dble(ntemp-1) * (itemp-1) + temp_low
    enddo
  else
    ln_high = log(temp_high)
    ln_low = log(temp_low)
    do itemp = 1, ntemp
      temp(itemp) = exp((ln_high-ln_low)/dble(ntemp-1) * (itemp-1) + ln_low)
    enddo
  endif

  ! density
  if (associated(EOSGasDensityPtr,EOSGasDensityConstant)) then
    eos_density_name = 'Constant'
  else if (associated(EOSGasDensityPtr,EOSGasDensityIdeal)) then
    eos_density_name = 'Ideal Gas Law'
  else if (associated(EOSGasDensityPtr,EOSGasDensityRKS)) then
    eos_density_name = 'Redlich-Kwong-Soave'
  else if (associated(EOSGasDensityPtr,EOSGasDensityPRMethane)) then
    eos_density_name = 'Peng-Robinson Methane'
  else 
    eos_density_name = 'Unknown'
  endif

  ! energy
  if (associated(EOSGasEnergyPtr,EOSGasEnergyConstant)) then
    eos_energy_name = 'Constant'
  else if (associated(EOSGasEnergyPtr,EOSGasEnergyIdeal)) then
    eos_energy_name = 'Ideal Gas Law'
  else if (associated(EOSGasEnergyPtr,EOSGasEnergyIdealMethane)) then
    eos_energy_name = 'Ideal Gas Law - Methane'
  else
    eos_energy_name = 'Unknown'
  endif

  ! viscosity
  if (associated(EOSGasViscosityPtr,EOSGasViscosityConstant)) then
    eos_viscosity_name = 'Constant'
  else if (associated(EOSGasViscosityPtr,EOSGasViscosity1)) then
    eos_viscosity_name = 'Default'
  else
    eos_viscosity_name = 'Unknown'
  endif

  ! Henry's constant
  if (associated(EOSGasHenryPtr,EOSGasHenryConstant)) then
    eos_henry_name = 'Constant'
  else if (associated(EOSGasHenryPtr,EOSGasHenry_air_noderiv)) then
    eos_henry_name = 'Default'
  else
    eos_henry_name = 'Unknown'
  endif

  ! saturation pressure
  eos_saturation_pressure_name = 'Unknown'

  do itemp = 1, ntemp
    ! Need saturation pressure to calculate air partial pressure
    call EOSWaterSaturationPressure(temp(itemp),saturation_pressure,ierr)
    saturation_pressure_array(itemp) = saturation_pressure
    call EOSGasHenryPtr(temp(itemp),saturation_pressure,henry(itemp))
    do ipres = 1, npres
      call EOSGasDensityPtr(temp(itemp),pres(ipres), &
                            density_kg(ipres,itemp), &
                            dum1,dum2,ierr)
      call EOSGasEnergyPtr(temp(itemp),pres(ipres),&
                           enthalpy(ipres,itemp),dum1,dum2, &
                           internal_energy(ipres,itemp),dum3,dum4,ierr)
      air_pressure = pres(ipres) - saturation_pressure
      if (air_pressure > 0.d0) then
        call EOSGasViscosityPtr(temp(itemp),air_pressure,pres(ipres), &
                                density_kg(ipres,itemp), &
                                viscosity(ipres,itemp),ierr)
      else
        viscosity(ipres,itemp) = NaN
      endif
    enddo
  enddo

100 format(100es16.8)
  if (len_trim(filename) == 0) filename = 'eos_gas_test.txt'
  open(unit=IUNIT_TEMP,file=filename)
  header = 'T[C], P[Pa], &
    &Density (' // trim(eos_density_name) // ') [kg/m^3], &
    &Enthalpy (' // trim(eos_energy_name) // ') [J/kmol], &
    &Internal Energy (' // trim(eos_energy_name) // ') [J/kmol], &
    &Viscosity (' // trim(eos_viscosity_name) // ') [Pa-s], &
    &Saturation Pressure (' // trim(eos_saturation_pressure_name) // ") [Pa], &
    &Henry's Constant (" // trim(eos_henry_name) // ') [-]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp, npres
  do itemp = 1, ntemp
    do ipres = 1, npres
      write(IUNIT_TEMP,100) temp(itemp), pres(ipres), &
            density_kg(ipres,itemp), enthalpy(ipres,itemp), &
            internal_energy(ipres,itemp), viscosity(ipres,itemp), &
            saturation_pressure_array(itemp), henry(itemp)
    enddo
  enddo
  close(IUNIT_TEMP)

  deallocate(temp)
  deallocate(pres)
  deallocate(density_kg)
  deallocate(enthalpy)
  deallocate(internal_energy)
  deallocate(viscosity)
  deallocate(saturation_pressure_array)
  deallocate(henry)

end subroutine EOSGasTest

! ************************************************************************** !

end module EOS_Gas_module
