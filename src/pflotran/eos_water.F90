module EOS_Water_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "finclude/petscsys.h"

  ! module variables
  PetscReal :: constant_density
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity

  ! In order to support generic EOS subroutines, we need the following:
  ! 1. An interface declaration that defines the argument list (best to have 
  !    "Dummy" appended.
  ! 2. A procedure pointer that is initially set to null.  This pointer is
  !    pointed to the appropriate subroutine later on (e.g. EOSWaterInit())
  ! 3. An interface for derivative/non-derivative versions

  ! procedure pointer declarations
  procedure(EOSWaterViscosityDummy), pointer :: EOSWaterViscosityPtr => null()
  procedure(EOSWaterSatPressDummy), pointer :: &
    EOSWaterSaturationPressurePtr => null()
  procedure(EOSWaterDensityEnthalpyDummy), pointer :: &
    EOSWaterDensityEnthalpyPtr => null()
  procedure(EOSWaterSteamDenEnthDummy), pointer :: &
    EOSWaterSteamDensityEnthalpyPtr => null()
  
  ! interface blocks
  interface
    subroutine EOSWaterViscosityDummy(T, P, PS, dPS_dT, VW, &
                                      calculate_derivatives, &
                                      dWV_dT, dWV_dP, dWV_dPS, ierr)
      implicit none
      PetscReal, intent(in) :: T, P, PS, dPS_dT
      PetscReal, intent(out) :: VW
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dWV_dT, dWV_dP, dWV_dPS
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterViscosityDummy
    subroutine EOSWaterSatPressDummy(T, PS, calculate_derivatives, &
                                               dPS_dT, ierr)
      implicit none
      PetscReal, intent(in) :: T
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: PS, dPS_dT
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterSatPressDummy
    subroutine EOSWaterDensityEnthalpyDummy(t,p,dw,dwmol,hw, &
                                            calculate_derivatives, &
                                            dwp,dwt,hwp,hwt,scale,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw,dwmol,dwp,dwt
      PetscReal, intent(out) :: hw,hwp,hwt
      PetscReal, intent(in) :: scale
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterDensityEnthalpyDummy
    subroutine EOSWaterSteamDenEnthDummy(t,p,pa,dg,dgmol,hg, &
                                         calculate_derivatives, &
                                         dgp,dgt,hgp,hgt,scale,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p,pa
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dg,dgmol,dgp,dgt
      PetscReal, intent(out) :: hg,hgp,hgt
      PetscReal, intent(in) :: scale
      PetscErrorCode, intent(out) :: ierr 
    end subroutine EOSWaterSteamDenEnthDummy
  end interface
  
  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  interface EOSWaterViscosity
    procedure EOSWaterViscosityNoDerive
    procedure EOSWaterViscosityDerive
  end interface
  interface EOSWaterSaturationPressure
    procedure EOSWaterSatPresNoDerive
    procedure EOSWaterSatPresDerive
  end interface
  interface EOSWaterDensityEnthalpy
    procedure EOSWaterDenEnthNoDerive
    procedure EOSWaterDenEnthDerive
  end interface
  interface EOSWaterSteamDensityEnthalpy
    procedure EOSWaterSteamDenEnthNoDerive
    procedure EOSWaterSteamDenEnthDerive
  end interface

  ! the "public" definition that makes subroutines visible outside.
  public :: EOSWaterInit, &
            EOSWaterVerify, &
            EOSWaterViscosity, &
            EOSWaterSaturationPressure, &
            EOSWaterDensityEnthalpy, &
            EOSWaterSteamDensityEnthalpy, &
            EOSWaterDuanMixture, &
            EOSWaterDensityNaCl, &
            EOSWaterViscosityNaCl, &
            EOSWaterdensity, &
            EOSWaterInternalEnergyIce, &
            EOSWaterDensityIce, &
            EOSWaterSaturationTemperature
  
  public :: EOSWaterSetDensityIFC67, &
            EOSWaterSetEnthalpyIFC67, &
            EOSWaterSetDensityConstant, &
            EOSWaterSetEnthalpyConstant, &
            EOSWaterSetViscosityConstant
 
  contains
  
! ************************************************************************** !
subroutine EOSWaterInit()

  implicit none
  
  constant_density = -999.d0
  constant_viscosity = -999.d0
  constant_enthalpy = -999.d0
  
  EOSWaterDensityEnthalpyPtr => EOSWaterDensityEnthalpyIFC67
  EOSWaterViscosityPtr => EOSWaterViscosity1
  EOSWaterSaturationPressurePtr => EOSWaterSaturationPressure1
  EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpy1
  
end subroutine EOSWaterInit

! ************************************************************************** !
subroutine EOSWaterVerify(ierr)

  implicit none
  
  PetscErrorCode :: ierr
  
  ierr = 0
  
  if ((associated(EOSWaterDensityEnthalpyPtr, &
                  EOSWaterDensityEnthalpyConstant) .and. &
       (constant_density < -998.d0 .or. constant_enthalpy < -998.d0)) .or. &
      (associated(EOSWaterDensityEnthalpyPtr, &
                  EOSWaterDensityEnthalpyIFC67) .and. &
       (constant_density > -998.d0 .or. constant_enthalpy > -998.d0))) then
    ierr = 1
  endif
  
  if ((associated(EOSWaterViscosityPtr, &
                  EOSWaterViscosityConstant) .and. &
       constant_viscosity < -998.d0) .or. &
      (associated(EOSWaterViscosityPtr, &
                  EOSWaterViscosity1) .and. &
       constant_viscosity > -998.d0)) then
    ierr = 1
  endif
  
end subroutine EOSWaterVerify

! ************************************************************************** !
subroutine EOSWaterSetDensityIFC67()

  implicit none
  
  EOSWaterDensityEnthalpyPtr => EOSWaterDensityEnthalpyIFC67
  
end subroutine EOSWaterSetDensityIFC67

! ************************************************************************** !
subroutine EOSWaterSetEnthalpyIFC67()

  implicit none
  
  EOSWaterDensityEnthalpyPtr => EOSWaterDensityEnthalpyIFC67
  
end subroutine EOSWaterSetEnthalpyIFC67

! ************************************************************************** !
subroutine EOSWaterSetDensityConstant(density)

  implicit none
  
  PetscReal :: density
  
  constant_density = density  
  EOSWaterDensityEnthalpyPtr => EOSWaterDensityEnthalpyConstant
  
end subroutine EOSWaterSetDensityConstant

! ************************************************************************** !
subroutine EOSWaterSetEnthalpyConstant(enthalpy)

  implicit none
  
  PetscReal :: enthalpy
  
  constant_enthalpy = enthalpy  
  EOSWaterDensityEnthalpyPtr => EOSWaterDensityEnthalpyConstant
  
end subroutine EOSWaterSetEnthalpyConstant

! ************************************************************************** !
subroutine EOSWaterSetViscosityConstant(viscosity)

  implicit none
  
  PetscReal :: viscosity
  
  constant_viscosity = viscosity  
  EOSWaterViscosityPtr => EOSWaterViscosityConstant
  
end subroutine EOSWaterSetViscosityConstant

! ************************************************************************** !
subroutine EOSWaterViscosityNoDerive(T, P, PS, VW, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(out) :: VW ! water viscosity
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dPS_dT ! derivative of PS with respect to temp
  PetscReal :: dum1, dum2, dum3
  
  dPS_dT = 0.d0
  call EOSWaterViscosityPtr(T, P, PS, dPS_dT, VW, PETSC_FALSE, dum1, dum2, &
                            dum3, ierr)
  
end subroutine EOSWaterViscosityNoDerive

! ************************************************************************** !
subroutine EOSWaterViscosityDerive(T, P, PS, dPS_dT, VW, dWV_dT, dWV_dP, &
                                   dWV_dPS, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dWV_dT, dWV_dP, dWV_dPS ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  call EOSWaterViscosityPtr(T, P, PS, dPS_dT, VW, PETSC_TRUE, &
                            dWV_dT, dWV_dP, dWV_dPS, ierr)
  
end subroutine EOSWaterViscosityDerive

! ************************************************************************** !
subroutine EOSWaterViscosity1(T, P, PS, dPS_dT, VW, calculate_derivatives, &
                              dWV_dT, dWV_dP, dWV_dPS, ierr)

! Calculates the viscosity of water and derivatives as a function of 
! temperature, pressure, and saturation pressure.

  implicit none
    
  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(out) :: VW ! water viscosity
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dWV_dT, dWV_dP, dWV_dPS ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: EX, PHI, AM, pwr, aln10
  
  EX  = 247.8d0/(T+133.15d0)
  PHI = 1.0467d0*(T-31.85d0)
  AM  = 1.d0+PHI*(P-PS)*1.d-11
  pwr = 10.d0**EX
  VW = 1.d-7*AM*241.4d0*pwr
    
  if (calculate_derivatives) then
    aln10 = log(10.d0)
    dWV_dT = VW/AM*1.d-11* &
            ! dAM_PHI_dT       dAM_PS_dT
            (1.0467d0*(P-PS) - PHI*dPS_dT) - &
            ! dpwr_EX_dT
            VW*aln10*247.8d0/(T+133.15d0)**2
    dWV_dP = VW/AM*PHI*1.d-11
    dWV_dPS = -1.d0 * dWV_dP
  else
    dWV_dT = -999.d0
    dWV_dP = -999.d0
    dWV_dPS = -999.d0
  endif
  ierr = 0
 
end subroutine EOSWaterViscosity1

! ************************************************************************** !
subroutine EOSWaterViscosityConstant(T, P, PS, dPS_dT, VW, &
                                     calculate_derivatives, &
                                     dWV_dT, dWV_dP, dWV_dPS, ierr)

! Calculates the viscosity of water and derivatives as a function of 
! temperature, pressure, and saturation pressure.

  implicit none
    
  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(out) :: VW ! water viscosity
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dWV_dT, dWV_dP, dWV_dPS ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  VW = constant_viscosity
  
  dWV_dT = 0.d0
  dWV_dP = 0.d0
  dWV_dPS = 0.d0
  
end subroutine EOSWaterViscosityConstant

! ************************************************************************** !
subroutine EOSWaterSatPresNoDerive(T, PS, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(out) :: PS ! Saturation pres. and derivative
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dummy
  
  call EOSWaterSaturationPressurePtr(T, PS, PETSC_FALSE, dummy, ierr)
  
end subroutine EOSWaterSatPresNoDerive

! ************************************************************************** !
subroutine EOSWaterSatPresDerive(T, PS, dPS_dT, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(out) :: ierr
  
  call EOSWaterSaturationPressurePtr(T, PS, PETSC_TRUE, dPS_dT, ierr)
  
end subroutine EOSWaterSatPresDerive

! ************************************************************************** !
subroutine EOSWaterSaturationPressure1(T, PS, calculate_derivatives, dPS_dT, &
                                       ierr)

    implicit none

    PetscReal, intent(in) :: T ! temperature
    PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
    PetscBool, intent(in) :: calculate_derivatives
    PetscErrorCode, intent(out) :: ierr
  
    PetscReal, save, dimension(9) :: A(9)
    PetscReal :: TC, SC, PCAP, E1, E2
    PetscReal :: one_m_tc, one_m_tc_sq, E2_bottom
    PetscReal :: dTC_dT, dSC_dTC, dE1_dTC, dE2_dTC, dPC_dSC, dPC_dTC
    
    DATA A/ &
      -7.691234564d0,-2.608023696d1,-1.681706546d2,6.423285504d1, &
      -1.189646225d2,4.167117320d0,2.097506760E1,1.d9,6.d0/
   
    if (T .GT. 500.d0) then
      ierr = 1
      return
    end if
    TC = (T+273.15d0)/647.3d0
    one_m_tc = 1.d0-TC
    one_m_tc_sq = one_m_tc*one_m_tc
    SC = A(1)*one_m_tc+A(2)*one_m_tc_sq+A(3)*one_m_tc**3.d0+ &
         A(4)*one_m_tc**4.d0+A(5)*one_m_tc**5.d0
    E1 = TC*(1.d0+A(6)*one_m_tc+A(7)*one_m_tc_sq)
    E2_bottom = A(8)*one_m_tc_sq+A(9)
    E2 = one_m_tc/E2_bottom
    PCAP = EXP(SC/E1-E2)
   
    PS = PCAP*2.212d7
    
    if (calculate_derivatives) then
      dTC_dT = 1.d0/647.3d0
      dSC_dTC = -A(1)-2.d0*A(2)*one_m_tc-3.d0*A(3)*one_m_tc_sq- &
                4.d0*A(4)*one_m_tc**3.-5.d0*A(5)*one_m_tc**4.
      dE1_dTC = (1.d0+A(6)*one_m_tc+A(7)*one_m_tc_sq)+ &
                TC*(-A(6)-2.d0*A(7)*one_m_tc)
      dE2_dTC = -1.d0/E2_bottom+one_m_tc/(E2_bottom*E2_bottom)*2.d0*one_m_tc
      dPC_dTC = (-SC/(E1*E1)*dE1_dTC-dE2_dTC)*PCAP
      dPC_dSC = 1.d0/E1*PCAP
      dPS_dT = (dPC_dSC*dSC_dTC+dPC_dTC)*dTC_dT*2.212d7
    else
      dPS_dT = -999.d0
    endif
    ierr = 0

end subroutine EOSWaterSaturationPressure1  

! ************************************************************************** !
subroutine EOSWaterDenEnthNoDerive(t,p,dw,dwmol,hw,scale,ierr)

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscReal, intent(out) :: dw,dwmol
  PetscReal, intent(out) :: hw
  PetscReal, intent(in) :: scale
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2, dum3, dum4
  
  call EOSWaterDensityEnthalpyPtr(t,p,dw,dwmol,hw,PETSC_FALSE, &
                                  dum1,dum2,dum3,dum4,scale,ierr)
  
end subroutine EOSWaterDenEnthNoDerive

! ************************************************************************** !
subroutine EOSWaterDenEnthDerive(t,p,dw,dwmol,hw, &
                                 dwp,dwt,hwp,hwt,scale,ierr)
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscReal, intent(in) :: scale
  PetscErrorCode, intent(out) :: ierr
  
  call EOSWaterDensityEnthalpyPtr(t,p,dw,dwmol,hw,PETSC_TRUE, &
                                  dwp,dwt,hwp,hwt,scale,ierr)
  
end subroutine EOSWaterDenEnthDerive

! ************************************************************************** !
subroutine EOSWaterDensityEnthalpyIFC67(t,p,dw,dwmol,hw, &
                                        calculate_derivatives, &
                                        dwp,dwt,hwp,hwt,scale,ierr)

!  This subroutine calculates water and steam-gas mixture properties.
!  The water and steam properties are valid in the range of:
!
!            0 < p < 165.4 * 10^5 pascals (165.4 bars)
!            0 < t < 800 centigrade (1073.15 Kelvin)
!
!  The properties cover densities, enthalpies, internal energies,
!  and partial derivatives of these quanties with respect to
!  pressure and temperature.
!
!  For saturated fluid, it will also calculate water saturation
!  temperature for the specified pressure using Newton-Raphson and
!  the derivative dts/dp (=tsp) or Ps for a given temperature.
!
!  Ref.: International Formulation Committee of the Sixth International
!       Conference on Properties of Steam (1967).

  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscReal, intent(in) :: scale
  PetscErrorCode, intent(out) :: ierr
  
  PetscInt :: i
    
  PetscReal, save :: aa(0:22)
  PetscReal, save :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
  
  PetscReal :: beta,beta2x,beta4,theta,utheta,theta2x,theta18,theta20
  PetscReal :: xx,yy,zz
  PetscReal :: u0,u1,u2,u3,u4,u5,u6,u7,u8,u9
  PetscReal :: tempreal
  PetscReal :: v0_1, v1_1, v2_1, v3_1, v4_1
  PetscReal :: v1_2, v2_2, v3_2, v4_2, v20_2, v40_2
  PetscReal :: v1_3, v2_3, v3_3, v4_3
  PetscReal :: v1_4, v2_4, v3_4
  PetscReal :: v1_5, v2_5
  PetscReal :: v1_6
  PetscReal :: term1,term2,term2t,term3,term3t,term3p,term4,term4t,term4p, &
               term5,term5t,term5p,term6,term6t,term6p,term7,term7t,term7p
  PetscReal :: dv2t,dv2p,dv3t
  PetscReal :: vr,ypt,yptt,zpt,zpp,vrpt,vrpp,cnv
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol,vc1molh
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: four = 4.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: seven = 7.d0
  PetscReal, parameter :: eight = 8.d0
  PetscReal, parameter :: nine = 9.d0
  PetscReal, parameter :: ten = 10.d0
  
!   save aa,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
    
!   save
  
  data aa/ &
!-----data aa0,aa1,aa2,aa3/
        6.824687741d03,-5.422063673d02,-2.096666205d04, 3.941286787d04, &
!-----data aa4,aa5,aa6,aa7/
        -6.733277739d04, 9.902381028d04,-1.093911774d05, 8.590841667d04, &
!-----data aa8,aa9,aa10,aa11/
        -4.511168742d04, 1.418138926d04,-2.017271113d03, 7.982692717d00, &
!-----data aa12,aa13,aa14,aa15/
        -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, 2.421647003d02, &
!-----data aa16,aa17,aa18,aa19/
        1.269716088d-10,2.074838328d-7, 2.174020350d-8, 1.105710498d-9, &
!-----data aa20,aa21,aa22/    
        1.293441934d01, 1.308119072d-5, 6.047626338d-14/

  data a1,a2,a3,a4/ &
  8.438375405d-1, 5.362162162d-4, 1.720000000d00, 7.342278489d-2/
  data a5,a6,a7,a8/ &
  4.975858870d-2, 6.537154300d-1, 1.150000000d-6, 1.510800000d-5/
  data a9,a10,a11,a12/ &
  1.418800000d-1, 7.002753165d00, 2.995284926d-4, 2.040000000d-1/
    
  ierr = 0

  tc1 = 647.3d0    ! K
  pc1 = 2.212d7    ! Pa 
  vc1 = 0.00317d0  ! m^3/kg
  utc1 = one/tc1   ! 1/C
  upc1 = one/pc1   ! 1/Pa
  vc1mol = vc1*FMWH2O ! m^3/kmol
    
  theta = (t+273.15d0)*utc1
  theta2x = theta*theta
  theta18 = theta**18
  theta20 = theta18*theta2x
    
  beta = p*upc1
  beta2x = beta*beta
  beta4  = beta2x*beta2x
    
  yy = one-a1*theta2x-a2*theta**(-6)
  xx = a3*yy*yy-two*(a4*theta-a5*beta)
    
!   Note: xx may become negative near the critical point-pcl.
  if (xx.gt.zero) then
    xx = sqrt(xx)
  else
    write(*,*) 'Warning: negative term in density (wateos): ','t= ',t,' p= ',p,' xx= ',xx
    ierr = 1
    xx = 1.e-6               !set arbitrarily
  end if
  zz = yy + xx                                     
  u0 = -five/17.d0
  u1 = aa(11)*a5*zz**u0
  u2 = one/(a8+theta**11)
  u3 = aa(17)+(two*aa(18)+three*aa(19)*beta)*beta
  u4 = one/(a7+theta18*theta)
  u5 = (a10+beta)**(-4)
  u6 = a11-three*u5
  u7 = aa(20)*theta18*(a9+theta2x)
  u8 = aa(15)*(a6-theta)**9
    
  vr = u1+aa(12)+theta*(aa(13)+aa(14)*theta)+u8*(a6-theta) &
        +aa(16)*u4-u2*u3-u6*u7+(three*aa(21)*(a12-theta) &
        +four*aa(22)*beta/theta20)*beta2x
    
  dwmol = one/(vr*vc1mol) ! kmol/m^3
  dw = one/(vr*vc1) ! kg/m^3

  ! ypt used for enthalpy even if derivative not calculated
  ypt = six*a2*theta**(-7)-two*a1*theta
  
  !---calculate derivatives for water density
  if (calculate_derivatives) then
    zpt = ypt+(a3*yy*ypt-a4)/xx
    zpp = a5/xx
    u9 = u0*u1/zz
    vrpt = u9*zpt+aa(13)+two*aa(14)*theta-ten*u8 &
        -19.d0*aa(16)*u4*u4*theta18+11.d0*u2*u2*u3*theta**10 &
        -aa(20)*u6*(18.d0*a9*theta18+20.d0*theta20)/theta &
        -(three*aa(21)+80.d0*aa(22)*beta/(theta20*theta))*beta2x
    
    vrpp = u9*zpp-u2*(two*aa(18)+six*aa(19)*beta)-12.d0*u7*u5/ &
        (a10+beta)+(six*aa(21)*(a12-theta)+12.d0*aa(22)*beta/ &
        theta20)*beta
    
    cnv = -one/(vc1mol*vr*vr)
    dwt = cnv*vrpt*utc1 ! kmol/m^3/C
    dwp = cnv*vrpp*upc1 ! kmol/m^3/Pa
  else
    dwt = -999.d0
    dwp = -999.d0
  endif

!   print *,'water_eos: ',p,t,dwp,cnv,vrpp,upc1
    
!---compute enthalpy internal energy and derivatives for water
  utheta = one/theta
  term1 = aa(0)*theta
  term2 = -aa(1)
  ! term2t is part of the derivative calc., but left here to avoid
  ! recomputing the expensive do loop
  term2t = zero
  do i = 3,10
    tempreal = dfloat(i-2)*aa(i)*theta**(i-1)
    term2t = term2t+tempreal*utheta*dfloat(i-1)
    term2 = term2+tempreal                            
  end do
    
  ! "v" section 1
  v0_1 = u1/a5
  v2_1 = 17.d0*(zz/29.d0-yy/12.d0)+five*theta*ypt/12.d0
  v3_1 = a4*theta-(a3-one)*theta*yy*ypt
  v1_1 = zz*v2_1+v3_1
  term3 = v0_1*v1_1
  
  ! block 1 removed from here

  ! "v" section 2
  v1_2 = nine*theta+a6
  v20_2 = (a6-theta)
  v2_2 = v20_2**9
  v3_2 = a7+20.d0*theta**19
  v40_2 = a7+theta**19
  v4_2 = one/(v40_2*v40_2)
  ! term4p is a derivative, but left due to dependency in term4
  term4p = aa(12)-aa(14)*theta2x+aa(15)*v1_2*v2_2+aa(16)*v3_2*v4_2
  term4 = term4p*beta
  
  ! block 2 removed from here
    
  ! "v" section 3
  v1_3 = beta*(aa(17)+aa(18)*beta+aa(19)*beta2x)
  v2_3 = 12.d0*theta**11+a8
  v4_3 = one/(a8+theta**11)
  v3_3 = v4_3*v4_3
  term5 = v1_3*v2_3*v3_3
  
  ! block 3 removed from here
    
  ! "v" section 4
  v1_4 = (a10+beta)**(-3)+a11*beta
  v3_4 = (17.d0*a9+19.d0*theta2x)
  v2_4 = aa(20)*theta18*v3_4       
  term6 = v1_4*v2_4

  ! block 4 removed from here
    
  ! "v" section 5
  v1_5 = 21.d0*aa(22)/theta20*beta4
  v2_5 = aa(21)*a12*beta2x*beta
  term7 = v1_5+v2_5  
    
  ! The scale (usually 1e-6) should be applied below to pc1 to convert
  ! from Pa -> MPa.  Applying it to critical volume does the same, but is
  ! misleading - geh
  vc1molh = vc1mol*scale
    
  ! "v" section 6
  v1_6 = pc1*vc1molh ! MPa * m^3/kmol = MJ/kmol assuming scale = 1e-6
  hw = (term1-term2+term3+term4-term5+term6+term7)*v1_6
    
  if (calculate_derivatives) then

    ! block 1
    yptt = -two*a1-42.d0*a2/theta**8
    dv2t = 17.d0*(zpt/29.d0-ypt/12.d0)+five/12.d0*(ypt+theta*yptt) 
    dv3t = a4-(a3-one)*(theta*yy*yptt+yy*ypt+theta*ypt*ypt)
    dv2p = 17.d0*zpp/29.d0
    v4_1 = five*v1_1/(17.d0*zz)       
    term3t = v0_1*(zz*dv2t+(v2_1-v4_1)*zpt+dv3t)
    term3p = v0_1*(zz*dv2p+(v2_1-v4_1)*zpp)
  
    ! block 2
    term4t = (-two*aa(14)*theta+nine*aa(15)*(v2_2-v1_2*v2_2/v20_2) &
             +38.d0*theta18*aa(16)*(ten*v4_2-v3_2*v4_2/v40_2))*beta

    ! block 3
    term5p = v3_3*v2_3*(aa(17)+two*aa(18)*beta+three*aa(19)*beta2x)
    term5t = v1_3*(132.d0*v3_3*theta**10-22.d0*v2_3*v3_3*v4_3*theta**10)
  
    ! block 4
    term6p = v2_4*(a11-three*(a10+beta)**(-4))
    term6t = v1_4*aa(20)*theta18*(18.d0*v3_4*utheta+38.d0*theta)
    
    ! block 5
    term7p = beta2x*(three*aa(21)*a12+84.d0*aa(22)*beta/theta20)
    term7t = -420.d0*aa(22)*beta4/(theta20*theta)

    hwp = (term3p+term4p-term5p+term6p+term7p)*vc1molh
    hwt = (aa(0)-term2t+term3t+term4t-term5t+term6t+term7t)*v1_6*utc1
  else
    hwp = -999.d0
    hwt = -999.d0
  endif
    
end subroutine EOSWaterDensityEnthalpyIFC67

! ************************************************************************** !
subroutine EOSWaterDensityEnthalpyConstant(t,p,dw,dwmol,hw, &
                                           calculate_derivatives, &
                                           dwp,dwt,hwp,hwt,scale,ierr)
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscReal, intent(in) :: scale
  PetscErrorCode, intent(out) :: ierr
  
  dw = constant_density ! kg/m^3
  dwmol = dw*FMWH2O ! kmol/m^3
  
  dwp = 0.d0
  dwt = 0.d0

  hw = constant_enthalpy ! MJ/kmol
  
  hwp = 0.d0
  hwt = 0.d0
  
end subroutine EOSWaterDensityEnthalpyConstant

! ************************************************************************** !
subroutine EOSWaterSteamDenEnthNoDerive(t,p,pa,dg,dgmol,hg,scale,ierr)

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: p,pa   ! Pressure in Pascals.
  PetscReal, intent(out) :: dg,dgmol
  PetscReal, intent(out) :: hg
  PetscReal, intent(in) :: scale
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2, dum3, dum4
  
  call EOSWaterSteamDensityEnthalpyPtr(t,p,pa,dg,dgmol,hg,PETSC_FALSE, &
                                       dum1,dum2,dum3,dum4,scale,ierr)
  
end subroutine EOSWaterSteamDenEnthNoDerive

! ************************************************************************** !
subroutine EOSWaterSteamDenEnthDerive(t,p,pa,dg,dgmol,hg, &
                                      dgp,dgt,hgp,hgt,scale,ierr)
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: p,pa   ! Pressure in Pascals.
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscReal, intent(in) :: scale
  PetscErrorCode, intent(out) :: ierr
  
  call EOSWaterSteamDensityEnthalpyPtr(t,p,pa,dg,dgmol,hg,PETSC_TRUE, &
                                       dgp,dgt,hgp,hgt,scale,ierr)
  
end subroutine EOSWaterSteamDenEnthDerive

! ************************************************************************** !
subroutine EOSWaterSteamDensityEnthalpy1(t,p,pa,dg,dgmol,hg, &
                                         calculate_derivatives, &
                                         dgp,dgt,hgp,hgt,scale,ierr)
! t/C  p/Pa dgmol/(mol/m^3)  h/MJ/mol
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: p,pa   ! Pressure in Pascals.
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscReal, intent(in) :: scale
  PetscErrorCode, intent(out) :: ierr
  
  PetscInt, save :: n(8),ll(8),x(8,2),z(8,3)
  PetscReal, save :: xi1,xl0,xl1,xl2
  PetscReal, save :: b(8,2),bb(0:9,0:6)
  PetscReal :: sumbx(8),sumbxt(8)
  
  PetscInt :: i,j
  PetscReal, save :: delp,delt
  PetscReal :: beta,betap,betal,betalp,betalp1,betal1,ubeta,theta, &
            thetap,utheta
  PetscReal :: xx,xxp
  PetscReal :: f1,f2,fim1,fim2,sum,sumt,sum1,sum1t,sum1p, &
            sum2,sum2t
  PetscReal :: u1,u1t,u1p,u2,u2p,u2t,u3,u3p,u3t,v1,v1t
  PetscReal :: term,term1,term1t,term1p,term2,term2t,term2p, &
            term3,term3t,term3p,term4,term4t,term4p,term5,term5t,term5p
  PetscReal :: hr,hrp,hrt,hrpt,hrpp
  PetscReal :: vr,vrpt,vrpp
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol,vc1molh
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: four = 4.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: seven = 7.d0
  PetscReal, parameter :: eight = 8.d0
  PetscReal, parameter :: nine = 9.d0
  PetscReal, parameter :: ten = 10.d0
 
  data delt,delp/1.d-6,1.d-6/
    
  data n/2,3,2,2,3,2,2,2/
  data ll/0,0,0,0,0,1,1,2/
  data x/0,0,0,0,0,14,19,54, &
          0,0,0,0,0, 0, 0,27/
  data z/13,18,18,25,32,12,24,24, &
          3, 2,10,14,28,11,18,14, &
          0, 1, 0, 0,24, 0, 0, 0/
    
  data b/7.6333333333d-1,0.d0,0.d0,0.d0,0.d0,4.006073948d-1, &
          8.636081627d-2,-8.532322921d-1, &
          0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,3.460208861d-1/

!---sub-region 2
  data bb/ &
!---bb(0-9,0)
          1.683599274d1,  8*0.d0        , 1.936587558d2, &
!---bb(0-9,1)
          2.856067796d01, 6.670375918d-2, 8.390104328d-2, &
          4.520918904d-1,-5.975336707d-1, 5.958051609d-1, &
          1.190610271d-1, 1.683998803d-1, 6.552390126d-3, &
          -1.388522425d+3, &
!---bb(0-9,2)
          -5.438923329d01, 1.388983801d0 , 2.614670893d-2, &
          1.069036614d-1,-8.847535804d-2,-5.159303373d-1, &
          -9.867174132d-2,-5.809438001d-2, 5.710218649d-4, &
          4.126607219d03, &
!---bb(0-9,3)
          4.330662834d-1, 0.d0          ,-3.373439453d-2, &
          0.d0          , 0.d0          , 2.075021122d-1, &
          0.d0          , 0.d0          , 0.d0          , &
          -6.508211677d03, &
!---bb(0-9,4)
          -6.547711697d-1, 8*0.d0        , 5.745984054d03, &
!---bb(0-9,5)
          8.565182058d-2, 8*0.d0        ,-2.693088365d03, &
!---bb(0-9,6)
                          9*0.d0        , 5.235718623d02/
    
  data xi1/4.260321148d0/
  data xl0,xl1,xl2/15.74373327d0,-34.17061978d0,19.31380707d0/
    
  ierr = 0
  
  tc1 = 647.3d0    
  pc1 = 2.212d7     
  vc1 = 0.00317d0
  utc1 = one/tc1
  upc1 = one/pc1
  vc1mol = vc1*18.0153d0
  vc1molh = vc1mol*scale   

  theta  = (t+273.15d0)*utc1
  beta   = (p-pa)*upc1
  ubeta  = one/beta
  utheta = one/theta
  xx = exp(b(1,1)*(one-theta))
    
!---compute steam density and derivatives
    
  term1 = xi1*theta*ubeta
  term1t = xi1*ubeta
  term1p = -term1*ubeta        
    
  do i = 1,8
    sum = zero
    sumt = zero
    do j = 1,n(i)
      u1 = bb(i,j)*xx**z(i,j)     
      sumt = sumt-b(1,1)*z(i,j)*u1                     
      sum = sum + u1                     
    end do
    sumbx(i) = sum
    sumbxt(i) = sumt
  end do
    
  term2  = sumbx(1)+beta*(two*sumbx(2)+beta*(three*sumbx(3) &
            +beta*(four*sumbx(4)+five*beta*sumbx(5))))
  term2t = sumbxt(1)+beta*(two*sumbxt(2)+beta*(three*sumbxt(3) &
            +beta*(four*sumbxt(4)+five*beta*sumbxt(5))))
  term2p = two*sumbx(2)+beta*(six*sumbx(3)+beta*(12.d0*sumbx(4) &
            +20.d0*beta*sumbx(5)))
    
  term3  = zero
  term3p = zero
  term3t = zero
    
  do i = 6,8
    fim1 = dfloat(i-1)
    fim2 = fim1-one    
    
    sum2 = zero
    sum2t = zero
    do j = 1,ll(i)
      u1 = b(i,j)*xx**x(i,j)
      sum2t = sum2t-b(1,1)*x(i,j)*u1
      sum2 = sum2 + u1
    end do
    
    f1 = fim2*beta**(1-i)
    f2 = beta**(2-i)
    u1 = one/(f2+sum2)**2
    term = u1*f1*sumbx(i)
    term3 = term3+term      
    u2 = two*term/(f2+sum2)
    term3t = term3t+u1*(f1*sumbxt(i)+u2*sum2t)    
    term3p = term3p-u1*(f1*fim1*sumbx(i)+u2*fim2*f2)*ubeta
  end do
    
  term4 = bb(9,0)
  term4t = zero  
  do i = 1,6
    u1 = bb(9,i)*xx**i 
    term4t = term4t-b(1,1)*dfloat(i)*u1
    term4 = term4 + u1             
  end do
    
  betal = xl0+theta*(xl1+xl2*theta)
  betalp = xl1+ two*xl2*theta
  u1 = 11.d0*(beta/betal)**10
  term4t = u1*(term4t-10.d0*betalp*term4/betal)
  term4 = term4*u1                  
  term4p = 10.d0*term4*ubeta             
    
  vr = term1-term2-term3+term4
  vrpt = term1t-term2t-term3t+term4t
  vrpp = term1p-term2p-term3p+term4p
    
  u1 = -one/(vc1mol*vr)
  dgmol = -u1
  dg = one/(vc1*vr)
  
  if (calculate_derivatives) then
    u1 = u1/vr
    dgt = u1*vrpt*utc1
    dgp = u1*vrpp*upc1
  else
    dgt = -999.d0
    dgp = -999.d0
  endif
    
!---compute steam enthalpy, internal energy, and derivatives
  thetap = theta+delt 
  betap  = beta +delp 
  xxp    = exp(b(1,1)*(one-thetap))
    
  term1  = bb(0,0)*theta
  term1t = bb(0,0)*thetap
    
  term2  = zero
  term2t = zero
  do i = 1,5
    u1 = bb(0,i)*(dfloat(i)-two)
    term2t = term2t+u1*thetap**(i-1)
    term2  = term2+u1*theta**(i-1)
  end do
    
  term3  = zero
  term3t = zero
  term3p = zero
  u1 = one
  u1p = one
  do i = 1,5
    u1 = u1*beta
    u1p = u1p*betap
    
    sum = zero
    sumt = zero
    do j = 1,n(i)
      sumt = sumt+bb(i,j)*(one+z(i,j)*b(1,1)*thetap)*xxp**z(i,j)
      sum  = sum+bb(i,j)*(one+z(i,j)*b(1,1)*theta)*xx**z(i,j)
    end do
    
    term3t = term3t+u1*sumt
    term3p = term3p+u1p*sum
    term3  = term3+u1*sum
  end do
    
  term4  = zero
  term4t = zero
  term4p = zero
  do i = 6,8
    
    sum1  = zero
    sum2  = zero
    sum1t = zero
    sum2t = zero
    
    do j = 1,ll(i)
      u1 = b(i,j)*xxp**x(i,j)
      sum1t = sum1t+x(i,j)*u1               
      sum2t = sum2t+u1                 
      u1 = b(i,j)*xx**x(i,j)
      sum1 = sum1+x(i,j)*u1               
      sum2 = sum2+u1                 
    end do
    
    u1 = one/(beta**(2-i)+sum2)
    u2 = one-b(1,1)*theta*sum1*u1                  
    u3 = b(1,1)*theta
    
    u1t = one/(beta**(2-i)+sum2t)
    u2t = one-b(1,1)*thetap*sum1t*u1t                 
    u3t = b(1,1)*thetap
    
    u1p = one/(betap**(2-i)+sum2)
    u2p = one-b(1,1)*theta*sum1*u1p                   
    u3p = u3          
    
    sum1 = zero
    sum1t = zero
    sum1p = zero
    do j = 1,n(i)
      sum1t = sum1t + bb(i,j)*xxp**z(i,j)*(u2t+z(i,j)*u3t)
      sum1p = sum1p + bb(i,j)*xx **z(i,j)*(u2p+z(i,j)*u3p)
      sum1  = sum1  + bb(i,j)*xx **z(i,j)*(u2 +z(i,j)*u3)
    end do
    
    term4t = term4t+sum1t*u1t
    term4p = term4p+sum1p*u1p
    term4  = term4+sum1*u1
  end do
    
  u1 = ten*betalp/betal
  term5 = (one+theta*u1)*bb(9,0)
    
  betal1  = xl0+thetap*(xl1+xl2*thetap)
  betalp1 = xl1+two*xl2*thetap
  u1t     = ten*betalp1/betal1
  term5t  = (one+thetap*u1t)*bb(9,0)
    
  do i = 1,6
    v1     = one+theta*(u1+dfloat(i)*b(1,1))
    v1t    = one+thetap*(u1t+dfloat(i)*b(1,1))
    term5t = v1t*bb(9,i)*xxp**i+term5t                        
    term5  =  v1*bb(9,i)*xx **i+term5                         
  end do
    
  term5  = term5*beta*(beta/betal)**10
  term5t = term5t*beta*(beta/betal1)**10
  term5p = term5*(betap/beta)**11        
    
  hr   = term1 -term2 -term3 -term4 +term5
  hrt  = term1t-term2t-term3t-term4t+term5t
  hrp  = term1 -term2 -term3p-term4p+term5p
  hrpt = (hrt-hr)/delt
  hrpp = (hrp-hr)/delp
    
  v1 = pc1*vc1molh
  hg = hr*v1
  
  if (calculate_derivatives) then
    hgt = hrpt*v1*utc1
    hgp = hrpp*vc1molh
  else
    hgt = -999.d0
    hgp = -999.d0
  endif

end subroutine EOSWaterSteamDensityEnthalpy1

! ************************************************************************** !
subroutine EOSWaterDuanMixture(t,p,xmol,y_nacl,avgmw,dw_kg,denmix)

! Duan et al. (2008) Energy and Fuels, v 22, 1666-1674.

  implicit none

  PetscReal :: t,tk,p,xco2,xmol,x1,y_nacl,vphi_a1,vphi_a2,vphi,denmix,pw_kg,dw_kg,avgmw

  PetscReal :: fmwh2o = 18.01534d0
  PetscReal :: fmwco2 = 44.0098d0
  PetscReal :: fmwnacl = 58.44277d0

  !duan mixing **************************
  tk = t + 273.15D0; xco2 = xmol;
  call EOSWaterDensityNaCl(t, p*1.D-6, 0.D0, pw_kg)
  pw_kg = pw_kg*1.D3;
  x1 = 1.D0-xco2;
  vphi_a1 = (0.3838402D-3*tk - 0.5595385D0)*tk + 0.30429268D3 + &
            (-0.72044305D5 + 0.63003388D7/tk)/tk;  
  vphi_a2 = (-0.57709332D-5*tk + 0.82764653D-2)*tk - 0.43813556D1 + &
            (0.10144907D4 - 0.86777045D5/tk)/tk;  
  vphi = (1.D0 + vphi_a1 + vphi_a2*p*1.D-6)*(fmwh2o*1.D-3/pw_kg); 
  vphi = x1*((1.D0 - y_nacl)*fmwh2o + y_nacl*fmwnacl)*1.D-3/dw_kg + xco2*vphi;
  denmix = (x1*((1.D0 - y_nacl)*fmwh2o + y_nacl*fmwnacl) + xco2*fmwco2)*1.D-3/vphi;
  denmix = denmix/avgmw
  
end subroutine EOSWaterDuanMixture

! ************************************************************************** !
subroutine EOSWaterDensityNaCl(t,p,xnacl,dnacl)

  ! density: Batzle & Wang (1992)

  implicit none

  PetscReal :: rw0,dnacl,t,p,xnacl
  PetscReal :: rw_mol, hw
  PetscReal :: dum1,dum2,dum3,dum4
  PetscErrorCode  :: ierr
  !units: t [C], p [MPa], xnacl [mass fraction NaCl], dnacl [g/cm^3]

  !rw0 = 1.d0 + 1.d-6*(-80.d0*t - 3.3d0*t**2 + 0.00175d0*t**3 &
  !      + 489.d0*p - 2.d0*t*p + 0.016d0*t**2*p - 1.3d-5*t**3*p &
  !      - 0.333d0*p**2 - 0.002d0*t*p**2)
  call EOSWaterDensityEnthalpyPtr(t,p*1D6,rw0,rw_mol,hw,PETSC_FALSE, &
                                  dum1,dum2,dum3,dum4,1.D-6,ierr)
  rw0=rw0*1.d-3
  dnacl = rw0 + xnacl*(0.668d0 + 0.44d0*xnacl &
          + 1.d-6*(300.d0*p - 2400.d0*p*xnacl + t*(80.d0 &
          + 3.d0*t - 3300.d0*xnacl - 13.d0*p + 47.d0*p*xnacl)))

end subroutine EOSWaterDensityNaCl

! ************************************************************************** !
subroutine EOSWaterViscosityNaCl (t,p,xnacl,visnacl)

  !viscosity: Kestin et al. (1981)

  implicit none

  !units: t [C], p [MPa], xnacl [mass fraction NaCl], visnacl [Pa s]

  PetscReal, save :: a1,a2,a3,b1,b2,b3,c1,c2,c3,c4,wnacl
  PetscReal :: ak,bk,ck
  PetscReal :: beta,betap,betas,betaw
  PetscReal :: t,tt,p,mnacl,xnacl,visnacl,fac,mu0,ms

  data a1,a2,a3 / 3.324d-2, 3.624d-3, -1.879d-4 /
  data b1,b2,b3 / -3.96d-2, 1.02d-2, -7.02d-4 /
  data c1,c2,c3,c4 / 1.2378d0, -1.303d-3, 3.06d-6, 2.55d-8 /

  data wnacl / 58.44277d-3 / ! (kg/mol NaCl)

  !convert pressure to GPa
  p = p*1.d-3

  mnacl = xnacl/(1.d0-xnacl)/wnacl

  tt = 20.d0-t
  ck = (c1 + (c2 + (c3+c4*tt)*tt)*tt)*tt/(96.d0+t)
  ak = (a1 + (a2 + a3*mnacl)*mnacl)*mnacl
  bk = (b1 + (b2 + b3*mnacl)*mnacl)*mnacl

  ms = 6.044d0 + (2.8d-3 + 3.6d-5*t)*t
  fac = mnacl/ms
  betaw = -1.297d0 + (5.74d-2 + (-6.97d-4 + (4.47d-6 - 1.05d-8*t)*t)*t)*t
  betas = 0.545d0 + 2.8d-3 * t - betaw
  betap = (2.5d0 + (-2.d0 + 0.5d0*fac)*fac)*fac
  beta = betas*betap + betaw

  mu0 = 1001.74d-6 * 10.d0**(ak + ck*(bk + 1.d0))

  visnacl = mu0*(1.d0 + beta*p)

end subroutine EOSWaterViscosityNaCl


! ************************************************************************** !
subroutine EOSWaterSaturationTemperature(ts,ps,t_ps,ts_guess,ierr)

  !  This function calculates saturation temperature for a given Ps c  
  !  Ref.: International Formulation Committee of the Sixth International
  !       Conference on Properties of Steam (1967).

  !    ps  = saturation pressure (pascals)
  !    ts  = saturation temperature (deg. C)
  !    tsp = estimated ts on entry and dT/dps on return

  implicit none

  PetscReal :: ts
  PetscReal :: ps
  PetscReal :: t_ps
  PetscReal :: ts_guess
  PetscErrorCode :: ierr
  

  PetscReal, parameter :: epsilon = 1.d-10
  PetscReal, parameter :: tc1 = 647.3d0
  PetscReal, parameter :: pc1 = 2.212d7  
      
  PetscReal :: theta, beta, u1, err
  PetscReal :: t1num, t1nump
  PetscReal :: t1, t1den, t1denp, term1, term1p, t2, term2, term2p
  PetscReal :: f, fp
  PetscReal :: kn(9)
  PetscInt :: itr

  data kn / -7.691234564d0,-2.608023696d1,-1.681706546d2, &
            6.423285504d1,-1.189646225d2, 4.167117320d0, &
            2.097506760d1, 1.d9         , 6.d0/

  ierr = 0

!geh  if(ipvtab.eq.0 .or. tsp.gt.369.d0) then

!-------newton-raphson iteration for calculating ts by analytical funcs

!-------on entry, ts_guess = estimated ts
  theta = (ts_guess+273.15d0)/tc1
  beta  = ps/pc1

  u1  = 1.d0-theta
  itr = 0

  do         
    itr = itr + 1

    t1num  = u1*(kn(1)+u1*(kn(2)+u1*(kn(3)+u1*(kn(4)+u1*kn(5)))))
    t1nump = u1*(2.d0*kn(2)+u1*(3.d0*kn(3)+u1*(4.d0*kn(4)+5.d0* &
             u1*kn(5))))+kn(1)
    t1     = 1.d0+u1*(kn(6)+kn(7)*u1)
    t1den  = 1.d0/(theta*t1)
    t1denp = theta*(kn(6)+2.d0*kn(7)*u1)-t1
    term1  = t1num*t1den
    term1p = (t1nump-term1*t1denp)*t1den

    t2     = 1.d0/(kn(8)*u1*u1+kn(9))
    term2  = u1*t2
    term2p = t2*(1.d0-2.d0*kn(8)*u1*u1*t2)
    f      = exp(term1-term2)-beta
    fp     = (f+beta)*(term1p-term2p)
    err    = f/fp
    u1     = u1-err
    theta  = 1.d0-u1
    
    if (dabs(err) <= epsilon .or. itr >= 20) exit

  enddo

  ts = theta*tc1-273.15d0

!-------Note-(dbeta/dtheta) = -fp ; tsp = dT/dps
  t_ps = -tc1/(pc1*fp)

end subroutine EOSWaterSaturationTemperature

! ************************************************************************** !
subroutine EOSWaterDensityIce(T, P, den_ice, dden_ice_dT, dden_ice_dP)
  ! Subroutine to calculate the density of ice at given temperature
  ! and pressure
  ! T is in deg C, P is in Pa, density is in kmol/m3
  implicit none
  
  PetscReal :: T
  PetscReal :: P
  PetscReal :: den_ice
  PetscReal :: dden_ice_dT, dden_ice_dP 
  PetscInt :: ierr
  PetscReal, parameter :: P_ref = 1.d5
  PetscReal, parameter :: alpha = 3.3d-10
  PetscReal, parameter :: beta = 1.53d-4

  den_ice = 5.09424d1*(1.d0 + alpha*(P - P_ref) - beta*(T)) !in Kmol/m3
  dden_ice_dT = 5.09424d1*(-beta)
  dden_ice_dP = 5.09424d1*alpha
  
end subroutine EOSWaterDensityIce

! ************************************************************************** !
subroutine EOSWaterInternalEnergyIce(T, u_ice, du_ice_dT)
  ! Subroutine to calculate the internal energy of ice at given temperature and 
  ! pressure
  ! T is in deg C, internal energy is in J/mol
  implicit none

  PetscReal :: T
  PetscReal :: u_ice
  PetscReal :: du_ice_dT
  PetscInt :: ierr
  PetscReal, parameter :: a = -10.6644d0
  PetscReal, parameter :: b = 0.1698d0
  PetscReal, parameter :: c = 198148.d0
  PetscReal, parameter :: T_ref = 273.15d0

  ! from Maier-Kelly type fit (integrated tref to t)
  ! in J/mol

  u_ice = a*(T) + b/2.d0*((T + T_ref)**(2.d0) - T_ref**(2.d0)) + &
          c*(1.d0/T_ref - 1.d0/(T + T_ref))
  u_ice = u_ice - HEAT_OF_FUSION*FMWH2O*1.d-3   ! kJ/kmol
  du_ice_dT = a + b*(T + T_ref) + c/((T + T_ref)**(2.d0)) !kJ/kmol/K
  
end subroutine EOSWaterInternalEnergyIce

! ************************************************************************** !
subroutine EOSWaterDensityEnthalpy2(T, P, den_water_kg, den_water_kmol, &
                                    h_MJ_kmol, &
                                    calculate_derivatives, dden_water_dp, &
                                    dden_water_dt, dh_dp, dh_dt, ierr)
! wateos_simple: Simple water equation of state from Scott Painter
! author: Satish Karra, LANL
! date: 02/1/12
! T in C, P in Pa
  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: den_water_kg, den_water_kmol, h_MJ_kmol
  PetscReal, intent(out) :: dden_water_dp, dden_water_dt
  PetscReal, intent(out) :: dh_dp, dh_dt
  
  PetscErrorCode, intent(out) :: ierr  
  
  PetscReal, parameter :: a = 999.915d0
  PetscReal, parameter :: b = 0.0416516d0
  PetscReal, parameter :: c = 0.0100836d0
  PetscReal, parameter :: d = 0.000206355
  PetscReal, parameter :: alpha = 5.0d-10     ! in Pa^(-1)
  PetscReal, parameter :: T_ref = 273.15d0    ! in K
  PetscReal, parameter :: P_ref = 1.0d5       ! in Pa
 
  PetscReal :: den_w_one_bar, T_K
  PetscReal :: u_J_mol, u_J_kg, h_J_kg
  PetscReal :: du_dt

  ierr = 0

  ! Density of water
  T_K = T + T_ref    ! convert to Kelvin
  den_w_one_bar = a + b*(T_K - T_ref) + c*(T_K - T_ref)**(2.d0) + &
                  d*(T_K - T_ref)**(3.d0)
  den_water_kg = den_w_one_bar*(1 + alpha*(P - P_ref))
  den_water_kmol = den_water_kg/FMWH2O     ! in mol
  
  ! Internal energy
  u_J_mol = 76.0d0*(T_K - T_ref)        ! in J/mol
  u_J_kg = 4.217*1.0d3*(T_K - T_ref)    ! in J/kg
  h_J_kg = u_J_kg + P/den_water_kg    ! in J/kg
  h_MJ_kmol = h_J_kg*FMWH2O*1.d-6     ! in MJ/kmol
  
  if (calculate_derivatives) then
    ! Derivatives of density
    dden_water_dp = 1/FMWH2O*den_w_one_bar*alpha    ! in Kmol/Pa
    dden_water_dt = 1/FMWH2O*(1 + alpha*(P - P_ref))*(b + 2.d0*c*(T_K - T_ref) + &
                              3.d0*d*(T_K - T_ref)**(2.d0))      ! in Kmol/K
                            
    ! Derivatives of enthalpy
    dh_dp = FMWH2O*1.d-6/den_water_kg   ! in MJ/kmol/Pa
    du_dt = 4.217*1.d3                  ! in J/kg/K
    dh_dt = FMWH2O*1.d-6*(du_dt + P*(-1.d0/den_water_kg**(2.d0))* &
                          dden_water_dt*FMWH2O)    ! in MJ/kmol/K
  else
    dden_water_dp = -999.d0
    dden_water_dp = -999.d0
    dh_dp = -999.d0
    du_dt = -999.d0
    dh_dt = -999.d0
  endif
   
end subroutine EOSWaterDensityEnthalpy2

! ************************************************************************** !
subroutine EOSWaterdensity (tc,p,d)
  !-units: tc - degrees c, p - pascals, d - kg/m^3
  !-compute density and internal energy of liquid water as
  ! function of temperature and pressure, using the steam
  ! table equations as given by the international formulation
  ! committee (1967).
  implicit none

  PetscReal :: tc, p, d, par1, par2, par3, par4, par5, pnmr, tkr, &
            v, vmkr, y, z, zp, zero = 0.d0
  PetscReal :: a(23),sa(12)

  data a/ &
    6.824687741d3,  -5.422063673d2, -2.096666205d4,  3.941286787d4,&
    -6.733277739d4,   9.902381028d4, -1.093911774d5,  8.590841667d4,&
    -4.511168742d4,   1.418138926d4, -2.017271113d3,  7.982692717d0,&
    -2.616571843d-2,  1.522411790d-3, 2.284279054d-2, 2.421647003d2,&
    1.269716088d-10, 2.074838328d-7, 2.174020350d-8, 1.105710498d-9,&
    1.293441934d1,   1.308119072d-5, 6.047626338d-14/
  data sa/ &
    8.438375405d-1, 5.362162162d-4, 1.7200000d0,    7.342278489d-2,&
    4.975858870d-2, 6.537154300d-1, 1.150d-6,       1.51080d-5,&
    1.41880d-1,     7.002753165d0,  2.995284926d-4, 2.040d-1/

  tkr = (tc + 273.15d0)/647.3d0
  pnmr = p/2.212d7   ! *1.d+5    mss nov, 97
  y = 1.d0 - sa(1)*tkr*tkr - sa(2)/tkr**6
  zp = (sa(3)*y*y - 2.d0*sa(4)*tkr + 2.d0*sa(5)*pnmr)
  if (zp.lt.zero) goto 1
  z = y + sqrt(zp)
  par1 = a(12)*sa(5)/z**(5.d0/17.d0)
  par2 = a(13) + a(14)*tkr + a(15)*tkr*tkr + a(16)*(sa(6) - tkr)**10 + &
    a(17)/(sa(7) + tkr**19)
  par3 = (a(18) + 2.d0*a(19)*pnmr + 3.d0*a(20)*pnmr*pnmr)/(sa(8) + tkr**11)
  par4 = a(21)*tkr**18*(sa(9) + tkr*tkr)*(-3.d0/(sa(10) + pnmr)**4 + sa(11))
  par5 = 3.d0*a(22)*(sa(12) - tkr)*pnmr*pnmr + 4.d0*a(23)/tkr**20*pnmr**3
  vmkr = par1 + par2 - par3 - par4 + par5
  v = vmkr*3.17d-3
  d = 1.d0/v
    
  if (d.le.zero) then
    print *,'stop! negative density  = ',d    
    stop
  end if
  return
    
1 continue
  print 2,tc
2 format (' Temperature = ',d12.5,'  Out of Range in DENSITY')

end subroutine EOSWaterdensity

end module EOS_Water_module
