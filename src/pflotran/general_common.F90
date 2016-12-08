module General_Common_module

  use General_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

#define CONVECTION
#define DIFFUSION
#define LIQUID_DIFFUSION
#define CONDUCTION
  
!#define DEBUG_GENERAL_FILEOUTPUT
!#define DEBUG_FLUXES  

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

#ifdef DEBUG_GENERAL_FILEOUTPUT
  PetscInt, parameter :: debug_unit = 87
  PetscInt, parameter :: debug_info_unit = 86
  character(len=MAXWORDLENGTH) :: debug_filename
  PetscInt :: debug_flag = 0
  PetscInt :: debug_iteration_count
  PetscInt :: debug_timestep_cut_count
  PetscInt :: debug_timestep_count
#endif

  public :: GeneralAccumulation, &
            GeneralFlux, &
            GeneralBCFlux, &
            GeneralSrcSink, &
            GeneralAccumDerivative, &
            GeneralFluxDerivative, &
            GeneralBCFluxDerivative, &
            GeneralSrcSinkDerivative

contains

! ************************************************************************** !

subroutine GeneralAccumulation(gen_auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res,Jac, &
                               analytical_derivatives,debug_cell)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none

  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_cell
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: porosity
  PetscReal :: volume_over_dt
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use gen_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = gen_auxvar%effective_porosity
  
  ! accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] * 
    !                           xmol[kmol comp/kmol phase]
    do icomp = 1, option%nflowspec
#ifdef DEBUG_GENERAL
      ! for debug version, aux var entries are initialized to NaNs.  even if
      ! saturation is zero, density may be a NaN.  So the conditional prevents
      ! this calculation.  For non-debug, aux var entries are initialized to
      ! 0.d0
      if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
      Res(icomp) = Res(icomp) + gen_auxvar%sat(iphase) * &
                                gen_auxvar%den(iphase) * &
                                gen_auxvar%xmol(icomp,iphase)
#ifdef DEBUG_GENERAL
      endif
#endif
    enddo
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * volume_over_dt

  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
#ifdef DEBUG_GENERAL
    ! for debug version, aux var entries are initialized to NaNs.  even if
    ! saturation is zero, density may be a NaN.  So the conditional prevents
    ! this calculation.  For non-debug, aux var entries are initialized to
    ! 0.d0
    if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
    Res(energy_id) = Res(energy_id) + gen_auxvar%sat(iphase) * &
                                      gen_auxvar%den(iphase) * &
                                      gen_auxvar%U(iphase)
#ifdef DEBUG_GENERAL
    endif
#endif
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * gen_auxvar%temp) * volume_over_dt
  
  if (analytical_derivatives) then
    Jac = 0.d0
    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        ! satl = 1
        ! ----------
        ! Water Equation
        ! por * satl * denl * Xwl
        ! ---
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Xwl + 
        ! por * ddenl_dpl * Xwl
        Jac(1,1) = &
          gen_auxvar%d%por_pl * gen_auxvar%den(1) * gen_auxvar%xmol(1,1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%xmol(1,1)
        ! w/respect to air mole fraction
        ! liquid phase density is indepenent of air mole fraction
        ! por * denl * dXwl_dXal
        ! Xwl = 1. - Xal
        ! dXwl_dXal = -1.
        Jac(1,2) = porosity * gen_auxvar%den(1) * (-1.d0)
        ! w/repect to temperature
        ! por * ddenl_dT * Xwl
        Jac(1,3) = porosity * gen_auxvar%d%denl_T * gen_auxvar%xmol(1,1)
        ! ----------
        ! Air Equation
        ! por * satl * denl * Xal
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Xal + 
        ! por * ddenl_dpl * Xal
        Jac(2,1) = &
          gen_auxvar%d%por_pl * gen_auxvar%den(1) * gen_auxvar%xmol(2,1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%xmol(2,1)
        ! w/respect to air mole fraction
        Jac(2,2) = porosity * gen_auxvar%den(1)
        ! w/repect to temperature
        ! por * ddenl_dT * Xwl
        Jac(2,3) = porosity * gen_auxvar%d%denl_T * gen_auxvar%xmol(2,1)
        ! ----------
        ! Energy Equation
        ! por * satl * denl * Ul + (1-por) * dens * Cp * T
        ! w/respect to liquid pressure
        ! dpor_dpl * denl * Ul + 
        ! por * ddenl_dpl * Ul + 
        ! por * denl * dUl_dpl + 
        ! -dpor_dpl * dens * Cp * T
        Jac(3,1) = &
          gen_auxvar%d%por_pl * gen_auxvar%den(1) * gen_auxvar%U(1) + &
          porosity * gen_auxvar%d%denl_pl * gen_auxvar%U(1) + &
          porosity * gen_auxvar%den(1) * gen_auxvar%d%Ul_pl + &
          (-1.d0) * gen_auxvar%d%por_pl * &
            material_auxvar%soil_particle_density * &
            soil_heat_capacity * gen_auxvar%temp
        ! w/respect to air mole fraction
        Jac(3,2) = 0.d0
        ! w/respect to temperature
        Jac(3,3) = &
          porosity * gen_auxvar%den(1) * gen_auxvar%d%Ul_T + &
          (1.d0 - porosity) * material_auxvar%soil_particle_density * &
            soil_heat_capacity
      case(GAS_STATE)
      case(TWO_PHASE_STATE)
!        if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
    end select
    Jac = Jac * volume_over_dt
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,7es24.15)') 'accum:', Res
  endif
#endif                    

end subroutine GeneralAccumulation

! ************************************************************************** !

subroutine GeneralFlux(gen_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       sir_up, &
                       thermal_conductivity_up, &
                       gen_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area, dist, general_parameter, &
                       option,v_darcy,Res,Jup,Jdn, &
                       analytical_derivatives, &
                       debug_connection)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  use Fracture_module
  use Klinkenberg_module
  
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  
  PetscReal :: xmol(option%nflowspec)
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass
  PetscReal :: delta_X_whatever
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  
  PetscReal :: dummy_dperm_up, dummy_dperm_dn
  PetscReal :: temp_perm_up, temp_perm_dn

  PetscReal :: dden_up, dden_dn
  PetscReal :: dden_dden_kg_up, dden_dden_kg_dn
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  PetscReal :: dmobility_dpup, dmobility_dpdn
  PetscReal :: dmobility_dsatup, dmobility_dsatdn
  PetscReal :: dmobility_dTup, dmobility_dTdn
  PetscReal :: dmole_flux_dpup, dmole_flux_dpdn
  PetscReal :: dmole_flux_dTup, dmole_flux_dTdn
  PetscReal :: dv_darcy_dpup, dv_darcy_dpdn
  PetscReal :: dv_darcy_dTup, dv_darcy_dTdn
  PetscReal :: duH_dpup, duH_dpdn
  PetscReal :: duH_dTup, duH_dTdn
  PetscReal :: dxmol_up, dxmol_dn
   
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_up%fracture)) then
    call FracturePermEvaluate(material_auxvar_up,perm_up,temp_perm_up, &
                              dummy_dperm_up,dist)
    perm_up = temp_perm_up
  endif
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif
  
  if (associated(klinkenberg)) then
    perm_ave_over_dist(1) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
    temp_perm_up = klinkenberg%Evaluate(perm_up, &
                                         gen_auxvar_up%pres(option%gas_phase))
    temp_perm_dn = klinkenberg%Evaluate(perm_dn, &
                                         gen_auxvar_dn%pres(option%gas_phase))
    perm_ave_over_dist(2) = (temp_perm_up * temp_perm_dn) / &
                            (dist_up*temp_perm_dn + dist_dn*temp_perm_up)
  else
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
  endif
      
  Res = 0.d0
  
  v_darcy = 0.d0
#ifdef DEBUG_FLUXES  
  adv_flux = 0.d0
  diff_flux = 0.d0
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

#ifdef CONVECTION
  do iphase = 1, option%nphase
 
    if (gen_auxvar_up%mobility(iphase) + &
        gen_auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif

    density_kg_ave = GeneralAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           gen_auxvar_up%den_kg, &
                                           gen_auxvar_dn%den_kg, &
                                           dden_up,dden_dn)
    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = gen_auxvar_up%pres(iphase) - &
                     gen_auxvar_dn%pres(iphase) + &
                     gravity_term

#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_dphi(iphase) = delta_pressure
#endif

    if (delta_pressure >= 0.D0) then
      mobility = gen_auxvar_up%mobility(iphase)
      xmol(:) = gen_auxvar_up%xmol(:,iphase)
      H_ave = gen_auxvar_up%H(iphase)
      uH = H_ave
    else
      mobility = gen_auxvar_dn%mobility(iphase)
      xmol(:) = gen_auxvar_dn%xmol(:,iphase)
      H_ave = gen_auxvar_dn%H(iphase)
      uH = H_ave
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = GeneralAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          gen_auxvar_up%den, &
                                          gen_auxvar_dn%den, &
                                          dden_up,dden_dn)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave
      ! Res[kmol total/sec]
      do icomp = 1, option%nflowspec
        ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
        !                      xmol[kmol comp/kmol phase]
        Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      enddo
#ifdef DEBUG_FLUXES  
      do icomp = 1, option%nflowspec
        adv_flux(icomp) = adv_flux(icomp) + mole_flux * xmol(icomp)
      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      do icomp = 1, option%nflowspec
        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
#endif
      Res(energy_id) = Res(energy_id) + mole_flux * uH
#ifdef DEBUG_FLUXES  
      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_dphi(iphase) = delta_pressure
      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
#endif
    endif                   

  enddo
! CONVECTION
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'delta pressure :', debug_dphi(:)
    write(debug_unit,'(a,7es24.15)') 'adv flux (liquid):', debug_flux(:,1)
    write(debug_unit,'(a,7es24.15)') 'adv flux (gas):', debug_flux(:,2)
  endif
  debug_flux = 0.d0
#endif                    

#ifdef DIFFUSION
  ! add in gas component diffusion in gas and liquid phases
  do iphase = 1, option%nphase
    
#ifndef LIQUID_DIFFUSION  
    if (iphase == LIQUID_PHASE) cycle
#endif    
    
    sat_up = gen_auxvar_up%sat(iphase)
    sat_dn = gen_auxvar_dn%sat(iphase)
    !geh: changed to .and. -> .or.
    if (sqrt(sat_up*sat_dn) < eps) cycle
    if (sat_up > eps .or. sat_dn > eps) then
      ! for now, if liquid state neighboring gas, we allow for minute
      ! diffusion in liquid phase.
      if (iphase == option%liquid_phase) then
        if ((sat_up > eps .or. sat_dn > eps)) then
          sat_up = max(sat_up,eps)
          sat_dn = max(sat_dn,eps)
        endif
      endif
      if (general_harmonic_diff_density) then
        den_up = gen_auxvar_up%den(iphase)
        den_dn = gen_auxvar_dn%den(iphase)
      else
        ! we use upstream weighting when iphase is not equal, otherwise
        ! arithmetic with 50/50 weighting
        den_up = GeneralAverageDensity(iphase, &
                                       global_auxvar_up%istate, &
                                       global_auxvar_dn%istate, &
                                       gen_auxvar_up%den, &
                                       gen_auxvar_dn%den, &
                                       dden_up,dden_dn)
        ! by setting both equal, we avoid the harmonic weighting below
        den_dn = den_up
      endif
      stpd_up = sat_up*material_auxvar_up%tortuosity* &
                gen_auxvar_up%effective_porosity*den_up
      stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
                gen_auxvar_dn%effective_porosity*den_dn
      if (general_diffuse_xmol) then ! delta of mole fraction
        delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                     gen_auxvar_dn%xmol(air_comp_id,iphase)
        delta_X_whatever = delta_xmol
      else ! delta of mass fraction
        xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
        xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
        xmass_air_up = xmol_air_up*fmw_comp(2) / &
                   (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
        xmass_air_dn = xmol_air_dn*fmw_comp(2) / &
                   (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
        delta_xmass = xmass_air_up - xmass_air_dn
        delta_X_whatever = delta_xmass
      endif
      ! units = [mole/m^4 bulk]
      stpd_ave_over_dist = (stpd_up*stpd_dn)/(stpd_up*dist_dn+stpd_dn*dist_up)
      ! need to account for multiple phases
      tempreal = 1.d0
      ! Eq. 1.9b.  The gas density is added below
      if (general_temp_dep_gas_air_diff .and. &
          iphase == option%gas_phase) then
        temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
        pressure_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                              gen_auxvar_dn%pres(iphase))
        tempreal = ((temp_ave+273.15d0)/273.15d0)**1.8d0 * &
                    101325.d0 / pressure_ave
      endif
      ! units = mole/sec
      mole_flux = stpd_ave_over_dist * tempreal * &
                  general_parameter%diffusion_coefficient(iphase) * &
                  delta_X_whatever * area
      Res(wat_comp_id) = Res(wat_comp_id) - mole_flux
      Res(air_comp_id) = Res(air_comp_id) + mole_flux
#ifdef DEBUG_FLUXES  
      diff_flux(wat_comp_id) = diff_flux(wat_comp_id) - mole_flux
      diff_flux(air_comp_id) = diff_flux(air_comp_id) + mole_flux      
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_flux(wat_comp_id,iphase) = debug_flux(wat_comp_id,iphase) - mole_flux 
      debug_flux(air_comp_id,iphase) = debug_flux(air_comp_id,iphase) + mole_flux 
#endif
    endif
  enddo
! DIFFUSION
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
  k_eff_up = thermal_conductivity_up(1) + &
             sqrt(gen_auxvar_up%sat(option%liquid_phase)) * &
             (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  k_eff_dn = thermal_conductivity_dn(1) + &
             sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
             (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  else
    k_eff_ave = 0.d0
  endif
  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = gen_auxvar_up%temp - gen_auxvar_dn%temp
  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux
! CONDUCTION
#endif

  if (analytical_derivatives) then
    Jup = 0.d0
    Jdn = 0.d0
    
    do iphase = 1, option%nphase
 
      if (gen_auxvar_up%mobility(iphase) + &
          gen_auxvar_dn%mobility(iphase) < eps) then
        cycle
      endif

      density_kg_ave = GeneralAverageDensity(iphase, &
                                             global_auxvar_up%istate, &
                                             global_auxvar_dn%istate, &
                                             gen_auxvar_up%den_kg, &
                                             gen_auxvar_dn%den_kg, &
                                             dden_dden_kg_up,dden_dden_kg_dn)
      gravity_term = density_kg_ave * dist_gravity
      delta_pressure = gen_auxvar_up%pres(iphase) - &
                       gen_auxvar_dn%pres(iphase) + &
                       gravity_term
                       
      ddelta_pressure_dpup = 1.d0 + dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_pl*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_pl*FMWH2O)
      ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_pl*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_pl*FMWH2O)
      ddelta_pressure_dTup = dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_T*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_T*FMWH2O)
      ddelta_pressure_dTdn = dist_gravity * &
              (dden_dden_kg_up * gen_auxvar_up%d%denl_T*FMWH2O + &
               dden_dden_kg_dn * gen_auxvar_dn%d%denl_T*FMWH2O)

      if (delta_pressure >= 0.D0) then
        mobility = gen_auxvar_up%mobility(iphase)
        xmol(:) = gen_auxvar_up%xmol(:,iphase)
        H_ave = gen_auxvar_up%H(iphase)
        uH = H_ave
        
        duH_dpup = gen_auxvar_up%d%Hl_pl
        duH_dpdn = 0.d0
        duH_dTup = gen_auxvar_up%d%Hl_T
        duH_dTdn = 0.d0
        dxmol_up = 1.d0
        dxmol_dn = 0.d0
        dmobility_dpup = gen_auxvar_up%d%mobilityl_pl
        dmobility_dsatup = gen_auxvar_up%d%mobilityl_satg
        dmobility_dTup = gen_auxvar_up%d%mobilityl_T
        dmobility_dpdn = 0.d0
        dmobility_dsatdn = 0.d0
        dmobility_dTdn = 0.d0        
      else
        mobility = gen_auxvar_dn%mobility(iphase)
        xmol(:) = gen_auxvar_dn%xmol(:,iphase)
        H_ave = gen_auxvar_dn%H(iphase)
        uH = H_ave

        duH_dpup = 0.d0
        duH_dTup = 0.d0
        duH_dTdn = gen_auxvar_dn%d%Hl_T
        dxmol_up = 0.d0
        dxmol_dn = 1.d0
        dmobility_dpup = 0.d0
        dmobility_dsatup = 0.d0
        dmobility_dTup = 0.d0
        dmobility_dpdn = gen_auxvar_dn%d%mobilityl_pl
        dmobility_dsatdn = gen_auxvar_dn%d%mobilityl_satg
        dmobility_dTdn = gen_auxvar_dn%d%mobilityl_T
      endif      

      if (mobility > floweps) then
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
        
        tempreal = perm_ave_over_dist(iphase)
        dv_darcy_dpup = tempreal * &
          (dmobility_dpup * delta_pressure + mobility * ddelta_pressure_dpup)
        dv_darcy_dTup = tempreal * &
          (dmobility_dTup * delta_pressure + mobility * ddelta_pressure_dTup)
        dv_darcy_dpdn = tempreal * &
          (dmobility_dpdn * delta_pressure + mobility * ddelta_pressure_dpdn)
        dv_darcy_dTdn = tempreal * &
          (dmobility_dTdn * delta_pressure + mobility * ddelta_pressure_dTdn)
        
        density_ave = GeneralAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            gen_auxvar_up%den, &
                                            gen_auxvar_dn%den, &
                                            dden_up,dden_dn)
        ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy(iphase) * area  
        ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
        !                             density_ave[kmol phase/m^3 phase]        
        mole_flux = q*density_ave
        ! Res[kmol total/sec]
!        do icomp = 1, option%nflowspec
!          ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
!          !                      xmol[kmol comp/kmol phase]
!          Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
!        enddo
!        Res(energy_id) = Res(energy_id) + mole_flux * uH

        select case(global_auxvar_up%istate)
          case(LIQUID_STATE)
            dmole_flux_dpup = &
              (dv_darcy_dpup * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_pl + &
                 dden_dn * gen_auxvar_dn%d%denl_pl))
            dmole_flux_dTup = &
              (dv_darcy_dTup * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_T + &
                 dden_dn * gen_auxvar_dn%d%denl_T))
            do icomp = 1, option%nflowspec
              Jup(icomp,1) = Jup(icomp,1) + dmole_flux_dpup * xmol(icomp)
              Jup(icomp,2) = Jup(icomp,2) + mole_flux * dxmol_up
              Jup(icomp,3) = Jup(icomp,3) + dmole_flux_dTup * xmol(icomp)
            enddo
            Jup(energy_id,1) = Jup(energy_id,1) + &
              (dmole_flux_dpup * uH + mole_flux * duH_dpup)
            Jup(energy_id,3) = Jup(energy_id,3) + &
              (dmole_flux_dTup * uH + mole_flux * duH_dTup)
          case(GAS_STATE)
          case(TWO_PHASE_STATE)
        end select
        select case(global_auxvar_dn%istate)
          case(LIQUID_STATE)
            dmole_flux_dpdn = &
              (dv_darcy_dpdn * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_pl + &
                 dden_dn * gen_auxvar_dn%d%denl_pl)) 
            dmole_flux_dTdn = &
              (dv_darcy_dTdn * density_ave + &
                v_darcy(iphase) * &
                (dden_up * gen_auxvar_up%d%denl_T + &
                 dden_dn * gen_auxvar_dn%d%denl_T))
            do icomp = 1, option%nflowspec
              Jdn(icomp,1) = Jdn(icomp,1) + dmole_flux_dpdn * xmol(icomp)
              Jdn(icomp,2) = Jdn(icomp,2) + mole_flux * dxmol_dn
              Jdn(icomp,3) = Jdn(icomp,3) + dmole_flux_dTdn * xmol(icomp)
            enddo
            Jdn(energy_id,1) = Jdn(energy_id,1) + &
              (dmole_flux_dpdn * uH + mole_flux * duH_dpdn)
            Jdn(energy_id,3) = Jdn(energy_id,3) + &
              (dmole_flux_dTdn * uH + mole_flux * duH_dTdn)
          case(GAS_STATE)
          case(TWO_PHASE_STATE)
        end select        
      endif                   
    enddo  
    Jup = Jup * area
    Jdn = Jdn * area

    ! add in gas component diffusion in gas and liquid phases
    do iphase = 1, option%nphase
    
      if (iphase == LIQUID_PHASE) cycle
    
      sat_up = gen_auxvar_up%sat(iphase)
      sat_dn = gen_auxvar_dn%sat(iphase)
      !geh: changed to .and. -> .or.
      if (sqrt(sat_up*sat_dn) < eps) cycle
      if (sat_up > eps .or. sat_dn > eps) then
        ! for now, if liquid state neighboring gas, we allow for minute
        ! diffusion in liquid phase.
        if (iphase == option%liquid_phase) then
          if ((sat_up > eps .or. sat_dn > eps)) then
            sat_up = max(sat_up,eps)
            sat_dn = max(sat_dn,eps)
          endif
        endif
        if (general_harmonic_diff_density) then
          den_up = gen_auxvar_up%den(iphase)
          den_dn = gen_auxvar_dn%den(iphase)
        else
          ! we use upstream weighting when iphase is not equal, otherwise
          ! arithmetic with 50/50 weighting
          den_up = GeneralAverageDensity(iphase, &
                                         global_auxvar_up%istate, &
                                         global_auxvar_dn%istate, &
                                         gen_auxvar_up%den, &
                                         gen_auxvar_dn%den, &
                                         dden_up,dden_dn)
          ! by setting both equal, we avoid the harmonic weighting below
          den_dn = den_up
        endif
        stpd_up = sat_up*material_auxvar_up%tortuosity* &
                  gen_auxvar_up%effective_porosity*den_up
        stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
                  gen_auxvar_dn%effective_porosity*den_dn
        if (general_diffuse_xmol) then
          delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                       gen_auxvar_dn%xmol(air_comp_id,iphase)
          delta_X_whatever = delta_xmol
        else
          xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
          xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
          xmass_air_up = xmol_air_up*fmw_comp(2) / &
                     (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
          xmass_air_dn = xmol_air_dn*fmw_comp(2) / &
                     (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
          delta_xmass = xmass_air_up - xmass_air_dn
          delta_X_whatever = delta_xmass
        endif
        ! units = [mole/m^4 bulk]
        stpd_ave_over_dist = (stpd_up*stpd_dn)/(stpd_up*dist_dn+stpd_dn*dist_up)
        ! need to account for multiple phases
        tempreal = 1.d0
        ! Eq. 1.9b.  The gas density is added below
        if (general_temp_dep_gas_air_diff .and. &
            iphase == option%gas_phase) then
          temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
          pressure_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                                gen_auxvar_dn%pres(iphase))
          tempreal = ((temp_ave+273.15d0)/273.15d0)**1.8d0 * &
                      101325.d0 / pressure_ave
        endif
        ! units = mole/sec
        mole_flux = stpd_ave_over_dist * tempreal * &
                    general_parameter%diffusion_coefficient(iphase) * &
                    delta_X_whatever * area
        Res(wat_comp_id) = Res(wat_comp_id) - mole_flux
        Res(air_comp_id) = Res(air_comp_id) + mole_flux
      endif
    enddo
  ! DIFFUSION
  endif
  
#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
!    write(*,'(a,7es12.4)') 'in: ', adv_flux(:)*dist(1), diff_flux(:)*dist(1)
    write(*,'('' phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(2), gen_auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(2), gen_auxvar_dn%sat(2)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,2), gen_auxvar_dn%xmol(1,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,2)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,2), gen_auxvar_dn%xmol(2,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,2)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,2) + heat_flux)*1.d6
    write(*,'('' phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(1), gen_auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(1), gen_auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,1), gen_auxvar_dn%xmol(1,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,1)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,1), gen_auxvar_dn%xmol(2,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,1)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,1) + heat_flux)*1.d6
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'dif flux (liquid):', debug_flux(:,1)
    write(debug_unit,'(a,7es24.15)') 'dif flux (gas):', debug_flux(:,2)
  endif
#endif

end subroutine GeneralFlux

! ************************************************************************** !

subroutine GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         gen_auxvar_up,global_auxvar_up, &
                         gen_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         sir_dn, &
                         thermal_conductivity_dn, &
                         area,dist,general_parameter, &
                         option,v_darcy,Res,debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module                              
  use Material_Aux_class
  use Fracture_module
  use Klinkenberg_module
  
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy(option%nphase), area
  type(general_parameter_type) :: general_parameter
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)
  PetscBool :: debug_connection
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  PetscInt :: bc_type
  PetscReal :: xmol(option%nflowspec)  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  PetscReal :: boundary_pressure
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass  
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: tempreal
  PetscReal :: delta_X_whatever

  PetscReal :: dden_dn, dden_up

  PetscInt :: idof
  PetscBool :: neumann_bc_present
  
  PetscReal :: temp_perm_dn
  PetscReal :: dummy_dperm_dn
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  Res = 0.d0
  v_darcy = 0.d0  
#ifdef DEBUG_FLUXES    
  adv_flux = 0.d0
  diff_flux = 0.d0
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

  neumann_bc_present = PETSC_FALSE
  
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif  
  
  if (associated(klinkenberg)) then
    perm_dn_adj(1) = perm_dn
                                          
    perm_dn_adj(2) = klinkenberg%Evaluate(perm_dn, &
                                          gen_auxvar_dn%pres(option%gas_phase))
  else
    perm_dn_adj(:) = perm_dn
  endif
  
#ifdef CONVECTION  
  do iphase = 1, option%nphase
 
    bc_type = ibndtype(iphase)
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(GENERAL_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(GENERAL_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
          
        ! using residual saturation cannot be correct! - geh
        ! reusing sir_dn for bounary auxvar
#define BAD_MOVE1 ! this works
#ifndef BAD_MOVE1       
        if (gen_auxvar_up%sat(iphase) > sir_dn(iphase) .or. &
            gen_auxvar_dn%sat(iphase) > sir_dn(iphase)) then
#endif
          boundary_pressure = gen_auxvar_up%pres(iphase)
          if (iphase == LIQUID_PHASE .and. &
              global_auxvar_up%istate == GAS_STATE) then
            ! the idea here is to accommodate a free surface boundary
            ! face.  this will not work for an interior grid cell as
            ! there should be capillary pressure in force.
            boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
          endif
          density_kg_ave = GeneralAverageDensity(iphase, &
                                                 global_auxvar_up%istate, &
                                                 global_auxvar_dn%istate, &
                                                 gen_auxvar_up%den_kg, &
                                                 gen_auxvar_dn%den_kg, &
                                                 dden_up, dden_dn)
          gravity_term = density_kg_ave * dist_gravity
          delta_pressure = boundary_pressure - &
                           gen_auxvar_dn%pres(iphase) + &
                           gravity_term

#ifdef DEBUG_GENERAL_FILEOUTPUT
          debug_dphi(iphase) = delta_pressure
#endif

          if (bc_type == SEEPAGE_BC .or. &
              bc_type == CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                gen_auxvar_up%pres(iphase) - &
                 option%reference_pressure < eps) then
              delta_pressure = 0.d0
            endif
          endif
            
          if (delta_pressure >= 0.D0) then
            mobility = gen_auxvar_up%mobility(iphase)
            xmol(:) = gen_auxvar_up%xmol(:,iphase)
            uH = gen_auxvar_up%H(iphase)
          else
            mobility = gen_auxvar_dn%mobility(iphase)
            xmol(:) = gen_auxvar_dn%xmol(:,iphase)
            uH = gen_auxvar_dn%H(iphase)
          endif      

          if (mobility > floweps) then
            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
            v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
            ! only need average density if velocity > 0.
            density_ave = GeneralAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                gen_auxvar_up%den, &
                                                gen_auxvar_dn%den, &
                                                dden_up,dden_dn)
          endif
#ifndef BAD_MOVE1        
        endif ! sat > eps
#endif

      case(NEUMANN_BC)
        select case(iphase)
          case(LIQUID_PHASE)
            idof = auxvar_mapping(GENERAL_LIQUID_FLUX_INDEX)
          case(GAS_PHASE)
            idof = auxvar_mapping(GENERAL_GAS_FLUX_INDEX)
        end select
        
        neumann_bc_present = PETSC_TRUE
        xmol = 0.d0
        xmol(iphase) = 1.d0
        if (dabs(auxvars(idof)) > floweps) then
          v_darcy(iphase) = auxvars(idof)
          if (v_darcy(iphase) > 0.d0) then 
            density_ave = gen_auxvar_up%den(iphase)
            uH = gen_auxvar_up%H(iphase)
          else 
            density_ave = gen_auxvar_dn%den(iphase)
            uH = gen_auxvar_dn%H(iphase)
          endif 
        endif
      case default
        option%io_buffer = &
          'Boundary condition type not recognized in GeneralBCFlux phase loop.'
        call printErrMsg(option)
    end select

    if (dabs(v_darcy(iphase)) > 0.d0) then
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                              density_ave[kmol phase/m^3 phase]
      if (density_ave < 1.d-40) then
        option%io_buffer = 'Zero density in GeneralBCFlux()'
        call printErrMsgByRank(option)
      endif
      mole_flux = q*density_ave       
      ! Res[kmol total/sec]
      do icomp = 1, option%nflowspec
        ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
        !                      xmol[kmol comp/mol phase]
        Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      enddo
#ifdef DEBUG_FLUXES  
      do icomp = 1, option%nflowspec
        adv_flux(icomp,iphase) = adv_flux(icomp,iphase) + mole_flux * xmol(icomp)
      enddo
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      do icomp = 1, option%nflowspec
        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
      enddo
#endif
      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH ! H_ave
#ifdef DEBUG_FLUXES  
      adv_flux(energy_id,iphase) = adv_flux(energy_id,iphase) + mole_flux * uH
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
#endif
    endif
  enddo
! CONVECTION
#endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then 
    write(debug_unit,'(a,7es24.15)') 'bc delta pressure :', debug_dphi(:)  
    write(debug_unit,'(a,7es24.15)') 'bc adv flux (liquid):', debug_flux(:,1)
    write(debug_unit,'(a,7es24.15)') 'bc adv flux (gas):', debug_flux(:,2)
  endif
  debug_flux = 0.d0
#endif  

#ifdef DIFFUSION
  ! add in gas component diffusion in gas and liquid phases
  do iphase = 1, option%nphase
  
#ifdef LIQUID_DIFFUSION    
!    if (neumann_bc_present) cycle
    if (ibndtype(iphase) == NEUMANN_BC) cycle
#else
    if (iphase == LIQUID_PHASE) cycle
#endif
    
    ! diffusion all depends upon the downwind cell.  phase diffusion only
    ! occurs if a phase exists in both auxvars (boundary and internal) or
    ! a liquid phase exists in the internal cell. so, one could say that
    ! liquid diffusion always exists as the internal cell has a liquid phase,
    ! but gas phase diffusion only occurs if the internal cell has a gas
    ! phase.
    if (gen_auxvar_dn%sat(iphase) > eps) then
      sat_dn = gen_auxvar_dn%sat(iphase)
      if (general_harmonic_diff_density) then
        den_dn = gen_auxvar_dn%den(iphase)
      else
        ! we use upstream weighting when iphase is not equal, otherwise
        ! arithmetic with 50/50 weighting
        den_dn = GeneralAverageDensity(iphase, &
                                       global_auxvar_up%istate, &
                                       global_auxvar_dn%istate, &
                                       gen_auxvar_up%den, &
                                       gen_auxvar_dn%den, &
                                       dden_up,dden_dn)
      endif
      ! units = [mole/m^4 bulk]
      stpd_ave_over_dist = sat_dn*material_auxvar_dn%tortuosity * &
                           gen_auxvar_dn%effective_porosity * &
                           den_dn / dist(0)
      if (general_diffuse_xmol) then
        delta_xmol = gen_auxvar_up%xmol(air_comp_id,iphase) - &
                     gen_auxvar_dn%xmol(air_comp_id,iphase)
        delta_X_whatever = delta_xmol
      else
        xmol_air_up = gen_auxvar_up%xmol(air_comp_id,iphase)
        xmol_air_dn = gen_auxvar_dn%xmol(air_comp_id,iphase)
        xmass_air_up = xmol_air_up*fmw_comp(2) / &
                   (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
        xmass_air_dn = xmol_air_dn*fmw_comp(2) / &
                   (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
        delta_xmass = xmass_air_up - xmass_air_dn
        delta_X_whatever = delta_xmass
      endif
      ! need to account for multiple phases
      ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
      tempreal = 1.d0
      ! Eq. 1.9b.  The gas density is added below
      if (general_temp_dep_gas_air_diff .and. &
          iphase == option%gas_phase) then
        temp_ave = 0.5d0*(gen_auxvar_up%temp+gen_auxvar_dn%temp)
        pres_ave = 0.5d0*(gen_auxvar_up%pres(iphase)+ &
                          gen_auxvar_dn%pres(iphase))
        tempreal = ((temp_ave+273.15d0)/273.15d0)**1.8d0 * &
                    101325.d0 / pres_ave
      endif
      ! units = mole/sec
      mole_flux = stpd_ave_over_dist * tempreal * &
                  general_parameter%diffusion_coefficient(iphase) * &
                  delta_X_whatever * area
      Res(wat_comp_id) = Res(wat_comp_id) - mole_flux
      Res(air_comp_id) = Res(air_comp_id) + mole_flux
#ifdef DEBUG_FLUXES  
      ! equal but opposite
      diff_flux(wat_comp_id,iphase) = diff_flux(wat_comp_id,iphase) - mole_flux
      diff_flux(air_comp_id,iphase) = diff_flux(air_comp_id,iphase) + mole_flux
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
      debug_flux(wat_comp_id,iphase) = debug_flux(wat_comp_id,iphase) - mole_flux
      debug_flux(air_comp_id,iphase) = debug_flux(air_comp_id,iphase) + mole_flux
#endif
    endif
  enddo
! DIFFUSION
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  heat_flux = 0.d0
  select case (ibndtype(GENERAL_ENERGY_EQUATION_INDEX))
    case (DIRICHLET_BC)
      ! based on Somerton et al., 1974:
      ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
      k_eff_dn = thermal_conductivity_dn(1) + &
                 sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
                 (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
      ! units:
      ! k_eff = W/K/m/m = J/s/K/m/m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = J/s
      k_eff_ave = k_eff_dn / dist(0)
      delta_temp = gen_auxvar_up%temp - gen_auxvar_dn%temp
      heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! convert W -> MW
    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = auxvars(auxvar_mapping(GENERAL_ENERGY_FLUX_INDEX)) * area

    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'GeneralBCFlux heat conduction loop.'
      call printErrMsg(option)
  end select
  Res(energy_id) = Res(energy_id) + heat_flux ! MW
! CONDUCTION
#endif

#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
!    write(*,'(a,7es12.4)') 'in: ', adv_flux(:)*dist(1), diff_flux(:)*dist(1)
    write(*,'('' phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(2), gen_auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(2), gen_auxvar_dn%sat(2)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,2), gen_auxvar_dn%xmol(1,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,2)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,2)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,2), gen_auxvar_dn%xmol(2,2)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,2)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,2) + heat_flux)*1.d6
    write(*,'('' phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') gen_auxvar_up%pres(1), gen_auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') gen_auxvar_up%sat(1), gen_auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(1,1), gen_auxvar_dn%xmol(1,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(1,1)
    write(*,'(''  air --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2,1)
    write(*,'(''   xmol      :'',2es12.4)') gen_auxvar_up%xmol(2,1), gen_auxvar_dn%xmol(2,1)
    write(*,'(''   diff flux :'',es12.4)') diff_flux(2,1)
    write(*,'(''  heat flux  :'',es12.4)') (adv_flux(3,1) + heat_flux)*1.d6
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'bc dif flux (liquid):', debug_flux(:,1)*dist(3)
    write(debug_unit,'(a,7es24.15)') 'bc dif flux (gas):', debug_flux(:,2)*dist(3)
  endif
#endif
  
end subroutine GeneralBCFlux

! ************************************************************************** !

subroutine GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                          gen_auxvar,global_auxvar,ss_flow_vol_flux, &
                          scale,Res,debug_cell)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
  PetscBool :: debug_cell
      
  PetscReal :: qsrc_mol
  PetscReal :: enthalpy, internal_energy
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: icomp
  PetscErrorCode :: ierr

  Res = 0.d0
  do icomp = 1, option%nflowspec
    qsrc_mol = 0.d0
    select case(flow_src_sink_type)
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(icomp)/fmw_comp(icomp) ! kg/sec -> kmol/sec
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_mol = qsrc(icomp)/fmw_comp(icomp)*scale 
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec
        qsrc_mol = qsrc(icomp)*gen_auxvar%den(icomp) ! den = kmol/m^3
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_mol = qsrc(icomp)*gen_auxvar%den(icomp)*scale 
    end select
    ! icomp here is really iphase
    ss_flow_vol_flux(icomp) = qsrc_mol/gen_auxvar%den(icomp)
    Res(icomp) = qsrc_mol
  enddo
  if (dabs(qsrc(TWO_INTEGER)) < 1.d-40 .and. &
      qsrc(ONE_INTEGER) < 0.d0) then ! extraction only
    ! Res(1) holds qsrc_mol for water.  If the src/sink value for air is zero,
    ! remove/add the equivalent mole fraction of air in the liquid phase.
    qsrc_mol = Res(ONE_INTEGER)*gen_auxvar%xmol(TWO_INTEGER,ONE_INTEGER)
    Res(TWO_INTEGER) = qsrc_mol
    ss_flow_vol_flux(TWO_INTEGER) = qsrc_mol/gen_auxvar%den(TWO_INTEGER)
  endif
  ! energy units: MJ/sec
  if (size(qsrc) == THREE_INTEGER) then
    if (dabs(qsrc(THREE_INTEGER)) < 1.d-40) then
      cell_pressure = &
        maxval(gen_auxvar%pres(option%liquid_phase:option%gas_phase))
      if (dabs(qsrc(ONE_INTEGER)) > 0.d0) then
        call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,enthalpy,ierr)
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
        ! enthalpy units: MJ/kmol                       ! water component mass
        Res(option%energy_id) = Res(option%energy_id) + Res(ONE_INTEGER) * &
                                                        enthalpy
      endif
      if (dabs(qsrc(TWO_INTEGER)) > 0.d0) then
        ! this is pure air, we use the enthalpy of air, NOT the air/water
        ! mixture in gas
        ! air enthalpy is only a function of temperature and the 
        dummy_pressure = 0.d0
        call EOSGasEnergy(gen_auxvar%temp,dummy_pressure, &
                          enthalpy,internal_energy,ierr)
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> MJ/kmol                                  
        ! enthalpy units: MJ/kmol                       ! air component mass
        Res(option%energy_id) = Res(option%energy_id) + Res(TWO_INTEGER) * &
                                                        enthalpy
      endif
    else
      Res(option%energy_id) = qsrc(THREE_INTEGER)*scale ! MJ/s
    endif
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'src/sink:', Res(1)-Res(2),Res(12:3)
  endif
#endif   
  
end subroutine GeneralSrcSink

! ************************************************************************** !

subroutine GeneralAccumDerivative(gen_auxvar,global_auxvar,material_auxvar, &
                                  soil_heat_capacity,option,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none

  type(general_auxvar_type) :: gen_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: jac(option%nflowdof,option%nflowdof)
  PetscReal :: jac_pert(option%nflowdof,option%nflowdof)
  PetscInt :: idof, irow

!geh:print *, 'GeneralAccumDerivative'

  call GeneralAccumulation(gen_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar,soil_heat_capacity,option, &
                           res,jac,general_analytical_derivatives, &
                           PETSC_FALSE)
                           
  do idof = 1, option%nflowdof
    call GeneralAccumulation(gen_auxvar(idof), &
                             global_auxvar, &
                             material_auxvar,soil_heat_capacity, &
                             option,res_pert,jac_pert,PETSC_FALSE,PETSC_FALSE)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar(idof)%pert
!geh:print *, irow, idof, J(irow,idof), gen_auxvar(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_analytical_derivatives) then
    J = jac
  endif

  if (general_isothermal) then
    J(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    J(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    J(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    J(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,10es24.15)') 'accum deriv:', J
  endif
#endif

end subroutine GeneralAccumDerivative

! ************************************************************************** !

subroutine GeneralFluxDerivative(gen_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 sir_up, &
                                 thermal_conductivity_up, &
                                 gen_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 sir_dn, &
                                 thermal_conductivity_dn, &
                                 area, dist, &
                                 general_parameter, &
                                 option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(general_auxvar_type) :: gen_auxvar_up(0:), gen_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_up(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_dn(option%nflowdof,option%nflowdof)
  PetscReal :: Jdummy(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0
  
!geh:print *, 'GeneralFluxDerivative'
  option%iflag = -2
  call GeneralFlux(gen_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up,sir_up, &
                   thermal_conductivity_up, &
                   gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn,sir_dn, &
                   thermal_conductivity_dn, &
                   area,dist,general_parameter, &
                   option,v_darcy,res,Janal_up,Janal_dn,&
                   general_analytical_derivatives,PETSC_FALSE)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_auxvar_up(idof),global_auxvar_up, &
                     material_auxvar_up,sir_up, &
                     thermal_conductivity_up, &
                     gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn,sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res_pert,Jdummy,Jdummy, &
                     PETSC_FALSE,PETSC_FALSE)
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_up(idof)%pert
!geh:print *, 'up: ', irow, idof, Jup(irow,idof), gen_auxvar_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralFlux(gen_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                     material_auxvar_up,sir_up, &
                     thermal_conductivity_up, &
                     gen_auxvar_dn(idof),global_auxvar_dn, &
                     material_auxvar_dn,sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res_pert,Jdummy,Jdummy, &
                     PETSC_FALSE,PETSC_FALSE)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_dn(idof)%pert
!geh:print *, 'dn: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_isothermal) then
    Jup(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jup(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
    Jdn(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jup(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jup(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
    Jdn(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,20es24.15)') 'flux deriv:', Jup, Jdn
  endif
#endif
  
end subroutine GeneralFluxDerivative

! ************************************************************************** !

subroutine GeneralBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   gen_auxvar_up, &
                                   global_auxvar_up, &
                                   gen_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   sir_dn, &
                                   thermal_conductivity_dn, &
                                   area,dist,general_parameter, &
                                   option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module 
  use Material_Aux_class
  
  implicit none

  PetscReal :: auxvars(:) ! from aux_real_var array
  type(general_auxvar_type) :: gen_auxvar_up, gen_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(general_parameter_type) :: general_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jdn = 0.d0
!geh:print *, 'GeneralBCFluxDerivative'

  option%iflag = -2
  call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     gen_auxvar_up,global_auxvar_up, &
                     gen_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,general_parameter, &
                     option,v_darcy,res,PETSC_FALSE)                     
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                       gen_auxvar_up,global_auxvar_up, &
                       gen_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area,dist,general_parameter, &
                       option,v_darcy,res_pert,PETSC_FALSE)   
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvar_dn(idof)%pert
!print *, 'bc: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (general_isothermal) then
    Jdn(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jdn(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jdn(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,10es24.15)') 'bc flux deriv:', Jdn
  endif
#endif
  
end subroutine GeneralBCFluxDerivative

! ************************************************************************** !

subroutine GeneralSrcSinkDerivative(option,qsrc,flow_src_sink_type, &
                                    gen_auxvars,global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(general_auxvar_type) :: gen_auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  option%iflag = -3
  call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                      gen_auxvars(ZERO_INTEGER),global_auxvar,dummy_real, &
                      scale,res,PETSC_FALSE)
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                        gen_auxvars(idof),global_auxvar,dummy_real, &
                        scale,res_pert,PETSC_FALSE)            
    do irow = 1, option%nflowdof
      Jac(irow,idof) = (res_pert(irow)-res(irow))/gen_auxvars(idof)%pert
    enddo !irow
  enddo ! idof
  
  if (general_isothermal) then
    Jac(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jac(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
  if (general_no_air) then
    Jac(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
    Jac(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
  endif  
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,20es24.15)') 'src/sink deriv:', Jac
  endif
#endif

end subroutine GeneralSrcSinkDerivative

! ************************************************************************** !

function GeneralAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn,dden_up,dden_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/14
  ! 

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)
  PetscReal :: dden_up, dden_dn

  PetscReal :: GeneralAverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0
  if (iphase == LIQUID_PHASE) then
    if (istate_up == GAS_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0
    else if (istate_dn == GAS_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == LIQUID_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0      
    else if (istate_dn == LIQUID_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0      
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0      
    endif
  endif

end function GeneralAverageDensity

end module General_Common_module
