module General_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  ! thermodynamic state of fluid ids
  PetscInt, parameter, public :: NULL_STATE = 0
  PetscInt, parameter, public :: LIQUID_STATE = 1
  PetscInt, parameter, public :: GAS_STATE = 2
  PetscInt, parameter, public :: TWO_PHASE_STATE = 3
  PetscInt, parameter, public :: ANY_STATE = 4
  
  PetscInt, parameter, public :: PREV_TS = 1
  PetscInt, parameter, public :: PREV_IT = 2

  PetscInt, parameter, public :: GENERAL_LIQUID_PRESSURE_DOF = 1
  PetscInt, parameter, public :: GENERAL_GAS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: GENERAL_AIR_PRESSURE_DOF = 2
  PetscInt, parameter, public :: GENERAL_GAS_SATURATION_DOF = 3
  PetscInt, parameter, public :: GENERAL_LIQUID_FLUX_DOF = 1
  PetscInt, parameter, public :: GENERAL_GAS_FLUX_DOF = 1
  
  PetscInt, parameter, public :: GENERAL_LIQUID_STATE_TEMPERATURE_DOF = 3
  PetscInt, parameter, public :: GENERAL_GAS_STATE_TEMPERATURE_DOF = 3
  PetscInt, parameter, public :: GENERAL_2PHASE_STATE_TEMPERATURE_DOF = 2
  
  PetscInt, parameter, public :: GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF = 2
  
  PetscInt, parameter, public :: GENERAL_LIQUID_CONDUCTANCE_DOF = -1
  PetscInt, parameter, public :: GENERAL_GAS_CONDUCTANCE_DOF = -2
  PetscInt, parameter, public :: GENERAL_FLUX_DOF = 4

  type, public :: general_auxvar_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: den(:)    ! (iphase)
    PetscReal, pointer :: den_kg(:) ! (iphase)
    PetscReal :: temp
    PetscReal, pointer :: xmol(:,:) ! (icomp,iphase)
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
!    PetscReal, pointer :: dsat_dp(:,:)
!    PetscReal, pointer :: dden_dp(:,:)
!    PetscReal, pointer :: dsat_dt(:)
!    PetscReal, pointer :: dden_dt(:)
    PetscReal, pointer :: kvr(:)
    PetscReal :: pert
!    PetscReal, pointer :: dkvr_dp(:)
  end type general_auxvar_type
  
  type, public :: general_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
    PetscReal, pointer :: thermal_conductivity(:) ! (iphase)
  end type general_parameter_type
  
  type, public :: general_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(general_parameter_type), pointer :: general_parameter
    type(general_auxvar_type), pointer :: auxvars(:,:)
    type(general_auxvar_type), pointer :: auxvars_bc(:)
    type(general_auxvar_type), pointer :: auxvars_ss(:)
  end type general_type

  interface GeneralAuxVarDestroy
    module procedure GeneralAuxVarSingleDestroy
    module procedure GeneralAuxVarArray1Destroy
    module procedure GeneralAuxVarArray2Destroy
  end interface GeneralAuxVarDestroy
  
  interface GeneralOutputAuxVars
    module procedure GeneralOutputAuxVars1
    module procedure GeneralOutputAuxVars2
  end interface GeneralOutputAuxVars
  
  public :: GeneralAuxCreate, GeneralAuxDestroy, &
            GeneralAuxVarCompute, GeneralAuxVarInit, &
            GeneralAuxVarCopy, GeneralAuxVarDestroy, &
            GeneralAuxVarStrip, GeneralAuxVarUpdateState, &
            GeneralPrintAuxVars, GeneralOutputAuxVars

contains

! ************************************************************************** !

function GeneralAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
    
  type(general_type), pointer :: GeneralAuxCreate
  
  type(general_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  allocate(aux%general_parameter)
  allocate(aux%general_parameter%diffusion_coefficient(option%nphase))
  aux%general_parameter%diffusion_coefficient(LIQUID_PHASE) = 1.d-9
  aux%general_parameter%diffusion_coefficient(GAS_PHASE) = 1.d-5
  allocate(aux%general_parameter%thermal_conductivity(option%nphase))
  aux%general_parameter%thermal_conductivity = 0.d0

  GeneralAuxCreate => aux
  
end function GeneralAuxCreate

! ************************************************************************** !

subroutine GeneralAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%istate_store = NULL_STATE
  allocate(auxvar%pres(option%nphase+FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den(option%nphase))
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(option%nphase))
  auxvar%den_kg = 0.d0
  ! keep at 25 C.
  auxvar%temp = 25.d0
  allocate(auxvar%xmol(option%nflowspec,option%nphase))
  auxvar%xmol = 0.d0
  allocate(auxvar%H(option%nphase))
  auxvar%H = 0.d0
  allocate(auxvar%U(option%nphase))
  auxvar%U = 0.d0
  allocate(auxvar%kvr(option%nphase))
  auxvar%kvr = 0.d0
  
  auxvar%pert = 0.d0
  
end subroutine GeneralAuxVarInit

! ************************************************************************** !

subroutine GeneralAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate_store = auxvar%istate_store
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%xmol = auxvar%xmol
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%kvr = auxvar%kvr
  auxvar2%pert = auxvar%pert

end subroutine GeneralAuxVarCopy

! ************************************************************************** !

subroutine GeneralAuxVarCompute(x,gen_auxvar,global_auxvar,material_auxvar, &
                                saturation_function,ghosted_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module
  use Global_Aux_module
  use Gas_EOS_module
  use EOS_Water_module
  use Saturation_Function_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: ghosted_id

  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: cell_pressure
  PetscReal :: den_wat_vap, den_kg_wat_vap, h_wat_vap
  PetscReal :: den_air, h_air
  PetscReal :: den_gp, den_gt, hgp, hgt, dgp, dgt, u
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas
  PetscReal :: krl, visl, dkrl_Se
  PetscReal :: krg, visg, dkrg_Se
  PetscReal :: K_H_tilde
  PetscReal :: guess, dummy
  PetscInt :: apid, cpid, vpid, spid
  character(len=8) :: state_char
  PetscErrorCode :: ierr

  ! from init.F90
!  option%nphase = 2
!  option%liquid_phase = 1  ! liquid_pressure
!  option%gas_phase = 2     ! gas_pressure

!  option%air_pressure_id = 3
!  option%capillary_pressure_id = 4
!  option%vapor_pressure_id = 5

!  option%water_id = 1
!  option%air_id = 2
!  option%energy_id = 3

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  !geh gen_auxvar%temp = 0.d0
#ifdef DEBUG_GENERAL  
  gen_auxvar%H = -999.d0
  gen_auxvar%U = -999.d0
  gen_auxvar%pres = -999.d0
  gen_auxvar%sat = -999.d0
  gen_auxvar%den = -999.d0
  gen_auxvar%den_kg = -999.d0
  gen_auxvar%xmol = -999.d0
  select case(global_auxvar%istate)
    case(1)
      state_char = 'L'
    case(2)
      state_char = 'G'
    case(3)
      state_char = '2P'
  end select
#else
  gen_auxvar%H = 0.d0
  gen_auxvar%U = 0.d0
  gen_auxvar%pres = 0.d0
  gen_auxvar%sat = 0.d0
  gen_auxvar%den = 0.d0
  gen_auxvar%den_kg = 0.d0
  gen_auxvar%xmol = 0.d0
#endif  
  gen_auxvar%kvr = 0.d0

#if 0
  if (option%iflag >= 1) then
    if (option%iflag == 1) then
      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
        ghosted_id, x(1:3), trim(state_char)
    else
!      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
!        -1*ghosted_id, x(1:3), trim(state_char)
    endif
  endif
#endif
  
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      gen_auxvar%pres(lid) = x(GENERAL_LIQUID_PRESSURE_DOF)
      gen_auxvar%xmol(acid,lid) = x(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF)
      gen_auxvar%temp = x(GENERAL_LIQUID_STATE_TEMPERATURE_DOF)

      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(:,gid) = 0.d0
      gen_auxvar%sat(lid) = 1.d0
      gen_auxvar%sat(gid) = 0.d0

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                      gen_auxvar%pres(spid),ierr)
      !geh: Henry_air_xxx returns K_H in units of Pa, but I am not confident
      !     that K_H is truly K_H_tilde (i.e. p_g * K_H).
      call Henry_air_noderiv(dummy,gen_auxvar%temp, &
                             gen_auxvar%pres(spid),K_H_tilde)
      gen_auxvar%pres(gid) = gen_auxvar%pres(lid)
      gen_auxvar%pres(apid) = K_H_tilde*gen_auxvar%xmol(acid,lid)
      ! need vpres for liq -> 2ph check
      gen_auxvar%pres(vpid) = gen_auxvar%pres(lid) - gen_auxvar%pres(apid)
      gen_auxvar%pres(cpid) = 0.d0
      
    case(GAS_STATE)
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%pres(apid) = x(GENERAL_AIR_PRESSURE_DOF)
      gen_auxvar%temp = x(GENERAL_GAS_STATE_TEMPERATURE_DOF)

      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 1.d0
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / &
                                   gen_auxvar%pres(gid)
      gen_auxvar%xmol(:,lid) = 0.d0
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)
      ! we have to have a liquid pressure to counter a neighboring 
      ! liquid pressure.  Set to gas pressure.
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid)
      gen_auxvar%pres(cpid) = 0.d0
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                      gen_auxvar%pres(spid),ierr)
      
    case(TWO_PHASE_STATE)
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%pres(apid) = x(GENERAL_AIR_PRESSURE_DOF)
      gen_auxvar%sat(gid) = x(GENERAL_GAS_SATURATION_DOF)
      
      gen_auxvar%sat(lid) = 1.d0 - gen_auxvar%sat(gid)
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)
      
      gen_auxvar%pres(spid) = gen_auxvar%pres(vpid)
      guess = gen_auxvar%temp
      call EOSWaterSaturationTemperature(gen_auxvar%temp, &
                                         gen_auxvar%pres(spid),dummy, &
                                         guess,ierr)
      
      call SatFuncGetCapillaryPressure(gen_auxvar%pres(cpid), &
                                       gen_auxvar%sat(lid), &
                                       saturation_function,option)      

      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - &
                              gen_auxvar%pres(cpid)

      call Henry_air_noderiv(dummy,gen_auxvar%temp, &
                             gen_auxvar%pres(spid),K_H_tilde)
      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / &
                                   gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)

  end select

  cell_pressure = max(gen_auxvar%pres(lid),gen_auxvar%pres(gid))

  ! ALWAYS UPDATE THERMODYNAMIC PROPERTIES FOR BOTH PHASES!!!
  ! Liquid phase thermodynamic properties
  ! must use cell_pressure as the pressure, not %pres(lid)
  call EOSWaterDensityEnthalpy(gen_auxvar%temp,cell_pressure, &
                               gen_auxvar%den_kg(lid),gen_auxvar%den(lid), &
                               gen_auxvar%H(lid),option%scale,ierr)

  ! MJ/kmol comp
  gen_auxvar%U(lid) = gen_auxvar%H(lid) - &
                       ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                       (cell_pressure / gen_auxvar%den(lid) * &
                        option%scale)

  ! Gas phase thermodynamic properties
  call ideal_gaseos_noderiv(gen_auxvar%pres(apid),gen_auxvar%temp, &
                            option%scale,den_air,h_air,u)
!  call steameos(gen_auxvar%temp,gen_auxvar%pres(gid), &
!                gen_auxvar%pres(apid),den_kg_wat_vap,den_wat_vap,dgp,dgt, &
!                h_wat_vap,hgp,hgt,option%scale,ierr) 
!  call EOSWaterSteamDensityEnthalpy(gen_auxvar%temp,gen_auxvar%pres(gid), &
  call EOSWaterSteamDensityEnthalpy(gen_auxvar%temp,cell_pressure, &
                                    gen_auxvar%pres(apid),den_kg_wat_vap, &
                                    den_wat_vap,h_wat_vap,option%scale,ierr)
  
  gen_auxvar%den(gid) = den_wat_vap + den_air
  gen_auxvar%den_kg(gid) = den_kg_wat_vap + den_air*FMWAIR
  ! if xmol not set for gas phase, set based on densities
  if (gen_auxvar%xmol(acid,gid) < 1.d-40) then
    xmol_air_in_gas = den_air / gen_auxvar%den(gid)
    xmol_water_in_gas = 1.d0 - gen_auxvar%xmol(acid,gid)
  else
    xmol_air_in_gas = gen_auxvar%xmol(acid,gid)
    xmol_water_in_gas = gen_auxvar%xmol(wid,gid)
  endif
  ! MJ/kmol
  gen_auxvar%H(gid) = xmol_water_in_gas*h_wat_vap + &
                      xmol_air_in_gas*h_air
  gen_auxvar%U(gid) = gen_auxvar%H(gid) - &
                       ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
!                       (gen_auxvar%pres(gid) / gen_auxvar%den(gid) * &
                       (cell_pressure / gen_auxvar%den(gid) * &
                        option%scale)

  if (global_auxvar%istate == LIQUID_STATE .or. &
      global_auxvar%istate == TWO_PHASE_STATE) then
    ! this does not need to be calculated for LIQUID_STATE (=1)
    call SatFuncGetRelPermFromSat(gen_auxvar%sat(lid),krl,dkrl_Se, &
                                  saturation_function,lid,PETSC_FALSE,option)
!    call EOSWaterViscosity(gen_auxvar%temp,gen_auxvar%pres(lid), &
    call EOSWaterViscosity(gen_auxvar%temp,cell_pressure, &
                           gen_auxvar%pres(spid),visl,ierr)
    gen_auxvar%kvr(lid) = krl/visl
  endif

  if (global_auxvar%istate == GAS_STATE .or. &
      global_auxvar%istate == TWO_PHASE_STATE) then
    ! this does not need to be calculated for GAS_STATE (=1)
    call SatFuncGetRelPermFromSat(gen_auxvar%sat(gid),krg,dkrg_Se, &
                                  saturation_function,gid,PETSC_FALSE,option)
    call visgas_noderiv(gen_auxvar%temp,gen_auxvar%pres(apid), &
!                        gen_auxvar%pres(gid),den_air,visg)
                        cell_pressure,den_air,visg)
    gen_auxvar%kvr(gid) = krg/visg
  endif

#if 0
  if (option%iflag == 1) then
    write(*,'(a,i3,6f13.4,a3)') 'i/l/g/a/c/v/s: ', &
      ghosted_id, gen_auxvar%pres(1:5), gen_auxvar%sat(1), trim(state_char)
  endif
#endif

end subroutine GeneralAuxVarCompute

! ************************************************************************** !

subroutine GeneralAuxVarUpdateState(x,gen_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    saturation_function,ghosted_id, &
                                    option)
  ! 
  ! GeneralUpdateState: Updates the state and swaps primary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/25/11
  ! 

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use Gas_EOS_module
  use Saturation_Function_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  PetscInt :: ghosted_id
  type(saturation_function_type) :: saturation_function
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscReal, parameter :: epsilon = 1.d-8
  PetscReal :: x(option%nflowdof)
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: dummy, guess
  PetscBool :: flag
  character(len=MAXSTRINGLENGTH) :: state_change_string
  PetscErrorCode :: ierr

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  flag = PETSC_FALSE
  
  gen_auxvar%istate_store(PREV_IT) = global_auxvar%istate
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      if (gen_auxvar%pres(vpid) <= gen_auxvar%pres(spid)) then
#ifdef DEBUG_GENERAL
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,ghosted_id, &
                                 'Before Update',option)
        if (option%iflag == 1) then
          write(state_change_string,'(''Liquid -> 2 Phase at Cell '',i5)') &
            ghosted_id
        else
          write(state_change_string,'(''Liquid -> 2 Phase at Boundary Face '', &
                                    & i5)') ghosted_id
        endif
#endif      
        global_auxvar%istate = TWO_PHASE_STATE
        x(GENERAL_GAS_PRESSURE_DOF) = &
          gen_auxvar%pres(lid) * (1.d0 + epsilon)
!        x(GENERAL_AIR_PRESSURE_DOF) = &
!          gen_auxvar%pres(apid) * (1.d0 + epsilon)
        x(GENERAL_AIR_PRESSURE_DOF) = &
          ! pa = pg - ps
          x(GENERAL_GAS_PRESSURE_DOF) - gen_auxvar%pres(spid)
        x(GENERAL_GAS_SATURATION_DOF) = epsilon
        flag = PETSC_TRUE
      endif
    case(GAS_STATE)
      if (gen_auxvar%pres(vpid) >= gen_auxvar%pres(spid)) then
#ifdef DEBUG_GENERAL
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,ghosted_id, &
                                 'Before Update',option)
        if (option%iflag == 1) then
          write(state_change_string,'(''Gas -> 2 Phase at Cell '',i5)') &
            ghosted_id
        else
          write(state_change_string,'(''Gas -> 2 Phase at Boundary Face '', &
                                    & i5)') ghosted_id
        endif
#endif      
        global_auxvar%istate = TWO_PHASE_STATE
        ! first two primary dependent variables do not change
        x(GENERAL_GAS_SATURATION_DOF) = 1.d0 - epsilon
        flag = PETSC_TRUE
      endif
    case(TWO_PHASE_STATE)
      if (gen_auxvar%sat(gid) < 0.d0) then
#ifdef DEBUG_GENERAL
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,ghosted_id, &
                                 'Before Update',option)
        if (option%iflag == 1) then
          write(state_change_string,'(''2 Phase -> Liquid at Cell '',i5)') &
            ghosted_id
        else
          write(state_change_string,'(''2 Phase -> Liquid at Boundary Face '', &
                                    & i5)') ghosted_id
        endif
#endif      
        ! convert to liquid state
        global_auxvar%istate = LIQUID_STATE
        x(GENERAL_LIQUID_PRESSURE_DOF) = &
          gen_auxvar%pres(gid) * (1.d0 + epsilon)
        x(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF) = &
          gen_auxvar%xmol(acid,lid) * (1+epsilon)
        x(GENERAL_LIQUID_STATE_TEMPERATURE_DOF) = &
          gen_auxvar%temp * (1.d0 - epsilon)
        flag = PETSC_TRUE
      else if (gen_auxvar%sat(gid) > 1.d0) then
#ifdef DEBUG_GENERAL
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,ghosted_id, &
                                 'Before Update',option)
        if (option%iflag == 1) then
          write(state_change_string,'(''2 Phase -> Gas at Cell '',i5)') &
            ghosted_id
        else
          write(state_change_string,'(''2 Phase -> Gas at Boundary Face '', &
                                    & i5)') ghosted_id
        endif
#endif      
        ! convert to gas state
        global_auxvar%istate = GAS_STATE
        ! first two primary dependent variables do not change
        x(GENERAL_GAS_STATE_TEMPERATURE_DOF) = &
          gen_auxvar%temp * (1.d0 + epsilon)
        flag = PETSC_TRUE
      endif
  end select
  
  if (flag) then
    call GeneralAuxVarCompute(x,gen_auxvar, global_auxvar,material_auxvar, &
                              saturation_function,ghosted_id,option)
#ifdef DEBUG_GENERAL
    call printMsg(option,state_change_string)
    call GeneralPrintAuxVars(gen_auxvar,global_auxvar,ghosted_id, &
                             'After Update',option)
#endif
    option%variables_swapped = PETSC_TRUE
  endif

end subroutine GeneralAuxVarUpdateState

! ************************************************************************** !

subroutine GeneralPrintAuxVars(general_auxvar,global_auxvar,ghosted_id, &
                               string,option)
  ! 
  ! Prints out the contents of an auxvar
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: ghosted_id
  character(len=*) :: string
  type(option_type) :: option

  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '             cell id: ', ghosted_id
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      print *, ' Thermodynamic state: Liquid phase'
    case(GAS_STATE)
      print *, ' Thermodynamic state: Gas phase'
    case(TWO_PHASE_STATE)
      print *, ' Thermodynamic state: Two phase'
  end select
  print *, '     liquid pressure: ', general_auxvar%pres(lid)
  print *, '        gas pressure: ', general_auxvar%pres(gid)
  print *, '        air pressure: ', general_auxvar%pres(apid)
  print *, '  capillary pressure: ', general_auxvar%pres(cpid)
  print *, '      vapor pressure: ', general_auxvar%pres(vpid)
  print *, ' saturation pressure: ', general_auxvar%pres(spid)
  print *, '   liquid saturation: ', general_auxvar%sat(lid)
  print *, '      gas saturation: ', general_auxvar%sat(gid)
  print *, 'liquid density [mol]: ', general_auxvar%den(lid)
  print *, '   gas density [mol]: ', general_auxvar%den(gid)
  print *, ' liquid density [kg]: ', general_auxvar%den_kg(lid)
  print *, '    gas density [kg]: ', general_auxvar%den_kg(gid)
  print *, '     temperature [C]: ', general_auxvar%temp
  print *, '  liquid H [MJ/kmol]: ', general_auxvar%H(lid)
  print *, '     gas H [MJ/kmol]: ', general_auxvar%H(gid)
  print *, '  liquid U [MJ/kmol]: ', general_auxvar%U(lid)
  print *, '     gas U [MJ/kmol]: ', general_auxvar%U(gid)
  print *, ' X (water in liquid): ', general_auxvar%xmol(lid,lid)
  print *, '   X (air in liquid): ', general_auxvar%xmol(gid,lid)
  print *, '    X (water in gas): ', general_auxvar%xmol(lid,gid)
  print *, '      X (air in gas): ', general_auxvar%xmol(gid,gid)
  print *, '          liquid kvr: ', general_auxvar%kvr(lid)
  print *, '             gas kvr: ', general_auxvar%kvr(gid)
  print *, '--------------------------------------------------------'

end subroutine GeneralPrintAuxVars

! ************************************************************************** !

subroutine GeneralOutputAuxVars1(general_auxvar,global_auxvar,ghosted_id, &
                                string,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: ghosted_id
  character(len=*) :: string
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2
  PetscInt :: apid, cpid, vpid
  PetscInt :: gid, lid, acid, wid, eid

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  write(string2,*) ghosted_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  open(unit=86,file=string2)

  write(86,*) '--------------------------------------------------------'
  write(86,*) trim(string)
  write(86,*) '             cell id: ', ghosted_id
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      write(86,*) ' Thermodynamic state: Liquid phase'
    case(GAS_STATE)
      write(86,*) ' Thermodynamic state: Gas phase'
    case(TWO_PHASE_STATE)
      write(86,*) ' Thermodynamic state: Two phase'
  end select
  write(86,*) '     liquid pressure: ', general_auxvar%pres(lid)
  write(86,*) '        gas pressure: ', general_auxvar%pres(gid)
  write(86,*) '        air pressure: ', general_auxvar%pres(apid)
  write(86,*) '  capillary pressure: ', general_auxvar%pres(cpid)
  write(86,*) '      vapor pressure: ', general_auxvar%pres(vpid)
  write(86,*) '     temperature [C]: ', general_auxvar%temp
  write(86,*) '   liquid saturation: ', general_auxvar%sat(lid)
  write(86,*) '      gas saturation: ', general_auxvar%sat(gid)
  write(86,*) 'liquid density [mol]: ', general_auxvar%den(lid)
  write(86,*) ' liquid density [kg]: ', general_auxvar%den_kg(lid)
  write(86,*) '   gas density [mol]: ', general_auxvar%den(gid)
  write(86,*) '    gas density [kg]: ', general_auxvar%den_kg(gid)
  write(86,*) ' X (water in liquid): ', general_auxvar%xmol(lid,lid)
  write(86,*) '   X (air in liquid): ', general_auxvar%xmol(gid,lid)
  write(86,*) '    X (water in gas): ', general_auxvar%xmol(lid,gid)
  write(86,*) '      X (air in gas): ', general_auxvar%xmol(gid,gid)
  write(86,*) '  liquid H [MJ/kmol]: ', general_auxvar%H(lid)
  write(86,*) '     gas H [MJ/kmol]: ', general_auxvar%H(gid)
  write(86,*) '  liquid U [MJ/kmol]: ', general_auxvar%U(lid)
  write(86,*) '     gas U [MJ/kmol]: ', general_auxvar%U(gid)
  write(86,*) '          liquid kvr: ', general_auxvar%kvr(lid)
  write(86,*) '             gas kvr: ', general_auxvar%kvr(gid)
  write(86,*) '...'
  write(86,*) general_auxvar%pres(lid)
  write(86,*) general_auxvar%pres(gid)
  write(86,*) general_auxvar%pres(apid)
  write(86,*) general_auxvar%pres(cpid)
  write(86,*) general_auxvar%pres(vpid)
  write(86,*) general_auxvar%temp
  write(86,*) general_auxvar%sat(lid)
  write(86,*) general_auxvar%sat(gid)
  write(86,*) general_auxvar%den(lid)
  write(86,*) general_auxvar%den_kg(lid)
  write(86,*) general_auxvar%den(gid)
  write(86,*) general_auxvar%den_kg(gid)
  write(86,*) general_auxvar%xmol(lid,lid)
  write(86,*) general_auxvar%xmol(gid,lid)
  write(86,*) general_auxvar%xmol(lid,gid)
  write(86,*) general_auxvar%xmol(gid,gid)
  write(86,*) general_auxvar%H(lid)
  write(86,*) general_auxvar%H(gid)
  write(86,*) general_auxvar%U(lid)
  write(86,*) general_auxvar%U(gid)
  write(86,*) ''
  write(86,*) general_auxvar%kvr(lid)
  write(86,*) general_auxvar%kvr(gid)
  write(86,*) '--------------------------------------------------------'
  
  close(86)

end subroutine GeneralOutputAuxVars1

! ************************************************************************** !

subroutine GeneralOutputAuxVars2(general_auxvars,global_auxvars,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvars(0:,:)
  type(global_auxvar_type) :: global_auxvars(:)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: apid, cpid, vpid
  PetscInt :: gid, lid, acid, wid, eid
  PetscInt :: i, n, idof

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  string = 'general_auxvar.txt'
  open(unit=86,file=string)
  
  n = size(global_auxvars)

100 format(a,100('','',i9))
  
  write(86,'(a,100('','',i9))') '             cell id: ', &
    ((i,i=1,n),idof=0,3)
  write(86,'(a,100('','',i2))') '                idof: ', &
    ((idof,i=1,n),idof=0,3)
  write(86,'(a,100('','',i2))') '               state: ', &
    (global_auxvars(i)%istate,i=1,n)
  write(86,100) '     liquid pressure: ', &
    ((general_auxvars(idof,i)%pres(lid),i=1,n),idof=0,3)
  write(86,100) '        gas pressure: ', &
    ((general_auxvars(idof,i)%pres(gid),i=1,n),idof=0,3)
  write(86,100) '        air pressure: ', &
    ((general_auxvars(idof,i)%pres(apid),i=1,n),idof=0,3)
  write(86,100) '  capillary pressure: ', &
    ((general_auxvars(idof,i)%pres(cpid),i=1,n),idof=0,3)
  write(86,100) '      vapor pressure: ', &
    ((general_auxvars(idof,i)%pres(vpid),i=1,n),idof=0,3)
  write(86,100) '     temperature [C]: ', &
    ((general_auxvars(idof,i)%temp,i=1,n),idof=0,3)
  write(86,100) '   liquid saturation: ', &
    ((general_auxvars(idof,i)%sat(lid),i=1,n),idof=0,3)
  write(86,100) '      gas saturation: ', &
    ((general_auxvars(idof,i)%sat(gid),i=1,n),idof=0,3)
  write(86,100) 'liquid density [mol]: ', &
    ((general_auxvars(idof,i)%den(lid),i=1,n),idof=0,3)
  write(86,100) ' liquid density [kg]: ', &
    ((general_auxvars(idof,i)%den_kg(lid),i=1,n),idof=0,3)
  write(86,100) '   gas density [mol]: ', &
    ((general_auxvars(idof,i)%den(gid),i=1,n),idof=0,3)
  write(86,100) '    gas density [kg]: ', &
    ((general_auxvars(idof,i)%den_kg(gid),i=1,n),idof=0,3)
  write(86,100) ' X (water in liquid): ', &
    ((general_auxvars(idof,i)%xmol(lid,lid),i=1,n),idof=0,3)
  write(86,100) '   X (air in liquid): ', &
    ((general_auxvars(idof,i)%xmol(gid,lid),i=1,n),idof=0,3)
  write(86,100) '    X (water in gas): ', &
    ((general_auxvars(idof,i)%xmol(lid,gid),i=1,n),idof=0,3)
  write(86,100) '      X (air in gas): ', &
    ((general_auxvars(idof,i)%xmol(gid,gid),i=1,n),idof=0,3)
  write(86,100) '  liquid H [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%H(lid),i=1,n),idof=0,3)
  write(86,100) '     gas H [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%H(gid),i=1,n),idof=0,3)
  write(86,100) '  liquid U [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%U(lid),i=1,n),idof=0,3)
  write(86,100) '     gas U [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%U(gid),i=1,n),idof=0,3)
  write(86,*)
  write(86,100) '          liquid kvr: ', &
    ((general_auxvars(idof,i)%kvr(lid),i=1,n),idof=0,3)
  write(86,100) '             gas kvr: ', &
    ((general_auxvars(idof,i)%kvr(gid),i=1,n),idof=0,3)
  
  close(86)

end subroutine GeneralOutputAuxVars2

! ************************************************************************** !

subroutine GeneralAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(general_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call GeneralAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)  

end subroutine GeneralAuxVarSingleDestroy

! ************************************************************************** !

subroutine GeneralAuxVarArray1Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(general_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call GeneralAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine GeneralAuxVarArray1Destroy

! ************************************************************************** !

subroutine GeneralAuxVarArray2Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(general_auxvar_type), pointer :: auxvars(:,:)
  
  PetscInt :: iaux, idof
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call GeneralAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine GeneralAuxVarArray2Destroy

! ************************************************************************** !

subroutine GeneralAuxVarStrip(auxvar)
  ! 
  ! GeneralAuxVarDestroy: Deallocates a general auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  implicit none

  type(general_auxvar_type) :: auxvar
  
  if (associated(auxvar%pres)) deallocate(auxvar%pres)
  nullify(auxvar%pres)
  if (associated(auxvar%sat)) deallocate(auxvar%sat)
  nullify(auxvar%sat)
  if (associated(auxvar%den)) deallocate(auxvar%den)
  nullify(auxvar%den)
  if (associated(auxvar%den_kg)) deallocate(auxvar%den_kg)
  nullify(auxvar%den_kg)
  if (associated(auxvar%xmol)) deallocate(auxvar%xmol)
  nullify(auxvar%xmol)
  if (associated(auxvar%H)) deallocate(auxvar%H)
  nullify(auxvar%H)
  if (associated(auxvar%U)) deallocate(auxvar%U)
  nullify(auxvar%U)
  if (associated(auxvar%kvr)) deallocate(auxvar%kvr)
  nullify(auxvar%kvr)
  
end subroutine GeneralAuxVarStrip

! ************************************************************************** !

subroutine GeneralAuxDestroy(aux)
  ! 
  ! Deallocates a general auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  implicit none

  type(general_type), pointer :: aux
  PetscInt :: iaux, idof
  
  if (.not.associated(aux)) return
  
  call GeneralAuxVarDestroy(aux%auxvars)
  call GeneralAuxVarDestroy(aux%auxvars_bc)
  call GeneralAuxVarDestroy(aux%auxvars_ss)

  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

  if (associated(aux%general_parameter)) then
    if (associated(aux%general_parameter%diffusion_coefficient)) &
      deallocate(aux%general_parameter%diffusion_coefficient)
    nullify(aux%general_parameter%diffusion_coefficient)
    if (associated(aux%general_parameter%thermal_conductivity)) &
      deallocate(aux%general_parameter%thermal_conductivity)
    nullify(aux%general_parameter%thermal_conductivity)
    deallocate(aux%general_parameter)
  endif
  nullify(aux%general_parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GeneralAuxDestroy

end module General_Aux_module
