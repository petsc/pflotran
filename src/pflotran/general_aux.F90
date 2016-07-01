module General_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: window_epsilon = 1.d-4
  PetscReal, public :: fmw_comp(2) = [FMWH2O,FMWAIR]
  PetscReal, public :: general_max_pressure_change = 5.d4
  PetscInt, public :: general_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: general_damping_factor = 0.6d0
  PetscReal, public :: general_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: general_itol_scaled_res = 1.d-5
  PetscReal, public :: general_tough2_itol_scaled_res_e1(3,3) = 1.d-5
  PetscReal, public :: general_tough2_itol_scaled_res_e2 = 1.d0
  PetscBool, public :: general_tough2_conv_criteria = PETSC_FALSE
  PetscInt, public :: general_debug_cell_id = UNINITIALIZED_INTEGER
#if defined(MATCH_TOUGH2)
  PetscBool, public :: general_temp_dep_gas_air_diff = PETSC_FALSE
  PetscBool, public :: general_diffuse_xmol = PETSC_FALSE
  PetscBool, public :: general_harmonic_diff_density = PETSC_TRUE
#else
  PetscBool, public :: general_diffuse_xmol = PETSC_TRUE
  PetscBool, public :: general_temp_dep_gas_air_diff = PETSC_TRUE
  PetscBool, public :: general_harmonic_diff_density = PETSC_TRUE
#endif

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
  PetscInt, parameter, public :: GENERAL_2PH_STATE_AIR_PRESSURE_DOF = 3
  PetscInt, parameter, public :: GENERAL_GAS_STATE_AIR_PRESSURE_DOF = 2
  PetscInt, parameter, public :: GENERAL_GAS_SATURATION_DOF = 2
  
  PetscInt, parameter, public :: GENERAL_ENERGY_DOF = 3
  PetscInt, parameter, public :: GENERAL_LIQUID_STATE_X_MOLE_DOF = 2
  
  PetscInt, parameter, public :: GENERAL_STATE_INDEX = 1
  PetscInt, parameter, public :: GENERAL_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: GENERAL_GAS_EQUATION_INDEX = 2
  PetscInt, parameter, public :: GENERAL_ENERGY_EQUATION_INDEX = 3
  
  PetscInt, parameter, public :: GENERAL_LIQUID_PRESSURE_INDEX = 2
  PetscInt, parameter, public :: GENERAL_GAS_PRESSURE_INDEX = 3
  PetscInt, parameter, public :: GENERAL_AIR_PRESSURE_INDEX = 4
  PetscInt, parameter, public :: GENERAL_MOLE_FRACTION_INDEX = 5
  PetscInt, parameter, public :: GENERAL_TEMPERATURE_INDEX = 6
  PetscInt, parameter, public :: GENERAL_GAS_SATURATION_INDEX = 7
  PetscInt, parameter, public :: GENERAL_LIQUID_FLUX_INDEX = 8
  PetscInt, parameter, public :: GENERAL_GAS_FLUX_INDEX = 9
  PetscInt, parameter, public :: GENERAL_ENERGY_FLUX_INDEX = 10
  PetscInt, parameter, public :: GENERAL_LIQUID_CONDUCTANCE_INDEX = 11
  PetscInt, parameter, public :: GENERAL_GAS_CONDUCTANCE_INDEX = 12
  PetscInt, parameter, public :: GENERAL_MAX_INDEX = 13
  
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_BOUNDARY = 2

  PetscReal, parameter, public :: general_pressure_scale = 1.d0
  
  ! these variables, which are global to general, can be modified
  PetscInt, public :: dof_to_primary_variable(3,3)
  PetscInt, public :: general_2ph_energy_dof = GENERAL_TEMPERATURE_INDEX
  PetscBool, public :: general_isothermal = PETSC_FALSE
  PetscBool, public :: general_no_air = PETSC_FALSE
  
  type, public :: general_auxvar_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal :: temp
    PetscReal, pointer :: xmol(:,:) ! (icomp,iphase)
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
!    PetscReal, pointer :: dsat_dp(:,:)
!    PetscReal, pointer :: dden_dp(:,:)
!    PetscReal, pointer :: dsat_dt(:)
!    PetscReal, pointer :: dden_dt(:)
    PetscReal, pointer :: mobility(:) ! relative perm / kinematic viscosity
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: pert
!    PetscReal, pointer :: dmobility_dp(:)
  end type general_auxvar_type
  
  type, public :: general_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
  end type general_parameter_type
  
  type, public :: general_type
    PetscInt :: n_inactive_rows
    PetscInt, pointer :: inactive_rows_local(:), inactive_rows_local_ghosted(:)
    PetscInt, pointer :: row_zeroing_array(:)

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
  
  public :: GeneralAuxCreate, &
            GeneralAuxDestroy, &
            GeneralAuxSetEnergyDOF, &
            GeneralAuxVarCompute, &
            GeneralAuxVarInit, &
            GeneralAuxVarCopy, &
            GeneralAuxVarDestroy, &
            GeneralAuxVarStrip, &
            GeneralAuxVarUpdateState, &
            GeneralAuxVarPerturb, &
            GeneralPrintAuxVars, &
            GeneralOutputAuxVars

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

  dof_to_primary_variable(1:3,1:3) = &
    reshape([GENERAL_LIQUID_PRESSURE_INDEX, GENERAL_MOLE_FRACTION_INDEX, &
             GENERAL_TEMPERATURE_INDEX, &
             GENERAL_GAS_PRESSURE_INDEX, GENERAL_AIR_PRESSURE_INDEX, &
             GENERAL_TEMPERATURE_INDEX, &
             GENERAL_GAS_PRESSURE_INDEX, GENERAL_GAS_SATURATION_INDEX, &
             general_2ph_energy_dof],shape(dof_to_primary_variable))
  
  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_inactive_rows = 0
  nullify(aux%inactive_rows_local)
  nullify(aux%inactive_rows_local_ghosted)
  nullify(aux%row_zeroing_array)

  allocate(aux%general_parameter)
  allocate(aux%general_parameter%diffusion_coefficient(option%nphase))
  !geh: there is no point in setting default lquid diffusion coeffcient values 
  !     here as they will be overwritten by the fluid property defaults.
  aux%general_parameter%diffusion_coefficient(LIQUID_PHASE) = &
                                                           UNINITIALIZED_DOUBLE
  aux%general_parameter%diffusion_coefficient(GAS_PHASE) = 2.13d-5
  aux%general_parameter%newton_inf_scaled_res_tol = 1.d-50
  aux%general_parameter%check_post_converged = PETSC_FALSE

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
  auxvar%temp = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%pert = 0.d0
  
  allocate(auxvar%pres(option%nphase+FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den(option%nphase))
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(option%nphase))
  auxvar%den_kg = 0.d0
  allocate(auxvar%xmol(option%nflowspec,option%nphase))
  auxvar%xmol = 0.d0
  allocate(auxvar%H(option%nphase))
  auxvar%H = 0.d0
  allocate(auxvar%U(option%nphase))
  auxvar%U = 0.d0
  allocate(auxvar%mobility(option%nphase))
  auxvar%mobility = 0.d0
  
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
  auxvar2%mobility = auxvar%mobility
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%pert = auxvar%pert

end subroutine GeneralAuxVarCopy

! ************************************************************************** !

subroutine GeneralAuxSetEnergyDOF(energy_keyword,option)
  ! 
  ! Sets the two phase primary dependent variable for energy based on user 
  ! input.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/25/14
  ! 
  use Option_module
  use String_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: energy_keyword
  type(option_type) :: option

  call StringToUpper(energy_keyword)
  select case(energy_keyword)
    case('TEMPERATURE')
      general_2ph_energy_dof = GENERAL_TEMPERATURE_INDEX
    case('AIR_PRESSURE')
      general_2ph_energy_dof = GENERAL_AIR_PRESSURE_INDEX
    case default
      option%io_buffer = 'Energy Keyword: ' // trim(energy_keyword) // &
                          ' not recognized in General Mode'    
      call printErrMsg(option)
  end select

end subroutine GeneralAuxSetEnergyDOF

! ************************************************************************** !

subroutine GeneralAuxVarCompute(x,gen_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Creep_Closure_module
  use Fracture_module
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id

  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: u_water_vapor, h_water_vapor
  PetscReal :: den_air, h_air, u_air
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas
  PetscReal :: krl, visl, dkrl_Se
  PetscReal :: krg, visg, dkrg_Se
  PetscReal :: K_H_tilde
  PetscReal :: guess, dummy
  PetscInt :: apid, cpid, vpid, spid
  PetscReal :: NaN
  PetscReal :: creep_closure_time
  PetscReal :: xmass_air_in_gas
  PetscReal :: xmass_air_in_liquid
  PetscReal :: Ugas_J_kg, Hgas_J_kg
  PetscReal :: Uair_J_kg, Hair_J_kg
  PetscReal :: Uvapor_J_kg, Hvapor_J_kg
  PetscReal :: Hg_mixture_fractioned  
  PetscReal :: aux(1)
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
  
#ifdef DEBUG_GENERAL  
  ! create a NaN
  NaN = 0.d0
  NaN = 1.d0/NaN
  NaN = 0.d0*NaN
  gen_auxvar%H = NaN
  gen_auxvar%U = NaN
  gen_auxvar%pres = NaN
  gen_auxvar%sat = NaN
  gen_auxvar%den = NaN
  gen_auxvar%den_kg = NaN
  gen_auxvar%xmol = NaN
  gen_auxvar%effective_porosity = NaN
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      state_char = 'L'
    case(GAS_STATE)
      state_char = 'G'
    case(TWO_PHASE_STATE)
      state_char = '2P'
  end select
#else
  !geh: do not initialize gen_auxvar%temp a the previous value is used as the
  !     initial guess for two phase.
  gen_auxvar%H = 0.d0
  gen_auxvar%U = 0.d0
  gen_auxvar%pres = 0.d0
  gen_auxvar%sat = 0.d0
  gen_auxvar%den = 0.d0
  gen_auxvar%den_kg = 0.d0
  gen_auxvar%xmol = 0.d0
  gen_auxvar%effective_porosity = 0.d0
#endif  
  gen_auxvar%mobility = 0.d0

#if 0
  if (option%iflag >= GENERAL_UPDATE_FOR_ACCUM) then
    if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
        natural_id, x(1:3), trim(state_char)
    else
!      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
!        -1*natural_id, x(1:3), trim(state_char)
    endif
  endif
#endif
  
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      gen_auxvar%pres(lid) = x(GENERAL_LIQUID_PRESSURE_DOF)
      gen_auxvar%xmol(acid,lid) = x(GENERAL_LIQUID_STATE_X_MOLE_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      ! with the gas state, we must calculate the mole fraction of air in 
      ! in the liquid phase, even though the liquid phase does not exist
      ! due to air diffusion between neighboring GAS and LIQUID cells (this
      ! is more of an issue on a boundary condition).  this is not 
      ! necessary for water since we do not calculate water diffusion 
      ! explicitly.  set mole fractions to zero in gas phase.
      gen_auxvar%xmol(:,gid) = 0.d0
      gen_auxvar%sat(lid) = 1.d0
      gen_auxvar%sat(gid) = 0.d0

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                      gen_auxvar%pres(spid),ierr)
      !geh: Henry_air_xxx returns K_H in units of Pa, but I am not confident
      !     that K_H is truly K_H_tilde (i.e. p_g * K_H).
      call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid),K_H_tilde)
      gen_auxvar%pres(gid) = max(gen_auxvar%pres(lid),gen_auxvar%pres(spid))
      gen_auxvar%pres(apid) = K_H_tilde*gen_auxvar%xmol(acid,lid)
      ! need vpres for liq -> 2ph check
      
      ! at this point, if the liquid pressure is negative, we have to go to 
      ! two phase or everything blows up:
      if (gen_auxvar%pres(gid) <= 0.d0) then
        write(option%io_buffer,'(''Negative gas pressure at cell '', &
          & i8,'' in GeneralAuxVarCompute(LIQUID_STATE).  Attempting bailout.'')') &
          natural_id
!        call printErrMsgByRank(option)
        call printMsgByRank(option)
        ! set vapor pressure to just under saturation pressure
        gen_auxvar%pres(vpid) = 0.5d0*gen_auxvar%pres(spid)
        ! set gas pressure to vapor pressure + air pressure
        gen_auxvar%pres(gid) = gen_auxvar%pres(vpid) + gen_auxvar%pres(apid)
        ! capillary pressure won't matter here.        
      else
        gen_auxvar%pres(vpid) = gen_auxvar%pres(lid) - gen_auxvar%pres(apid)
      endif
      gen_auxvar%pres(cpid) = 0.d0
      
    case(GAS_STATE)
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%pres(apid) = x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 1.d0
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / &
                                   gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
      ! need to set mole fractions in liquid phase in equilibrium with
      ! water saturated with air in order to accommodate air diffusion between
      ! GAS_STATE cell and TWO_PHASE/LIQUID_STATE cells as air should still
      ! diffuse through the liquid phase.
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                      gen_auxvar%pres(spid),ierr)
      call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid),K_H_tilde)
      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      ! set water mole fraction to zero as there is no water in liquid phase
      gen_auxvar%xmol(wid,lid) = 0.d0
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

       
      ! we have to have a liquid pressure to counter a neighboring 
      ! liquid pressure.  Set to gas pressure.
!      gen_auxvar%pres(lid) = gen_auxvar%pres(gid)
!      gen_auxvar%pres(cpid) = 0.d0

      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid),gen_auxvar%pres(cpid), &
                               option)                             
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - &
                             gen_auxvar%pres(cpid)
      
    case(TWO_PHASE_STATE)
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(gid) = x(GENERAL_GAS_SATURATION_DOF)
      if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
        gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
        call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                        gen_auxvar%pres(spid),ierr)
        gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
        gen_auxvar%pres(apid) = gen_auxvar%pres(gid) - gen_auxvar%pres(vpid)
      
      else
        gen_auxvar%pres(apid) = x(GENERAL_ENERGY_DOF)
        gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)
      
        gen_auxvar%pres(spid) = gen_auxvar%pres(vpid)
        guess = gen_auxvar%temp
        call EOSWaterSaturationTemperature(gen_auxvar%temp, &
                                           gen_auxvar%pres(spid),dummy, &
                                           guess,ierr)
      endif

      gen_auxvar%sat(lid) = 1.d0 - gen_auxvar%sat(gid)
      
      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid),gen_auxvar%pres(cpid), &
                               option)                             
!      gen_auxvar%pres(cpid) = 0.d0
 
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - &
                              gen_auxvar%pres(cpid)

      call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid),K_H_tilde)
      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      ! immiscible.
!      gen_auxvar%xmol(acid,lid) = 1.d-10
      
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / &
                                   gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)

    case default
      write(option%io_buffer,*) global_auxvar%istate
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in GeneralAuxVarCompute.'
      call printErrMsgByRank(option)

  end select

  cell_pressure = max(gen_auxvar%pres(lid),gen_auxvar%pres(gid), &
                      gen_auxvar%pres(spid))
        
  ! calculate effective porosity as a function of pressure
  if (option%iflag /= GENERAL_UPDATE_FOR_BOUNDARY) then
    gen_auxvar%effective_porosity = material_auxvar%porosity_base
#if 0
!geh this code is no longer valid
    if (associated(material_auxvar%fracture) .and. & 
      material_auxvar%fracture%setup) then
      ! The initiating pressure and maximum pressure must be calculated
      ! before fracture function applies - Heeho
      call FractureInitialSetup(material_auxvar,cell_pressure)
    endif
    if (soil_compressibility_index > 0 .and. &
      material_auxvar%setup_reference_pressure) then
      call MaterialReferencePressureSetup(material_auxvar,cell_pressure)
    endif
#endif
    if (associated(creep_closure)) then
      if (creep_closure%imat == material_auxvar%id) then
        ! option%time here is the t time, not t + dt time.
        creep_closure_time = option%time
        if (option%iflag /= GENERAL_UPDATE_FOR_FIXED_ACCUM) then
          creep_closure_time = creep_closure_time + option%flow_dt
        endif
        gen_auxvar%effective_porosity = &
          creep_closure%Evaluate(creep_closure_time,cell_pressure)
      else if (associated(material_auxvar%fracture)) then
          call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dummy)
      else if (soil_compressibility_index > 0) then
          call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dummy)
      endif
    else if (associated(material_auxvar%fracture)) then
      call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dummy)
    else if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dummy)
    endif
    if (option%iflag /= GENERAL_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = gen_auxvar%effective_porosity
    endif
  endif

  ! ALWAYS UPDATE THERMODYNAMIC PROPERTIES FOR BOTH PHASES!!!

  ! Liquid phase thermodynamic properties
  ! must use cell_pressure as the pressure, not %pres(lid)
  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(gen_auxvar%temp,cell_pressure, &
                         gen_auxvar%den_kg(lid),gen_auxvar%den(lid),ierr)
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(gen_auxvar%temp,celL_pressure,aux, &
                            gen_auxvar%den_kg(lid),gen_auxvar%den(lid),ierr)
  endif
  call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,gen_auxvar%H(lid),ierr)
  gen_auxvar%H(lid) = gen_auxvar%H(lid) * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp
  gen_auxvar%U(lid) = gen_auxvar%H(lid) - &
                       ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                       (cell_pressure / gen_auxvar%den(lid) * &
                        1.d-6)

  ! Gas phase thermodynamic properties
  ! we cannot use %pres(vpid) as vapor pressre in the liquid phase, since
  ! it can go negative
  if (global_auxvar%istate /= LIQUID_STATE) then
    if (global_auxvar%istate == GAS_STATE) then
      water_vapor_pressure = gen_auxvar%pres(vpid)
    else
      water_vapor_pressure = gen_auxvar%pres(spid)
    endif
    call EOSGasDensityEnergy(gen_auxvar%temp,gen_auxvar%pres(apid),den_air, &
                             h_air,u_air,ierr)
    h_air = h_air * 1.d-6 ! J/kmol -> MJ/kmol
    u_air = u_air * 1.d-6 ! J/kmol -> MJ/kmol
    call EOSWaterSteamDensityEnthalpy(gen_auxvar%temp,water_vapor_pressure, &
                                      den_kg_water_vapor,den_water_vapor, &
                                      h_water_vapor,ierr)
    u_water_vapor = h_water_vapor - &
                    ! Pa / kmol/m^3 = J/kmol
                    water_vapor_pressure / den_water_vapor
    h_water_vapor = h_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    u_water_vapor = u_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    gen_auxvar%den(gid) = den_water_vapor + den_air
    gen_auxvar%den_kg(gid) = den_kg_water_vapor + den_air*fmw_comp(gid)
    ! if xmol not set for gas phase, as is the case for LIQUID_STATE, 
    ! set based on densities
!    if (gen_auxvar%xmol(acid,gid) < 1.d-40) then
!      xmol_air_in_gas = den_air / gen_auxvar%den(gid)
!      xmol_water_in_gas = 1.d0 - gen_auxvar%xmol(acid,gid)
!    else
      xmol_air_in_gas = gen_auxvar%xmol(acid,gid)
      xmol_water_in_gas = gen_auxvar%xmol(wid,gid)
!    endif
      
#ifdef DEBUG_GENERAL  
    xmass_air_in_gas = xmol_air_in_gas*fmw_comp(gid) / &
                       (xmol_water_in_gas*FMWH2O + &
                        xmol_air_in_gas*fmw_comp(gid))
    Hair_J_kg = h_air*1.d6/fmw_comp(gid)
    Uair_J_kg = u_air*1.d6/fmw_comp(gid)
    Hvapor_J_kg = h_water_vapor*1.d6/FMWH2O
    Uvapor_J_kg = u_water_vapor*1.d6/FMWH2O
    Ugas_J_kg = xmass_air_in_gas*Uair_J_kg + &
                (1.d0-xmass_air_in_gas)*Uvapor_J_kg
    Hgas_J_kg = Ugas_J_kg + &
                gen_auxvar%pres(gid)/gen_auxvar%den_kg(gid)
#endif
    
    ! MJ/kmol
    gen_auxvar%U(gid) = xmol_water_in_gas * u_water_vapor + &
                        xmol_air_in_gas * u_air
    Hg_mixture_fractioned = xmol_water_in_gas*h_water_vapor + &
                            xmol_air_in_gas*h_air
    gen_auxvar%H(gid) = gen_auxvar%U(gid) + &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        gen_auxvar%pres(gid)/gen_auxvar%den(gid) * 1.d-6
  endif ! istate /= LIQUID_STATE
  
  if (global_auxvar%istate == LIQUID_STATE .or. &
      global_auxvar%istate == TWO_PHASE_STATE) then
    ! this does not need to be calculated for LIQUID_STATE (=1)
!    call SatFuncGetLiqRelPermFromSat(gen_auxvar%sat(lid),krl,dkrl_Se, &
!                                     saturation_function,lid,PETSC_FALSE,option)
    call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(gen_auxvar%sat(lid),krl,dkrl_Se,option)                            
                               
    ! use cell_pressure; cell_pressure - psat calculated internally
    if (.not.option%flow%density_depends_on_salinity) then
      call EOSWaterViscosity(gen_auxvar%temp,cell_pressure, &
                             gen_auxvar%pres(spid),visl,ierr)
    else
      aux(1) = global_auxvar%m_nacl(1)
      call EOSWaterViscosityExt(gen_auxvar%temp,cell_pressure, &
                                gen_auxvar%pres(spid),aux,visl,ierr)
    endif
    gen_auxvar%mobility(lid) = krl/visl
  endif

  if (global_auxvar%istate == GAS_STATE .or. &
      global_auxvar%istate == TWO_PHASE_STATE) then
    ! this does not need to be calculated for GAS_STATE (=1)
!    call SatFuncGetGasRelPermFromSat(gen_auxvar%sat(lid),krg, &
!                                     saturation_function,option)
    call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(gen_auxvar%sat(lid),krg,dkrg_Se,option)                            
    ! STOMP uses separate functions for calculating viscosity of vapor and
    ! and air (WATGSV,AIRGSV) and then uses GASVIS to calculate mixture 
    ! viscosity.
    call EOSGasViscosity(gen_auxvar%temp,gen_auxvar%pres(apid), &
                         gen_auxvar%pres(gid),den_air,visg,ierr)
    gen_auxvar%mobility(gid) = krg/visg
  endif

#if 0
  if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
    write(*,'(a,i3,2f13.4,es13.6,f13.4,a3)') 'i/l/g/x[a/l]/t: ', &
      natural_id, gen_auxvar%pres(1:2), gen_auxvar%xmol(acid,lid), &
      gen_auxvar%temp, trim(state_char)
    if (natural_id == 5) then
#if 0
    write(*,'(a,i3,8f13.4,a3)') 'i/l/g/a/c/v/s/t: ', &
      natural_id, gen_auxvar%pres(1:6), gen_auxvar%sat(1), gen_auxvar%temp, &
      trim(state_char)
#endif
#if 0
    if (gen_auxvar%sat(2) > 0.d0) then
      write(*,'(a,7es13.6)') 'kmol/kmol/kmol/MJ/MJ/MJ: ', &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(1,2),  &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(2,2),  &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(1,2) + &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(2,2),  &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1),  &
        gen_auxvar%sat(2)*gen_auxvar%den(2)*gen_auxvar%U(2),  &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1) + &
        gen_auxvar%sat(2)*gen_auxvar%den(2)*gen_auxvar%U(2)
    else
      write(*,'(a,7es13.6)') 'kmol/kmol/kmol/MJ/MJ/MJ: ', &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1), &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1), &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1) + &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1), &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1), 0.d0, &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1)
    endif
#endif
    endif
  endif
#endif

end subroutine GeneralAuxVarCompute

! ************************************************************************** !

subroutine GeneralAuxVarUpdateState(x,gen_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    characteristic_curves,natural_id, &
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
!  use Gas_EOS_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  class(characteristic_curves_type) :: characteristic_curves
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

! based on min_pressure in CheckPre set to zero
!  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal, parameter :: epsilon = 0.d0
  PetscReal :: liquid_epsilon
  PetscReal :: two_phase_epsilon
  PetscReal :: x(option%nflowdof)
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: dummy, guess
  PetscReal :: n_air, n_air_in_gas, n_air_in_liquid, RT, K_H_tilde, theta
  PetscBool :: flag
  character(len=MAXSTRINGLENGTH) :: state_change_string, string
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
      ! scaling by window_epsilon forces vapor pressure to enter two phase
      ! region a finite amount before phase change can occur
      if (gen_auxvar%pres(vpid) <= &
          gen_auxvar%pres(spid)*(1.d0-window_epsilon)) then
!#ifdef DEBUG_GENERAL
#ifdef DEBUG_GENERAL_INFO
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,'Before Update',option)
#endif
        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''Liquid -> 2 Phase at Cell '',i8)') &
            natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
            '(''Liquid -> 2 Phase at Cell (due to perturbation) '',i8)') &
            natural_id
        else
          write(state_change_string,'(''Liquid -> 2 Phase at Boundary Face '', &
                                    & i8)') natural_id
        endif
!#endif      
        global_auxvar%istate = TWO_PHASE_STATE
        liquid_epsilon = epsilon
!TOUGH2 approach
        ! geh: just testing something here
        x(GENERAL_GAS_PRESSURE_DOF) = max(gen_auxvar%pres(gid), &
                                          gen_auxvar%pres(spid))* &
                                         (1.d0 + liquid_epsilon)
!        x(GENERAL_GAS_PRESSURE_DOF) = max(gen_auxvar%pres(gid), &
!                                          gen_auxvar%pres(spid)+ &
!                                          gen_auxvar%pres(apid))* &
!                                         (1.d0 + liquid_epsilon)
        if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
          ! do nothing as the energy dof has not changed
          if (.not.general_isothermal) then
            X(GENERAL_ENERGY_DOF) = X(GENERAL_ENERGY_DOF) * &
                                   (1.d0 + liquid_epsilon)
          endif
        else
          ! pa = pg - ps
          x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
            x(GENERAL_GAS_PRESSURE_DOF) - gen_auxvar%pres(spid)
          if (x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) <= 0.d0) then
            write(string,*) natural_id
            option%io_buffer = 'Negative air pressure during state change ' // &
              'at ' // trim(adjustl(string))
  !          call printErrMsg(option)
            call printMsgByRank(option)
            x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
              0.01d0*x(GENERAL_GAS_PRESSURE_DOF)
          endif
        endif
        x(GENERAL_GAS_SATURATION_DOF) = 1.d-6
!        x(GENERAL_GAS_SATURATION_DOF) = liquid_epsilon
        flag = PETSC_TRUE
      endif
    case(GAS_STATE)
      ! scaling by window_epsilon forces vapor pressure to enter two phase
      ! region a finite amount before phase change can occur
      if (gen_auxvar%pres(vpid) >= &
          gen_auxvar%pres(spid)*(1.d0+window_epsilon)) then
!#ifdef DEBUG_GENERAL
#ifdef DEBUG_GENERAL_INFO
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,'Before Update',option)
#endif
        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''Gas -> 2 Phase at Cell '',i8)') &
            natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
            '(''Gas -> 2 Phase at Cell (due to perturbation) '',i8)') &
            natural_id
        else
          write(state_change_string,'(''Gas -> 2 Phase at Boundary Face '', &
                                    & i8)') natural_id
        endif
!#endif      
        global_auxvar%istate = TWO_PHASE_STATE
        ! first two primary dependent variables do not change
        if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
          if (.not.general_isothermal) then
            X(GENERAL_ENERGY_DOF) = X(GENERAL_ENERGY_DOF) * &
                                    (1.d0 - epsilon)
          endif
        else
          X(GENERAL_ENERGY_DOF) = gen_auxvar%pres(apid)*(1.d0-epsilon)
        endif
!        x(GENERAL_GAS_SATURATION_DOF) = 1.d0 - epsilon
        x(GENERAL_GAS_SATURATION_DOF) = 1.d0 - 1.d-6
        flag = PETSC_TRUE
      endif
    case(TWO_PHASE_STATE)
      if (gen_auxvar%sat(gid) < 0.d0) then
!#ifdef DEBUG_GENERAL
#ifdef DEBUG_GENERAL_INFO
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,'Before Update',option)
#endif
        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''2 Phase -> Liquid at Cell '',i8)') &
            natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
            '(''2 Phase -> Liquid at Cell (due to perturbation) '',i8)') &
            natural_id        
        else
          write(state_change_string,'(''2 Phase -> Liquid at Boundary Face '', &
                                    & i8)') natural_id
        endif
!#endif      
        two_phase_epsilon = epsilon
        ! convert to liquid state
        global_auxvar%istate = LIQUID_STATE
        x(GENERAL_LIQUID_PRESSURE_DOF) = &
          gen_auxvar%pres(gid) * (1.d0 + two_phase_epsilon) ! 4.94500E+01, 8800 NI, 10 cuts
!          gen_auxvar%pres(gid) ! 4.95500E+01, 9119 NI, 7 cuts
        x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
          gen_auxvar%xmol(acid,lid) * (1.d0 - two_phase_epsilon) ! 4.94500E+01, 8800 NI, 10 cuts
!          gen_auxvar%xmol(acid,lid) ! 4.95298E+01, 10355 NI, 6 cuts
        if (x(GENERAL_LIQUID_STATE_X_MOLE_DOF) <= 0.d0) then
          write(string,*) natural_id
          option%io_buffer = 'Negative air mole fraction during state change ' // &
            'at ' // trim(adjustl(string))
          call printMsgByRank(option)
          x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = two_phase_epsilon
        endif
        if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
          if (.not.general_isothermal) then
            X(GENERAL_ENERGY_DOF) = X(GENERAL_ENERGY_DOF) * &
                                    (1.d0-two_phase_epsilon)
          endif
        else
          X(GENERAL_ENERGY_DOF) = gen_auxvar%temp*(1.d0-two_phase_epsilon)
        endif
        flag = PETSC_TRUE
      else if (gen_auxvar%sat(gid) > 1.d0) then
!#ifdef DEBUG_GENERAL
#ifdef DEBUG_GENERAL_INFO
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,'Before Update',option)
#endif
        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''2 Phase -> Gas at Cell '',i8)') &
            natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
            '(''2 Phase -> Gas at Cell (due to perturbation) '',i8)') &
            natural_id       
        else
          write(state_change_string,'(''2 Phase -> Gas at Boundary Face '', &
                                    & i8)') natural_id
        endif
!#endif      
        two_phase_epsilon = epsilon !
        ! convert to gas state
        global_auxvar%istate = GAS_STATE
        if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
          ! first two primary dependent variables do not change
          x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
            gen_auxvar%pres(apid) * (1.d0 + two_phase_epsilon)
          if (.not.general_isothermal) then
            X(GENERAL_ENERGY_DOF) = X(GENERAL_ENERGY_DOF) * &
                                    (1.d0+two_phase_epsilon)
          endif
        else
          X(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = gen_auxvar%pres(apid)* &
                                                  (1.d0+two_phase_epsilon)
          X(GENERAL_ENERGY_DOF) = gen_auxvar%temp*(1.d0+two_phase_epsilon)
        endif
        flag = PETSC_TRUE
      endif
  end select
  
  if (flag) then
    call GeneralAuxVarCompute(x,gen_auxvar, global_auxvar,material_auxvar, &
                              characteristic_curves,natural_id,option)
!#ifdef DEBUG_GENERAL
    state_change_string = 'State Transition: ' // trim(state_change_string)
    call printMsgByRank(option,state_change_string)
#ifdef DEBUG_GENERAL_INFO
    call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                             natural_id,'After Update',option)
#endif
!#endif
  endif

end subroutine GeneralAuxVarUpdateState

! ************************************************************************** !

subroutine GeneralAuxVarPerturb(gen_auxvar,global_auxvar, &
                                material_auxvar, &
                                characteristic_curves,natural_id, &
                                option)
  ! 
  ! Calculates auxiliary variables for perturbed system
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(general_auxvar_type) :: gen_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
     
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: tempreal
!#define LEGACY_PERTURBATION
#ifdef LEGACY_PERTURBATION
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
#else
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
#endif
  PetscReal, parameter :: min_mole_fraction_pert = 1.d-12
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof

#ifdef DEBUG_GENERAL
  character(len=MAXWORDLENGTH) :: word
  type(global_auxvar_type) :: global_auxvar_debug
  type(general_auxvar_type) :: general_auxvar_debug
  call GlobalAuxVarInit(global_auxvar_debug,option)
  call GeneralAuxVarInit(general_auxvar_debug,option)
#endif

  select case(global_auxvar%istate)
    case(LIQUID_STATE)
       x(GENERAL_LIQUID_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
       x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
       x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp
#ifdef LEGACY_PERTURBATION
       ! if the liquid state, the liquid pressure will always be greater
       ! than zero.
       pert(GENERAL_LIQUID_PRESSURE_DOF) = &
         max(perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF), &
             perturbation_tolerance)
       ! if the air mole fraction perturbation is too small, the derivatives
       ! can be poor.
       pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
         -1.d0*max(perturbation_tolerance*x(GENERAL_LIQUID_STATE_X_MOLE_DOF), &
                   min_mole_fraction_pert)
       pert(GENERAL_ENERGY_DOF) = &
         -1.d0*perturbation_tolerance*x(GENERAL_ENERGY_DOF)
#else
       pert(GENERAL_LIQUID_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF) + &
         min_perturbation
       if (x(GENERAL_LIQUID_STATE_X_MOLE_DOF) > &
           1.d3 * perturbation_tolerance) then
         pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = -1.d0 * perturbation_tolerance
       else
         pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = perturbation_tolerance
       endif
       pert(GENERAL_ENERGY_DOF) = -1.d0 * &
         (perturbation_tolerance*x(GENERAL_ENERGY_DOF) + min_perturbation)
#endif
    case(GAS_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
       x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(GENERAL_ENERGY_DOF) = gen_auxvar(ZERO_INTEGER)%temp
#ifdef LEGACY_PERTURBATION
       ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a)
       ! and p(v) >= Psat (i.e. an increase in p(v)) results in two phase.
       pert(GENERAL_GAS_PRESSURE_DOF) = &
         -1.d0*perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)
       ! perturb air pressure towards gas pressure unless the perturbed
       ! air pressure exceeds the gas pressure
       if (x(GENERAL_GAS_PRESSURE_DOF) - &
           x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > &
           perturbation_tolerance* &
           x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)) then 
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
           perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
       else
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
           -1.d0*perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
       endif
       pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF)
#else
       ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a)
       ! and p(v) >= Psat (i.e. an increase in p(v)) results in two phase.
       pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 * &
         (perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF) + min_perturbation)
       ! perturb air pressure towards gas pressure unless the perturbed
       ! air pressure exceeds the gas pressure
       tempreal = perturbation_tolerance* &
                  x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) + min_perturbation
       if (x(GENERAL_GAS_PRESSURE_DOF) - &
           x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > tempreal) then
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = tempreal
       else
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
       endif
       pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF) + min_perturbation
#endif
    case(TWO_PHASE_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
!       x(GENERAL_AIR_PRESSURE_DOF) = &
!         gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
#ifdef LEGACY_PERTURBATION
       if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)
       else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > &
             perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
             perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
             -1.d0*perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)
         endif
       endif
       pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)
       ! always perturb toward 0.5
       if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then 
         pert(GENERAL_GAS_SATURATION_DOF) = &
           -1.d0*perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)
       else
         pert(GENERAL_GAS_SATURATION_DOF) = &
           perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)
       endif
#else
       if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation
       else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         tempreal = perturbation_tolerance* &
                    x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) + min_perturbation
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > tempreal) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = tempreal
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
         endif
       endif
       pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
       if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then 
         pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
       else
         pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
       endif
#endif
  end select
  
  ! GENERAL_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = GENERAL_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    gen_auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call GeneralAuxVarCompute(x_pert,gen_auxvar(idof),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id,option)
#ifdef DEBUG_GENERAL
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_debug,option)
    call GeneralAuxVarCopy(gen_auxvar(idof),general_auxvar_debug,option)
    call GeneralAuxVarUpdateState(x_pert,general_auxvar_debug, &
                                  global_auxvar_debug, &
                                  material_auxvar, &
                                  characteristic_curves, &
                                  natural_id,option)
    if (global_auxvar%istate /= global_auxvar_debug%istate) then
      write(option%io_buffer, &
            &'(''Change in state due to perturbation: '',i3,'' -> '',i3, &
            &'' at cell '',i3,'' for dof '',i3)') &
        global_auxvar%istate, global_auxvar_debug%istate, natural_id, idof
      call printMsg(option)
      write(option%io_buffer,'(''orig: '',6es17.8)') x(1:3)
      call printMsg(option)
      write(option%io_buffer,'(''pert: '',6es17.8)') x_pert_save(1:3)
      call printMsg(option)
    endif
#endif

  enddo

  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert / general_pressure_scale
    case(TWO_PHASE_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / general_pressure_scale
      if (general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
        gen_auxvar(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)%pert = &
          gen_auxvar(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)%pert / &
          general_pressure_scale
      endif
    case(GAS_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / general_pressure_scale
      gen_auxvar(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)%pert / general_pressure_scale
  end select
  
#ifdef DEBUG_GENERAL
  call GlobalAuxVarStrip(global_auxvar_debug)
  call GeneralAuxVarStrip(general_auxvar_debug)
#endif
 
end subroutine GeneralAuxVarPerturb

! ************************************************************************** !

subroutine GeneralPrintAuxVars(general_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,option)
  ! 
  ! Prints out the contents of an auxvar
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      print *, '     Thermodynamic state: Liquid phase'
      liquid_density = general_auxvar%den(lid)
      liquid_energy = general_auxvar%U(lid)
      liquid_saturation = general_auxvar%sat(lid)
    case(GAS_STATE)
      print *, '     Thermodynamic state: Gas phase'
      gas_density = general_auxvar%den(gid)
      gas_energy = general_auxvar%U(gid)
      gas_saturation = general_auxvar%sat(gid)
    case(TWO_PHASE_STATE)
      print *, '     Thermodynamic state: Two phase'
      liquid_density = general_auxvar%den(lid)
      gas_density = general_auxvar%den(gid)
      liquid_energy = general_auxvar%U(lid)
      gas_energy = general_auxvar%U(gid)
      liquid_saturation = general_auxvar%sat(lid)
      gas_saturation = general_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*general_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*general_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 general_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*general_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*general_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              general_auxvar%effective_porosity*material_auxvar%volume
  print *, 'tot liq comp mass [kmol]: ', liquid_mass
  print *, 'tot gas comp mass [kmol]: ', gas_mass
  print *, '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                         gas_mass*gas_energy
  print *, '         liquid pressure: ', general_auxvar%pres(lid)
  print *, '            gas pressure: ', general_auxvar%pres(gid)
  print *, '            air pressure: ', general_auxvar%pres(apid)
  print *, '      capillary pressure: ', general_auxvar%pres(cpid)
  print *, '          vapor pressure: ', general_auxvar%pres(vpid)
  print *, '     saturation pressure: ', general_auxvar%pres(spid)
  print *, '       liquid saturation: ', general_auxvar%sat(lid)
  print *, '          gas saturation: ', general_auxvar%sat(gid)
  print *, '   liquid density [kmol]: ', general_auxvar%den(lid)
  print *, '      gas density [kmol]: ', general_auxvar%den(gid)
  print *, '     liquid density [kg]: ', general_auxvar%den_kg(lid)
  print *, '        gas density [kg]: ', general_auxvar%den_kg(gid)
  print *, '         temperature [C]: ', general_auxvar%temp
  print *, '      liquid H [MJ/kmol]: ', general_auxvar%H(lid)
  print *, '         gas H [MJ/kmol]: ', general_auxvar%H(gid)
  print *, '      liquid U [MJ/kmol]: ', general_auxvar%U(lid)
  print *, '         gas U [MJ/kmol]: ', general_auxvar%U(gid)
  print *, '     X (water in liquid): ', general_auxvar%xmol(lid,lid)
  print *, '       X (air in liquid): ', general_auxvar%xmol(gid,lid)
  print *, '        X (water in gas): ', general_auxvar%xmol(lid,gid)
  print *, '          X (air in gas): ', general_auxvar%xmol(gid,gid)
  print *, '         liquid mobility: ', general_auxvar%mobility(lid)
  print *, '            gas mobility: ', general_auxvar%mobility(gid)
  print *, '      effective porosity: ', general_auxvar%effective_porosity
  print *, '--------------------------------------------------------'

end subroutine GeneralPrintAuxVars

! ************************************************************************** !

subroutine GeneralOutputAuxVars1(general_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=86,file=string2,position='append')
  else
    open(unit=86,file=string2)
  endif

  write(86,*) '--------------------------------------------------------'
  write(86,*) trim(string)
  write(86,*) '             cell id: ', natural_id
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      write(86,*) ' Thermodynamic state: Liquid phase'
      liquid_density = general_auxvar%den(lid)
      liquid_energy = general_auxvar%U(lid)
      liquid_saturation = general_auxvar%sat(lid)
    case(GAS_STATE)
      write(86,*) ' Thermodynamic state: Gas phase'
      gas_density = general_auxvar%den(gid)
      gas_energy = general_auxvar%U(gid)
      gas_saturation = general_auxvar%sat(gid)
    case(TWO_PHASE_STATE)
      write(86,*) ' Thermodynamic state: Two phase'
      liquid_density = general_auxvar%den(lid)
      gas_density = general_auxvar%den(gid)
      liquid_energy = general_auxvar%U(lid)
      gas_energy = general_auxvar%U(gid)
      liquid_saturation = general_auxvar%sat(lid)
      gas_saturation = general_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*general_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*general_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 general_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*general_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*general_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              general_auxvar%effective_porosity*material_auxvar%volume
  write(86,*) 'tot liq comp mass [kmol]: ', liquid_mass
  write(86,*) 'tot gas comp mass [kmol]: ', gas_mass
  write(86,*) '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                            gas_mass*gas_energy
  write(86,*) '         liquid pressure: ', general_auxvar%pres(lid)
  write(86,*) '            gas pressure: ', general_auxvar%pres(gid)
  write(86,*) '            air pressure: ', general_auxvar%pres(apid)
  write(86,*) '      capillary pressure: ', general_auxvar%pres(cpid)
  write(86,*) '          vapor pressure: ', general_auxvar%pres(vpid)
  write(86,*) '     saturation pressure: ', general_auxvar%pres(spid)
  write(86,*) '         temperature [C]: ', general_auxvar%temp
  write(86,*) '       liquid saturation: ', general_auxvar%sat(lid)
  write(86,*) '          gas saturation: ', general_auxvar%sat(gid)
  write(86,*) '   liquid density [kmol]: ', general_auxvar%den(lid)
  write(86,*) '     liquid density [kg]: ', general_auxvar%den_kg(lid)
  write(86,*) '      gas density [kmol]: ', general_auxvar%den(gid)
  write(86,*) '        gas density [kg]: ', general_auxvar%den_kg(gid)
  write(86,*) '     X (water in liquid): ', general_auxvar%xmol(lid,lid)
  write(86,*) '       X (air in liquid): ', general_auxvar%xmol(gid,lid)
  write(86,*) '        X (water in gas): ', general_auxvar%xmol(lid,gid)
  write(86,*) '          X (air in gas): ', general_auxvar%xmol(gid,gid)
  write(86,*) '      liquid H [MJ/kmol]: ', general_auxvar%H(lid)
  write(86,*) '         gas H [MJ/kmol]: ', general_auxvar%H(gid)
  write(86,*) '      liquid U [MJ/kmol]: ', general_auxvar%U(lid)
  write(86,*) '         gas U [MJ/kmol]: ', general_auxvar%U(gid)
  write(86,*) '         liquid mobility: ', general_auxvar%mobility(lid)
  write(86,*) '            gas mobility: ', general_auxvar%mobility(gid)
  write(86,*) '      effective porosity: ', general_auxvar%effective_porosity
  write(86,*) '...'
  write(86,*) liquid_mass
  write(86,*) gas_mass
  write(86,*) liquid_mass*general_auxvar%U(lid) + &
              gas_mass*general_auxvar%U(gid)
  write(86,*) general_auxvar%pres(lid)
  write(86,*) general_auxvar%pres(gid)
  write(86,*) general_auxvar%pres(apid)
  write(86,*) general_auxvar%pres(cpid)
  write(86,*) general_auxvar%pres(vpid)
  write(86,*) general_auxvar%pres(spid)
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
  write(86,*) general_auxvar%mobility(lid)
  write(86,*) general_auxvar%mobility(gid)
  write(86,*) general_auxvar%effective_porosity
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
  write(86,100) '      liquid pressure: ', &
    ((general_auxvars(idof,i)%pres(lid),i=1,n),idof=0,3)
  write(86,100) '         gas pressure: ', &
    ((general_auxvars(idof,i)%pres(gid),i=1,n),idof=0,3)
  write(86,100) '         air pressure: ', &
    ((general_auxvars(idof,i)%pres(apid),i=1,n),idof=0,3)
  write(86,100) '   capillary pressure: ', &
    ((general_auxvars(idof,i)%pres(cpid),i=1,n),idof=0,3)
  write(86,100) '       vapor pressure: ', &
    ((general_auxvars(idof,i)%pres(vpid),i=1,n),idof=0,3)
  write(86,100) '      temperature [C]: ', &
    ((general_auxvars(idof,i)%temp,i=1,n),idof=0,3)
  write(86,100) '    liquid saturation: ', &
    ((general_auxvars(idof,i)%sat(lid),i=1,n),idof=0,3)
  write(86,100) '       gas saturation: ', &
    ((general_auxvars(idof,i)%sat(gid),i=1,n),idof=0,3)
  write(86,100) 'liquid density [kmol]: ', &
    ((general_auxvars(idof,i)%den(lid),i=1,n),idof=0,3)
  write(86,100) '  liquid density [kg]: ', &
    ((general_auxvars(idof,i)%den_kg(lid),i=1,n),idof=0,3)
  write(86,100) '   gas density [kmol]: ', &
    ((general_auxvars(idof,i)%den(gid),i=1,n),idof=0,3)
  write(86,100) '     gas density [kg]: ', &
    ((general_auxvars(idof,i)%den_kg(gid),i=1,n),idof=0,3)
  write(86,100) '  X (water in liquid): ', &
    ((general_auxvars(idof,i)%xmol(lid,lid),i=1,n),idof=0,3)
  write(86,100) '    X (air in liquid): ', &
    ((general_auxvars(idof,i)%xmol(gid,lid),i=1,n),idof=0,3)
  write(86,100) '     X (water in gas): ', &
    ((general_auxvars(idof,i)%xmol(lid,gid),i=1,n),idof=0,3)
  write(86,100) '       X (air in gas): ', &
    ((general_auxvars(idof,i)%xmol(gid,gid),i=1,n),idof=0,3)
  write(86,100) '   liquid H [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%H(lid),i=1,n),idof=0,3)
  write(86,100) '      gas H [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%H(gid),i=1,n),idof=0,3)
  write(86,100) '   liquid U [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%U(lid),i=1,n),idof=0,3)
  write(86,100) '      gas U [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%U(gid),i=1,n),idof=0,3)
  write(86,*)
  write(86,100) '      liquid mobility: ', &
    ((general_auxvars(idof,i)%mobility(lid),i=1,n),idof=0,3)
  write(86,100) '         gas mobility: ', &
    ((general_auxvars(idof,i)%mobility(gid),i=1,n),idof=0,3)
  write(86,100) '   effective porosity: ', &
    ((general_auxvars(idof,i)%effective_porosity,i=1,n),idof=0,3)
  
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
  use Utility_module, only : DeallocateArray

  implicit none

  type(general_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)  
  call DeallocateArray(auxvar%sat)  
  call DeallocateArray(auxvar%den)  
  call DeallocateArray(auxvar%den_kg)  
  call DeallocateArray(auxvar%xmol)  
  call DeallocateArray(auxvar%H)  
  call DeallocateArray(auxvar%U)  
  call DeallocateArray(auxvar%mobility)  
  
end subroutine GeneralAuxVarStrip

! ************************************************************************** !

subroutine GeneralAuxDestroy(aux)
  ! 
  ! Deallocates a general auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(general_type), pointer :: aux
  PetscInt :: iaux, idof
  
  if (.not.associated(aux)) return
  
  call GeneralAuxVarDestroy(aux%auxvars)
  call GeneralAuxVarDestroy(aux%auxvars_bc)
  call GeneralAuxVarDestroy(aux%auxvars_ss)

  call DeallocateArray(aux%inactive_rows_local)
  call DeallocateArray(aux%inactive_rows_local_ghosted)
  call DeallocateArray(aux%row_zeroing_array)

  if (associated(aux%general_parameter)) then
    call DeallocateArray(aux%general_parameter%diffusion_coefficient)
    deallocate(aux%general_parameter)
  endif
  nullify(aux%general_parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GeneralAuxDestroy

end module General_Aux_module
