module General_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  ! thermodynamic state of fluid ids
  PetscInt, parameter, public :: LIQUID_STATE = 1
  PetscInt, parameter, public :: GAS_STATE = 2
  PetscInt, parameter, public :: TWO_PHASE_STATE = 3

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

    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(general_parameter_type), pointer :: general_parameter
    type(general_auxvar_type), pointer :: aux_vars(:,:)
    type(general_auxvar_type), pointer :: aux_vars_bc(:)
    type(general_auxvar_type), pointer :: aux_vars_ss(:)
  end type general_type

  interface GeneralAuxVarDestroy
    module procedure GeneralAuxVarSingleDestroy
    module procedure GeneralAuxVarArray1Destroy
    module procedure GeneralAuxVarArray2Destroy
  end interface GeneralAuxVarDestroy
  
  public :: GeneralAuxCreate, GeneralAuxDestroy, &
            GeneralAuxVarCompute, GeneralAuxVarInit, &
            GeneralAuxVarCopy, GeneralAuxVarDestroy, &
            GeneralAuxVarStrip, GeneralAuxVarUpdateState

contains


! ************************************************************************** !
!
! GeneralAuxCreate: Allocate and initialize auxiliary object
! author: Glenn Hammond
! date: 03/07/11
!
! ************************************************************************** !
function GeneralAuxCreate(option)

  use Option_module

  implicit none

  type(option_type) :: option
    
  type(general_type), pointer :: GeneralAuxCreate
  
  type(general_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  allocate(aux%general_parameter)
  allocate(aux%general_parameter%diffusion_coefficient(option%nphase))
  aux%general_parameter%diffusion_coefficient = 1.d-9
  allocate(aux%general_parameter%thermal_conductivity(option%nphase))
  aux%general_parameter%thermal_conductivity = 0.d0

  GeneralAuxCreate => aux
  
end function GeneralAuxCreate

! ************************************************************************** !
!
! GeneralAuxVarInit: Initialize auxiliary object
! author: Glenn Hammond
! date: 03/07/11
!
! ************************************************************************** !
subroutine GeneralAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: aux_var
  type(option_type) :: option

  allocate(aux_var%pres(option%nphase+THREE_INTEGER))
  aux_var%pres = 0.d0
  allocate(aux_var%sat(option%nphase))
  aux_var%sat = 0.d0
  allocate(aux_var%den(option%nphase))
  aux_var%den = 0.d0
  allocate(aux_var%den_kg(option%nphase))
  aux_var%den_kg = 0.d0
  ! keep at 25 C.
  aux_var%temp = 25.d0
  allocate(aux_var%xmol(option%nflowspec,option%nphase))
  aux_var%xmol = 0.d0
  allocate(aux_var%H(option%nphase))
  aux_var%H = 0.d0
  allocate(aux_var%U(option%nphase))
  aux_var%U = 0.d0
  allocate(aux_var%kvr(option%nphase))
  aux_var%kvr = 0.d0
  
  aux_var%pert = 0.d0
  
end subroutine GeneralAuxVarInit

! ************************************************************************** !
!
! GeneralAuxVarCopy: Copies an auxiliary variable
! author: Glenn Hammond
! date: 03/07/11
!
! ************************************************************************** !  
subroutine GeneralAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%pres = aux_var%pres
  aux_var2%temp = aux_var%temp
  aux_var2%sat = aux_var%sat
  aux_var2%den = aux_var%den
  aux_var2%den_kg = aux_var%den_kg
  aux_var2%xmol = aux_var%xmol
  aux_var2%H = aux_var%H
  aux_var2%U = aux_var%U
  aux_var2%kvr = aux_var%kvr
  aux_var2%pert = aux_var%pert

end subroutine GeneralAuxVarCopy
  
! ************************************************************************** !
!
! GeneralAuxVarCompute: Computes auxiliary variables for each grid cell
! author: Glenn Hammond
! date: 03/07/11
!
! ************************************************************************** !
subroutine GeneralAuxVarCompute(x,gen_aux_var, global_aux_var,&
                                saturation_function,por,perm,option)

  use Option_module
  use Global_Aux_module
  use water_eos_module
  use Gas_Eos_module
  use Saturation_Function_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(general_auxvar_type) :: gen_aux_var
  type(global_auxvar_type) :: global_aux_var

  PetscReal :: por, perm
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: den_wat_vap, den_kg_wat_vap, h_wat_vap
  PetscReal :: den_air, h_air
  PetscReal :: den_gp, den_gt, hgp, hgt, dgp, dgt, u
  PetscReal :: krl, visl, dkrl_Se
  PetscReal :: krg, visg, dkrg_Se
  PetscReal :: K_H_tilde, P_sat
  PetscReal :: guess, dummy
  PetscInt :: apid, cpid, vpid
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

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  !geh gen_aux_var%temp = 0.d0
#ifdef DEBUG_GENERAL  
  gen_aux_var%H = -999.d0
  gen_aux_var%U = -999.d0
  gen_aux_var%kvr = -999.d0
  gen_aux_var%pres = -999.d0
  gen_aux_var%sat = -999.d0
  gen_aux_var%den = -999.d0
  gen_aux_var%den_kg = -999.d0
  gen_aux_var%xmol = -999.d0
#else
  gen_aux_var%H = 0.d0
  gen_aux_var%U = 0.d0
  gen_aux_var%kvr = 0.d0
  gen_aux_var%pres = 0.d0
  gen_aux_var%sat = 0.d0
  gen_aux_var%den = 0.d0
  gen_aux_var%den_kg = 0.d0
  gen_aux_var%xmol = 0.d0
#endif  
  
  select case(global_aux_var%istate)
    case(LIQUID_STATE)
      gen_aux_var%pres(lid) = x(GENERAL_LIQUID_PRESSURE_DOF)
      gen_aux_var%xmol(acid,lid) = x(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF)
      gen_aux_var%temp = x(GENERAL_LIQUID_STATE_TEMPERATURE_DOF)

      gen_aux_var%xmol(wid,lid) = 1.d0 - gen_aux_var%xmol(acid,lid)
      gen_aux_var%xmol(:,gid) = 0.d0
      gen_aux_var%sat(lid) = 1.d0
      gen_aux_var%sat(gid) = 0.d0

      call psat(gen_aux_var%temp,P_sat,ierr)
      call Henry_air_noderiv(dummy,gen_aux_var%temp, &
                             P_sat,K_H_tilde)
      gen_aux_var%pres(gid) = gen_aux_var%pres(lid)
      gen_aux_var%pres(apid) = K_H_tilde*gen_aux_var%xmol(acid,lid)
      ! need vpres for liq -> 2ph check
      gen_aux_var%pres(vpid) = gen_aux_var%pres(lid) - gen_aux_var%pres(apid)
      gen_aux_var%pres(cpid) = 0.d0
      
    case(GAS_STATE)
      gen_aux_var%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_aux_var%pres(apid) = x(GENERAL_AIR_PRESSURE_DOF)
      gen_aux_var%temp = x(GENERAL_GAS_STATE_TEMPERATURE_DOF)

      gen_aux_var%sat(lid) = 0.d0
      gen_aux_var%sat(gid) = 1.d0
      gen_aux_var%xmol(acid,gid) = gen_aux_var%pres(apid) / gen_aux_var%pres(gid)
      gen_aux_var%xmol(:,lid) = 0.d0
      gen_aux_var%xmol(wid,gid) = 1.d0 - gen_aux_var%xmol(acid,gid)
      gen_aux_var%pres(vpid) = gen_aux_var%pres(gid) - gen_aux_var%pres(apid)
      gen_aux_var%pres(lid) = 0.d0
      
    case(TWO_PHASE_STATE)
      gen_aux_var%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_aux_var%pres(apid) = x(GENERAL_AIR_PRESSURE_DOF)
      gen_aux_var%sat(gid) = x(GENERAL_GAS_SATURATION_DOF)
      
      gen_aux_var%sat(lid) = 1.d0 - gen_aux_var%sat(gid)
      gen_aux_var%pres(vpid) = gen_aux_var%pres(gid) - gen_aux_var%pres(apid)
      
      P_sat = gen_aux_var%pres(vpid)
      guess = gen_aux_var%temp
      call Tsat(gen_aux_var%temp,P_sat,dummy,guess,ierr)
      
      call SatFuncGetCapillaryPressure(gen_aux_var%pres(cpid), &
                                       gen_aux_var%sat(lid), &
                                       saturation_function,option)      

      gen_aux_var%pres(lid) = gen_aux_var%pres(gid) - gen_aux_var%pres(cpid)

      call Henry_air_noderiv(dummy,gen_aux_var%temp, &
                             P_sat,K_H_tilde)
      gen_aux_var%xmol(acid,lid) = gen_aux_var%pres(apid) / K_H_tilde
      gen_aux_var%xmol(wid,lid) = 1.d0 - gen_aux_var%xmol(acid,lid)
      gen_aux_var%xmol(acid,gid) = gen_aux_var%pres(apid) / gen_aux_var%pres(gid)
      gen_aux_var%xmol(wid,gid) = 1.d0 - gen_aux_var%xmol(acid,gid)

  end select


  if (global_aux_var%istate == LIQUID_STATE .or. &
      global_aux_var%istate == TWO_PHASE_STATE) then
    call wateos_noderiv(gen_aux_var%temp,gen_aux_var%pres(lid), &
                        gen_aux_var%den_kg(lid),gen_aux_var%den(lid), &
                        gen_aux_var%H(lid),option%scale,ierr)

    ! MJ/kmol comp
    gen_aux_var%U(lid) = gen_aux_var%H(lid) - &
                         ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                         (gen_aux_var%pres(lid) / gen_aux_var%den(lid) * &
                          option%scale)
    ! this does not need to be calculated for LIQUID_STATE (=1)                     
    call SatFuncGetRelPermFromSat(gen_aux_var%sat(lid),krl,dkrl_Se, &
                                  saturation_function,lid,PETSC_FALSE,option)
    call visw_noderiv(gen_aux_var%temp,gen_aux_var%pres(lid), &
                      P_sat,visl,ierr)
    gen_aux_var%kvr(lid) = krl/visl
  endif

  if (global_aux_var%istate == GAS_STATE .or. &
      global_aux_var%istate == TWO_PHASE_STATE) then
    call ideal_gaseos_noderiv(gen_aux_var%pres(apid),gen_aux_var%temp, &
                              option%scale,den_air,h_air,u)
    call steameos(gen_aux_var%temp,gen_aux_var%pres(gid), &
                  gen_aux_var%pres(apid),den_kg_wat_vap,den_wat_vap,dgp,dgt, &
                  h_wat_vap,hgp,hgt,option%scale,ierr)      
    
    gen_aux_var%den(gid) = den_wat_vap + den_air
    gen_aux_var%den_kg(gid) = den_kg_wat_vap + den_air*FMWAIR
    ! MJ/kmol
    gen_aux_var%H(gid) = gen_aux_var%xmol(wid,gid)*h_wat_vap + &
                         gen_aux_var%xmol(acid,gid)*h_air
    gen_aux_var%U(gid) = gen_aux_var%H(gid) - &
                         ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                         (gen_aux_var%pres(gid) / gen_aux_var%den(gid) * &
                          option%scale)

    ! this does not need to be calculated for GAS_STATE (=1)
    call SatFuncGetRelPermFromSat(gen_aux_var%sat(gid),krg,dkrg_Se, &
                                  saturation_function,gid,PETSC_FALSE,option)
    call visgas_noderiv(gen_aux_var%temp,gen_aux_var%pres(apid), &
                        gen_aux_var%pres(gid),den_air,visg)
    gen_aux_var%kvr(gid) = krg/visg
  endif

end subroutine GeneralAuxVarCompute


! ************************************************************************** !
!
! GeneralUpdateState: Updates the state and swaps primary variables
! author: Glenn Hammond
! date: 05/25/11
!
! ************************************************************************** !
subroutine GeneralAuxVarUpdateState(x,gen_aux_var,global_aux_var, &
                                    saturation_function,por,perm,ghosted_id, &
                                    option)

  use Option_module
  use Global_Aux_module
  use water_eos_module
  use Gas_Eos_module
  use Saturation_Function_module
  
  implicit none

  type(option_type) :: option
  PetscInt :: ghosted_id
  type(saturation_function_type) :: saturation_function
  type(general_auxvar_type) :: gen_aux_var
  type(global_auxvar_type) :: global_aux_var

  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal :: x(option%nflowdof)
  PetscReal :: por, perm
  PetscInt :: apid, cpid, vpid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: dummy, guess
  PetscReal :: P_sat
  PetscBool :: flag
  PetscErrorCode :: ierr

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  flag = PETSC_FALSE
  
  select case(global_aux_var%istate)
    case(LIQUID_STATE)
      call psat(gen_aux_var%temp,P_sat,ierr)
      if (gen_aux_var%pres(vpid) <= P_sat) then
        global_aux_var%istate = TWO_PHASE_STATE
!geh        x(GENERAL_GAS_PRESSURE_DOF) = gen_aux_var%pres(vpid)
!geh        x(GENERAL_AIR_PRESSURE_DOF) = epsilon
        ! vapor pressure needs to be >= epsilon
        x(GENERAL_GAS_PRESSURE_DOF) = gen_aux_var%pres(apid) + P_sat
        x(GENERAL_AIR_PRESSURE_DOF) = gen_aux_var%pres(apid)
        x(GENERAL_GAS_SATURATION_DOF) = epsilon
        flag = PETSC_TRUE
#ifdef DEBUG_GENERAL
        write(option%io_buffer,'(''Liquid -> 2 Phase at Cell '',i11)') ghosted_id
        call printMsg(option)
#endif        
      endif
    case(GAS_STATE)
      call psat(gen_aux_var%temp,P_sat,ierr)
      if (gen_aux_var%pres(vpid) >= P_sat) then
        global_aux_var%istate = TWO_PHASE_STATE
!geh        x(GENERAL_GAS_PRESSURE_DOF) = gen_aux_var%pres(vpid)
!        x(GENERAL_GAS_PRESSURE_DOF) = gen_aux_var%pres(vpid)+gen_aux_var%pres(apid)
        ! first two primary dep vars do not change
        !x(GENERAL_GAS_PRESSURE_DOF) = gen_aux_var%pres(gid)
        !x(GENERAL_AIR_PRESSURE_DOF) = gen_aux_var%pres(apid)
        x(GENERAL_GAS_SATURATION_DOF) = 1.d0 - epsilon
        flag = PETSC_TRUE
#ifdef DEBUG_GENERAL
        write(option%io_buffer,'(''Gas -> 2 Phase at Cell '',i11)') ghosted_id
        call printMsg(option)
#endif        
      endif
    case(TWO_PHASE_STATE)
      if (gen_aux_var%sat(gid) < 0.d0) then
        ! convert to liquid state
        global_aux_var%istate = LIQUID_STATE
!        x(GENERAL_LIQUID_PRESSURE_DOF) = (1.d0+epsilon)* &
!                                         gen_aux_var%pres(lid)
        x(GENERAL_LIQUID_PRESSURE_DOF) = gen_aux_var%pres(gid)
        x(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF) = &
          gen_aux_var%xmol(acid,lid)
        x(GENERAL_LIQUID_STATE_TEMPERATURE_DOF) = gen_aux_var%temp
        flag = PETSC_TRUE
#ifdef DEBUG_GENERAL
        write(option%io_buffer,'(''2 Phase -> Liquid at Cell '',i11)') ghosted_id
        call printMsg(option)
#endif        
      else if (gen_aux_var%sat(gid) > 1.d0) then
        ! convert to gas state
        global_aux_var%istate = GAS_STATE
        ! first two pdv do not change
        !x(GENERAL_GAS_PRESSURE_DOF) = (1.d0-epsilon)* &
        !                              gen_aux_var%pres(vpid)
        !x(GENERAL_AIR_PRESSURE_DOF) = gen_aux_var%pres(apid)
        x(GENERAL_GAS_STATE_TEMPERATURE_DOF) = gen_aux_var%temp
        flag = PETSC_TRUE
#ifdef DEBUG_GENERAL
        write(option%io_buffer,'(''2 Phase -> Gas at Cell '',i11)') ghosted_id
        call printMsg(option)
#endif        
      endif
  end select
  
  if (flag) then
#ifdef DEBUG_GENERAL
    call GeneralPrintAuxVars(gen_aux_var,global_aux_var,ghosted_id, &
                             'Before Update',option)
#endif
    call GeneralAuxVarCompute(x,gen_aux_var, global_aux_var,&
                              saturation_function,por,perm,option)
#ifdef DEBUG_GENERAL
    call GeneralPrintAuxVars(gen_aux_var,global_aux_var,ghosted_id, &
                             'After Update',option)
#endif
    option%variables_swapped = PETSC_TRUE
  endif

end subroutine GeneralAuxVarUpdateState

! ************************************************************************** !
!
! GeneralPrintAuxVars: Prints out the contents of an auxvar
! author: Glenn Hammond
! date: 02/18/13
!
! ************************************************************************** !
subroutine GeneralPrintAuxVars(general_auxvar,global_auxvar,ghosted_id, &
                               string,option)

  use Global_Aux_module
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: ghosted_id
  character(len=*) :: string
  type(option_type) :: option

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
!
! GeneralAuxVarSingleDestroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine GeneralAuxVarSingleDestroy(aux_var)

  implicit none

  type(general_auxvar_type), pointer :: aux_var
  
  if (associated(aux_var)) then
    call GeneralAuxVarStrip(aux_var)
    deallocate(aux_var)
  endif
  nullify(aux_var)  

end subroutine GeneralAuxVarSingleDestroy
  
! ************************************************************************** !
!
! GeneralAuxVarArray1Destroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine GeneralAuxVarArray1Destroy(aux_vars)

  implicit none

  type(general_auxvar_type), pointer :: aux_vars(:)
  
  PetscInt :: iaux
  
  if (associated(aux_vars)) then
    do iaux = 1, size(aux_vars)
      call GeneralAuxVarStrip(aux_vars(iaux))
    enddo  
    deallocate(aux_vars)
  endif
  nullify(aux_vars)  

end subroutine GeneralAuxVarArray1Destroy

! ************************************************************************** !
!
! GeneralAuxVarArray2Destroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine GeneralAuxVarArray2Destroy(aux_vars)

  implicit none

  type(general_auxvar_type), pointer :: aux_vars(:,:)
  
  PetscInt :: iaux, idof
  
  if (associated(aux_vars)) then
    do iaux = 1, size(aux_vars,2)
      do idof = 1, size(aux_vars,1)
        call GeneralAuxVarStrip(aux_vars(idof-1,iaux))
      enddo
    enddo  
    deallocate(aux_vars)
  endif
  nullify(aux_vars)  

end subroutine GeneralAuxVarArray2Destroy

! ************************************************************************** !
!
! GeneralAuxVarDestroy: Deallocates a general auxiliary object
! author: Glenn Hammond
! date: 03/07/11
!
! ************************************************************************** !
subroutine GeneralAuxVarStrip(aux_var)

  implicit none

  type(general_auxvar_type) :: aux_var
  
  if (associated(aux_var%pres)) deallocate(aux_var%pres)
  nullify(aux_var%pres)
  if (associated(aux_var%sat)) deallocate(aux_var%sat)
  nullify(aux_var%sat)
  if (associated(aux_var%den)) deallocate(aux_var%den)
  nullify(aux_var%den)
  if (associated(aux_var%den_kg)) deallocate(aux_var%den_kg)
  nullify(aux_var%den_kg)
  if (associated(aux_var%xmol)) deallocate(aux_var%xmol)
  nullify(aux_var%xmol)
  if (associated(aux_var%H)) deallocate(aux_var%H)
  nullify(aux_var%H)
  if (associated(aux_var%U)) deallocate(aux_var%U)
  nullify(aux_var%U)
  if (associated(aux_var%kvr)) deallocate(aux_var%kvr)
  nullify(aux_var%kvr)
  
end subroutine GeneralAuxVarStrip

! ************************************************************************** !
!
! GeneralAuxDestroy: Deallocates a general auxiliary object
! author: Glenn Hammond
! date: 03/07/11
!
! ************************************************************************** !
subroutine GeneralAuxDestroy(aux)

  implicit none

  type(general_type), pointer :: aux
  PetscInt :: iaux, idof
  
  if (.not.associated(aux)) return
  
  call GeneralAuxVarDestroy(aux%aux_vars)
  call GeneralAuxVarDestroy(aux%aux_vars_bc)
  call GeneralAuxVarDestroy(aux%aux_vars_ss)

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
