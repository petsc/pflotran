module General_Derivative_module

  
  use PFLOTRAN_Constants_module
  use Option_module
  use General_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: GeneralDerivative

contains
  
! ************************************************************************** !

subroutine GeneralDerivative(option)

  use Characteristic_Curves_module
  use Option_module

  implicit none

  type(option_type), pointer :: option

  type(general_auxvar_type), pointer :: general_auxvar(:)
  type(global_auxvar_type), pointer :: global_auxvar(:)
  class(material_auxvar_type), pointer :: material_auxvar(:)
  class(characteristic_curves_type), pointer :: characteristic_curves
  class(sat_func_VG_type), pointer :: sf
  class(rpf_Mualem_VG_liq_type), pointer :: rpf_liq
  class(rpf_Mualem_VG_gas_type), pointer :: rpf_gas
  PetscReal :: xx(3), pert(3), xx_pert(3)
  PetscInt :: natural_id = 1
  PetscBool :: analytical_derivative = PETSC_TRUE
  character(len=MAXSTRINGLENGTH) :: strings(3,3)
  PetscInt :: i
  
  strings(1,1) = 'Liquid Pressure'
  strings(2,1) = 'Air Mole Fraction in Liquid'
  strings(3,1) = 'Temperature'
  strings(1,2) = 'Gas Pressure'
  strings(2,2) = 'Air Pressure'
  strings(3,2) = 'Temperature'
  strings(1,3) = 'Gas Pressure'
  strings(2,3) = 'Gas Saturation'
  strings(3,3) = 'Temperature'
  
  call GeneralDerivativeSetFlowMode(option)

  allocate(general_auxvar(0:3))
  allocate(global_auxvar(0:3))
  allocate(material_auxvar(0:3))
  
  do i = 0, 3
    call GeneralAuxVarInit(general_auxvar(i),analytical_derivative,option)
    call GlobalAuxVarInit(global_auxvar(i),option)
    call MaterialAuxVarInit(material_auxvar(i),option)
    material_auxvar(i)%porosity_base = 0.25d0
    material_auxvar(i)%volume = 1.d0
    global_auxvar(i)%istate = TWO_PHASE_STATE
  enddo
  
  characteristic_curves => CharacteristicCurvesCreate()
  sf => SF_VG_Create()
  rpf_liq => RPF_Mualem_VG_Liq_Create()
  rpf_gas => RPF_Mualem_VG_Gas_Create()
  sf%m = 0.5d0
  sf%alpha = 1.d-4
  sf%Sr = 0.d0
  sf%pcmax = 1.d6
  characteristic_curves%saturation_function => sf
  rpf_liq%m = 0.5d0
  rpf_liq%Sr = 0.d0
  characteristic_curves%liq_rel_perm_function => rpf_liq
  rpf_gas%m = 0.5d0
  rpf_gas%Sr = 0.d0
  rpf_gas%Srg = 1.d-40
  characteristic_curves%gas_rel_perm_function => rpf_gas
  
  xx(1) = 1.d6
  xx(2) = 0.5d0
  xx(3) = 30.d0
  option%iflag = GENERAL_UPDATE_FOR_ACCUM
  call GeneralAuxVarCompute(xx,general_auxvar(0),global_auxvar(0), &
                            material_auxvar(0),characteristic_curves, &
                            natural_id,option)  
  call GeneralPrintAuxVars(general_auxvar(0),global_auxvar(0),material_auxvar(0), &
                           natural_id,'',option)
  do i = 1, 3
    option%iflag = GENERAL_UPDATE_FOR_DERIVATIVE
    xx_pert = xx
    pert(i) = 1.d-6 * xx(i)
    xx_pert(i) = xx(i) + pert(i)
    call GeneralAuxVarCompute(xx_pert,general_auxvar(i),global_auxvar(i), &
                              material_auxvar(i),characteristic_curves, &
                              natural_id,option)    
    call GeneralAuxVarDiff(i,general_auxvar(0),global_auxvar(0), &
                           material_auxvar(0), &
                           general_auxvar(i),global_auxvar(i), &
                           material_auxvar(i), &
                           pert(i), &
                           strings(i,global_auxvar(i)%istate), &
                           analytical_derivative,option) 
  enddo
  
  do i = 0, 3
    call GeneralAuxVarStrip(general_auxvar(i))
    call GlobalAuxVarStrip(global_auxvar(i))
    call MaterialAuxVarStrip(material_auxvar(i))
  enddo
  deallocate(general_auxvar)
  deallocate(global_auxvar)
  deallocate(material_auxvar)
  call CharacteristicCurvesDestroy(characteristic_curves)
  
end subroutine GeneralDerivative

! ************************************************************************** !

subroutine GeneralDerivativeSetFlowMode(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  
  option%iflowmode = G_MODE
  option%nphase = 2
  option%liquid_phase = 1  ! liquid_pressure
  option%gas_phase = 2     ! gas_pressure

  option%air_pressure_id = 3
  option%capillary_pressure_id = 4
  option%vapor_pressure_id = 5
  option%saturation_pressure_id = 6

  option%water_id = 1
  option%air_id = 2
  option%energy_id = 3

  option%nflowdof = 3
  option%nflowspec = 2
  option%use_isothermal = PETSC_FALSE
      
end subroutine GeneralDerivativeSetFlowMode

! ************************************************************************** !

subroutine GeneralAuxVarDiff(idof,general_auxvar,global_auxvar, &
                             material_auxvar, &
                             general_auxvar_pert,global_auxvar_pert, &
                             material_auxvar_pert, &
                             pert,string,compare_analytical_derivative, &
                             option)

  use Option_module
  use General_Aux_module
  use Global_Aux_module
  use Material_Aux_class  

  implicit none
  
  type(option_type) :: option
  PetscInt :: idof
  type(general_auxvar_type) :: general_auxvar, general_auxvar_pert
  type(global_auxvar_type) :: global_auxvar, global_auxvar_pert
  class(material_auxvar_type) :: material_auxvar, material_auxvar_pert
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: compare_analytical_derivative
  PetscReal :: pert

  
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation
  PetscReal :: liquid_mass_pert, gas_mass_pert
  PetscReal :: liquid_density_pert, gas_density_pert
  PetscReal :: liquid_energy_pert, gas_energy_pert
  PetscReal :: liquid_saturation_pert, gas_saturation_pert
  
  PetscReal :: dpl 
  PetscReal :: dpg 
  PetscReal :: dpa 
  PetscReal :: dpc 
  PetscReal :: dpv 
  PetscReal :: dps 
  PetscReal :: dsatl
  PetscReal :: dsatg
  PetscReal :: ddenl  
  PetscReal :: ddeng  
  PetscReal :: dUl 
  PetscReal :: dHl  
  PetscReal :: dUg  
  PetscReal :: dHg  
  PetscReal :: dpsat  
  PetscReal :: dmobilityl  
  PetscReal :: dmobilityg  
  
  PetscReal, parameter :: uninitialized_value = 0.d0
  
  dpl = uninitialized_value
  dpg = uninitialized_value
  dpa = uninitialized_value
  dpc = uninitialized_value
  dpv = uninitialized_value
  dps = uninitialized_value
  dsatl = uninitialized_value
  dsatg = uninitialized_value
  ddenl = uninitialized_value
  ddeng = uninitialized_value
  dUl = uninitialized_value
  dHl = uninitialized_value
  dUg = uninitialized_value
  dHg = uninitialized_value
  dpsat = uninitialized_value
  dmobilityl = uninitialized_value
  dmobilityg = uninitialized_value

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
    
  if (compare_analytical_derivative) then
    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        select case(idof)
          case(1)
          case(2)
          case(3)
        end select
      case(GAS_STATE)
        select case(idof)
          case(1)
          case(2)
          case(3)
        end select
      case(TWO_PHASE_STATE)
        select case(idof)
          case(1) ! pg
            dpl = 1.d0 ! pl = pg - pc
            dpg = 1.d0 ! pg = pg
            dpa = 1.d0 ! pa = pg - pv
            dpv = 0.d0
            dps = 0.d0
            dsatl = 0.d0
            dsatg = 0.d0
            ddenl = general_auxvar%d%denl_pl
            ddeng = general_auxvar%d%deng_pg
            dUl = general_auxvar%d%Ul_pl
            dHl = general_auxvar%d%Hl_pl
            dUg = general_auxvar%d%Ug_pg
            dHg = general_auxvar%d%Hg_pg
            dmobilityl = general_auxvar%d%mobilityl_pl
            dmobilityg = general_auxvar%d%mobilityg_pg
          case(2)
            dpl = -1.d0*general_auxvar%d%pc_satg ! pl = pg - pc
            dpg = 1.d0 ! pg = pg
            dpa = 1.d0 ! pa = pg - pv
            dpc = general_auxvar%d%pc_satg
            dpv = 0.d0 ! pv = pg - pa
            dps = 0.d0
            dsatl = -1.d0
            dsatg = 1.d0
            ddenl = 0.d0
            dmobilityl = general_auxvar%d%mobilityl_satg
            dmobilityg = general_auxvar%d%mobilityg_satg
          case(3)
            dpl = 0.d0 ! pl = pg - pc
            dpg = 0.d0 ! pg = pg
            dpa = -1.d0*general_auxvar%d%psat_dT ! pa = pg - pv
            dpv = general_auxvar%d%psat_dT
            dps = general_auxvar%d%psat_dT
            dsatl = 0.d0
            dsatg = 0.d0            
            ddenl = general_auxvar%d%Hg_T
            ddeng = general_auxvar%d%Hg_T
            dUl = general_auxvar%d%Ul_T
            dHl = general_auxvar%d%Hl_T
            dUg = general_auxvar%d%Ug_T
            dHg = general_auxvar%d%Hg_T
            dmobilityl = general_auxvar%d%mobilityl_T
            dmobilityg = general_auxvar%d%mobilityg_T
        end select
      end select
    endif

  print *, '--------------------------------------------------------'
  print *, 'Derivative with respect to ' // trim(string)
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
  select case(global_auxvar_pert%istate)
    case(LIQUID_STATE)
      print *, '     Thermodynamic state (pert): Liquid phase'
      liquid_density_pert = general_auxvar_pert%den(lid)
      liquid_energy_pert = general_auxvar_pert%U(lid)
      liquid_saturation_pert = general_auxvar_pert%sat(lid)
    case(GAS_STATE)
      print *, '     Thermodynamic state (pert): Gas phase'
      gas_density_pert = general_auxvar_pert%den(gid)
      gas_energy_pert = general_auxvar_pert%U(gid)
      gas_saturation_pert = general_auxvar_pert%sat(gid)
    case(TWO_PHASE_STATE)
      print *, '     Thermodynamic state (pert): Two phase'
      liquid_density_pert = general_auxvar_pert%den(lid)
      gas_density_pert = general_auxvar_pert%den(gid)
      liquid_energy_pert = general_auxvar_pert%U(lid)
      gas_energy_pert = general_auxvar_pert%U(gid)
      liquid_saturation_pert = general_auxvar_pert%sat(lid)
      gas_saturation_pert = general_auxvar_pert%sat(gid)
  end select  
  liquid_mass_pert = (liquid_density_pert*general_auxvar_pert%xmol(lid,lid)* & 
                 liquid_saturation_pert+ &
                 gas_density_pert*general_auxvar_pert%xmol(lid,gid)* & 
                 gas_saturation_pert)* & 
                 general_auxvar_pert%effective_porosity*material_auxvar_pert%volume
  gas_mass_pert = (liquid_density_pert*general_auxvar_pert%xmol(gid,lid)* & 
              liquid_saturation_pert+ &
              gas_density_pert*general_auxvar_pert%xmol(gid,gid)* & 
              gas_saturation_pert)* & 
              general_auxvar_pert%effective_porosity*material_auxvar_pert%volume 
100 format(a,100('','',es13.5))  
  write(*,100) 'tot liq comp mass [kmol]: ', (liquid_mass_pert-liquid_mass)/pert
  write(*,100) 'tot gas comp mass [kmol]: ', (gas_mass_pert-gas_mass)/pert
  write(*,100) '             energy [MJ]: ', ((liquid_mass_pert*liquid_energy_pert + &
                                           gas_mass_pert*gas_energy_pert)- &
                                          (liquid_mass*liquid_energy + &
                                           gas_mass*gas_energy))/pert
  write(*,100) '         liquid pressure: ', (general_auxvar_pert%pres(lid)-general_auxvar%pres(lid))/pert,dpl
  write(*,100) '            gas pressure: ', (general_auxvar_pert%pres(gid)-general_auxvar%pres(gid))/pert,dpg
  write(*,100) '            air pressure: ', (general_auxvar_pert%pres(apid)-general_auxvar%pres(apid))/pert,dpa
  write(*,100) '      capillary pressure: ', (general_auxvar_pert%pres(cpid)-general_auxvar%pres(cpid))/pert,dpc
  write(*,100) '          vapor pressure: ', (general_auxvar_pert%pres(vpid)-general_auxvar%pres(vpid))/pert,dpv
  write(*,100) '     saturation pressure: ', (general_auxvar_pert%pres(spid)-general_auxvar%pres(spid))/pert,dps
  write(*,100) '       liquid saturation: ', (general_auxvar_pert%sat(lid)-general_auxvar%sat(lid))/pert
  write(*,100) '          gas saturation: ', (general_auxvar_pert%sat(gid)-general_auxvar%sat(gid))/pert
  write(*,100) '   liquid density [kmol]: ', (general_auxvar_pert%den(lid)-general_auxvar%den(lid))/pert,ddenl
  write(*,100) '      gas density [kmol]: ', (general_auxvar_pert%den(gid)-general_auxvar%den(gid))/pert,ddeng
  write(*,100) '     liquid density [kg]: ', (general_auxvar_pert%den_kg(lid)-general_auxvar%den_kg(lid))/pert
  write(*,100) '        gas density [kg]: ', (general_auxvar_pert%den_kg(gid)-general_auxvar%den_kg(gid))/pert
  write(*,100) '         temperature [C]: ', (general_auxvar_pert%temp-general_auxvar%temp)/pert
  write(*,100) '      liquid H [MJ/kmol]: ', (general_auxvar_pert%H(lid)-general_auxvar%H(lid))/pert,dHl
  write(*,100) '         gas H [MJ/kmol]: ', (general_auxvar_pert%H(gid)-general_auxvar%H(gid))/pert,dHg
  write(*,100) '      liquid U [MJ/kmol]: ', (general_auxvar_pert%U(lid)-general_auxvar%U(lid))/pert,dUl
  write(*,100) '         gas U [MJ/kmol]: ', (general_auxvar_pert%U(gid)-general_auxvar%U(gid))/pert,dUg
  write(*,100) '     X (water in liquid): ', (general_auxvar_pert%xmol(lid,lid)-general_auxvar%xmol(lid,lid))/pert
  write(*,100) '       X (air in liquid): ', (general_auxvar_pert%xmol(gid,lid)-general_auxvar%xmol(gid,lid))/pert
  write(*,100) '        X (water in gas): ', (general_auxvar_pert%xmol(lid,gid)-general_auxvar%xmol(lid,gid))/pert
  write(*,100) '          X (air in gas): ', (general_auxvar_pert%xmol(gid,gid)-general_auxvar%xmol(gid,gid))/pert
  write(*,100) '         liquid mobility: ', (general_auxvar_pert%mobility(lid)-general_auxvar%mobility(lid))/pert,dmobilityl
  write(*,100) '            gas mobility: ', (general_auxvar_pert%mobility(gid)-general_auxvar%mobility(gid))/pert,dmobilityg
  write(*,100) '      effective porosity: ', (general_auxvar_pert%effective_porosity-general_auxvar%effective_porosity)/pert
  write(*,100) '--------------------------------------------------------'  
  
end subroutine GeneralAuxVarDiff

end module General_Derivative_module
