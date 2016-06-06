module Well_WaterInjector_class

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_energy_type) :: well_water_injector_type
    ! .................
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintWatInj
    !procedure, public :: ConnInit => WellTOilImsConnInit !should go in parent classes 
    procedure, public  :: PrintOutputHeader => WellWatInjPrintOutputHeader
    procedure, public :: VarsExplUpdate => WellWatInjVarsExplUpdate
    procedure, public :: LimitCheck => WellWatInjLimitCheck
    procedure, public :: ConnDenUpdate => WellWatInjConnDenUpdate
    procedure, public :: ConnMob => WellWatInjConnMob
  end type  well_water_injector_type

  !public :: CreateTOilImsWell

contains

! ************************************************************************** !

subroutine PrintWatInj(this)

  implicit none

  class(well_water_injector_type) :: this

  write(*,*) "Well PrintWatInj Printing message"

end subroutine PrintWatInj

! ************************************************************************** !

subroutine WellWatInjPrintOutputHeader(this,output_option,file_unit)
  ! 
  ! Write header for well_TOilIms output file
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Output_Aux_module

  implicit none

  class(well_water_injector_type) :: this
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: file_unit

  !character(len=MAXWORDLENGTH) :: tunit

  !tunit = trim(output_option%tunit)

  write(*,*) "Well PrintOutputHeaderWellWatInj to be extended"
  

end subroutine WellWatInjPrintOutputHeader


! ************************************************************************** !

subroutine WellWatInjVarsExplUpdate(this,grid,option)
  !
  ! Explicit update of well variable for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/04/2016
  !

  use Grid_module
  use Option_module
  use EOS_Water_module

  class(well_water_injector_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  PetscReal :: enth_src_h2o, dw_h2o_kg,dw_h2o_mol
  PetscInt :: ierr

  !write(*,"('WellWatInjVarsExplUpdate d11 before = ',e10.4)"), this%flow_energy_auxvarsy(0,1)%den(1)
  !write(*,"('WellWatInjVarsExplUpdate d12 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%den(2) 
  !write(*,"('WellWatInjVarsExplUpdate p11 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%pres(1) 
  !write(*,"('WellWatInjVarsExplUpdate t1 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%temp 

  !write(*,"('WellWatInjVarsExplUpdate rate = ',e10.4)"), &
  !        this%flow_condition%flow_well%rate%dataset%rarray(1)
  !write(*,"('WellWatInjVarsExplUpdate temp = ',e10.4)"), &
  !        this%flow_condition%flow_well%temperature%dataset%rarray(1) 

  !write(*,"('WellWatInjVarsExplUpdate press = ',e10.4)"), &
  !        this%flow_condition%flow_well%pressure%dataset%rarray(1) 

  ! note: cntrl_var is used instead of this%cntrl_var because the
  ! well control variable can change during the well press computation
  ! this%cntrl_var is the initial value 

  if(this%connection_set%num_connections == 0 ) return

  this%tw_ref = this%flow_condition%flow_well%temperature%dataset%rarray(1)

  select case(this%spec%cntrl_var) 
    case(CNTRL_VAR_BHP)
      this%pw_ref = this%flow_condition%flow_well%pressure%dataset%rarray(1)
      call this%QPhase(grid,option%liquid_phase,option)      
      call this%MRPhase(grid,option%liquid_phase,option)

    case(CNTRL_VAR_MASS_RATE)      
      call this%PressRef(grid,option%liquid_phase,option)
      this%mr_fld(option%liquid_phase) = &
                  this%flow_condition%flow_well%rate%dataset%rarray(1)
      call this%QPhase(grid,option%liquid_phase,option)
 
    case(CNTRL_VAR_VOL_RATE)
      this%q_fld(option%liquid_phase) = &
                this%flow_condition%flow_well%rate%dataset%rarray(1)
      call this%PressRef(grid,option%liquid_phase,option)
      call this%MRPhase(grid,option%liquid_phase,option)
  end select
  ! should not be required, they are initialised to zero 
  !this%q_fld(option%oil_phase) = 0.0d0
  !this%mr_fld(option%oil_phase) = 0.0d0
  !call EOSWaterDensityEnthalpy(this%tw_ref,this%pw_ref,dw_h2o_kg,dw_h2o_mol, &
  !                             enth_src_h2o,ierr)

  call EOSWaterDensity(this%tw_ref,this%pw_ref, &
                       dw_h2o_kg,dw_h2o_mol,ierr) 
  !call EOSWaterEnthalpy(toil_auxvar%temp,cell_pressure,toil_auxvar%H(lid),ierr)

  !can load modlar density as well - 
  !might be worth computing the enthalpy as well?? Required in the Res comp later..
  this%dw_kg_ref(option%liquid_phase) = dw_h2o_kg
  ! should not be required, this is already initialised to zero
  !this%dw_kg_ref(option%oil_phase) = 0.0d0


end subroutine WellWatInjVarsExplUpdate

! ************************************************************************** !

subroutine WellWatInjLimitCheck(this,pass)
  !
  ! Perform limit check for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !

  implicit none

  class(well_water_injector_type) :: this
  PetscBool :: pass

  PetscInt :: cntrl_var_tmp 
  PetscReal :: press_max

  pass = PETSC_TRUE
  !if no changes cntrl_var mantains its initial value
  cntrl_var_tmp = this%spec%cntrl_var

  select case(this%spec%cntrl_var)
 
    case(CNTRL_VAR_MASS_RATE)
      if(this%spec%lmt_var(LMT_PRESS_MAX)) then
        press_max = this%flow_condition%flow_well%pressure%dataset%rarray(1) 
        if(this%pw_ref > press_max) then
          print *, "water_injector control switch: " // &
                   "MASS_RATE -> BHP, pw_ref/pw_max= ",&
                    & this%pw_ref, press_max
          ! updates after check pw_ref after check
          this%pw_ref = press_max 
          cntrl_var_tmp = CNTRL_VAR_BHP
          pass = PETSC_FALSE
        end if
      end if 

  end select

  ! update control variable
  this%spec%cntrl_var = cntrl_var_tmp

end subroutine WellWatInjLimitCheck

! ************************************************************************** !

subroutine WellWatInjConnDenUpdate(this,grid,ss_fluxes,option)
  !
  ! Compute connection densities for a water injector 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !
  use Grid_module
  use Option_module

  use EOS_Water_module

  implicit none

  class(well_water_injector_type) :: this
  type(grid_type), pointer :: grid !not currently used
  PetscReal :: ss_fluxes(:,:)      !not currently used
  type(option_type) :: option

  PetscReal :: dw_kg_injw, dw_mol_injw   
  PetscInt :: ierr

  if( this%connection_set%num_connections == 0 ) return

  call EOSWaterDensity(this%tw_ref,this%pw_ref, &
                       dw_kg_injw,dw_mol_injw,ierr) 

  ! all well conns are assigned the same density
  ! a more accurate iterative model could be implemented 
  ! (e.g. hydrostatic coupler)

  this%conn_den_kg = dw_kg_injw !conn_densities is an array

end subroutine WellWatInjConnDenUpdate

! ************************************************************************** !

function WellWatInjConnMob(this,mobility,iphase)
  !  
  ! Compute well connection mobility for mphase mode 
  ! For injection equal to total mobility (summ of all phase mobilities) of
  ! the perforated grid block
  !
  ! example of a function common to all injectors 
  ! (should create an injector class)
  ! otherwise use select case a move this to flow
  !
  ! for efficieny - ConnMob can be extended to well_xxx_mode. 
  ! Instead of performing a loop, a simplue sum can be used 
  ! E.g. WellWatInjConnMob = mobility(liquid_phase) + mobility(oil_phase)
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 1/07/2015
  !

  implicit none

  class(well_water_injector_type) :: this
  PetscInt :: iphase  !not required for this WellWatInjConnMob extension
  PetscReal :: mobility(:)

  PetscReal :: WellWatInjConnMob

  PetscInt :: iph

  WellWatInjConnMob = 0.0d0
  do iph=1,size(mobility)
    WellWatInjConnMob = WellWatInjConnMob + mobility(iph)
  end do


end function WellWatInjConnMob
!*****************************************************************************!


end module Well_WaterInjector_class




