module Well_FlowEnergy_class

  use PFLOTRAN_Constants_module
  use Well_Flow_class
  use AuxVars_FlowEnergy_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_type) :: well_flow_energy_type
      PetscReal :: tw_ref                        ! [°C] well temperature at reference elevation
      PetscReal, pointer :: ent_ref(:)           ! ent_ref(iphase) MJ/kmol, well fluid enthalpy of iphase
      class(auxvar_flow_energy_type), pointer :: flow_energy_auxvars(:,:)
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintFlowEnergy
    procedure, public :: ExplUpdate => FlowEnergyExplUpdate !could move this to flow
    procedure, public :: VarsExplUpdate => FlowEnergyVarsExplUpdate
    !procedure, public :: QPhase => FlowEnergyQPhase
    procedure, public :: ConnMob => WellFlowEnergyConnMob
    !procedure, public :: Init => WellAuxVarBaseInit
    !procedure, public :: Read => WellAuxVarBaseRead
    !procedure, public :: WellAuxVarClear => WellAuxVarBaseClear
    !procedure, public :: WellInit => WellBaseInit
    !procedure, public :: UpdateConnFactor
    !procedure, public :: Output
    !procedure  WellConnInit ! init all vars related to well connections
    !procedure  :: InitWellZRefCntrlConn
    !procedure  :: WellConnSort
  end type  well_flow_energy_type

  public :: WellFlowEnergyInit

contains

! ************************************************************************** !

subroutine PrintFlowEnergy(this)

  implicit none

  class(well_flow_energy_type) :: this

  write(*,*) "Well FlowEnergy Printing message"

end subroutine PrintFlowEnergy

! ************************************************************************** !

subroutine WellFlowEnergyInit(this,option)
  ! 
  ! Initializes variables/objects in flow and energy well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 5/20/2016
  ! 

  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(option_type) :: option

  this%tw_ref = 0.0d0; 

  allocate( this%ent_ref(option%nphase) );
  this%ent_ref = 0.0d0;

  nullify(this%flow_energy_auxvars);

end subroutine WellFlowEnergyInit

! ************************************************************************** !

subroutine FlowEnergyExplUpdate(this,grid,option)
  ! 
  ! - Update FlowEnergy well vars
  ! - Perform a limit on well checks 
  ! - Update well control variable in case of switch when a limit is reached
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 6/03/2016
  ! 

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  PetscBool :: pass

  write(*,"('FlowEnergyExpl d11 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%den(1)
  write(*,"('FlowEnergyExpl d12 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%den(2) 
  write(*,"('FlowEnergyExpl p11 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%pres(1) 
  write(*,"('FlowEnergyExpl t1 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%temp 

  if(this%connection_set%num_connections == 0 ) return

  pass = PETSC_FALSE

  !cntrl_var = this%cntrl_var ! initialise well control variable
  do
    if(pass) exit ! the well limits are satisfied

    call this%VarsExplUpdate(grid,option)

    ! NOW IMPLEMENT CHECK
    !call this%ExplIPRMphase(grid,cur_connection_set, &
    !                 flow_condition,auxvars, pw_ref,tw_ref, q_liq,q_gas, &
    !                 m_liq, m_gas, dw_ref, cntrl_var,ivar)

    !print *, "pw_ref = ", pw_ref," ivar = ", ivar
    !pass = PETSC_TRUE ! at the moment no checks
    ! at the moment only checks for gas producer 
    !call this%WellMphaseCheck(flow_condition,pw_ref,q_liq,q_gas,m_liq,m_gas, &
    !                          dw_ref,cntrl_var,ivar,pass)

    ! call this%CheckLimits(pass,cntrl_var,pw_ref,q_liq,q_gas)
    ! during the well checks limits, cntrl_var,pw_ref,q_liq,q_gas can change    

    ! end IPR computation
 
    ! well check - is pw admissible? volumtric rates needed for VFPs
    ! if updates might change the well control variable, those repeating 
    ! the previous operations

  !if well check ok, ends IPR iterative computation
  end do

  ! update well variables 
  !if(wellvar_update) then
  !  this%pw_ref = pw_ref
  !  this%tw_ref = tw_ref
  !  this%q_fld(LIQUID_PHASE) = q_liq
  !  this%q_fld(GAS_PHASE) = q_gas
  !  this%mr_fld(LIQUID_PHASE) = m_liq
  !  this%mr_fld(GAS_PHASE) = m_gas
  !  this%dw_ref(LIQUID_PHASE) = dw_ref(LIQUID_PHASE)
  !  this%dw_ref(GAS_PHASE) = dw_ref(GAS_PHASE)
  !  this%cntrl_var = cntrl_var
  !end if



end subroutine FlowEnergyExplUpdate

! ************************************************************************** !

subroutine FlowEnergyVarsExplUpdate(this,grid,option)

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  print *, "FlowEnergyVarsExplUpdate must be extended"
  stop

end subroutine FlowEnergyVarsExplUpdate


!*****************************************************************************!

function WellFlowEnergyConnMob(this,mobility,iphase)

  implicit none

  class(well_flow_energy_type) :: this
  PetscInt :: iphase  
  PetscReal :: mobility(:)

  PetscReal :: WellFlowEnergyConnMob

  print *, "WellFlowEnergyConnMob must be extended"
  stop

end function WellFlowEnergyConnMob
!*****************************************************************************!


end module Well_FlowEnergy_class


