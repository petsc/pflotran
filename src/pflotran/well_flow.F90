module Well_Flow_class

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use AuxVars_Flow_module
  use Condition_module
  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_base_type) :: well_flow_type
    PetscReal :: pw_ref                        ! [Pa] well pressure at reference elevation
    PetscReal, pointer :: dw_ref(:)            ! dw_ref(iphase) [kg/m3] well fluid density of iphase at reference elevation
    PetscReal, pointer :: q_fld(:)             ! q_fld(iphase)  [m3/s] well fluid flow rates of iphase
    PetscReal, pointer :: mr_fld(:)            ! mr_fld(iphase) [kg/s] well fluid mass rates of iphase
    PetscReal, pointer :: conn_h(:)            ! connection hydrostatic pressure corrections
    class(auxvar_flow_type), pointer :: flow_auxvars(:,:) !pointer to flow auxvars
    type(flow_condition_type), pointer :: flow_condition ! pointer to flow_condition associated with the well
    !PetscReal, pointer :: conn_mobs(:,:)       ! well connection mobilities ! TO REMOVE - computed when needed flight
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintFlow
    procedure, public :: ConnInit => WellFlowConnInit
    procedure, public :: ExplUpdate => FlowExplUpdate
    procedure, public :: VarsExplUpdate => FlowVarsExplUpdate
    procedure, public :: ConnMob => WellFlowConnMob
    procedure, public :: PressRef => FlowPressRef  
    procedure, public :: QPhase => FlowQPhase
    procedure, public :: MRPhase => FlowMRPhase
    procedure, public :: LimitCheck => WellFlowLimitCheck
    !------------------------------------------------
    !procedure, public :: Init => WellAuxVarBaseInit
    !procedure, public :: Read => WellAuxVarBaseRead
    !procedure, public :: WellAuxVarClear => WellAuxVarBaseClear
    !procedure, public :: WellInit => WellBaseInit
    !procedure, public :: UpdateConnFactor
    !procedure, public :: Output
    !procedure  WellConnInit ! init all vars related to well connections
    !procedure  :: InitWellZRefCntrlConn
    !procedure  :: WellConnSort
  end type  well_flow_type

  public :: WellFlowInit, WellFlowConnInit

contains

! ************************************************************************** !

subroutine PrintFlow(this)

  implicit none

  class(well_flow_type) :: this

  write(*,*) "Well Flow Printing message"

end subroutine PrintFlow

! ************************************************************************** !

subroutine WellFlowInit(this,option)
  ! 
  ! Initializes variables/objects in flow well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module

  implicit none

  class(well_flow_type) :: this

  type(option_type) :: option

  this%pw_ref = 0.0d0;

  allocate( this%dw_ref(option%nphase) );
  this%dw_ref = 0.0d0;
  allocate( this%q_fld(option%nphase) );
  this%q_fld = 0.0d0;
  allocate( this%mr_fld(option%nphase) );
  this%mr_fld = 0.0d0;

  nullify(this%flow_auxvars) 
 
end subroutine WellFlowInit

! ************************************************************************** !

subroutine WellFlowConnInit(this,num_connections,option)
  ! 
  ! Allocate and initilize well_base connections arrays
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/20/2016
  !

  use Option_module

  implicit none

  class(well_flow_type) :: this
  PetscInt, intent(in) :: num_connections 
  type(option_type) :: option  

  call WellBaseConnInit(this,num_connections,option);

  nullify(this%conn_h);
  allocate(this%conn_h(num_connections));
  this%conn_h = 0.0d0; 

end subroutine WellFlowConnInit

! ************************************************************************** !

subroutine FlowExplUpdate(this,grid,option)
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

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  PetscBool :: pass

  write(*,"('FlowExpl d11 before = ',e10.4)"), this%flow_auxvars(0,1)%den(1)
  write(*,"('FlowExpl d12 before = ',e10.4)"), this%flow_auxvars(0,1)%den(2) 
  write(*,"('FlowExpl p11 before = ',e10.4)"), this%flow_auxvars(0,1)%pres(1) 
  !write(*,"('FlowExpl t1 before = ',e10.4)"), this%flow_auxvars(0,1)%temp 

  if(this%connection_set%num_connections == 0 ) return

  pass = PETSC_FALSE

  !cntrl_var = this%cntrl_var ! initialise well control variable
  do
    if(pass) exit ! the well limits are satisfied

    call this%VarsExplUpdate(grid,option)

    call this%LimitCheck(pass)
    ! NOW IMPLEMENT CHECK

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

  !print *, "Well => FlowExplUpdate must be extended"
  !stop  

end subroutine FlowExplUpdate

! ************************************************************************** !

subroutine FlowVarsExplUpdate(this,grid,option)

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  print *, "FlowVarsExplUpdate must be extended"
  stop

end subroutine FlowVarsExplUpdate

!*****************************************************************************!
subroutine FlowPressRef(this,grid,phase,option)
  !  
  ! Compute well p_ref given the volumetric rate of one phase 
  ! IPR sign convention: rate > 0 for fluid being produced
  ! Tested for production well only but should work also for injectors
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 25/06/2015
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option
  
  !PetscInt :: phase, cntrl_var, ivar
  !type(connection_set_type), pointer :: connection_set
  !type(mphase_auxvar_type), pointer :: auxvars(:)

  PetscReal :: rate 

  PetscReal :: press_div_loc, press_div
  PetscReal :: conn_loc, conn_tot
  !PetscReal :: p_well_lc, p_well
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob
 
  if(this%connection_set%num_connections > 0 ) then
    
    rate =  this%flow_condition%flow_well%rate%dataset%rarray(1)
 
    ! divisor computation
    press_div_loc = 0.0d0
    do iconn = 1, this%connection_set%num_connections
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      select case(this%spec%cntrl_var)
        case(CNTRL_VAR_VOL_RATE)
          press_div_loc = press_div_loc + this%conn_factors(iconn) * mob  
        case(CNTRL_VAR_MASS_RATE)
          press_div_loc = press_div_loc + this%conn_factors(iconn) * mob * &
                      this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase)
      end select 
    end do

    call MPI_ALLREDUCE(press_div_loc, press_div, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

    ! compute local well connection flux contribution
    conn_loc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! vol_rate > 0 for fluid entering the well  
        select case(this%spec%cntrl_var)
          case(CNTRL_VAR_VOL_RATE)
            conn_loc = conn_loc + this%conn_factors(iconn) * mob * &
                 ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                   this%conn_h(iconn) )
          case(CNTRL_VAR_MASS_RATE)
            conn_loc = conn_loc + this%conn_factors(iconn) * mob * &
                 this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase) * &
                 ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                   this%conn_h(iconn) )
        end select 
      end if
    end do ! end loop on well connections
    ! connection fluxes contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(conn_loc, conn_tot, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if
 
  ! rate can be volume or mass rate depending on the control variable
  this%pw_ref = (conn_tot - rate) / press_div
  

end subroutine FlowPressRef


! ************************************************************************** !

subroutine FlowQPhase(this,grid,phase,option)
  !  
  ! Compute well volumetric rate for one phase given p_ref 
  ! IPR sign convention: vol_rate > 0 for fluid being produced 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 25/06/2015
  !

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option

  PetscReal :: vol_rate_lc, vol_rate 
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob

  vol_rate = 0.0d0
  if(this%connection_set%num_connections > 0 ) then

    vol_rate_lc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! vol_rate > 0 for fluid entering the well  
        vol_rate_lc = vol_rate_lc + this%conn_factors(iconn) * mob * &
                  ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &                         
                    this%pw_ref - this%conn_h(iconn) )
      end if
    end do ! end loop on well connections
    ! vol_rate contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(vol_rate_lc, vol_rate, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if

  this%q_fld(phase) = vol_rate

end subroutine FlowQPhase

!*****************************************************************************!

subroutine FlowMRPhase(this,grid,phase,option)
  !  
  ! Compute well mass rate for one phase given p_ref 
  ! IPR sign convention: mass_rate > 0 for fluid being produced 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 11/07/2015
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option

  !type(connection_set_type), pointer :: connection_set

  PetscReal :: mass_rate_lc, mass_rate 
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob

  mass_rate = 0.0d0
  if(this%connection_set%num_connections > 0 ) then

    mass_rate_lc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! mass_rate > 0 for fluid entering the well  
        mass_rate_lc = mass_rate_lc + this%conn_factors(iconn) * mob * &
                  this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase) * &
                  (this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                   this%pw_ref - this%conn_h(iconn) )
      end if
    end do ! end loop on well connections
    ! mass_rate contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(mass_rate_lc, mass_rate, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if

  this%mr_fld(phase) = mass_rate

end subroutine FlowMRPhase

!*****************************************************************************!

subroutine WellFlowLimitCheck(this,pass)
  ! 
  !
  ! Perform limit check for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !

  implicit none

  class(well_flow_type) :: this
  PetscBool :: pass

  print *, "WellFlowLimitCheck must be extended"
  stop  

end subroutine WellFlowLimitCheck

!*****************************************************************************!

function WellFlowConnMob(this,mobility,iphase)

  implicit none

  class(well_flow_type) :: this 
  PetscInt :: iphase  !not required for this WellWatInjConnMob extension
  PetscReal :: mobility(:)

  PetscReal :: WellFlowConnMob

  print *, "WellFlowConnMob must be extended"
  stop

end function WellFlowConnMob

!*****************************************************************************!

end module Well_Flow_class
