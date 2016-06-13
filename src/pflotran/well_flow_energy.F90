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
    !procedure, public :: ExplUpdate => FlowEnergyExplUpdate !could move this to flow
    procedure, public :: VarsExplUpdate => FlowEnergyVarsExplUpdate
    !procedure, public :: QPhase => FlowEnergyQPhase
    !procedure, public :: ConnMob => WellFlowEnergyConnMob
    procedure, public :: ExplJDerivative => WellFlowEnergyExplJDerivative
    procedure, public :: AverageTemp => WellFlowEnergyAverageTemp
    !------------------------------------------------------------
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
subroutine WellFlowEnergyExplJDerivative(this,iconn,ghosted_id,isothermal, &
                                         energy_equation_index,option,Jac)
  ! 
  ! Computes the well derivatives terms for the jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 6/06/16
  ! 

  use Option_module
  !use Condition_module

  implicit none

  class(well_flow_energy_type) :: this
  PetscInt :: iconn 
  PetscInt :: ghosted_id
  PetscBool :: isothermal
  PetscInt :: energy_equation_index
  type(option_type) :: option
  PetscReal :: Jac(option%nflowdof,option%nflowdof)

  !type(flow_toil_ims_condition_type), pointer :: src_sink_condition
  !type(toil_ims_auxvar_type) :: toil_auxvar(0:)
  !class(auxvar_toil_ims_type) :: toil_auxvar(0:)
  !type(auxvar_toil_ims_type) :: toil_auxvar(0:)
  !type(global_auxvar_type) :: global_auxvar
  !PetscReal :: scale
  
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  option%iflag = -3

  !call TOilImsSrcSink(option,src_sink_condition,toil_auxvar(ZERO_INTEGER), &
  !                        global_auxvar,dummy_real,scale,Res)

#ifdef WELL_DEBUG
  write(*,"('ExplJDerivative p011 = ',e16.10)"), this%flow_energy_auxvars(0,1)%pres(1)
  write(*,"('ExplJDerivative p111 = ',e16.10)"), this%flow_energy_auxvars(1,1)%pres(1)
  write(*,"('ExplJDerivative p211 = ',e16.10)"), this%flow_energy_auxvars(2,1)%pres(1)
  write(*,"('ExplJDerivative p311 = ',e16.10)"), this%flow_energy_auxvars(3,1)%pres(1) 
#endif


  call this%ExplRes(iconn,dummy_real,isothermal,ghosted_id,ZERO_INTEGER,&
                    option,res)

  ! downgradient derivatives
  do idof = 1, option%nflowdof

    !call TOilImsSrcSink(option,src_sink_condition,toil_auxvar(idof), &
    !                    global_auxvar,dummy_real,scale,res_pert)
    call this%ExplRes(iconn,dummy_real,isothermal,ghosted_id,idof, &
                      option,res_pert)

    do irow = 1, option%nflowdof
      !Jac(irow,idof) = (res_pert(irow)-res(irow))/toil_auxvar(idof)%pert
      Jac(irow,idof) = (res_pert(irow)-res(irow)) / &
                         this%flow_energy_auxvars(idof,ghosted_id)%pert
    enddo !irow
  enddo ! idof
  
  if (isothermal) then
    !Jac(TOIL_IMS_ENERGY_EQUATION_INDEX,:) = 0.d0
    !Jac(:,TOIL_IMS_ENERGY_EQUATION_INDEX) = 0.d0
    Jac(energy_equation_index,:) = 0.d0
    Jac(:,energy_equation_index) = 0.d0
  endif   

    !Jac =  0.0d0

end subroutine WellFlowEnergyExplJDerivative

! ************************************************************************** !

subroutine WellFlowEnergyAverageTemp(this,grid,ss_fluxes,option)
  !
  ! Compute connection densities for producers
  ! to be overwritten for injectors 
  ! Compute well fluid average density from the well segment phase densities 
  ! using grid-blocks/well fluxes as weights 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/12/2016
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid 
  PetscReal :: ss_fluxes(:,:)      
  type(option_type) :: option

  PetscReal :: q_sum
  PetscReal :: q_ph(option%nphase)

  PetscInt :: iconn, local_id, ghosted_id, ierr
  PetscInt :: i_ph

  PetscReal :: q_sum_lc, q_sum_well 
  PetscReal :: temp_q_lc, temp_q_well
  PetscReal :: temp_lc, temp_well
 
  q_sum_lc = 0.0d0
  q_sum_well = 0.0d0
  temp_q_lc = 0.0d0
  temp_q_well = 0.0d0
  temp_lc = 0.0d0
  temp_well = 0.0d0

  do iconn = 1,this%connection_set%num_connections
    local_id = this%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    ! need to change signs because in the PFLOTRAN convention
    ! producing wells have ss_fluxes < 0
    q_sum = 0.0d0
    do i_ph = 1,option%nphase
      q_ph(i_ph) = 0.0d0
      if( ss_fluxes(i_ph,iconn) < 0.0d0) then
        !should not have ss_fluxes > 0.0d0, this reverse flow situation
        ! should be detected in Res computation 
        q_ph(i_ph) = -1.0d0 * ss_fluxes(i_ph,iconn)
      end if
      q_sum = q_sum + q_ph(i_ph)
    end do

    q_sum_lc = q_sum_lc + q_sum     
    temp_q_lc = temp_q_lc + &
       q_sum * this%flow_energy_auxvars(ZERO_INTEGER,ghosted_id)%temp 
    temp_lc = temp_lc + &
       this%flow_energy_auxvars(ZERO_INTEGER,ghosted_id)%temp
  end do

  call MPI_ALLREDUCE(q_sum_lc, q_sum_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)

  call MPI_ALLREDUCE(temp_q_lc, temp_q_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)

  call MPI_ALLREDUCE(temp_lc, temp_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)
  
  if ( q_sum_well > wfloweps ) then
    this%tw_ref = temp_q_well / q_sum_well
  else
    this%tw_ref = temp_well / dble(this%connection_set%num_connections)
  end if


end subroutine WellFlowEnergyAverageTemp


! ************************************************************************** !

subroutine FlowEnergyVarsExplUpdate(this,grid,ss_fluxes,option)

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  PetscReal :: ss_fluxes(:,:)
  type(option_type) :: option

  print *, "FlowEnergyVarsExplUpdate must be extended"
  stop

end subroutine FlowEnergyVarsExplUpdate


!*****************************************************************************!

!function WellFlowEnergyConnMob(this,mobility,iphase)
!
!  implicit none
!
!  class(well_flow_energy_type) :: this
!  PetscInt :: iphase  
!  PetscReal :: mobility(:)
!
!  PetscReal :: WellFlowEnergyConnMob
!
!  print *, "WellFlowEnergyConnMob must be extended"
!  stop
!
!end function WellFlowEnergyConnMob
!*****************************************************************************!


end module Well_FlowEnergy_class


