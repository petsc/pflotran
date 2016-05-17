module Well_Base_class

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: well_base_type 
    !PetscInt :: id                            ! well id is not needed - use coupler id
    !character(len=MAXWORDLENGTH) :: name      ! well name is not needed - use coupler name
    PetscMPIInt :: comm                        ! well group 
    PetscMPIInt :: group                       ! well communicator
    PetscMPIInt :: cntr_rank                   ! rank where the controlling connection is located
    PetscInt, pointer  :: w_rank_conn(:)       ! number of well conns in each well rank
    PetscInt, pointer  :: disp_rank_conn(:)    ! conns stride for each well rank
    PetscReal :: z_pw_ref                      ! well elevation where the reference pressure is defined
    PetscInt  :: iwconn_ref                    ! index of reference connection in w_conn_z(:)
    PetscInt :: well_num_conns                 ! number of connection for the entire well
    PetscBool :: cntrl_conn                    ! true if the local well sgment coontains the control connection  
    PetscInt  :: cntrl_lcell_id                ! if(cntrl_conn) = local cell id of control connection, otherwise -999  
    PetscReal, pointer :: conn_factors(:)      ! well connection factors
    PetscInt, pointer  :: conn_drill_dir(:)    ! connection drilling directions
    PetscInt, pointer  :: conn_status(:)       ! connection status (can be open, closed or auto)
    PetscInt, pointer  :: w_conn_order(:)      ! connection order for ascending elevation
    PetscInt, pointer  :: conn_l2w(:)          ! map the list of local well conns to well conns 
    PetscReal, pointer :: w_conn_z(:)          ! all well connection elevations ordered by ascending z
    !well_flow_type
      !PetscReal :: pw_ref                        ! [Pa] well pressure at reference elevation
      !PetscReal, pointer :: dw_ref(:)            ! dw_ref(iphase) [kg/m3] well fluid density of iphase at reference elevation
      !PetscReal, pointer :: q_fld(:)             ! q_fld(iphase)  [m3/s] well fluid flow rates of iphase
      !PetscReal, pointer :: mr_fld(:)            ! mr_fld(iphase) [kg/s] well fluid mass rates of iphase
      !PetscReal, pointer :: conn_h(:)            ! connection hydrostatic pressure corrections
      !PetscReal, pointer :: conn_mobs(:,:)       ! well connection mobilities ! TO REMOVE - computed when needed flight
    !well_flow_type end
    !well_flow_energy_type
      !PetscReal :: tw_ref                        ! [°C] well temperature at reference elevation
      !!PetscReal, pointer :: ent_ref(:)          ! MJ/kmol 
    !well_flow_energy_type
    class(well_spec_base_type), pointer :: spec  !well_spec pointer
  contains  ! add here type-bound procedure 
    !procedure, public :: Init => WellAuxVarBaseInit
    !procedure, public :: Read => WellAuxVarBaseRead
    !procedure, public :: WellAuxVarClear => WellAuxVarBaseClear
    !procedure, public :: WellInit => WellBaseInit
    !procedure, public :: UpdateConnFactor
    !procedure, public :: Output
    !procedure  WellConnInit ! init all vars related to well connections
    !procedure  :: InitWellZRefCntrlConn
    !procedure  :: WellConnSort
  end type  well_base_type


end module Well_Base_class


