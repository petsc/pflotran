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
    class(well_spec_base_type), pointer :: spec  !well_spec pointer
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintBase
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

  public :: WellBaseInit

contains

! ************************************************************************** !

subroutine PrintBase(this)

  implicit none

  class(well_base_type) :: this

  write(*,*) "Well Base Printing message"

end subroutine PrintBase

! ************************************************************************** !

subroutine WellBaseInit(this)
  ! 
  ! Initializes variables/objects in base well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  implicit none

  class(well_base_type) :: this

  this%comm=0;                         
  this%group=0;
  this%cntr_rank=0;                       
  this%z_pw_ref = 0.0;
  this%iwconn_ref = -111
  this%well_num_conns = 0;
  this%cntrl_conn = PETSC_FALSE 
  this%cntrl_lcell_id = -999

  nullify(this%w_rank_conn);   
  nullify(this%disp_rank_conn);
  nullify(this%conn_factors);  
  nullify(this%conn_drill_dir); 
  nullify(this%conn_status); 
  nullify(this%w_conn_order); 
  nullify(this%conn_l2w); 
  nullify(this%w_conn_z);
  nullify(this%spec);

end subroutine WellBaseInit

! ************************************************************************** !

end module Well_Base_class


