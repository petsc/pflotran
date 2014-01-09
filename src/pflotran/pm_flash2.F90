module PM_Flash2_class

  use PM_Base_class
  use Realization_class
  use Communicator_Base_module
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(pm_base_type) :: pm_flash2_type
    class(realization_type), pointer :: realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Init => PMFlash2Init
    procedure, public :: PMFlash2SetRealization
    procedure, public :: InitializeRun => PMFlash2InitializeRun
    procedure, public :: FinalizeRun => PMFlash2FinalizeRun
    procedure, public :: InitializeTimestep => PMFlash2InitializeTimestep
    procedure, public :: FinalizeTimestep => PMFlash2FinalizeTimeStep
    procedure, public :: Residual => PMFlash2Residual
    procedure, public :: Jacobian => PMFlash2Jacobian
    procedure, public :: UpdateTimestep => PMFlash2UpdateTimestep
    procedure, public :: PreSolve => PMFlash2PreSolve
    procedure, public :: PostSolve => PMFlash2PostSolve
    procedure, public :: AcceptSolution => PMFlash2AcceptSolution
#if 0
    procedure, public :: CheckUpdatePre => PMFlash2CheckUpdatePre
    procedure, public :: CheckUpdatePost => PMFlash2CheckUpdatePost
#endif
    procedure, public :: TimeCut => PMFlash2TimeCut
    procedure, public :: UpdateSolution => PMFlash2UpdateSolution
    procedure, public :: MaxChange => PMFlash2MaxChange
    procedure, public :: ComputeMassBalance => PMFlash2ComputeMassBalance
    procedure, public :: Checkpoint => PMFlash2Checkpoint    
    procedure, public :: Restart => PMFlash2Restart  
    procedure, public :: Destroy => PMFlash2Destroy
  end type pm_flash2_type
  
  public :: PMFlash2Create
  
contains

! ************************************************************************** !

function PMFlash2Create()
  ! 
  ! Creates Flash2 process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_flash2_type), pointer :: PMFlash2Create

  class(pm_flash2_type), pointer :: flash2_pm
  
#ifdef PM_FLAHS2_DEBUG
  print *, 'PMFlash2Create()'
#endif  

  allocate(flash2_pm)
  nullify(flash2_pm%option)
  nullify(flash2_pm%output_option)
  nullify(flash2_pm%realization)
  nullify(flash2_pm%comm1)

  call PMBaseCreate(flash2_pm)
  flash2_pm%name = 'PMFlash2'

  PMFlash2Create => flash2_pm
  
end function PMFlash2Create

! ************************************************************************** !

subroutine PMFlash2Init(this)
  ! 
  ! Initializes variables associated with Richard
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

#ifndef SIMPLIFY  
  use Discretization_module
  use Structured_Communicator_class
  use Unstructured_Communicator_class
  use Grid_module 
#endif

  implicit none
  
  class(pm_flash2_type) :: this

#ifdef PM_FLAHS2_DEBUG
  call printMsg(this%option,'PMFlash2%Init()')
#endif
  
#ifndef SIMPLIFY  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      this%comm1 => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select
  call this%comm1%SetDM(this%realization%discretization%dm_1dof)
#endif

end subroutine PMFlash2Init

! ************************************************************************** !

subroutine PMFlash2SetRealization(this,realization)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Realization_class
  use Grid_module

  implicit none
  
  class(pm_flash2_type) :: this
  class(realization_type), pointer :: realization

#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%SetRealization()')
#endif
  
  this%realization => realization
  this%realization_base => realization

  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then 
    this%solution_vec = realization%field%flow_xx_faces
    this%residual_vec = realization%field%flow_r_faces
  else
    this%solution_vec = realization%field%flow_xx
    this%residual_vec = realization%field%flow_r
  endif
  
end subroutine PMFlash2SetRealization

! ************************************************************************** !

subroutine PMFlash2InitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2InitializeTimestep
  use Global_module
  
  implicit none
  
  class(pm_flash2_type) :: this

#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%InitializeTimestep()')
#endif

  this%option%flow_dt = this%option%dt

#ifndef SIMPLIFY  
  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%porosity_loc, &
                              this%realization%field%porosity_loc)
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," FLASH2 FLOW ",62("="))')
  endif
  
  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
  endif  
  
  call Flash2InitializeTimestep(this%realization)
  
end subroutine PMFlash2InitializeTimestep

! ************************************************************************** !

subroutine PMFlash2PreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Global_module

  implicit none
  
  class(pm_flash2_type) :: this
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%PreSolve()')
#endif

end subroutine PMFlash2PreSolve

! ************************************************************************** !

subroutine PMFlash2PostSolve(this)
  ! 
  ! PMFlash2UpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Global_module

  implicit none
  
  class(pm_flash2_type) :: this
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%PostSolve()')
#endif
  
end subroutine PMFlash2PostSolve

! ************************************************************************** !

subroutine PMFlash2FinalizeTimestep(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2MaxChange
  use Global_module

  implicit none
  
  class(pm_flash2_type) :: this
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%FinalizeTimestep()')
#endif
  
  if (this%option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
  endif
  
  call Flash2MaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          this%option%dpmax,this%option%dtmpmax,this%option%dcmax, &
          this%option%dsmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
      this%option%dpmax,this%option%dtmpmax,this%option%dcmax, &
      this%option%dsmax
  endif  
  
end subroutine PMFlash2FinalizeTimestep

! ************************************************************************** !

function PMFlash2AcceptSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_flash2_type) :: this
  
  PetscBool :: PMFlash2AcceptSolution
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%AcceptSolution()')
#endif
  ! do nothing
  PMFlash2AcceptSolution = PETSC_TRUE
  
end function PMFlash2AcceptSolution

! ************************************************************************** !

subroutine PMFlash2UpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_flash2_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
  PetscReal :: uc
  PetscReal :: uus
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%UpdateTimestep()')
#endif
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%option%dpmxe/(this%option%dpmax+0.1)
      utmp = this%option%dtmpmxe/(this%option%dtmpmax+1.d-5)
      uc = this%option%dcmxe/(this%option%dcmax+1.d-6)
      uus= this%option%dsmxe/(this%option%dsmax+1.d-6)
      ut = min(up,utmp,uc,uus)
    endif
    dtt = fac * dt * (1.d0 + ut)
  else
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    dt_tfac = tfac(ifac) * dt

    fac = 0.5d0
    up = this%option%dpmxe/(this%option%dpmax+0.1)
    dt_p = fac * dt * (1.d0 + up)

    dtt = min(dt_tfac,dt_p)
  endif
  
  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
      
  dt = dtt
  
end subroutine PMFlash2UpdateTimestep

! ************************************************************************** !

recursive subroutine PMFlash2InitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2UpdateSolution

  implicit none
  
  class(pm_flash2_type) :: this
  
 
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%InitializeRun()')
#endif
  
  ! restart
  ! init to steady state
  
#if 0  
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif  
  call Flash2UpdateSolution(this%realization)
#endif  
    
end subroutine PMFlash2InitializeRun

! ************************************************************************** !

recursive subroutine PMFlash2FinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_flash2_type) :: this
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%FinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMFlash2FinalizeRun

! ************************************************************************** !

subroutine PMFlash2Residual(this,snes,xx,r,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2Residual

  implicit none
  
  class(pm_flash2_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%Residual()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call Flash2ResidualMFDLP(snes,xx,r,this%realization,ierr)
    case default
      call Flash2Residual(snes,xx,r,this%realization,ierr)
  end select

end subroutine PMFlash2Residual

! ************************************************************************** !

subroutine PMFlash2Jacobian(this,snes,xx,A,B,flag,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2Jacobian

  implicit none
  
  class(pm_flash2_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%Jacobian()')
#endif
  
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
!      call Flash2JacobianMFDLP(snes,xx,A,B,flag,this%realization,ierr)
    case default
      call Flash2Jacobian(snes,xx,A,B,flag,this%realization,ierr)
  end select

end subroutine PMFlash2Jacobian
    
#if 0

! ************************************************************************** !

subroutine PMFlash2CheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2CheckUpdatePre

  implicit none
  
  class(pm_flash2_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%CheckUpdatePre()')
#endif
  
#ifndef SIMPLIFY  
  call Flash2CheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
#endif

end subroutine PMFlash2CheckUpdatePre

! ************************************************************************** !

subroutine PMFlash2CheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2CheckUpdatePost

  implicit none
  
  class(pm_flash2_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%CheckUpdatePost()')
#endif
  
#ifndef SIMPLIFY  
  call Flash2CheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)
#endif

end subroutine PMFlash2CheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMFlash2TimeCut(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2TimeCut

  implicit none
  
  class(pm_flash2_type) :: this
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%TimeCut()')
#endif
  
  this%option%flow_dt = this%option%dt

  call Flash2TimeCut(this%realization)

end subroutine PMFlash2TimeCut

! ************************************************************************** !

subroutine PMFlash2UpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2UpdateSolution
  use Condition_module

  implicit none
  
  class(pm_flash2_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%UpdateSolution()')
#endif

  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%realization%flow_conditions, &
                           this%realization%option, &
                           this%realization%option%time)
  ! right now, RealizUpdateAllCouplerAuxVars only updates flow
  call RealizUpdateAllCouplerAuxVars(this%realization,force_update_flag)
  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif  
  ! end from RealizationUpdate()
  call Flash2UpdateSolution(this%realization)

end subroutine PMFlash2UpdateSolution     

! ************************************************************************** !

subroutine PMFlash2MaxChange(this)
  ! 
  ! Not needed given PMFlash2MaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2MaxChange

  implicit none
  
  class(pm_flash2_type) :: this
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%MaxChange()')
#endif

  call Flash2MaxChange(this%realization)

end subroutine PMFlash2MaxChange

! ************************************************************************** !

subroutine PMFlash2ComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  !use Flash2_module, only : Flash2ComputeMassBalance

  implicit none
  
  class(pm_flash2_type) :: this
  PetscReal :: mass_balance_array(:)
  
#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2%ComputeMassBalance()')
#endif

#ifndef SIMPLIFY
  !geh: currently does not include "trapped" mass
  !call Flash2ComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMFlash2ComputeMassBalance

! ************************************************************************** !

subroutine PMFlash2Checkpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Flash2 PM
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Checkpoint_module

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_flash2_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModel(viewer,this%realization) 
  
end subroutine PMFlash2Checkpoint

! ************************************************************************** !

subroutine PMFlash2Restart(this,viewer)
  ! 
  ! Restarts data associated with Flash2 PM
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Checkpoint_module
  use Flash2_module, only : Flash2UpdateAuxVars

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_flash2_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModel(viewer,this%realization)
  call Flash2UpdateAuxVars(this%realization)
  call this%UpdateSolution()
  
end subroutine PMFlash2Restart

! ************************************************************************** !

subroutine PMFlash2Destroy(this)
  ! 
  ! Destroys Flash2 process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2Destroy

  implicit none
  
  class(pm_flash2_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_FLAHS2_DEBUG  
  call printMsg(this%option,'PMFlash2Destroy()')
#endif

#ifndef SIMPLIFY 
  call Flash2Destroy(this%realization)
#endif
  call this%comm1%Destroy()
  
end subroutine PMFlash2Destroy
  
end module PM_Flash2_class
