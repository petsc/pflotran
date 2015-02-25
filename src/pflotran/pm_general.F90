module PM_General_class

  use PM_Base_class
  use PM_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(pm_subsurface_type) :: pm_general_type
    PetscReal :: dPmax
    PetscReal :: dTmax
    PetscReal :: dXmax
    PetscReal :: dSmax
    PetscReal :: dPmax_allowable
    PetscReal :: dTmax_allowable
    PetscReal :: dXmax_allowable
    PetscReal :: dSmax_allowable
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
  contains
    procedure, public :: SetupSolvers => PMGeneralSetupSolvers
    procedure, public :: InitializeRun => PMGeneralInitializeRun
    procedure, public :: InitializeTimestep => PMGeneralInitializeTimestep
    procedure, public :: Residual => PMGeneralResidual
    procedure, public :: Jacobian => PMGeneralJacobian
    procedure, public :: UpdateTimestep => PMGeneralUpdateTimestep
    procedure, public :: PreSolve => PMGeneralPreSolve
    procedure, public :: PostSolve => PMGeneralPostSolve
    procedure, public :: CheckUpdatePre => PMGeneralCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMGeneralCheckUpdatePost
    procedure, public :: TimeCut => PMGeneralTimeCut
    procedure, public :: UpdateSolution => PMGeneralUpdateSolution
    procedure, public :: UpdateAuxvars => PMGeneralUpdateAuxvars
    procedure, public :: MaxChange => PMGeneralMaxChange
    procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    procedure, public :: Checkpoint => PMGeneralCheckpoint
    procedure, public :: Restart => PMGeneralRestart
    procedure, public :: Destroy => PMGeneralDestroy
  end type pm_general_type
  
  public :: PMGeneralCreate
  
contains

! ************************************************************************** !

function PMGeneralCreate()
  ! 
  ! Creates General process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_general_type), pointer :: PMGeneralCreate

  class(pm_general_type), pointer :: general_pm
  
#ifdef PM_GENERAL_DEBUG  
  print *, 'PMGeneralCreate()'
#endif  

  allocate(general_pm)

  general_pm%dPmax = 0.d0
  general_pm%dTmax = 0.d0
  general_pm%dXmax = 0.d0
  general_pm%dSmax = 0.d0
  general_pm%dPmax_allowable = 5.d5
  general_pm%dTmax_allowable = 5.d0
  general_pm%dXmax_allowable = 0.5d0
  general_pm%dSmax_allowable = 1.d0
  allocate(general_pm%max_change_ivar(6))
  general_pm%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                                LIQUID_MOLE_FRACTION, TEMPERATURE, &
                                GAS_SATURATION]
  allocate(general_pm%max_change_isubvar(6))
                                   ! UNINITIALIZED_INTEGER avoids zeroing of 
                                   ! pressures not represented in phase
                                       ! 2 = air in xmol(air,liquid)
  general_pm%max_change_isubvar = [0,0,0,2,0,0]
  
  call PMSubsurfaceCreate(general_pm)
  general_pm%name = 'PMGeneral'

  PMGeneralCreate => general_pm
  
end function PMGeneralCreate

! ************************************************************************** !

subroutine PMGeneralSetupSolvers(this,solver)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/14

  use General_module, only : GeneralCheckUpdatePre, GeneralCheckUpdatePost
  use Solver_module
  
  implicit none
  
  class(pm_general_type) :: this
  type(solver_type) :: solver
  
  SNESLineSearch :: linesearch
  PetscErrorCode :: ierr
  
  call PMSubsurfaceSetupSolvers(this,solver)

  call SNESGetLineSearch(solver%snes, linesearch, ierr);CHKERRQ(ierr)
  call SNESLineSearchSetPreCheck(linesearch, &
                                 GeneralCheckUpdatePre, &
                                 this%realization,ierr);CHKERRQ(ierr)
  if (solver%check_post_convergence) then
    call SNESLineSearchSetPostCheck(linesearch, &
                                    GeneralCheckUpdatePost, &
                                    this%realization,ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMGeneralSetupSolvers

! ************************************************************************** !

recursive subroutine PMGeneralInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 

  use Realization_Base_class
  
  implicit none
  
  class(pm_general_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,SIX_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 6
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo
  
  ! call parent implementation
  call PMSubsurfaceInitializeRun(this)

end subroutine PMGeneralInitializeRun

! ************************************************************************** !

subroutine PMGeneralInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_general_type) :: this

  call PMSubsurfaceInitializeTimestepA(this)                                 
!geh:remove   everywhere                                
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY,0)
                                 
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," GENERAL FLOW ",64("="))')
  endif
  
  call GeneralInitializeTimestep(this%realization)
  call PMSubsurfaceInitializeTimestepB(this)                                 
  
end subroutine PMGeneralInitializeTimestep

! ************************************************************************** !

subroutine PMGeneralPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

end subroutine PMGeneralPreSolve

! ************************************************************************** !

subroutine PMGeneralPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

end subroutine PMGeneralPostSolve

! ************************************************************************** !

subroutine PMGeneralUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ut, ux, us, umin
  PetscReal :: dtt
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%UpdateTimestep()')
#endif
  
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    up = this%dPmax_allowable/(this%dPmax+0.1)
    ut = this%dTmax_allowable/(this%dTmax+1.d-5)
    ux = this%dXmax_allowable/(this%dXmax+1.d-5)
    us = this%dSmax_allowable/(this%dSmax+1.d-5)
    umin = min(up,ut,ux,us)
  endif
  ifac = max(min(num_newton_iterations,size(tfac)),1)
  dtt = fac * dt * (1.d0 + umin)
  dt = min(dtt,tfac(ifac)*dt,dt_max)
  dt = max(dt,dt_min)

end subroutine PMGeneralUpdateTimestep

! ************************************************************************** !

subroutine PMGeneralResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralResidual

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceUpdatePropertiesNI(this)
  call GeneralResidual(snes,xx,r,this%realization,ierr)

end subroutine PMGeneralResidual

! ************************************************************************** !

subroutine PMGeneralJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralJacobian

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call GeneralJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMGeneralJacobian

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralCheckUpdatePre

  implicit none
  
  class(pm_general_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call GeneralCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

end subroutine PMGeneralCheckUpdatePre

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                    P1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralCheckUpdatePost

  implicit none
  
  class(pm_general_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
  call GeneralCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)

end subroutine PMGeneralCheckUpdatePost

! ************************************************************************** !

subroutine PMGeneralTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralTimeCut

  implicit none
  
  class(pm_general_type) :: this
  
  call PMSubsurfaceTimeCut(this)
  call GeneralTimeCut(this%realization)

end subroutine PMGeneralTimeCut

! ************************************************************************** !

subroutine PMGeneralUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralUpdateSolution, &
                             GeneralMapBCAuxvarsToGlobal

  implicit none
  
  class(pm_general_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call GeneralUpdateSolution(this%realization)
  call GeneralMapBCAuxvarsToGlobal(this%realization)

end subroutine PMGeneralUpdateSolution     

! ************************************************************************** !

subroutine PMGeneralUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14
  use General_module, only : GeneralUpdateAuxVars

  implicit none
  
  class(pm_general_type) :: this

  call GeneralUpdateAuxVars(this%realization,PETSC_FALSE)

end subroutine PMGeneralUpdateAuxvars   

! ************************************************************************** !

subroutine PMGeneralMaxChange(this)
  ! 
  ! Not needed given GeneralMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_Base_class
  use Realization_class
  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  use General_Aux_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
                               TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_general_type) :: this
  
  class(realization_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(6)
  PetscReal :: max_change_global(6)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  
  PetscErrorCode :: ierr
  
  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global = 0.d0
  max_change_local = 0.d0
  do i = 1, 6
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo
    max_change_local(i) = max_change
    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(max_change_local,max_change_global,SIX_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4,&
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4, &
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif
  ! ignore air pressure as it jumps during phase change
  this%dPmax = maxval(max_change_global(1:2))
  this%dXmax = max_change_global(4)
  this%dTmax = max_change_global(5)
  this%dSmax = max_change_global(6)
  
end subroutine PMGeneralMaxChange

! ************************************************************************** !

subroutine PMGeneralComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralComputeMassBalance

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call GeneralComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMGeneralComputeMassBalance

! ************************************************************************** !

subroutine PMGeneralCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call GlobalGetAuxVarVecLoc(this%realization, &
                             this%realization%field%iphas_loc, &
                             STATE,ZERO_INTEGER)
  call PMSubsurfaceCheckpoint(this,viewer)
  
end subroutine PMGeneralCheckpoint

! ************************************************************************** !

subroutine PMGeneralRestart(this,viewer)
  ! 
  ! Restarts data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceRestart(this,viewer)
  call GlobalSetAuxVarVecLoc(this%realization, &
                             this%realization%field%iphas_loc, &
                             STATE,ZERO_INTEGER)
  
end subroutine PMGeneralRestart
! ************************************************************************** !

subroutine PMGeneralDestroy(this)
  ! 
  ! Destroys General process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralDestroy

  implicit none
  
  class(pm_general_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  ! preserve this ordering
  call GeneralDestroy(this%realization)
  call PMSubsurfaceDestroy(this)
  
end subroutine PMGeneralDestroy
  
end module PM_General_class
