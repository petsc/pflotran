module PM_Richards_class

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

  type, public, extends(pm_subsurface_type) :: pm_richards_type
  contains
    procedure, public :: InitializeTimestep => PMRichardsInitializeTimestep
    procedure, public :: Residual => PMRichardsResidual
    procedure, public :: Jacobian => PMRichardsJacobian
    procedure, public :: UpdateTimestep => PMRichardsUpdateTimestep
    procedure, public :: PreSolve => PMRichardsPreSolve
    procedure, public :: PostSolve => PMRichardsPostSolve
    procedure, public :: CheckUpdatePre => PMRichardsCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRichardsCheckUpdatePost
    procedure, public :: TimeCut => PMRichardsTimeCut
    procedure, public :: UpdateSolution => PMRichardsUpdateSolution
    procedure, public :: UpdateAuxvars => PMRichardsUpdateAuxvars
    procedure, public :: MaxChange => PMRichardsMaxChange
    procedure, public :: ComputeMassBalance => PMRichardsComputeMassBalance
    procedure, public :: Destroy => PMRichardsDestroy
  end type pm_richards_type
  
  public :: PMRichardsCreate
  
contains

! ************************************************************************** !

function PMRichardsCreate()
  ! 
  ! Creates Richards process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_richards_type), pointer :: PMRichardsCreate

  class(pm_richards_type), pointer :: richards_pm
  
  allocate(richards_pm)
  call PMSubsurfaceCreate(richards_pm)
  richards_pm%name = 'PMRichards'

  PMRichardsCreate => richards_pm
  
end function PMRichardsCreate

! ************************************************************************** !

subroutine PMRichardsInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsInitializeTimestep
  
  implicit none
  
  class(pm_richards_type) :: this

  call PMSubsurfaceInitializeTimestepA(this)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," RICHARDS FLOW ",63("="))')
  endif
  
  call RichardsInitializeTimestep(this%realization)
  call PMSubsurfaceInitializeTimestepB(this)
  
end subroutine PMRichardsInitializeTimestep

! ************************************************************************** !

subroutine PMRichardsPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none
  
  class(pm_richards_type) :: this

end subroutine PMRichardsPreSolve

! ************************************************************************** !

subroutine PMRichardsPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none
  
  class(pm_richards_type) :: this
  
end subroutine PMRichardsPostSolve

! ************************************************************************** !

subroutine PMRichardsUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_richards_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%option%dpmxe/(this%option%dpmax+0.1)
      ut = up
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
  
end subroutine PMRichardsUpdateTimestep

! ************************************************************************** !

subroutine PMRichardsResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsResidual

  implicit none
  
  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceUpdatePropertiesNI(this)
  call RichardsResidual(snes,xx,r,this%realization,ierr)

end subroutine PMRichardsResidual

! ************************************************************************** !

subroutine PMRichardsJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsJacobian

  implicit none
  
  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call RichardsJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMRichardsJacobian

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsCheckUpdatePre

  implicit none
  
  class(pm_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call RichardsCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

end subroutine PMRichardsCheckUpdatePre

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsCheckUpdatePost

  implicit none
  
  class(pm_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
  call RichardsCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)

end subroutine PMRichardsCheckUpdatePost

! ************************************************************************** !

subroutine PMRichardsTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsTimeCut

  implicit none
  
  class(pm_richards_type) :: this
  
  call PMSubsurfaceTimeCut(this)
  call RichardsTimeCut(this%realization)

end subroutine PMRichardsTimeCut

! ************************************************************************** !

subroutine PMRichardsUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsUpdateSolution, &
                              RichardsUpdateSurfacePress

  implicit none
  
  class(pm_richards_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call RichardsUpdateSolution(this%realization)
  if (this%option%surf_flow_on) &
    call RichardsUpdateSurfacePress(this%realization)

end subroutine PMRichardsUpdateSolution     

! ************************************************************************** !

subroutine PMRichardsUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Richards_module, only : RichardsUpdateAuxVars
  
  implicit none
  
  class(pm_richards_type) :: this

  call RichardsUpdateAuxVars(this%realization)

end subroutine PMRichardsUpdateAuxvars   

! ************************************************************************** !

subroutine PMRichardsMaxChange(this)
  ! 
  ! Not needed given RichardsMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsMaxChange

  implicit none
  
  class(pm_richards_type) :: this
  
  call RichardsMaxChange(this%realization)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4)') this%option%dpmax
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') &
      this%option%dpmax
  endif    

end subroutine PMRichardsMaxChange

! ************************************************************************** !

subroutine PMRichardsComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsComputeMassBalance

  implicit none
  
  class(pm_richards_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call RichardsComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMRichardsComputeMassBalance

! ************************************************************************** !

subroutine PMRichardsDestroy(this)
  ! 
  ! Destroys Richards process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsDestroy

  implicit none
  
  class(pm_richards_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call RichardsDestroy(this%realization)
  call PMSubsurfaceDestroy(this)
  
end subroutine PMRichardsDestroy
  
end module PM_Richards_class
