module PM_Flash2_class

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

  type, public, extends(pm_subsurface_type) :: pm_flash2_type
  contains
    procedure, public :: InitializeTimestep => PMFlash2InitializeTimestep
    procedure, public :: Residual => PMFlash2Residual
    procedure, public :: Jacobian => PMFlash2Jacobian
    procedure, public :: UpdateTimestep => PMFlash2UpdateTimestep
    procedure, public :: PreSolve => PMFlash2PreSolve
    procedure, public :: PostSolve => PMFlash2PostSolve
#if 0
    procedure, public :: CheckUpdatePre => PMFlash2CheckUpdatePre
    procedure, public :: CheckUpdatePost => PMFlash2CheckUpdatePost
#endif
    procedure, public :: TimeCut => PMFlash2TimeCut
    procedure, public :: UpdateSolution => PMFlash2UpdateSolution
    procedure, public :: UpdateAuxvars => PMFlash2UpdateAuxvars
    procedure, public :: MaxChange => PMFlash2MaxChange
    procedure, public :: ComputeMassBalance => PMFlash2ComputeMassBalance
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
  
  allocate(flash2_pm)
  call PMSubsurfaceCreate(flash2_pm)
  flash2_pm%name = 'PMFlash2'

  PMFlash2Create => flash2_pm
  
end function PMFlash2Create

! ************************************************************************** !

subroutine PMFlash2InitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2InitializeTimestep
  
  implicit none
  
  class(pm_flash2_type) :: this

  call PMSubsurfaceInitializeTimestep(this)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," FLASH2 FLOW ",62("="))')
  endif
  
  call Flash2InitializeTimestep(this%realization)
  
end subroutine PMFlash2InitializeTimestep

! ************************************************************************** !

subroutine PMFlash2PreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_flash2_type) :: this
  
end subroutine PMFlash2PreSolve

! ************************************************************************** !

subroutine PMFlash2PostSolve(this)
  ! 
  ! PMFlash2UpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_flash2_type) :: this
  
end subroutine PMFlash2PostSolve

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
  
  call Flash2Residual(snes,xx,r,this%realization,ierr)

end subroutine PMFlash2Residual

! ************************************************************************** !

subroutine PMFlash2Jacobian(this,snes,xx,A,B,ierr)
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
  PetscErrorCode :: ierr
  
  call Flash2Jacobian(snes,xx,A,B,this%realization,ierr)

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
  
  call Flash2CheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

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
  
  call Flash2CheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                             P1_changed,this%realization,ierr)

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
  
  call PMSubsurfaceTimeCut(this)
  call Flash2TimeCut(this%realization)

end subroutine PMFlash2TimeCut

! ************************************************************************** !

subroutine PMFlash2UpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Flash2_module, only : Flash2UpdateSolution

  implicit none
  
  class(pm_flash2_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call Flash2UpdateSolution(this%realization)

end subroutine PMFlash2UpdateSolution     

! ************************************************************************** !

subroutine PMFlash2UpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Flash2_module, only : Flash2UpdateAuxVars

  implicit none
  
  class(pm_flash2_type) :: this

  call Flash2UpdateAuxVars(this%realization)

end subroutine PMFlash2UpdateAuxvars   

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
  
  !geh: currently does not include "trapped" mass
  !call Flash2ComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMFlash2ComputeMassBalance

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

  ! preserve this ordering
  call Flash2Destroy(this%realization)
  call PMSubsurfaceDestroy(this)
  
end subroutine PMFlash2Destroy
  
end module PM_Flash2_class
