module PM_Miscible_class

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

  type, public, extends(pm_subsurface_type) :: pm_miscible_type
  contains
    procedure, public :: InitializeTimestep => PMMiscibleInitializeTimestep
    procedure, public :: Residual => PMMiscibleResidual
    procedure, public :: Jacobian => PMMiscibleJacobian
    procedure, public :: UpdateTimestep => PMMiscibleUpdateTimestep
    procedure, public :: PreSolve => PMMisciblePreSolve
    procedure, public :: PostSolve => PMMisciblePostSolve
#if 0
    procedure, public :: CheckUpdatePre => PMMiscibleCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMMiscibleCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMMiscibleTimeCut
    procedure, public :: UpdateSolution => PMMiscibleUpdateSolution
    procedure, public :: UpdateAuxvars => PMMiscibleUpdateAuxvars
    procedure, public :: MaxChange => PMMiscibleMaxChange
    procedure, public :: ComputeMassBalance => PMMiscibleComputeMassBalance
    procedure, public :: Destroy => PMMiscibleDestroy
  end type pm_miscible_type
  
  public :: PMMiscibleCreate
  
contains

! ************************************************************************** !

function PMMiscibleCreate()
  ! 
  ! Creates Miscible process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type), pointer :: PMMiscibleCreate

  class(pm_miscible_type), pointer :: miscible_pm
  
#ifdef PM__DEBUG  
  print *, 'PMMiscibleCreate()'
#endif  

  allocate(miscible_pm)

  call PMSubsurfaceCreate(miscible_pm)
  miscible_pm%name = 'PMMiscible'

  PMMiscibleCreate => miscible_pm
  
end function PMMiscibleCreate

! ************************************************************************** !

subroutine PMMiscibleInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleInitializeTimestep
  
  implicit none
  
  class(pm_miscible_type) :: this

  call PMSubsurfaceInitializeTimestep(this)         

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," MISCIBLE FLOW ",62("="))')
  endif
  
  call MiscibleInitializeTimestep(this%realization)
  
end subroutine PMMiscibleInitializeTimestep

! ************************************************************************** !

subroutine PMMisciblePreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
  
end subroutine PMMisciblePreSolve

! ************************************************************************** !

subroutine PMMisciblePostSolve(this)
  ! 
  ! PMMiscibleUpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
  
end subroutine PMMisciblePostSolve

! ************************************************************************** !

subroutine PMMiscibleUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
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
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%UpdateTimestep()')
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
  
end subroutine PMMiscibleUpdateTimestep

! ************************************************************************** !

subroutine PMMiscibleResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleResidual

  implicit none
  
  class(pm_miscible_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call MiscibleResidual(snes,xx,r,this%realization,ierr)

end subroutine PMMiscibleResidual

! ************************************************************************** !

subroutine PMMiscibleJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleJacobian

  implicit none
  
  class(pm_miscible_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call MiscibleJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMMiscibleJacobian
    
#if 0

! ************************************************************************** !

subroutine PMMiscibleCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleCheckUpdatePre

  implicit none
  
  class(pm_miscible_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call MiscibleCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

end subroutine PMMiscibleCheckUpdatePre

! ************************************************************************** !

subroutine PMMiscibleCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleCheckUpdatePost

  implicit none
  
  class(pm_miscible_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
  call MiscibleCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)

end subroutine PMMiscibleCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMMiscibleTimeCut(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleTimeCut

  implicit none
  
  class(pm_miscible_type) :: this
  
  call PMSubsurfaceTimeCut(this)
  call MiscibleTimeCut(this%realization)

end subroutine PMMiscibleTimeCut

! ************************************************************************** !

subroutine PMMiscibleUpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleUpdateSolution

  implicit none
  
  class(pm_miscible_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call MiscibleUpdateSolution(this%realization)

end subroutine PMMiscibleUpdateSolution     

! ************************************************************************** !

subroutine PMMiscibleUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Miscible_module, only : MiscibleUpdateAuxVars
  
  implicit none
  
  class(pm_miscible_type) :: this

  call MiscibleUpdateAuxVars(this%realization)

end subroutine PMMiscibleUpdateAuxvars   

! ************************************************************************** !

subroutine PMMiscibleMaxChange(this)
  ! 
  ! Not needed given PMMiscibleMaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleMaxChange

  implicit none
  
  class(pm_miscible_type) :: this
  
  call MiscibleMaxChange(this%realization)
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

end subroutine PMMiscibleMaxChange

! ************************************************************************** !

subroutine PMMiscibleComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleComputeMassBalance

  implicit none
  
  class(pm_miscible_type) :: this
  PetscReal :: mass_balance_array(:)
  
  !geh: currently does not include "trapped" mass
  !call MiscibleComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMMiscibleComputeMassBalance

! ************************************************************************** !

subroutine PMMiscibleDestroy(this)
  ! 
  ! Destroys Miscible process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleDestroy

  implicit none
  
  class(pm_miscible_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call MiscibleDestroy(this%realization)
  call PMSubsurfaceDestroy(this)
  
end subroutine PMMiscibleDestroy
  
end module PM_Miscible_class
