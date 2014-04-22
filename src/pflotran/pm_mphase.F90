module PM_Mphase_class

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

  type, public, extends(pm_subsurface_type) :: pm_mphase_type
  contains
    procedure, public :: InitializeTimestep => PMMphaseInitializeTimestep
    procedure, public :: Residual => PMMphaseResidual
    procedure, public :: Jacobian => PMMphaseJacobian
    procedure, public :: UpdateTimestep => PMMphaseUpdateTimestep
    procedure, public :: PreSolve => PMMphasePreSolve
    procedure, public :: PostSolve => PMMphasePostSolve
#if 0
    procedure, public :: CheckUpdatePre => PMMphaseCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMMphaseCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMMphaseTimeCut
    procedure, public :: UpdateSolution => PMMphaseUpdateSolution
    procedure, public :: UpdateAuxvars => PMMphaseUpdateAuxvars
    procedure, public :: MaxChange => PMMphaseMaxChange
    procedure, public :: ComputeMassBalance => PMMphaseComputeMassBalance
    procedure, public :: Destroy => PMMphaseDestroy
  end type pm_mphase_type
  
  public :: PMMphaseCreate
  
contains

! ************************************************************************** !

function PMMphaseCreate()
  ! 
  ! Creates Mphase process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_mphase_type), pointer :: PMMphaseCreate

  class(pm_mphase_type), pointer :: mphase_pm

  allocate(mphase_pm)

  call PMSubsurfaceCreate(mphase_pm)
  mphase_pm%name = 'PMMphase'

  PMMphaseCreate => mphase_pm
  
end function PMMphaseCreate

! ************************************************************************** !

subroutine PMMphaseInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseInitializeTimestep
  
  implicit none
  
  class(pm_mphase_type) :: this

  call PMSubsurfaceInitializeTimestep(this)         

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," MPHASE FLOW ",62("="))')
  endif
  
  call MphaseInitializeTimestep(this%realization)
  
end subroutine PMMphaseInitializeTimestep

! ************************************************************************** !

subroutine PMMphasePreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none
  
  class(pm_mphase_type) :: this

end subroutine PMMphasePreSolve

! ************************************************************************** !

subroutine PMMphasePostSolve(this)
  ! 
  ! PMMphaseUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none
  
  class(pm_mphase_type) :: this
  
end subroutine PMMphasePostSolve

! ************************************************************************** !

subroutine PMMphaseUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_mphase_type) :: this
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
  
end subroutine PMMphaseUpdateTimestep

! ************************************************************************** !

subroutine PMMphaseResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseResidual

  implicit none
  
  class(pm_mphase_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call MphaseResidual(snes,xx,r,this%realization,ierr)

end subroutine PMMphaseResidual

! ************************************************************************** !

subroutine PMMphaseJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseJacobian

  implicit none
  
  class(pm_mphase_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call MphaseJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMMphaseJacobian
    
#if 0

! ************************************************************************** !

subroutine PMMphaseCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseCheckUpdatePre

  implicit none
  
  class(pm_mphase_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call MphaseCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

end subroutine PMMphaseCheckUpdatePre

! ************************************************************************** !

subroutine PMMphaseCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseCheckUpdatePost

  implicit none
  
  class(pm_mphase_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
  call MphaseCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)

end subroutine PMMphaseCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMMphaseTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseTimeCut

  implicit none
  
  class(pm_mphase_type) :: this
  
  call PMSubsurfaceTimeCut(this)
  call MphaseTimeCut(this%realization)

end subroutine PMMphaseTimeCut

! ************************************************************************** !

subroutine PMMphaseUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseUpdateSolution

  implicit none
  
  class(pm_mphase_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call MphaseUpdateSolution(this%realization)

end subroutine PMMphaseUpdateSolution     

! ************************************************************************** !

subroutine PMMphaseUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Mphase_module, only : MphaseUpdateAuxVars
  
  implicit none
  
  class(pm_mphase_type) :: this

  call MphaseUpdateAuxVars(this%realization)

end subroutine PMMphaseUpdateAuxvars   

! ************************************************************************** !

subroutine PMMphaseMaxChange(this)
  ! 
  ! Not needed given MphaseMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseMaxChange

  implicit none
  
  class(pm_mphase_type) :: this
  
  call MphaseMaxChange(this%realization)
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

end subroutine PMMphaseMaxChange

! ************************************************************************** !

subroutine PMMphaseComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseComputeMassBalance

  implicit none
  
  class(pm_mphase_type) :: this
  PetscReal :: mass_balance_array(:)
  
  !geh: currently does not include "trapped" mass
  !call MphaseComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMMphaseComputeMassBalance

! ************************************************************************** !

subroutine PMMphaseDestroy(this)
  ! 
  ! Destroys Mphase process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseDestroy

  implicit none
  
  class(pm_mphase_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call MphaseDestroy(this%realization)
  call PMSubsurfaceDestroy(this)
  
end subroutine PMMphaseDestroy
  
end module PM_Mphase_class
