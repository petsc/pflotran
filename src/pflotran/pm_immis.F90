module PM_Immis_class

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

  type, public, extends(pm_subsurface_type) :: pm_immis_type
  contains
    procedure, public :: InitializeTimestep => PMImmisInitializeTimestep
    procedure, public :: Residual => PMImmisResidual
    procedure, public :: Jacobian => PMImmisJacobian
    procedure, public :: UpdateTimestep => PMImmisUpdateTimestep
    procedure, public :: PreSolve => PMImmisPreSolve
    procedure, public :: PostSolve => PMImmisPostSolve
#if 0
    procedure, public :: CheckUpdatePre => PMImmisCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMImmisCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMImmisTimeCut
    procedure, public :: UpdateSolution => PMImmisUpdateSolution
    procedure, public :: UpdateAuxvars => PMImmisUpdateAuxvars
    procedure, public :: MaxChange => PMImmisMaxChange
    procedure, public :: ComputeMassBalance => PMImmisComputeMassBalance
    procedure, public :: Destroy => PMImmisDestroy
  end type pm_immis_type
  
  public :: PMImmisCreate
  
contains

! ************************************************************************** !

function PMImmisCreate()
  ! 
  ! Creates Immiscible process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type), pointer :: PMImmisCreate

  class(pm_immis_type), pointer :: immis_pm
  
  allocate(immis_pm)

  call PMSubsurfaceCreate(immis_pm)
  immis_pm%name = 'PMImmis'

  PMImmisCreate => immis_pm
  
end function PMImmisCreate

! ************************************************************************** !

subroutine PMImmisInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisInitializeTimestep
  
  implicit none
  
  class(pm_immis_type) :: this

  call PMSubsurfaceInitializeTimestep(this)         

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," IMMISCIBLE FLOW ",62("="))')
  endif
  
  call ImmisInitializeTimestep(this%realization)
  
end subroutine PMImmisInitializeTimestep

! ************************************************************************** !

subroutine PMImmisPreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
end subroutine PMImmisPreSolve

! ************************************************************************** !

subroutine PMImmisPostSolve(this)
  ! 
  ! PMImmisUpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
end subroutine PMImmisPostSolve

! ************************************************************************** !

subroutine PMImmisUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
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
  
end subroutine PMImmisUpdateTimestep

! ************************************************************************** !

subroutine PMImmisResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisResidual

  implicit none
  
  class(pm_immis_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call ImmisResidual(snes,xx,r,this%realization,ierr)

end subroutine PMImmisResidual

! ************************************************************************** !

subroutine PMImmisJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisJacobian

  implicit none
  
  class(pm_immis_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call ImmisJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMImmisJacobian
    
#if 0

! ************************************************************************** !

subroutine PMImmisCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisCheckUpdatePre

  implicit none
  
  class(pm_immis_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call ImmisCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)

end subroutine PMImmisCheckUpdatePre

! ************************************************************************** !

subroutine PMImmisCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisCheckUpdatePost

  implicit none
  
  class(pm_immis_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
  call ImmisCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)

end subroutine PMImmisCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMImmisTimeCut(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisTimeCut

  implicit none
  
  class(pm_immis_type) :: this
  
  call PMSubsurfaceTimeCut(this)
  call ImmisTimeCut(this%realization)

end subroutine PMImmisTimeCut

! ************************************************************************** !

subroutine PMImmisUpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisUpdateSolution

  implicit none
  
  class(pm_immis_type) :: this
  
  call PMSubsurfaceUpdateSolution(this)
  call ImmisUpdateSolution(this%realization)

end subroutine PMImmisUpdateSolution     

! ************************************************************************** !

subroutine PMImmisUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Immis_module, only : ImmisUpdateAuxVars
    
  implicit none
  
  class(pm_immis_type) :: this

  call ImmisUpdateAuxVars(this%realization)

end subroutine PMImmisUpdateAuxvars   

! ************************************************************************** !

subroutine PMImmisMaxChange(this)
  ! 
  ! Not needed given PMImmisMaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisMaxChange

  implicit none
  
  class(pm_immis_type) :: this
  
  call ImmisMaxChange(this%realization)
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

end subroutine PMImmisMaxChange

! ************************************************************************** !

subroutine PMImmisComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisComputeMassBalance

  implicit none
  
  class(pm_immis_type) :: this
  PetscReal :: mass_balance_array(:)
  
  !geh: currently does not include "trapped" mass
  !call ImmisComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMImmisComputeMassBalance

! ************************************************************************** !

recursive subroutine PMImmisFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMImmisFinalizeRun

! ************************************************************************** !

subroutine PMImmisDestroy(this)
  ! 
  ! Destroys Immiscible process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisDestroy

  implicit none
  
  class(pm_immis_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call ImmisDestroy(this%realization)
  call PMSubsurfaceDestroy(this)
  
end subroutine PMImmisDestroy
  
end module PM_Immis_class
