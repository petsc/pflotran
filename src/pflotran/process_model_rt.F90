module Process_Model_RT_class

  use Process_Model_Base_class
  use Reactive_Transport_module
  use Realization_class
  
  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(process_model_base_type) :: process_model_rt_type
    class(realization_type), pointer :: realization
  contains
    procedure, public :: Init => PMRTInit
    procedure, public :: InitializeTimeStep => PMRTInitializeTimestep
    procedure, public :: Residual => PMRTResidual
    procedure, public :: Jacobian => PMRTJacobian
    procedure, public :: CheckUpdatePre => PMRTCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRTCheckUpdatePost
    procedure, public :: TimeCut => PMRTTimeCut
    procedure, public :: UpdateSolution => PMRTUpdateSolution
    procedure, public :: MaxChange => PMRTMaxChange
    procedure, public :: ComputeMassBalance => PMRTComputeMassBalance
    procedure, public :: Destroy => PMRTDestroy
  end type process_model_rt_type
  
contains

! ************************************************************************** !
!
! PMRTInit: Initializes variables associated with Richard
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTInit(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
#if 0  
  call printMsg(option,"  Beginning setup of TRAN SNES ")
    
  call SolverCreateSNES(tran_solver,option%mycomm)  
  call SNESSetOptionsPrefix(tran_solver%snes, "tran_",ierr)
  call SolverCheckCommandLine(tran_solver)
      
    if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
    if (tran_solver%Jpre_mat_type == '') then
      if (tran_solver%J_mat_type /= MATMFFD) then
        tran_solver%Jpre_mat_type = tran_solver%J_mat_type
      else
        tran_solver%Jpre_mat_type = MATBAIJ
      endif
    endif
    call DiscretizationCreateJacobian(discretization,NTRANDOF, &
                                      tran_solver%Jpre_mat_type, &
                                      tran_solver%Jpre,option)
  else
    tran_solver%J_mat_type = MATAIJ
    tran_solver%Jpre_mat_type = MATAIJ

    call DiscretizationCreateJacobian(discretization,ONEDOF, &
                                      tran_solver%Jpre_mat_type, &
                                      tran_solver%Jpre,option)
  endif

  if (tran_solver%J_mat_type /= MATMFFD) then
    tran_solver%J = tran_solver%Jpre
  endif
    
  call MatSetOptionsPrefix(tran_solver%Jpre,"tran_",ierr)
    
  if (tran_solver%use_galerkin_mg) then
    call DiscretizationCreateInterpolation(discretization,NTRANDOF, &
                                            tran_solver%interpolation, &
                                            tran_solver%galerkin_mg_levels_x, &
                                            tran_solver%galerkin_mg_levels_y, &
                                            tran_solver%galerkin_mg_levels_z, &
                                            option)
  endif

  if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then

    call SNESSetFunction(tran_solver%snes,field%tran_r,RTResidual,&
                          realization,ierr)

    if (tran_solver%J_mat_type == MATMFFD) then
      call MatCreateSNESMF(tran_solver%snes,tran_solver%J,ierr)
    endif
      
    call SNESSetJacobian(tran_solver%snes,tran_solver%J,tran_solver%Jpre, &
                          RTJacobian,realization,ierr)

    ! this could be changed in the future if there is a way to ensure that the linesearch
    ! update does not perturb concentrations negative.
    call SNESGetSNESLineSearch(tran_solver%snes, linesearch, ierr)
    call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr)
      
    if (option%use_mc) then
      call SNESLineSearchSetPostCheck(linesearch, &
                                      SecondaryRTUpdateIterate, &
                                      realization,ierr)      
    endif
      
    ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
    if (option%verbosity >= 1) then
      string = '-tran_snes_view'
      call PetscOptionsInsertString(string, ierr)
    endif

  endif

  ! ensure setting of SNES options since they set KSP and PC options too
  call SolverSetSNESOptions(tran_solver)

  option%io_buffer = 'Solver: ' // trim(tran_solver%ksp_type)
  call printMsg(option)
  option%io_buffer = 'Preconditioner: ' // trim(tran_solver%pc_type)
  call printMsg(option)

  if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then

    ! shell for custom convergence test.  The default SNES convergence test  
    ! is call within this function. 
    tran_stepper%convergence_context => &
      ConvergenceContextCreate(tran_solver,option,grid)
    call SNESSetConvergenceTest(tran_solver%snes,ConvergenceTest, &
                                tran_stepper%convergence_context, &
                                PETSC_NULL_FUNCTION,ierr) 

    ! this update check must be in place, otherwise reactive transport is likely
    ! to fail
    if (associated(realization%reaction)) then
      if (realization%reaction%check_update) then
        call SNESGetSNESLineSearch(tran_solver%snes, linesearch, ierr)
        call SNESLineSearchSetPreCheck(linesearch,RTCheckUpdate, &
                                        realization,ierr)
      endif
    endif
  endif
    
  call printMsg(option,"  Finished setting up TRAN SNES ")  
#endif  
  
  call RTSetup(this%realization)
  
end subroutine PMRTInit

! ************************************************************************** !
!
! PMRTInitializeTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTInitializeTimestep(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTSetup(this%realization)

end subroutine PMRTInitializeTimestep 

! ************************************************************************** !
!
! PMRTResidual: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTResidual(this,snes,xx,r,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call RTResidual(snes,xx,r,this%realization,ierr)
  
end subroutine PMRTResidual

! ************************************************************************** !
!
! PMRTJacobian: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTJacobian(this,snes,xx,A,B,flag,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
  call RTJacobian(snes,xx,A,B,flag,this%realization,ierr)
  
end subroutine PMRTJacobian
    
! ************************************************************************** !
!
! PMRTCheckUpdatePre: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTCheckUpdatePre(this,line_search,P,dP,changed,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call RTCheckUpdate(line_search,P,dP,changed,this%realization,ierr)
  
end subroutine PMRTCheckUpdatePre
    
! ************************************************************************** !
!
! PMRTCheckUpdatePost: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  implicit none
  
  class(process_model_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
!  call RTCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
!                               P1_changed,this%realization,ierr)

end subroutine PMRTCheckUpdatePost
  
! ************************************************************************** !
!
! PMRTTimeCut: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTTimeCut(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTTimeCut(this%realization)

end subroutine PMRTTimeCut
    
! ************************************************************************** !
!
! PMRTUpdateSolution: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTUpdateSolution(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTUpdateSolution(this%realization)

end subroutine PMRTUpdateSolution     

! ************************************************************************** !
!
! PMRTMaxChange: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTMaxChange(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTMaxChange(this%realization)

end subroutine PMRTMaxChange
    
! ************************************************************************** !
!
! PMRTComputeMassBalance: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTComputeMassBalance(this,mass_balance_array)

  implicit none
  
  class(process_model_rt_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call RTComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMRTComputeMassBalance

! ************************************************************************** !
!
! PMRTDestroy: Destroys RT process model
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTDestroy(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTDestroy(this%realization)
  
end subroutine PMRTDestroy
  
end module Process_Model_RT_class
