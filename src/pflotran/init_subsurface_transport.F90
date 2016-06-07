module Init_Subsurface_Tran_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: InitSubsurfTranSetupRealization, &
            InitSubsurfTranSetupSolvers
  
contains

! ************************************************************************** !

subroutine InitSubsurfTranSetupRealization(realization)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_Subsurface_class
  use Option_module
  
  use Reactive_Transport_module
  use Global_module
  use Condition_Control_module
  use Variables_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  
  option => realization%option
  
  call RTSetup(realization)

  ! initialize densities and saturations
  if (option%nflowdof == 0) then
    call GlobalSetAuxVarScalar(realization,option%reference_pressure, &
                                LIQUID_PRESSURE)
    call GlobalSetAuxVarScalar(realization,option%reference_temperature, &
                                TEMPERATURE)
    call GlobalSetAuxVarScalar(realization,option%reference_saturation, &
                                LIQUID_SATURATION)
    call GlobalSetAuxVarScalar(realization,option%reference_water_density, &
                                LIQUID_DENSITY)
  else
    call GlobalUpdateAuxVars(realization,TIME_T,0.d0)
    call GlobalWeightAuxVars(realization,0.d0)
  endif

  ! initial concentrations must be assigned after densities are set !!!
  call CondControlAssignTranInitCond(realization)
  
end subroutine InitSubsurfTranSetupRealization

! ************************************************************************** !

subroutine InitSubsurfTranSetupSolvers(realization,convergence_context,solver)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_Subsurface_class
  use Option_module
  use Init_Common_module

  use Reactive_Transport_module
  use Secondary_Continuum_module
  use Solver_module
  use Convergence_module
  use Discretization_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscpc.h"
  
  class(realization_subsurface_type) :: realization
  type(convergence_context_type), pointer :: convergence_context
  type(solver_type), pointer :: solver
  
  type(option_type), pointer :: option
  SNESLineSearch :: linesearch
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  option => realization%option
  
  call printMsg(option,"  Beginning setup of TRAN SNES ")
    
  call SolverCreateSNES(solver,option%mycomm)  
  call SNESSetOptionsPrefix(solver%snes, "tran_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)
    
  if (option%transport%reactive_transport_coupling == GLOBAL_IMPLICIT) then
    if (solver%Jpre_mat_type == '') then
      if (solver%J_mat_type /= MATMFFD) then
        solver%Jpre_mat_type = solver%J_mat_type
      else
        solver%Jpre_mat_type = MATBAIJ
      endif
    endif
    call DiscretizationCreateJacobian(realization%discretization,NTRANDOF, &
                                      solver%Jpre_mat_type, &
                                      solver%Jpre,option)
  else
    solver%J_mat_type = MATAIJ
    solver%Jpre_mat_type = MATAIJ

    call DiscretizationCreateJacobian(realization%discretization,ONEDOF, &
                                      solver%Jpre_mat_type, &
                                      solver%Jpre,option)
  endif

  if (solver%J_mat_type /= MATMFFD) then
    solver%J = solver%Jpre
  endif
    
  call MatSetOptionsPrefix(solver%Jpre,"tran_",ierr);CHKERRQ(ierr)
    
  if (solver%use_galerkin_mg) then
    call DiscretizationCreateInterpolation(realization%discretization,NTRANDOF, &
                                            solver%interpolation, &
                                            solver%galerkin_mg_levels_x, &
                                            solver%galerkin_mg_levels_y, &
                                            solver%galerkin_mg_levels_z, &
                                            option)
  endif

  if (option%transport%reactive_transport_coupling == GLOBAL_IMPLICIT) then

    if (solver%J_mat_type == MATMFFD) then
      call MatCreateSNESMF(solver%snes,solver%J, &
                            ierr);CHKERRQ(ierr)
    endif
      
    ! this could be changed in the future if there is a way to ensure that 
    ! the linesearch update does not perturb concentrations negative.
    call SNESGetLineSearch(solver%snes, linesearch, ierr);CHKERRQ(ierr)
    call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC,  &
                                ierr);CHKERRQ(ierr)
      
    if (option%use_mc) then
      call SNESLineSearchSetPostCheck(linesearch, &
                                      SecondaryRTUpdateIterate, &
                                      realization,ierr);CHKERRQ(ierr)
    endif
      
    ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
    if (option%verbosity >= 2) then
      string = '-tran_snes_view'
      call PetscOptionsInsertString(PETSC_NULL_OBJECT, &
                                    string, ierr);CHKERRQ(ierr)
    endif

  endif

  ! ensure setting of SNES options since they set KSP and PC options too
  call SolverSetSNESOptions(solver)

  option%io_buffer = 'Solver: ' // trim(solver%ksp_type)
  call printMsg(option)
  option%io_buffer = 'Preconditioner: ' // trim(solver%pc_type)
  call printMsg(option)

  if (option%transport%reactive_transport_coupling == GLOBAL_IMPLICIT) then

    ! shell for custom convergence test.  The default SNES convergence test  
    ! is call within this function. 
    !TODO(geh): free this convergence context somewhere!
    option%io_buffer = 'DEALLOCATE TRANSPORT CONVERGENCE CONTEXT somewhere!!!'
    convergence_context => ConvergenceContextCreate(solver,option, &
                                                    realization%patch%grid)
    call SNESSetConvergenceTest(solver%snes,ConvergenceTest, &
                                convergence_context, &
                                PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)
  endif
    
  call printMsg(option,"  Finished setting up TRAN SNES ")
  
end subroutine InitSubsurfTranSetupSolvers

end module Init_Subsurface_Tran_module
