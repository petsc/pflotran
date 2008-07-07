module Mass_Balance_module

  implicit none

  private

#include "definitions.h"
  
  public :: MassBalanceCreate, MassBalanceUpdate
  
contains

! ************************************************************************** !
!
! MassBalanceCreate: 
! author: Glenn Hammond
! date: 06/18/08
!
! ************************************************************************** !
subroutine MassBalanceCreate(realization)

  use Realization_module
  use Option_module
  use Field_module
  
  implicit none
  
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  type(realization_type) :: realization
  type(field_type), pointer :: field
  PetscErrorCode :: ierr
  
  field => realization%field
  
  if (realization%option%nflowdof > 0) then
    call VecDuplicate(field%flow_xx,field%flow_ts_mass_balance,ierr)
    call VecDuplicate(field%flow_xx,field%flow_total_mass_balance,ierr)
  endif
  
  if (realization%option%nflowdof > 0) then
    call VecDuplicate(field%tran_xx,field%tran_ts_mass_balance,ierr)
    call VecDuplicate(field%tran_xx,field%tran_total_mass_balance,ierr)
  endif

end subroutine MassBalanceCreate

! ************************************************************************** !
!
! MassBalanceUpdate: 
! author: Glenn Hammond
! date: 06/18/08
!
! ************************************************************************** !
subroutine MassBalanceUpdate(realization,flow_solver,tran_solver)

  use Realization_module
  use Option_module
  use Field_module
  use Solver_module
  use Logging_module
  
  use Richards_module, only : RichardsResidualToMass
  
  implicit none
  
#include "include/finclude/petscdef.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscsnes.h"

  type(realization_type) :: realization
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  call PetscLogEventBegin(logging%event_mass_balance, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                                
  option => realization%option
  field => realization%field
  
  if (option%nflowdof > 0) then
    call SNESComputeFunction(flow_solver%snes,field%flow_xx, &
                             field%flow_ts_mass_balance,ierr)
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        call RichardsResidualToMass(realization)      
        call VecAXPY(field%tran_total_mass_balance, &
                     field%tran_ts_mass_balance,1.d0,ierr)
      case(RICHARDS_LITE_MODE)
        call VecAXPY(field%tran_total_mass_balance, &
                     field%tran_ts_mass_balance,1.d0,ierr)
      case(MPH_MODE)
    end select
  endif
  
  if (option%ntrandof > 0) then
    call SNESComputeFunction(tran_solver%snes,field%tran_xx, &
                             field%tran_ts_mass_balance,ierr)
    call VecScale(field%tran_ts_mass_balance,option%tran_dt,ierr)                     
    call VecAXPY(field%tran_total_mass_balance, &
                 field%tran_ts_mass_balance,1.d0,ierr)
  endif

  call PetscLogEventEnd(logging%event_mass_balance, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  
end subroutine MassBalanceUpdate

end module Mass_Balance_module
