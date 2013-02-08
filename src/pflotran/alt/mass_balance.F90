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
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

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
  
  use TH_module, only : THResidualToMass
  use THC_module, only : THCResidualToMass
  use THMC_module, only : THMCResidualToMass

  
  implicit none
  
#include "finclude/petscdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  call PetscLogEventBegin(logging%event_mass_balance,ierr)
                                
  option => realization%option
  field => realization%field
  
  if (option%nflowdof > 0) then
    call SNESComputeFunction(flow_solver%snes,field%flow_xx, &
                             field%flow_ts_mass_balance,ierr)
    select case(option%iflowmode)
      case(TH_MODE)
        call THResidualToMass(realization)      
        call VecAXPY(field%flow_total_mass_balance, &
                     1.d0,field%flow_ts_mass_balance,ierr)
      case(THC_MODE)
        call THCResidualToMass(realization)      
        call VecAXPY(field%flow_total_mass_balance, &
                     1.d0,field%flow_ts_mass_balance,ierr)
      case(THMC_MODE)
        call THMCResidualToMass(realization)      
        call VecAXPY(field%flow_total_mass_balance, &
                     1.d0,field%flow_ts_mass_balance,ierr)

      case(RICHARDS_MODE,G_MODE)
        call VecAXPY(field%flow_total_mass_balance, &
                     1.d0,field%flow_ts_mass_balance,ierr)
      case(MPH_MODE)
      case(IMS_MODE)
      case(MIS_MODE)
      case(FLASH2_MODE)
    end select
  endif
  
  if (option%ntrandof > 0) then
    call SNESComputeFunction(tran_solver%snes,field%tran_xx, &
                             field%tran_ts_mass_balance,ierr)
    call VecScale(field%tran_ts_mass_balance,option%tran_dt,ierr)                     
    call VecAXPY(field%tran_total_mass_balance, &
                 1.d0,field%tran_ts_mass_balance,ierr)
  endif

  call PetscLogEventEnd(logging%event_mass_balance,ierr)
  
end subroutine MassBalanceUpdate

end module Mass_Balance_module
