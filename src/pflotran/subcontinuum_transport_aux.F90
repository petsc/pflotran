module Subcontinuum_Transport_Aux_module
  ! 11/01/2010, Jitu: Replicating reactive_transport_aux implementations with
  ! slight modifications for subcontinuum problem. May simply and merge them
  ! together at a later stage if appropriate
  ! 
  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  use Matrix_Block_Aux_module
  use Reactive_Transport_module

  implicit none
  
  private 

#include "definitions.h"
 
  type, public :: subcontinuum_transport_type
    PetscInt :: num_aux, num_aux_bc
    PetscInt, pointer :: zero_rows_local(:)
    PetscInt :: n_zero_rows
    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    type(reactive_transport_param_type), pointer :: rt_parameter
    type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
    type(reactive_transport_auxvar_type), pointer :: aux_vars_bc(:)
  end type subcontinuum_transport_type
 
  type, public :: subcontinuum_transport_type1
    PetscInt :: num_subgrids
    type(subcontinuum_transport_type), pointer :: st_type(:)
  end type subcontinuum_transport_type1
  
  type, public :: subcontinuum_transport_type2
    PetscInt :: num_subcontinuum
    type(subcontinuum_transport_type1), pointer :: st_type1(:)
  end type subcontinuum_transport_type2

  public :: STAuxCreate, STAuxDestroy, &
            STAuxVarInit, STAuxVarCopy, STAuxVarDestroy
            
contains


! ************************************************************************ !
!
! STAuxCreate: Allocate and initialize subcontinuum auxilliary object
! author: Jitendra Kumar
! date: 11/01/2010 
!
! Duplicated from RTAuxCreate
! ************************************************************************ !
function STAuxCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(subcontinuum_transport_type), pointer :: STAuxCreate
  
  type(subcontinuum_transport_type), pointer :: aux

  allocate(aux)  
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE

  allocate(aux%rt_parameter)
  allocate(aux%rt_parameter%diffusion_coefficient(option%nphase))
  aux%rt_parameter%diffusion_coefficient = 0.d0
  aux%rt_parameter%dispersivity = 0.d0
  aux%rt_parameter%ncomp = 0
  aux%rt_parameter%naqcomp = 0
  aux%rt_parameter%nimcomp = 0
  aux%rt_parameter%ncoll = 0
  aux%rt_parameter%ncollcomp = 0
  aux%rt_parameter%offset_aq = 0
  aux%rt_parameter%offset_coll = 0
  aux%rt_parameter%offset_collcomp = 0
  nullify(aux%rt_parameter%pri_spec_to_coll_spec)
  nullify(aux%rt_parameter%coll_spec_to_pri_spec)
#ifdef OS_STATISTICS
  aux%rt_parameter%newton_call_count = 0
  aux%rt_parameter%sum_newton_call_count = 0.d0
  aux%rt_parameter%newton_iterations = 0
  aux%rt_parameter%sum_newton_iterations = 0.d0
  aux%rt_parameter%max_newton_iterations = 0
  aux%rt_parameter%overall_max_newton_iterations = 0
#endif   
  STAuxCreate => aux
  
end function STAuxCreate

! ************************************************************************ !
!
! STAuxVarInit: Initialize subcontinuum auxilliary object
! author: Jitendra Kumar 
! date: 11/01/2010 
!
! Duplicated from RTAuxVarInit
! ************************************************************************ !
subroutine STAuxVarInit(aux_var,reaction,option)

  use Option_module
  use Reaction_Aux_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  type(reaction_type) :: reaction
  type(option_type) :: option  
  
  call RTAuxVarInit(aux_var, reaction, option)

end subroutine STAuxVarInit

! ************************************************************************ !
!
! STAuxVarCopy: Copys an auxilliary object
! author: Jitendra Kumar 
! date: 11/01/2010
! 
! Duplicated from RTAuxVarCopy
! ************************************************************************ !
subroutine STAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option  
  
  call RTAuxVarCopy(aux_var, aux_var2, option)

end subroutine STAuxVarCopy


! ************************************************************************ !
!
! STAuxVarDestroy: Deallocates a reactive transport auxilliary object
! author: Jitendra Kumar 
! date: 11/01/2010
!
! Duplicated from RTAuxVarDestroy
! ************************************************************************ !
subroutine STAuxVarDestroy(aux_var)

  implicit none

  type(reactive_transport_auxvar_type) :: aux_var
  
  call RTAuxVarDestroy(aux_var)
  
end subroutine STAuxVarDestroy

! ************************************************************************ !
!
! STAuxDestroy: Deallocates a reactive transport auxilliary object
! author: Jitendra Kumar 
! date: 11/01/2010 
!
! Duplicated from RTAuxDestroy
! ************************************************************************ !
subroutine STAuxDestroy(aux)

  implicit none

  type(subcontinuum_transport_type), pointer :: aux
  
  call RTAuxDestroy(aux)
    
end subroutine STAuxDestroy

end module Subcontinuum_Transport_Aux_module
