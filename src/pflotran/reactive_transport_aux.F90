module Reactive_Transport_Aux_module

  implicit none
  
  private 

#include "definitions.h"
  
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petsclog.h"

  type, public :: reactive_transport_auxvar_type
    PetscReal, pointer :: total(:)
    PetscReal, pointer :: dtotal(:,:)
  end type reactive_transport_auxvar_type
  
  type, public :: reactive_transport_aux_type
    PetscInt :: num_aux, num_aux_bc
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscInt :: n_zero_rows
    type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
    type(reactive_transport_auxvar_type), pointer :: aux_vars_bc(:)
  end type reactive_transport_aux_type

  public :: ReactiveTransportAuxCreate, ReactiveTransportAuxDestroy, &
            RTAuxVarCompute, RTAuxVarInit

contains


! ************************************************************************** !
!
! ReactiveTransportAuxCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
function ReactiveTransportAuxCreate()

  use Option_module

  implicit none
  
  type(reactive_transport_aux_type), pointer :: ReactiveTransportAuxCreate
  
  type(reactive_transport_aux_type), pointer :: aux

  allocate(aux)  
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  ReactiveTransportAuxCreate => aux
  
end function ReactiveTransportAuxCreate

! ************************************************************************** !
!
! RTAuxVarInit: Initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  type(option_type) :: option  
  
  allocate(aux_var%total(option%ncomp))
  aux_var%total = 0.d0
  allocate(aux_var%dtotal(option%ncomp,option%ncomp))
  aux_var%dtotal = 0.d0
  
end subroutine RTAuxVarInit

! ************************************************************************** !
!
! RTAuxVarCompute: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxVarCompute(x,aux_var,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  PetscReal :: x(option%ncomp)  
  type(reactive_transport_auxvar_type) :: aux_var

  ! update totals  
  aux_var%total(1:option%ncomp) = x(1:option%ncomp)
  ! add in other later
  
end subroutine RTAuxVarCompute

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a reactive transport auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(reactive_transport_auxvar_type) :: aux_var
  
  if (associated(aux_var%total)) deallocate(aux_var%total)
  nullify(aux_var%total)
  if (associated(aux_var%dtotal))deallocate(aux_var%dtotal)
  nullify(aux_var%dtotal)

end subroutine AuxVarDestroy

! ************************************************************************** !
!
! ReactiveTransportAuxDestroy: Deallocates a reactive transport auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine ReactiveTransportAuxDestroy(aux)

  implicit none

  type(reactive_transport_aux_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    call AuxVarDestroy(aux%aux_vars(iaux))
  enddo  
  do iaux = 1, aux%num_aux_bc
    call AuxVarDestroy(aux%aux_vars_bc(iaux))
  enddo  
  
  if (associated(aux%aux_vars)) deallocate(aux%aux_vars)
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) deallocate(aux%aux_vars_bc)
  nullify(aux%aux_vars_bc)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
    
end subroutine ReactiveTransportAuxDestroy

end module Reactive_Transport_Aux_module
