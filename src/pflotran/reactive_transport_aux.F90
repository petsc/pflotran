module Reactive_Transport_Aux_module

  ! this module cannot depend on any other modules beside Option_module

  implicit none
  
  private 

#include "definitions.h"
 
  type, public :: reactive_transport_auxvar_type
    ! phase dependent totals
    PetscReal, pointer :: total(:,:)
    PetscReal, pointer :: dtotal(:,:,:)
    ! aqueous species
    PetscReal, pointer :: primary_spec(:)
    ! aqueous complexes
    PetscReal, pointer :: secondary_spec(:)
    ! sorption reactions
    PetscReal, pointer :: kinsurfcmplx_spec(:)
    PetscReal, pointer :: kinionx_molfrac(:)
    PetscReal, pointer :: eqsurfcmplx_spec(:)
    PetscReal, pointer :: eqionx_molfrac(:)
    ! mineral reactions
    PetscReal, pointer :: mnrl_volfrac(:)
    PetscReal, pointer :: mnrl_area0(:)
  end type reactive_transport_auxvar_type
  
  type, public :: reactive_transport_type
    PetscInt :: num_aux, num_aux_bc
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscInt :: n_zero_rows
    type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
    type(reactive_transport_auxvar_type), pointer :: aux_vars_bc(:)
  end type reactive_transport_type

  public :: RTAuxCreate, RTAuxDestroy, &
            RTAuxVarInit

contains


! ************************************************************************** !
!
! RTAuxCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
function RTAuxCreate()

  use Option_module

  implicit none
  
  type(reactive_transport_type), pointer :: RTAuxCreate
  
  type(reactive_transport_type), pointer :: aux

  allocate(aux)  
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  RTAuxCreate => aux
  
end function RTAuxCreate

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
  
  allocate(aux_var%total(option%ncomp,option%nphase))
  aux_var%total = 0.d0
  allocate(aux_var%dtotal(option%ncomp,option%ncomp,option%nphase))
  aux_var%dtotal = 0.d0
  allocate(aux_var%primary_spec(option%ncomp))
  aux_var%primary_spec = 0.d0
  allocate(aux_var%secondary_spec(option%ncmplx))
  aux_var%secondary_spec = 0.d0
  
end subroutine RTAuxVarInit

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
! RTAuxDestroy: Deallocates a reactive transport auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxDestroy(aux)

  implicit none

  type(reactive_transport_type), pointer :: aux
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
    
end subroutine RTAuxDestroy

end module Reactive_Transport_Aux_module
