module Global_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: global_auxvar_type
    PetscReal, pointer :: pres(:)
    PetscReal, pointer :: pres_store(:,:)
    PetscReal, pointer :: temp(:)
    PetscReal, pointer :: temp_store(:,:)
    PetscReal, pointer :: sat(:)
    PetscReal, pointer :: sat_store(:,:)
    PetscReal, pointer :: den(:)
    PetscReal, pointer :: den_kg(:)
    PetscReal, pointer :: den_store(:,:)
    PetscReal, pointer :: den_kg_store(:,:)
    PetscReal, pointer :: fugacoeff(:)
    PetscReal, pointer :: fugacoeff_store(:,:)
    PetscReal, pointer :: mass_balance(:)
    PetscReal, pointer :: mass_balance_delta(:)
  end type global_auxvar_type
  
  type, public :: global_type
    PetscInt :: num_aux, num_aux_bc
    type(global_auxvar_type), pointer :: aux_vars(:)
    type(global_auxvar_type), pointer :: aux_vars_bc(:)
  end type global_type

  public :: GlobalAuxCreate, GlobalAuxDestroy, &
            GlobalAuxVarInit, GlobalAuxVarCopy, &
            GlobalAuxVarDestroy

contains


! ************************************************************************** !
!
! GlobalAuxCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
function GlobalAuxCreate()

  use Option_module

  implicit none
  
  type(global_type), pointer :: GlobalAuxCreate
  
  type(global_type), pointer :: aux

  allocate(aux) 
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)

  GlobalAuxCreate => aux
  
end function GlobalAuxCreate

! ************************************************************************** !
!
! GlobalAuxVarInit: Initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine GlobalAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: aux_var
  type(option_type) :: option
  
  allocate(aux_var%pres(option%nphase))
  aux_var%pres = 0.d0
  allocate(aux_var%temp(ONE_INTEGER))
  aux_var%temp = 0.d0
  allocate(aux_var%sat(option%nphase))
  aux_var%sat = 0.d0
  allocate(aux_var%den(option%nphase))
  aux_var%den = 0.d0
  allocate(aux_var%den_kg(option%nphase))
  aux_var%den_kg = 0.d0
  allocate(aux_var%sat_store(option%nphase,TWO_INTEGER))
  aux_var%sat_store = 0.d0
  allocate(aux_var%den_kg_store(option%nphase,TWO_INTEGER))
  aux_var%den_kg_store = 0.d0

  if(option%iflowmode == IMS_MODE)then
    allocate(aux_var%pres_store(option%nphase,TWO_INTEGER))
    aux_var%pres_store = 0.d0
    allocate(aux_var%temp_store(ONE_INTEGER,TWO_INTEGER))
    aux_var%temp_store = 0.d0
    allocate(aux_var%fugacoeff(ONE_INTEGER))
    aux_var%fugacoeff = 1.d0
    allocate(aux_var%fugacoeff_store(ONE_INTEGER,TWO_INTEGER))
    aux_var%fugacoeff_store = 1.d0    
    allocate(aux_var%den_store(option%nphase,TWO_INTEGER))
    aux_var%den_store = 0.d0
  else
    nullify(aux_var%pres_store)
    nullify(aux_var%temp_store)
    nullify(aux_var%fugacoeff)
    nullify(aux_var%fugacoeff_store)
    nullify(aux_var%den_store)
  endif

  if(option%iflowmode == MPH_MODE)then
    allocate(aux_var%pres_store(option%nphase,TWO_INTEGER))
    aux_var%pres_store = 0.d0
    allocate(aux_var%temp_store(ONE_INTEGER,TWO_INTEGER))
    aux_var%temp_store = 0.d0
    allocate(aux_var%fugacoeff(ONE_INTEGER))
    aux_var%fugacoeff = 1.d0
    allocate(aux_var%fugacoeff_store(ONE_INTEGER,TWO_INTEGER))
    aux_var%fugacoeff_store = 1.d0
    allocate(aux_var%den_store(option%nphase,TWO_INTEGER))
    aux_var%den_store = 0.d0
  else
    nullify(aux_var%pres_store)
    nullify(aux_var%temp_store)
    nullify(aux_var%fugacoeff)
    nullify(aux_var%fugacoeff_store)
    nullify(aux_var%den_store)
  endif


  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(aux_var%mass_balance(option%nphase))
    aux_var%mass_balance = 0.d0
    allocate(aux_var%mass_balance_delta(option%nphase))
    aux_var%mass_balance_delta = 0.d0
  else
    nullify(aux_var%mass_balance)
    nullify(aux_var%mass_balance_delta)
  endif

end subroutine GlobalAuxVarInit

! ************************************************************************** !
!
! GlobalAuxVarCopy: Copies an auxilliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine GlobalAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%pres = aux_var%pres
  aux_var2%temp = aux_var%temp
  aux_var2%sat = aux_var%sat
  aux_var2%den = aux_var%den
  aux_var2%den_kg = aux_var%den_kg

  aux_var2%sat_store = aux_var%sat_store
  aux_var2%den_kg_store = aux_var%den_kg_store
  
  if (associated(aux_var%fugacoeff) .and. &
      associated(aux_var2%fugacoeff)) then
    aux_var2%fugacoeff = aux_var%fugacoeff  
  endif
  if (associated(aux_var%pres_store) .and. &
      associated(aux_var2%pres_store)) then
    aux_var2%pres_store = aux_var%pres_store  
  endif
  if (associated(aux_var%den_store) .and. &
      associated(aux_var2%den_store)) then
    aux_var2%pres_store = aux_var%den_store  
  endif
  if (associated(aux_var%temp_store) .and. &
      associated(aux_var2%temp_store)) then
    aux_var2%temp_store = aux_var%temp_store  
  endif
  if (associated(aux_var%fugacoeff_store) .and. &
      associated(aux_var2%fugacoeff_store)) then
    aux_var2%fugacoeff_store = aux_var%fugacoeff_store  
  endif

  
  if (associated(aux_var%mass_balance) .and. &
      associated(aux_var2%mass_balance)) then
    aux_var2%mass_balance = aux_var%mass_balance
    aux_var2%mass_balance_delta = aux_var%mass_balance_delta
  endif

end subroutine GlobalAuxVarCopy
  
! ************************************************************************** !
!
! GlobalAuxVarDestroy: Deallocates a mode auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine GlobalAuxVarDestroy(aux_var)

  implicit none

  type(global_auxvar_type) :: aux_var
  
  if (associated(aux_var%pres)) deallocate(aux_var%pres)
  nullify(aux_var%pres)
  if (associated(aux_var%temp)) deallocate(aux_var%temp)
  nullify(aux_var%temp)
  if (associated(aux_var%sat)) deallocate(aux_var%sat)
  nullify(aux_var%sat)
  if (associated(aux_var%den)) deallocate(aux_var%den)
  nullify(aux_var%den)
  if (associated(aux_var%fugacoeff)) deallocate(aux_var%fugacoeff)
  nullify(aux_var%fugacoeff)
  if (associated(aux_var%den_kg)) deallocate(aux_var%den_kg)
  nullify(aux_var%den_kg)
  
  if (associated(aux_var%pres_store)) deallocate(aux_var%pres_store)
  nullify(aux_var%pres_store)
  if (associated(aux_var%temp_store)) deallocate(aux_var%temp_store)
  nullify(aux_var%temp_store)
  if (associated(aux_var%fugacoeff_store)) deallocate(aux_var%fugacoeff_store)
  nullify(aux_var%fugacoeff_store)
  if (associated(aux_var%sat_store)) deallocate(aux_var%sat_store)
  nullify(aux_var%sat_store)
  if (associated(aux_var%den_kg_store)) deallocate(aux_var%den_kg_store)
  nullify(aux_var%den_kg_store)
  
  if (associated(aux_var%mass_balance)) deallocate(aux_var%mass_balance)
  nullify(aux_var%mass_balance)
  if (associated(aux_var%mass_balance_delta)) deallocate(aux_var%mass_balance_delta)
  nullify(aux_var%mass_balance_delta)

end subroutine GlobalAuxVarDestroy

! ************************************************************************** !
!
! GlobalAuxDestroy: Deallocates a mode auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine GlobalAuxDestroy(aux)

  implicit none

  type(global_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call GlobalAuxVarDestroy(aux%aux_vars(iaux))
    enddo  
    deallocate(aux%aux_vars)
  endif
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call GlobalAuxVarDestroy(aux%aux_vars_bc(iaux))
    enddo  
    deallocate(aux%aux_vars_bc)
  endif
  nullify(aux%aux_vars_bc)
    
end subroutine GlobalAuxDestroy

end module Global_Aux_module
