module Global_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: global_auxvar_type
    PetscInt :: istate
    PetscReal, pointer :: pres(:)
    PetscReal, pointer :: pres_store(:,:)
    PetscReal, pointer :: temp(:)
    PetscReal, pointer :: temp_store(:,:)
    PetscReal, pointer :: sat(:)
    PetscReal, pointer :: sat_store(:,:)
    PetscReal, pointer :: den(:)  ! kmol/m^3
    PetscReal, pointer :: den_kg(:) ! kg/m^3
    PetscReal, pointer :: den_store(:,:)
    PetscReal, pointer :: den_kg_store(:,:)
    PetscReal, pointer :: fugacoeff(:)
    PetscReal, pointer :: fugacoeff_store(:,:)
    PetscReal, pointer :: m_nacl(:)
    PetscReal, pointer :: xmass(:)
    PetscReal, pointer :: mass_balance(:,:) ! kg
    PetscReal, pointer :: mass_balance_delta(:,:) ! kmol
    PetscReal, pointer :: reaction_rate(:)
    PetscReal, pointer :: reaction_rate_store(:)
!   PetscReal, pointer :: reaction_rate_store(:,:)
    PetscReal, pointer :: displacement(:)
    PetscReal, pointer :: dphi(:,:)
    PetscReal :: scco2_eq_logK ! SC CO2
  end type global_auxvar_type
  
  type, public :: global_type
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(global_auxvar_type), pointer :: aux_vars(:)
    type(global_auxvar_type), pointer :: aux_vars_bc(:)
    type(global_auxvar_type), pointer :: aux_vars_ss(:)
  end type global_type
  
  interface GlobalAuxVarDestroy
    module procedure GlobalAuxVarSingleDestroy
    module procedure GlobalAuxVarArrayDestroy
  end interface GlobalAuxVarDestroy
  
  public :: GlobalAuxCreate, GlobalAuxDestroy, &
            GlobalAuxVarInit, GlobalAuxVarCopy, &
            GlobalAuxVarDestroy, GlobalAuxVarStrip

contains


! ************************************************************************** !
!
! GlobalAuxCreate: Allocate and initialize auxiliary object
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
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)

  GlobalAuxCreate => aux
  
end function GlobalAuxCreate

! ************************************************************************** !
!
! GlobalAuxVarInit: Initialize auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine GlobalAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: aux_var
  type(option_type) :: option
  
  aux_var%istate = 0

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
  allocate(aux_var%displacement(THREE_INTEGER))
  aux_var%displacement = 0.d0
  allocate(aux_var%dphi(option%nphase,THREE_INTEGER))
  aux_var%dphi = 0.d0

  aux_var%scco2_eq_logK = 0.d0

  select case(option%iflowmode)
    case(IMS_MODE, MPH_MODE, FLASH2_MODE)
      allocate(aux_var%xmass(option%nphase))
      aux_var%xmass = 1.d0
      allocate(aux_var%pres_store(option%nphase,TWO_INTEGER))
      aux_var%pres_store = option%reference_pressure
      allocate(aux_var%temp_store(ONE_INTEGER,TWO_INTEGER))
      aux_var%temp_store = 0.d0
      allocate(aux_var%fugacoeff(ONE_INTEGER))
      aux_var%fugacoeff = 1.d0
      allocate(aux_var%fugacoeff_store(ONE_INTEGER,TWO_INTEGER))
      aux_var%fugacoeff_store = 1.d0    
      allocate(aux_var%den_store(option%nphase,TWO_INTEGER))
      aux_var%den_store = 0.d0
      allocate(aux_var%m_nacl(TWO_INTEGER))
      aux_var%m_nacl = option%m_nacl
      allocate(aux_var%reaction_rate(option%nflowspec))
      aux_var%reaction_rate = 0.d0
      allocate(aux_var%reaction_rate_store(option%nflowspec))
      aux_var%reaction_rate_store = 0.d0
    ! allocate(aux_var%reaction_rate_store(option%nflowspec,TWO_INTEGER))
    ! aux_var%reaction_rate_store = 0.d0
    case(TH_MODE,THC_MODE,THMC_MODE)
    ! allocate(aux_var%xmass(option%nphase))
    ! aux_var%xmass = 1.d0
      allocate(aux_var%pres_store(option%nphase,TWO_INTEGER))
      aux_var%pres_store = 0.d0
      allocate(aux_var%temp_store(ONE_INTEGER,TWO_INTEGER))
      aux_var%temp_store = 0.d0
    ! allocate(aux_var%fugacoeff(ONE_INTEGER))
    ! aux_var%fugacoeff = 1.d0
    ! allocate(aux_var%fugacoeff_store(ONE_INTEGER,TWO_INTEGER))
    ! aux_var%fugacoeff_store = 1.d0    
      allocate(aux_var%den_store(option%nphase,TWO_INTEGER))
      aux_var%den_store = 0.d0
    ! allocate(aux_var%m_nacl(TWO_INTEGER))
    ! aux_var%m_nacl = option%m_nacl
    ! allocate(aux_var%reaction_rate(option%nflowspec))
    ! aux_var%reaction_rate = 0.d0
    ! allocate(aux_var%reaction_rate_store(option%nflowspec))
    ! aux_var%reaction_rate_store = 0.d0
    ! allocate(aux_var%reaction_rate_store(option%nflowspec,TWO_INTEGER))
    ! aux_var%reaction_rate_store = 0.d0
      nullify(aux_var%xmass)
      nullify(aux_var%fugacoeff)
      nullify(aux_var%fugacoeff_store)
      nullify(aux_var%m_nacl)
      nullify(aux_var%reaction_rate)
      nullify(aux_var%reaction_rate_store)  
    case (G_MODE)
      nullify(aux_var%xmass)
      nullify(aux_var%pres_store)
      nullify(aux_var%temp_store)
      nullify(aux_var%fugacoeff)
      nullify(aux_var%fugacoeff_store)
      nullify(aux_var%den_store)
      nullify(aux_var%m_nacl)
      nullify(aux_var%reaction_rate)
      nullify(aux_var%reaction_rate_store)  
    case default
      nullify(aux_var%xmass)
      nullify(aux_var%pres_store)
      nullify(aux_var%temp_store)
      nullify(aux_var%fugacoeff)
      nullify(aux_var%fugacoeff_store)
      nullify(aux_var%den_store)
      nullify(aux_var%m_nacl)
      nullify(aux_var%reaction_rate)
      nullify(aux_var%reaction_rate_store)
  end select
  
  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(aux_var%mass_balance(option%nflowspec,option%nphase))
    aux_var%mass_balance = 0.d0
    allocate(aux_var%mass_balance_delta(option%nflowspec,option%nphase))
    aux_var%mass_balance_delta = 0.d0
  else
    nullify(aux_var%mass_balance)
    nullify(aux_var%mass_balance_delta)
  endif
  
end subroutine GlobalAuxVarInit

! ************************************************************************** !
!
! GlobalAuxVarCopy: Copies an auxiliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine GlobalAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%istate = aux_var%istate
  aux_var2%pres = aux_var%pres
  aux_var2%temp = aux_var%temp
  aux_var2%sat = aux_var%sat
  aux_var2%den = aux_var%den
  aux_var2%den_kg = aux_var%den_kg
  aux_var2%sat_store = aux_var%sat_store
  aux_var2%den_kg_store = aux_var%den_kg_store
  aux_var2%displacement = aux_var%displacement
!  aux_var2%dphi = aux_var%dphi
  
  if (associated(aux_var%reaction_rate) .and. &
      associated(aux_var2%reaction_rate)) then
    aux_var2%reaction_rate = aux_var%reaction_rate
  endif
  
  if (associated(aux_var%m_nacl) .and. &
      associated(aux_var2%m_nacl)) then
    aux_var2%m_nacl = aux_var%m_nacl
  endif
  
  if (associated(aux_var%reaction_rate) .and. &
      associated(aux_var2%reaction_rate)) then
  endif
  
  
  if (associated(aux_var%fugacoeff) .and. &
      associated(aux_var2%fugacoeff)) then
    aux_var2%fugacoeff = aux_var%fugacoeff  
  endif
  if (associated(aux_var%xmass) .and. &
      associated(aux_var2%xmass)) then
    aux_var2%xmass = aux_var%xmass  
  endif
  if (associated(aux_var%pres_store) .and. &
      associated(aux_var2%pres_store)) then
    aux_var2%pres_store = aux_var%pres_store  
  endif
  if (associated(aux_var%den_store) .and. &
      associated(aux_var2%den_store)) then
    aux_var2%den_store = aux_var%den_store  
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
! GlobalAuxVarSingleDestroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine GlobalAuxVarSingleDestroy(aux_var)

  implicit none

  type(global_auxvar_type), pointer :: aux_var
  
  if (associated(aux_var)) then
    call GlobalAuxVarStrip(aux_var)
    deallocate(aux_var)
  endif
  nullify(aux_var)

end subroutine GlobalAuxVarSingleDestroy
  
! ************************************************************************** !
!
! GlobalAuxVarArrayDestroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine GlobalAuxVarArrayDestroy(aux_vars)

  implicit none

  type(global_auxvar_type), pointer :: aux_vars(:)
  
  PetscInt :: iaux
  
  if (associated(aux_vars)) then
    do iaux = 1, size(aux_vars)
      call GlobalAuxVarStrip(aux_vars(iaux))
    enddo  
    deallocate(aux_vars)
  endif
  nullify(aux_vars)

end subroutine GlobalAuxVarArrayDestroy
  
! ************************************************************************** !
!
! GlobalAuxVarStrip: Deallocates all members of single auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine GlobalAuxVarStrip(aux_var)

  use Utility_module, only: DeallocateArray

  implicit none

  type(global_auxvar_type) :: aux_var
  
  call DeallocateArray(aux_var%pres)
  call DeallocateArray(aux_var%temp)
  call DeallocateArray(aux_var%sat)
  call DeallocateArray(aux_var%den)
  call DeallocateArray(aux_var%fugacoeff)
  call DeallocateArray(aux_var%den_kg)
  call DeallocateArray(aux_var%m_nacl)
  call DeallocateArray(aux_var%xmass)
  call DeallocateArray(aux_var%reaction_rate)
  call DeallocateArray(aux_var%displacement)
  call DeallocateArray(aux_var%dphi)

  call DeallocateArray(aux_var%pres_store)
  call DeallocateArray(aux_var%temp_store)
  call DeallocateArray(aux_var%fugacoeff_store)
  call DeallocateArray(aux_var%sat_store)
  call DeallocateArray(aux_var%den_kg_store)
  
  call DeallocateArray(aux_var%mass_balance)
  call DeallocateArray(aux_var%mass_balance_delta)

end subroutine GlobalAuxVarStrip

! ************************************************************************** !
!
! GlobalAuxDestroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine GlobalAuxDestroy(aux)

  implicit none

  type(global_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  call GlobalAuxVarDestroy(aux%aux_vars)
  call GlobalAuxVarDestroy(aux%aux_vars_bc)
  call GlobalAuxVarDestroy(aux%aux_vars_ss)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GlobalAuxDestroy

end module Global_Aux_module
