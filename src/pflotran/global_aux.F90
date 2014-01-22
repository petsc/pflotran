module Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

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
    PetscReal, pointer :: dphi(:,:)
    PetscReal :: scco2_eq_logK ! SC CO2
  end type global_auxvar_type
  
  type, public :: global_type
    PetscReal :: time_t, time_tpdt
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(global_auxvar_type), pointer :: auxvars(:)
    type(global_auxvar_type), pointer :: auxvars_bc(:)
    type(global_auxvar_type), pointer :: auxvars_ss(:)
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

function GlobalAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(global_type), pointer :: GlobalAuxCreate
  
  type(global_type), pointer :: aux

  allocate(aux) 
  aux%time_t = 0.d0
  aux%time_tpdt = 0.d0
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  GlobalAuxCreate => aux
  
end function GlobalAuxCreate

! ************************************************************************** !

subroutine GlobalAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%istate = 0

  allocate(auxvar%pres(option%nphase))
  auxvar%pres = 0.d0
  allocate(auxvar%temp(ONE_INTEGER))
  auxvar%temp = 0.d0
  allocate(auxvar%sat(option%nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den(option%nphase))
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(option%nphase))
  auxvar%den_kg = 0.d0
  allocate(auxvar%sat_store(option%nphase,TWO_INTEGER))
  auxvar%sat_store = 0.d0
  allocate(auxvar%den_kg_store(option%nphase,TWO_INTEGER))
  auxvar%den_kg_store = 0.d0
  allocate(auxvar%dphi(option%nphase,THREE_INTEGER))
  auxvar%dphi = 0.d0

  auxvar%scco2_eq_logK = 0.d0

  select case(option%iflowmode)
    case(IMS_MODE, MPH_MODE, FLASH2_MODE)
      allocate(auxvar%xmass(option%nphase))
      auxvar%xmass = 1.d0
      allocate(auxvar%pres_store(option%nphase,TWO_INTEGER))
      auxvar%pres_store = option%reference_pressure
      allocate(auxvar%temp_store(ONE_INTEGER,TWO_INTEGER))
      auxvar%temp_store = 0.d0
      allocate(auxvar%fugacoeff(ONE_INTEGER))
      auxvar%fugacoeff = 1.d0
      allocate(auxvar%fugacoeff_store(ONE_INTEGER,TWO_INTEGER))
      auxvar%fugacoeff_store = 1.d0    
      allocate(auxvar%den_store(option%nphase,TWO_INTEGER))
      auxvar%den_store = 0.d0
      allocate(auxvar%m_nacl(TWO_INTEGER))
      auxvar%m_nacl = option%m_nacl
      allocate(auxvar%reaction_rate(option%nflowspec))
      auxvar%reaction_rate = 0.d0
      allocate(auxvar%reaction_rate_store(option%nflowspec))
      auxvar%reaction_rate_store = 0.d0
    ! allocate(auxvar%reaction_rate_store(option%nflowspec,TWO_INTEGER))
    ! auxvar%reaction_rate_store = 0.d0
    case(TH_MODE,THC_MODE)
    ! allocate(auxvar%xmass(option%nphase))
    ! auxvar%xmass = 1.d0
      allocate(auxvar%pres_store(option%nphase,TWO_INTEGER))
      auxvar%pres_store = 0.d0
      allocate(auxvar%temp_store(ONE_INTEGER,TWO_INTEGER))
      auxvar%temp_store = 0.d0
    ! allocate(auxvar%fugacoeff(ONE_INTEGER))
    ! auxvar%fugacoeff = 1.d0
    ! allocate(auxvar%fugacoeff_store(ONE_INTEGER,TWO_INTEGER))
    ! auxvar%fugacoeff_store = 1.d0    
      allocate(auxvar%den_store(option%nphase,TWO_INTEGER))
      auxvar%den_store = 0.d0
    ! allocate(auxvar%m_nacl(TWO_INTEGER))
    ! auxvar%m_nacl = option%m_nacl
    ! allocate(auxvar%reaction_rate(option%nflowspec))
    ! auxvar%reaction_rate = 0.d0
    ! allocate(auxvar%reaction_rate_store(option%nflowspec))
    ! auxvar%reaction_rate_store = 0.d0
    ! allocate(auxvar%reaction_rate_store(option%nflowspec,TWO_INTEGER))
    ! auxvar%reaction_rate_store = 0.d0
      nullify(auxvar%xmass)
      nullify(auxvar%fugacoeff)
      nullify(auxvar%fugacoeff_store)
      nullify(auxvar%m_nacl)
      nullify(auxvar%reaction_rate)
      nullify(auxvar%reaction_rate_store)  
    case (G_MODE)
      nullify(auxvar%xmass)
      nullify(auxvar%pres_store)
      nullify(auxvar%temp_store)
      nullify(auxvar%fugacoeff)
      nullify(auxvar%fugacoeff_store)
      nullify(auxvar%den_store)
      nullify(auxvar%m_nacl)
      nullify(auxvar%reaction_rate)
      nullify(auxvar%reaction_rate_store)  
    case default
      nullify(auxvar%xmass)
      nullify(auxvar%pres_store)
      nullify(auxvar%temp_store)
      nullify(auxvar%fugacoeff)
      nullify(auxvar%fugacoeff_store)
      nullify(auxvar%den_store)
      nullify(auxvar%m_nacl)
      nullify(auxvar%reaction_rate)
      nullify(auxvar%reaction_rate_store)
  end select
  
  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(auxvar%mass_balance(option%nflowspec,option%nphase))
    auxvar%mass_balance = 0.d0
    allocate(auxvar%mass_balance_delta(option%nflowspec,option%nphase))
    auxvar%mass_balance_delta = 0.d0
  else
    nullify(auxvar%mass_balance)
    nullify(auxvar%mass_balance_delta)
  endif
  
end subroutine GlobalAuxVarInit

! ************************************************************************** !

subroutine GlobalAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate = auxvar%istate
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%sat_store = auxvar%sat_store
  auxvar2%den_kg_store = auxvar%den_kg_store
!  auxvar2%dphi = auxvar%dphi
  
  if (associated(auxvar%reaction_rate) .and. &
      associated(auxvar2%reaction_rate)) then
    auxvar2%reaction_rate = auxvar%reaction_rate
  endif
  
  if (associated(auxvar%m_nacl) .and. &
      associated(auxvar2%m_nacl)) then
    auxvar2%m_nacl = auxvar%m_nacl
  endif
  
  if (associated(auxvar%reaction_rate) .and. &
      associated(auxvar2%reaction_rate)) then
  endif
  
  
  if (associated(auxvar%fugacoeff) .and. &
      associated(auxvar2%fugacoeff)) then
    auxvar2%fugacoeff = auxvar%fugacoeff  
  endif
  if (associated(auxvar%xmass) .and. &
      associated(auxvar2%xmass)) then
    auxvar2%xmass = auxvar%xmass  
  endif
  if (associated(auxvar%pres_store) .and. &
      associated(auxvar2%pres_store)) then
    auxvar2%pres_store = auxvar%pres_store  
  endif
  if (associated(auxvar%den_store) .and. &
      associated(auxvar2%den_store)) then
    auxvar2%den_store = auxvar%den_store  
  endif
  if (associated(auxvar%temp_store) .and. &
      associated(auxvar2%temp_store)) then
    auxvar2%temp_store = auxvar%temp_store  
  endif
  if (associated(auxvar%fugacoeff_store) .and. &
      associated(auxvar2%fugacoeff_store)) then
    auxvar2%fugacoeff_store = auxvar%fugacoeff_store  
  endif

  if (associated(auxvar%mass_balance) .and. &
      associated(auxvar2%mass_balance)) then
    auxvar2%mass_balance = auxvar%mass_balance
    auxvar2%mass_balance_delta = auxvar%mass_balance_delta
  endif

end subroutine GlobalAuxVarCopy

! ************************************************************************** !

subroutine GlobalAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(global_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call GlobalAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine GlobalAuxVarSingleDestroy

! ************************************************************************** !

subroutine GlobalAuxVarArrayDestroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(global_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call GlobalAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine GlobalAuxVarArrayDestroy

! ************************************************************************** !

subroutine GlobalAuxVarStrip(auxvar)
  ! 
  ! Deallocates all members of single auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(global_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)
  call DeallocateArray(auxvar%temp)
  call DeallocateArray(auxvar%sat)
  call DeallocateArray(auxvar%den)
  call DeallocateArray(auxvar%fugacoeff)
  call DeallocateArray(auxvar%den_kg)
  call DeallocateArray(auxvar%m_nacl)
  call DeallocateArray(auxvar%xmass)
  call DeallocateArray(auxvar%reaction_rate)
  call DeallocateArray(auxvar%dphi)

  call DeallocateArray(auxvar%pres_store)
  call DeallocateArray(auxvar%temp_store)
  call DeallocateArray(auxvar%fugacoeff_store)
  call DeallocateArray(auxvar%sat_store)
  call DeallocateArray(auxvar%den_kg_store)
  
  call DeallocateArray(auxvar%mass_balance)
  call DeallocateArray(auxvar%mass_balance_delta)

end subroutine GlobalAuxVarStrip

! ************************************************************************** !

subroutine GlobalAuxDestroy(aux)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(global_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  call GlobalAuxVarDestroy(aux%auxvars)
  call GlobalAuxVarDestroy(aux%auxvars_bc)
  call GlobalAuxVarDestroy(aux%auxvars_ss)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GlobalAuxDestroy

end module Global_Aux_module
