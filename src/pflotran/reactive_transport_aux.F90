module Reactive_Transport_Aux_module

  ! this module cannot depend on any other modules beside Option_module

  implicit none
  
  private 

#include "definitions.h"
 
  type, public :: reactive_transport_auxvar_type
    ! molality
    PetscReal, pointer :: pri_molal(:)     ! mol/kg water
    ! phase dependent totals
    PetscReal, pointer :: total(:,:)       ! mol solute/L water
    PetscReal, pointer :: dtotal(:,:,:)    ! kg water/m^3 water
    ! sorbed totals
    PetscReal, pointer :: total_sorb(:)    ! mol/m^3 bulk
    PetscReal, pointer :: dtotal_sorb(:,:) ! kg water/m^3 bulk
    ! aqueous species
    ! aqueous complexes
    PetscReal, pointer :: sec_molal(:)
    PetscReal, pointer :: gas_molal(:)
    ! sorption reactions
    ! PetscReal, pointer :: kinsurfcmplx_spec(:)
    ! PetscReal, pointer :: kinionx_molfrac(:)
    PetscReal, pointer :: eqsurfcmplx_conc(:)
    PetscReal, pointer :: eqsurfcmplx_freesite_conc(:)
!   PetscReal, pointer :: eqsurf_site_density(:)
    PetscReal, pointer :: eqionx_ref_cation_sorbed_conc(:)
    PetscReal, pointer :: eqionx_conc(:,:)
!   PetscReal, pointer :: eqionx_cec(:)
    ! PetscReal, pointer :: eqionx_molfrac(:)
    ! mineral reactions
    PetscReal, pointer :: mnrl_volfrac0(:)
    PetscReal, pointer :: mnrl_volfrac(:)
    PetscReal, pointer :: mnrl_area0(:)
    PetscReal, pointer :: mnrl_area(:)
    PetscReal, pointer :: mnrl_rate(:)
    ! activity coefficients
    PetscReal :: act_h2o
    PetscReal, pointer :: pri_act_coef(:)
    PetscReal, pointer :: sec_act_coef(:)
    
    PetscReal, pointer :: mass_balance(:,:)
    PetscReal, pointer :: mass_balance_delta(:,:)
    
    PetscReal, pointer :: kinmr_total_sorb(:,:)
    PetscReal, pointer :: kinmr_total_sorb_prev(:,:)
    
  end type reactive_transport_auxvar_type
  
  type, public :: reactive_transport_param_type
    PetscReal :: dispersivity
    PetscReal, pointer :: diffusion_coefficient(:)
  end type reactive_transport_param_type

  type, public :: reactive_transport_type
    PetscInt :: num_aux, num_aux_bc
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscInt :: n_zero_rows
    PetscTruth :: aux_vars_up_to_date
    PetscTruth :: inactive_cells_exist
    type(reactive_transport_param_type), pointer :: rt_parameter
    type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
    type(reactive_transport_auxvar_type), pointer :: aux_vars_bc(:)
  end type reactive_transport_type
  
  public :: RTAuxCreate, RTAuxDestroy, &
            RTAuxVarInit, RTAuxVarCopy, RTAuxVarDestroy
            
contains


! ************************************************************************** !
!
! RTAuxCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
function RTAuxCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
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
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE

  allocate(aux%rt_parameter)
  allocate(aux%rt_parameter%diffusion_coefficient(option%nphase))
  aux%rt_parameter%diffusion_coefficient = 0.d0
  aux%rt_parameter%dispersivity = 0.d0

  RTAuxCreate => aux
  
end function RTAuxCreate

! ************************************************************************** !
!
! RTAuxVarInit: Initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxVarInit(aux_var,reaction,option)

  use Option_module
  use Reaction_Aux_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  type(reaction_type) :: reaction
  type(option_type) :: option  
  
  allocate(aux_var%pri_molal(reaction%ncomp))
  aux_var%pri_molal = 0.d0
  allocate(aux_var%total(reaction%ncomp,option%nphase))
  aux_var%total = 0.d0
  allocate(aux_var%dtotal(reaction%ncomp,reaction%ncomp,option%nphase))
  aux_var%dtotal = 0.d0
  
  if (reaction%neqcmplx > 0) then
    allocate(aux_var%sec_molal(reaction%neqcmplx))
    aux_var%sec_molal = 0.d0
  else
    nullify(aux_var%sec_molal)
  endif
  
  if (reaction%ngas > 0) then
    allocate(aux_var%gas_molal(reaction%ngas))
    aux_var%gas_molal = 0.d0
  else
    nullify(aux_var%gas_molal)
  endif

  
  if (reaction%nsorb > 0) then  
    allocate(aux_var%total_sorb(reaction%ncomp))
    aux_var%total_sorb = 0.d0
    allocate(aux_var%dtotal_sorb(reaction%ncomp,reaction%ncomp))
    aux_var%dtotal_sorb = 0.d0
  else
    nullify(aux_var%total_sorb)
    nullify(aux_var%dtotal_sorb)
  endif    
  
  if (reaction%neqsurfcmplxrxn > 0) then
    allocate(aux_var%eqsurfcmplx_conc(reaction%neqsurfcmplx))
    aux_var%eqsurfcmplx_conc = 0.d0
    allocate(aux_var%eqsurfcmplx_freesite_conc(reaction%neqsurfcmplxrxn))
    aux_var%eqsurfcmplx_freesite_conc = 1.d-9 ! initialize to guess
!   allocate(aux_var%eqsurf_site_density(reaction%neqsurfcmplxrxn))
!   aux_var%eqsurf_site_density = 0.d0
  else
    nullify(aux_var%eqsurfcmplx_conc)
    nullify(aux_var%eqsurfcmplx_freesite_conc)
!   nullify(aux_var%eqsurf_site_density)
  endif
  
  if (reaction%neqionxrxn > 0) then
    allocate(aux_var%eqionx_ref_cation_sorbed_conc(reaction%neqionxrxn))
    aux_var%eqionx_ref_cation_sorbed_conc = 1.d-9 ! initialize to guess
    allocate(aux_var%eqionx_conc(reaction%neqionxcation,reaction%neqionxrxn))
    aux_var%eqionx_conc = 1.d-9
!   allocate(aux_var%eqionx_cec(reaction%neqionxcation))
!   aux_var%eqionx_cec = 0.d0
  else
    nullify(aux_var%eqionx_ref_cation_sorbed_conc)
    nullify(aux_var%eqionx_conc)
!   nullify(aux_var%eqionx_cec)
  endif
  
  if (reaction%nkinmnrl > 0) then
    allocate(aux_var%mnrl_volfrac0(reaction%nkinmnrl))
    aux_var%mnrl_volfrac0 = 0.d0
    allocate(aux_var%mnrl_volfrac(reaction%nkinmnrl))
    aux_var%mnrl_volfrac = 0.d0
    allocate(aux_var%mnrl_area0(reaction%nkinmnrl))
    aux_var%mnrl_area0 = 0.d0
    allocate(aux_var%mnrl_area(reaction%nkinmnrl))
    aux_var%mnrl_area = 0.d0
    allocate(aux_var%mnrl_rate(reaction%nkinmnrl))
    aux_var%mnrl_rate = 0.d0
  else
    nullify(aux_var%mnrl_volfrac0)
    nullify(aux_var%mnrl_volfrac)
    nullify(aux_var%mnrl_area0)
    nullify(aux_var%mnrl_area)
    nullify(aux_var%mnrl_rate)
  endif
  
  aux_var%act_h2o = 1.d0
  allocate(aux_var%pri_act_coef(reaction%ncomp))
  aux_var%pri_act_coef = 1.d0
  if (reaction%neqcmplx > 0) then
    allocate(aux_var%sec_act_coef(reaction%neqcmplx))
    aux_var%sec_act_coef = 1.d0
  else
    nullify(aux_var%sec_act_coef)
  endif
  
  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(aux_var%mass_balance(reaction%ncomp,option%nphase))
    aux_var%mass_balance = 0.d0
    allocate(aux_var%mass_balance_delta(reaction%ncomp,option%nphase))
    aux_var%mass_balance_delta = 0.d0
  else
    nullify(aux_var%mass_balance)
    nullify(aux_var%mass_balance_delta)
  endif
  
  if (reaction%nkinmr_rates > 0) then
    allocate(aux_var%kinmr_total_sorb(reaction%ncomp,reaction%nkinmr_rates))
    aux_var%kinmr_total_sorb = 0.d0
    allocate(aux_var%kinmr_total_sorb_prev(reaction%ncomp,reaction%nkinmr_rates))
    aux_var%kinmr_total_sorb_prev = 0.d0
  else
    nullify(aux_var%kinmr_total_sorb)
    nullify(aux_var%kinmr_total_sorb_prev)
  endif
  
end subroutine RTAuxVarInit

! ************************************************************************** !
!
! RTAuxVarCopy: Copys an auxilliary object
! author: Glenn Hammond
! date: 09/05/08
!
! ************************************************************************** !
subroutine RTAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option  
  
  aux_var%pri_molal = aux_var2%pri_molal
  aux_var%total = aux_var2%total
  aux_var%dtotal = aux_var2%dtotal
  
  if (associated(aux_var%sec_molal)) &
    aux_var%sec_molal = aux_var2%sec_molal
  if (associated(aux_var%total_sorb)) then  
    aux_var%total_sorb = aux_var2%total_sorb
    aux_var%dtotal_sorb = aux_var2%dtotal_sorb
  endif
  
  if (associated(aux_var%gas_molal)) &
    aux_var%gas_molal = aux_var2%gas_molal
  
  if (associated(aux_var%eqsurfcmplx_conc)) then
    aux_var%eqsurfcmplx_conc = aux_var2%eqsurfcmplx_conc
    aux_var%eqsurfcmplx_freesite_conc = aux_var2%eqsurfcmplx_freesite_conc
  endif
  
  if (associated(aux_var%eqionx_ref_cation_sorbed_conc)) then
    aux_var%eqionx_ref_cation_sorbed_conc = aux_var2%eqionx_ref_cation_sorbed_conc
    aux_var%eqionx_conc = aux_var2%eqionx_conc
  endif  
  
  if (associated(aux_var%mnrl_volfrac)) then
    aux_var%mnrl_volfrac0 = aux_var2%mnrl_volfrac0
    aux_var%mnrl_volfrac = aux_var2%mnrl_volfrac
    aux_var%mnrl_area0 = aux_var2%mnrl_area0
    aux_var%mnrl_area = aux_var2%mnrl_area
    aux_var%mnrl_rate = aux_var2%mnrl_rate
  endif
  
  aux_var%act_h2o = aux_var2%act_h2o
  aux_var%pri_act_coef = aux_var2%pri_act_coef
  if (associated(aux_var%sec_act_coef)) &
    aux_var%sec_act_coef = aux_var2%sec_act_coef

  if (associated(aux_var%mass_balance) .and. &
      associated(aux_var2%mass_balance)) then
    aux_var%mass_balance = aux_var2%mass_balance
    aux_var%mass_balance_delta = aux_var2%mass_balance_delta
  endif

  if (associated(aux_var%kinmr_total_sorb) .and. &
      associated(aux_var2%kinmr_total_sorb)) then
    aux_var%kinmr_total_sorb = aux_var2%kinmr_total_sorb
    aux_var%kinmr_total_sorb_prev = aux_var2%kinmr_total_sorb_prev
  endif

end subroutine RTAuxVarCopy


! ************************************************************************** !
!
! RTAuxVarDestroy: Deallocates a reactive transport auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxVarDestroy(aux_var)

  implicit none

  type(reactive_transport_auxvar_type) :: aux_var
  
  if (associated(aux_var%pri_molal)) deallocate(aux_var%pri_molal)
  nullify(aux_var%pri_molal)
  
  if (associated(aux_var%total)) deallocate(aux_var%total)
  nullify(aux_var%total)
  if (associated(aux_var%dtotal))deallocate(aux_var%dtotal)
  nullify(aux_var%dtotal)
  
  if (associated(aux_var%sec_molal))deallocate(aux_var%sec_molal)
  nullify(aux_var%sec_molal)
  
  if (associated(aux_var%gas_molal))deallocate(aux_var%gas_molal)
  nullify(aux_var%gas_molal)
  
  if (associated(aux_var%total_sorb)) deallocate(aux_var%total_sorb)
  nullify(aux_var%total_sorb)
  if (associated(aux_var%dtotal_sorb))deallocate(aux_var%dtotal_sorb)
  nullify(aux_var%dtotal_sorb)
  
  if (associated(aux_var%eqionx_ref_cation_sorbed_conc)) deallocate(aux_var%eqionx_ref_cation_sorbed_conc)
  nullify(aux_var%eqionx_ref_cation_sorbed_conc)
  if (associated(aux_var%eqionx_conc)) deallocate(aux_var%eqionx_conc)
  nullify(aux_var%eqionx_conc)
  
  if (associated(aux_var%mnrl_volfrac0))deallocate(aux_var%mnrl_volfrac0)
  nullify(aux_var%mnrl_volfrac0)
  if (associated(aux_var%mnrl_volfrac))deallocate(aux_var%mnrl_volfrac)
  nullify(aux_var%mnrl_volfrac)
  if (associated(aux_var%mnrl_area0))deallocate(aux_var%mnrl_area0)
  nullify(aux_var%mnrl_area0)
  if (associated(aux_var%mnrl_area))deallocate(aux_var%mnrl_area)
  nullify(aux_var%mnrl_area)
  if (associated(aux_var%mnrl_rate))deallocate(aux_var%mnrl_rate)
  nullify(aux_var%mnrl_rate)
  
  if (associated(aux_var%pri_act_coef))deallocate(aux_var%pri_act_coef)
  nullify(aux_var%pri_act_coef)
  if (associated(aux_var%sec_act_coef))deallocate(aux_var%sec_act_coef)
  nullify(aux_var%sec_act_coef)

  if (associated(aux_var%mass_balance)) deallocate(aux_var%mass_balance)
  nullify(aux_var%mass_balance)
  if (associated(aux_var%mass_balance_delta)) deallocate(aux_var%mass_balance_delta)
  nullify(aux_var%mass_balance_delta)

  if (associated(aux_var%kinmr_total_sorb)) deallocate(aux_var%kinmr_total_sorb)
  nullify(aux_var%kinmr_total_sorb)
  if (associated(aux_var%kinmr_total_sorb_prev)) deallocate(aux_var%kinmr_total_sorb_prev)
  nullify(aux_var%kinmr_total_sorb_prev)

end subroutine RTAuxVarDestroy

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
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call RTAuxVarDestroy(aux%aux_vars(iaux))
    enddo  
    deallocate(aux%aux_vars)
  endif
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call RTAuxVarDestroy(aux%aux_vars_bc(iaux))
    enddo  
    deallocate(aux%aux_vars_bc)
  endif
  nullify(aux%aux_vars_bc)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%rt_parameter)) then
    deallocate(aux%rt_parameter%diffusion_coefficient)
    nullify(aux%rt_parameter%diffusion_coefficient)
    deallocate(aux%rt_parameter)
  endif
  nullify(aux%rt_parameter)
    
end subroutine RTAuxDestroy

end module Reactive_Transport_Aux_module
