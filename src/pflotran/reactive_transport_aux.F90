module Reactive_Transport_Aux_module

  ! this module cannot depend on any other modules beside Option_module

  implicit none
  
  private 

#include "definitions.h"
 
  type, public :: reactive_transport_auxvar_type
    PetscReal, pointer :: den(:) ! kg water / m^3 water
    PetscReal, pointer :: sat(:)
    PetscReal :: temp
    ! molality
    PetscReal, pointer :: primary_molal(:) ! kg solute / L water
    ! phase dependent totals
    PetscReal, pointer :: total(:,:) ! mol solute / L water
    PetscReal, pointer :: dtotal(:,:,:) ! kg water / m^3 water
    ! sorbed totals
    PetscReal, pointer :: total_sorb(:) ! mol solute / L water
    PetscReal, pointer :: dtotal_sorb(:,:) ! kg water / m^3 water
    ! aqueous species
    PetscReal, pointer :: primary_spec(:) ! mol solute / L water
    ! aqueous complexes
    PetscReal, pointer :: secondary_spec(:)
    ! sorption reactions
    ! PetscReal, pointer :: kinsurfcmplx_spec(:)
    ! PetscReal, pointer :: kinionx_molfrac(:)
    PetscReal, pointer :: eqsurfcmplx_spec(:)
    PetscReal, pointer :: eqsurfcmplx_freesite_conc(:)
    PetscReal, pointer :: eqionx_ref_cation_sorbed_conc(:)
    ! PetscReal, pointer :: eqionx_molfrac(:)
    ! mineral reactions
    PetscReal, pointer :: mnrl_volfrac(:)
    PetscReal, pointer :: mnrl_area0(:)
    PetscReal, pointer :: mnrl_rate(:)
    ! activity coefficients
    PetscReal :: act_h2o
    PetscReal, pointer :: pri_act_coef(:)
    PetscReal, pointer :: sec_act_coef(:)
  end type reactive_transport_auxvar_type
  
  type, public :: reactive_transport_type
    PetscInt :: num_aux, num_aux_bc
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscInt :: n_zero_rows
    PetscTruth :: aux_vars_up_to_date
    PetscTruth :: inactive_cells_exist
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
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE

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
  
  aux_var%temp = option%tref
  allocate(aux_var%den(option%nphase))
  aux_var%den = 0.d0
  allocate(aux_var%sat(option%nphase))
  aux_var%sat = 0.d0
  allocate(aux_var%primary_molal(reaction%ncomp))
  aux_var%primary_molal = 0.d0
  allocate(aux_var%total(reaction%ncomp,option%nphase))
  aux_var%total = 0.d0
  allocate(aux_var%dtotal(reaction%ncomp,reaction%ncomp,option%nphase))
  aux_var%dtotal = 0.d0
  
  allocate(aux_var%primary_spec(reaction%ncomp))
  aux_var%primary_spec = 0.d0
  
  if (reaction%neqcmplx > 0) then
    allocate(aux_var%secondary_spec(reaction%neqcmplx))
    aux_var%secondary_spec = 0.d0
  else
    nullify(aux_var%secondary_spec)
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
    allocate(aux_var%eqsurfcmplx_spec(reaction%neqsurfcmplx))
    aux_var%eqsurfcmplx_spec = 0.d0
    allocate(aux_var%eqsurfcmplx_freesite_conc(reaction%neqsurfcmplxrxn))
    aux_var%eqsurfcmplx_freesite_conc = 1.d-9 ! initialize to guess
  else
    nullify(aux_var%eqsurfcmplx_spec)
    nullify(aux_var%eqsurfcmplx_freesite_conc)
  endif
  
  if (reaction%neqionxrxn > 0) then
    allocate(aux_var%eqionx_ref_cation_sorbed_conc(reaction%neqionxrxn))
    aux_var%eqionx_ref_cation_sorbed_conc = 1.d-9 ! initialize to guess
  else
    nullify(aux_var%eqionx_ref_cation_sorbed_conc)
  endif
  
  if (reaction%nmnrl > 0) then
    allocate(aux_var%mnrl_volfrac(reaction%nmnrl))
    aux_var%mnrl_volfrac = 0.d0
    allocate(aux_var%mnrl_area0(reaction%nmnrl))
    aux_var%mnrl_area0 = 1.d0 ! Hardwired for now - geh
    allocate(aux_var%mnrl_rate(reaction%nmnrl))
    aux_var%mnrl_rate = 0.d0
  else
    nullify(aux_var%mnrl_volfrac)
    nullify(aux_var%mnrl_area0)
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
  
  aux_var%temp = aux_var2%temp
  aux_var%den = aux_var2%den
  aux_var%sat = aux_var2%sat
  aux_var%primary_molal = aux_var2%primary_molal
  aux_var%total = aux_var2%total
  aux_var%dtotal = aux_var2%dtotal
  aux_var%primary_spec = aux_var2%primary_spec
  if (associated(aux_var%secondary_spec)) &
    aux_var%secondary_spec = aux_var2%secondary_spec
  if (associated(aux_var%total_sorb)) then  
    aux_var%total_sorb = aux_var2%total_sorb
    aux_var%dtotal_sorb = aux_var2%dtotal_sorb
  endif
  if (associated(aux_var%eqsurfcmplx_spec)) then
    aux_var%eqsurfcmplx_spec = aux_var2%eqsurfcmplx_spec
    aux_var%eqsurfcmplx_freesite_conc = aux_var2%eqsurfcmplx_freesite_conc
  endif
  if (associated(aux_var%eqionx_ref_cation_sorbed_conc)) then
    aux_var%eqionx_ref_cation_sorbed_conc = aux_var2%eqionx_ref_cation_sorbed_conc
  endif  
  if (associated(aux_var%mnrl_volfrac)) then
    aux_var%mnrl_volfrac = aux_var2%mnrl_volfrac
    aux_var%mnrl_area0 = aux_var2%mnrl_area0
    aux_var%mnrl_rate = aux_var2%mnrl_rate
  endif
  aux_var%act_h2o = aux_var2%act_h2o
  aux_var%pri_act_coef = aux_var2%pri_act_coef
  if (associated(aux_var%sec_act_coef)) &
    aux_var%sec_act_coef = aux_var2%sec_act_coef

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
  
  if (associated(aux_var%den)) deallocate(aux_var%den)
  nullify(aux_var%den)
  if (associated(aux_var%sat)) deallocate(aux_var%sat)
  nullify(aux_var%sat)
  if (associated(aux_var%primary_molal)) deallocate(aux_var%primary_molal)
  nullify(aux_var%primary_molal)
  if (associated(aux_var%total)) deallocate(aux_var%total)
  nullify(aux_var%total)
  if (associated(aux_var%dtotal))deallocate(aux_var%dtotal)
  nullify(aux_var%dtotal)
  if (associated(aux_var%primary_spec))deallocate(aux_var%primary_spec)
  nullify(aux_var%primary_spec)
  if (associated(aux_var%secondary_spec))deallocate(aux_var%secondary_spec)
  nullify(aux_var%secondary_spec)
  if (associated(aux_var%total_sorb)) deallocate(aux_var%total_sorb)
  nullify(aux_var%total_sorb)
  if (associated(aux_var%dtotal_sorb))deallocate(aux_var%dtotal_sorb)
  nullify(aux_var%dtotal_sorb)
  if (associated(aux_var%eqionx_ref_cation_sorbed_conc)) deallocate(aux_var%eqionx_ref_cation_sorbed_conc)
  nullify(aux_var%eqionx_ref_cation_sorbed_conc)  
  if (associated(aux_var%mnrl_volfrac))deallocate(aux_var%mnrl_volfrac)
  nullify(aux_var%mnrl_volfrac)
  if (associated(aux_var%mnrl_area0))deallocate(aux_var%mnrl_area0)
  nullify(aux_var%mnrl_area0)
  if (associated(aux_var%mnrl_rate))deallocate(aux_var%mnrl_rate)
  nullify(aux_var%mnrl_rate)
  if (associated(aux_var%pri_act_coef))deallocate(aux_var%pri_act_coef)
  nullify(aux_var%pri_act_coef)
  if (associated(aux_var%sec_act_coef))deallocate(aux_var%sec_act_coef)
  nullify(aux_var%sec_act_coef)

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
  
  do iaux = 1, aux%num_aux
    call RTAuxVarDestroy(aux%aux_vars(iaux))
  enddo  
  do iaux = 1, aux%num_aux_bc
    call RTAuxVarDestroy(aux%aux_vars_bc(iaux))
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
