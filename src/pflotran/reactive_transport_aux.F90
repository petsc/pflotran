module Reactive_Transport_Aux_module

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
#ifdef REVISED_TRANSPORT
  use Matrix_Block_Aux_module
#endif  

  implicit none
  
  private 

#include "definitions.h"
 
  type, public :: reactive_transport_auxvar_type
    ! molality
    PetscReal, pointer :: pri_molal(:)     ! mol/kg water
    
    ! phase dependent totals
    PetscReal, pointer :: total(:,:)       ! mol solute/L water
#ifdef REVISED_TRANSPORT    
    type(matrix_block_auxvar_type), pointer :: aqueous
#else
    PetscReal, pointer :: dtotal(:,:,:)    ! kg water/m^3 water
#endif    
    ! sorbed totals
    PetscReal, pointer :: total_sorb_eq(:)    ! mol/m^3 bulk
    PetscReal, pointer :: dtotal_sorb_eq(:,:) ! kg water/m^3 bulk
    
    ! aqueous species
    ! aqueous complexes
    PetscReal, pointer :: sec_molal(:)
    PetscReal, pointer :: gas_molal(:)
    
    ! sorption reactions
    ! PetscReal, pointer :: kinionx_molfrac(:)
    PetscReal, pointer :: kinsrfcplx_conc(:) ! S_{i\alpha}^k
    PetscReal, pointer :: kinsrfcplx_conc_kp1(:) ! S_{i\alpha}^k+1
    PetscReal, pointer :: kinsrfcplx_free_site_conc(:)  ! S_\alpha
    PetscReal, pointer :: eqsrfcplx_conc(:)
    PetscReal, pointer :: eqsrfcplx_free_site_conc(:)
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
!   PetscReal :: act_h2o
    PetscReal, pointer :: pri_act_coef(:)
    PetscReal, pointer :: sec_act_coef(:)

    PetscReal :: ln_act_h2o
    
    PetscReal, pointer :: mass_balance(:,:)
    PetscReal, pointer :: mass_balance_delta(:,:)
    
    PetscReal, pointer :: kinmr_total_sorb(:,:)

#ifdef REVISED_TRANSPORT
    type(colloid_auxvar_type), pointer :: colloid
#endif
    
  end type reactive_transport_auxvar_type

  type, public :: reactive_transport_param_type
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: nimcomp
    PetscInt :: ncolcomp
    PetscInt :: offset_coll_sorb
    PetscReal :: dispersivity
    PetscReal, pointer :: diffusion_coefficient(:)
  end type reactive_transport_param_type

  ! Colloids
  type, public :: colloid_auxvar_type
    PetscReal, pointer :: total(:)
    type(matrix_block_auxvar_type), pointer :: dRj_dSic
    type(matrix_block_auxvar_type), pointer :: dRic_dCj
    type(matrix_block_auxvar_type), pointer :: dRic_dSic
  end type colloid_auxvar_type
  
  type, public :: colloid_param_type
    PetscInt :: num_colloids
    PetscInt :: num_colloid_comp
  end type colloid_param_type

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
  aux%rt_parameter%ncomp = 0
  aux%rt_parameter%naqcomp = 0
  aux%rt_parameter%nimcomp = 0
  aux%rt_parameter%ncolcomp = 0
  aux%rt_parameter%offset_coll_sorb = 0
  
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
  
  allocate(aux_var%pri_molal(reaction%naqcomp))
  aux_var%pri_molal = 0.d0

  allocate(aux_var%total(reaction%naqcomp,option%nphase))
  aux_var%total = 0.d0
#ifdef REVISED_TRANSPORT 
  aux_var%aqueous => MatrixBlockAuxVarCreate(option)
  call MatrixBlockAuxVarInit(aux_var%aqueous,reaction%naqcomp, &
                             reaction%naqcomp,option%nphase,option)
#else  
  allocate(aux_var%dtotal(reaction%naqcomp,reaction%naqcomp,option%nphase))
  aux_var%dtotal = 0.d0
#endif
  
  if (reaction%neqcplx > 0) then
    allocate(aux_var%sec_molal(reaction%neqcplx))
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

  if (reaction%neqsorb > 0) then  
    allocate(aux_var%total_sorb_eq(reaction%ncomp))
    aux_var%total_sorb_eq = 0.d0
    if (reaction%kinmr_nrate <= 0) then
      allocate(aux_var%dtotal_sorb_eq(reaction%ncomp,reaction%ncomp))
      aux_var%dtotal_sorb_eq = 0.d0
    else
      nullify(aux_var%dtotal_sorb_eq)
    endif
  else
    nullify(aux_var%total_sorb_eq)
    nullify(aux_var%dtotal_sorb_eq)
  endif    
  
  if (reaction%neqsrfcplxrxn > 0) then
    allocate(aux_var%eqsrfcplx_conc(reaction%neqsrfcplx))
    aux_var%eqsrfcplx_conc = 0.d0
    
    allocate(aux_var%eqsrfcplx_free_site_conc(reaction%neqsrfcplxrxn))
    aux_var%eqsrfcplx_free_site_conc = 1.d-9 ! initialize to guess
    
!   allocate(aux_var%eqsurf_site_density(reaction%neqsrfcplxrxn))
!   aux_var%eqsurf_site_density = 0.d0
  else
    nullify(aux_var%eqsrfcplx_conc)
    nullify(aux_var%eqsrfcplx_free_site_conc)
!   nullify(aux_var%eqsurf_site_density)
  endif
  
  if (reaction%nkinsrfcplxrxn > 0) then
    allocate(aux_var%kinsrfcplx_conc(reaction%nkinsrfcplx))
    aux_var%kinsrfcplx_conc = 0.d0

    allocate(aux_var%kinsrfcplx_conc_kp1(reaction%nkinsrfcplx))
    aux_var%kinsrfcplx_conc_kp1 = 0.d0
    
    allocate(aux_var%kinsrfcplx_free_site_conc(reaction%nkinsrfcplxrxn))
    aux_var%kinsrfcplx_free_site_conc = 0.d0 ! initialize to guess
    
!   allocate(aux_var%kinsurf_site_density(reaction%nkinsrfcplxrxn))
!   aux_var%kinsurf_site_density = 0.d0
  else
    nullify(aux_var%kinsrfcplx_conc)
    nullify(aux_var%kinsrfcplx_conc_kp1)
    nullify(aux_var%kinsrfcplx_free_site_conc)
!   nullify(aux_var%kinsurf_site_density)
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
  
  allocate(aux_var%pri_act_coef(reaction%ncomp))
  aux_var%pri_act_coef = 1.d0
  if (reaction%neqcplx > 0) then
    allocate(aux_var%sec_act_coef(reaction%neqcplx))
    aux_var%sec_act_coef = 1.d0
  else
    nullify(aux_var%sec_act_coef)
  endif

! initialize ln activity H2O
  aux_var%ln_act_h2o = 0.d0
  
  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(aux_var%mass_balance(reaction%ncomp,option%nphase))
    aux_var%mass_balance = 0.d0
    allocate(aux_var%mass_balance_delta(reaction%ncomp,option%nphase))
    aux_var%mass_balance_delta = 0.d0
  else
    nullify(aux_var%mass_balance)
    nullify(aux_var%mass_balance_delta)
  endif
  
  if (reaction%kinmr_nrate > 0) then
    allocate(aux_var%kinmr_total_sorb(reaction%ncomp,reaction%kinmr_nrate))
    aux_var%kinmr_total_sorb = 0.d0
  else
    nullify(aux_var%kinmr_total_sorb)
  endif

#ifdef REVISED_TRANSPORT
  if (reaction%ncolcomp > 0) then
    allocate(aux_var%colloid)
    allocate(aux_var%colloid%total(reaction%ncolcomp))
    ! dRj/dSic
    aux_var%colloid%dRj_dSic => MatrixBlockAuxVarCreate(option)
    call MatrixBlockAuxVarInit(aux_var%colloid%dRj_dSic,reaction%naqcomp, &
                               reaction%ncolcomp,ONE_INTEGER,option)
    ! dRic/dCj
    aux_var%colloid%dRic_dCj => MatrixBlockAuxVarCreate(option)
    call MatrixBlockAuxVarInit(aux_var%colloid%dRic_dCj,reaction%ncolcomp, &
                               reaction%naqcomp,ONE_INTEGER,option)
    ! dRic/dSic
    aux_var%colloid%dRic_dSic => MatrixBlockAuxVarCreate(option)
    call MatrixBlockAuxVarInit(aux_var%colloid%dRic_dSic,reaction%ncolcomp, &
                               reaction%ncolcomp,ONE_INTEGER,option)
  else
    nullify(aux_var%colloid)
  endif
#endif
  
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

#ifdef REVISED_TRANSPORT 
  call MatrixBlockAuxVarCopy(aux_var%aqueous,aux_var2%aqueous,option)
#else  
  aux_var%dtotal = aux_var2%dtotal
#endif
  
  if (associated(aux_var%sec_molal)) &
    aux_var%sec_molal = aux_var2%sec_molal
  if (associated(aux_var%total_sorb_eq)) then  
    aux_var%total_sorb_eq = aux_var2%total_sorb_eq
    aux_var%dtotal_sorb_eq = aux_var2%dtotal_sorb_eq
  endif
  
  if (associated(aux_var%gas_molal)) &
    aux_var%gas_molal = aux_var2%gas_molal
  
  if (associated(aux_var%eqsrfcplx_conc)) then
    aux_var%eqsrfcplx_conc = aux_var2%eqsrfcplx_conc
    aux_var%eqsrfcplx_free_site_conc = aux_var2%eqsrfcplx_free_site_conc
  endif
  
  if (associated(aux_var%kinsrfcplx_conc)) then
    aux_var%kinsrfcplx_conc = aux_var2%kinsrfcplx_conc
    aux_var%kinsrfcplx_conc_kp1 = aux_var2%kinsrfcplx_conc_kp1
    aux_var%kinsrfcplx_free_site_conc = aux_var2%kinsrfcplx_free_site_conc
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
  
  aux_var%ln_act_h2o = aux_var2%ln_act_h2o
  
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
  endif

#ifdef REVISED_TRANSPORT 
  if (associated(aux_var%colloid)) then
    aux_var%colloid%total = aux_var2%colloid%total
    ! dRj/dCj
    call MatrixBlockAuxVarCopy(aux_var%colloid%dRj_dSic, &
                               aux_var2%colloid%dRj_dSic,option)
    ! dRic/dCj
    call MatrixBlockAuxVarCopy(aux_var%colloid%dRic_dCj, &
                               aux_var2%colloid%dRic_dCj,option)
    ! dRic/dSic
    call MatrixBlockAuxVarCopy(aux_var%colloid%dRic_dSic, &
                               aux_var2%colloid%dRic_dSic,option)
  endif
#endif

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

#ifdef REVISED_TRANSPORT
  call MatrixBlockAuxVarDestroy(aux_var%aqueous)
#else  
  if (associated(aux_var%dtotal))deallocate(aux_var%dtotal)
  nullify(aux_var%dtotal)
#endif  
  
  if (associated(aux_var%sec_molal))deallocate(aux_var%sec_molal)
  nullify(aux_var%sec_molal)
  
  if (associated(aux_var%gas_molal))deallocate(aux_var%gas_molal)
  nullify(aux_var%gas_molal)
  
  if (associated(aux_var%total_sorb_eq)) deallocate(aux_var%total_sorb_eq)
  nullify(aux_var%total_sorb_eq)
  if (associated(aux_var%dtotal_sorb_eq))deallocate(aux_var%dtotal_sorb_eq)
  nullify(aux_var%dtotal_sorb_eq)

  if (associated(aux_var%eqsrfcplx_conc)) deallocate(aux_var%eqsrfcplx_conc)
  nullify(aux_var%eqsrfcplx_conc)
  if (associated(aux_var%eqsrfcplx_free_site_conc)) &
    deallocate(aux_var%eqsrfcplx_free_site_conc)
  nullify(aux_var%eqsrfcplx_free_site_conc)
  
  if (associated(aux_var%kinsrfcplx_conc)) deallocate(aux_var%kinsrfcplx_conc)
  nullify(aux_var%kinsrfcplx_conc)
  
  if (associated(aux_var%kinsrfcplx_conc_kp1)) deallocate(aux_var%kinsrfcplx_conc_kp1)
  nullify(aux_var%kinsrfcplx_conc_kp1)
  
  if (associated(aux_var%kinsrfcplx_free_site_conc)) &
    deallocate(aux_var%kinsrfcplx_free_site_conc)
  nullify(aux_var%kinsrfcplx_free_site_conc)
  
  if (associated(aux_var%eqionx_ref_cation_sorbed_conc)) &
    deallocate(aux_var%eqionx_ref_cation_sorbed_conc)
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
  
#ifdef REVISED_TRANSPORT
  if (associated(aux_var%colloid)) then
    if (associated(aux_var%colloid%total)) deallocate(aux_var%colloid%total)
    nullify(aux_var%colloid%total)
    ! dRj/dSic
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRj_dSic)
    ! dRic/dCj
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRic_dCj)
    ! dRic/dSic
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRic_dSic)
  endif
#endif    
  
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
