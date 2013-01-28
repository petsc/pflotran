module Reactive_Transport_Aux_module

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  use Matrix_Block_Aux_module

  implicit none
  
  private 

#include "definitions.h"
 
  type, public :: reactive_transport_auxvar_type
    ! molality
    PetscReal, pointer :: pri_molal(:)     ! mol/kg water
    
    ! phase dependent totals
    PetscReal, pointer :: total(:,:)       ! mol solute/L water
    type(matrix_block_auxvar_type), pointer :: aqueous

    ! sorbed totals
    PetscReal, pointer :: total_sorb_eq(:)    ! mol/m^3 bulk
    PetscReal, pointer :: dtotal_sorb_eq(:,:) ! kg water/m^3 bulk
    
    ! aqueous species
    ! aqueous complexes
    PetscReal, pointer :: sec_molal(:)
    PetscReal, pointer :: gas_molal(:)
    
    ! sorption reactions
    PetscReal, pointer :: srfcplxrxn_free_site_conc(:)
    PetscReal, pointer :: kinsrfcplx_conc(:,:) ! S_{i\alpha}^k
    PetscReal, pointer :: kinsrfcplx_conc_kp1(:,:) ! S_{i\alpha}^k+1
    PetscReal, pointer :: kinsrfcplx_free_site_conc(:)  ! S_\alpha
    PetscReal, pointer :: eqsrfcplx_conc(:)
    PetscReal, pointer :: eqionx_ref_cation_sorbed_conc(:)
    PetscReal, pointer :: eqionx_conc(:,:)
    
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
    
    PetscReal, pointer :: kinmr_total_sorb(:,:,:)

    type(colloid_auxvar_type), pointer :: colloid
    
    ! immobile species such as biomass
    PetscReal, pointer :: immobile(:)
    
  end type reactive_transport_auxvar_type

  type, public :: reactive_transport_param_type
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: nimcomp
    PetscInt :: ncoll
    PetscInt :: ncollcomp
    PetscInt :: offset_aqueous
    PetscInt :: offset_colloid
    PetscInt :: offset_collcomp
    PetscInt :: offset_immobile
    PetscInt, pointer :: pri_spec_to_coll_spec(:)
    PetscInt, pointer :: coll_spec_to_pri_spec(:)
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
#ifdef OS_STATISTICS
! use PetscReal for large counts
    PetscInt :: newton_call_count
    PetscReal :: sum_newton_call_count
    PetscInt :: newton_iterations
    PetscReal :: sum_newton_iterations
    PetscInt :: max_newton_iterations
    PetscInt :: overall_max_newton_iterations
#endif    

  end type reactive_transport_param_type

  ! Colloids
  type, public :: colloid_auxvar_type
    PetscReal, pointer :: conc_mob(:) ! mol/L water
    PetscReal, pointer :: conc_imb(:) ! mol/m^3 bulk
    PetscReal, pointer :: total_eq_mob(:) ! mol/L water
    PetscReal, pointer :: total_kin(:)
    type(matrix_block_auxvar_type), pointer :: dRj_dCj
    type(matrix_block_auxvar_type), pointer :: dRj_dSic
    type(matrix_block_auxvar_type), pointer :: dRic_dCj
    type(matrix_block_auxvar_type), pointer :: dRic_dSic
  end type colloid_auxvar_type
  
  type, public :: colloid_param_type
    PetscInt :: num_colloids
    PetscInt :: num_colloid_comp
  end type colloid_param_type

  type, public :: reactive_transport_type
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscInt :: n_zero_rows
    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    type(reactive_transport_param_type), pointer :: rt_parameter
    type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
    type(reactive_transport_auxvar_type), pointer :: aux_vars_bc(:)
    type(reactive_transport_auxvar_type), pointer :: aux_vars_ss(:)
  end type reactive_transport_type

  interface RTAuxVarDestroy
    module procedure RTAuxVarSingleDestroy
    module procedure RTAuxVarArrayDestroy
  end interface RTAuxVarDestroy
  
  public :: RTAuxCreate, RTAuxDestroy, &
            RTAuxVarInit, RTAuxVarCopy, RTAuxVarDestroy, &
            RTAuxVarStrip
            
contains


! ************************************************************************** !
!
! RTAuxCreate: Allocate and initialize auxiliary object
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
  aux%num_aux = 0      ! number of rt_auxvars objects for local and ghosted cells
  aux%num_aux_bc = 0   ! number of rt_auxvars objects for boundary connections
  aux%num_aux_ss = 0   ! number of rt_auxvars objects for source/sinks
  nullify(aux%aux_vars)      ! rt_auxvars for local and ghosted grid cells
  nullify(aux%aux_vars_bc)   ! rt_auxvars for boundary connections
  nullify(aux%aux_vars_ss)   ! rt_auxvars for source/sinks
  aux%n_zero_rows = 0    ! number of zeroed rows in Jacobian for inactive cells
  nullify(aux%zero_rows_local)  ! ids of zero rows in local, non-ghosted numbering
  nullify(aux%zero_rows_local_ghosted) ! ids of zero rows in ghosted numbering
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE

  allocate(aux%rt_parameter)
  allocate(aux%rt_parameter%diffusion_coefficient(option%nphase))
  allocate(aux%rt_parameter%diffusion_activation_energy(option%nphase))
  aux%rt_parameter%diffusion_coefficient = 1.d-9
  aux%rt_parameter%diffusion_activation_energy = 0.d0
  aux%rt_parameter%ncomp = 0
  aux%rt_parameter%naqcomp = 0
  aux%rt_parameter%nimcomp = 0
  aux%rt_parameter%ncoll = 0
  aux%rt_parameter%ncollcomp = 0
  aux%rt_parameter%offset_aqueous = 0
  aux%rt_parameter%offset_colloid = 0
  aux%rt_parameter%offset_collcomp = 0
  aux%rt_parameter%offset_immobile = 0
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
  RTAuxCreate => aux
  
end function RTAuxCreate

! ************************************************************************** !
!
! RTAuxVarInit: Initialize auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxVarInit(aux_var,reaction,option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type
  use Surface_Complexation_Aux_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(surface_complexation_type), pointer :: surface_complexation
  
  surface_complexation => reaction%surface_complexation
  
  allocate(aux_var%pri_molal(reaction%naqcomp))
  aux_var%pri_molal = 0.d0

  allocate(aux_var%total(reaction%naqcomp,option%nphase))
  aux_var%total = 0.d0
  aux_var%aqueous => MatrixBlockAuxVarCreate(option)
  call MatrixBlockAuxVarInit(aux_var%aqueous,reaction%naqcomp, &
                             reaction%naqcomp,option%nphase,option)
  
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
    allocate(aux_var%total_sorb_eq(reaction%naqcomp))
    aux_var%total_sorb_eq = 0.d0
    allocate(aux_var%dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp))
    aux_var%dtotal_sorb_eq = 0.d0
  else
    nullify(aux_var%total_sorb_eq)
    nullify(aux_var%dtotal_sorb_eq)
  endif    
  
  ! surface complexation
  nullify(aux_var%eqsrfcplx_conc)
  nullify(aux_var%srfcplxrxn_free_site_conc)
  nullify(aux_var%kinsrfcplx_conc)
  nullify(aux_var%kinsrfcplx_conc_kp1)
  nullify(aux_var%kinsrfcplx_free_site_conc)
  nullify(aux_var%kinmr_total_sorb)
  if (surface_complexation%nsrfcplxrxn > 0) then
    allocate(aux_var%srfcplxrxn_free_site_conc(surface_complexation%nsrfcplxrxn))
    aux_var%srfcplxrxn_free_site_conc = 1.d-9 ! initialize to guess
    if (surface_complexation%neqsrfcplxrxn > 0) then
      allocate(aux_var%eqsrfcplx_conc(surface_complexation%nsrfcplx))
      aux_var%eqsrfcplx_conc = 0.d0
    endif
    if (surface_complexation%nkinsrfcplxrxn > 0) then
      !geh: currently hardwired to only 1 reaction
      allocate(aux_var%kinsrfcplx_conc(surface_complexation%nkinsrfcplx,1))
      aux_var%kinsrfcplx_conc = 0.d0

      allocate(aux_var%kinsrfcplx_conc_kp1(surface_complexation%nkinsrfcplx,1))
      aux_var%kinsrfcplx_conc_kp1 = 0.d0
    endif
    if (surface_complexation%nkinmrsrfcplxrxn > 0) then
      ! the zeroth entry here stores the equilibrium concentration used in the 
      ! update
      ! the zeroth entry of kinmr_nrate holds the maximum number of rates
      ! prescribed in a multirate reaction...required for appropriate sizing
      allocate(aux_var%kinmr_total_sorb(reaction%naqcomp, &
                                        0:surface_complexation%kinmr_nrate(0), &
                                        surface_complexation%nkinmrsrfcplxrxn))
      aux_var%kinmr_total_sorb = 0.d0
    endif
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
  
  if (associated(reaction%mineral)) then
    if (reaction%mineral%nkinmnrl > 0) then
      allocate(aux_var%mnrl_volfrac0(reaction%mineral%nkinmnrl))
      aux_var%mnrl_volfrac0 = 0.d0
      allocate(aux_var%mnrl_volfrac(reaction%mineral%nkinmnrl))
      aux_var%mnrl_volfrac = 0.d0
      allocate(aux_var%mnrl_area0(reaction%mineral%nkinmnrl))
      aux_var%mnrl_area0 = 0.d0
      allocate(aux_var%mnrl_area(reaction%mineral%nkinmnrl))
      aux_var%mnrl_area = 0.d0
      allocate(aux_var%mnrl_rate(reaction%mineral%nkinmnrl))
      aux_var%mnrl_rate = 0.d0
    else
      nullify(aux_var%mnrl_volfrac0)
      nullify(aux_var%mnrl_volfrac)
      nullify(aux_var%mnrl_area0)
      nullify(aux_var%mnrl_area)
      nullify(aux_var%mnrl_rate)
    endif
  endif
  
  allocate(aux_var%pri_act_coef(reaction%naqcomp))
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
  
  if (reaction%ncollcomp > 0) then
    allocate(aux_var%colloid)
    allocate(aux_var%colloid%conc_mob(reaction%ncoll))
    allocate(aux_var%colloid%conc_imb(reaction%ncoll))
    allocate(aux_var%colloid%total_eq_mob(reaction%ncollcomp))
    allocate(aux_var%colloid%total_kin(reaction%ncollcomp))
    ! dRj/dCj
    aux_var%colloid%dRj_dCj => MatrixBlockAuxVarCreate(option)
    call MatrixBlockAuxVarInit(aux_var%colloid%dRj_dCj,reaction%naqcomp, &
                               reaction%naqcomp,ONE_INTEGER,option)
    ! dRj/dSic
    aux_var%colloid%dRj_dSic => MatrixBlockAuxVarCreate(option)
    call MatrixBlockAuxVarInit(aux_var%colloid%dRj_dSic,reaction%naqcomp, &
                               reaction%ncollcomp,ONE_INTEGER,option)
    ! dRic/dCj
    aux_var%colloid%dRic_dCj => MatrixBlockAuxVarCreate(option)
    call MatrixBlockAuxVarInit(aux_var%colloid%dRic_dCj,reaction%ncollcomp, &
                               reaction%naqcomp,ONE_INTEGER,option)
    ! dRic/dSic
    aux_var%colloid%dRic_dSic => MatrixBlockAuxVarCreate(option)
    call MatrixBlockAuxVarInit(aux_var%colloid%dRic_dSic,reaction%ncollcomp, &
                               reaction%ncollcomp,ONE_INTEGER,option)
  else
    nullify(aux_var%colloid)
  endif
  
  if (reaction%nimcomp > 0) then
    allocate(aux_var%immobile(reaction%nimcomp))
    aux_var%immobile = 0.d0
  else
    nullify(aux_var%immobile)
  endif
  
end subroutine RTAuxVarInit

! ************************************************************************** !
!
! RTAuxVarCopy: Copys an auxiliary object
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

  call MatrixBlockAuxVarCopy(aux_var%aqueous,aux_var2%aqueous,option)
  
  if (associated(aux_var%sec_molal)) &
    aux_var%sec_molal = aux_var2%sec_molal
  if (associated(aux_var%total_sorb_eq)) then  
    aux_var%total_sorb_eq = aux_var2%total_sorb_eq
  endif
  if (associated(aux_var%dtotal_sorb_eq)) then  
    aux_var%dtotal_sorb_eq = aux_var2%dtotal_sorb_eq
  endif
  
  if (associated(aux_var%gas_molal)) &
    aux_var%gas_molal = aux_var2%gas_molal
  
  if (associated(aux_var%srfcplxrxn_free_site_conc)) then
    aux_var%srfcplxrxn_free_site_conc = aux_var2%srfcplxrxn_free_site_conc
  endif
  
  if (associated(aux_var%eqsrfcplx_conc)) then
    aux_var%eqsrfcplx_conc = aux_var2%eqsrfcplx_conc
  endif
  
  if (associated(aux_var%kinsrfcplx_conc)) then
    aux_var%kinsrfcplx_conc = aux_var2%kinsrfcplx_conc
    aux_var%kinsrfcplx_conc_kp1 = aux_var2%kinsrfcplx_conc_kp1
    aux_var%kinsrfcplx_free_site_conc = aux_var2%kinsrfcplx_free_site_conc
  endif
  
  if (associated(aux_var%eqionx_ref_cation_sorbed_conc)) then
    aux_var%eqionx_ref_cation_sorbed_conc = &
      aux_var2%eqionx_ref_cation_sorbed_conc
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

  if (associated(aux_var%mass_balance)) then
    aux_var%mass_balance = aux_var2%mass_balance
    aux_var%mass_balance_delta = aux_var2%mass_balance_delta
  endif

  if (associated(aux_var%kinmr_total_sorb)) then
    aux_var%kinmr_total_sorb = aux_var2%kinmr_total_sorb
  endif

  if (associated(aux_var%colloid)) then
    aux_var%colloid%conc_mob = aux_var2%colloid%conc_mob
    aux_var%colloid%conc_imb = aux_var2%colloid%conc_imb
    aux_var%colloid%total_eq_mob = aux_var2%colloid%total_eq_mob
    aux_var%colloid%total_kin = aux_var2%colloid%total_kin
    ! dRj/dCj
    call MatrixBlockAuxVarCopy(aux_var%colloid%dRj_dCj, &
                               aux_var2%colloid%dRj_dCj,option)
    ! dRj/dSic
    call MatrixBlockAuxVarCopy(aux_var%colloid%dRj_dSic, &
                               aux_var2%colloid%dRj_dSic,option)
    ! dRic/dCj
    call MatrixBlockAuxVarCopy(aux_var%colloid%dRic_dCj, &
                               aux_var2%colloid%dRic_dCj,option)
    ! dRic/dSic
    call MatrixBlockAuxVarCopy(aux_var%colloid%dRic_dSic, &
                               aux_var2%colloid%dRic_dSic,option)
  endif

  if (associated(aux_var%immobile)) then
    aux_var%immobile = aux_var2%immobile
  endif
  
end subroutine RTAuxVarCopy

! ************************************************************************** !
!
! RTAuxVarSingleDestroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine RTAuxVarSingleDestroy(aux_var)

  implicit none

  type(reactive_transport_auxvar_type), pointer :: aux_var
  
  if (associated(aux_var)) then
    call RTAuxVarStrip(aux_var)
    deallocate(aux_var)
  endif
  nullify(aux_var)  

end subroutine RTAuxVarSingleDestroy
  
! ************************************************************************** !
!
! RTAuxVarArrayDestroy: Deallocates a mode auxiliary object
! author: Glenn Hammond
! date: 01/10/12
!
! ************************************************************************** !
subroutine RTAuxVarArrayDestroy(aux_vars)

  implicit none

  type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
  
  PetscInt :: iaux
  
  if (associated(aux_vars)) then
    do iaux = 1, size(aux_vars)
      call RTAuxVarStrip(aux_vars(iaux))
    enddo  
    deallocate(aux_vars)
  endif
  nullify(aux_vars)

end subroutine RTAuxVarArrayDestroy
  
! ************************************************************************** !
!
! RTAuxVarStrip: Deallocates all members of single auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxVarStrip(aux_var)

  use Utility_module, only: DeallocateArray

  implicit none

  type(reactive_transport_auxvar_type) :: aux_var
  
  call DeallocateArray(aux_var%pri_molal)
  call DeallocateArray(aux_var%total)
  
  call MatrixBlockAuxVarDestroy(aux_var%aqueous)

  call DeallocateArray(aux_var%sec_molal)
  call DeallocateArray(aux_var%gas_molal)
  call DeallocateArray(aux_var%total_sorb_eq)
  call DeallocateArray(aux_var%dtotal_sorb_eq)
  
  call DeallocateArray(aux_var%eqsrfcplx_conc)
  call DeallocateArray(aux_var%srfcplxrxn_free_site_conc)
  call DeallocateArray(aux_var%kinsrfcplx_conc)
  call DeallocateArray(aux_var%kinsrfcplx_conc_kp1)
  call DeallocateArray(aux_var%kinsrfcplx_free_site_conc)
  
  call DeallocateArray(aux_var%eqionx_ref_cation_sorbed_conc)
  call DeallocateArray(aux_var%eqionx_conc)
  
  call DeallocateArray(aux_var%mnrl_volfrac0)
  call DeallocateArray(aux_var%mnrl_volfrac)
  call DeallocateArray(aux_var%mnrl_area0)
  call DeallocateArray(aux_var%mnrl_area)
  call DeallocateArray(aux_var%mnrl_rate)
  
  call DeallocateArray(aux_var%pri_act_coef)
  call DeallocateArray(aux_var%sec_act_coef)
  
  call DeallocateArray(aux_var%mass_balance)
  call DeallocateArray(aux_var%mass_balance_delta)
  
  call DeallocateArray(aux_var%kinmr_total_sorb)
  
  if (associated(aux_var%colloid)) then
    call DeallocateArray(aux_var%colloid%conc_mob)
    call DeallocateArray(aux_var%colloid%conc_imb)
    call DeallocateArray(aux_var%colloid%total_eq_mob)
    call DeallocateArray(aux_var%colloid%total_kin)
    ! dRj/dCj
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRj_dCj)
    ! dRj/dSic
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRj_dSic)
    ! dRic/dCj
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRic_dCj)
    ! dRic/dSic
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRic_dSic)
    deallocate(aux_var%colloid)
    nullify(aux_var%colloid)
  endif
  
  call DeallocateArray(aux_var%immobile)
  
end subroutine RTAuxVarStrip

! ************************************************************************** !
!
! RTAuxDestroy: Deallocates a reactive transport auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxDestroy(aux)

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(reactive_transport_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  call RTAuxVarDestroy(aux%aux_vars)
  call RTAuxVarDestroy(aux%aux_vars_bc)
  call RTAuxVarDestroy(aux%aux_vars_ss)
  call DeallocateArray(aux%zero_rows_local)
  call DeallocateArray(aux%zero_rows_local_ghosted)

  if (associated(aux%rt_parameter)) then
    call DeallocateArray(aux%rt_parameter%diffusion_coefficient)
    call DeallocateArray(aux%rt_parameter%diffusion_activation_energy)
    call DeallocateArray(aux%rt_parameter%pri_spec_to_coll_spec)
    call DeallocateArray(aux%rt_parameter%coll_spec_to_pri_spec)
    deallocate(aux%rt_parameter)
  endif
  nullify(aux%rt_parameter)

  deallocate(aux)
  nullify(aux)

  end subroutine RTAuxDestroy

end module Reactive_Transport_Aux_module
