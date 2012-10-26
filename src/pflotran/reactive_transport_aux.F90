module Reactive_Transport_Aux_module

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  use Matrix_Block_Aux_module
  use Secondary_Continuum_module

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
    
  end type reactive_transport_auxvar_type

  ! START CHUNKED!!!!!
  type, public :: react_tran_auxvar_chunk_type
  
    PetscReal, pointer :: den(:,:,:)
    PetscReal, pointer :: temp(:,:,:)
    PetscReal, pointer :: sat(:,:,:)
    PetscReal, pointer :: vol(:,:)
    PetscReal, pointer :: por(:,:)
    
#ifdef CHUAN_CO2
    PetscReal, pointer :: pres(:,:,:)
    PetscReal, pointer :: xmass(:,:,:)
    PetscReal, pointer :: fugacoeff(:,:,:)
#endif    
  
    ! molality
    PetscReal, pointer :: pri_molal(:,:,:)     ! mol/kg water
    PetscReal, pointer :: ln_pri_molal(:,:,:)
    
    ! phase dependent totals
    PetscReal, pointer :: total(:,:,:,:)       ! mol solute/L water
    PetscReal, pointer :: dtotal(:,:,:,:,:)

    ! sorbed totals
    PetscReal, pointer :: total_sorb_eq(:,:,:)    ! mol/m^3 bulk
    PetscReal, pointer :: dtotal_sorb_eq(:,:,:,:) ! kg water/m^3 bulk
    
    ! aqueous species
    ! aqueous complexes
    PetscReal, pointer :: sec_molal(:,:,:)
    PetscReal, pointer :: gas_molal(:,:,:)
    
    PetscReal, pointer :: eqsrfcplx_conc(:,:,:)
    PetscReal, pointer :: eqsrfcplx_free_site_conc(:,:,:)

    ! mineral reactions
!    PetscReal, pointer :: mnrl_volfrac0(:,:,:)
    PetscReal, pointer :: mnrl_volfrac(:,:,:)
!    PetscReal, pointer :: mnrl_area0(:,:,:)
    PetscReal, pointer :: mnrl_area(:,:,:)
    PetscReal, pointer :: mnrl_rate(:,:,:)
    
    ! activity coefficients
!   PetscReal :: act_h2o
    PetscReal, pointer :: pri_act_coef(:,:,:)
    PetscReal, pointer :: sec_act_coef(:,:,:)
    
    PetscReal, pointer :: ln_act_h2o(:,:)

  end type react_tran_auxvar_chunk_type
  ! END CHUNKED!!!!!

  type, public :: reactive_transport_param_type
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: nimcomp
    PetscInt :: ncoll
    PetscInt :: ncollcomp
    PetscInt :: offset_aq
    PetscInt :: offset_coll
    PetscInt :: offset_collcomp
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
    type(sec_transport_type), pointer :: sec_transport_vars(:)
#ifdef CHUNK
    type(react_tran_auxvar_chunk_type), pointer :: aux_var_chunk
#endif
  end type reactive_transport_type

  interface RTAuxVarDestroy
    module procedure RTAuxVarSingleDestroy
    module procedure RTAuxVarArrayDestroy
  end interface RTAuxVarDestroy
  
  public :: RTAuxCreate, RTAuxDestroy, &
            RTAuxVarInit, RTAuxVarCopy, RTAuxVarDestroy, &
            RTAuxVarChunkDestroy, RTAuxVarStrip, &
            RTSecTransportAuxVarCompute
            
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
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)
#ifdef CHUNK
  nullify(aux%aux_var_chunk)
#endif  
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
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
  nullify(aux%sec_transport_vars)
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

  if (associated(aux_var%mass_balance) .and. &
      associated(aux_var2%mass_balance)) then
    aux_var%mass_balance = aux_var2%mass_balance
    aux_var%mass_balance_delta = aux_var2%mass_balance_delta
  endif

  if (associated(aux_var%kinmr_total_sorb) .and. &
      associated(aux_var2%kinmr_total_sorb)) then
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

end subroutine RTAuxVarCopy

! ************************************************************************** !
! 
! RTSecTransportAuxVarCompute: Computes secondary auxillary variables in each
!                              grid cell for transport only
! author: Satish Karra
! Date: 10/8/12
!
! ************************************************************************** !
subroutine RTSecTransportAuxVarCompute(sec_transport_vars,aux_var, &
                                       global_aux_var,reaction, &
                                       diffusion_coefficient,porosity, &
                                       option)

  use Option_module 
  use Global_Aux_module
  use Reaction_Aux_module  
  
  implicit none
  
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(sec_transport_vars%ncells)
  PetscReal :: coeff_diag(sec_transport_vars%ncells)
  PetscReal :: coeff_right(sec_transport_vars%ncells)
  PetscReal :: rhs(sec_transport_vars%ncells)
  PetscReal :: sec_conc(sec_transport_vars%ncells)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, diffusion_coefficient, porosity
  PetscReal :: conc_primary_node
  PetscReal :: m
  PetscReal :: kin_mnrl_rate
  PetscReal :: mnrl_area, mnrl_molar_vol
  PetscReal :: equil_conc, Im(sec_transport_vars%ncells)
  PetscReal :: sec_mnrl_volfrac(sec_transport_vars%ncells)
  PetscInt :: sec_zeta(sec_transport_vars%ncells)
  PetscReal :: diag_react, rhs_react
  PetscReal, parameter :: rgas = 8.3144621d-3
  PetscReal :: arrhenius_factor

  
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  diag_react = 0.d0
  rhs_react = 0.d0
  Im = 0.d0
  
  if (reaction%naqcomp > 1 .or. reaction%mineral%nkinmnrl > 1) then
    option%io_buffer = 'Currently only single component system with ' // &
                       'multiple continuum is implemented'
    call printErrMsg(option)
  endif

  conc_primary_node = aux_var%total(1,1)                             ! in mol/L 
  sec_mnrl_volfrac = sec_transport_vars%sec_mnrl_volfrac             ! dimensionless
  mnrl_area = sec_transport_vars%sec_mnrl_area                       ! in 1/cm
  
  if (reaction%mineral%nkinmnrl > 0) then
    kin_mnrl_rate = reaction%mineral%kinmnrl_rate(1)                 ! in mol/cm^2/s
    ! Arrhenius factor
    arrhenius_factor = 1.d0
    if (reaction%mineral%kinmnrl_activation_energy(1) > 0.d0) then
      arrhenius_factor = exp(reaction%mineral%kinmnrl_activation_energy(1)/rgas &
          *(1.d0/(25.d0+273.15d0)-1.d0/(global_aux_var%temp(1)+273.15d0)))
    endif    
    kin_mnrl_rate = kin_mnrl_rate*arrhenius_factor
    equil_conc = (10.d0)**(reaction%mineral%mnrl_logK(1))            ! in mol/L
    mnrl_molar_vol = reaction%mineral%kinmnrl_molar_vol(1)           ! in m^3/mol
    diag_react = kin_mnrl_rate/equil_conc*mnrl_area*option%tran_dt/porosity*1.d3
    rhs_react = diag_react*equil_conc                                ! in mol/L
  endif
  
  alpha = diffusion_coefficient*option%tran_dt 
  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0 &
                    + diag_react*sec_zeta(i)
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0 &
                  + diag_react*sec_zeta(1)
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0 + diag_react*sec_zeta(ngcells)
                        
  do i = 1, ngcells
    rhs(i) = sec_transport_vars%sec_conc(i) + rhs_react*sec_zeta(i) ! secondary continuum values from previous time step
  enddo
  
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 conc_primary_node 
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate concentration in the secondary continuum
  sec_conc(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_conc(i) = (rhs(i) - coeff_right(i)*sec_conc(i+1))/coeff_diag(i)
  enddo

 ! print *,'conc_dcdm= ',(sec_conc(i),i=1,ngcells)
 
   do i = 1, ngcells
    Im(i) = kin_mnrl_rate*mnrl_area*(sec_conc(i)/equil_conc - 1.d0) ! in mol/cm^3/s
    if (Im(i) > 0.d0) then 
      sec_mnrl_volfrac(i) = sec_mnrl_volfrac(i) + option%tran_dt*1.d6* &
                            mnrl_molar_vol*Im(i)
      sec_zeta(i) = 1
    else
      if (sec_mnrl_volfrac(i) > 0.d0) then
        sec_mnrl_volfrac(i) = sec_mnrl_volfrac(i) + option%tran_dt*1.d6* &
                              mnrl_molar_vol*Im(i)
        sec_zeta(i) = 1
      else
        Im(i) = 0.d0
        sec_zeta(i) = 0
      endif
    endif
    if (sec_mnrl_volfrac(i) < 0.d0) then
      sec_mnrl_volfrac(i) = 0.d0
      sec_zeta(i) = 0
    endif
  enddo

  sec_transport_vars%sec_conc = sec_conc
  sec_transport_vars%sec_mnrl_volfrac = sec_mnrl_volfrac
  sec_transport_vars%sec_zeta = sec_zeta

end subroutine RTSecTransportAuxVarCompute


! ************************************************************************** !
!
! RTAuxVarChunkDestroy: Deallocates a reactive transport auxilliary object
! author: Glenn Hammond
! date: 01/31/11
!
! ************************************************************************** !
subroutine RTAuxVarChunkDestroy(auxvar)

  implicit none

  type(react_tran_auxvar_chunk_type), pointer :: auxvar
  
    ! for global auxvar
  if (associated(auxvar%den)) deallocate(auxvar%den)
  nullify(auxvar%den)
  if (associated(auxvar%temp)) deallocate(auxvar%temp)
  nullify(auxvar%temp)
  if (associated(auxvar%sat)) deallocate(auxvar%sat)
  nullify(auxvar%sat)
  if (associated(auxvar%vol)) deallocate(auxvar%vol)
  nullify(auxvar%vol)
  if (associated(auxvar%por)) deallocate(auxvar%por)
  nullify(auxvar%por)

#ifdef CHUAN_CO2
  if (associated(auxvar%pres)) deallocate(auxvar%pres)
  nullify(auxvar%pres)
  if (associated(auxvar%xmass)) deallocate(auxvar%xmass)
  nullify(auxvar%xmass)
  if (associated(auxvar%fugacoeff)) deallocate(auxvar%fugacoeff)
  nullify(auxvar%fugacoeff)
#endif
  
  if (associated(auxvar%pri_molal)) deallocate(auxvar%pri_molal)
  nullify(auxvar%pri_molal)

  if (associated(auxvar%total)) deallocate(auxvar%total)
  nullify(auxvar%total)

  if (associated(auxvar%total)) deallocate(auxvar%dtotal)
  nullify(auxvar%dtotal)

  if (associated(auxvar%sec_molal))deallocate(auxvar%sec_molal)
  nullify(auxvar%sec_molal)
  
  if (associated(auxvar%gas_molal))deallocate(auxvar%gas_molal)
  nullify(auxvar%gas_molal)
  
  if (associated(auxvar%total_sorb_eq)) deallocate(auxvar%total_sorb_eq)
  nullify(auxvar%total_sorb_eq)
  if (associated(auxvar%dtotal_sorb_eq))deallocate(auxvar%dtotal_sorb_eq)
  nullify(auxvar%dtotal_sorb_eq)

  if (associated(auxvar%eqsrfcplx_conc)) deallocate(auxvar%eqsrfcplx_conc)
  nullify(auxvar%eqsrfcplx_conc)
  if (associated(auxvar%eqsrfcplx_free_site_conc)) &
    deallocate(auxvar%eqsrfcplx_free_site_conc)
  nullify(auxvar%eqsrfcplx_free_site_conc)
  
  if (associated(auxvar%mnrl_volfrac))deallocate(auxvar%mnrl_volfrac)
  nullify(auxvar%mnrl_volfrac)
  if (associated(auxvar%mnrl_area))deallocate(auxvar%mnrl_area)
  nullify(auxvar%mnrl_area)
  if (associated(auxvar%mnrl_rate))deallocate(auxvar%mnrl_rate)
  nullify(auxvar%mnrl_rate)
  
  if (associated(auxvar%pri_act_coef))deallocate(auxvar%pri_act_coef)
  nullify(auxvar%pri_act_coef)
  if (associated(auxvar%sec_act_coef))deallocate(auxvar%sec_act_coef)
  nullify(auxvar%sec_act_coef)

  if (associated(auxvar%ln_act_h2o))deallocate(auxvar%ln_act_h2o)
  nullify(auxvar%ln_act_h2o)
  
  deallocate(auxvar)
  nullify(auxvar)

end subroutine RTAuxVarChunkDestroy

! ************************************************************************** !
!
! RTAuxVarSingleDestroy: Deallocates a mode auxilliary object
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
! RTAuxVarArrayDestroy: Deallocates a mode auxilliary object
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
! RTAuxVarStrip: Deallocates all members of single auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTAuxVarStrip(aux_var)

  implicit none

  type(reactive_transport_auxvar_type) :: aux_var
  
  if (associated(aux_var%pri_molal)) deallocate(aux_var%pri_molal)
  nullify(aux_var%pri_molal)

  if (associated(aux_var%total)) deallocate(aux_var%total)
  nullify(aux_var%total)

  call MatrixBlockAuxVarDestroy(aux_var%aqueous)
  
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
  if (associated(aux_var%srfcplxrxn_free_site_conc)) &
    deallocate(aux_var%srfcplxrxn_free_site_conc)
  nullify(aux_var%srfcplxrxn_free_site_conc)
  
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
  
  if (associated(aux_var%colloid)) then
    if (associated(aux_var%colloid%conc_mob)) deallocate(aux_var%colloid%conc_mob)
    nullify(aux_var%colloid%conc_mob)
    if (associated(aux_var%colloid%conc_imb)) deallocate(aux_var%colloid%conc_imb)
    nullify(aux_var%colloid%conc_imb)
    if (associated(aux_var%colloid%total_eq_mob)) deallocate(aux_var%colloid%total_eq_mob)
    nullify(aux_var%colloid%total_eq_mob)
    if (associated(aux_var%colloid%total_kin)) deallocate(aux_var%colloid%total_kin)
    nullify(aux_var%colloid%total_kin)
    ! dRj/dCj
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRj_dCj)
    ! dRj/dSic
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRj_dSic)
    ! dRic/dCj
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRic_dCj)
    ! dRic/dSic
    call MatrixBlockAuxVarDestroy(aux_var%colloid%dRic_dSic)
  endif
  
end subroutine RTAuxVarStrip

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
  
  call RTAuxVarDestroy(aux%aux_vars)
  call RTAuxVarDestroy(aux%aux_vars_bc)
  call RTAuxVarDestroy(aux%aux_vars_ss)
#ifdef CHUNK
  if (associated(aux%aux_var_chunk)) then
    call RTAuxVarChunkDestroy(aux%aux_var_chunk)
  endif
#endif  
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%rt_parameter)) then
    if (associated(aux%rt_parameter%diffusion_coefficient)) &
      deallocate(aux%rt_parameter%diffusion_coefficient)
    nullify(aux%rt_parameter%diffusion_coefficient)
    if (associated(aux%rt_parameter%diffusion_activation_energy)) &
      deallocate(aux%rt_parameter%diffusion_activation_energy)
    nullify(aux%rt_parameter%diffusion_activation_energy)
    if (associated(aux%rt_parameter%pri_spec_to_coll_spec)) &
      deallocate(aux%rt_parameter%pri_spec_to_coll_spec)
    nullify(aux%rt_parameter%pri_spec_to_coll_spec)
    if (associated(aux%rt_parameter%coll_spec_to_pri_spec)) &
      deallocate(aux%rt_parameter%coll_spec_to_pri_spec)
    nullify(aux%rt_parameter%coll_spec_to_pri_spec)
    deallocate(aux%rt_parameter)
  endif
  nullify(aux%rt_parameter)
  
  if (associated(aux%sec_transport_vars)) deallocate (aux%sec_transport_vars)
  nullify (aux%sec_transport_vars)
  
  deallocate(aux)
  nullify(aux)  

  end subroutine RTAuxDestroy

end module Reactive_Transport_Aux_module
