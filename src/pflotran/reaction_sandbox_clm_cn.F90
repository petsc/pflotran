module Reaction_Sandbox_CLM_CN_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clm_cn_type
    PetscInt :: nrxn
    PetscInt :: npool
    PetscReal, pointer :: CN_ratio(:)
    PetscReal, pointer :: rate_constant(:)
    PetscReal, pointer :: respiration_fraction(:)
    PetscReal, pointer :: inhibition_constant(:)
    PetscInt, pointer :: upstream_pool_id(:)
    PetscInt, pointer :: downstream_pool_id(:)
    PetscInt, pointer :: inhibitor_id(:)
    PetscInt, pointer :: pool_id_to_species_id(:,:)
    PetscInt :: C_species_id
    PetscInt :: N_species_id
    type(pool_type), pointer :: pools
    type(clm_cn_reaction_type), pointer :: reactions
  contains
    procedure, public :: Init => CLM_CN_Init
    procedure, public :: ReadInput => CLM_CN_Read
    procedure, public :: SkipInput => CLM_CN_ReadSkipBlock
    procedure, public :: Evaluate => CLM_CN_React
    procedure, public :: Destroy => CLM_CN_Destroy
  end type reaction_sandbox_clm_cn_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: CN_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: clm_cn_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    character(len=MAXWORDLENGTH) :: downstream_pool_name
    character(len=MAXWORDLENGTH) :: inhibitor_name
    PetscReal :: rate_constant
    PetscReal :: respiration_fraction
    PetscReal :: inhibition_constant
    type(clm_cn_reaction_type), pointer :: next
  end type clm_cn_reaction_type
  
  public :: CLM_CN_Create

contains

! ************************************************************************** !
!
! RSandboxInit: Initializes reaction sandbox at beginning of simulation
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
function CLM_CN_Create()

  implicit none
  
  type(reaction_sandbox_clm_cn_type), pointer :: CLM_CN_Create
  
  allocate(CLM_CN_Create)
  CLM_CN_Create%nrxn = 0
  CLM_CN_Create%npool = 0
  nullify(CLM_CN_Create%CN_ratio)
  nullify(CLM_CN_Create%rate_constant)
  nullify(CLM_CN_Create%respiration_fraction)
  nullify(CLM_CN_Create%inhibition_constant)
  nullify(CLM_CN_Create%pool_id_to_species_id)
  nullify(CLM_CN_Create%upstream_pool_id)
  nullify(CLM_CN_Create%downstream_pool_id)
  nullify(CLM_CN_Create%inhibitor_id)
  CLM_CN_Create%C_species_id = 0
  CLM_CN_Create%N_species_id = 0
  nullify(CLM_CN_Create%next)
  nullify(CLM_CN_Create%pools)
  nullify(CLM_CN_Create%reactions)

end function CLM_CN_Create

! ************************************************************************** !
!
! CLM_CN_Init: Initializes reaction sandbox at beginning of simulation
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
subroutine CLM_CN_Init(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type
  use Option_module
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction  
  
  call CLM_CN_Map(this,reaction,option)

end subroutine CLM_CN_Init

! ************************************************************************** !
!
! CLM_CN_Read: Reads input deck for reaction sandbox parameters
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
subroutine CLM_CN_Read(this,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: temp_real
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(clm_cn_reaction_type), pointer :: new_reaction, prev_reaction
  
  nullify(new_pool)
  nullify(prev_pool)
  
  nullify(new_reaction)
  nullify(prev_reaction)
  
  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CN')
    call StringToUpper(word)   

    select case(trim(word))
      case('POOLS')
        do
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   

          allocate(new_pool)
          new_pool%name = ''
          new_pool%CN_ratio = 0.d0
          nullify(new_pool%next)

          call InputReadWord(input,option,new_pool%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'pool name', &
            'CHEMISTRY,REACTION_SANDBOX,CLM_CN,POOLS')
          call InputReadDouble(input,option,new_pool%CN_ratio)
          if (InputError(input)) then
            new_pool%CN_ratio = -999.d0
          endif
          if (associated(this%pools)) then
            prev_pool%next => new_pool
          else
            this%pools => new_pool
          endif
          prev_pool => new_pool
          nullify(new_pool)
        enddo
      case('REACTION')
      
        allocate(new_reaction)
        new_reaction%upstream_pool_name = ''
        new_reaction%downstream_pool_name = ''
        new_reaction%inhibitor_name = ''
        new_reaction%rate_constant = -999.d0
        new_reaction%respiration_fraction = -999.d0
        new_reaction%inhibition_constant = 0.d0
        
        do 
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,REACTION_SANDBOX,CLM_CN,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_CN,REACTION')
            case('DOWNSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%downstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_CN,REACTION')
            case('RATE_CONSTANT')
              call InputReadDouble(input,option,new_reaction%rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_CN,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLM_CN RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                new_reaction%rate_constant = new_reaction%rate_constant * &
                  UnitsConvertToInternal(word,option)
              endif
            case('RESPIRATION_FRACTION')
              call InputReadDouble(input,option,new_reaction%respiration_fraction)
              call InputErrorMsg(input,option,'respiration fraction', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_CN,REACTION')
            case('INHIBITION')
              call InputReadWord(input,option,new_reaction%inhibitor_name, &
                                 PETSC_TRUE)
              call InputErrorMsg(input,option,'inhibitor name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_CN,REACTION')
              call InputReadDouble(input,option,new_reaction%inhibition_constant)
              call InputErrorMsg(input,option,'inhibition constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_CN,REACTION')
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
              call printErrMsg(option)
          end select
        enddo
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine CLM_CN_Read

! ************************************************************************** !
!
! CLM_CN_ReadSkipBlock: Intelligently skips over block
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
subroutine CLM_CN_ReadSkipBlock(this,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word

  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CN Skip Block')
    call StringToUpper(word)   

    select case(trim(word))
      case('POOLS','REACTION')
        call InputSkipToEnd(input,option,word)
    end select
  enddo
  
end subroutine CLM_CN_ReadSkipBlock

! ************************************************************************** !
!
! CLM_CN_Map: Maps coefficients to primary dependent variables
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
subroutine CLM_CN_Map(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type
  use Option_module
  use String_module
  use Immobile_Aux_module
  
  implicit none

  class(reaction_sandbox_clm_cn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH), allocatable :: pool_names(:)
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt :: icount

  type(pool_type), pointer :: cur_pool
  type(clm_cn_reaction_type), pointer :: cur_rxn
  
  ! count # pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    cur_pool => cur_pool%next
  enddo
  this%npool = icount
  
  ! count # reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    cur_rxn => cur_rxn%next
  enddo
  this%nrxn = icount
  
  ! allocate and initialize arrays
  allocate(this%CN_ratio(this%npool))
  allocate(this%pool_id_to_species_id(0,2:this%npool))
  allocate(this%upstream_pool_id(this%nrxn))
  allocate(this%downstream_pool_id(this%nrxn))
  allocate(this%inhibitor_id(this%nrxn))
  allocate(this%rate_constant(this%nrxn))
  allocate(this%respiration_fraction(this%nrxn))
  allocate(this%inhibition_constant(this%nrxn))
  this%CN_ratio = 0.d0
  this%pool_id_to_species_id = 0
  this%upstream_pool_id = 0
  this%downstream_pool_id = 0
  this%inhibitor_id = 0
  this%rate_constant = 0.d0
  this%respiration_fraction = 0.d0
  this%inhibition_constant = 0.d0
  
  ! temporary array for mapping pools in reactions
  allocate(pool_names(this%npool))
  
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    this%CN_ratio(icount) = cur_pool%CN_ratio
    pool_names(icount) = cur_pool%name
    if (cur_pool%CN_ratio < 0.d0) then
      ! Since no CN ratio provided, must provide two species with the
      ! same name as the pool with C or N appended.
      word = trim(cur_pool%name) // 'C'
      this%pool_id_to_species_id(1,icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      word = trim(cur_pool%name) // 'N'
      this%pool_id_to_species_id(2,icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      this%pool_id_to_species_id(0,icount) = 2
      if (minval(this%pool_id_to_species_id) <= 0) then
        option%io_buffer = 'For CLM_CN pools with no CN ratio defined, ' // &
          'the user must define two immobile species with the same name ' // &
          'as the pool with "C" or "N" appended, respectively.'
        call printErrMsg(option)
      endif
    else ! only one species (e.g. SOMX)
      this%pool_id_to_species_id(1,icount) = &
        GetImmobileSpeciesIDFromName(cur_pool%name,reaction%immobile, &
                                     PETSC_TRUE,option)
      this%pool_id_to_species_id(0,icount) = 1
    endif
    cur_pool => cur_pool%next
  enddo
  
  ! map C and N species (solid phase for now)
  word = 'C'
  this%C_species_id = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                   PETSC_TRUE,option)
  word = 'N'
  this%N_species_id = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                   PETSC_TRUE,option)
  
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    this%upstream_pool_id(icount) = &
      StringFindEntryInList(cur_rxn%upstream_pool_name,pool_names)
    if (len_trim(cur_rxn%downstream_pool_name) > 0) then
      this%downstream_pool_id(icount) = &
        StringFindEntryInList(cur_rxn%downstream_pool_name,pool_names)
    endif
    if (len_trim(cur_rxn%inhibitor_name) > 0) then
      this%inhibitor_id(icount) = &
        GetImmobileSpeciesIDFromName(cur_rxn%inhibitor_name,reaction%immobile, &
                                     PETSC_TRUE,option)
    endif
    this%rate_constant(icount) = cur_rxn%rate_constant
    this%respiration_fraction(icount) = cur_rxn%respiration_fraction
    this%inhibition_constant(icount) = cur_rxn%inhibition_constant
    cur_rxn => cur_rxn%next
  enddo 
  
  deallocate(pool_names)
  
end subroutine CLM_CN_Map

! ************************************************************************** !
!
! CLM_CN_React: Evaluates reaction storing residual and/or Jacobian
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
subroutine CLM_CN_React(this,Res,Jac,compute_derivative,rt_auxvar, &
                        global_auxvar,porosity,volume,reaction,option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscInt, parameter :: iphase = 1
  PetscInt :: ipool_up, ipool_down
  PetscInt :: ispec_up, ispec_down
  PetscInt :: ispecC_up, ispecN_up
  PetscInt :: ires_up, ires_dn, ires_c, ires_n
  PetscInt :: iresC_up, iresN_up
  PetscReal :: drate, rate_const, rate
  PetscInt :: irxn
  
  ! inhibition variables
  PetscReal :: F_t
  PetscReal :: F_theta
  PetscReal :: constant_inhibition
  PetscReal :: temp_K
  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0

  PetscReal :: CN_ratio_up
  PetscReal :: resp_frac
  PetscReal :: stoich_N
  PetscReal :: stoich_C
  PetscReal :: stoich_downstream_pool
  PetscReal :: stoich_upstream_pool
  PetscReal :: stoich_upstreamC_pool, stoich_upstreamN_pool
  
  PetscReal :: inhibition_conc, species_inhibition, d_species_inhibition
  PetscInt :: inhibitor_id, ires_inhibitor
  PetscReal :: drate_inhibited
  PetscReal :: temp_real, u, d
  
  ! inhibition due to temperature
  ! Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))
  temp_K = global_auxvar%temp(1) + 273.15d0
  F_t = exp(308.56d0*(one_over_71_02 - 1.d0/(temp_K - 227.13d0)))
  
  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = log(theta_min/global_auxvar%sat(1)) * one_over_log_theta_min 
  
  constant_inhibition = F_t * F_theta
  
  ! indices for C and N species
  ires_c = reaction%offset_immobile + this%C_species_id
  ires_n = reaction%offset_immobile + this%N_species_id

    ! Litter pools
  do irxn = 1, this%nrxn
  
    ipool_up = this%upstream_pool_id(irxn)
    ispec_up = 0 ! this serves as a flag regardless of whether it is set later
    
    rate_const = this%rate_constant(irxn)*constant_inhibition
    resp_frac = this%respiration_fraction(irxn)
    
    ! inhibition by limitting species
    if (this%inhibitor_id(irxn) > 0) then
      inhibitor_id = this%inhibitor_id(irxn) 
      ires_inhibitor = reaction%offset_immobile + inhibitor_id
      inhibition_conc = rt_auxvar%immobile(inhibitor_id)
      temp_real = inhibition_conc + this%inhibition_constant(irxn)
      species_inhibition = inhibition_conc / temp_real
      d_species_inhibition = this%inhibition_constant(irxn) / &
                             (temp_real * temp_real)
    else 
      inhibitor_id = 0
      species_inhibition = 1.d0
    endif
    
    ! contributions for pools of carbon/nitrogen
    if (this%pool_id_to_species_id(0,ipool_up) == 2) then
      ! upstream pool is Litter pool with two species (C,N)
      ispecC_up = this%pool_id_to_species_id(1,ipool_up)
      ispecN_up = this%pool_id_to_species_id(2,ipool_up)
      CN_ratio_up = rt_auxvar%immobile(ispecC_up) / &
                    rt_auxvar%immobile(ispecN_up)
      ! right now, rate calculated based on carbon
      rate = rate_const * rt_auxvar%immobile(ispecC_up) * species_inhibition
      stoich_upstreamC_pool = -1.d0
      stoich_upstreamN_pool = stoich_upstreamC_pool / CN_ratio_up
    else
      ! upstream pool is an SOM pool with one species
      ispec_up = this%pool_id_to_species_id(1,ipool_up)
      CN_ratio_up = this%CN_ratio(ipool_up)
      rate = rate_const * rt_auxvar%immobile(ispec_up) * species_inhibition
      stoich_upstream_pool = -1.d0
    endif

    ipool_down = this%upstream_pool_id(irxn)
    if (ipool_down > 0) then
      ! downstream pool is always SOM
      ispec_down = this%pool_id_to_species_id(1,ipool_down)
      d = twelve_over_14/this%CN_ratio(ipool_down)
    else
      ispec_down = 0
      d = 0.d0
    endif

    ! calculation of n (stoichiometry for N in reaction)
    ! n = u - (1-f) * d
    u = twelve_over_14/CN_ratio_up
    stoich_N = u - (1.d0 - resp_frac) * d
    stoich_C = resp_frac

    ! calculation of residual
    
    ! carbon
    Res(ires_c) = Res(ires_c) - stoich_C * rate
    ! nitrogen
    Res(ires_n) = Res(ires_n) - stoich_N * rate

    if (ispec_up == 0) then ! C and N separate
      ! C species in upstream pool
      iresC_up = reaction%offset_immobile + ispecC_up
      Res(iresC_up) = Res(iresC_up) - stoich_upstreamC_pool * rate
      ! N species in upstream pool
      iresN_up = reaction%offset_immobile + ispecN_up
      Res(iresN_up) = Res(iresN_up) - stoich_upstreamN_pool * rate
    else
      ! upstream pool
      ires_up = reaction%offset_immobile + ispec_up
      Res(ires_up) = Res(ires_up) - stoich_upstream_pool * rate
    endif
    
    if (ispec_down > 0) then
      ! downstream pool
      stoich_downstream_pool = 1.d0 - resp_frac
      ires_dn = reaction%offset_immobile + ispec_down
      Res(ires_dn) = Res(ires_dn) - stoich_downstream_pool * rate
    endif
    
    if (compute_derivative) then
      drate = rate_const * species_inhibition
      drate_inhibited = rate / species_inhibition * d_species_inhibition
      ! upstream pool
      if (ispec_up == 0) then ! C and N separate
        Jac(iresC_up,iresC_up) = Jac(iresC_up,iresC_up) + &
          stoich_upstreamC_pool * drate
        Jac(iresN_up,iresN_up) = Jac(iresN_up,iresN_up) + &
          stoich_upstreamN_pool * drate
        if (inhibitor_id > 0) then
          Jac(iresC_up,ires_inhibitor) = Jac(iresC_up,ires_inhibitor) + &
            stoich_upstreamC_pool * drate_inhibited
          Jac(iresN_up,ires_inhibitor) = Jac(iresN_up,ires_inhibitor) + &
            stoich_upstreamN_pool * drate_inhibited
        endif
        ! set index of Jacobian column for downstream pool, C and N.
        ires_up = iresC_up
      else
        Jac(ires_up,ires_up) = Jac(ires_up,ires_up) + &
          stoich_upstream_pool * drate
        if (inhibitor_id > 0) then
          Jac(ires_up,ires_inhibitor) = Jac(ires_up,ires_inhibitor) + &
            stoich_upstream_pool * drate_inhibited
        endif
      endif
      ! downstream pool
      if (ispec_down > 0) then
        Jac(ires_dn,ires_up) = Jac(ires_dn,ires_up) + &
          stoich_downstream_pool * drate
        if (inhibitor_id > 0) then
          Jac(ires_dn,ires_inhibitor) = Jac(ires_dn,ires_inhibitor) + &
            stoich_downstream_pool * drate_inhibited
        endif
      endif
      ! carbon
      Jac(ires_c,ires_up) = Jac(ires_c,ires_up) + & stoich_C * drate
      ! nitrogen
      Jac(ires_n,ires_up) = Jac(ires_n,ires_up) + stoich_N * drate
      if (inhibitor_id > 0) then
        Jac(ires_c,ires_inhibitor) = Jac(ires_c,ires_inhibitor) + & 
          stoich_C * drate_inhibited
        Jac(ires_n,ires_inhibitor) = Jac(ires_n,ires_inhibitor) + &
          stoich_N * drate_inhibited
      endif
    endif
  enddo
  
end subroutine CLM_CN_React

! ************************************************************************** !
!
! CLM_CN_Destroy: Destroys allocatable or pointer objects created in this 
!                 module
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
subroutine CLM_CN_Destroy(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(clm_cn_reaction_type), pointer :: cur_reaction, prev_reaction
  
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    prev_pool => cur_pool
    cur_pool => cur_pool%next
    deallocate(prev_pool)
    nullify(prev_pool)
  enddo
  
  cur_reaction => this%reactions
  do
    if (.not.associated(cur_reaction)) exit
    prev_reaction => cur_reaction
    cur_reaction => cur_reaction%next
    deallocate(prev_reaction)
    nullify(prev_reaction)
  enddo
  
  call DeallocateArray(this%CN_ratio)
  call DeallocateArray(this%rate_constant)
  call DeallocateArray(this%respiration_fraction)
  call DeallocateArray(this%upstream_pool_id)
  call DeallocateArray(this%downstream_pool_id)
  call DeallocateArray(this%pool_id_to_species_id)
  
end subroutine CLM_CN_Destroy

end module Reaction_Sandbox_CLM_CN_class
