module Reaction_Sandbox_CLM_CN_class
#include "finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clm_cn_type
    PetscInt :: nlitter
    character(len=MAXWORDLENGTH), pointer :: litter_names(:)
    PetscReal, pointer :: litter_coefs(:,:)
    PetscReal, pointer :: litter_rate_const(:)
    PetscInt, pointer :: litter_species_ids(:,:)
    PetscInt, pointer :: litter_upstream_pool(:)
    PetscInt :: nsom
    character(len=MAXWORDLENGTH), pointer :: som_names(:)
    PetscReal, pointer :: som_coefs(:,:)
    PetscReal, pointer :: som_rate_const(:)
    PetscInt, pointer :: som_species_ids(:,:)
    PetscInt, pointer :: som_upstream_pool(:)
    character(len=MAXSTRINGLENGTH), pointer :: reaction_strings(:)
  contains
    procedure, public :: Init => CLM_CN_Init
    procedure, public :: ReadInput => CLM_CN_Read
    procedure, public :: Evaluate => CLM_CN_React
    procedure, public :: Destroy => CLM_CN_Destroy
  end type reaction_sandbox_clm_cn_type
  
  public :: CLM_CN_Create

contains

! ************************************************************************** !

function CLM_CN_Create()
  ! 
  ! RSandboxInit: Initializes reaction sandbox at beginning of simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/08/12
  ! 

  implicit none
  
  type(reaction_sandbox_clm_cn_type), pointer :: CLM_CN_Create
  
  allocate(CLM_CN_Create)
  CLM_CN_Create%nlitter = 0
  CLM_CN_Create%nsom = 0
  nullify(CLM_CN_Create%litter_names)
  nullify(CLM_CN_Create%som_names)
  nullify(CLM_CN_Create%litter_coefs)
  nullify(CLM_CN_Create%som_coefs)
  nullify(CLM_CN_Create%litter_species_ids)
  nullify(CLM_CN_Create%litter_upstream_pool)
  nullify(CLM_CN_Create%som_species_ids)
  nullify(CLM_CN_Create%litter_rate_const)
  nullify(CLM_CN_Create%som_rate_const)
  nullify(CLM_CN_Create%som_upstream_pool)
  nullify(CLM_CN_Create%next)
  
end function CLM_CN_Create

! ************************************************************************** !

subroutine CLM_CN_Init(this,reaction,option)
  ! 
  ! Initializes reaction sandbox at beginning of simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/13
  ! 

  use Reaction_Aux_module
  use Option_module
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction  
  
  call CLM_CN_Map(this,reaction,option)

end subroutine CLM_CN_Init

! ************************************************************************** !

subroutine CLM_CN_Map(this,reaction,option)
  ! 
  ! Maps coefficients to primary dependent variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/13
  ! 

  use Reaction_Aux_module
  use Option_module
  use String_module
  use Input_Aux_module
  use Database_Aux_module
  use Immobile_Aux_module
  
  implicit none

  class(reaction_sandbox_clm_cn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  type(database_rxn_type), pointer ::dbase_rxn
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  character(len=MAXSTRINGLENGTH) :: temp_string
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt :: ilit, isom, i, lit_max_spec, som_max_spec, pool_id
  PetscErrorCode :: ierr

  ! determine how many SOM and Litters
  ilit = 0
  isom = 0
  lit_max_spec = 0
  som_max_spec = 0
  do i = 1, size(this%reaction_strings)
    strings => StringSplit(this%reaction_strings(i),';')
    dbase_rxn => &
      DatabaseRxnCreateFromRxnString(strings(1), &
                                     reaction%naqcomp, &
                                     reaction%offset_aqueous, &
                                     reaction%primary_species_names, &
                                     reaction%nimcomp, &
                                     reaction%offset_immobile, &
                                     reaction%immobile%names, &
                                     option)
    if (StringStartsWith(this%reaction_strings(i),'SOM')) then
      som_max_spec = max(som_max_spec,dbase_rxn%nspec)
      isom = isom + 1
    else if (StringStartsWith(this%reaction_strings(i),'Lit')) then
      lit_max_spec = max(lit_max_spec,dbase_rxn%nspec)
      ilit = ilit + 1
    else
      option%io_buffer = 'Unrecognized reaction string for CLM_CN.'
      call printErrMsg(option)
    endif
    call DatabaseRxnDestroy(dbase_rxn)    
  enddo
  this%nsom = isom
  this%nlitter = ilit
  
  allocate(this%litter_coefs(lit_max_spec,this%nlitter))
  allocate(this%litter_species_ids(0:lit_max_spec,this%nlitter))
  allocate(this%litter_rate_const(this%nlitter))
  allocate(this%litter_upstream_pool(this%nlitter))
  this%litter_coefs = 0.d0
  this%litter_species_ids = 0
  this%litter_rate_const = 0.d0

  allocate(this%som_coefs(som_max_spec,this%nsom))
  allocate(this%som_species_ids(0:som_max_spec,this%nsom))
  allocate(this%som_rate_const(this%nsom))
  allocate(this%som_upstream_pool(this%nsom))
  this%som_coefs = 0.d0
  this%som_species_ids = 0
  this%som_rate_const = 0.d0
  
  ilit = 0
  isom = 0
  do i = 1, size(this%reaction_strings)
    strings => StringSplit(this%reaction_strings(i),';')
    ! find upstream pool.  The first species in the reaction string
    temp_string = strings(1)
    ierr = 0
    do ! need to loop inorder to remove first stoichiometry if it exists
      call InputReadWord(temp_string,word,PETSC_TRUE,ierr)
      ! no need for an error message since the string would have to be
      ! legitimate to pass the database rxn creation above.
      ! PETSC_FALSE instucts not to throw an error, but to return
      ! the unset value (= -999)
      pool_id = GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                             PETSC_FALSE,option)
      if (ierr /= 0 .or. pool_id > 0) exit
    enddo
    ! process full reaction
    dbase_rxn => &
      DatabaseRxnCreateFromRxnString(strings(1), &
                                     reaction%naqcomp, &
                                     reaction%offset_aqueous, &
                                     reaction%primary_species_names, &
                                     reaction%nimcomp, &
                                     reaction%offset_immobile, &
                                     reaction%immobile%names, &
                                     option)
    if (StringStartsWith(strings(1),'Lit')) then
      ilit = ilit + 1
      this%litter_coefs(1:dbase_rxn%nspec,ilit) = dbase_rxn%stoich(:)
      this%litter_species_ids(0,ilit) = dbase_rxn%nspec
      this%litter_species_ids(1:dbase_rxn%nspec,ilit) = &
        dbase_rxn%spec_ids(:) - reaction%offset_immobile
      this%litter_upstream_pool(ilit) = pool_id
      call InputReadDouble(strings(2),option,this%litter_rate_const(ilit),ierr)
    else if (StringStartsWith(strings(1),'SOM')) then
      isom = isom + 1
      this%som_coefs(1:dbase_rxn%nspec,isom) = dbase_rxn%stoich(:)
      this%som_species_ids(0,isom) = dbase_rxn%nspec
      this%som_species_ids(1:dbase_rxn%nspec,isom) = &
        dbase_rxn%spec_ids(:) - reaction%offset_immobile
      this%som_upstream_pool(isom) = pool_id
      call InputReadDouble(strings(2),option,this%som_rate_const(isom),ierr)
    endif
    ! as long as InputReadDouble is the last call in the conditional above
    ! both branches can use thir error statement.
    if (ierr /= 0) then
      option%io_buffer = 'Rate not found in CLM_CN reacton string (' // &
        trim(adjustl(this%reaction_strings(i))) // ').'
      call printErrMsg(option)
    endif
    deallocate(strings)
    nullify(strings)
    call DatabaseRxnDestroy(dbase_rxn)
  enddo
  
end subroutine CLM_CN_Map

! ************************************************************************** !

subroutine CLM_CN_Read(this,input,option)
  ! 
  ! Reads input deck for reaction sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/13
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word

  PetscInt :: icoef, icount
  character(len=MAXWORDLENGTH) :: names(TWELVE_INTEGER)
  PetscReal :: litter_coefs(THREE_INTEGER,TWELVE_INTEGER)
  PetscReal :: som_coefs(THREE_INTEGER,TWELVE_INTEGER)

  character(len=MAXSTRINGLENGTH) :: strings(100)

  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REACTION_SANDBOX')
    call StringToUpper(word)   

    select case(trim(word))
      case('REACTIONS')
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   
          icount = icount + 1
          strings(icount) = adjustl(input%buf)
        enddo
#if 0    
      case('LITTER_NAMES')
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   
          icount = icount + 1
          call InputReadWord(input,option,names(icount),PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,CLM_CN,LITTERNAMES')
        enddo
        ! deallocate just incase already allocated by mistake
        call DeallocateArray(this%litter_names)
        allocate(this%litter_names(icount))
        this%litter_names = names(1:icount)
      case('SOM_NAMES')
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   
          icount = icount + 1
          call InputReadWord(input,option,som_names(icount),PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,CLM_CN,SOM_NAMES')
        enddo
        call DeallocateArray(this%som_names)
        allocate(this%som_names(icount))
        this%som_names = names(1:icount)
      case('LITTER_COEFFICIENTS')
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   
          icount = icount + 1
          do icoef = 1, 3
            call InputReadWord(input,option,litter_coefs(icoef,icount), &
                               PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword', &
                               'CHEMISTRY,CLM_CN,LITTER_COEFFICIENTS')
          enddo
        enddo
        call DeallocateArray(this%litter_coefs)
        allocate(this%litter_coefs(icount))
        this%litter_coefs = litter_coefs(1:icount)
      case('SOM_COEFFICIENTS')
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   
          icount = icount + 1
          do icoef = 1, 3
            call InputReadWord(input,option,som_coefs(icoef,icount), &
                               PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword', &
                               'CHEMISTRY,CLM_CN,SOM_COEFFICIENTS')
          enddo
        enddo 
        call DeallocateArray(this%som_coefs)
        allocate(this%som_coefs(icount))
        this%som_coefs = som_coefs(1:icount)      
#endif        
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
#if 0  
  ! error checking
  if (size(this%litter_names) /= size(this%litter_coefs,2)) then
    option%io_buffer = &
      'Number of Liter names does not match number of coefficients.'
    call printErrMsg(option)
  else
    this%nlitter = this%litter_names
  endif
  if (size(this%som_names) /= size(this%som_coefs,2)) then
    option%io_buffer = &
      'Number of SOM names does not match number of coefficients.'
    call printErrMsg(option)
  else
    this%nsom = this%som_names
  endif
#endif  

  if (icount > 0) then
    allocate(this%reaction_strings(icount))
    this%reaction_strings = strings(1:icount)
  endif
  
end subroutine CLM_CN_Read

! ************************************************************************** !

subroutine CLM_CN_React(this,Res,Jac,compute_derivative,rt_auxvar, &
                        global_auxvar,porosity,volume,reaction,option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/13
  ! 

  use Option_module
  use Reaction_Aux_module
  
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
  PetscInt :: ilit, isom, iimmobile, idof, i, icomp
  PetscReal :: drate, rate_const, rate
  
  ! inhibition variables
  PetscReal :: F_t
  PetscReal :: F_theta
  PetscReal :: inhibition
  PetscReal :: temp_K
  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1

  
  ! inhibition due to temperature
  ! Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))
  temp_K = global_auxvar%temp + 273.15d0
  F_t = exp(308.56d0*(one_over_71_02 - 1.d0/(temp_K - 227.13d0)))
  
  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = log(theta_min/global_auxvar%sat(1)) * one_over_log_theta_min 
  
  inhibition = F_t * F_theta
  
  ! Litter pools
  do ilit = 1, this%nlitter
    ! first species will be the used in the rate
    iimmobile = this%litter_upstream_pool(ilit)
    idof = reaction%offset_immobile + iimmobile
    rate_const = this%litter_rate_const(ilit)*inhibition
    rate = rate_const * rt_auxvar%immobile(iimmobile)
    do i = 1, this%litter_species_ids(0,ilit)
      icomp = reaction%offset_immobile + this%litter_species_ids(i,ilit) 
      Res(icomp) = Res(icomp) - this%litter_coefs(i,ilit)*rate
    enddo
    if (compute_derivative) then
      drate = rate_const
      do i = 1, this%litter_species_ids(0,ilit)
        icomp = reaction%offset_immobile + this%litter_species_ids(i,ilit) 
        Jac(icomp,idof) = Jac(icomp,idof) + this%litter_coefs(i,ilit)*drate
      enddo
    endif
  enddo
  
  ! SOM pools
  do isom = 1, this%nsom
    ! first species will be the used in the rate
    iimmobile = this%som_upstream_pool(isom)
    idof = reaction%offset_immobile + iimmobile
    rate_const = this%som_rate_const(isom)*inhibition
    rate = rate_const * rt_auxvar%immobile(iimmobile)
    do i = 1, this%som_species_ids(0,isom)
      icomp = reaction%offset_immobile + this%som_species_ids(i,isom) 
      Res(icomp) = Res(icomp) - this%som_coefs(i,isom)*rate
    enddo
    if (compute_derivative) then
      drate = rate_const
      do i = 1, this%som_species_ids(0,isom)
        icomp = reaction%offset_immobile + this%som_species_ids(i,isom) 
        Jac(icomp,idof) = Jac(icomp,idof) + this%som_coefs(i,isom)*drate
      enddo
    endif
  enddo
  
end subroutine CLM_CN_React

! ************************************************************************** !

subroutine CLM_CN_Destroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/13
  ! 

  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  
end subroutine CLM_CN_Destroy

end module Reaction_Sandbox_CLM_CN_class
