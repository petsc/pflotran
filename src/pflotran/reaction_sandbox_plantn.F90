module Reaction_Sandbox_PlantN_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
! ------------------------------------------------------------------------------
! Description
! to handle plant N uptake with
! 1) Monod type downregulation N/(hs + N)
! 2) Cut off downregulation 1 if N >= N1, 0 if N <= N0, 1 - [1 - (x/d)^2]^2
!     with x = (N - N0)/(N1 - N0)
! 3) inhibition of NH3 on NO3- uptake (assuming plant take NH3 preferentially)
! Author: Guoping Tang
! Date:   07/08/14 
! -----------------------------------------------------------------------------

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_plantn_type
    PetscReal :: rate
    PetscReal :: half_saturation_nh3
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh3_no3
    PetscReal :: residual_conc_nh3
    PetscReal :: residual_conc_no3
    PetscReal :: downreg_no3_0 
    PetscReal :: downreg_no3_1
    PetscReal :: downreg_nh3_0
    PetscReal :: downreg_nh3_1

    PetscInt :: ispec_nh3
    PetscInt :: ispec_no3
    PetscInt :: ispec_plantn
    PetscInt :: ispec_nh4in
    PetscInt :: ispec_no3in
    PetscInt :: ispec_plantndemand

  contains
    procedure, public :: ReadInput => PlantNRead
    procedure, public :: Setup => PlantNSetup
    procedure, public :: Evaluate => PlantNReact
    procedure, public :: Destroy => PlantNDestroy
  end type reaction_sandbox_plantn_type

  public :: PlantNCreate

contains

! **************************************************************************** !
!
! PlantNCreate: Allocates plantn reaction sandbox object.
!
! **************************************************************************** !
function PlantNCreate()

  implicit none
  
  class(reaction_sandbox_plantn_type), pointer :: PlantNCreate

  allocate(PlantNCreate)
  PlantNCreate%rate = 0.d0
  PlantNCreate%half_saturation_nh3 = 1.d-15
  PlantNCreate%half_saturation_no3 = 1.d-15
  PlantNCreate%inhibition_nh3_no3  = 1.d-15
  PlantNCreate%residual_conc_nh3  = 1.d-20
  PlantNCreate%residual_conc_no3  = 1.d-20
  PlantNCreate%downreg_no3_0 = -1.0d-9 
  PlantNCreate%downreg_no3_1 = 1.0d-7
  PlantNCreate%downreg_nh3_0 = -1.0d-9 
  PlantNCreate%downreg_nh3_1 = 1.0d-7
  PlantNCreate%ispec_nh3 = -1
  PlantNCreate%ispec_no3 = -1
  PlantNCreate%ispec_plantn = -1
  PlantNCreate%ispec_nh4in = -1
  PlantNCreate%ispec_no3in = -1
  PlantNCreate%ispec_plantndemand = -1

  nullify(PlantNCreate%next)  
      
end function PlantNCreate

! **************************************************************************** !
!
! PlantNRead: Reads input deck for plantn reaction parameters
!
! **************************************************************************** !
subroutine PlantNRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_plantn_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,PLANTN')
    call StringToUpper(word)   

    select case(trim(word))
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate)
        call InputErrorMsg(input,option,'rate constant', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('HALF_SATURATION_NH3')
        call InputReadDouble(input,option,this%half_saturation_nh3)
        call InputErrorMsg(input,option,'half saturation NH3-', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('HALF_SATURATION_NO3')
        call InputReadDouble(input,option,this%half_saturation_no3)
        call InputErrorMsg(input,option,'half saturation NO3-', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('NH3_INHIBITION_NO3')
        call InputReadDouble(input,option,this%inhibition_nh3_no3)
        call InputErrorMsg(input,option,'NH3 inhibition on NO3-', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('RESIDUAL_CONC_NH3')
        call InputReadDouble(input,option,this%residual_conc_nh3)
        call InputErrorMsg(input,option,'residual concentration nh3', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('RESIDUAL_CONC_NO3')
        call InputReadDouble(input,option,this%residual_conc_no3)
        call InputErrorMsg(input,option,'residual concentration no3-', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('DOWNREGULATE_NH3')
        call InputReadDouble(input,option,this%downreg_nh3_0)
        call InputErrorMsg(input,option,'downreg_nh3_0', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
        call InputReadDouble(input,option,this%downreg_nh3_1)
        call InputErrorMsg(input,option,'downreg_nh3_1', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
        if (this%downreg_nh3_0 > this%downreg_nh3_1) then
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN,' // &
            'NH4+ down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('DOWNREGULATE_NO3')
        call InputReadDouble(input,option,this%downreg_no3_0)
        call InputErrorMsg(input,option,'downreg_no3_0', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
        call InputReadDouble(input,option,this%downreg_no3_1)
        call InputErrorMsg(input,option,'downreg_no3_1', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
        if (this%downreg_no3_0 > this%downreg_no3_1) then
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN,' // &
            'NO3- down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine PlantNRead

! **************************************************************************** !
!
! PlantNSetup: Sets up the plantn reaction with parameters
!
! **************************************************************************** !
subroutine PlantNSetup(this,reaction,option)

  use Reaction_Aux_module
  use Option_module
  use Immobile_Aux_module

  implicit none
  
  class(reaction_sandbox_plantn_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
  word = 'NH4+'
  this%ispec_nh3 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  if (this%ispec_nh3 < 0) then
    word = 'NH3(aq)'
    this%ispec_nh3 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE, &
      option)
  endif
  
  if (this%ispec_nh3 < 0) then
    option%io_buffer = 'Neither NH4+ nor NH3(aq) is specified in the input' // &
      'file for PlantN sandbox!'
    call printErrMsg(option)
  endif

  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  word = 'PlantN'
  this%ispec_plantn = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
    PETSC_FALSE,option)

  if (this%ispec_plantn < 0) then
    option%io_buffer = 'PlantN is specified in the input file!'
    call printErrMsg(option)
  endif

  word = 'Ain'
  this%ispec_nh4in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
    PETSC_FALSE,option)

  word = 'Tin'
  this%ispec_no3in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
    PETSC_FALSE,option)

  word = 'Plantndemand'
  this%ispec_plantndemand = GetImmobileSpeciesIDFromName(word, &
    reaction%immobile, PETSC_FALSE,option)

end subroutine PlantNSetup

! **************************************************************************** !
!
! PlantNReact: Evaluates reaction storing residual and/or Jacobian
!
! **************************************************************************** !
subroutine PlantNReact(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                       global_auxvar,material_auxvar,reaction,option)

  use Option_module
  use Reaction_Aux_module
  use Immobile_Aux_module
  use Material_Aux_class, only : material_auxvar_type

  implicit none

  class(reaction_sandbox_plantn_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscBool :: compute_derivative

  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: rate, drate, concN, rate0
  PetscReal :: volume, porosity
  PetscInt :: local_id
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  PetscInt, parameter :: iphase = 1
  PetscInt :: ires_nh3, ires_no3, ires_plantn
  PetscInt :: ires_nh4in, ires_no3in
  PetscInt :: ires_plantndemand

  PetscReal :: c_nh3         ! concentration (mole/L)
  PetscReal :: ac_nh3        ! activity coefficient
  PetscReal :: f_nh3         ! nh3 / (half_saturation + nh3)
  PetscReal :: d_nh3         ! half_saturation / (half_saturation + nh3)^2
  PetscReal :: f_nh3_inhibit ! inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: d_nh3_inhibit ! d inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: ac_no3        ! activity coefficient
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: temp_real

  PetscReal :: rate_nh3
  PetscReal :: rate_no3
  PetscReal :: drate_nh3_dnh3
  PetscReal :: drate_no3_dno3
  PetscReal :: drate_no3_dnh3
  PetscReal :: c_plantn, c_plantno3, c_plantnh3, c_plantndemand
  PetscReal :: xxx, delta, regulator, dregulator

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  ires_nh4in = this%ispec_nh4in + reaction%offset_immobile
  ires_no3in = this%ispec_no3in + reaction%offset_immobile
  ires_plantndemand = this%ispec_plantndemand + reaction%offset_immobile

  if (this%ispec_plantn < 0) then
    option%io_buffer = 'PlantN is not specified in the input file!'
    call printErrMsg(option)
  else
    ires_plantn = this%ispec_plantn + reaction%offset_immobile
  endif

  if (this%ispec_nh3 > 0) then
    c_nh3     = rt_auxvar%pri_molal(this%ispec_nh3)
    ac_nh3    = rt_auxvar%pri_act_coef(this%ispec_nh3)
    temp_real = (c_nh3 - this%residual_conc_nh3) * ac_nh3 &
              + this%half_saturation_nh3
    f_nh3     = (c_nh3 - this%residual_conc_nh3) * ac_nh3 / temp_real 
    d_nh3     = ac_nh3 * this%half_saturation_nh3 / temp_real / temp_real

    if (this%downreg_nh3_0 > 0.0d0) then
      ! additional down regulation for plant NH4+ uptake
      if (c_nh3 <= this%downreg_nh3_0) then
        regulator = 0.0d0
        dregulator = 0.0d0
      elseif (c_nh3 >= this%downreg_nh3_1 .or. &
              this%downreg_nh3_1 - this%downreg_nh3_0 <= 1.0d-20) then
        regulator = 1.0d0
        dregulator = 0.0d0
      else
        xxx = c_nh3 - this%downreg_nh3_0
        delta = this%downreg_nh3_1 - this%downreg_nh3_0
        regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                   / delta / delta
      endif
    
      ! rate = rate_orginal * regulator
      ! drate = drate_original * regulator + rate_orginal * dregulator
      d_nh3 = d_nh3 * regulator + f_nh3 * dregulator

      f_nh3 = f_nh3 * regulator

    endif

    ires_nh3 = this%ispec_nh3
  endif

  if (this%ispec_no3 > 0) then
    c_no3     = rt_auxvar%pri_molal(this%ispec_no3)
    ac_no3    = rt_auxvar%pri_act_coef(this%ispec_no3)
    temp_real = c_no3 -this%residual_conc_no3 + this%half_saturation_no3
    f_no3 = (c_no3 - this%residual_conc_no3) * ac_no3 / temp_real
    d_no3 = ac_no3 * this%half_saturation_no3 / temp_real / temp_real

    if (this%downreg_no3_0 > 0.0d0) then
      ! additional down regulation for plant NO3- uptake
      if (c_no3 <= this%downreg_no3_0) then
        regulator = 0.0d0
        dregulator = 0.0d0
      elseif (c_no3 >= this%downreg_no3_1 .or. &
              this%downreg_no3_1 - this%downreg_no3_0 <= 1.0d-20) then
        regulator = 1.0d0
        dregulator = 0.0d0
      else
        xxx = c_no3 - this%downreg_no3_0
        delta = this%downreg_no3_1 - this%downreg_no3_0
        regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                   / delta / delta
      endif

      ! rate = rate_orginal * regulator
      ! drate = drate_original * regulator + rate_orginal * dregulator
      d_no3 = d_no3 * regulator + f_no3 * dregulator
      f_no3 = f_no3 * regulator

    endif

    ires_no3 = this%ispec_no3

    if (this%ispec_nh3 > 0) then
      temp_real = this%inhibition_nh3_no3 + c_nh3 * ac_nh3
      f_nh3_inhibit = this%inhibition_nh3_no3/temp_real
      d_nh3_inhibit = -1.0d0 * this%inhibition_nh3_no3 * ac_nh3 &
                    / temp_real / temp_real
    else
      f_nh3_inhibit = 1.0d0
      d_nh3_inhibit = 0.0d0
    endif
  endif

  rate_nh3 = this%rate

  if (this%ispec_plantndemand > 0) then
    Residual(ires_plantndemand) = Residual(ires_plantndemand) - rate_nh3
  endif

  if (this%ispec_nh3 > 0) then

    drate_nh3_dnh3 = rate_nh3 * d_nh3
    rate_nh3 = rate_nh3 * f_nh3

    Residual(ires_nh3) = Residual(ires_nh3) + rate_nh3
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nh3

    if (this%ispec_nh4in > 0) then
      Residual(ires_nh4in) = Residual(ires_nh4in) - rate_nh3
    endif

    if (compute_derivative) then
      Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_nh3_dnh3

      Jacobian(ires_plantn,ires_nh3) = Jacobian(ires_plantn,ires_nh3) &
                                     - drate_nh3_dnh3

      if (this%ispec_nh4in > 0) then
        Jacobian(ires_nh4in,ires_nh3) = Jacobian(ires_nh4in,ires_nh3) &
                                      - drate_nh3_dnh3
      endif
    endif
  endif

  if (this%ispec_no3 > 0) then
    rate_no3 = this%rate * f_nh3_inhibit * f_no3
    drate_no3_dno3 = this%rate * f_nh3_inhibit * d_no3
    drate_no3_dnh3 = this%rate * d_nh3_inhibit * f_no3

    Residual(ires_no3) = Residual(ires_no3) + rate_no3
    Residual(ires_plantn) = Residual(ires_plantn) - rate_no3

    if (this%ispec_no3in > 0) then
      Residual(ires_no3in) = Residual(ires_no3in) - rate_no3
    endif

    if (compute_derivative) then
      Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + drate_no3_dno3

      Jacobian(ires_plantn,ires_no3) = Jacobian(ires_plantn,ires_no3) &
                                     - drate_no3_dno3

      if (this%ispec_no3in > 0) then
        Jacobian(ires_no3in,ires_no3) = Jacobian(ires_no3in,ires_no3) &
                                      - drate_no3_dno3
      endif

      Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) + drate_no3_dnh3

      Jacobian(ires_plantn,ires_nh3) = Jacobian(ires_plantn,ires_nh3) &
                                     - drate_no3_dnh3

      if (this%ispec_no3in > 0) then
        Jacobian(ires_no3in,ires_nh3) = Jacobian(ires_no3in,ires_nh3) &
                                      - drate_no3_dnh3
      endif

      if (PETSC_FALSE) then
        c_plantn = rt_auxvar%immobile(this%ispec_plantn)
        write(*, *) c_nh3, c_no3, this%rate, rate_nh3, rate_no3, c_plantn
      endif
    endif
  endif

end subroutine PlantNReact

! **************************************************************************** !
!
! PlantNDestroy: Destroys allocatable or pointer objects created in this module
!
! **************************************************************************** !
subroutine PlantNDestroy(this)

  implicit none
  
  class(reaction_sandbox_plantn_type) :: this  

end subroutine PlantNDestroy

end module Reaction_Sandbox_PlantN_class
