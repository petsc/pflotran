module Reaction_Sandbox_Cyber_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"


  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_cyber_type
    PetscInt :: nh4_id
    PetscInt :: o2_id
    PetscInt :: no3_id
    PetscInt :: no2_id
    PetscInt :: n2_id
    PetscInt :: doc_id
    PetscInt :: biomass_id
    PetscInt :: co2_id
    PetscReal :: f1
    PetscReal :: f2
    PetscReal :: f3
    PetscReal :: f_act   ! fraction of active biomass
    PetscReal :: k_deg   ! biomass degradation rate
    PetscReal :: k1      ! nitrate rate constant
    PetscReal :: k2      ! nitrite rate constant
    PetscReal :: k3      ! oxygen rate constant
    PetscReal :: Kd1
    PetscReal :: Ka1
    PetscReal :: Kd2
    PetscReal :: Ka2
    PetscReal :: Kd3
    PetscReal :: Ka3
    PetscReal :: stoich_1_doc
    PetscReal :: stoich_1_nh4
    PetscReal :: stoich_1_no3
    PetscReal :: stoich_1_no2
    PetscReal :: stoich_1_co2
    PetscReal :: stoich_1_biomass
    PetscReal :: stoich_2_doc
    PetscReal :: stoich_2_nh4
    PetscReal :: stoich_2_no2
    PetscReal :: stoich_2_n2
    PetscReal :: stoich_2_co2
    PetscReal :: stoich_2_biomass
    PetscReal :: stoich_3_doc
    PetscReal :: stoich_3_nh4
    PetscReal :: stoich_3_o2
    PetscReal :: stoich_3_co2
    PetscReal :: stoich_3_biomass
    PetscInt :: nrxn
    PetscInt, pointer :: nrow(:)
    PetscInt, pointer :: ncol(:)
    PetscInt, pointer :: irow(:,:)
    PetscInt, pointer :: icol(:,:)
    PetscReal, pointer :: stoich_row(:,:)
  contains
    procedure, public :: ReadInput => CyberRead
    procedure, public :: Setup => CyberSetup
    procedure, public :: Evaluate => CyberReact
    procedure, public :: Destroy => CyberDestroy
  end type reaction_sandbox_cyber_type
  
  public :: CyberCreate

contains

! ************************************************************************** !

function CyberCreate()
  ! 
  ! Allocates PNNL N reaction object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  implicit none
  
  class(reaction_sandbox_cyber_type), pointer :: CyberCreate

  allocate(CyberCreate)
  CyberCreate%o2_id = UNINITIALIZED_INTEGER
  CyberCreate%nh4_id = UNINITIALIZED_INTEGER
  CyberCreate%no3_id = UNINITIALIZED_INTEGER
  CyberCreate%no2_id = UNINITIALIZED_INTEGER
  CyberCreate%n2_id = UNINITIALIZED_INTEGER
  CyberCreate%doc_id = UNINITIALIZED_INTEGER
  CyberCreate%biomass_id = UNINITIALIZED_INTEGER
  CyberCreate%co2_id = UNINITIALIZED_INTEGER
  CyberCreate%f1 = UNINITIALIZED_DOUBLE  
  CyberCreate%f2 = UNINITIALIZED_DOUBLE  
  CyberCreate%f3 = UNINITIALIZED_DOUBLE  
  CyberCreate%f_act = UNINITIALIZED_DOUBLE  
  CyberCreate%k_deg = UNINITIALIZED_DOUBLE  
  CyberCreate%k1 = UNINITIALIZED_DOUBLE  
  CyberCreate%k2 = UNINITIALIZED_DOUBLE  
  CyberCreate%k3 = UNINITIALIZED_DOUBLE  
  CyberCreate%Kd1 = UNINITIALIZED_DOUBLE  
  CyberCreate%Ka1 = UNINITIALIZED_DOUBLE  
  CyberCreate%Kd2 = UNINITIALIZED_DOUBLE  
  CyberCreate%Ka2 = UNINITIALIZED_DOUBLE  
  CyberCreate%Kd3 = UNINITIALIZED_DOUBLE  
  CyberCreate%Ka3 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_1_doc = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_1_nh4 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_1_no3 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_1_no2 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_1_co2 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_1_biomass = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_2_doc = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_2_nh4 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_2_no2 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_2_n2 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_2_co2 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_2_biomass = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_3_doc = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_3_nh4 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_3_o2 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_3_co2 = UNINITIALIZED_DOUBLE  
  CyberCreate%stoich_3_biomass = UNINITIALIZED_DOUBLE  
  CyberCreate%nrxn = UNINITIALIZED_INTEGER
  nullify(CyberCreate%nrow)
  nullify(CyberCreate%ncol)
  nullify(CyberCreate%irow)
  nullify(CyberCreate%icol)
  nullify(CyberCreate%stoich_row)

  nullify(CyberCreate%next)  
      
end function CyberCreate

! ************************************************************************** !

subroutine CyberRead(this,input,option)
  ! 
  ! Reads input deck for PNNL N reaction parameters (if any)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_cyber_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units, units
  PetscReal :: units_conversion
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CYBERNETIC')
    call StringToUpper(word)   

    select case(trim(word))
      case('F1')
        call InputReadDouble(input,option,this%f1)  
        call InputErrorMsg(input,option,'f1', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
      case('F2')
        call InputReadDouble(input,option,this%f2)  
        call InputErrorMsg(input,option,'f2', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
      case('F3')
        call InputReadDouble(input,option,this%f3)  
        call InputErrorMsg(input,option,'f3', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
      case('K1','K_NO3-')
        call InputReadDouble(input,option,this%k1)  
        call InputErrorMsg(input,option,'k1', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = '1/sec'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%k1 = this%k1 * units_conversion
        endif
      case('K2','K_NO2-')
        call InputReadDouble(input,option,this%k2)  
        call InputErrorMsg(input,option,'k2', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = '1/sec'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%k2 = this%k2 * units_conversion
        endif
      case('K3','K_O2(aq)')
        call InputReadDouble(input,option,this%k3)  
        call InputErrorMsg(input,option,'k3', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = '1/sec'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%k3 = this%k3 * units_conversion
        endif    
      case('KA1','KA_NO3-')
        call InputReadDouble(input,option,this%Ka1)  
        call InputErrorMsg(input,option,'Ka1', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = 'M'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%Ka1 = this%Ka1 * units_conversion
        endif    
      case('KA2','KA_NO2-')
        call InputReadDouble(input,option,this%Ka2)  
        call InputErrorMsg(input,option,'Ka2', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = 'M'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%Ka2 = this%Ka2 * units_conversion
        endif    
      case('KA3','KA_O2(aq)')
        call InputReadDouble(input,option,this%Ka3)  
        call InputErrorMsg(input,option,'Ka3', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = 'M'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%Ka3 = this%Ka3 * units_conversion
        endif    
      case('KD1','KD_NO3-')
        call InputReadDouble(input,option,this%Kd1)  
        call InputErrorMsg(input,option,'Kd1', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = 'M'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%Kd1 = this%Kd1 * units_conversion
        endif    
      case('KD2','KD_NO2-')
        call InputReadDouble(input,option,this%Kd2)  
        call InputErrorMsg(input,option,'Kd2', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = 'M'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%Kd2 = this%Kd2 * units_conversion
        endif    
      case('KD3','KD_O2(aq)')
        call InputReadDouble(input,option,this%Kd3)  
        call InputErrorMsg(input,option,'Kd3', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = 'M'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%Kd3 = this%Kd3 * units_conversion
        endif   
      case('KDEG')
        call InputReadDouble(input,option,this%k_deg)  
        call InputErrorMsg(input,option,'kdeg', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = '1/sec'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%k_deg = this%k_deg * units_conversion
        endif    
      case('F_ACT')
        call InputReadDouble(input,option,this%f_act)  
        call InputErrorMsg(input,option,'f_act', &
                           'CHEMISTRY,REACTION_SANDBOX_CYBERNETIC')      
        call InputReadWord(input,option,units,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = '1/sec'
          units_conversion = UnitsConvertToInternal(units,internal_units,option)
          this%f_act = this%f_act * units_conversion
        endif    
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,CYBERNETIC',option)
    end select
  enddo
  
end subroutine CyberRead

! ************************************************************************** !

subroutine CyberSetup(this,reaction,option)
  ! 
  ! Sets up the PNNL N reaction either with parameters either
  ! read from the input deck or hardwired.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_cyber_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: irxn
  
  PetscReal, parameter :: per_day_to_per_sec = 1.d0 / 24.d0 / 3600.d0

  word = 'NH4+'
  this%nh4_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%o2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO3-'
  this%no3_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO2-'
  this%no2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'N2(aq)'
  this%n2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'CH2O(aq)'
  this%doc_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'C5H7O2N(aq)'
  this%biomass_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
!    GetImmobileSpeciesIDFromName(word,reaction%immobile,option) + reaction%offset_immobile
  word = 'CO2(aq)'
  this%co2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  
  ! constants based on Hyun's writeup on 6/21/16 entitled "Mini-cybernetic 
  ! model of batch denitrification process"
  if (Uninitialized(this%f1)) this%f1 = 0.497d0
  if (Uninitialized(this%f2)) this%f2 = 0.999d0
  if (Uninitialized(this%f3)) this%f3 = 0.066d0
  if (Uninitialized(this%k1)) this%k1 = 25.04d0 * per_day_to_per_sec
  if (Uninitialized(this%k2)) this%k2 = 17.82d0 * per_day_to_per_sec
  if (Uninitialized(this%k3)) this%k3 = 75.12d0 * per_day_to_per_sec
  if (Uninitialized(this%Kd1)) this%Kd1 = 0.25d-3
  if (Uninitialized(this%Kd2)) this%Kd2 = 0.25d-3
  if (Uninitialized(this%Kd3)) this%Kd3 = 0.25d-3
  if (Uninitialized(this%Ka1)) this%Ka1 = 0.001d-3
  if (Uninitialized(this%Ka2)) this%Ka2 = 0.004d-3
  if (Uninitialized(this%Ka3)) this%Ka3 = 0.001d-3
  if (Uninitialized(this%f_act)) this%f_act = 0.126d0
  if (Uninitialized(this%k_deg)) this%k_deg = 0.532d0 * per_day_to_per_sec

! uncomment these to zero out reactions
!  this%k1 = 0.d0
!  this%k2 = 0.d0
!  this%k3 = 0.d0
!  this%Kd1 = 0.d0
!  this%Ka1 = 0.d0
!  this%f_act = 1.d40
!  this%k_deg = 0.d0

  this%stoich_1_doc = -1.d0
  this%stoich_1_nh4 = -0.2d0*(1.d0-this%f1)
  this%stoich_1_no3 = -2.d0*this%f1
  this%stoich_1_no2 = 2.d0*this%f1
  this%stoich_1_co2 = this%f1
  this%stoich_1_biomass = 0.2d0*(1.d0-this%f1)

  this%stoich_2_doc = -1.d0
  this%stoich_2_nh4 = -0.2d0*(1.d0-this%f2)
  this%stoich_2_no2 = -4.d0/3.d0*this%f2
  this%stoich_2_n2 = 2.d0/3.d0*this%f2
  this%stoich_2_co2 = this%f2
  this%stoich_2_biomass = 0.2d0*(1.d0-this%f2)

  this%stoich_3_doc = -1.d0
  this%stoich_3_nh4 = -0.2d0*(1.d0-this%f3)
  this%stoich_3_o2 = -1.d0*this%f3
  this%stoich_3_co2 = this%f3
  this%stoich_3_biomass = 0.2d0*(1.d0-this%f3)

  if (this%k3 > 1.d-40) then
  this%nrxn = 3
  else
    this%nrxn = 2
  endif
  allocate(this%nrow(this%nrxn))
  this%nrow = UNINITIALIZED_INTEGER
  allocate(this%ncol(this%nrxn))
  this%ncol = UNINITIALIZED_INTEGER
  allocate(this%irow(6,this%nrxn))
  this%irow = UNINITIALIZED_INTEGER
  allocate(this%icol(4,this%nrxn))
  this%icol = UNINITIALIZED_INTEGER
  allocate(this%stoich_row(6,this%nrxn))
  this%stoich_row = UNINITIALIZED_DOUBLE
  
  ! NO3- -> NO2-
  irxn = 1
  this%nrow(irxn) = 6
  this%irow(1,irxn) = this%doc_id
  this%irow(2,irxn) = this%nh4_id
  this%irow(3,irxn) = this%no3_id
  this%irow(4,irxn) = this%no2_id
  this%irow(5,irxn) = this%co2_id
  this%irow(6,irxn) = this%biomass_id
  this%stoich_row(1,irxn) = this%stoich_1_doc
  this%stoich_row(2,irxn) = this%stoich_1_nh4
  this%stoich_row(3,irxn) = this%stoich_1_no3
  this%stoich_row(4,irxn) = this%stoich_1_no2
  this%stoich_row(5,irxn) = this%stoich_1_co2
  this%stoich_row(6,irxn) = this%stoich_1_biomass
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%doc_id
  this%icol(2,irxn) = this%no3_id
  this%icol(3,irxn) = this%no2_id
  this%icol(4,irxn) = this%o2_id
  ! NO2- -> N2
  irxn = 2
  this%nrow(irxn) = 6
  this%irow(1,irxn) = this%doc_id
  this%irow(2,irxn) = this%nh4_id
  this%irow(3,irxn) = this%no2_id
  this%irow(4,irxn) = this%n2_id
  this%irow(5,irxn) = this%co2_id
  this%irow(6,irxn) = this%biomass_id
  this%stoich_row(1,irxn) = this%stoich_2_doc
  this%stoich_row(2,irxn) = this%stoich_2_nh4
  this%stoich_row(3,irxn) = this%stoich_2_no2
  this%stoich_row(4,irxn) = this%stoich_2_n2
  this%stoich_row(5,irxn) = this%stoich_2_co2
  this%stoich_row(6,irxn) = this%stoich_2_biomass
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%doc_id
  this%icol(2,irxn) = this%no3_id
  this%icol(3,irxn) = this%no2_id
  this%icol(4,irxn) = this%o2_id
  if (this%nrxn > 2) then
  ! O2 -> H2O + CO2
  irxn = 3 
  this%nrow(irxn) = 5
  this%irow(1,irxn) = this%doc_id
  this%irow(2,irxn) = this%nh4_id
  this%irow(3,irxn) = this%o2_id
  this%irow(4,irxn) = this%co2_id
  this%irow(5,irxn) = this%biomass_id
  this%stoich_row(1,irxn) = this%stoich_3_doc
  this%stoich_row(2,irxn) = this%stoich_3_nh4
  this%stoich_row(3,irxn) = this%stoich_3_o2
  this%stoich_row(4,irxn) = this%stoich_3_co2
  this%stoich_row(5,irxn) = this%stoich_3_biomass
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%doc_id
  this%icol(2,irxn) = this%no3_id
  this%icol(3,irxn) = this%no2_id
  this%icol(4,irxn) = this%o2_id
  endif
  
end subroutine CyberSetup

! ************************************************************************** !

subroutine CyberReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_cyber_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  PetscReal :: kg_water
  
  PetscInt :: i, j, irxn

  PetscReal :: Co2, Cno3, Cno2, Cn2, Cdoc, Cco2, Cnh4, X
  PetscReal :: r1docmonod, r1docmonod_denom
  PetscReal :: r2docmonod, r2docmonod_denom
  PetscReal :: r3docmonod, r3docmonod_denom
  PetscReal :: r1no3monod, r1no3monod_denom
  PetscReal :: r2no2monod, r2no2monod_denom
  PetscReal :: r3o2monod, r3o2monod_denom
  PetscReal :: r1kin, r2kin, r3kin
  PetscReal :: sumkin, sumkinsq
  PetscReal :: u1, u2, u3
  PetscReal :: dr1kin_ddoc, dr1kin_dno3
  PetscReal :: dr2kin_ddoc, dr2kin_dno2
  PetscReal :: dr3kin_ddoc, dr3kin_do2
  PetscReal :: du_denom_dr
  PetscReal :: du1_dr1kin, du1_dr2kin, du1_dr3kin
  PetscReal :: du2_dr1kin, du2_dr2kin, du2_dr3kin
  PetscReal :: du3_dr1kin, du3_dr2kin, du3_dr3kin
  PetscReal :: dr1_ddoc, dr1_dno3, dr1_dno2, dr1_do2
  PetscReal :: dr2_ddoc, dr2_dno3, dr2_dno2, dr2_do2
  PetscReal :: dr3_ddoc, dr3_dno3, dr3_dno2, dr3_do2
  PetscReal :: du1_ddoc, du1_dno3, du1_dno2, du1_do2
  PetscReal :: du2_ddoc, du2_dno3, du2_dno2, du2_do2
  PetscReal :: du3_ddoc, du3_dno3, du3_dno2, du3_do2
  PetscReal :: molality_to_molarity

  PetscReal :: rate(3), derivative_col(6,3)
  
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! m^3 -> L
  kg_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
             global_auxvar%den_kg(iphase)*material_auxvar%volume

  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
    
  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    option%io_buffer = 'Activity coefficients not currently supported in &
      &CyberReact().'
    call printErrMsg(option)
  endif
  
  ! concentrations are molarities [M]
  Co2 = rt_auxvar%pri_molal(this%o2_id)* &
        rt_auxvar%pri_act_coef(this%o2_id)*molality_to_molarity
  Cnh4 = rt_auxvar%pri_molal(this%nh4_id)* &
         rt_auxvar%pri_act_coef(this%nh4_id)*molality_to_molarity
  Cno3 = rt_auxvar%pri_molal(this%no3_id)* &
         rt_auxvar%pri_act_coef(this%no3_id)*molality_to_molarity
  Cno2 = rt_auxvar%pri_molal(this%no2_id)* &
         rt_auxvar%pri_act_coef(this%no2_id)*molality_to_molarity
  Cn2 = rt_auxvar%pri_molal(this%n2_id)* &
        rt_auxvar%pri_act_coef(this%n2_id)*molality_to_molarity
  Cdoc = rt_auxvar%pri_molal(this%doc_id)* &
         rt_auxvar%pri_act_coef(this%doc_id)*molality_to_molarity
!  X = rt_auxvar%immobile(this%biomass_id-reaction%offset_immobile)
  X = rt_auxvar%pri_molal(this%biomass_id)* &
      rt_auxvar%pri_act_coef(this%biomass_id)*molality_to_molarity
  
  ! NO3- -> NO2-
  r1docmonod_denom = Cdoc+this%Kd1
  r1docmonod = Cdoc/r1docmonod_denom
  r1no3monod_denom = Cno3+this%Ka1
  r1no3monod = Cno3/r1no3monod_denom
  r1kin = this%k1*r1docmonod*r1no3monod
                 
  ! NO2- -> N2
  r2docmonod_denom = Cdoc+this%Kd2 
  r2docmonod = Cdoc/r2docmonod_denom
  r2no2monod_denom = Cno2+this%Ka2
  r2no2monod = Cno2/r2no2monod_denom
  r2kin = this%k2*r2docmonod*r2no2monod
                 
  ! O2 -> 
  r3docmonod_denom = Cdoc+this%Kd3
  r3docmonod = Cdoc/r3docmonod_denom
  r3o2monod_denom = Co2+this%Ka3
  r3o2monod = Co2/r3o2monod_denom
  r3kin = this%k3*r3docmonod*r3o2monod
                 
  sumkin = r1kin + r2kin + r3kin
  sumkinsq = sumkin * sumkin
                 
  u1 = 0.d0
  if (r1kin > 0.d0) u1 = r1kin/sumkin
  u2 = 0.d0
  if (r2kin > 0.d0) u2 = r2kin/sumkin
  u3 = 0.d0
  if (r3kin > 0.d0) u3 = r3kin/sumkin
  
  rate(1) = u1*r1kin  ! mol/mol biomass/sec
  rate(2) = u2*r2kin
  rate(3) = u3*r3kin
  
  do irxn = 1, this%nrxn
    do i = 1, this%nrow(irxn)
      ! mol/sec
      ! X is in [M]
      Residual(this%irow(i,irxn)) = Residual(this%irow(i,irxn)) - &
        this%stoich_row(i,irxn) * rate(irxn) * X * &
        ! if biomass is aqueous multiply by L_water
        ! if biomass is immobile multiply by material_auxvar%volume
        L_water
    enddo
  enddo
  
  ! decay of biomass
  ! if biomass is aqueous multiply by L_water
  ! if biomass is immobile multiply by material_auxvar%volume
  Residual(this%biomass_id) = Residual(this%biomass_id) + &
                              this%k_deg * X * L_water

  ! production of doc by biomass decay
  ! note the addition
  ! mol/sec
  Residual(this%doc_id) = Residual(this%doc_id) - &
                          this%k_deg/this%f_act * X * L_water
                 
  if (compute_derivative) then
  
    dr1kin_ddoc = this%k1 * &
                  (r1docmonod/Cdoc - r1docmonod/r1docmonod_denom) * &
                  r1no3monod * &
                  rt_auxvar%pri_act_coef(this%doc_id)
    dr1kin_dno3 = this%k1 * &
                  r1docmonod * &
                  (r1no3monod/Cno3 - r1no3monod/r1no3monod_denom) * &
                  rt_auxvar%pri_act_coef(this%no3_id)
    dr2kin_ddoc = this%k2 * &
                  (r2docmonod/Cdoc - r2docmonod/r2docmonod_denom) * &
                  r2no2monod * &
                  rt_auxvar%pri_act_coef(this%doc_id)
    dr2kin_dno2 = this%k2 * &
                  r2docmonod * &
                  (r2no2monod/Cno2 - r2no2monod/r2no2monod_denom) * &
                  rt_auxvar%pri_act_coef(this%no2_id)
    dr3kin_ddoc = this%k3 * &
                  (r3docmonod/Cdoc - r3docmonod/r3docmonod_denom) * &
                  r3o2monod * &
                  rt_auxvar%pri_act_coef(this%doc_id)
    dr3kin_do2 = this%k3 * &
                  r3docmonod * &
                  (r3o2monod/Co2 - r3o2monod/r3o2monod_denom) * &
                  rt_auxvar%pri_act_coef(this%o2_id)
                
    du_denom_dr = -1.d0/sumkinsq
    du1_dr1kin = 1.d0/sumkin + r1kin*du_denom_dr
    du1_dr2kin = r1kin*du_denom_dr
    du1_dr3kin = r1kin*du_denom_dr
    du2_dr1kin = r2kin*du_denom_dr
    du2_dr2kin = 1.d0/sumkin + r2kin*du_denom_dr
    du2_dr3kin = r2kin*du_denom_dr
    du3_dr1kin = r3kin*du_denom_dr
    du3_dr2kin = r3kin*du_denom_dr
    du3_dr3kin = 1.d0/sumkin + r3kin*du_denom_dr

    du1_ddoc = du1_dr1kin * dr1kin_ddoc + du1_dr2kin * dr2kin_ddoc + &
               du1_dr3kin * dr3kin_ddoc
    du1_dno3 = du1_dr1kin * dr1kin_dno3
    du1_dno2 = du1_dr2kin * dr2kin_dno2
    du1_do2 = du1_dr3kin * dr3kin_do2
    du2_ddoc = du2_dr1kin * dr1kin_ddoc + du2_dr2kin * dr2kin_ddoc + &
               du2_dr3kin * dr3kin_ddoc
    du2_dno3 = du2_dr1kin * dr1kin_dno3
    du2_dno2 = du2_dr2kin * dr2kin_dno2
    du2_do2 = du2_dr3kin * dr3kin_do2
    du3_ddoc = du3_dr1kin * dr1kin_ddoc + du3_dr2kin * dr2kin_ddoc + &
               du3_dr3kin * dr3kin_ddoc
    du3_dno3 = du3_dr1kin * dr1kin_dno3
    du3_dno2 = du3_dr2kin * dr2kin_dno2
    du3_do2 = du3_dr3kin * dr3kin_do2

    dr1_ddoc = dr1kin_ddoc * u1 + du1_ddoc * r1kin
    dr1_dno3 = dr1kin_dno3 * u1 + du1_dno3 * r1kin
    dr1_dno2 = du1_dno2 * r1kin
    dr1_do2 = du1_do2 * r1kin
    dr2_ddoc = dr2kin_ddoc * u2 + du2_ddoc * r2kin
    dr2_dno3 = du2_dno3 * r2kin
    dr2_dno2 = dr2kin_dno2 * u2 + du2_dno2 * r2kin
    dr2_do2 = du2_do2 * r2kin
    dr3_ddoc = dr3kin_ddoc * u3 + du3_ddoc * r3kin
    dr3_dno3 = du3_dno3 * r3kin
    dr3_dno2 = du3_dno2 * r3kin
    dr3_do2 = dr3kin_do2 * u3 + du3_do2 * r3kin
  
    irxn = 1
    derivative_col(1,irxn) = dr1_ddoc
    derivative_col(2,irxn) = dr1_dno3
    derivative_col(3,irxn) = dr1_dno2
    derivative_col(4,irxn) = dr1_do2  
    irxn = 2
    derivative_col(1,irxn) = dr2_ddoc
    derivative_col(2,irxn) = dr2_dno3
    derivative_col(3,irxn) = dr2_dno2
    derivative_col(4,irxn) = dr2_do2  
    irxn = 3
    derivative_col(1,irxn) = dr3_ddoc
    derivative_col(2,irxn) = dr3_dno3
    derivative_col(3,irxn) = dr3_dno2
    derivative_col(4,irxn) = dr3_do2      
    
    ! fill the Jacobian
    ! units = kg water/sec. Multiply by kg_water
    do irxn = 1, this%nrxn
      do j = 1, this%ncol(irxn)
        do i = 1, this%nrow(irxn)
          Jacobian(this%irow(i,irxn),this%icol(j,irxn)) = &
            Jacobian(this%irow(i,irxn),this%icol(j,irxn)) - &
            this%stoich_row(i,irxn) * derivative_col(j,irxn) * X * kg_water
        enddo
      enddo
      ! if biomass is aqueous, units = kg water/sec. Multiply by kg_water
      ! if biomass is immobile, units = m^3 bulk/sec. Multiply by 
      !   material_auxvar%volume
      do i = 1, this%nrow(irxn)
        Jacobian(this%irow(i,irxn),this%biomass_id) = &
          Jacobian(this%irow(i,irxn),this%biomass_id) - &
           this%stoich_row(i,irxn) * rate(irxn) * kg_water
      enddo
    enddo 

    ! decay of biomass
    ! if biomass is aqueous, units = kg water/sec. Multiply by kg_water
    ! if biomass is immobile, units = m^3 bulk/sec. Multiply by 
    !   material_auxvar%volume
    Jacobian(this%biomass_id,this%biomass_id) = &
      Jacobian(this%biomass_id,this%biomass_id) + &
      this%k_deg * kg_water

    ! production of doc by biomass decay
    ! units = kg water/sec. Multiply by kg_water
    Jacobian(this%doc_id,this%biomass_id) = &
      Jacobian(this%doc_id,this%biomass_id) - &
      this%k_deg/this%f_act * kg_water
    
  endif
  
end subroutine CyberReact

! ************************************************************************** !

subroutine CyberDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_cyber_type) :: this  

  call DeallocateArray(this%nrow)
  call DeallocateArray(this%irow)
  call DeallocateArray(this%icol)
  call DeallocateArray(this%stoich_row)

end subroutine CyberDestroy

end module Reaction_Sandbox_Cyber_class
