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
    PetscReal :: f_act
    PetscReal :: k_deg
    PetscReal :: k1
    PetscReal :: k2
    PetscReal :: k3
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
  
  PetscReal, parameter :: STOICH_CH2O = 1.d0

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
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CYBERNETIC')
    call StringToUpper(word)   

    select case(trim(word))
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
  
  this%f1 = 0.497d0
  this%k1 = 25.04d0 * per_day_to_per_sec
  this%Kd1 = 0.25d-3
  this%Ka1 = 1.d-6
  this%stoich_1_nh4 = 0.2d0*(1.d0-this%f1)
  this%stoich_1_no3 = 2.d0*this%f1
  this%stoich_1_no2 = 2.d0*this%f1
  this%stoich_1_co2 = this%f1
  this%stoich_1_biomass = 0.2d0*(1.d0-this%f1)
  
  this%f2 = 0.999d0
  this%k2 = 17.82d0 * per_day_to_per_sec
  this%Kd2 = 0.25d-3
  this%Ka2 = 4.d-6
  this%stoich_2_nh4 = 0.2d0*(1.d0-this%f2)
  this%stoich_2_no2 = 4.d0/3.d0*this%f2
  this%stoich_2_n2 = 2.d0/3.d0*this%f2
  this%stoich_2_co2 = this%f2
  this%stoich_2_biomass = 0.2d0*(1.d0-this%f2)
  
  this%f3 = 0.066d0
  this%k3 = 75.12d0 * per_day_to_per_sec
  this%Kd3 = 0.25d-3
  this%Ka3 = 1.d-6
  this%stoich_3_nh4 = 0.2d0*(1.d0-this%f3)
  this%stoich_3_o2 = 2.d0*this%f3
  this%stoich_3_co2 = this%f3
  this%stoich_3_biomass = 0.2d0*(1.d0-this%f3)

  this%nrxn = 3
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
  this%stoich_row(1,irxn) = STOICH_CH2O
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
  this%stoich_row(1,irxn) = STOICH_CH2O
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
  ! O2 -> H2O + CO2
  irxn = 3 
  this%nrow(irxn) = 5
  this%irow(1,irxn) = this%doc_id
  this%irow(2,irxn) = this%nh4_id
  this%irow(3,irxn) = this%o2_id
  this%irow(4,irxn) = this%co2_id
  this%irow(5,irxn) = this%biomass_id
  this%stoich_row(1,irxn) = STOICH_CH2O
  this%stoich_row(2,irxn) = this%stoich_3_nh4
  this%stoich_row(3,irxn) = this%stoich_3_o2
  this%stoich_row(4,irxn) = this%stoich_3_co2
  this%stoich_row(5,irxn) = this%stoich_3_biomass
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%doc_id
  this%icol(2,irxn) = this%no3_id
  this%icol(3,irxn) = this%no2_id
  this%icol(4,irxn) = this%o2_id
  
  this%f_act = 0.126d0
  this%k_deg = 0.532d0 * per_day_to_per_sec
      
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

  PetscReal :: rate(3), derivative_col(6,3)
  
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
    
  Co2 = rt_auxvar%pri_molal(this%o2_id)*rt_auxvar%pri_act_coef(this%o2_id)
  Cnh4 = rt_auxvar%pri_molal(this%o2_id)*rt_auxvar%pri_act_coef(this%nh4_id)
  Cno3 = rt_auxvar%pri_molal(this%no3_id)*rt_auxvar%pri_act_coef(this%no3_id)
  Cno2 = rt_auxvar%pri_molal(this%no2_id)*rt_auxvar%pri_act_coef(this%no2_id)
  Cn2 = rt_auxvar%pri_molal(this%n2_id)*rt_auxvar%pri_act_coef(this%n2_id)
  Cdoc = rt_auxvar%pri_molal(this%doc_id)*rt_auxvar%pri_act_coef(this%doc_id)
!  X = rt_auxvar%immobile(this%biomass_id-reaction%offset_immobile)
  X = rt_auxvar%pri_molal(this%biomass_id)
  
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
  r2no2monod = Cno3/r2no2monod_denom
  r2kin = this%k2*r2docmonod*r2no2monod
                 
  ! O2 -> 
  r3docmonod_denom = Cdoc+this%Kd3
  r3docmonod = Cdoc/r3docmonod_denom
  r3o2monod_denom = Co2+this%Ka3
  r3o2monod = Co2/r3o2monod_denom
  r3kin = this%k3*r3docmonod*r3o2monod
                 
  sumkin = r1kin + r2kin + r3kin
  sumkinsq = sumkin * sumkin
                 
  u1 = r1kin/sumkin
  u2 = r2kin/sumkin
  u3 = r3kin/sumkin
  
  rate(1) = u1*r1kin
  rate(2) = u2*r2kin
  rate(3) = u3*r3kin
  
  do irxn = 1, this%nrxn
    do i = 1, this%nrow(irxn)
      Residual(this%irow(i,irxn)) = Residual(this%irow(i,irxn)) + &
        this%stoich_row(i,irxn) * rate(irxn) * X
    enddo
  enddo
  
  ! decay of doc
  Residual(this%doc_id) = Residual(this%doc_id) + &
                          -1.d0 * this%k_deg/this%f_act * X
  
  ! decay of biomass
  Residual(this%biomass_id) = Residual(this%biomass_id) + &
                              -1.d0 * this%k_deg * X
                 
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
    do irxn = 1, this%nrxn
      do j = 1, this%ncol(irxn)
        do i = 1, this%nrow(irxn)
          Jacobian(this%irow(i,irxn),this%icol(j,irxn)) = &
            Jacobian(this%irow(i,irxn),this%icol(j,irxn)) + &
            this%stoich_row(i,irxn) * derivative_col(j,irxn) * X
        enddo
      enddo
      do i = 1, this%nrow(irxn)
        Jacobian(this%irow(i,irxn),this%biomass_id) = &
          Jacobian(this%irow(i,irxn),this%biomass_id) + &
           this%stoich_row(i,irxn) * rate(irxn)
      enddo
    enddo 

    ! decay of doc
    Jacobian(this%doc_id,this%doc_id) = &
      Jacobian(this%doc_id,this%doc_id) + &
      -1.d0 * this%k_deg/this%f_act

    ! decay of biomass
    Jacobian(this%biomass_id,this%biomass_id) = &
      Jacobian(this%biomass_id,this%biomass_id) + &
      -1.d0 * this%k_deg
    
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
