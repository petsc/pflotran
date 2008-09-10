module Reaction_module
  
  implicit none
  
  private 

#include "definitions.h"
  
  type, public :: aq_species_type
    character(len=MAXNAMELENGTH) :: spec_name
    PetscReal :: a0
    PetscReal :: molar_weight
    PetscReal :: valence
    type(equilibrium_rxn_type), pointer :: eqrxn
    type(aq_species_type), pointer :: next
  end type aq_species_type

  type, public :: gas_species_type
    character(len=MAXNAMELENGTH) :: spec_name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    type(equilibrium_rxn_type), pointer :: eqrxn
    type(gas_species_type), pointer :: next    
  end type gas_species_type

  type, public :: equilibrium_rxn_type
    PetscInt :: nspec
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK
  end type equilibrium_rxn_type

  type, public :: kinetic_rxn_type
    PetscInt :: nspec
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal :: logK
    PetscReal :: rate_forward
    PetscReal :: rate_reverse
  end type kinetic_rxn_type

  type, public :: mineral_type
    character(len=MAXNAMELENGTH) :: mnrl_name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    type(transition_state_rxn_type), pointer :: tstrxn
    type(mineral_type), pointer :: next
  end type mineral_type

  type, public :: transition_state_rxn_type
    PetscInt :: nspec
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: stoich(:)
    PetscInt :: nspec_primary_prefactor
    character(len=MAXNAMELENGTH), pointer :: spec_name_primary_prefactor(:)
    PetscInt, pointer :: spec_ids_primary_prefactor(:)
    PetscReal, pointer :: stoich_primary_prefactor(:)
    PetscInt :: nspec_secondary_prefactor
    character(len=MAXNAMELENGTH), pointer :: spec_name_secondary_prefactor(:)
    PetscInt, pointer :: spec_ids_secondary_prefactor(:)
    PetscReal, pointer :: stoich_secondary_prefactor(:)
    PetscReal :: affinity_factor_sigma
    PetscReal :: affinity_factor_beta
    PetscReal :: logK
    PetscReal :: rate
    PetscReal :: area0
  end type transition_state_rxn_type
  
  type, public :: ion_exchange_rxn_type
    character(len=MAXNAMELENGTH) :: reference_cation_name
    character(len=MAXNAMELENGTH) :: mnrl_name
    PetscInt :: ncation
    character(len=MAXNAMELENGTH), pointer :: cation_name(:)
    PetscInt, pointer :: cation_ids(:)
    PetscReal, pointer :: k(:)  ! selectivity coefficient
    PetscReal :: CEC
    type (ion_exchange_rxn_type), pointer :: next
  end type ion_exchange_rxn_type

  type, public :: surface_complexation_rxn_type
    character(len=MAXNAMELENGTH) :: surfcmplx_name
    PetscInt :: ncomp
    character(len=MAXNAMELENGTH), pointer :: spec_name(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: stoich(:)
    PetscReal, pointer :: free_site_stoich
    PetscReal :: logK
    PetscReal :: valence
    type (surface_complexation_rxn_type), pointer :: next
  end type surface_complexation_rxn_type
  
  type, public :: reaction_type
    type(aq_species_type), pointer :: primary_species_list
    type(aq_species_type), pointer :: secondary_species_list
    type(gas_species_type), pointer :: gas_species_list
    type(mineral_type), pointer :: mineral_list
    type(ion_exchange_rxn_type), pointer :: ion_exchange_list
    type(surface_complexation_rxn_type), pointer :: surface_complex_list
    ! compressed arrays for efficient computation
    ! primary aqueous complexes
    PetscReal, pointer :: primary_spec_molwt(:)
    PetscReal, pointer :: primary_spec_a0(:)
    PetscReal, pointer :: primary_spec_Z(:)
    ! aqueous complexes
    PetscInt :: neqcmplx
    PetscInt, pointer :: eqcmplxspecid(:,:)   ! (0:ncomp in rxn)
    PetscReal, pointer :: eqcmplxstoich(:,:)
    PetscReal, pointer :: eqcmplx_a0(:)  ! DeBye-Huckel constant
    PetscReal, pointer :: eqcmplx_K(:)
    ! ionx exchange reactions
    PetscInt, pointer :: eqionx_ncation(:)
    PetscReal, pointer :: eqionx_CEC(:)
    PetscReal, pointer :: eqionx_k(:,:)
    PetscInt, pointer :: eqionx_cationid(:)
    PetscInt, pointer :: eqionx_rxn_offset(:)
    PetscInt, pointer :: kinionx_ncation(:)
    PetscReal, pointer :: kinionx_CEC(:)
    PetscReal, pointer :: kinionx_k(:,:)
    PetscInt, pointer :: kinionx_cationid(:)
    PetscInt, pointer :: kinionx_rxn_offset(:)
    ! surface complexation reactions
    PetscInt, pointer :: eqsurfcmplxspecid(:,:)
    PetscReal, pointer :: eqsurfcmplxstoich(:,:)
    PetscReal, pointer :: eqsurfcmplx_freesite_stoich(:,:)
    PetscReal, pointer :: eqsurfcmplx_K(:)
    PetscReal, pointer :: eqsurfcmplx_Z(:)  ! valence
    PetscInt, pointer :: kinsurfcmplxspecid(:,:)
    PetscReal, pointer :: kinsurfcmplxstoich(:,:)
    PetscReal, pointer :: kinsurfcmplx_freesite_stoich(:,:)
    PetscReal, pointer :: kinsurfcmplx_K(:)
    PetscReal, pointer :: kinsurfcmplx_Z(:)  ! valence
    ! mineral reactions
      ! for saturation states
    PetscInt, pointer :: mnrlspecid(:,:)
    PetscReal, pointer :: mnrlstoich(:,:)
    PetscReal, pointer :: mnrl_K(:)
    PetscReal, pointer :: mnrl_molar_vol(:)
      ! for kinetic reactions
    PetscInt :: nkinmnrl
    PetscInt, pointer :: kinmnrlspecid(:,:)
    PetscReal, pointer :: kinmnrlstoich(:,:)
    PetscReal, pointer :: kinmnrl_K(:)
    PetscReal, pointer :: kinmnrl_rate(:,:)
    PetscInt, pointer :: kinmnrl_num_prefactors(:)
    PetscInt, pointer :: kinmnrl_pri_prefactor_id(:,:,:)
    PetscReal, pointer :: kinmnrl_pri_pref_alpha_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_pri_pref_beta_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_pri_pref_atten_coef(:,:,:)
    PetscInt, pointer :: kinmnrl_sec_prefactor_id(:,:,:)
    PetscReal, pointer :: kinmnrl_sec_pref_alpha_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_sec_pref_beta_stoich(:,:,:)
    PetscReal, pointer :: kinmnrl_sec_pref_atten_coef(:,:,:)
    PetscReal, pointer :: kinmnrl_Tempkin_const(:)
    PetscReal, pointer :: kinmnrl_affinity_power(:)
  end type reaction_type

  public :: ReactionCreate, &
            ReactionRead, &
            GetPrimarySpeciesCount, &
            GetPrimarySpeciesNames, &
            GetSecondarySpeciesCount, &
            GetSecondarySpeciesNames, &
            GetMineralCount, &
            GetMineralNames, &
            CarbonateTestProblemCreate

contains

! ************************************************************************** !
!
! CarbonateTestProblemCreate: Creates a carbonate test problem for reactive
!                             transport
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
function CarbonateTestProblemCreate(option)

  use Option_module
  
  type(reaction_type), pointer :: CarbonateTestProblemCreate
  type(option_type) :: option

  type(reaction_type), pointer :: reaction
  PetscInt :: icomp, irxn
  
  reaction => ReactionCreate()

  
  ! Assumes primary components
  ! 1 H+
  ! 2 HCO3-
  ! 3 Ca+2
  
  ! aqueous complexes
  ! CO2(aq) (combined with H2CO3(aq)
  ! CO3-2
  ! CaCO3(aq)
  
  ! minerals
  ! CaCO3(s)
  
  reaction%neqcmplx = 5
  allocate(reaction%eqcmplxspecid(0:option%ncomp,reaction%neqcmplx))
  reaction%eqcmplxspecid = 0
  allocate(reaction%eqcmplxstoich(option%ncomp,reaction%neqcmplx))
  reaction%eqcmplxstoich = 0.d0
  allocate(reaction%eqcmplx_a0(reaction%neqcmplx))
  reaction%eqcmplx_a0 = 0.d0
  allocate(reaction%eqcmplx_K(reaction%neqcmplx))
  reaction%eqcmplx_K = 0.d0
  
  ! CO2(aq)
  irxn = 1
  reaction%eqcmplxspecid(0,irxn) = 2
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  reaction%eqcmplxstoich(1,irxn) = 1.d0 ! H+
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplx_K(irxn) = -6.3447d0
  
  ! CO3-2
  irxn = 2
  reaction%eqcmplxspecid(0,irxn) = 2
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  reaction%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplx_K(irxn) = 10.3288d0
  
  ! CaCO3(aq)
  irxn = 3
  reaction%eqcmplxspecid(0,irxn) = 3
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  reaction%eqcmplxspecid(3,irxn) = 3    ! Ca+2
  reaction%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplxstoich(3,irxn) = 1.d0 ! Ca+2
  reaction%eqcmplx_K(irxn) = 7.0017d0

  ! CaHCO3-
  irxn = 4
  reaction%eqcmplxspecid(0,irxn) = 2
  reaction%eqcmplxspecid(1,irxn) = 2    ! HCO3-
  reaction%eqcmplxspecid(2,irxn) = 3    ! Ca+2
  reaction%eqcmplxstoich(1,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! Ca+2
  reaction%eqcmplx_K(irxn) = -1.0467d0

  ! OH-
  irxn = 5
  reaction%eqcmplxspecid(0,irxn) = 1
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  reaction%eqcmplx_K(irxn) = 13.9951
  
  reaction%nkinmnrl = 1
  allocate(reaction%kinmnrlspecid(0:option%ncomp,reaction%nkinmnrl))
  allocate(reaction%kinmnrlstoich(option%ncomp,reaction%nkinmnrl))
  allocate(reaction%kinmnrl_K(reaction%nkinmnrl))
  allocate(reaction%kinmnrl_rate(1,reaction%nkinmnrl))
  allocate(reaction%mnrl_molar_vol(reaction%nkinmnrl))
  allocate(reaction%kinmnrl_num_prefactors(reaction%nkinmnrl))
  
  ! CaCO3(s)
  irxn = 1
  reaction%kinmnrlspecid(0,irxn) = 3
  reaction%kinmnrlspecid(1,irxn) = 1    ! H+
  reaction%kinmnrlspecid(2,irxn) = 2    ! HCO3-
  reaction%kinmnrlspecid(3,irxn) = 3    ! Ca+2
  reaction%kinmnrlstoich(1,irxn) = -1.d0 ! H+
  reaction%kinmnrlstoich(2,irxn) = 1.d0 ! HCO3-
  reaction%kinmnrlstoich(3,irxn) = 1.d0 ! Ca+2
  reaction%kinmnrl_K(irxn) = 1.8487d0
  reaction%kinmnrl_rate(1,irxn) = 1.d-6
  reaction%mnrl_molar_vol(irxn) = 36.9340d0/1.d6  ! based on 36.934 cm^3/mol
  reaction%kinmnrl_num_prefactors(irxn) = 0
  
  CarbonateTestProblemCreate => reaction
     
end function CarbonateTestProblemCreate

! ************************************************************************** !
!
! ReactionCreate: Allocate and initialize reaction object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function ReactionCreate()

  use Option_module

  implicit none
  
  type(reaction_type), pointer :: ReactionCreate
  
  type(reaction_type), pointer :: reaction

  allocate(reaction)  

  nullify(reaction%primary_species_list)
  nullify(reaction%secondary_species_list)
  nullify(reaction%gas_species_list)
  nullify(reaction%mineral_list)
  nullify(reaction%ion_exchange_list)
  nullify(reaction%surface_complex_list)
  
  nullify(reaction%primary_spec_molwt)
  nullify(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_Z)
  
  reaction%neqcmplx = 0
  nullify(reaction%eqcmplxspecid)
  nullify(reaction%eqcmplxstoich)
  nullify(reaction%eqcmplx_a0)
  nullify(reaction%eqcmplx_K)
  
  nullify(reaction%eqionx_ncation)
  nullify(reaction%eqionx_CEC)
  nullify(reaction%eqionx_k)
  nullify(reaction%eqionx_cationid)
  nullify(reaction%eqionx_rxn_offset)
  
  nullify(reaction%kinionx_ncation)
  nullify(reaction%kinionx_CEC)
  nullify(reaction%kinionx_k)
  nullify(reaction%kinionx_cationid)
  nullify(reaction%kinionx_rxn_offset)
  
  nullify(reaction%eqsurfcmplxspecid)
  nullify(reaction%eqsurfcmplxstoich)
  nullify(reaction%eqsurfcmplx_freesite_stoich)
  nullify(reaction%eqsurfcmplx_K)
  nullify(reaction%eqsurfcmplx_Z)
  
  nullify(reaction%kinsurfcmplxspecid)
  nullify(reaction%kinsurfcmplxstoich)
  nullify(reaction%kinsurfcmplx_freesite_stoich)
  nullify(reaction%kinsurfcmplx_K)
  nullify(reaction%kinsurfcmplx_Z)

  nullify(reaction%mnrlspecid)
  nullify(reaction%mnrlstoich)
  nullify(reaction%mnrl_K)
  nullify(reaction%mnrl_molar_vol)
  
  reaction%nkinmnrl = 0  
  nullify(reaction%kinmnrlspecid)
  nullify(reaction%kinmnrlstoich)
  nullify(reaction%kinmnrl_K)
  nullify(reaction%kinmnrl_rate)
  nullify(reaction%kinmnrl_num_prefactors)
  nullify(reaction%kinmnrl_pri_prefactor_id)
  nullify(reaction%kinmnrl_pri_pref_alpha_stoich)
  nullify(reaction%kinmnrl_pri_pref_beta_stoich)
  nullify(reaction%kinmnrl_pri_pref_atten_coef)
  nullify(reaction%kinmnrl_sec_prefactor_id)
  nullify(reaction%kinmnrl_sec_pref_alpha_stoich)
  nullify(reaction%kinmnrl_sec_pref_beta_stoich)
  nullify(reaction%kinmnrl_sec_pref_atten_coef)
  nullify(reaction%kinmnrl_Tempkin_const)
  nullify(reaction%kinmnrl_affinity_power)

  ReactionCreate => reaction
  
end function ReactionCreate

! ************************************************************************** !
!
! AqueousSpeciesCreate: Allocate and initialize an aqueous species object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function AqueousSpeciesCreate()

  use Option_module

  implicit none
  
  type(aq_species_type), pointer :: AqueousSpeciesCreate
  
  type(aq_species_type), pointer :: species

  allocate(species)  
  species%spec_name = ''
  species%a0 = 0.d0
  species%molar_weight = 0.d0
  species%valence = 0.d0
  nullify(species%eqrxn)
  nullify(species%next)

  AqueousSpeciesCreate => species
  
end function AqueousSpeciesCreate

! ************************************************************************** !
!
! GasSpeciesCreate: Allocate and initialize a gas species object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function GasSpeciesCreate()

  use Option_module

  implicit none
  
  type(gas_species_type), pointer :: GasSpeciesCreate
  
  type(gas_species_type), pointer :: species

  allocate(species)  
  species%spec_name = ''
  species%molar_volume = 0.d0
  species%molar_weight = 0.d0
  nullify(species%eqrxn)
  nullify(species%next)

  GasSpeciesCreate => species
  
end function GasSpeciesCreate

! ************************************************************************** !
!
! MineralCreate: Allocate and initialize a mineral object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function MineralCreate()

  use Option_module

  implicit none
  
  type(mineral_type), pointer :: MineralCreate
  
  type(mineral_type), pointer :: mineral

  allocate(mineral)  
  mineral%mnrl_name = ''
  mineral%molar_volume = 0.d0
  mineral%molar_weight = 0.d0
  nullify(mineral%tstrxn)
  nullify(mineral%next)
  
  MineralCreate => mineral
  
end function MineralCreate


! ************************************************************************** !
!
! GetPrimarySpeciesNames: Returns the names of primary species in an array
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetPrimarySpeciesNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetPrimarySpeciesNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = GetPrimarySpeciesCount(reaction)
  allocate(names(count))
  
  count = 1
  species => reaction%primary_species_list
  do
    if (.not.associated(species)) exit
    names(count) = species%spec_name
    count = count + 1
    species => species%next
  enddo

  GetPrimarySpeciesNames => names
  
end function GetPrimarySpeciesNames

! ************************************************************************** !
!
! GetPrimarySpeciesCount: Returns the number of primary species
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetPrimarySpeciesCount(reaction)

  implicit none
  
  integer :: GetPrimarySpeciesCount
  type(reaction_type) :: reaction

  type(aq_species_type), pointer :: species

  GetPrimarySpeciesCount = 0
  species => reaction%primary_species_list
  do
    if (.not.associated(species)) exit
    GetPrimarySpeciesCount = GetPrimarySpeciesCount + 1
    species => species%next
  enddo

end function GetPrimarySpeciesCount

! ************************************************************************** !
!
! GetSecondarySpeciesNames: Returns the names of secondary species in an array
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetSecondarySpeciesNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetSecondarySpeciesNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = GetSecondarySpeciesCount(reaction)
  allocate(names(count))
  
  count = 1
  species => reaction%secondary_species_list
  do
    if (.not.associated(species)) exit
    names(count) = species%spec_name
    count = count + 1
    species => species%next
  enddo

  GetSecondarySpeciesNames => names
  
end function GetSecondarySpeciesNames

! ************************************************************************** !
!
! GetSecondarySpeciesCount: Returns the number of secondary species
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetSecondarySpeciesCount(reaction)

  implicit none
  
  integer :: GetSecondarySpeciesCount
  type(reaction_type) :: reaction

  type(aq_species_type), pointer :: species

  GetSecondarySpeciesCount = 0
  species => reaction%secondary_species_list
  do
    if (.not.associated(species)) exit
    GetSecondarySpeciesCount = GetSecondarySpeciesCount + 1
    species => species%next
  enddo

end function GetSecondarySpeciesCount

! ************************************************************************** !
!
! GetMineralNames: Returns the names of minerals in an array
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
function GetMineralNames(reaction)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetMineralNames(:)
  type(reaction_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(mineral_type), pointer :: mineral

  count = GetMineralCount(reaction)
  allocate(names(count))
  
  count = 1
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    names(count) = mineral%mnrl_name
    count = count + 1
    mineral => mineral%next
  enddo

  GetMineralNames => names
  
end function GetMineralNames

! ************************************************************************** !
!
! GetMineralCount: Returns the number of primary species
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetMineralCount(reaction)

  implicit none
  
  integer :: GetMineralCount
  type(reaction_type) :: reaction

  type(mineral_type), pointer :: mineral

  GetMineralCount = 0
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    GetMineralCount = GetMineralCount + 1
    mineral => mineral%next
  enddo

end function GetMineralCount

! ************************************************************************** !
!
! ReactionRead: Reads chemical species
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
subroutine ReactionRead(reaction,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  type(aq_species_type), pointer :: species, prev_species
  type(gas_species_type), pointer :: gas, prev_gas
  type(mineral_type), pointer :: mineral, prev_mineral
  PetscInt :: length
  PetscErrorCode :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY', ierr)
    call fiWordToUpper(word)   

    select case(trim(word))
    
      case('PRIMARY_SPECIES')
        nullify(prev_species)
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          species => AqueousSpeciesCreate()
          call fiReadWord(string,species%spec_name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,PRIMARY_SPECIES', ierr)    
          if (.not.associated(reaction%primary_species_list)) &
            reaction%primary_species_list => species
          if (associated(prev_species)) prev_species%next => species
          prev_species => species
          nullify(species)
        enddo
      case('SECONDARY_SPECIES')
        nullify(prev_species)
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          species => AqueousSpeciesCreate()
          call fiReadWord(string,species%spec_name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,PRIMARY_SPECIES', ierr)    
          if (.not.associated(reaction%secondary_species_list)) &
            reaction%secondary_species_list => species
          if (associated(prev_species)) prev_species%next => species
          prev_species => species
          nullify(species)
        enddo
      case('GAS_SPECIES')
        nullify(prev_gas)
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          gas => GasSpeciesCreate()
          call fiReadWord(string,gas%spec_name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,GAS_SPECIES', ierr)    
          if (.not.associated(reaction%gas_species_list)) &
            reaction%gas_species_list => gas
          if (associated(prev_gas)) prev_gas%next => gas
          prev_gas => gas
          nullify(gas)
        enddo
      case('MINERALS')
        nullify(prev_mineral)
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          mineral => MineralCreate()
          call fiReadWord(string,mineral%mnrl_name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,MINERAL', ierr)    
          if (.not.associated(reaction%mineral_list)) &
            reaction%mineral_list => mineral
          if (associated(prev_mineral)) prev_mineral%next => mineral
          prev_mineral => mineral
          nullify(mineral)
        enddo
      case('END')
        exit
      case('RUN_CARBONATE')
        ! skip        
      case default
        call printErrMsg(option,'CHEMISTRY keyword: '//trim(word)//' not recognized')
    end select
  enddo
 
end subroutine ReactionRead

! ************************************************************************** !
!
! AqueousSpeciesDestroy: Deallocates an aqueous species
! author: Glenn Hammond
! date: 05/29/07
!
! ************************************************************************** !
subroutine AqueousSpeciesDestroy(species)

  implicit none
    
  type(aq_species_type), pointer :: species

  if (associated(species%eqrxn)) call EquilibriumRxnDestroy(species%eqrxn)
  deallocate(species)  
  nullify(species)

end subroutine AqueousSpeciesDestroy

! ************************************************************************** !
!
! GasSpeciesDestroy: Deallocates a gas species
! author: Glenn Hammond
! date: 05/29/07
!
! ************************************************************************** !
subroutine GasSpeciesDestroy(species)

  implicit none
    
  type(gas_species_type), pointer :: species

  if (associated(species%eqrxn)) call EquilibriumRxnDestroy(species%eqrxn)
  deallocate(species)  
  nullify(species)

end subroutine GasSpeciesDestroy

! ************************************************************************** !
!
! MineralDestroy: Deallocates a mineral
! author: Glenn Hammond
! date: 05/29/07
!
! ************************************************************************** !
subroutine MineralDestroy(mineral)

  implicit none
    
  type(mineral_type), pointer :: mineral

  if (associated(mineral%tstrxn)) call TransitionStateRxnDestroy(mineral%tstrxn)
  deallocate(mineral)  
  nullify(mineral)

end subroutine MineralDestroy

! ************************************************************************** !
!
! EquilibriumRxnDestroy: Deallocates an equilibrium reaction
! author: Glenn Hammond
! date: 05/29/07
!
! ************************************************************************** !
subroutine EquilibriumRxnDestroy(eqrxn)

  implicit none
    
  type(equilibrium_rxn_type), pointer :: eqrxn

  if (.not.associated(eqrxn)) return
  
  if (associated(eqrxn%spec_name)) deallocate(eqrxn%spec_name)
  nullify(eqrxn%spec_name)
  if (associated(eqrxn%spec_ids)) deallocate(eqrxn%spec_ids)
  nullify(eqrxn%spec_ids)
  if (associated(eqrxn%stoich)) deallocate(eqrxn%stoich)
  nullify(eqrxn%stoich)

  deallocate(eqrxn)  
  nullify(eqrxn)

end subroutine EquilibriumRxnDestroy

! ************************************************************************** !
!
! TransitionStateRxnDestroy: Deallocates a transition state reaction
! author: Glenn Hammond
! date: 05/29/07
!
! ************************************************************************** !
subroutine TransitionStateRxnDestroy(tstrxn)

  implicit none
    
  type(transition_state_rxn_type), pointer :: tstrxn

  if (.not.associated(tstrxn)) return
  
  if (associated(tstrxn%spec_name)) deallocate(tstrxn%spec_name)
  nullify(tstrxn%spec_name)
  if (associated(tstrxn%spec_ids)) deallocate(tstrxn%spec_ids)
  nullify(tstrxn%spec_ids)
  if (associated(tstrxn%stoich)) deallocate(tstrxn%stoich)
  nullify(tstrxn%stoich)
  if (associated(tstrxn%spec_name_primary_prefactor)) deallocate(tstrxn%spec_name_primary_prefactor)
  nullify(tstrxn%spec_name_primary_prefactor)
  if (associated(tstrxn%spec_ids_primary_prefactor)) deallocate(tstrxn%spec_ids_primary_prefactor)
  nullify(tstrxn%spec_ids_primary_prefactor)
  if (associated(tstrxn%stoich_primary_prefactor)) deallocate(tstrxn%stoich_primary_prefactor)
  nullify(tstrxn%stoich_primary_prefactor)
  if (associated(tstrxn%spec_name_secondary_prefactor)) deallocate(tstrxn%spec_name_secondary_prefactor)
  nullify(tstrxn%spec_name_secondary_prefactor)
  if (associated(tstrxn%spec_ids_secondary_prefactor)) deallocate(tstrxn%spec_ids_secondary_prefactor)
  nullify(tstrxn%spec_ids_secondary_prefactor)
  if (associated(tstrxn%stoich_secondary_prefactor)) deallocate(tstrxn%stoich_secondary_prefactor)
  nullify(tstrxn%stoich_secondary_prefactor)

  deallocate(tstrxn)  
  nullify(tstrxn)

end subroutine TransitionStateRxnDestroy

! ************************************************************************** !
!
! IonExchangeRxnDestroy: Deallocates an ion exchange reaction
! author: Glenn Hammond
! date: 05/29/07
!
! ************************************************************************** !
subroutine IonExchangeRxnDestroy(ionxrxn)

  implicit none
    
  type(ion_exchange_rxn_type), pointer :: ionxrxn

  if (.not.associated(ionxrxn)) return
  
  if (associated(ionxrxn%cation_name)) deallocate(ionxrxn%cation_name)
  nullify(ionxrxn%cation_name)
  if (associated(ionxrxn%cation_ids)) deallocate(ionxrxn%cation_ids)
  nullify(ionxrxn%cation_ids)
  if (associated(ionxrxn%k)) deallocate(ionxrxn%k)
  nullify(ionxrxn%k)

  deallocate(ionxrxn)  
  nullify(ionxrxn)

end subroutine IonExchangeRxnDestroy

! ************************************************************************** !
!
! SurfaceComplexationRxnDestroy: Deallocates a surface complexation reaction
! author: Glenn Hammond
! date: 05/29/07
!
! ************************************************************************** !
subroutine SurfaceComplexationRxnDestroy(surfcplxrxn)

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: surfcplxrxn

  if (.not.associated(surfcplxrxn)) return
  
  if (associated(surfcplxrxn%spec_name)) deallocate(surfcplxrxn%spec_name)
  nullify(surfcplxrxn%spec_name)
  if (associated(surfcplxrxn%spec_ids)) deallocate(surfcplxrxn%spec_ids)
  nullify(surfcplxrxn%spec_ids)
  if (associated(surfcplxrxn%stoich)) deallocate(surfcplxrxn%stoich)
  nullify(surfcplxrxn%stoich)

  deallocate(surfcplxrxn)  
  nullify(surfcplxrxn)

end subroutine SurfaceComplexationRxnDestroy

! ************************************************************************** !
!
! ReactionDestroy: Deallocates a reaction object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine ReactionDestroy(reaction)

  implicit none

  type(reaction_type), pointer :: reaction
  
  type(aq_species_type), pointer :: aq_species, prev_aq_species
  type(gas_species_type), pointer :: gas_species, prev_gas_species
  type(mineral_type), pointer :: mineral, prev_mineral
  type(ion_exchange_rxn_type), pointer :: ionxrxn, prev_ionxrxn
  type(surface_complexation_rxn_type), pointer :: surfcplxrxn, prev_surfcplxrxn

  if (.not.associated(reaction)) return 

  ! primary species
  aq_species => reaction%primary_species_list
  do
    if (.not.associated(aq_species)) exit
    prev_aq_species => aq_species
    aq_species => aq_species%next
    call AqueousSpeciesDestroy(prev_aq_species)
  enddo  
  nullify(reaction%primary_species_list)

  ! secondary species
  aq_species => reaction%secondary_species_list
  do
    if (.not.associated(aq_species)) exit
    prev_aq_species => aq_species
    aq_species => aq_species%next
    call AqueousSpeciesDestroy(prev_aq_species)
  enddo  
  nullify(reaction%secondary_species_list)

  ! gas species
  gas_species => reaction%gas_species_list
  do
    if (.not.associated(gas_species)) exit
    prev_gas_species => gas_species
    gas_species => gas_species%next
    call GasSpeciesDestroy(prev_gas_species)
  enddo  
  nullify(reaction%gas_species_list)
  
  ! mineral species
  mineral => reaction%mineral_list
  do
    if (.not.associated(mineral)) exit
    prev_mineral => mineral
    mineral => mineral%next
    call MineralDestroy(prev_mineral)
  enddo    
  nullify(reaction%mineral_list)
  
  ! ionx exchange reactions
  ionxrxn => reaction%ion_exchange_list
  do
    if (.not.associated(ionxrxn)) exit
    prev_ionxrxn => ionxrxn
    ionxrxn => ionxrxn%next
    call IonExchangeRxnDestroy(prev_ionxrxn)
  enddo    
  nullify(reaction%ion_exchange_list)

  ! surface complexation reactions
  surfcplxrxn => reaction%surface_complex_list
  do
    if (.not.associated(surfcplxrxn)) exit
    prev_surfcplxrxn => surfcplxrxn
    surfcplxrxn => surfcplxrxn%next
    call SurfaceComplexationRxnDestroy(prev_surfcplxrxn)
  enddo    
  nullify(reaction%surface_complex_list)
  
  if (associated(reaction%primary_spec_molwt)) deallocate(reaction%primary_spec_molwt)
  nullify(reaction%primary_spec_molwt)
  if (associated(reaction%primary_spec_a0)) deallocate(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_a0)
  if (associated(reaction%primary_spec_Z)) deallocate(reaction%primary_spec_Z)
  nullify(reaction%primary_spec_Z)
  
  if (associated(reaction%eqcmplxspecid)) deallocate(reaction%eqcmplxspecid)
  nullify(reaction%eqcmplxspecid)
  if (associated(reaction%eqcmplxstoich)) deallocate(reaction%eqcmplxstoich)
  nullify(reaction%eqcmplxstoich)
  if (associated(reaction%eqcmplx_a0)) deallocate(reaction%eqcmplx_a0)
  nullify(reaction%eqcmplx_a0)
  if (associated(reaction%eqcmplx_K)) deallocate(reaction%eqcmplx_K)
  nullify(reaction%eqcmplx_K)
  
  if (associated(reaction%eqionx_ncation)) deallocate(reaction%eqionx_ncation)
  nullify(reaction%eqionx_ncation)
  if (associated(reaction%eqionx_CEC)) deallocate(reaction%eqionx_CEC)
  nullify(reaction%eqionx_CEC)
  if (associated(reaction%eqionx_k)) deallocate(reaction%eqionx_k)
  nullify(reaction%eqionx_k)
  if (associated(reaction%eqionx_cationid)) deallocate(reaction%eqionx_cationid)
  nullify(reaction%eqionx_cationid)
  if (associated(reaction%eqionx_cationid)) deallocate(reaction%eqionx_cationid)
  nullify(reaction%eqionx_rxn_offset)
  
  if (associated(reaction%kinionx_ncation)) deallocate(reaction%kinionx_ncation)
  nullify(reaction%kinionx_ncation)
  if (associated(reaction%kinionx_CEC)) deallocate(reaction%kinionx_CEC)
  nullify(reaction%kinionx_CEC)
  if (associated(reaction%kinionx_k)) deallocate(reaction%kinionx_k)
  nullify(reaction%kinionx_k)
  if (associated(reaction%kinionx_cationid)) deallocate(reaction%kinionx_cationid)
  nullify(reaction%kinionx_cationid)
  if (associated(reaction%kinionx_rxn_offset)) deallocate(reaction%kinionx_rxn_offset)
  nullify(reaction%kinionx_rxn_offset)
  
  if (associated(reaction%eqsurfcmplxspecid)) deallocate(reaction%eqsurfcmplxspecid)
  nullify(reaction%eqsurfcmplxspecid)
  if (associated(reaction%eqsurfcmplxstoich)) deallocate(reaction%eqsurfcmplxstoich)
  nullify(reaction%eqsurfcmplxstoich)
  if (associated(reaction%eqsurfcmplx_freesite_stoich)) deallocate(reaction%eqsurfcmplx_freesite_stoich)
  nullify(reaction%eqsurfcmplx_freesite_stoich)
  if (associated(reaction%eqsurfcmplx_K)) deallocate(reaction%eqsurfcmplx_K)
  nullify(reaction%eqsurfcmplx_K)
  if (associated(reaction%eqsurfcmplx_Z)) deallocate(reaction%eqsurfcmplx_Z)
  nullify(reaction%eqsurfcmplx_Z)
  
  if (associated(reaction%kinsurfcmplxspecid)) deallocate(reaction%kinsurfcmplxspecid)
  nullify(reaction%kinsurfcmplxspecid)
  if (associated(reaction%kinsurfcmplxstoich)) deallocate(reaction%kinsurfcmplxstoich)
  nullify(reaction%kinsurfcmplxstoich)
  if (associated(reaction%kinsurfcmplx_freesite_stoich)) deallocate(reaction%kinsurfcmplx_freesite_stoich)
  nullify(reaction%kinsurfcmplx_freesite_stoich)
  if (associated(reaction%kinsurfcmplx_K)) deallocate(reaction%kinsurfcmplx_K)
  nullify(reaction%kinsurfcmplx_K)
  if (associated(reaction%kinsurfcmplx_Z)) deallocate(reaction%kinsurfcmplx_Z)
  nullify(reaction%kinsurfcmplx_Z)
  
  if (associated(reaction%mnrlspecid)) deallocate(reaction%mnrlspecid)
  nullify(reaction%mnrlspecid)
  if (associated(reaction%mnrlstoich)) deallocate(reaction%mnrlstoich)
  nullify(reaction%mnrlstoich)
  if (associated(reaction%mnrl_K)) deallocate(reaction%mnrl_K)
  nullify(reaction%mnrl_K)
  if (associated(reaction%mnrl_molar_vol)) deallocate(reaction%mnrl_molar_vol)
  nullify(reaction%mnrl_molar_vol)
  
  if (associated(reaction%kinmnrlspecid)) deallocate(reaction%kinmnrlspecid)
  nullify(reaction%kinmnrlspecid)
  if (associated(reaction%kinmnrlstoich)) deallocate(reaction%kinmnrlstoich)
  nullify(reaction%kinmnrlstoich)
  if (associated(reaction%kinmnrl_K)) deallocate(reaction%kinmnrl_K)
  nullify(reaction%kinmnrl_K)
  if (associated(reaction%kinmnrl_rate)) deallocate(reaction%kinmnrl_rate)
  nullify(reaction%kinmnrl_rate)
  if (associated(reaction%kinmnrl_num_prefactors)) deallocate(reaction%kinmnrl_num_prefactors)
  nullify(reaction%kinmnrl_num_prefactors)
  if (associated(reaction%kinmnrl_pri_prefactor_id)) deallocate(reaction%kinmnrl_pri_prefactor_id)
  nullify(reaction%kinmnrl_pri_prefactor_id)
  if (associated(reaction%kinmnrl_pri_pref_alpha_stoich)) deallocate(reaction%kinmnrl_pri_pref_alpha_stoich)
  nullify(reaction%kinmnrl_pri_pref_alpha_stoich)
  if (associated(reaction%kinmnrl_pri_pref_beta_stoich)) deallocate(reaction%kinmnrl_pri_pref_beta_stoich)
  nullify(reaction%kinmnrl_pri_pref_beta_stoich)
  if (associated(reaction%kinmnrl_pri_pref_atten_coef)) deallocate(reaction%kinmnrl_pri_pref_atten_coef)
  nullify(reaction%kinmnrl_pri_pref_atten_coef)
  if (associated(reaction%kinmnrl_sec_prefactor_id)) deallocate(reaction%kinmnrl_sec_prefactor_id)
  nullify(reaction%kinmnrl_sec_prefactor_id)
  if (associated(reaction%kinmnrl_sec_pref_alpha_stoich)) deallocate(reaction%kinmnrl_sec_pref_alpha_stoich)
  nullify(reaction%kinmnrl_sec_pref_alpha_stoich)
  if (associated(reaction%kinmnrl_sec_pref_beta_stoich)) deallocate(reaction%kinmnrl_sec_pref_beta_stoich)
  nullify(reaction%kinmnrl_sec_pref_beta_stoich)
  if (associated(reaction%kinmnrl_sec_pref_atten_coef)) deallocate(reaction%kinmnrl_sec_pref_atten_coef)
  nullify(reaction%kinmnrl_sec_pref_atten_coef)
  if (associated(reaction%kinmnrl_Tempkin_const)) deallocate(reaction%kinmnrl_Tempkin_const)
  nullify(reaction%kinmnrl_Tempkin_const)
  if (associated(reaction%kinmnrl_affinity_power)) deallocate(reaction%kinmnrl_affinity_power)
  nullify(reaction%kinmnrl_affinity_power)

end subroutine ReactionDestroy

end module Reaction_module
