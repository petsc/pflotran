module Chemistry_module
  
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
    PetscReal, pointer :: eqcmplx_a0(:)
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
    PetscReal, pointer :: kinmnrl_rate(:)
    PetscInt, pointer :: kinmnrl_pri_prefactor_id(:,:)
    PetscInt, pointer :: kinmnrl_sec_prefactor_id(:,:)
    PetscInt, pointer :: kinmnrl_pri_prefactor_stoich(:,:)
    PetscInt, pointer :: kinmnrl_sec_prefactor_stoich(:,:)
    PetscReal, pointer :: kinmnrl_affinity_factor_sigma(:)
    PetscReal, pointer :: kinmnrl_affinity_factor_beta(:)
  end type reaction_type

  public :: ChemistryCreate, &
            ChemistryRead, &
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

  type(reaction_type), pointer :: chemistry
  PetscInt :: icomp, irxn
  
  chemistry => ChemistryCreate()

  
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
  
  chemistry%neqcmplx = 5
  allocate(chemistry%eqcmplxspecid(0:option%ncomp,chemistry%neqcmplx))
  chemistry%eqcmplxspecid = 0
  allocate(chemistry%eqcmplxstoich(option%ncomp,chemistry%neqcmplx))
  chemistry%eqcmplxstoich = 0.d0
  allocate(chemistry%eqcmplx_a0(chemistry%neqcmplx))
  chemistry%eqcmplx_a0 = 0.d0
  allocate(chemistry%eqcmplx_K(chemistry%neqcmplx))
  chemistry%eqcmplx_K = 0.d0
  
  ! CO2(aq)
  irxn = 1
  chemistry%eqcmplxspecid(0,irxn) = 2
  chemistry%eqcmplxspecid(1,irxn) = 1    ! H+
  chemistry%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  chemistry%eqcmplxstoich(1,irxn) = 1.d0 ! H+
  chemistry%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  chemistry%eqcmplx_K(irxn) = -6.3447d0
  
  ! CO3-2
  irxn = 2
  chemistry%eqcmplxspecid(0,irxn) = 2
  chemistry%eqcmplxspecid(1,irxn) = 1    ! H+
  chemistry%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  chemistry%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  chemistry%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  chemistry%eqcmplx_K(irxn) = 10.3288d0
  
  ! CaCO3(aq)
  irxn = 3
  chemistry%eqcmplxspecid(0,irxn) = 3
  chemistry%eqcmplxspecid(1,irxn) = 1    ! H+
  chemistry%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  chemistry%eqcmplxspecid(3,irxn) = 3    ! Ca+2
  chemistry%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  chemistry%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  chemistry%eqcmplxstoich(3,irxn) = 1.d0 ! Ca+2
  chemistry%eqcmplx_K(irxn) = 7.0017d0

  ! CaHCO3-
  irxn = 4
  chemistry%eqcmplxspecid(0,irxn) = 2
  chemistry%eqcmplxspecid(1,irxn) = 2    ! HCO3-
  chemistry%eqcmplxspecid(2,irxn) = 3    ! Ca+2
  chemistry%eqcmplxstoich(1,irxn) = 1.d0 ! HCO3-
  chemistry%eqcmplxstoich(2,irxn) = 1.d0 ! Ca+2
  chemistry%eqcmplx_K(irxn) = -1.0467d0

  ! OH-
  irxn = 5
  chemistry%eqcmplxspecid(0,irxn) = 1
  chemistry%eqcmplxspecid(1,irxn) = 1    ! H+
  chemistry%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  chemistry%eqcmplx_K(irxn) = 13.9951
  
  chemistry%nkinmnrl = 1
  allocate(chemistry%kinmnrlspecid(0:option%ncomp,chemistry%nkinmnrl))
  allocate(chemistry%kinmnrlstoich(option%ncomp,chemistry%nkinmnrl))
  allocate(chemistry%kinmnrl_K(chemistry%nkinmnrl))
  allocate(chemistry%kinmnrl_rate(chemistry%nkinmnrl))
  allocate(chemistry%mnrl_molar_vol(chemistry%nkinmnrl))
  
  ! CaCO3(s)
  irxn = 1
  chemistry%kinmnrlspecid(0,irxn) = 3
  chemistry%kinmnrlspecid(1,irxn) = 1    ! H+
  chemistry%kinmnrlspecid(2,irxn) = 2    ! HCO3-
  chemistry%kinmnrlspecid(3,irxn) = 3    ! Ca+2
  chemistry%kinmnrlstoich(1,irxn) = -1.d0 ! H+
  chemistry%kinmnrlstoich(2,irxn) = 1.d0 ! HCO3-
  chemistry%kinmnrlstoich(3,irxn) = 1.d0 ! Ca+2
  chemistry%kinmnrl_K(irxn) = 1.8487d0
  chemistry%kinmnrl_rate(irxn) = 1.d-6
  chemistry%mnrl_molar_vol(irxn) = 36.9340d0/1.d6  ! based on 36.934 cm^3/mol
  
  CarbonateTestProblemCreate => chemistry
     
end function CarbonateTestProblemCreate

! ************************************************************************** !
!
! ChemistryCreate: Allocate and initialize chemistry object
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
function ChemistryCreate()

  use Option_module

  implicit none
  
  type(reaction_type), pointer :: ChemistryCreate
  
  type(reaction_type), pointer :: chemistry

  allocate(chemistry)  

  nullify(chemistry%primary_species_list)
  nullify(chemistry%secondary_species_list)
  nullify(chemistry%gas_species_list)
  nullify(chemistry%mineral_list)
  nullify(chemistry%ion_exchange_list)
  nullify(chemistry%surface_complex_list)
  
  nullify(chemistry%primary_spec_molwt)
  nullify(chemistry%primary_spec_a0)
  nullify(chemistry%primary_spec_Z)
  
  chemistry%neqcmplx = 0
  nullify(chemistry%eqcmplxspecid)
  nullify(chemistry%eqcmplxstoich)
  nullify(chemistry%eqcmplx_a0)
  nullify(chemistry%eqcmplx_K)
  
  nullify(chemistry%eqionx_ncation)
  nullify(chemistry%eqionx_CEC)
  nullify(chemistry%eqionx_k)
  nullify(chemistry%eqionx_cationid)
  nullify(chemistry%eqionx_rxn_offset)
  
  nullify(chemistry%kinionx_ncation)
  nullify(chemistry%kinionx_CEC)
  nullify(chemistry%kinionx_k)
  nullify(chemistry%kinionx_cationid)
  nullify(chemistry%kinionx_rxn_offset)
  
  nullify(chemistry%eqsurfcmplxspecid)
  nullify(chemistry%eqsurfcmplxstoich)
  nullify(chemistry%eqsurfcmplx_freesite_stoich)
  nullify(chemistry%eqsurfcmplx_K)
  nullify(chemistry%eqsurfcmplx_Z)
  
  nullify(chemistry%kinsurfcmplxspecid)
  nullify(chemistry%kinsurfcmplxstoich)
  nullify(chemistry%kinsurfcmplx_freesite_stoich)
  nullify(chemistry%kinsurfcmplx_K)
  nullify(chemistry%kinsurfcmplx_Z)

  nullify(chemistry%mnrlspecid)
  nullify(chemistry%mnrlstoich)
  nullify(chemistry%mnrl_K)
  nullify(chemistry%mnrl_molar_vol)
  
  chemistry%nkinmnrl = 0  
  nullify(chemistry%kinmnrlspecid)
  nullify(chemistry%kinmnrlstoich)
  nullify(chemistry%kinmnrl_K)
  nullify(chemistry%kinmnrl_rate)
  nullify(chemistry%kinmnrl_pri_prefactor_id)
  nullify(chemistry%kinmnrl_sec_prefactor_id)
  nullify(chemistry%kinmnrl_pri_prefactor_stoich)
  nullify(chemistry%kinmnrl_sec_prefactor_stoich)
  nullify(chemistry%kinmnrl_affinity_factor_sigma)
  nullify(chemistry%kinmnrl_affinity_factor_beta)

  ChemistryCreate => chemistry
  
end function ChemistryCreate

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
function GetPrimarySpeciesNames(chemistry)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetPrimarySpeciesNames(:)
  type(reaction_type) :: chemistry

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = GetPrimarySpeciesCount(chemistry)
  allocate(names(count))
  
  count = 1
  species => chemistry%primary_species_list
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
function GetPrimarySpeciesCount(chemistry)

  implicit none
  
  integer :: GetPrimarySpeciesCount
  type(reaction_type) :: chemistry

  type(aq_species_type), pointer :: species

  GetPrimarySpeciesCount = 0
  species => chemistry%primary_species_list
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
function GetSecondarySpeciesNames(chemistry)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetSecondarySpeciesNames(:)
  type(reaction_type) :: chemistry

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = GetSecondarySpeciesCount(chemistry)
  allocate(names(count))
  
  count = 1
  species => chemistry%secondary_species_list
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
function GetSecondarySpeciesCount(chemistry)

  implicit none
  
  integer :: GetSecondarySpeciesCount
  type(reaction_type) :: chemistry

  type(aq_species_type), pointer :: species

  GetSecondarySpeciesCount = 0
  species => chemistry%secondary_species_list
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
function GetMineralNames(chemistry)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: GetMineralNames(:)
  type(reaction_type) :: chemistry

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(mineral_type), pointer :: mineral

  count = GetMineralCount(chemistry)
  allocate(names(count))
  
  count = 1
  mineral => chemistry%mineral_list
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
function GetMineralCount(chemistry)

  implicit none
  
  integer :: GetMineralCount
  type(reaction_type) :: chemistry

  type(mineral_type), pointer :: mineral

  GetMineralCount = 0
  mineral => chemistry%mineral_list
  do
    if (.not.associated(mineral)) exit
    GetMineralCount = GetMineralCount + 1
    mineral => mineral%next
  enddo

end function GetMineralCount

! ************************************************************************** !
!
! ChemistryRead: Reads chemical species
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
subroutine ChemistryRead(chemistry,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(reaction_type) :: chemistry
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
          if (.not.associated(chemistry%primary_species_list)) &
            chemistry%primary_species_list => species
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
          if (.not.associated(chemistry%secondary_species_list)) &
            chemistry%secondary_species_list => species
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
          if (.not.associated(chemistry%gas_species_list)) &
            chemistry%gas_species_list => gas
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
          if (.not.associated(chemistry%mineral_list)) &
            chemistry%mineral_list => mineral
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
 
end subroutine ChemistryRead

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
! ChemistryDestroy: Deallocates a chemistry object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine ChemistryDestroy(chemistry)

  implicit none

  type(reaction_type), pointer :: chemistry
  
  type(aq_species_type), pointer :: aq_species, prev_aq_species
  type(gas_species_type), pointer :: gas_species, prev_gas_species
  type(mineral_type), pointer :: mineral, prev_mineral
  type(ion_exchange_rxn_type), pointer :: ionxrxn, prev_ionxrxn
  type(surface_complexation_rxn_type), pointer :: surfcplxrxn, prev_surfcplxrxn

  if (.not.associated(chemistry)) return 

  ! primary species
  aq_species => chemistry%primary_species_list
  do
    if (.not.associated(aq_species)) exit
    prev_aq_species => aq_species
    aq_species => aq_species%next
    call AqueousSpeciesDestroy(prev_aq_species)
  enddo  
  nullify(chemistry%primary_species_list)

  ! secondary species
  aq_species => chemistry%secondary_species_list
  do
    if (.not.associated(aq_species)) exit
    prev_aq_species => aq_species
    aq_species => aq_species%next
    call AqueousSpeciesDestroy(prev_aq_species)
  enddo  
  nullify(chemistry%secondary_species_list)

  ! gas species
  gas_species => chemistry%gas_species_list
  do
    if (.not.associated(gas_species)) exit
    prev_gas_species => gas_species
    gas_species => gas_species%next
    call GasSpeciesDestroy(prev_gas_species)
  enddo  
  nullify(chemistry%gas_species_list)
  
  ! mineral species
  mineral => chemistry%mineral_list
  do
    if (.not.associated(mineral)) exit
    prev_mineral => mineral
    mineral => mineral%next
    call MineralDestroy(prev_mineral)
  enddo    
  nullify(chemistry%mineral_list)
  
  ! ionx exchange reactions
  ionxrxn => chemistry%ion_exchange_list
  do
    if (.not.associated(ionxrxn)) exit
    prev_ionxrxn => ionxrxn
    ionxrxn => ionxrxn%next
    call IonExchangeRxnDestroy(prev_ionxrxn)
  enddo    
  nullify(chemistry%ion_exchange_list)

  ! surface complexation reactions
  surfcplxrxn => chemistry%surface_complex_list
  do
    if (.not.associated(surfcplxrxn)) exit
    prev_surfcplxrxn => surfcplxrxn
    surfcplxrxn => surfcplxrxn%next
    call SurfaceComplexationRxnDestroy(prev_surfcplxrxn)
  enddo    
  nullify(chemistry%surface_complex_list)
  
  if (associated(chemistry%primary_spec_molwt)) deallocate(chemistry%primary_spec_molwt)
  nullify(chemistry%primary_spec_molwt)
  if (associated(chemistry%primary_spec_a0)) deallocate(chemistry%primary_spec_a0)
  nullify(chemistry%primary_spec_a0)
  if (associated(chemistry%primary_spec_Z)) deallocate(chemistry%primary_spec_Z)
  nullify(chemistry%primary_spec_Z)
  
  if (associated(chemistry%eqcmplxspecid)) deallocate(chemistry%eqcmplxspecid)
  nullify(chemistry%eqcmplxspecid)
  if (associated(chemistry%eqcmplxstoich)) deallocate(chemistry%eqcmplxstoich)
  nullify(chemistry%eqcmplxstoich)
  if (associated(chemistry%eqcmplx_a0)) deallocate(chemistry%eqcmplx_a0)
  nullify(chemistry%eqcmplx_a0)
  if (associated(chemistry%eqcmplx_K)) deallocate(chemistry%eqcmplx_K)
  nullify(chemistry%eqcmplx_K)
  
  if (associated(chemistry%eqionx_ncation)) deallocate(chemistry%eqionx_ncation)
  nullify(chemistry%eqionx_ncation)
  if (associated(chemistry%eqionx_CEC)) deallocate(chemistry%eqionx_CEC)
  nullify(chemistry%eqionx_CEC)
  if (associated(chemistry%eqionx_k)) deallocate(chemistry%eqionx_k)
  nullify(chemistry%eqionx_k)
  if (associated(chemistry%eqionx_cationid)) deallocate(chemistry%eqionx_cationid)
  nullify(chemistry%eqionx_cationid)
  if (associated(chemistry%eqionx_cationid)) deallocate(chemistry%eqionx_cationid)
  nullify(chemistry%eqionx_rxn_offset)
  
  if (associated(chemistry%kinionx_ncation)) deallocate(chemistry%kinionx_ncation)
  nullify(chemistry%kinionx_ncation)
  if (associated(chemistry%kinionx_CEC)) deallocate(chemistry%kinionx_CEC)
  nullify(chemistry%kinionx_CEC)
  if (associated(chemistry%kinionx_k)) deallocate(chemistry%kinionx_k)
  nullify(chemistry%kinionx_k)
  if (associated(chemistry%kinionx_cationid)) deallocate(chemistry%kinionx_cationid)
  nullify(chemistry%kinionx_cationid)
  if (associated(chemistry%kinionx_rxn_offset)) deallocate(chemistry%kinionx_rxn_offset)
  nullify(chemistry%kinionx_rxn_offset)
  
  if (associated(chemistry%eqsurfcmplxspecid)) deallocate(chemistry%eqsurfcmplxspecid)
  nullify(chemistry%eqsurfcmplxspecid)
  if (associated(chemistry%eqsurfcmplxstoich)) deallocate(chemistry%eqsurfcmplxstoich)
  nullify(chemistry%eqsurfcmplxstoich)
  if (associated(chemistry%eqsurfcmplx_freesite_stoich)) deallocate(chemistry%eqsurfcmplx_freesite_stoich)
  nullify(chemistry%eqsurfcmplx_freesite_stoich)
  if (associated(chemistry%eqsurfcmplx_K)) deallocate(chemistry%eqsurfcmplx_K)
  nullify(chemistry%eqsurfcmplx_K)
  if (associated(chemistry%eqsurfcmplx_Z)) deallocate(chemistry%eqsurfcmplx_Z)
  nullify(chemistry%eqsurfcmplx_Z)
  
  if (associated(chemistry%kinsurfcmplxspecid)) deallocate(chemistry%kinsurfcmplxspecid)
  nullify(chemistry%kinsurfcmplxspecid)
  if (associated(chemistry%kinsurfcmplxstoich)) deallocate(chemistry%kinsurfcmplxstoich)
  nullify(chemistry%kinsurfcmplxstoich)
  if (associated(chemistry%kinsurfcmplx_freesite_stoich)) deallocate(chemistry%kinsurfcmplx_freesite_stoich)
  nullify(chemistry%kinsurfcmplx_freesite_stoich)
  if (associated(chemistry%kinsurfcmplx_K)) deallocate(chemistry%kinsurfcmplx_K)
  nullify(chemistry%kinsurfcmplx_K)
  if (associated(chemistry%kinsurfcmplx_Z)) deallocate(chemistry%kinsurfcmplx_Z)
  nullify(chemistry%kinsurfcmplx_Z)
  
  if (associated(chemistry%mnrlspecid)) deallocate(chemistry%mnrlspecid)
  nullify(chemistry%mnrlspecid)
  if (associated(chemistry%mnrlstoich)) deallocate(chemistry%mnrlstoich)
  nullify(chemistry%mnrlstoich)
  if (associated(chemistry%mnrl_K)) deallocate(chemistry%mnrl_K)
  nullify(chemistry%mnrl_K)
  if (associated(chemistry%mnrl_molar_vol)) deallocate(chemistry%mnrl_molar_vol)
  nullify(chemistry%mnrl_molar_vol)
  
  if (associated(chemistry%kinmnrlspecid)) deallocate(chemistry%kinmnrlspecid)
  nullify(chemistry%kinmnrlspecid)
  if (associated(chemistry%kinmnrlstoich)) deallocate(chemistry%kinmnrlstoich)
  nullify(chemistry%kinmnrlstoich)
  if (associated(chemistry%kinmnrl_K)) deallocate(chemistry%kinmnrl_K)
  nullify(chemistry%kinmnrl_K)
  if (associated(chemistry%kinmnrl_rate)) deallocate(chemistry%kinmnrl_rate)
  nullify(chemistry%kinmnrl_rate)
  if (associated(chemistry%kinmnrl_pri_prefactor_id)) deallocate(chemistry%kinmnrl_sec_prefactor_id)
  nullify(chemistry%kinmnrl_pri_prefactor_id)
  if (associated(chemistry%kinmnrl_sec_prefactor_id)) deallocate(chemistry%kinmnrl_sec_prefactor_id)
  nullify(chemistry%kinmnrl_sec_prefactor_id)
  if (associated(chemistry%kinmnrl_pri_prefactor_stoich)) deallocate(chemistry%kinmnrl_pri_prefactor_stoich)
  nullify(chemistry%kinmnrl_pri_prefactor_stoich)
  if (associated(chemistry%kinmnrl_sec_prefactor_stoich)) deallocate(chemistry%kinmnrl_sec_prefactor_stoich)
  nullify(chemistry%kinmnrl_sec_prefactor_stoich)
  if (associated(chemistry%kinmnrl_affinity_factor_sigma)) deallocate(chemistry%kinmnrl_affinity_factor_sigma)
  nullify(chemistry%kinmnrl_affinity_factor_sigma)
  if (associated(chemistry%kinmnrl_affinity_factor_beta)) deallocate(chemistry%kinmnrl_affinity_factor_beta)
  nullify(chemistry%kinmnrl_affinity_factor_beta)

end subroutine ChemistryDestroy

end module Chemistry_module
