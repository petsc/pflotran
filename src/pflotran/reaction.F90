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
  
  type, public :: ion_exchange_type
    character(len=MAXNAMELENGTH) :: reference_cation_name
    character(len=MAXNAMELENGTH) :: mnrl_name
    PetscInt :: ncation
    character(len=MAXNAMELENGTH), pointer :: cation_name(:)
    PetscInt, pointer :: cation_ids(:)
    PetscReal, pointer :: k(:)  ! selectivity coefficient
    PetscReal :: CEC
    type (ion_exchange_type), pointer :: next
  end type ion_exchange_type

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
    type(ion_exchange_type), pointer :: ion_exchange_list
    type(surface_complexation_rxn_type), pointer :: surface_complex_list
    ! compressed arrays for efficient computation
    ! primary aqueous complexes
    PetscReal, pointer :: primary_spec_molwt(:)
    PetscReal, pointer :: primary_spec_a0(:)
    PetscReal, pointer :: primary_spec_Z(:)
    ! aqueous complexes
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
      ! for kinetic reactions
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
            ChemistryRead

contains

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
          if (.not.associated(chemistry%primary_species_list)) &
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
      case default
        call printErrMsg(option,'CHEMISTRY keyword: '//trim(word)//' not recognized')
    end select
  enddo
 
end subroutine ChemistryRead

end module Chemistry_module
