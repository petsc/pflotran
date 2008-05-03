module Reaction_Aux_module

  ! this module cannot depend on any other modules beside Option_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: aq_species_type
    character(len=MAXNAMELENGTH) :: spec_name
    PetscReal :: a0
    PetscReal :: molar_weight
    PetscReal :: valence
    type(equilibrium_rxn_type), pointer :: eqrxn
  end type aq_species_type

  type, public :: gas_species_type
    character(len=MAXNAMELENGTH) :: spec_name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    type(equilibrium_rxn_type), pointer :: eqrxn
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
    PetscReal, allocatable :: primary_spec_molwt(:)
    PetscReal, allocatable :: primary_spec_a0(:)
    PetscReal, allocatable :: primary_spec_Z(:)
    ! aqueous complexes
    PetscInt, allocatable :: eqcmplxspecid(:,:)   ! (0:ncomp in rxn)
    PetscReal, allocatable :: eqcmplxstoich(:,:)
    PetscReal, allocatable :: eqcmplx_a0(:)
    PetscReal, allocatable :: eqcmplx_K(:)
    ! ionx exchange reactions
    PetscInt, allocatable :: eqionx_ncation(:)
    PetscReal, allocatable :: eqionx_CEC(:)
    PetscReal, allocatable :: eqionx_k(:,:)
    PetscInt, allocatable :: eqionx_cationid(:)
    PetscInt, allocatable :: eqionx_rxn_offset(:)
    PetscInt, allocatable :: kinionx_ncation(:)
    PetscReal, allocatable :: kinionx_CEC(:)
    PetscReal, allocatable :: kinionx_k(:,:)
    PetscInt, allocatable :: kinionx_cationid(:)
    PetscInt, allocatable :: kinionx_rxn_offset(:)
    ! surface complexation reactions
    PetscInt, allocatable :: eqsurfcmplxspecid(:,:)
    PetscReal, allocatable :: eqsurfcmplxstoich(:,:)
    PetscReal, allocatable :: eqsurfcmplx_freesite_stoich(:,:)
    PetscReal, allocatable :: eqsurfcmplx_K(:)
    PetscReal, allocatable :: eqsurfcmplx_Z(:)  ! valence
    PetscInt, allocatable :: kinsurfcmplxspecid(:,:)
    PetscReal, allocatable :: kinsurfcmplxstoich(:,:)
    PetscReal, allocatable :: kinsurfcmplx_freesite_stoich(:,:)
    PetscReal, allocatable :: kinsurfcmplx_K(:)
    PetscReal, allocatable :: kinsurfcmplx_Z(:)  ! valence
    ! mineral reactions
      ! for saturation states
    PetscInt, allocatable :: mnrlspecid(:,:)
    PetscReal, allocatable :: mnrlstoich(:,:)
    PetscReal, allocatable :: mnrl_K(:)
      ! for kinetic reactions
    PetscInt, allocatable :: kinmnrlspecid(:,:)
    PetscReal, allocatable :: kinmnrlstoich(:,:)
    PetscReal, allocatable :: kinmnrl_K(:)
    PetscReal, allocatable :: kinmnrl_rate(:)
    PetscInt, allocatable :: kinmnrl_pri_prefactor_id(:,:)
    PetscInt, allocatable :: kinmnrl_sec_prefactor_id(:,:)
    PetscInt, allocatable :: kinmnrl_pri_prefactor_stoich(:,:)
    PetscInt, allocatable :: kinmnrl_sec_prefactor_stoich(:,:)
    PetscReal, allocatable :: kinmnrl_affinity_factor_sigma(:)
    PetscReal, allocatable :: kinmnrl_affinity_factor_beta(:)
  end type reaction_type
  
contains


! ************************************************************************** !
!
! ReactionAuxCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 04/23/08
!
! ************************************************************************** !
function ReactionAuxCreate()

  use Option_module

  implicit none
  
  type(reaction_type), pointer :: ReactionAuxCreate
  
  type(reaction_type), pointer :: aux

  allocate(aux)  

  ReactionAuxCreate => aux
  
end function ReactionAuxCreate

! ************************************************************************** !
!
! ReactionAuxVarInit: Initialize auxilliary object
! author: Glenn Hammond
! date: 04/23/08
!
! ************************************************************************** !
subroutine ReactionAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(reaction_type) :: aux_var
  type(option_type) :: option  
  
  
end subroutine ReactionAuxVarInit

! ************************************************************************** !
!
! ReactionAuxVarCompute: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 04/23/08
!
! ************************************************************************** !
subroutine ReactionAuxVarCompute(x,aux_var,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  PetscReal :: x(option%ncomp)  
  type(reaction_type) :: aux_var
  
end subroutine ReactionAuxVarCompute

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a reaction auxilliary object
! author: Glenn Hammond
! date: 04/23/08
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(reaction_type) :: aux_var

end subroutine AuxVarDestroy

! ************************************************************************** !
!
! ReactionAuxDestroy: Deallocates a reaction auxilliary object
! author: Glenn Hammond
! date: 04/23/08
!
! ************************************************************************** !
subroutine ReactionAuxDestroy(aux)

  implicit none

  type(reaction_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return 
      
end subroutine ReactionAuxDestroy

end module Reaction_Aux_module
