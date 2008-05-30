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
