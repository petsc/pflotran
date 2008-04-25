module Reaction_Aux_module

  ! this module cannot depend on any other modules beside Option_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: equilibrium_rxn_type
    PetscInt :: ncomp
    character(len=MAXNAMELENGTH), pointer :: comp_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: comp_ids(:)
    PetscReal, pointer :: a0
    PetscReal, pointer :: z
    PetscReal, pointer :: logK
    type(equilibrium_rxn_type), pointer :: next
  end type equilibrium_rxn_type

  type, public :: kinetic_rxn_type
    PetscInt :: ncomp
    character(len=MAXNAMELENGTH), pointer :: comp_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: comp_ids(:)
    PetscReal, pointer :: logK
    PetscReal :: rate_forward
    PetscReal :: rate_reverse
    type(kinetic_rxn_type), pointer :: next
  end type kinetic_rxn_type

  type, public :: mineral_rxn_type
    character(len=MAXNAMELENGTH) :: mnrl_name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    PetscInt :: ncomp
    character(len=MAXNAMELENGTH), pointer :: comp_name(:)
    PetscInt, pointer :: comp_ids(:)
    PetscReal, pointer :: stoich(:)
    PetscReal :: logK
    PetscReal :: rate
    PetscReal :: area
    type(mineral_rxn_type), pointer :: next
  end type mineral_rxn_type
  
  type, public :: ion_exchange_rxn_type
    character(len=MAXNAMELENGTH) :: cation_name
    PetscInt :: ncomp
    character(len=MAXNAMELENGTH), pointer :: comp_name(:)
    PetscInt, pointer :: comp_ids(:)
    PetscReal, pointer :: z(:)
    PetscReal, pointer :: k(:)
    PetscReal, pointer :: X_(:)
    PetscReal, pointer :: prev_X_(:)
    PetscReal :: CEC
    type (ion_exchange_rxn_type), pointer :: next
  end type ion_exchange_rxn_type
 
  type, public :: surface_complexation_rxn_type
    character(len=MAXNAMELENGTH) :: surfcmplx_name
    PetscInt :: ncomp
    character(len=MAXNAMELENGTH), pointer :: comp_name(:)
    PetscInt, pointer :: comp_ids(:)
    PetscReal, pointer :: stoich(:)
    PetscReal :: charge
    type (surface_complexation_rxn_type), pointer :: next
  end type surface_complexation_rxn_type
  
  type, public :: reaction_type
    type(equilibrium_rxn_type), pointer :: equilibrium_list
    type(kinetic_rxn_type), pointer :: kinetic_list
    type(mineral_rxn_type), pointer :: mineral_list
    type(ion_exchange_rxn_type), pointer :: ion_exchange_list
    type(surface_complexation_rxn_type), pointer :: surface_complex_list
    ! compressed arrays for efficient computation
    ! equilibrium complexation
    PetscInt, allocatable :: eqspeccmpid(:,:)
    PetscReal, allocatable :: steqspec(:,:)
    PetscReal, allocatable :: Keqspec(:)
    ! kinetic reversible
    PetscInt, allocatable :: kincmpid(:,:)
    PetscReal, allocatable :: stkin(:,:)
    PetscReal, allocatable :: Kkin(:)
    PetscReal, allocatable :: ratefor(:)
    PetscReal, allocatable :: ratebac(:)
    ! kinetic mineral precipitation/dissolution
    PetscInt, allocatable :: kinmnrlcmpid(:,:)
    PetscReal, allocatable :: stkinmnrl(:,:)
    PetscReal, allocatable :: Kkinmnrl(:)
    PetscReal, allocatable :: kinmnrlrate(:)
    PetscReal, allocatable :: kinmnrlarea(:)
    ! ion exchange
    ! surface complexation
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
