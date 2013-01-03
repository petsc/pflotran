module Microbial_Aux_module
  
  use Database_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  PetscInt, parameter :: INHIBITION_THRESHOLD = 1
  PetscInt, parameter :: INHIBITION_THERMODYNAMIC = 2
  PetscInt, parameter :: INHIBITION_MONOD = 3
  PetscInt, parameter :: INHIBITION_INVERSE_MONOD = 4
  
  type, public :: biomass_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(biomass_species_type), pointer :: next    
  end type biomass_species_type
  
  type, public :: microbial_rxn_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: rate_constant
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn    
    type(monod_type), pointer :: monod
    type(inhibition_type), pointer :: inhibition
    type(biomass_type), pointer :: biomass
    type(microbial_rxn_type), pointer :: next
  end type microbial_rxn_type
  
  type, public :: monod_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: half_saturation_constant
    type(monod_type), pointer :: next
  end type monod_type

  type, public :: inhibition_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: inhibition_constant
    PetscReal :: concentration_threshold
    type(inhibition_type), pointer :: next
  end type inhibition_type

  type, public :: biomass_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: yield
  end type biomass_type
  
  type, public :: microbe_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
!TODO(geh): set up constraints for microbial life
!    character(len=MAXWORDLENGTH), pointer :: names(:)
!    PetscReal, pointer :: constraint_vol_frac(:)
!    PetscReal, pointer :: constraint_area(:)
!    PetscReal, pointer :: basis_vol_frac(:)
!    PetscReal, pointer :: basis_area(:)
!    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
!    PetscBool, pointer :: external_dataset(:)
  end type microbe_constraint_type
  
  type, public :: microbial_type

    PetscInt :: nrxn
    PetscInt :: nbiomass
    
    type(microbial_rxn_type), pointer :: microbial_rxn_list
    type(biomass_species_type), pointer :: biomass_list

    ! biomass species
    character(len=MAXWORDLENGTH), pointer :: biomass_names(:)
    PetscBool, pointer :: biomass_print(:)    

    ! microbial reactions
    PetscReal, pointer :: rate_constant(:)
    PetscReal, pointer :: stoich(:,:)
    PetscInt, pointer :: specid(:,:)
    PetscInt, pointer :: biomassid(:)
    PetscReal, pointer :: biomass_yield(:)
    PetscInt, pointer :: monodid(:,:)
    PetscInt, pointer :: inhibitionid(:,:)
    PetscInt, pointer :: monod_specid(:)
    PetscReal, pointer :: monod_K(:)
    PetscInt, pointer :: inhibition_specid(:)
    PetscReal, pointer :: inhibition_C(:)
    
  end type microbial_type

  public :: MicrobialCreate, &
            MicrobialRxnCreate, &
            MicrobialMonodCreate, &
            MicrobialInhibitionCreate, &
            MicrobialBiomassCreate, &
            MicrobialBiomassSpeciesCreate, &
            MicrobialGetMonodCount, &
            MicrobialGetInhibitionCount, &
            MicrobialGetBiomassCount, &
            MicrobialRxnDestroy, &
            MicrobialDestroy
             
contains

! ************************************************************************** !
!
! MicrobialCreate: Allocate and initialize microbial object
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
function MicrobialCreate()

  implicit none
  
  type(microbial_type), pointer :: MicrobialCreate
  
  type(microbial_type), pointer :: microbial

  allocate(microbial)  
    
  nullify(microbial%microbial_rxn_list)
  nullify(microbial%biomass_list)
    
  microbial%nrxn = 0
  microbial%nbiomass = 0

  nullify(microbial%rate_constant)
  nullify(microbial%stoich)
  nullify(microbial%specid)
  nullify(microbial%biomassid)
  nullify(microbial%biomass_yield)
  nullify(microbial%monodid)
  nullify(microbial%inhibitionid)
  nullify(microbial%monod_specid)
  nullify(microbial%monod_K)
  nullify(microbial%inhibition_specid)
  nullify(microbial%inhibition_C)
  
  MicrobialCreate => microbial
  
end function MicrobialCreate

! ************************************************************************** !
!
! MicrobialRxnCreate: Allocate and initialize a microbial object
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
function MicrobialRxnCreate()

  implicit none
  
  type(microbial_rxn_type), pointer :: MicrobialRxnCreate
  
  type(microbial_rxn_type), pointer :: microbial_rxn

  allocate(microbial_rxn)  
  microbial_rxn%id = 0
  microbial_rxn%itype = 0
  microbial_rxn%reaction = ''
  microbial_rxn%rate_constant = 0.d0
  microbial_rxn%print_me = PETSC_FALSE
  nullify(microbial_rxn%biomass)
  nullify(microbial_rxn%dbaserxn)
  nullify(microbial_rxn%monod)
  nullify(microbial_rxn%inhibition)
  nullify(microbial_rxn%next)
  
  MicrobialRxnCreate => microbial_rxn
  
end function MicrobialRxnCreate

! ************************************************************************** !
!
! MicrobialMonodCreate: Allocate and initialize a microbial monod object
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
function MicrobialMonodCreate()

  implicit none
  
  type(monod_type), pointer :: MicrobialMonodCreate
  
  type(monod_type), pointer :: monod

  allocate(monod)  
  monod%id = 0
  monod%species_name = ''
  monod%half_saturation_constant = 0.d0
  nullify(monod%next)
  
  MicrobialMonodCreate => monod
  
end function MicrobialMonodCreate

! ************************************************************************** !
!
! MicrobialInhibitionCreate: Allocate and initialize a microbial inhibition 
!                            object
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
function MicrobialInhibitionCreate()

  implicit none
  
  type(inhibition_type), pointer :: MicrobialInhibitionCreate
  
  type(inhibition_type), pointer :: inhibition
  
  allocate(inhibition)  
  inhibition%id = 0
  inhibition%itype = 0
  inhibition%species_name = ''
  inhibition%inhibition_constant = 0.d0
  inhibition%concentration_threshold = 0.d0
  nullify(inhibition%next)
  
  MicrobialInhibitionCreate => inhibition
  
end function MicrobialInhibitionCreate

! ************************************************************************** !
!
! MicrobialBiomassCreate: Allocate and initialize a microbial biomass object
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function MicrobialBiomassCreate()

  implicit none
  
  type(biomass_type), pointer :: MicrobialBiomassCreate
  
  type(biomass_type), pointer :: biomass

  allocate(biomass)  
  biomass%id = 0
  biomass%species_name = ''
  biomass%yield = 0.d0
  
  MicrobialBiomassCreate => biomass
  
end function MicrobialBiomassCreate

! ************************************************************************** !
!
! MicrobialBiomassSpeciesCreate: Allocate and initialize a biomass species 
!                                object
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function MicrobialBiomassSpeciesCreate()

  implicit none
  
  type(biomass_species_type), pointer :: MicrobialBiomassSpeciesCreate
  
  type(biomass_species_type), pointer :: species

  allocate(species)  
  species%id = 0
  species%name = ''
  species%molar_weight = 0.d0
  species%print_me = PETSC_FALSE
  nullify(species%next)

  MicrobialBiomassSpeciesCreate => species
  
end function MicrobialBiomassSpeciesCreate

! ************************************************************************** !
!
! MicrobialGetMonodCount: Counts number of monod expressions in
!                              microbial reaction
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
function MicrobialGetMonodCount(microbial_rxn)

  implicit none
  
  type(microbial_rxn_type) :: microbial_rxn
  
  PetscInt :: MicrobialGetMonodCount
  
  type(monod_type), pointer :: cur_monod
  PetscInt :: icount
  
  icount = 0
  cur_monod => microbial_rxn%monod
  do
    if (.not.associated(cur_monod)) exit
    icount = icount + 1
    cur_monod => cur_monod%next
  enddo
  
  MicrobialGetMonodCount = icount
  
end function MicrobialGetMonodCount

! ************************************************************************** !
!
! MicrobialGetInhibitionCount: Counts number of inhibiton expressions in
!                              microbial reaction
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
function MicrobialGetInhibitionCount(microbial_rxn)

  implicit none
  
  type(microbial_rxn_type) :: microbial_rxn
  
  PetscInt :: MicrobialGetInhibitionCount
  
  type(inhibition_type), pointer :: cur_inhibition
  PetscInt :: icount
  
  icount = 0
  cur_inhibition => microbial_rxn%inhibition
  do
    if (.not.associated(cur_inhibition)) exit
    icount = icount + 1
    cur_inhibition => cur_inhibition%next
  enddo
  
  MicrobialGetInhibitionCount = icount
  
end function MicrobialGetInhibitionCount

! ************************************************************************** !
!
! MicrobialGetBiomassCount: Returns the number of biomass species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
function MicrobialGetBiomassCount(microbial)

  implicit none
  
  PetscInt :: MicrobialGetBiomassCount
  type(microbial_type) :: microbial

  type(biomass_species_type), pointer :: biomass

  MicrobialGetBiomassCount = 0
  biomass => microbial%biomass_list
  do
    if (.not.associated(biomass)) exit
    MicrobialGetBiomassCount = MicrobialGetBiomassCount + 1
    biomass => biomass%next
  enddo

end function MicrobialGetBiomassCount

! ************************************************************************** !
!
! MicrobialRxnDestroy: Deallocates a microbial rxn object
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
subroutine MicrobialRxnDestroy(microbial)

  implicit none
    
  type(microbial_rxn_type), pointer :: microbial

  call DatabaseRxnDestroy(microbial%dbaserxn)
  call MicrobialMonodDestroy(microbial%monod)
  call MicrobialInhibitionDestroy(microbial%inhibition)
  call MicrobialBiomassDestroy(microbial%biomass)

  deallocate(microbial)  
  nullify(microbial)

end subroutine MicrobialRxnDestroy

! ************************************************************************** !
!
! MicrobialMonodDestroy: Deallocates a microbial monod object
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
recursive subroutine MicrobialMonodDestroy(monod)

  implicit none
    
  type(monod_type), pointer :: monod

  if (.not.associated(monod)) return
  
  call MicrobialMonodDestroy(monod%next)
  
  deallocate(monod)
  nullify(monod)
  
end subroutine MicrobialMonodDestroy

! ************************************************************************** !
!
! MicrobialInhibitionDestroy: Deallocates a microbial inhibition object
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
recursive subroutine MicrobialInhibitionDestroy(inhibition)

  implicit none
    
  type(inhibition_type), pointer :: inhibition

  if (.not. associated(inhibition)) return

  call MicrobialInhibitionDestroy(inhibition%next)
  
  deallocate(inhibition)
  nullify(inhibition)

end subroutine MicrobialInhibitionDestroy

! ************************************************************************** !
!
! MicrobialBiomassDestroy: Deallocates a microbial biomass object
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
subroutine MicrobialBiomassDestroy(biomass)

  implicit none
    
  type(biomass_type), pointer :: biomass

  if (.not.associated(biomass)) return
  
  deallocate(biomass)
  nullify(biomass)
  
end subroutine MicrobialBiomassDestroy

! ************************************************************************** !
!
! MicrobialBiomassSpeciesDestroy: Deallocates a biomass species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
subroutine MicrobialBiomassSpeciesDestroy(species)

  implicit none
    
  type(biomass_species_type), pointer :: species

  deallocate(species)  
  nullify(species)

end subroutine MicrobialBiomassSpeciesDestroy

! ************************************************************************** !
!
! MicrobialDestroy: Deallocates a microbial object
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine MicrobialDestroy(microbial)

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(microbial_type), pointer :: microbial
  
  type(microbial_rxn_type), pointer :: cur_microbial, prev_microbial
  type(biomass_species_type), pointer :: cur_biomass, prev_biomass

  if (.not.associated(microbial)) return
  
  ! microbial reactions
  cur_microbial => microbial%microbial_rxn_list
  do
    if (.not.associated(cur_microbial)) exit
    prev_microbial => cur_microbial
    cur_microbial => cur_microbial%next
    call MicrobialRxnDestroy(prev_microbial)
  enddo    
  nullify(microbial%microbial_rxn_list)
  
  ! biomass species
  cur_biomass => microbial%biomass_list
  do
    if (.not.associated(cur_biomass)) exit
    prev_biomass => cur_biomass
    cur_biomass => cur_biomass%next
    call MicrobialBiomassSpeciesDestroy(prev_biomass)
  enddo    
  nullify(microbial%biomass_list)
  
  call DeallocateArray(microbial%biomass_names)
  call DeallocateArray(microbial%biomass_print)
  
  call DeallocateArray(microbial%rate_constant)
  call DeallocateArray(microbial%stoich)
  call DeallocateArray(microbial%specid)
  call DeallocateArray(microbial%biomassid)
  call DeallocateArray(microbial%biomass_yield)
  call DeallocateArray(microbial%monodid)
  call DeallocateArray(microbial%inhibitionid)
  call DeallocateArray(microbial%monod_specid)
  call DeallocateArray(microbial%monod_K)
  call DeallocateArray(microbial%inhibition_specid)
  call DeallocateArray(microbial%inhibition_C)
  
  deallocate(microbial)
  nullify(microbial)

end subroutine MicrobialDestroy

end module Microbial_Aux_module
