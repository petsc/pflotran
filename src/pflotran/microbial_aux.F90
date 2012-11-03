module Microbial_Aux_module
  
  use Database_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  PetscInt, parameter :: INHIBITION_THRESHOLD = 1
  PetscInt, parameter :: INHIBITION_THERMODYNAMIC = 2
  PetscInt, parameter :: INHIBITION_MONOD = 3
  PetscInt, parameter :: INHIBITION_INVERSE_MONOD = 4

  type, public :: microbial_rxn_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: rate_constant
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn    
    type(monod_type), pointer :: monod
    type(inhibition_type), pointer :: inhibition
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
  
  type, public :: microbial_type

    PetscInt :: nrxn
    
    type(microbial_rxn_type), pointer :: microbial_rxn_list

    ! for saturation states
    PetscReal, pointer :: rate_constant(:)
    PetscReal, pointer :: stoich(:,:)
    PetscInt, pointer :: specid(:,:)
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
            MicrobialGetMonodCount, &
            MicrobialGetInhibitionCount, &
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
    
  microbial%nrxn = 0  

  nullify(microbial%rate_constant)
  nullify(microbial%stoich)
  nullify(microbial%specid)
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

  if (associated(monod%next)) then
    call MicrobialMonodDestroy(monod%next)
  endif
  
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

  if (associated(inhibition%next)) then
    call MicrobialInhibitionDestroy(inhibition%next)
  endif
  
  deallocate(inhibition)
  nullify(inhibition)

end subroutine MicrobialInhibitionDestroy

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
  
  call DeallocateArray(microbial%rate_constant)
  call DeallocateArray(microbial%stoich)
  call DeallocateArray(microbial%specid)
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
