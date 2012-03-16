module Database_Aux_module

  implicit none
  
  private
  
#include "definitions.h"

  type, public :: database_rxn_type
    PetscInt :: nspec
    character(len=MAXWORDLENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: spec_ids(:)
    PetscReal, pointer :: logK(:)
    PetscReal, pointer :: logKCoeff_hpt(:)
  end type database_rxn_type

  public :: BasisAlignSpeciesInRxn, &
            BasisSubSpeciesInGasOrSecRxn, &
            BasisSubSpeciesinMineralRxn, &
            DatabaseRxnCreate, &
            DatabaseRxnDestroy
            
contains

! ************************************************************************** !
!
! DatabaseRxnCreate: Allocate and initialize an equilibrium reaction
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
function DatabaseRxnCreate()

  implicit none
    
  type(database_rxn_type), pointer :: DatabaseRxnCreate

  type(database_rxn_type), pointer :: dbaserxn

  allocate(dbaserxn)
  dbaserxn%nspec = 0
  nullify(dbaserxn%spec_name)
  nullify(dbaserxn%stoich)
  nullify(dbaserxn%spec_ids)
  nullify(dbaserxn%logK)
  
  DatabaseRxnCreate => dbaserxn
  
end function DatabaseRxnCreate

! ************************************************************************** !
!
! DatabaseRxnDestroy: Deallocates a database reaction
! author: Glenn Hammond
! date: 05/29/08
!
! ************************************************************************** !
subroutine DatabaseRxnDestroy(dbaserxn)

  implicit none
    
  type(database_rxn_type), pointer :: dbaserxn

  if (.not.associated(dbaserxn)) return
  
  if (associated(dbaserxn%spec_name)) deallocate(dbaserxn%spec_name)
  nullify(dbaserxn%spec_name)
  if (associated(dbaserxn%spec_ids)) deallocate(dbaserxn%spec_ids)
  nullify(dbaserxn%spec_ids)
  if (associated(dbaserxn%stoich)) deallocate(dbaserxn%stoich)
  nullify(dbaserxn%stoich)
  if (associated(dbaserxn%logK)) deallocate(dbaserxn%logK)
  nullify(dbaserxn%logK)

  deallocate(dbaserxn)  
  nullify(dbaserxn)

end subroutine DatabaseRxnDestroy

! ************************************************************************** !
!
! BasisAlignSpeciesInRxn: Aligns the ordering of species in reaction with
!                         the current basis
! author: Glenn Hammond
! date: 10/07/08
!
! ************************************************************************** !
subroutine BasisAlignSpeciesInRxn(num_basis_species,basis_names, &
                                  num_rxn_species,rxn_species_names, &
                                  rxn_stoich,rxn_species_ids,species_name, &
                                  option)

  use Option_module
  use String_module
  
  implicit none
  
  PetscInt :: num_basis_species
  character(len=MAXWORDLENGTH) :: basis_names(num_basis_species), species_name
  PetscInt :: num_rxn_species
  character(len=MAXWORDLENGTH) :: rxn_species_names(num_rxn_species)
  PetscReal :: rxn_stoich(num_rxn_species)
  PetscInt :: rxn_species_ids(num_rxn_species)
  type(option_type) :: option

  PetscInt :: i_rxn_species
  PetscInt :: i_basis_species
  PetscReal :: stoich_new(num_basis_species)
  PetscBool :: found
  
  stoich_new = 0.d0
  do i_rxn_species = 1, num_rxn_species
    found = PETSC_FALSE
    do i_basis_species = 1, num_basis_species
      if (StringCompare(rxn_species_names(i_rxn_species), &
                            basis_names(i_basis_species), &
                            MAXWORDLENGTH)) then
        stoich_new(i_basis_species) = rxn_stoich(i_rxn_species)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = trim(rxn_species_names(i_rxn_species)) // &
               ' not found in basis (BasisAlignSpeciesInRxn) for species ' // &
               trim(species_name)
      call printErrMsg(option)
    endif
  enddo
  
  ! zero everthing out
  rxn_species_names = ''
  rxn_stoich = 0.d0
  rxn_species_ids = 0

  ! fill in
  i_rxn_species = 0
  do i_basis_species = 1, num_basis_species
    if (dabs(stoich_new(i_basis_species)) > 1.d-40) then
      i_rxn_species = i_rxn_species + 1
      rxn_species_names(i_rxn_species) = basis_names(i_basis_species)
      rxn_stoich(i_rxn_species) = stoich_new(i_basis_species)
      rxn_species_ids(i_rxn_species) = i_basis_species
    endif
  enddo
  
  if (i_rxn_species /= num_rxn_species) then
    write(option%io_buffer,*) &
                   'Number of reaction species does not match original:', &
                    i_rxn_species, num_rxn_species
    call printErrMsg(option)
  endif

end subroutine BasisAlignSpeciesInRxn 
     
! ************************************************************************** !
!
! BasisSubSpeciesInGasOrSecRxn: Swaps out a chemical species in a chemical  
!                                reaction, replacing it with the species in a 
!                                secondary reaction (swaps 1 into 2)
! author: Glenn Hammond
! date: 10/06/08
!
! ************************************************************************** !
subroutine BasisSubSpeciesInGasOrSecRxn(name1,dbaserxn1,dbaserxn2,scale)

  use String_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: name1
  type(database_rxn_type) :: dbaserxn1
  type(database_rxn_type) :: dbaserxn2
  PetscReal :: scale
  
  PetscInt :: i, j, tempcount, prevcount
  character(len=MAXWORDLENGTH) :: tempnames(20)
  PetscReal :: tempstoich(20)
  PetscBool :: found

  tempnames = ''
  tempstoich = 0.d0

  ! load species in reaction other than species 1 into new arrays
  scale = 1.d0
  tempcount = 0
  do i=1,dbaserxn2%nspec
    if (.not.StringCompare(name1, &
                             dbaserxn2%spec_name(i), &
                             MAXWORDLENGTH)) then
      tempcount = tempcount + 1
      tempnames(tempcount) = dbaserxn2%spec_name(i)
      tempstoich(tempcount) = dbaserxn2%stoich(i)
    else
      scale = dbaserxn2%stoich(i)
    endif
  enddo
  
  ! search for duplicate species and add stoichs or add new species
  ! if not duplicated
  do j=1,dbaserxn1%nspec
    found = PETSC_FALSE
    do i=1,tempcount
      if (StringCompare(tempnames(i), &
                          dbaserxn1%spec_name(j), &
                          MAXWORDLENGTH)) then
        tempstoich(i) = tempstoich(i) + scale*dbaserxn1%stoich(j)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      tempcount = tempcount + 1
      tempnames(tempcount) = dbaserxn1%spec_name(j)
      tempstoich(tempcount) = scale*dbaserxn1%stoich(j)
    endif
  enddo
  
  ! deallocate arrays
  deallocate(dbaserxn2%spec_name)
  deallocate(dbaserxn2%stoich)
  
  ! check for zero stoichiometries due to cancelation
  prevcount = tempcount
  tempcount = 0
  do i=1,prevcount
    if (dabs(tempstoich(i)) > 1.d-10) then
      tempcount = tempcount + 1
      tempnames(tempcount) = tempnames(i)
      tempstoich(tempcount) = tempstoich(i)
    endif
  enddo
  
  tempnames(tempcount+1:) = ''
  tempstoich(tempcount+1:) = 0.d0
  
  ! reallocate
  allocate(dbaserxn2%spec_name(tempcount))
  allocate(dbaserxn2%stoich(tempcount))
  
  ! fill arrays in dbaserxn
  dbaserxn2%nspec = tempcount
  do i=1,tempcount
    dbaserxn2%spec_name(i) = tempnames(i)
    dbaserxn2%stoich(i) = tempstoich(i)
  enddo
  
  !dbaserxn2%logK = dbaserxn2%logK + scale*dbaserxn1%logK

end subroutine BasisSubSpeciesInGasOrSecRxn

! ************************************************************************** !
!
! BasisSubSpeciesInMineralRxn: Swaps out a chemical species in a chemical  
!                                reaction, replacing it with the species in a 
!                                secondary reaction (swaps 1 into 2)
! author: Glenn Hammond
! date: 10/06/08
!
! ************************************************************************** !
subroutine BasisSubSpeciesInMineralRxn(name,sec_dbaserxn,mnrl_dbaserxn,scale)

  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(database_rxn_type) :: sec_dbaserxn
  type(database_rxn_type) :: mnrl_dbaserxn
  PetscReal :: scale
  
  PetscInt :: i, j, tempcount, prevcount
  character(len=MAXWORDLENGTH) :: tempnames(20)
  PetscReal :: tempstoich(20)
  PetscBool :: found

  tempnames = ''
  tempstoich = 0.d0
  
  ! load species in reaction other than species 1 into new arrays
  scale = 1.d0
  tempcount = 0
  do i=1,mnrl_dbaserxn%nspec
    if (.not.StringCompare(name, &
                           mnrl_dbaserxn%spec_name(i), &
                           MAXWORDLENGTH)) then
      tempcount = tempcount + 1
      tempnames(tempcount) = mnrl_dbaserxn%spec_name(i)
      tempstoich(tempcount) = mnrl_dbaserxn%stoich(i)
    else
      scale = mnrl_dbaserxn%stoich(i)
    endif
  enddo
  
  ! search for duplicate species and add stoichs or add new species
  ! if not duplicated
  do j=1,sec_dbaserxn%nspec
    found = PETSC_FALSE
    do i=1,mnrl_dbaserxn%nspec
      if (StringCompare(tempnames(i), &
                        sec_dbaserxn%spec_name(j), &
                        MAXWORDLENGTH)) then
        tempstoich(i) = tempstoich(i) + scale*sec_dbaserxn%stoich(j)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      tempcount = tempcount + 1
      tempnames(tempcount) = sec_dbaserxn%spec_name(j)
      tempstoich(tempcount) = scale*sec_dbaserxn%stoich(j)
    endif
  enddo
  
  ! deallocate arrays
  deallocate(mnrl_dbaserxn%spec_name)
  deallocate(mnrl_dbaserxn%stoich)
  
  ! check for zero stoichiometries due to cancelation
  prevcount = tempcount
  tempcount = 0
  do i=1,prevcount
    if (dabs(tempstoich(i)) > 1.d-10) then
      tempcount = tempcount + 1
      tempnames(tempcount) = tempnames(i)
      tempstoich(tempcount) = tempstoich(i)
    endif
  enddo

  tempnames(tempcount+1:) = ''
  tempstoich(tempcount+1:) = 0.d0
    
  ! reallocate
  allocate(mnrl_dbaserxn%spec_name(tempcount))
  allocate(mnrl_dbaserxn%stoich(tempcount))
  
  ! fill arrays in dbaserxn
  mnrl_dbaserxn%nspec = tempcount
  do i=1,tempcount
    mnrl_dbaserxn%spec_name(i) = tempnames(i)
    mnrl_dbaserxn%stoich(i) = tempstoich(i)
  enddo
  
   !mnrl_dbaserxn%logK = mnrl_dbaserxn%logK + scale*sec_dbaserxn%logK

end subroutine BasisSubSpeciesInMineralRxn

end module Database_Aux_module
