module Reaction_module

  use Reaction_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"
  
  public :: ReactionCreate, &
            ReactionRead, &
            CarbonateTestProblemCreate, &
            ReactionReadMineralRates, &
            ReactionReadSurfaceComplexes, &
            ReactionInitializeConstraint

contains

! ************************************************************************** !
!
! CarbonateTestProblemCreate: Creates a carbonate test problem for reactive
!                             transport
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine CarbonateTestProblemCreate(reaction,option)

  use Option_module
  
  type(reaction_type), pointer :: reaction
  type(option_type) :: option

  PetscInt :: icomp, irxn
  
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
  reaction%debyeA = 0.5114d0
  reaction%debyeB = 0.3288d0
  reaction%debyeBdot = 0.0410d0 
  reaction%num_dbase_temperatures = 1
  
  reaction%ncomp = option%ncomp
  allocate(reaction%primary_spec_Z(option%ncomp))
  reaction%primary_spec_Z(1) = 1.d0
  reaction%primary_spec_Z(2) = -1.d0
  reaction%primary_spec_Z(3) = 2.d0
  allocate(reaction%primary_spec_a0(option%ncomp))
  reaction%primary_spec_a0(1) = 9.d0
  reaction%primary_spec_a0(2) = 4.d0
  reaction%primary_spec_a0(3) = 6.d0
  
  reaction%neqcmplx = 5
  allocate(reaction%eqcmplxspecid(0:option%ncomp,reaction%neqcmplx))
  reaction%eqcmplxspecid = 0
  allocate(reaction%eqcmplxstoich(option%ncomp,reaction%neqcmplx))
  reaction%eqcmplxstoich = 0.d0
  allocate(reaction%eqcmplx_logK(reaction%neqcmplx))
  reaction%eqcmplx_logK = 0.d0
  allocate(reaction%eqcmplx_Z(reaction%neqcmplx))
  reaction%eqcmplx_Z = 0.d0
  allocate(reaction%eqcmplx_a0(reaction%neqcmplx))
  reaction%eqcmplx_a0 = 0.d0
  
  ! CO2(aq)
  irxn = 1
  reaction%eqcmplxspecid(0,irxn) = 2
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  reaction%eqcmplxstoich(1,irxn) = 1.d0 ! H+
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplx_logK(irxn) = -6.3447d0
  reaction%eqcmplx_Z(irxn) = 0.d0
  reaction%eqcmplx_a0(irxn) = 3.d0
  
  ! CO3-2
  irxn = 2
  reaction%eqcmplxspecid(0,irxn) = 2
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  reaction%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplx_logK(irxn) = 10.3288d0
  reaction%eqcmplx_z(irxn) = -2.d0
  reaction%eqcmplx_a0(irxn) = 4.5d0
  
  ! CaCO3(aq)
  irxn = 3
  reaction%eqcmplxspecid(0,irxn) = 3
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxspecid(2,irxn) = 2    ! HCO3-
  reaction%eqcmplxspecid(3,irxn) = 3    ! Ca+2
  reaction%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplxstoich(3,irxn) = 1.d0 ! Ca+2
  reaction%eqcmplx_logK(irxn) = 7.0017d0
  reaction%eqcmplx_z(irxn) = 0.d0
  reaction%eqcmplx_a0(irxn) = 3.d0

  ! CaHCO3+
  irxn = 4
  reaction%eqcmplxspecid(0,irxn) = 2
  reaction%eqcmplxspecid(1,irxn) = 2    ! HCO3-
  reaction%eqcmplxspecid(2,irxn) = 3    ! Ca+2
  reaction%eqcmplxstoich(1,irxn) = 1.d0 ! HCO3-
  reaction%eqcmplxstoich(2,irxn) = 1.d0 ! Ca+2
  reaction%eqcmplx_logK(irxn) = -1.0467d0
  reaction%eqcmplx_z(irxn) = 1.d0
  reaction%eqcmplx_a0(irxn) = 4.d0

  ! OH-
  irxn = 5
  reaction%eqcmplxspecid(0,irxn) = 1
  reaction%eqcmplxspecid(1,irxn) = 1    ! H+
  reaction%eqcmplxstoich(1,irxn) = -1.d0 ! H+
  reaction%eqcmplx_logK(irxn) = 13.9951
  reaction%eqcmplx_z(irxn) = -1.d0
  reaction%eqcmplx_a0(irxn) = 3.5d0
  
  reaction%nkinmnrl = 1
  allocate(reaction%kinmnrlspecid(0:option%ncomp,reaction%nkinmnrl))
  allocate(reaction%kinmnrlstoich(option%ncomp,reaction%nkinmnrl))
  allocate(reaction%kinmnrl_logK(reaction%nkinmnrl))
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
  reaction%kinmnrl_logK(irxn) = 1.8487d0
  reaction%kinmnrl_rate(1,irxn) = 1.d-6
  reaction%mnrl_molar_vol(irxn) = 36.9340d0/1.d6  ! based on 36.934 cm^3/mol
  reaction%kinmnrl_num_prefactors(irxn) = 0
  
end subroutine CarbonateTestProblemCreate

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
  
    call fiReadFlotranString(fid,string,ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY', ierr)
    call fiWordToUpper(word)   

    select case(trim(word))
    
      case('PRIMARY_SPECIES')
        nullify(prev_species)
        do
          call fiReadFlotranString(fid,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          species => AqueousSpeciesCreate()
          call fiReadWord(string,species%name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,PRIMARY_SPECIES', ierr)    
          if (.not.associated(reaction%primary_species_list)) then
            reaction%primary_species_list => species
            species%id = 1
          endif
          if (associated(prev_species)) then
            prev_species%next => species
            species%id = prev_species%id + 1
          endif
          prev_species => species
          nullify(species)
        enddo
      case('SECONDARY_SPECIES')
        nullify(prev_species)
        do
          call fiReadFlotranString(fid,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          species => AqueousSpeciesCreate()
          call fiReadWord(string,species%name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,PRIMARY_SPECIES', ierr)    
          if (.not.associated(reaction%secondary_species_list)) then
            reaction%secondary_species_list => species
            species%id = 1
          endif
          if (associated(prev_species)) then
            prev_species%next => species
            species%id = prev_species%id + 1            
          endif
          prev_species => species
          nullify(species)
        enddo
      case('GAS_SPECIES')
        nullify(prev_gas)
        do
          call fiReadFlotranString(fid,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          gas => GasSpeciesCreate()
          call fiReadWord(string,gas%name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,GAS_SPECIES', ierr)    
          if (.not.associated(reaction%gas_species_list)) then
            reaction%gas_species_list => gas
            gas%id = 1
          endif
          if (associated(prev_gas)) then
            prev_gas%next => gas
            gas%id = prev_gas%id + 1
          endif
          prev_gas => gas
          nullify(gas)
        enddo
      case('MINERALS')
        nullify(prev_mineral)
        do
          call fiReadFlotranString(fid,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
          mineral => MineralCreate()
          call fiReadWord(string,mineral%name,.true.,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,MINERAL', ierr)    
          if (.not.associated(reaction%mineral_list)) then
            reaction%mineral_list => mineral
            mineral%id = 1
          endif
          if (associated(prev_mineral)) then
            prev_mineral%next => mineral
            mineral%id = prev_mineral%id + 1
          endif
          prev_mineral => mineral
          nullify(mineral)
        enddo
      case('MINERAL_RATES')
        call fiSkipToEND(fid,option%myrank,word)
      case('SORPTION')
        do
          call fiReadFlotranString(fid,string,ierr)
          if (ierr /= 0) exit
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit

          call fiReadWord(string,word,.true.,ierr)
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SORPTION', ierr)
          call fiWordToUpper(word)   

          select case(trim(word))
            case('SURFACE_COMPLEXATION_RXN')
          !   call SurfaceComplexRXNRead
            case('ION_EXCHANGE_RXN')
          !   call IonExchangeRXNRead
            case('DISTRIBUTION_COEF')
          !   call DistributionCoefRead
          end select
        enddo
      case('DATABASE')
        call fiReadNChars(string,reaction%database_filename,MAXSTRINGLENGTH,.true.,ierr)  
        call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,DATABASE FILENAME', ierr)          
      case('END','/','.')
        exit
      case default
        call printErrMsg(option,'CHEMISTRY keyword: '//trim(word)//' not recognized')
    end select
  enddo
 
end subroutine ReactionRead

! ************************************************************************** !
!
! ReactionInitializeConstraint: Initializes constraints based on primary
!                               species in system
! author: Glenn Hammond
! date: 10/014/08
!
! ************************************************************************** !
subroutine ReactionInitializeConstraint(constraint_name, &
                                        aq_species_constraint, &
                                        mineral_constraint,option)
  use Option_module
  use Fileio_module
  
  implicit none
  
  character(len=MAXNAMELENGTH) :: constraint_name
  type(aq_species_constraint_type) :: aq_species_constraint
  type(mineral_constraint_type) :: mineral_constraint
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: found
  PetscInt :: icomp, jcomp
  
  ! aqueous species
  do icomp = 1, option%ncomp
    found = PETSC_FALSE
    do jcomp = 1, option%ncomp
      if (fiStringCompare(aq_species_constraint%names(icomp), &
                          option%comp_names(jcomp), &
                          MAXNAMELENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      string = 'Species ' // trim(aq_species_constraint%names(icomp)) // &
               'from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among primary species.'
      call printErrMsg(option,string)
    else
      aq_species_constraint%basis_conc(jcomp) = aq_species_constraint%conc(icomp)
    endif
  enddo
  
  do icomp = 1, option%nmnrl
    found = PETSC_FALSE
    do jcomp = 1, option%nmnrl
      if (fiStringCompare(mineral_constraint%names(icomp), &
                          option%mnrl_names(jcomp), &
                          MAXNAMELENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      string = 'Mineral ' // trim(mineral_constraint%names(icomp)) // &
               'from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among primary species.'
      call printErrMsg(option,string)
    else
      mineral_constraint%basis_conc(jcomp) = mineral_constraint%conc(icomp)
    endif  
  enddo

end subroutine ReactionInitializeConstraint

! ************************************************************************** !
!
! ReactionRead: Reads chemical species
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
subroutine ReactionReadMineralRates(reaction,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXNAMELENGTH) :: name
  
  type(mineral_type), pointer :: cur_mineral
  PetscErrorCode :: ierr

  cur_mineral => reaction%mineral_list
  do 
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -1*abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)
    if (ierr /= 0) exit

    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadNChars(string,name,MAXNAMELENGTH,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY', ierr)
    
    cur_mineral => reaction%mineral_list
    do 
      if (.not.associated(cur_mineral)) exit
      if (fiStringCompare(cur_mineral%name,name,MAXNAMELENGTH)) then
        if (.not.associated(cur_mineral%tstrxn)) then
          cur_mineral%tstrxn => TransitionStateTheoryRxnCreate()
        endif
        call fiReadDouble(string,cur_mineral%tstrxn%rate,ierr)
        cur_mineral%id = abs(cur_mineral%id)
        exit
      endif
      cur_mineral => cur_mineral%next
    enddo
    
  enddo
 
  cur_mineral => reaction%mineral_list
  do 
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%id < 0) then
      string = 'No rate provided in input file for mineral: ' // &
               trim(cur_mineral%name) // '.'
      call printErrMsg(option,string)
    endif
    cur_mineral => cur_mineral%next
  enddo
  
end subroutine ReactionReadMineralRates
! ************************************************************************** !
!
! ReactionRead: Reads chemical species
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
subroutine ReactionReadSurfaceComplexes(reaction,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXNAMELENGTH) :: name
  
  type(surface_complexation_rxn_type), pointer :: cur_srfcmplx
  PetscErrorCode :: ierr

  cur_srfcmplx => reaction%surface_complex_list
  do 
    if (.not.associated(cur_srfcmplx)) exit
    cur_srfcmplx%id = -1*abs(cur_srfcmplx%id)
    cur_srfcmplx => cur_srfcmplx%next
  enddo

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)
    if (ierr /= 0) exit

    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadNChars(string,name,MAXNAMELENGTH,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY', ierr)
    
    cur_srfcmplx => reaction%surface_complex_list
    do 
      if (.not.associated(cur_srfcmplx)) exit
      if (fiStringCompare(cur_srfcmplx%name,name,MAXNAMELENGTH)) then
!       if (.not.associated(cur_srfcmplx%tstrxn)) then
!         cur_srfcmplx%tstrxn => TransitionStateTheoryRxnCreate()
!       endif
!       call fiReadDouble(string,cur_srfcmplx%tstrxn%rate,ierr)
!       cur_srfcmplx%id = abs(cur_srfcmplx%id)
        exit
      endif
      cur_srfcmplx => cur_srfcmplx%next
    enddo
    
  enddo
 
  cur_srfcmplx => reaction%surface_complex_list
  do 
    if (.not.associated(cur_srfcmplx)) exit
    if (cur_srfcmplx%id < 0) then
      string = 'No surface complex site density provided for mineral: ' // &
               trim(cur_srfcmplx%name) // '.'
      call printErrMsg(option,string)
    endif
    cur_srfcmplx => cur_srfcmplx%next
  enddo
  
end subroutine ReactionReadSurfaceComplexes

end module Reaction_module
