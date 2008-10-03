module Database_module

  implicit none
  
  private
  
#include "definitions.h"

  public :: DatabaseRead, BasisInit
  
contains

! ************************************************************************** !
!
! DatabaseRead: Collects parameters from geochemical database
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
subroutine DatabaseRead(reaction,option)

  use Reaction_module
  use Option_module
  use Fileio_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_type), pointer :: cur_mineral
  type(surface_complexation_rxn_type), pointer :: cur_surfcplx
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXNAMELENGTH) :: name
  character(len=MAXNAMELENGTH) :: null_name
  
  PetscTruth :: flag, found
  PetscInt :: ispec, itemp
  PetscInt, parameter :: dbase_id = 86
  PetscInt :: iostat
  PetscInt :: num_nulls
  PetscErrorCode :: ierr
  
  ! negate ids for use as flags
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    cur_aq_spec%id = -cur_aq_spec%id
    cur_aq_spec => cur_aq_spec%next
  enddo
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    cur_aq_spec%id = -cur_aq_spec%id
    cur_aq_spec => cur_aq_spec%next
  enddo  
  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    cur_gas_spec%id = -cur_gas_spec%id
    cur_gas_spec => cur_gas_spec%next
  enddo  
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -cur_mineral%id
    cur_mineral => cur_mineral%next
  enddo
  cur_surfcplx => reaction%surface_complex_list
  do
    if (.not.associated(cur_surfcplx)) exit
    cur_surfcplx%id = -cur_surfcplx%id
    cur_surfcplx => cur_surfcplx%next
  enddo  
  
  reaction%database_filename = 'hanford.dat'
  open(unit=dbase_id,file=reaction%database_filename,status='old',iostat=ierr)
  if (ierr /= 0) then
    string = 'DATABASE File: ' // reaction%database_filename // ' not found.'
    call printErrMsg(option,string)  
  endif

  ! read temperatures
  call fiReadDBaseString(dbase_id,string,ierr)
  ! remove comment
  call fiReadDBaseName(dbase_id,string,name,.true.,ierr)
  call fiReadDBaseInt(dbase_id,string,reaction%num_dbase_temperatures,ierr)
  call fiErrorMsg(option%myrank,'Number of database temperatures','DATABASE',ierr)  
  allocate(reaction%dbase_temperatures(reaction%num_dbase_temperatures))
  reaction%dbase_temperatures = 0.d0  
  do itemp = 1, reaction%num_dbase_temperatures
    call fiReadDBaseDouble(dbase_id,string,reaction%dbase_temperatures(itemp),ierr)
    call fiErrorMsg(option%myrank,'Database temperatures','DATABASE',ierr)            
  enddo

  num_nulls = 0
  null_name = 'null'
  do ! loop over every entry in the database
    call fiReadDBaseString(dbase_id,string,ierr)
    call fiReadStringErrorMsg(option%myrank,'DATABASE',ierr)

    call fiReadDBaseName(dbase_id,string,name,.true.,ierr)
    ! 'null's mark the end of a section in the database.  We count these 
    ! to determine which species we are reading.
    ! --
    ! primary species
    ! null
    ! aq complexes
    ! null
    ! gases
    ! null
    ! minerals
    ! null
    ! surface complexes
    ! null
    ! --
    
    if (fiStringCompare(name,null_name,MAXNAMELENGTH)) then
      num_nulls = num_nulls + 1
      if (num_nulls >= 5) exit
      cycle
    endif
    
    select case(num_nulls)
      case(0,1) ! primary and secondary aq species
        cur_aq_spec => reaction%primary_species_list
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_aq_spec)) exit
          if (fiStringCompare(name,cur_aq_spec%name,MAXNAMELENGTH)) then
            found = PETSC_TRUE
            ! change negative id to positive, indicating it was found in database
            cur_aq_spec%id = abs(cur_aq_spec%id)
            exit
          endif
          cur_aq_spec => cur_aq_spec%next
        enddo
        if (.not.found) cur_aq_spec => reaction%secondary_species_list
        do
          if (found .or. .not.associated(cur_aq_spec)) exit
          if (fiStringCompare(name,cur_aq_spec%name,MAXNAMELENGTH)) then
            found = PETSC_TRUE          
          ! change negative id to positive, indicating it was found in database
            cur_aq_spec%id = abs(cur_aq_spec%id)
            exit
          endif
          cur_aq_spec => cur_aq_spec%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        if (num_nulls > 0) then ! secondary species in database
          ! create aqueous equilibrium reaction
          cur_aq_spec%eqrxn => EquilibriumRxnCreate()
          ! read the number of primary species in secondary rxn
          call fiReadDBaseInt(dbase_id,string,cur_aq_spec%eqrxn%nspec,ierr)
          call fiErrorMsg(option%myrank,'Number of species in aqueous complex', &
                          'DATABASE',ierr)  
          ! allocate arrays for rxn
          allocate(cur_aq_spec%eqrxn%spec_name(cur_aq_spec%eqrxn%nspec))
          cur_aq_spec%eqrxn%spec_name = ''
          allocate(cur_aq_spec%eqrxn%stoich(cur_aq_spec%eqrxn%nspec))
          cur_aq_spec%eqrxn%stoich = 0.d0
          allocate(cur_aq_spec%eqrxn%logK(reaction%num_dbase_temperatures))
          cur_aq_spec%eqrxn%logK = 0.d0
          ! read in species and stoichiometries
          do ispec = 1, cur_aq_spec%eqrxn%nspec
            call fiReadDBaseDouble(dbase_id,string,cur_aq_spec%eqrxn%stoich(ispec),ierr)
            call fiErrorMsg(option%myrank,'EQRXN species stoichiometry','DATABASE',ierr)            
            call fiReadDBaseName(dbase_id,string,cur_aq_spec%eqrxn%spec_name(ispec),.true.,ierr)
            call fiErrorMsg(option%myrank,'EQRXN species name','DATABASE',ierr)            
          enddo
          do itemp = 1, reaction%num_dbase_temperatures
            call fiReadDBaseDouble(dbase_id,string,cur_aq_spec%eqrxn%logK(itemp),ierr)
            call fiErrorMsg(option%myrank,'EQRXN logKs','DATABASE',ierr)            
          enddo
        endif 
        ! read the Debye-Huckel ion size parameter (a0)
        call fiReadDBaseDouble(dbase_id,string,cur_aq_spec%a0,ierr)
        call fiErrorMsg(option%myrank,'AQ Species a0','DATABASE',ierr)            
        ! read the valence
        call fiReadDBaseDouble(dbase_id,string,cur_aq_spec%Z,ierr)
        call fiErrorMsg(option%myrank,'AQ Species Z','DATABASE',ierr)            
        ! read the molar weight
        call fiReadDBaseDouble(dbase_id,string,cur_aq_spec%molar_weight,ierr)
        call fiErrorMsg(option%myrank,'AQ Species molar weight','DATABASE',ierr)
        
                    
      case(2) ! gas species
        cur_gas_spec => reaction%gas_species_list
        if (.not.associated(cur_gas_spec)) cycle
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_gas_spec)) exit
          if (fiStringCompare(name,cur_gas_spec%name,MAXNAMELENGTH)) then
            found = PETSC_TRUE          
          ! change negative id to positive, indicating it was found in database
            cur_gas_spec%id = abs(cur_gas_spec%id)
            exit
          endif
          cur_gas_spec => cur_gas_spec%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        ! read the molar volume
        call fiReadDBaseDouble(dbase_id,string,cur_gas_spec%molar_volume,ierr)
        call fiErrorMsg(option%myrank,'GAS molar volume','DATABASE',ierr)            
        ! create aqueous equilibrium reaction
        cur_gas_spec%eqrxn => EquilibriumRxnCreate()
        ! read the number of aqueous species in secondary rxn
        call fiReadDBaseInt(dbase_id,string,cur_gas_spec%eqrxn%nspec,ierr)
        call fiErrorMsg(option%myrank,'Number of species in gas reaction', &
                        'DATABASE',ierr)  
        ! allocate arrays for rxn
        allocate(cur_gas_spec%eqrxn%spec_name(cur_gas_spec%eqrxn%nspec))
        cur_gas_spec%eqrxn%spec_name = ''
        allocate(cur_gas_spec%eqrxn%stoich(cur_gas_spec%eqrxn%nspec))
        cur_gas_spec%eqrxn%stoich = 0.d0
        allocate(cur_gas_spec%eqrxn%logK(reaction%num_dbase_temperatures))
        cur_gas_spec%eqrxn%logK = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_gas_spec%eqrxn%nspec
          call fiReadDBaseDouble(dbase_id,string,cur_gas_spec%eqrxn%stoich(ispec),ierr)
          call fiErrorMsg(option%myrank,'GAS species stoichiometry','DATABASE',ierr)            
          call fiReadDBaseName(dbase_id,string,cur_gas_spec%eqrxn%spec_name(ispec),.true.,ierr)
          call fiErrorMsg(option%myrank,'GAS species name','DATABASE',ierr)            
        enddo
        do itemp = 1, reaction%num_dbase_temperatures
          call fiReadDBaseDouble(dbase_id,string,cur_gas_spec%eqrxn%logK(itemp),ierr)
          call fiErrorMsg(option%myrank,'GAS logKs','DATABASE',ierr)            
        enddo
        ! read the molar weight
        call fiReadDBaseDouble(dbase_id,string,cur_gas_spec%molar_weight,ierr)
        call fiErrorMsg(option%myrank,'GAS molar weight','DATABASE',ierr)     
        
               
      case(3) ! minerals
        cur_mineral => reaction%mineral_list
        if (.not.associated(cur_mineral)) cycle
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_mineral)) exit
          if (fiStringCompare(name,cur_mineral%name,MAXNAMELENGTH)) then
            found = PETSC_TRUE          
          ! change negative id to positive, indicating it was found in database
            cur_mineral%id = abs(cur_mineral%id)
            exit
          endif
          cur_mineral => cur_mineral%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        ! read the molar volume
        call fiReadDBaseDouble(dbase_id,string,cur_mineral%molar_volume,ierr)
        call fiErrorMsg(option%myrank,'MINERAL molar volume','DATABASE',ierr)            
        ! create mienral reaction
        cur_mineral%tstrxn => TransitionStateTheoryRxnCreate()
        ! read the number of aqueous species in mineral rxn
        call fiReadDBaseInt(dbase_id,string,cur_mineral%tstrxn%nspec,ierr)
        call fiErrorMsg(option%myrank,'Number of species in mineral reaction', &
                        'DATABASE',ierr)  
        ! allocate arrays for rxn
        allocate(cur_mineral%tstrxn%spec_name(cur_mineral%tstrxn%nspec))
        cur_mineral%tstrxn%spec_name = ''
        allocate(cur_mineral%tstrxn%stoich(cur_mineral%tstrxn%nspec))
        cur_mineral%tstrxn%stoich = 0.d0
        allocate(cur_mineral%tstrxn%logK(reaction%num_dbase_temperatures))
        cur_mineral%tstrxn%logK = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_mineral%tstrxn%nspec
          call fiReadDBaseDouble(dbase_id,string,cur_mineral%tstrxn%stoich(ispec),ierr)
          call fiErrorMsg(option%myrank,'MINERAL species stoichiometry','DATABASE',ierr)            
          call fiReadDBaseName(dbase_id,string,cur_mineral%tstrxn%spec_name(ispec),.true.,ierr)
          call fiErrorMsg(option%myrank,'MINERAL species name','DATABASE',ierr)            
        enddo
        do itemp = 1, reaction%num_dbase_temperatures
          call fiReadDBaseDouble(dbase_id,string,cur_mineral%tstrxn%logK(itemp),ierr)
          call fiErrorMsg(option%myrank,'MINERAL logKs','DATABASE',ierr)            
        enddo
        ! read the molar weight
        call fiReadDBaseDouble(dbase_id,string,cur_mineral%molar_weight,ierr)
        call fiErrorMsg(option%myrank,'MINERAL molar weight','DATABASE',ierr)            
        
        
      case(4) ! surface complexes
        cur_surfcplx => reaction%surface_complex_list
        if (.not.associated(cur_surfcplx)) cycle
        found = PETSC_FALSE
        do
          if (.not.associated(cur_surfcplx)) exit
          if (fiStringCompare(name,cur_surfcplx%name,MAXNAMELENGTH)) then
            found = PETSC_TRUE          
          ! change negative id to positive, indicating it was found in database
            cur_surfcplx%id = abs(cur_surfcplx%id)
            exit
          endif
          cur_surfcplx => cur_surfcplx%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        ! read the number of aqueous species in surface complexation rxn
        call fiReadDBaseInt(dbase_id,string,cur_surfcplx%nspec,ierr)
        call fiErrorMsg(option%myrank,'Number of species in surface complexation reaction', &
                        'DATABASE',ierr)  
        ! allocate arrays for rxn
        allocate(cur_surfcplx%spec_name(cur_surfcplx%nspec))
        cur_surfcplx%spec_name = ''
        allocate(cur_surfcplx%stoich(cur_surfcplx%nspec))
        cur_surfcplx%stoich = 0.d0
        allocate(cur_surfcplx%logK(reaction%num_dbase_temperatures))
        cur_surfcplx%logK = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_surfcplx%nspec
          call fiReadDBaseDouble(dbase_id,string,cur_surfcplx%stoich(ispec),ierr)
          call fiErrorMsg(option%myrank,'SURFACE COMPLEX species stoichiometry','DATABASE',ierr)            
          call fiReadDBaseName(dbase_id,string,cur_surfcplx%spec_name(ispec),.true.,ierr)
          call fiErrorMsg(option%myrank,'SURFACE COMPLEX species name','DATABASE',ierr)            
        enddo
        do itemp = 1, reaction%num_dbase_temperatures
          call fiReadDBaseDouble(dbase_id,string,cur_surfcplx%logK(itemp),ierr)
          call fiErrorMsg(option%myrank,'SURFACE COMPLEX logKs','DATABASE',ierr)            
        enddo
        ! skip molar weight
      
    end select
    
  enddo
  
  ! check that all species, etc. were read
  flag = PETSC_FALSE
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%id < 0) then
      flag = PETSC_TRUE
      string = 'Aqueous primary species (' // trim(cur_aq_spec%name) // &
               ') not found in database.'
      call printMsg(option,string)
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%id < 0) then
      flag = PETSC_TRUE
      string = 'Aqueous secondary species (' // trim(cur_aq_spec%name) // &
               ') not found in database.'
      call printMsg(option,string)
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo  
  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    if (cur_gas_spec%id < 0) then
      flag = PETSC_TRUE
      string = 'Gas species (' // trim(cur_gas_spec%name) // &
               ') not found in database.'
      call printMsg(option,string)
    endif
    cur_gas_spec => cur_gas_spec%next
  enddo  
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%id < 0) then
      flag = PETSC_TRUE
      string = 'Mineral (' // trim(cur_mineral%name) // &
               ') not found in database.'
      call printMsg(option,string)
    endif
    cur_mineral => cur_mineral%next
  enddo
  cur_surfcplx => reaction%surface_complex_list
  do
    if (.not.associated(cur_surfcplx)) exit
    if (cur_surfcplx%id < 0) then
      flag = PETSC_TRUE
      string = 'Surface species (' // trim(cur_surfcplx%name) // &
               ') not found in database.'
      call printMsg(option,string)
    endif
    cur_surfcplx => cur_surfcplx%next
  enddo  
  
  if (flag) call printErrMsg(option,'Species not found in database.')

  close(dbase_id)

end subroutine DatabaseRead

! ************************************************************************** !
!
! BasisInit: Initializes the basis for geochemistry
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
subroutine BasisInit(reaction,option)

  use Option_module
  use Reaction_module
  use Fileio_module

  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_type), pointer :: cur_mineral
  type(surface_complexation_rxn_type), pointer :: cur_surfcplx

  character(len=MAXNAMELENGTH), allocatable :: old_basis_names(:)
  character(len=MAXNAMELENGTH), allocatable :: new_basis_names(:)
  
  PetscInt :: ispec, itemp
  PetscInt :: icount_old, icount_new
  PetscInt :: i_old, i_new
  
  PetscTruth :: compute_new_basis
  PetscTruth :: found
  
  reaction%ncomp = GetPrimarySpeciesCount(reaction)
  reaction%neqcmplx = GetSecondarySpeciesCount(reaction)

  allocate(old_basis_names(reaction%ncomp+reaction%neqcmplx))
  allocate(new_basis_names(reaction%ncomp))
  
  call BasisPrint(reaction,'Initial Basis',option)
  
  icount_old = 1
  icount_new = 1
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    icount_new = icount_new + 1
    new_basis_names(icount_new) = cur_aq_spec%name
    if (.not.associated(cur_aq_spec%eqrxn)) then
      icount_old = icount_old + 1
      old_basis_names(icount_old) = cur_aq_spec%name
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (.not.associated(cur_aq_spec%eqrxn)) then
      icount_old = icount_old + 1
      old_basis_names(icount_old) = cur_aq_spec%name
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo
  
  ! check if basis needs to be swapped
  compute_new_basis = PETSC_FALSE
  do i_old = 1, icount_old
    found = PETSC_FALSE
    do i_new = 1, icount_new
      if (fiStringCompare(old_basis_names(i_old), &
                          new_basis_names(i_new),MAXNAMELENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      compute_new_basis = PETSC_TRUE
      exit
    endif
  enddo
  
  ! check for seconday aqueous and gas species in reactions and swap them
  ! out (e.g. HS- is a secondary aqueous complex in the database, but found 
  ! in many secondary aqueous and mineral reactions)
  

  call BasisPrint(reaction,'Final Basis',option)
  
end subroutine BasisInit

! ************************************************************************** !
!
! BasisPrint: Prints the basis
! author: Glenn Hammond
! date: 09/01/08
!
! ************************************************************************** !
subroutine BasisPrint(reaction,title,option)

  use Option_module
  use Reaction_module

  implicit none
  
  type(reaction_type) :: reaction
  character(len=*) :: title
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_type), pointer :: cur_mineral
  type(surface_complexation_rxn_type), pointer :: cur_surfcplx

  PetscInt :: ispec, itemp

  if (option%myrank == 0) then
    write(IUNIT2,*) trim(title)
    write(IUNIT2,*)
    write(IUNIT2,*)
    write(IUNIT2,*)

    write(IUNIT2,*) 'Primary Components'
    cur_aq_spec => reaction%primary_species_list
    do
      if (.not.associated(cur_aq_spec)) exit
      write(IUNIT2,*) cur_aq_spec%name
      write(IUNIT2,*) '  Charge: ', cur_aq_spec%Z
      write(IUNIT2,*) '  Molar Weight: ', cur_aq_spec%molar_weight
      write(IUNIT2,*) '  Debye-Huckel a0: ', cur_aq_spec%name
      if (associated(cur_aq_spec%eqrxn)) then
        write(IUNIT2,*) '  Equilibrium Aqueous Reaction: '
        write(IUNIT2,*) '    ', -1.d0, cur_aq_spec%name
        do ispec = 1, cur_aq_spec%eqrxn%nspec
          write(IUNIT2,*) '    ', cur_aq_spec%eqrxn%stoich(ispec), &
                          cur_aq_spec%eqrxn%spec_name(ispec)
        enddo
        write(IUNIT2,*) '    logK: ', (cur_aq_spec%eqrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(IUNIT2,*)
      cur_aq_spec => cur_aq_spec%next
    enddo
    
    write(IUNIT2,*)
    write(IUNIT2,*) 'Secondary Components'
    cur_aq_spec => reaction%secondary_species_list
    do
      if (.not.associated(cur_aq_spec)) exit
      write(IUNIT2,*) cur_aq_spec%name
      write(IUNIT2,*) '  Charge: ', cur_aq_spec%Z
      write(IUNIT2,*) '  Molar Weight: ', cur_aq_spec%molar_weight
      write(IUNIT2,*) '  Debye-Huckel a0: ', cur_aq_spec%name
      if (associated(cur_aq_spec%eqrxn)) then
        write(IUNIT2,*) '  Equilibrium Aqueous Reaction: '
        write(IUNIT2,*) '    ', -1.d0, cur_aq_spec%name
        do ispec = 1, cur_aq_spec%eqrxn%nspec
          write(IUNIT2,*) '    ', cur_aq_spec%eqrxn%stoich(ispec), &
                          cur_aq_spec%eqrxn%spec_name(ispec)
        enddo
        write(IUNIT2,*) '    logK: ', (cur_aq_spec%eqrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(IUNIT2,*)
      cur_aq_spec => cur_aq_spec%next
    enddo
    
    write(IUNIT2,*)
    write(IUNIT2,*) 'Gas Components'
    cur_gas_spec => reaction%gas_species_list
    do
      if (.not.associated(cur_gas_spec)) exit
      write(IUNIT2,*) cur_gas_spec%name
      write(IUNIT2,*) '  Molar Weight: ', cur_gas_spec%molar_weight
      if (associated(cur_gas_spec%eqrxn)) then
        write(IUNIT2,*) '  Gas Reaction: '
        write(IUNIT2,*) '    ', -1.d0, cur_gas_spec%name
        do ispec = 1, cur_gas_spec%eqrxn%nspec
          write(IUNIT2,*) '    ', cur_gas_spec%eqrxn%stoich(ispec), &
                          cur_gas_spec%eqrxn%spec_name(ispec)
        enddo
        write(IUNIT2,*) '    logK: ', (cur_gas_spec%eqrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(IUNIT2,*)
      cur_gas_spec => cur_gas_spec%next
    enddo
    
    write(IUNIT2,*)
    write(IUNIT2,*) 'Minerals'
    cur_mineral => reaction%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      write(IUNIT2,*) cur_mineral%name
      write(IUNIT2,*) '  Molar Weight: ', cur_mineral%molar_weight
      write(IUNIT2,*) '  Molar Volume: ', cur_mineral%molar_volume
      if (associated(cur_mineral%tstrxn)) then
        write(IUNIT2,*) '  Mineral Reaction: '
        write(IUNIT2,*) '    ', -1.d0, cur_mineral%name
        do ispec = 1, cur_mineral%tstrxn%nspec
          write(IUNIT2,*) '    ', cur_mineral%tstrxn%stoich(ispec), &
                          cur_mineral%tstrxn%spec_name(ispec)
        enddo
        write(IUNIT2,*) '    logK: ', (cur_mineral%tstrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(IUNIT2,*)
      cur_mineral => cur_mineral%next
    enddo
    
    write(IUNIT2,*)
    write(IUNIT2,*) 'Surface Complexes'
    cur_surfcplx => reaction%surface_complex_list
    do
      if (.not.associated(cur_surfcplx)) exit
      write(IUNIT2,*) cur_surfcplx%name
      write(IUNIT2,*) '  Charge: ', cur_surfcplx%Z
      write(IUNIT2,*) '  Surface Complex Reaction: '
      write(IUNIT2,*) '    ', -1.d0, cur_surfcplx%name
      do ispec = 1, cur_surfcplx%nspec
        write(IUNIT2,*) '    ', cur_surfcplx%stoich(ispec), &
                        cur_surfcplx%spec_name(ispec)
      enddo
      write(IUNIT2,*) '    logK: ', (cur_surfcplx%logK(itemp),itemp=1, &
                                     reaction%num_dbase_temperatures)
      write(IUNIT2,*)
      cur_surfcplx => cur_surfcplx%next
    enddo
    
  endif

end subroutine BasisPrint

end module Database_module
