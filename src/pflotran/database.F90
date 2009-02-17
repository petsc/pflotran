module Database_module

  use Reaction_module
  use Reaction_Aux_module

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

  use Option_module
  use Input_module
  use String_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec, cur_aq_spec2
  type(gas_species_type), pointer :: cur_gas_spec, cur_gas_spec2
  type(mineral_type), pointer :: cur_mineral, cur_mineral2
  type(surface_complexation_rxn_type), pointer :: cur_surfcplx_rxn
  type(surface_complex_type), pointer :: cur_surfcplx, cur_surfcplx2
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: null_name
  
  PetscTruth :: flag, found
  PetscInt :: ispec, itemp, i
  PetscReal :: stoich
  type(input_type), pointer :: input
  PetscInt :: iostat
  PetscInt :: num_nulls
  
  ! negate ids for use as flags
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    cur_aq_spec%id = -abs(cur_aq_spec%id)
    cur_aq_spec => cur_aq_spec%next
  enddo
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    cur_aq_spec%id = -abs(cur_aq_spec%id)
    cur_aq_spec => cur_aq_spec%next
  enddo  
  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    cur_gas_spec%id = -abs(cur_gas_spec%id)
    cur_gas_spec => cur_gas_spec%next
  enddo  
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo
  cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_surfcplx_rxn)) exit
    cur_surfcplx => cur_surfcplx_rxn%complex_list
    do  
      if (.not.associated(cur_surfcplx)) exit
      cur_surfcplx%id = -abs(cur_surfcplx%id)
      cur_surfcplx => cur_surfcplx%next
    enddo
    cur_surfcplx_rxn => cur_surfcplx_rxn%next
  enddo  
  
  input => InputCreate(IUNIT_TEMP,reaction%database_filename)

  ! read temperatures
  call InputReadFlotranString(input,option)
  ! remove comment
  call InputReadQuotedWord(input,option,name,PETSC_TRUE)
  call InputReadInt(input,option,reaction%num_dbase_temperatures)
  call InputErrorMsg(input,option,'Number of database temperatures','DATABASE')  
  allocate(reaction%dbase_temperatures(reaction%num_dbase_temperatures))
  reaction%dbase_temperatures = 0.d0  
  do itemp = 1, reaction%num_dbase_temperatures
    call InputReadDouble(input,option,reaction%dbase_temperatures(itemp))
    call InputErrorMsg(input,option,'Database temperatures','DATABASE')            
  enddo

  num_nulls = 0
  null_name = 'null'
  do ! loop over every entry in the database
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'DATABASE')

    call InputReadQuotedWord(input,option,name,PETSC_TRUE)
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
    
    if (StringCompare(name,null_name,MAXWORDLENGTH)) then
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
          if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
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
          if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
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
          if (.not.associated(cur_aq_spec%eqrxn)) &
            cur_aq_spec%eqrxn => EquilibriumRxnCreate()
          ! read the number of primary species in secondary rxn
          call InputReadInt(input,option,cur_aq_spec%eqrxn%nspec)
          call InputErrorMsg(input,option,'Number of species in aqueous complex', &
                          'DATABASE')  
          ! allocate arrays for rxn
          allocate(cur_aq_spec%eqrxn%spec_name(cur_aq_spec%eqrxn%nspec))
          cur_aq_spec%eqrxn%spec_name = ''
          allocate(cur_aq_spec%eqrxn%stoich(cur_aq_spec%eqrxn%nspec))
          cur_aq_spec%eqrxn%stoich = 0.d0
          allocate(cur_aq_spec%eqrxn%logK(reaction%num_dbase_temperatures))
          cur_aq_spec%eqrxn%logK = 0.d0
          ! read in species and stoichiometries
          do ispec = 1, cur_aq_spec%eqrxn%nspec
            call InputReadDouble(input,option,cur_aq_spec%eqrxn%stoich(ispec))
            call InputErrorMsg(input,option,'EQRXN species stoichiometry','DATABASE')            
            call InputReadQuotedWord(input,option,cur_aq_spec%eqrxn%spec_name(ispec),PETSC_TRUE)
            call InputErrorMsg(input,option,'EQRXN species name','DATABASE')            
          enddo
          do itemp = 1, reaction%num_dbase_temperatures
            call InputReadDouble(input,option,cur_aq_spec%eqrxn%logK(itemp))
            call InputErrorMsg(input,option,'EQRXN logKs','DATABASE')            
          enddo
        endif 
        ! read the Debye-Huckel ion size parameter (a0)
        call InputReadDouble(input,option,cur_aq_spec%a0)
        call InputErrorMsg(input,option,'AQ Species a0','DATABASE')            
        ! read the valence
        call InputReadDouble(input,option,cur_aq_spec%Z)
        call InputErrorMsg(input,option,'AQ Species Z','DATABASE')            
        ! read the molar weight
        call InputReadDouble(input,option,cur_aq_spec%molar_weight)
        call InputErrorMsg(input,option,'AQ Species molar weight','DATABASE')
        
                    
      case(2) ! gas species
        cur_gas_spec => reaction%gas_species_list
        if (.not.associated(cur_gas_spec)) cycle
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_gas_spec)) exit
          if (StringCompare(name,cur_gas_spec%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
          ! change negative id to positive, indicating it was found in database
            cur_gas_spec%id = abs(cur_gas_spec%id)
            exit
          endif
          cur_gas_spec => cur_gas_spec%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        ! read the molar volume
        call InputReadDouble(input,option,cur_gas_spec%molar_volume)
        call InputErrorMsg(input,option,'GAS molar volume','DATABASE')
        ! convert from cm^3/mol to m^3/mol
        cur_gas_spec%molar_volume = cur_gas_spec%molar_volume*1.d-6
        ! create aqueous equilibrium reaction
        if (.not.associated(cur_gas_spec%eqrxn)) &
          cur_gas_spec%eqrxn => EquilibriumRxnCreate()
        ! read the number of aqueous species in secondary rxn
        call InputReadInt(input,option,cur_gas_spec%eqrxn%nspec)
        call InputErrorMsg(input,option,'Number of species in gas reaction', &
                        'DATABASE')  
        ! allocate arrays for rxn
        allocate(cur_gas_spec%eqrxn%spec_name(cur_gas_spec%eqrxn%nspec))
        cur_gas_spec%eqrxn%spec_name = ''
        allocate(cur_gas_spec%eqrxn%stoich(cur_gas_spec%eqrxn%nspec))
        cur_gas_spec%eqrxn%stoich = 0.d0
        allocate(cur_gas_spec%eqrxn%logK(reaction%num_dbase_temperatures))
        cur_gas_spec%eqrxn%logK = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_gas_spec%eqrxn%nspec
          call InputReadDouble(input,option,cur_gas_spec%eqrxn%stoich(ispec))
          call InputErrorMsg(input,option,'GAS species stoichiometry','DATABASE')            
          call InputReadQuotedWord(input,option,cur_gas_spec%eqrxn%spec_name(ispec),PETSC_TRUE)
          call InputErrorMsg(input,option,'GAS species name','DATABASE')            
        enddo
        do itemp = 1, reaction%num_dbase_temperatures
          call InputReadDouble(input,option,cur_gas_spec%eqrxn%logK(itemp))
          call InputErrorMsg(input,option,'GAS logKs','DATABASE')            
        enddo
        ! read the molar weight
        call InputReadDouble(input,option,cur_gas_spec%molar_weight)
        call InputErrorMsg(input,option,'GAS molar weight','DATABASE')     
        
               
      case(3) ! minerals
        cur_mineral => reaction%mineral_list
        if (.not.associated(cur_mineral)) cycle
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_mineral)) exit
          if (StringCompare(name,cur_mineral%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
          ! change negative id to positive, indicating it was found in database
            cur_mineral%id = abs(cur_mineral%id)
            exit
          endif
          cur_mineral => cur_mineral%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        ! read the molar volume
        call InputReadDouble(input,option,cur_mineral%molar_volume)
        call InputErrorMsg(input,option,'MINERAL molar volume','DATABASE')            
        ! convert from cm^3/mol to m^3/mol
        cur_mineral%molar_volume = cur_mineral%molar_volume*1.d-6
        ! create mineral reaction
        if (.not.associated(cur_mineral%tstrxn)) &
          cur_mineral%tstrxn => TransitionStateTheoryRxnCreate()
        ! read the number of aqueous species in mineral rxn
        call InputReadInt(input,option,cur_mineral%tstrxn%nspec)
        call InputErrorMsg(input,option,'Number of species in mineral reaction', &
                        'DATABASE')  
        ! allocate arrays for rxn
        allocate(cur_mineral%tstrxn%spec_name(cur_mineral%tstrxn%nspec))
        cur_mineral%tstrxn%spec_name = ''
        allocate(cur_mineral%tstrxn%stoich(cur_mineral%tstrxn%nspec))
        cur_mineral%tstrxn%stoich = 0.d0
        allocate(cur_mineral%tstrxn%logK(reaction%num_dbase_temperatures))
        cur_mineral%tstrxn%logK = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_mineral%tstrxn%nspec
          call InputReadDouble(input,option,cur_mineral%tstrxn%stoich(ispec))
          call InputErrorMsg(input,option,'MINERAL species stoichiometry','DATABASE')            
          call InputReadQuotedWord(input,option,cur_mineral%tstrxn%spec_name(ispec),PETSC_TRUE)
          call InputErrorMsg(input,option,'MINERAL species name','DATABASE')            
        enddo
        do itemp = 1, reaction%num_dbase_temperatures
          call InputReadDouble(input,option,cur_mineral%tstrxn%logK(itemp))
          call InputErrorMsg(input,option,'MINERAL logKs','DATABASE')            
        enddo
        ! read the molar weight
        call InputReadDouble(input,option,cur_mineral%molar_weight)
        call InputErrorMsg(input,option,'MINERAL molar weight','DATABASE')            
        
        
      case(4) ! surface complexes
        cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
        found = PETSC_FALSE
        do
          if (.not.associated(cur_surfcplx_rxn)) exit
          cur_surfcplx => cur_surfcplx_rxn%complex_list
          do
            if (.not.associated(cur_surfcplx)) exit
            if (StringCompare(name,cur_surfcplx%name,MAXWORDLENGTH)) then
              found = PETSC_TRUE          
            ! change negative id to positive, indicating it was found in database
              cur_surfcplx%id = abs(cur_surfcplx%id)
              exit
            endif
            cur_surfcplx => cur_surfcplx%next
          enddo
          if (found) exit
          cur_surfcplx_rxn => cur_surfcplx_rxn%next
        enddo
        
        if (.not.found) cycle ! go to next line in database

        if (.not.associated(cur_surfcplx%eqrxn)) &
          cur_surfcplx%eqrxn => EquilibriumRxnCreate()
            
        ! read the number of aqueous species in surface complexation rxn
        call InputReadInt(input,option,cur_surfcplx%eqrxn%nspec)
        call InputErrorMsg(input,option,'Number of species in surface complexation reaction', &
                        'DATABASE')  
        ! decrement number of species since free site will not be included
        cur_surfcplx%eqrxn%nspec = cur_surfcplx%eqrxn%nspec - 1
        ! allocate arrays for rxn
        allocate(cur_surfcplx%eqrxn%spec_name(cur_surfcplx%eqrxn%nspec))
        cur_surfcplx%eqrxn%spec_name = ''
        allocate(cur_surfcplx%eqrxn%stoich(cur_surfcplx%eqrxn%nspec))
        cur_surfcplx%eqrxn%stoich = 0.d0
        allocate(cur_surfcplx%eqrxn%logK(reaction%num_dbase_temperatures))
        cur_surfcplx%eqrxn%logK = 0.d0
        ! read in species and stoichiometries
        ispec = 0
        found = PETSC_FALSE
        do i = 1, cur_surfcplx%eqrxn%nspec+1 ! recall that nspec was decremented above
          call InputReadDouble(input,option,stoich)
          call InputErrorMsg(input,option,'SURFACE COMPLEX species stoichiometry','DATABASE')            
          call InputReadQuotedWord(input,option,name,PETSC_TRUE)
          call InputErrorMsg(input,option,'SURFACE COMPLEX species name','DATABASE')            
          if (StringCompare(name,cur_surfcplx_rxn%free_site_name,MAXWORDLENGTH)) then
            found = PETSC_TRUE
            cur_surfcplx%free_site_stoich = stoich
          else
            ispec = ispec + 1
            cur_surfcplx%eqrxn%stoich(ispec) = stoich
            cur_surfcplx%eqrxn%spec_name(ispec) = name
          endif
        enddo
        if (.not.found) then
          option%io_buffer = 'Free site name: ' // &
                             trim(cur_surfcplx_rxn%free_site_name) // &
                             ' not found in surface complex:' // &
                             trim(cur_surfcplx%name)
          call printErrMsg(option)
        endif
        do itemp = 1, reaction%num_dbase_temperatures
          call InputReadDouble(input,option,cur_surfcplx%eqrxn%logK(itemp))
          call InputErrorMsg(input,option,'SURFACE COMPLEX logKs','DATABASE')            
        enddo
        ! read the valence
        call InputReadDouble(input,option,cur_surfcplx%Z)
        call InputErrorMsg(input,option,'Surface Complex Z','DATABASE')            

      
    end select
    
  enddo
  
  ! check for duplicate species
  flag = PETSC_FALSE
 
  ! aqueous primary species
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    
    ! aqueous primary species
    cur_aq_spec2 => cur_aq_spec%next
    do
      if (.not.associated(cur_aq_spec2)) exit
      if (cur_aq_spec%id /= cur_aq_spec2%id .and. &
          StringCompare(cur_aq_spec%name, &
                          cur_aq_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = &
                 'Aqueous primary species (' // trim(cur_aq_spec%name) // &
                 ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_aq_spec2 => cur_aq_spec2%next
    enddo

    cur_aq_spec2 => reaction%secondary_species_list
    do
      if (.not.associated(cur_aq_spec2)) exit
      if (StringCompare(cur_aq_spec%name, &
                          cur_aq_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous primary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated as secondary species in input file.'
        call printMsg(option)                          
      endif
      cur_aq_spec2 => cur_aq_spec2%next
    enddo

    cur_gas_spec2 => reaction%gas_species_list
    do
      if (.not.associated(cur_gas_spec2)) exit
      if (StringCompare(cur_aq_spec%name, &
                          cur_gas_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous primary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated as gas species in input file.'
        call printMsg(option)                          
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    cur_aq_spec => cur_aq_spec%next  
  enddo
  
  ! aqueous secondary species
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    
    ! already checked against primary
    ! aqueous secondary species
    cur_aq_spec2 => cur_aq_spec%next
    do
      if (.not.associated(cur_aq_spec2)) exit
      if (cur_aq_spec%id /= cur_aq_spec2%id .and. &
          StringCompare(cur_aq_spec%name, &
                          cur_aq_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous secondary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_aq_spec2 => cur_aq_spec2%next
    enddo

    cur_gas_spec2 => reaction%gas_species_list
    do
      if (.not.associated(cur_gas_spec2)) exit
      if (StringCompare(cur_aq_spec%name, &
                          cur_gas_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous secondary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated as gas species in input file.'
        call printMsg(option)                          
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    cur_aq_spec => cur_aq_spec%next  
  enddo
  
  ! gas species
  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    
    ! already checked against primary
    ! already checked against secondary
    ! gas species
    cur_gas_spec2 => cur_gas_spec%next
    do
      if (.not.associated(cur_gas_spec2)) exit
      if (cur_gas_spec%id /= cur_gas_spec2%id .and. &
          StringCompare(cur_aq_spec%name, &
                          cur_gas_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Gas species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    cur_aq_spec => cur_aq_spec%next  
  enddo
  
  ! minerals
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral2 => cur_mineral%next
    do
      if (.not.associated(cur_mineral2)) exit
      if (cur_mineral%id /= cur_mineral2%id .and. &
          StringCompare(cur_mineral%name, &
                          cur_mineral2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Mineral (' // &
                           trim(cur_mineral%name) // &
                           ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_mineral2 => cur_mineral2%next
    enddo
    cur_mineral => cur_mineral%next
  enddo
  
  ! surface complexes
  cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_surfcplx_rxn)) exit
    cur_surfcplx => cur_surfcplx_rxn%complex_list
    do
      if (.not.associated(cur_surfcplx)) exit
      cur_surfcplx2 => cur_surfcplx%next
      do
        if (.not.associated(cur_surfcplx2)) exit
        if (cur_surfcplx%id /= cur_surfcplx2%id .and. &
            StringCompare(cur_surfcplx%name, &
                            cur_surfcplx2%name,MAXWORDLENGTH)) then
          flag = PETSC_TRUE
          option%io_buffer = 'Surface complex (' // &
                             trim(cur_surfcplx2%name) // &
                      ') duplicated in input file surface complex reaction.'
          call printMsg(option)                          
        endif
        cur_surfcplx2 => cur_surfcplx2%next
      enddo
      cur_surfcplx => cur_surfcplx%next
    enddo
    cur_surfcplx_rxn => cur_surfcplx_rxn%next
  enddo  
  
  if (flag) call printErrMsg(option,'Species duplicated in input file.')

  ! check that all species, etc. were read
  flag = PETSC_FALSE
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = 'Aqueous primary species (' // &
               trim(cur_aq_spec%name) // &
               ') not found in database.'
      call printMsg(option)
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = &
               'Aqueous secondary species (' // trim(cur_aq_spec%name) // &
               ') not found in database.'
      call printMsg(option)
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo  
  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    if (cur_gas_spec%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = 'Gas species (' // trim(cur_gas_spec%name) // &
                         ') not found in database.'
      call printMsg(option)
    endif
    cur_gas_spec => cur_gas_spec%next
  enddo  
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = 'Mineral (' // trim(cur_mineral%name) // &
               ') not found in database.'
      call printMsg(option)
    endif
    cur_mineral => cur_mineral%next
  enddo
  cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_surfcplx_rxn)) exit
    cur_surfcplx => cur_surfcplx_rxn%complex_list
    do
      if (.not.associated(cur_surfcplx)) exit
      if (cur_surfcplx%id < 0) then
        flag = PETSC_TRUE
        option%io_buffer = 'Surface species (' // trim(cur_surfcplx%name) // &
                 ') not found in database.'
        call printMsg(option)
      endif
      cur_surfcplx => cur_surfcplx%next
    enddo  
    cur_surfcplx_rxn => cur_surfcplx_rxn%next
  enddo 
    
  if (flag) call printErrMsg(option,'Species not found in database.')

  call InputDestroy(input)

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
  use String_module
  use Utility_module

  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec
  type(aq_species_type), pointer :: cur_pri_aq_spec
  type(aq_species_type), pointer :: cur_sec_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_type), pointer :: cur_mineral
  type(aq_species_type), pointer :: cur_sec_aq_spec1
  type(aq_species_type), pointer :: cur_sec_aq_spec2
  type(gas_species_type), pointer :: cur_gas_spec1
  type(gas_species_type), pointer :: cur_gas_spec2
  type(surface_complexation_rxn_type), pointer :: cur_surfcplx_rxn
  type(surface_complex_type), pointer :: cur_surfcplx
  type(surface_complex_type), pointer :: cur_surfcplx2
  type(ion_exchange_rxn_type), pointer :: cur_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cur_cation
  
  character(len=MAXWORDLENGTH), allocatable :: old_basis_names(:)
  character(len=MAXWORDLENGTH), allocatable :: new_basis_names(:)

  character(len=MAXWORDLENGTH), parameter :: h2oname = 'H2O'
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt, parameter :: h2o_id = 1

  PetscReal :: logK(reaction%num_dbase_temperatures)
  PetscReal, allocatable :: transformation(:,:), old_basis(:,:), new_basis(:,:)
  PetscReal, allocatable :: stoich_new(:), stoich_prev(:), logKvector(:,:)
  PetscInt, allocatable :: indices(:)
  
  PetscReal, allocatable :: pri_matrix(:,:), sec_matrix(:,:)
  PetscReal, allocatable :: sec_matrix_inverse(:,:)
  PetscReal, allocatable :: stoich_matrix(:,:)
  PetscReal, allocatable :: unit_vector(:)
  character(len=MAXWORDLENGTH), allocatable :: pri_names(:)
  character(len=MAXWORDLENGTH), allocatable :: sec_names(:)
  character(len=MAXWORDLENGTH), allocatable :: gas_names(:)
  PetscReal, allocatable :: logKvector_swapped(:,:)
  
  PetscInt :: ispec, itemp
  PetscInt :: spec_id
  PetscInt :: ncomp_h2o, ncomp_secondary
  PetscInt :: icount_old, icount_new, icount, icount2
  PetscInt :: i, j, irow, icol
  PetscInt :: ipri_spec, isec_spec, imnrl, igas_spec, ikinmnrl
  PetscInt :: i_old, i_new
  PetscInt :: isurfcplx, irxn
  PetscInt :: ication
  PetscInt :: idum
  PetscReal :: temp_high, temp_low
  PetscInt :: itemp_high, itemp_low
  
  PetscTruth :: compute_new_basis
  PetscTruth :: found

! get database temperature based on REFERENCE_TEMPERATURE
  if (option%reference_temperature <= 0.01d0) then
    reaction%debyeA = 0.4939d0 
    reaction%debyeB = 0.3253d0 
    reaction%debyeBdot = 0.0374d0
  else if (option%reference_temperature > 0.d0 .and. &
           option%reference_temperature <= 25.d0) then
    temp_low = 0.d0
    temp_high = 25.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.5114d0,0.4939d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3288d0,0.3253d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0410d0,0.0374d0,reaction%debyeBdot)
  else if (option%reference_temperature > 25.d0 .and. &
           option%reference_temperature <= 60.d0) then
    temp_low = 25.d0
    temp_high = 60.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.5465d0,0.5114d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3346d0,0.3288d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0440d0,0.0410d0,reaction%debyeBdot)
  else if (option%reference_temperature > 60.d0 .and. &
           option%reference_temperature <= 100.d0) then
    temp_low = 60.d0
    temp_high = 100.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.5995d0,0.5465d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3421d0,0.3346d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0460d0,0.0440d0,reaction%debyeBdot)
  else if (option%reference_temperature > 100.d0 .and. &
           option%reference_temperature <= 150.d0) then
    temp_low = 100.d0
    temp_high = 150.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.6855d0,0.5995d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3525d0,0.3421d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0470d0,0.0460d0,reaction%debyeBdot)
  else if (option%reference_temperature > 150.d0 .and. &
           option%reference_temperature <= 200.d0) then
    temp_low = 150.d0
    temp_high = 200.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.7994d0,0.6855d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3639d0,0.3525d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0470d0,0.0470d0,reaction%debyeBdot)
  else if (option%reference_temperature > 200.d0 .and. &
           option%reference_temperature <= 250.d0) then
    temp_low = 200.d0
    temp_high = 250.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.9593d0,0.7994d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3766d0,0.3639d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0340d0,0.0470d0,reaction%debyeBdot)
  else if (option%reference_temperature > 250.d0 .and. &
           option%reference_temperature <= 300.d0) then
    temp_low = 250.d0
    temp_high = 300.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     1.2180d0,0.9593d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3925d0,0.3766d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0000d0,0.0340d0,reaction%debyeBdot)
  else if (option%reference_temperature > 300.d0 .and. &
           option%reference_temperature <= 350.d0) then
    temp_low = 300.d0
    temp_high = 350.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     1.2180d0,1.2180d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3925d0,0.3925d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0000d0,0.0000d0,reaction%debyeBdot)
  else if (option%reference_temperature > 350.d0) then
    reaction%debyeA = 1.2180d0 
    reaction%debyeB = 0.3925d0 
    reaction%debyeBdot = 0.0000d0
  endif

  if (option%reference_temperature <= reaction%dbase_temperatures(1)) then
    itemp_low = 1
    itemp_high = 1
    temp_low = reaction%dbase_temperatures(itemp_low)
    temp_high = reaction%dbase_temperatures(itemp_high)
  else if (option%reference_temperature >= &
           reaction%dbase_temperatures(reaction%num_dbase_temperatures)) then
    itemp_low = reaction%num_dbase_temperatures
    itemp_high = reaction%num_dbase_temperatures
    temp_low = reaction%dbase_temperatures(itemp_low)
    temp_high = reaction%dbase_temperatures(itemp_high)
  else
    do itemp = 1, reaction%num_dbase_temperatures-1
      itemp_low = itemp
      itemp_high = itemp+1
      temp_low = reaction%dbase_temperatures(itemp_low)
      temp_high = reaction%dbase_temperatures(itemp_high)
      if (option%reference_temperature > temp_low .and. &
          option%reference_temperature <= temp_high) then
        exit
      endif
    enddo
  endif

  reaction%ncomp = GetPrimarySpeciesCount(reaction)
  reaction%neqcmplx = GetSecondarySpeciesCount(reaction)
  reaction%ngas = GetGasCount(reaction)

  ! account for H2O in the basis by adding 1
  ncomp_h2o = reaction%ncomp+1
  
  allocate(old_basis_names(ncomp_h2o+reaction%neqcmplx))
  allocate(new_basis_names(ncomp_h2o))
  old_basis_names = ''
  new_basis_names = ''
  
  call BasisPrint(reaction,'Initial Basis',option)
  
  !--------------------------------------------
#if 1 

  ncomp_secondary = reaction%neqcmplx+reaction%ngas
  allocate(pri_matrix(ncomp_secondary,ncomp_h2o))
  pri_matrix = 0.d0
  allocate(pri_names(ncomp_h2o))
  pri_names = ''
  allocate(sec_matrix(ncomp_secondary,ncomp_secondary))
  sec_matrix = 0.d0
  allocate(sec_names(reaction%neqcmplx))
  sec_names = ''
  allocate(gas_names(reaction%ngas))
  gas_names = ''
  
  allocate(logKvector(reaction%num_dbase_temperatures,ncomp_secondary))
  logKvector = 0.d0
  
  ! fill in names
  icount = 1
  pri_names(icount) = h2oname
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    icount = icount + 1
    pri_names(icount) = cur_aq_spec%name
    cur_aq_spec => cur_aq_spec%next
  enddo
  icount = 0
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    icount = icount + 1
    sec_names(icount) = cur_aq_spec%name
    cur_aq_spec => cur_aq_spec%next
  enddo
  icount= 0
  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    icount = icount + 1
    gas_names(icount) = cur_gas_spec%name
    cur_gas_spec => cur_gas_spec%next
  enddo
  
  ! fill in matrices
  icount = 0
  cur_pri_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    if (associated(cur_pri_aq_spec%eqrxn)) then
      icount = icount + 1
      logKvector(:,icount) = cur_pri_aq_spec%eqrxn%logK
      i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                            cur_pri_aq_spec%name, &
                            cur_pri_aq_spec%name, &
                            pri_names,sec_names,gas_names)
      if (i < 0) then
        option%io_buffer = 'Primary species ' // &
                 trim(cur_pri_aq_spec%name) // &
                 ' found in secondary or gas list.'
        call printErrMsg(option)
      endif
      pri_matrix(icount,i) = -1.d0
      do ispec=1,cur_pri_aq_spec%eqrxn%nspec
        i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                              cur_pri_aq_spec%name, &
                              cur_pri_aq_spec%eqrxn%spec_name(ispec), &
                              pri_names,sec_names,gas_names)
        if (i > 0) then
          pri_matrix(icount,i) = cur_pri_aq_spec%eqrxn%stoich(ispec)
        else
          sec_matrix(icount,-i) = cur_pri_aq_spec%eqrxn%stoich(ispec)
        endif
      enddo
    endif
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo

  cur_sec_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_sec_aq_spec)) exit
    if (associated(cur_sec_aq_spec%eqrxn)) then
      icount = icount + 1
      logKvector(:,icount) = cur_sec_aq_spec%eqrxn%logK
      i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                            cur_sec_aq_spec%name, &
                            cur_sec_aq_spec%name, &
                            pri_names,sec_names,gas_names)
      if (i > 0) then
        option%io_buffer = 'Secondary aqueous species ' // &
                 trim(cur_sec_aq_spec%name) // &
                 ' found in primary species list.'
        call printErrMsg(option)
      endif
      sec_matrix(icount,-i) = -1.d0
      do ispec=1,cur_sec_aq_spec%eqrxn%nspec
        i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                              cur_sec_aq_spec%name, &
                              cur_sec_aq_spec%eqrxn%spec_name(ispec), &
                              pri_names,sec_names,gas_names)
        if (i > 0) then
          pri_matrix(icount,i) = cur_sec_aq_spec%eqrxn%stoich(ispec)
        else
          sec_matrix(icount,-i) = cur_sec_aq_spec%eqrxn%stoich(ispec)
        endif
      enddo
    endif
    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo

  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    if (associated(cur_gas_spec%eqrxn)) then
      icount = icount + 1
      logKvector(:,icount) = cur_gas_spec%eqrxn%logK
      i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                            cur_gas_spec%name, &
                            cur_gas_spec%name, &
                            pri_names,sec_names,gas_names)
      if (i > 0) then
        option%io_buffer = 'Gas species ' // &
                 trim(cur_gas_spec%name) // &
                 ' found in primary species list.'
        call printErrMsg(option)
      endif
      sec_matrix(icount,-i) = -1.d0
      do ispec=1,cur_gas_spec%eqrxn%nspec
        i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                              cur_gas_spec%name, &
                              cur_gas_spec%eqrxn%spec_name(ispec), &
                              pri_names,sec_names,gas_names)
        if (i > 0) then
          pri_matrix(icount,i) = cur_gas_spec%eqrxn%stoich(ispec)
        else
          sec_matrix(icount,-i) = cur_gas_spec%eqrxn%stoich(ispec)
        endif
      enddo
    endif
    cur_gas_spec => cur_gas_spec%next
  enddo

  ncomp_secondary = reaction%neqcmplx+reaction%ngas
  allocate(indices(ncomp_secondary))
  indices = 0
  allocate(unit_vector(ncomp_secondary))
  unit_vector = 0.d0
  allocate(sec_matrix_inverse(ncomp_secondary,ncomp_secondary))
  sec_matrix_inverse = 0.d0
 
  call ludcmp(sec_matrix,ncomp_secondary,indices,idum)
  do ispec = 1, ncomp_secondary
    unit_vector = 0.d0
    unit_vector(ispec) = 1.d0
    call lubksb(sec_matrix,ncomp_secondary,indices,unit_vector)
    sec_matrix_inverse(:,ispec) = unit_vector(:)
  enddo

  ! invert the secondary species matrix
  allocate(stoich_matrix(ncomp_secondary,ncomp_h2o))
  stoich_matrix = 0.d0
  do j = 1, ncomp_h2o
    do i = 1, ncomp_secondary
      do ispec = 1, ncomp_secondary
        stoich_matrix(i,j) = stoich_matrix(i,j) + &
          sec_matrix_inverse(i,ispec)*pri_matrix(ispec,j)
      enddo
    enddo
  enddo
  stoich_matrix = -1.d0*stoich_matrix

  allocate(logKvector_swapped(reaction%num_dbase_temperatures,ncomp_secondary))
  logKvector_swapped = 0.d0
  
  do j = 1, ncomp_secondary
    do i = 1, reaction%num_dbase_temperatures
      logKvector_swapped(i,j) = logKvector_swapped(i,j) - &
        dot_product(sec_matrix_inverse(j,1:ncomp_secondary), &
                    logKvector(i,1:ncomp_secondary))
    enddo
  enddo
    
  deallocate(pri_matrix)
  deallocate(sec_matrix)
  deallocate(indices)
  deallocate(unit_vector)
  deallocate(sec_matrix_inverse)
  deallocate(logKvector)
  
  cur_pri_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    if (associated(cur_pri_aq_spec%eqrxn)) then
      call EquilibriumRxnDestroy(cur_pri_aq_spec%eqrxn)
    endif
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo

  icount = 0
  cur_sec_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_sec_aq_spec)) exit
    icount = icount + 1
    ! destory old reaction
    call EquilibriumRxnDestroy(cur_sec_aq_spec%eqrxn)
    ! allocate new
    cur_sec_aq_spec%eqrxn => EquilibriumRxnCreate()

    ! count # of species in reaction
    icount2 = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        cur_sec_aq_spec%eqrxn%nspec = cur_sec_aq_spec%eqrxn%nspec + 1
      endif
    enddo
    
    allocate(cur_sec_aq_spec%eqrxn%stoich(cur_sec_aq_spec%eqrxn%nspec))
    cur_sec_aq_spec%eqrxn%stoich = 0.d0
    allocate(cur_sec_aq_spec%eqrxn%spec_name(cur_sec_aq_spec%eqrxn%nspec))
    cur_sec_aq_spec%eqrxn%spec_name = ''
    allocate(cur_sec_aq_spec%eqrxn%spec_ids(cur_sec_aq_spec%eqrxn%nspec))
    cur_sec_aq_spec%eqrxn%spec_ids = 0
    allocate(cur_sec_aq_spec%eqrxn%logK(reaction%num_dbase_temperatures))
    cur_sec_aq_spec%eqrxn%logK = 0.d0

    ispec = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        ispec = ispec + 1
        cur_sec_aq_spec%eqrxn%spec_name(ispec) = pri_names(icol)
        cur_sec_aq_spec%eqrxn%stoich(ispec) = stoich_matrix(icount,icol)
        cur_sec_aq_spec%eqrxn%spec_ids(ispec) = icol
      endif
    enddo

    cur_sec_aq_spec%eqrxn%logK = logKvector_swapped(:,icount)

    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo

  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    icount = icount + 1
    ! destory old reaction
    call EquilibriumRxnDestroy(cur_gas_spec%eqrxn)
    ! allocate new
    cur_gas_spec%eqrxn => EquilibriumRxnCreate()

    ! count # of species in reaction
    icount2 = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        cur_gas_spec%eqrxn%nspec = cur_gas_spec%eqrxn%nspec + 1
      endif
    enddo
    
    allocate(cur_gas_spec%eqrxn%stoich(cur_gas_spec%eqrxn%nspec))
    cur_gas_spec%eqrxn%stoich = 0.d0
    allocate(cur_gas_spec%eqrxn%spec_name(cur_gas_spec%eqrxn%nspec))
    cur_gas_spec%eqrxn%spec_name = ''
    allocate(cur_gas_spec%eqrxn%spec_ids(cur_gas_spec%eqrxn%nspec))
    cur_gas_spec%eqrxn%spec_ids = 0
    allocate(cur_gas_spec%eqrxn%logK(reaction%num_dbase_temperatures))
    cur_gas_spec%eqrxn%logK = 0.d0

    ispec = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        ispec = ispec + 1
        cur_gas_spec%eqrxn%spec_name(ispec) = pri_names(icol)
        cur_gas_spec%eqrxn%stoich(ispec) = stoich_matrix(icount,icol)
        cur_gas_spec%eqrxn%spec_ids(ispec) = icol
      endif
    enddo

    cur_gas_spec%eqrxn%logK = logKvector_swapped(:,icount)

    cur_gas_spec => cur_gas_spec%next
  enddo

  new_basis_names = pri_names

  deallocate(stoich_matrix)
  deallocate(logKvector_swapped)

  deallocate(pri_names)
  deallocate(sec_names)
  deallocate(gas_names)

  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_surfcplx_rxn)
  nullify(cur_surfcplx)
    
  ! first off, lets remove all the secondary gases from all other reactions
  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    
    ! gases in mineral reactions
    cur_mineral => reaction%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%tstrxn%nspec) exit
          if (StringCompare(cur_gas_spec%name, &
                              cur_mineral%tstrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn(cur_gas_spec%name, &
                                             cur_gas_spec%eqrxn, &
                                             cur_mineral%tstrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo
    nullify(cur_mineral)

    ! gases in surface complex reactions
    cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_surfcplx_rxn)) exit
      cur_surfcplx2 => cur_surfcplx_rxn%complex_list
      do
        if (.not.associated(cur_surfcplx2)) exit
        
        if (associated(cur_surfcplx2%eqrxn)) then
          ispec = 1
          do
            if (ispec > cur_surfcplx2%eqrxn%nspec) exit
            if (StringCompare(cur_gas_spec%name, &
                                cur_surfcplx2%eqrxn%spec_name(ispec), &
                                MAXWORDLENGTH)) then
              call BasisSubSpeciesInGasOrSecRxn(cur_gas_spec%name, &
                                                cur_gas_spec%eqrxn, &
                                                cur_surfcplx2%eqrxn)
              ispec = 0
            endif
            ispec = ispec + 1
          enddo
        endif
        cur_surfcplx2 => cur_surfcplx2%next
      enddo
      nullify(cur_surfcplx2)
      cur_surfcplx_rxn => cur_surfcplx_rxn%next
    enddo
    nullify(cur_surfcplx_rxn)

    cur_gas_spec => cur_gas_spec%next
  enddo

  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_surfcplx_rxn)
  nullify(cur_surfcplx)

  ! secondary aqueous species
  cur_sec_aq_spec => reaction%secondary_species_list
  do

    if (.not.associated(cur_sec_aq_spec)) exit
    
    ! secondary aqueous species in mineral reactions
    cur_mineral => reaction%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%tstrxn%nspec) exit
          if (StringCompare(cur_sec_aq_spec%name, &
                              cur_mineral%tstrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn(cur_sec_aq_spec%name, &
                                             cur_sec_aq_spec%eqrxn, &
                                             cur_mineral%tstrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo

    ! secondary aqueous species in surface complex reactions
    cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_surfcplx_rxn)) exit
      cur_surfcplx2 => cur_surfcplx_rxn%complex_list
      do
        if (.not.associated(cur_surfcplx2)) exit
        
        if (associated(cur_surfcplx2%eqrxn)) then
          ispec = 1
          do
            if (ispec > cur_surfcplx2%eqrxn%nspec) exit
            if (StringCompare(cur_sec_aq_spec%name, &
                                cur_surfcplx2%eqrxn%spec_name(ispec), &
                                MAXWORDLENGTH)) then
              call BasisSubSpeciesInGasOrSecRxn(cur_sec_aq_spec%name, &
                                                cur_sec_aq_spec%eqrxn, &
                                                cur_surfcplx2%eqrxn)
              ispec = 0
            endif
            ispec = ispec + 1
          enddo
        endif
        cur_surfcplx2 => cur_surfcplx2%next
      enddo
      nullify(cur_surfcplx2)
      cur_surfcplx_rxn => cur_surfcplx_rxn%next
    enddo
    nullify(cur_surfcplx_rxn)    
    
    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo
  
  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_surfcplx_rxn)
  nullify(cur_surfcplx)

#else
  
  !---------------------------------------------
  icount_old = 0
  icount_new = 0
  
  icount_new = icount_new + 1
  new_basis_names(icount_new) = h2oname
  icount_old = icount_old + 1
  old_basis_names(icount_old) = h2oname
  
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

  nullify(cur_aq_spec)
  nullify(cur_pri_aq_spec)
  nullify(cur_sec_aq_spec)
  nullify(cur_sec_aq_spec1)
  nullify(cur_sec_aq_spec2)
  nullify(cur_gas_spec)
  nullify(cur_gas_spec1)
  nullify(cur_gas_spec2)
  nullify(cur_mineral)
  nullify(cur_surfcplx_rxn)
  nullify(cur_surfcplx)
  nullify(cur_surfcplx2)
    
  ! first off, lets remove all the secondary gases from all other reactions
  cur_gas_spec1 => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec1)) exit
    
    if (.not.associated(cur_gas_spec1%eqrxn)) then
      cur_gas_spec1 => cur_gas_spec1%next
      cycle
    endif
    
    ! gases in primary aqueous reactions
    cur_pri_aq_spec => reaction%primary_species_list
    do
      if (.not.associated(cur_pri_aq_spec)) exit
      
      if (associated(cur_pri_aq_spec%eqrxn)) then
        ispec = 1
        do
          if (ispec > cur_pri_aq_spec%eqrxn%nspec) exit
          if (StringCompare(cur_gas_spec1%name, &
                              cur_pri_aq_spec%eqrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInGasOrSecRxn(cur_gas_spec1%name, &
                                              cur_gas_spec1%eqrxn, &
                                              cur_pri_aq_spec%eqrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_pri_aq_spec => cur_pri_aq_spec%next
    enddo
    nullify(cur_pri_aq_spec)

    ! gases in secondary aqueous reactions
    cur_sec_aq_spec2 => reaction%secondary_species_list
    do
      if (.not.associated(cur_sec_aq_spec2)) exit
      
      if (associated(cur_sec_aq_spec2%eqrxn)) then
        ispec = 1
        do
          if (ispec > cur_sec_aq_spec2%eqrxn%nspec) exit
          if (StringCompare(cur_gas_spec1%name, &
                              cur_sec_aq_spec2%eqrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInGasOrSecRxn(cur_gas_spec1%name, &
                                              cur_gas_spec1%eqrxn, &
                                              cur_sec_aq_spec2%eqrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_sec_aq_spec2 => cur_sec_aq_spec2%next
    enddo
    nullify(cur_sec_aq_spec2)

    cur_gas_spec1 => cur_gas_spec1%next
  enddo

  nullify(cur_aq_spec)
  nullify(cur_pri_aq_spec)
  nullify(cur_sec_aq_spec)
  nullify(cur_sec_aq_spec1)
  nullify(cur_sec_aq_spec2)
  nullify(cur_gas_spec)
  nullify(cur_gas_spec1)
  nullify(cur_gas_spec2)
  nullify(cur_mineral)
  nullify(cur_surfcplx_rxn)
  nullify(cur_surfcplx)
  nullify(cur_surfcplx2)

  ! check if basis needs to be swapped
  compute_new_basis = PETSC_FALSE
  do i_old = 1, icount_old
    found = PETSC_FALSE
    do i_new = 1, icount_new
      if (StringCompare(old_basis_names(i_old), &
                          new_basis_names(i_new),MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      compute_new_basis = PETSC_TRUE
      exit
    endif
  enddo

  allocate(stoich_prev(ncomp_h2o))
  stoich_prev = 0.d0
  allocate(stoich_new(ncomp_h2o))
  stoich_new = 0.d0

  if (compute_new_basis) then
  
    allocate(new_basis(ncomp_h2o,ncomp_h2o))
    new_basis = 0.d0
    allocate(transformation(ncomp_h2o,ncomp_h2o))
    transformation = 0.d0
    allocate(old_basis(ncomp_h2o,ncomp_h2o))
    old_basis = 0.d0
    allocate(logKvector(reaction%num_dbase_temperatures,ncomp_h2o))
    logKvector = 0.d0
    allocate(indices(ncomp_h2o))
    indices = 0
    
    ! account for H2O
    new_basis(1,1) = 1.d0

    ipri_spec = 2 ! since water is 1
    cur_pri_aq_spec => reaction%primary_species_list
    do
      if (.not.associated(cur_pri_aq_spec)) exit
      do irow = 1, ncomp_h2o
        if (StringCompare(cur_pri_aq_spec%name,new_basis_names(irow), &
                            MAXWORDLENGTH)) then
          if (associated(cur_pri_aq_spec%eqrxn)) then
            logKvector(:,ipri_spec) = &
              cur_pri_aq_spec%eqrxn%logK(:)
            do i=1,cur_pri_aq_spec%eqrxn%nspec
              found = PETSC_FALSE
              do icol = 1, icount_old
                if (StringCompare(cur_pri_aq_spec%eqrxn%spec_name(i), &
                                    old_basis_names(icol), &
                                    MAXWORDLENGTH)) then
                  new_basis(irow,icol) = cur_pri_aq_spec%eqrxn%stoich(i)
                  found = PETSC_TRUE
                  exit
                endif
              enddo
              if (.not.found) then
                string = 'One or more species not found in problem statement' // &
                         ' for species: ' // &
                         trim(cur_pri_aq_spec%name) // ' ('
                do j=1,cur_pri_aq_spec%eqrxn%nspec
                  found = PETSC_FALSE
                  do icol = 1, icount_old
                    if (StringCompare(cur_pri_aq_spec%eqrxn%spec_name(j), &
                                        old_basis_names(icol), &
                                      MAXWORDLENGTH)) then
                      found = PETSC_TRUE
                      exit
                    endif
                  enddo
                  if (.not.found) then
                    string = trim(string) // ' ' // &
                             trim(cur_pri_aq_spec%eqrxn%spec_name(j))
                  endif
                enddo
                string = trim(string) // ' )'
                call printErrMsg(option,string)             
              endif
            enddo
          else
            logKvector(:,ipri_spec) = 0.d0
            do icol = 1, icount_old
              if (StringCompare(new_basis_names(irow), &
                                  old_basis_names(icol), &
                                  MAXWORDLENGTH)) then
                new_basis(irow,icol) = 1.d0
                exit
              endif            
            enddo
          endif
        endif
      enddo
      cur_pri_aq_spec => cur_pri_aq_spec%next
      ipri_spec = ipri_spec + 1
    enddo
    
    old_basis = new_basis
    
    ! solve system of equations for the new basis and create a transformation
    ! matrix
    call ludcmp(new_basis,ncomp_h2o,indices,idum)
    do ispec = 1, ncomp_h2o
      stoich_prev = 0.d0
      stoich_prev(ispec) = 1.d0
      call lubksb(new_basis,ncomp_h2o,indices,stoich_prev)
      transformation(:,ispec) = stoich_prev(:)
    enddo
  
    ! at this point, any primary species with eqrxns can have them removed
    ! and secondary species without eqrxns need one added

    cur_pri_aq_spec => reaction%primary_species_list
    do
      if (.not.associated(cur_pri_aq_spec)) exit
      if (associated(cur_pri_aq_spec%eqrxn)) then
        call EquilibriumRxnDestroy(cur_pri_aq_spec%eqrxn)
      endif
      cur_pri_aq_spec => cur_pri_aq_spec%next
    enddo

    cur_sec_aq_spec => reaction%secondary_species_list
    do
      if (.not.associated(cur_sec_aq_spec)) exit
      if (.not.associated(cur_sec_aq_spec%eqrxn)) then
        cur_sec_aq_spec%eqrxn => EquilibriumRxnCreate()
        stoich_prev = 0.d0
        logK = 0.d0
        found = PETSC_FALSE
        do icol = 1, icount_old
          if (StringCompare(cur_sec_aq_spec%name, &
                              old_basis_names(icol),MAXWORDLENGTH)) then
            stoich_prev(icol) = 1.d0
            found = PETSC_TRUE
            exit
          endif
        enddo  
        if (.not.found) then
          string = 'Species ' // trim(cur_sec_aq_spec%name) // ', which is' // &
                   ' being swapped out of the basis was not found in' // &
                   ' original basis.'
          call printErrMsg(option,string)        
        endif   
        do icol = 1, ncomp_h2o
          stoich_new(icol) = &
            dot_product(transformation(1:ncomp_h2o,icol), &
                        stoich_prev(1:ncomp_h2o))
        enddo
        do i = 1, reaction%num_dbase_temperatures
          logK(i) = logK(i) - &
            dot_product(stoich_new(1:ncomp_h2o),logKvector(i,1:ncomp_h2o))
        enddo

        ! count # of species in reaction
        do icol = 1, ncomp_h2o
          if (dabs(stoich_new(icol)) > 1.d-40) &
            cur_sec_aq_spec%eqrxn%nspec = cur_sec_aq_spec%eqrxn%nspec + 1
        enddo
        allocate(cur_sec_aq_spec%eqrxn%stoich(cur_sec_aq_spec%eqrxn%nspec))
        cur_sec_aq_spec%eqrxn%stoich = 0.d0
        allocate(cur_sec_aq_spec%eqrxn%spec_name(cur_sec_aq_spec%eqrxn%nspec))
        cur_sec_aq_spec%eqrxn%spec_name = ''
        allocate(cur_sec_aq_spec%eqrxn%spec_ids(cur_sec_aq_spec%eqrxn%nspec))
        cur_sec_aq_spec%eqrxn%spec_ids = 0
        allocate(cur_sec_aq_spec%eqrxn%logK(reaction%num_dbase_temperatures))
        cur_sec_aq_spec%eqrxn%logK = 0.d0

        ispec = 0
        do icol = 1, ncomp_h2o
          if (dabs(stoich_new(icol)) > 1.d-40) then
            ispec = ispec + 1
            cur_sec_aq_spec%eqrxn%spec_name(ispec) = new_basis_names(icol)
            cur_sec_aq_spec%eqrxn%stoich(ispec) = stoich_new(icol)
            cur_sec_aq_spec%eqrxn%spec_ids(ispec) = icol
          endif
        enddo
        cur_sec_aq_spec%eqrxn%logK = logK
    
      endif
      cur_sec_aq_spec => cur_sec_aq_spec%next
    enddo
  endif

  ! now substitute in secondary aqueous species and gases

  ! check for secondary aqueous and gas species in reactions and swap them
  ! out (e.g. HS- is a secondary aqueous complex in the database, but found 
  ! in many secondary aqueous and mineral reactions)

  nullify(cur_aq_spec)
  nullify(cur_pri_aq_spec)
  nullify(cur_sec_aq_spec)
  nullify(cur_sec_aq_spec1)
  nullify(cur_sec_aq_spec2)
  nullify(cur_gas_spec)
  nullify(cur_gas_spec1)
  nullify(cur_gas_spec2)
  nullify(cur_mineral)
    
  ! first off, lets remove all the secondary gases from all other reactions
  cur_gas_spec1 => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec1)) exit
    
    if (.not.associated(cur_gas_spec1%eqrxn)) then
      cur_gas_spec1 => cur_gas_spec1%next
      cycle
    endif
    
    ! gases in gas reactions
    cur_gas_spec2 => reaction%gas_species_list
    do
      if (.not.associated(cur_gas_spec2)) exit
      
      if (associated(cur_gas_spec2%eqrxn)) then
        ispec = 1
        do
          if (ispec > cur_gas_spec2%eqrxn%nspec) exit
          if (StringCompare(cur_gas_spec1%name, &
                              cur_gas_spec2%eqrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInGasOrSecRxn(cur_gas_spec1%name, &
                                              cur_gas_spec1%eqrxn, &
                                              cur_gas_spec2%eqrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    nullify(cur_gas_spec2)

    ! gases in secondary aqueous reactions
    cur_sec_aq_spec2 => reaction%secondary_species_list
    do
      if (.not.associated(cur_sec_aq_spec2)) exit
      
      if (associated(cur_sec_aq_spec2%eqrxn)) then
        ispec = 1
        do
          if (ispec > cur_sec_aq_spec2%eqrxn%nspec) exit
          if (StringCompare(cur_gas_spec1%name, &
                              cur_sec_aq_spec2%eqrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInGasOrSecRxn(cur_gas_spec1%name, &
                                              cur_gas_spec1%eqrxn, &
                                              cur_sec_aq_spec2%eqrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_sec_aq_spec2 => cur_sec_aq_spec2%next
    enddo
    nullify(cur_sec_aq_spec2)

    ! gases in mineral reactions
    cur_mineral => reaction%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%tstrxn%nspec) exit
          if (StringCompare(cur_gas_spec1%name, &
                              cur_mineral%tstrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn(cur_gas_spec1%name, &
                                             cur_gas_spec1%eqrxn, &
                                             cur_mineral%tstrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo
    nullify(cur_mineral)

    ! gases in surface complex reactions
    cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_surfcplx_rxn)) exit
      cur_surfcplx2 => cur_surfcplx_rxn%complex_list
      do
        if (.not.associated(cur_surfcplx2)) exit
        
        if (associated(cur_surfcplx2%eqrxn)) then
          ispec = 1
          do
            if (ispec > cur_surfcplx2%eqrxn%nspec) exit
            if (StringCompare(cur_gas_spec1%name, &
                                cur_surfcplx2%eqrxn%spec_name(ispec), &
                                MAXWORDLENGTH)) then
              call BasisSubSpeciesInGasOrSecRxn(cur_gas_spec1%name, &
                                                cur_gas_spec1%eqrxn, &
                                                cur_surfcplx2%eqrxn)
              ispec = 0
            endif
            ispec = ispec + 1
          enddo
        endif
        cur_surfcplx2 => cur_surfcplx2%next
      enddo
      nullify(cur_surfcplx2)
      cur_surfcplx_rxn => cur_surfcplx_rxn%next
    enddo
    nullify(cur_surfcplx_rxn)

    cur_gas_spec1 => cur_gas_spec1%next
  enddo

  nullify(cur_aq_spec)
  nullify(cur_pri_aq_spec)
  nullify(cur_sec_aq_spec)
  nullify(cur_sec_aq_spec1)
  nullify(cur_sec_aq_spec2)
  nullify(cur_gas_spec)
  nullify(cur_gas_spec1)
  nullify(cur_gas_spec2)
  nullify(cur_mineral)
  nullify(cur_surfcplx_rxn)
  nullify(cur_surfcplx)
  nullify(cur_surfcplx2)

  ! secondary aqueous species
  cur_sec_aq_spec1 => reaction%secondary_species_list
  do

    if (.not.associated(cur_sec_aq_spec1)) exit
    
    if (.not.associated(cur_sec_aq_spec1%eqrxn)) then
      cur_sec_aq_spec1 => cur_sec_aq_spec1%next
      cycle
    endif
    
    ! secondary aqueous species in gas reactions
    cur_gas_spec2 => reaction%gas_species_list
    do
      if (.not.associated(cur_gas_spec2)) exit
      
      if (associated(cur_gas_spec2%eqrxn)) then
        ispec = 1
        do
          if (ispec > cur_gas_spec2%eqrxn%nspec) exit
          if (StringCompare(cur_sec_aq_spec1%name, &
                              cur_gas_spec2%eqrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInGasOrSecRxn(cur_sec_aq_spec1%name, &
                                              cur_sec_aq_spec1%eqrxn, &
                                              cur_gas_spec2%eqrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    nullify(cur_gas_spec2)

    ! secondary aqueous species in secondary aqueous reactions
    cur_sec_aq_spec2 => reaction%secondary_species_list
    do
      if (.not.associated(cur_sec_aq_spec2)) exit
      
      if (associated(cur_sec_aq_spec2%eqrxn)) then
        ispec = 1
        do
          if (ispec > cur_sec_aq_spec2%eqrxn%nspec) exit
          if (StringCompare(cur_sec_aq_spec1%name, &
                              cur_sec_aq_spec2%eqrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInGasOrSecRxn(cur_sec_aq_spec1%name, &
                                              cur_sec_aq_spec1%eqrxn, &
                                              cur_sec_aq_spec2%eqrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_sec_aq_spec2 => cur_sec_aq_spec2%next      
    enddo
    nullify(cur_sec_aq_spec2)

    ! secondary aqueous species in mineral reactions
    cur_mineral => reaction%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%tstrxn%nspec) exit
          if (StringCompare(cur_sec_aq_spec1%name, &
                              cur_mineral%tstrxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn(cur_sec_aq_spec1%name, &
                                             cur_sec_aq_spec1%eqrxn, &
                                             cur_mineral%tstrxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo

    ! secondary aqueous species in surface complex reactions
    cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_surfcplx_rxn)) exit
      cur_surfcplx2 => cur_surfcplx_rxn%complex_list
      do
        if (.not.associated(cur_surfcplx2)) exit
        
        if (associated(cur_surfcplx2%eqrxn)) then
          ispec = 1
          do
            if (ispec > cur_surfcplx2%eqrxn%nspec) exit
            if (StringCompare(cur_sec_aq_spec1%name, &
                                cur_surfcplx2%eqrxn%spec_name(ispec), &
                                MAXWORDLENGTH)) then
              call BasisSubSpeciesInGasOrSecRxn(cur_sec_aq_spec1%name, &
                                                cur_sec_aq_spec1%eqrxn, &
                                                cur_surfcplx2%eqrxn)
              ispec = 0
            endif
            ispec = ispec + 1
          enddo
        endif
        cur_surfcplx2 => cur_surfcplx2%next
      enddo
      nullify(cur_surfcplx2)
      cur_surfcplx_rxn => cur_surfcplx_rxn%next
    enddo
    nullify(cur_surfcplx_rxn)    
    
    cur_sec_aq_spec1 => cur_sec_aq_spec1%next
  enddo
  

  nullify(cur_aq_spec)
  nullify(cur_pri_aq_spec)
  nullify(cur_sec_aq_spec)
  nullify(cur_sec_aq_spec1)
  nullify(cur_sec_aq_spec2)
  nullify(cur_gas_spec)
  nullify(cur_gas_spec1)
  nullify(cur_gas_spec2)
  nullify(cur_mineral)
  nullify(cur_surfcplx_rxn)
  nullify(cur_surfcplx)
  nullify(cur_surfcplx2)  

  ! align all species in rxns with basis
  cur_sec_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_sec_aq_spec)) exit
    if (.not.associated(cur_sec_aq_spec%eqrxn%spec_ids)) then
      allocate(cur_sec_aq_spec%eqrxn%spec_ids(cur_sec_aq_spec%eqrxn%nspec))
      cur_sec_aq_spec%eqrxn%spec_ids = 0
    endif
    call BasisAlignSpeciesInRxn(ncomp_h2o,new_basis_names, &
                                cur_sec_aq_spec%eqrxn%nspec, &
                                cur_sec_aq_spec%eqrxn%spec_name, &
                                cur_sec_aq_spec%eqrxn%stoich, &
                                cur_sec_aq_spec%eqrxn%spec_ids,option)  
    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo

  cur_gas_spec => reaction%gas_species_list
  do
    if (.not.associated(cur_gas_spec)) exit
    if (.not.associated(cur_gas_spec%eqrxn%spec_ids)) then
      allocate(cur_gas_spec%eqrxn%spec_ids(cur_gas_spec%eqrxn%nspec))
      cur_gas_spec%eqrxn%spec_ids = 0
    endif
    call BasisAlignSpeciesInRxn(ncomp_h2o,new_basis_names, &
                                cur_gas_spec%eqrxn%nspec, &
                                cur_gas_spec%eqrxn%spec_name, &
                                cur_gas_spec%eqrxn%stoich, &
                                cur_gas_spec%eqrxn%spec_ids,option)     
    cur_gas_spec => cur_gas_spec%next
  enddo

#endif

  ! substitute new basis into mineral and surface complexation rxns,
  ! if necessary
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    if (.not.associated(cur_mineral%tstrxn%spec_ids)) then
      allocate(cur_mineral%tstrxn%spec_ids(cur_mineral%tstrxn%nspec))
      cur_mineral%tstrxn%spec_ids = 0
    endif    
    call BasisAlignSpeciesInRxn(ncomp_h2o,new_basis_names, &
                                cur_mineral%tstrxn%nspec, &
                                cur_mineral%tstrxn%spec_name, &
                                cur_mineral%tstrxn%stoich, &
                                cur_mineral%tstrxn%spec_ids,option)     
    cur_mineral => cur_mineral%next
  enddo  

  cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_surfcplx_rxn)) exit
    cur_surfcplx => cur_surfcplx_rxn%complex_list
    do
      if (.not.associated(cur_surfcplx)) exit
      if (.not.associated(cur_surfcplx%eqrxn%spec_ids)) then
        allocate(cur_surfcplx%eqrxn%spec_ids(cur_surfcplx%eqrxn%nspec))
        cur_surfcplx%eqrxn%spec_ids = 0
      endif
      call BasisAlignSpeciesInRxn(ncomp_h2o,new_basis_names, &
                                  cur_surfcplx%eqrxn%nspec, &
                                  cur_surfcplx%eqrxn%spec_name, &
                                  cur_surfcplx%eqrxn%stoich, &
                                  cur_surfcplx%eqrxn%spec_ids,option) 
      cur_surfcplx => cur_surfcplx%next
    enddo
    nullify(cur_surfcplx)
    cur_surfcplx_rxn => cur_surfcplx_rxn%next
  enddo
  nullify(cur_surfcplx_rxn)  

  ! fill reaction arrays, swapping if necessary
  if (associated(reaction%primary_species_names)) &
    deallocate(reaction%primary_species_names)
  allocate(reaction%primary_species_names(reaction%ncomp))
  reaction%primary_species_names = ''
  allocate(reaction%primary_species_print(reaction%ncomp))
  reaction%primary_species_print = PETSC_FALSE
  allocate(reaction%primary_spec_Z(reaction%ncomp))
  reaction%primary_spec_Z = 0.d0
  allocate(reaction%primary_spec_a0(reaction%ncomp))
  reaction%primary_spec_a0 = 0.d0
  
    ! pack in reaction arrays
  cur_pri_aq_spec => reaction%primary_species_list
  ispec = 1
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    reaction%primary_species_names(ispec) = &
      cur_pri_aq_spec%name
    reaction%primary_spec_Z(ispec) = &
      cur_pri_aq_spec%Z
    reaction%primary_spec_a0(ispec) = &
      cur_pri_aq_spec%a0
    reaction%primary_species_print(ispec) = cur_pri_aq_spec%print_me .or. &
                                            reaction%print_all_species
    ispec = ispec + 1
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo
  nullify(cur_pri_aq_spec)
  ispec = -1 ! to catch bugs
  
  ! secondary aqueous complexes
  reaction%neqcmplx = GetSecondarySpeciesCount(reaction)
  
  if (reaction%neqcmplx > 0) then
    allocate(reaction%secondary_species_names(reaction%neqcmplx))
    reaction%secondary_species_names = ''
    allocate(reaction%secondary_species_print(reaction%neqcmplx))
    reaction%secondary_species_print = PETSC_FALSE
    allocate(reaction%eqcmplx_basis_names(reaction%ncomp,reaction%neqcmplx))
    reaction%eqcmplx_basis_names = ''
    allocate(reaction%eqcmplxspecid(0:reaction%ncomp,reaction%neqcmplx))
    reaction%eqcmplxspecid = 0
    allocate(reaction%eqcmplxstoich(0:reaction%ncomp,reaction%neqcmplx))
    reaction%eqcmplxstoich = 0.d0
    allocate(reaction%eqcmplxh2oid(reaction%neqcmplx))
    reaction%eqcmplxh2oid = 0
    allocate(reaction%eqcmplxh2ostoich(reaction%neqcmplx))
    reaction%eqcmplxh2ostoich = 0.d0
    allocate(reaction%eqcmplx_logK(reaction%neqcmplx))
    reaction%eqcmplx_logK = 0.d0
    allocate(reaction%eqcmplx_logKcoef(reaction%num_dbase_temperatures,reaction%neqcmplx))
    reaction%eqcmplx_logKcoef = 0.d0
    allocate(reaction%eqcmplx_Z(reaction%neqcmplx))
    reaction%eqcmplx_Z = 0.d0
    allocate(reaction%eqcmplx_a0(reaction%neqcmplx))
    reaction%eqcmplx_a0 = 0.d0

    ! pack in reaction arrays
    cur_sec_aq_spec => reaction%secondary_species_list
    isec_spec = 1
    do
      if (.not.associated(cur_sec_aq_spec)) exit

      reaction%secondary_species_names(isec_spec) = &
        cur_sec_aq_spec%name
      reaction%secondary_species_print(isec_spec) = cur_sec_aq_spec%print_me .or. &
                                            reaction%print_all_species
      ispec = 0
      do i = 1, cur_sec_aq_spec%eqrxn%nspec
      
!       print *,'database: ',i,cur_sec_aq_spec%eqrxn%spec_name(i)
        
        if (cur_sec_aq_spec%eqrxn%spec_ids(i) /= h2o_id) then
          ispec = ispec + 1
          spec_id = cur_sec_aq_spec%eqrxn%spec_ids(i)
          if (spec_id > h2o_id) spec_id = spec_id - 1
          reaction%eqcmplxspecid(ispec,isec_spec) = spec_id
          reaction%eqcmplx_basis_names(ispec,isec_spec) = &
            cur_sec_aq_spec%eqrxn%spec_name(i)
          reaction%eqcmplxstoich(ispec,isec_spec) = &
            cur_sec_aq_spec%eqrxn%stoich(i)
            
        else ! fill in h2o id and stoich
          reaction%eqcmplxh2oid(isec_spec) = h2o_id
          reaction%eqcmplxh2ostoich(isec_spec) = &
            cur_sec_aq_spec%eqrxn%stoich(i)
        endif
      enddo
      reaction%eqcmplxspecid(0,isec_spec) = ispec
      reaction%eqcmplx_logKcoef(:,isec_spec) = &
        cur_sec_aq_spec%eqrxn%logK
      call Interpolate(temp_high,temp_low,option%reference_temperature, &
                       cur_sec_aq_spec%eqrxn%logK(itemp_high), &
                       cur_sec_aq_spec%eqrxn%logK(itemp_low), &
                       reaction%eqcmplx_logK(isec_spec))
!      reaction%eqcmplx_logK(isec_spec) = cur_sec_aq_spec%eqrxn%logK(option%itemp_ref)
      reaction%eqcmplx_Z(isec_spec) = cur_sec_aq_spec%Z
      reaction%eqcmplx_a0(isec_spec) = cur_sec_aq_spec%a0
  
      isec_spec = isec_spec + 1
      cur_sec_aq_spec => cur_sec_aq_spec%next
    enddo

  endif
  nullify(cur_sec_aq_spec)
  isec_spec = -1 ! to catch bugs

  ! gas complexes
  reaction%ngas = GetGasCount(reaction)
  
  if (reaction%ngas > 0) then
    allocate(reaction%gas_species_names(reaction%ngas))
    reaction%gas_species_names = ''
    allocate(reaction%gas_species_print(reaction%ngas))
    reaction%gas_species_print = PETSC_FALSE
    allocate(reaction%eqgasspecid(0:reaction%ncomp,reaction%ngas))
    reaction%eqgasspecid = 0
    allocate(reaction%eqgasstoich(0:reaction%ncomp,reaction%ngas))
    reaction%eqgasstoich = 0.d0
    allocate(reaction%eqgash2oid(reaction%ngas))
    reaction%eqgash2oid = 0
    allocate(reaction%eqgash2ostoich(reaction%ngas))
    reaction%eqgash2ostoich = 0.d0
    allocate(reaction%eqgas_logK(reaction%ngas))
    reaction%eqgas_logK = 0.d0
#if TEMP_DEPENDENT_LOGK
!Peter change here
    allocate(reaction%eqgas_logKcoef(FIVE_INTEGER,reaction%ngas))
    reaction%eqgas_logKcoef = 0.d0
#else
    allocate(reaction%eqgas_logKcoef(reaction%num_dbase_temperatures, &
                                     reaction%ngas))
    reaction%eqgas_logKcoef = 0.d0
#endif

    ! pack in reaction arrays
    cur_gas_spec => reaction%gas_species_list
    igas_spec = 1
    do
      if (.not.associated(cur_gas_spec)) exit

      reaction%gas_species_names(igas_spec) = &
        cur_gas_spec%name
      reaction%gas_species_print(igas_spec) = cur_gas_spec%print_me .or. &
                                            reaction%print_all_species
      ispec = 0
      do i = 1, cur_gas_spec%eqrxn%nspec
        if (cur_gas_spec%eqrxn%spec_ids(i) /= h2o_id) then
          ispec = ispec + 1
          spec_id = cur_gas_spec%eqrxn%spec_ids(i)
          if (spec_id > h2o_id) spec_id = spec_id - 1
          reaction%eqgasspecid(ispec,igas_spec) = spec_id
          reaction%eqgasstoich(ispec,igas_spec) = &
            cur_gas_spec%eqrxn%stoich(i)
            
        else ! fill in h2o id and stoich
          reaction%eqgash2oid(igas_spec) = h2o_id
          reaction%eqgash2ostoich(igas_spec) = &
            cur_gas_spec%eqrxn%stoich(i)
        endif
      enddo
      reaction%eqgasspecid(0,igas_spec) = ispec
#if TEMP_DEPENDENT_LOGK
!Peter change here
      call ReactionFitLogKCoef(reaction%eqgas_logKcoef(:,igas_spec),cur_gas_spec%eqrxn%logK, &
                               option,reaction)
      call ReactionInitializeLogK(reaction%eqgas_logKcoef(:,igas_spec), &
                                  cur_gas_spec%eqrxn%logK, &
                                  reaction%eqgas_logK(igas_spec), &
                                  option,reaction)
#else
      reaction%eqgas_logKcoef(:,igas_spec) = &
        cur_gas_spec%eqrxn%logK
      call Interpolate(temp_high,temp_low,option%reference_temperature, &
                       cur_gas_spec%eqrxn%logK(itemp_high), &
                       cur_gas_spec%eqrxn%logK(itemp_low), &
                       reaction%eqgas_logK(igas_spec))
!      reaction%eqgas_logK(igas_spec) = cur_gas_spec%eqrxn%logK(option%itemp_ref)
#endif      
  
      igas_spec = igas_spec + 1
      cur_gas_spec => cur_gas_spec%next
    enddo

  endif
  nullify(cur_gas_spec)
  igas_spec = -1 ! to catch bugs

  ! minerals
!  reaction%nmnrl = GetMineralCount(reaction)
  reaction%nkinmnrl = GetKineticMineralCount(reaction)

  if (reaction%nmnrl > 0) then
    allocate(reaction%mineral_names(reaction%nmnrl))
    reaction%mineral_names = ''
    allocate(reaction%mnrlspecid(0:reaction%ncomp,reaction%nmnrl))
    reaction%mnrlspecid = 0
    allocate(reaction%mnrlstoich(reaction%ncomp,reaction%nmnrl))
    reaction%mnrlstoich = 0.d0
    allocate(reaction%mnrlh2oid(reaction%nmnrl))
    reaction%mnrlh2oid = 0
    allocate(reaction%mnrlh2ostoich(reaction%nmnrl))
    reaction%mnrlh2ostoich = 0.d0
    allocate(reaction%mnrl_logK(reaction%nmnrl))
    reaction%mnrl_logK = 0.d0
    allocate(reaction%mnrl_logKcoef(reaction%num_dbase_temperatures, &
                                    reaction%nmnrl))
    reaction%mnrl_logKcoef = 0.d0

    allocate(reaction%kinmnrl_names(reaction%nkinmnrl))
    reaction%kinmnrl_names = ''
    allocate(reaction%kinmnrl_print(reaction%nkinmnrl))
    reaction%kinmnrl_print = PETSC_FALSE
    allocate(reaction%kinmnrlspecid(0:reaction%ncomp,reaction%nkinmnrl))
    reaction%kinmnrlspecid = 0
    allocate(reaction%kinmnrlstoich(reaction%ncomp,reaction%nkinmnrl))
    reaction%kinmnrlstoich = 0.d0
    allocate(reaction%kinmnrlh2oid(reaction%nkinmnrl))
    reaction%kinmnrlh2oid = 0
    allocate(reaction%kinmnrlh2ostoich(reaction%nkinmnrl))
    reaction%kinmnrlh2ostoich = 0.d0
    allocate(reaction%kinmnrl_logK(reaction%nkinmnrl))
    reaction%kinmnrl_logK = 0.d0
    allocate(reaction%kinmnrl_logKcoef(reaction%num_dbase_temperatures, &
                                       reaction%nkinmnrl))
    reaction%kinmnrl_logKcoef = 0.d0
    allocate(reaction%kinmnrl_rate(1,reaction%nkinmnrl))
    reaction%kinmnrl_rate = 0.d0
    allocate(reaction%kinmnrl_molar_vol(reaction%nkinmnrl))
    reaction%kinmnrl_molar_vol = 0.d0
    allocate(reaction%kinmnrl_num_prefactors(reaction%nkinmnrl))
    reaction%kinmnrl_num_prefactors = 0
    
    cur_mineral => reaction%mineral_list
    imnrl = 1
    ikinmnrl = 1
    do
      if (.not.associated(cur_mineral)) exit

      reaction%mineral_names(imnrl) = cur_mineral%name
      ispec = 0
      do i = 1, cur_mineral%tstrxn%nspec
        if (cur_mineral%tstrxn%spec_ids(i) /= h2o_id) then
          ispec = ispec + 1
          spec_id = cur_mineral%tstrxn%spec_ids(i)
          if (spec_id > h2o_id) spec_id = spec_id - 1
          reaction%mnrlspecid(ispec,imnrl) = spec_id
          reaction%mnrlstoich(ispec,imnrl) = &
            cur_mineral%tstrxn%stoich(i)
            
        else ! fill in h2o id and stoich
          reaction%mnrlh2oid(imnrl) = h2o_id
          reaction%mnrlh2ostoich(imnrl) = &
            cur_mineral%tstrxn%stoich(i)
        endif
      enddo
      reaction%mnrlspecid(0,imnrl) = ispec
      reaction%mnrl_logKcoef(:,imnrl) = &
        cur_mineral%tstrxn%logK
      call Interpolate(temp_high,temp_low,option%reference_temperature, &
                       cur_mineral%tstrxn%logK(itemp_high), &
                       cur_mineral%tstrxn%logK(itemp_low), &
                       reaction%mnrl_logK(imnrl))
!      reaction%mnrl_logK(imnrl) = cur_mineral%tstrxn%logK(option%itemp_ref)
  
      if (cur_mineral%itype == MINERAL_KINETIC) then
        reaction%kinmnrl_names(ikinmnrl) = reaction%mineral_names(imnrl)
        reaction%kinmnrl_print(ikinmnrl) = cur_mineral%print_me .or. &
                                            reaction%print_all_species
        reaction%kinmnrlspecid(:,ikinmnrl) = reaction%mnrlspecid(:,imnrl)
        reaction%kinmnrlstoich(:,ikinmnrl) = reaction%mnrlstoich(:,imnrl)
        reaction%kinmnrlh2oid(ikinmnrl) = reaction%mnrlh2oid(imnrl)
        reaction%kinmnrlh2ostoich(ikinmnrl) = reaction%mnrlh2ostoich(imnrl)
        reaction%kinmnrl_logK(ikinmnrl) = reaction%mnrl_logK(imnrl)
        reaction%kinmnrl_logKcoef(:,ikinmnrl) = reaction%mnrl_logKcoef(:,imnrl)
        reaction%kinmnrl_rate(1,ikinmnrl) = cur_mineral%tstrxn%rate
        reaction%kinmnrl_molar_vol(ikinmnrl) = cur_mineral%molar_volume
        ikinmnrl = ikinmnrl + 1
      endif

      cur_mineral => cur_mineral%next
      imnrl = imnrl + 1
    enddo
  endif
  
  if (reaction%neqsurfcmplx > 0) then
  
    ! determine max # complexes for a given site
    icount = 0
    cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_surfcplx_rxn)) exit
      isurfcplx = 0
      cur_surfcplx => cur_surfcplx_rxn%complex_list
      do
        if (.not.associated(cur_surfcplx)) exit
        isurfcplx = isurfcplx + 1
        cur_surfcplx => cur_surfcplx%next
      enddo
      if (isurfcplx > icount) icount = isurfcplx
      cur_surfcplx_rxn => cur_surfcplx_rxn%next
    enddo
    nullify(cur_surfcplx_rxn)  

    allocate(reaction%eqsurfcmplx_rxn_to_mineral(reaction%neqsurfcmplxrxn))
    reaction%eqsurfcmplx_rxn_to_mineral = 0
    allocate(reaction%eqsurfcmplx_rxn_to_complex(0:icount, &
                                                 reaction%neqsurfcmplxrxn))
    reaction%eqsurfcmplx_rxn_to_complex = 0
    allocate(reaction%surface_site_names(reaction%neqsurfcmplxrxn))
    reaction%surface_site_names = ''
    allocate(reaction%surface_site_print(reaction%neqsurfcmplxrxn))
    reaction%surface_site_print = PETSC_FALSE
    allocate(reaction%eqsurfcmplx_rxn_site_density(reaction%neqsurfcmplxrxn))
    reaction%eqsurfcmplx_rxn_site_density = 0.d0
    allocate(reaction%eqsurfcmplx_rxn_stoich_flag(reaction%neqsurfcmplxrxn))
    reaction%eqsurfcmplx_rxn_stoich_flag = PETSC_FALSE
    allocate(reaction%surface_complex_names(reaction%neqsurfcmplx))
    reaction%surface_complex_names = ''
    allocate(reaction%surface_complex_print(reaction%neqsurfcmplx))
    reaction%surface_complex_print = PETSC_FALSE
    allocate(reaction%eqsurfcmplxspecid(0:reaction%ncomp,reaction%neqsurfcmplx))
    reaction%eqsurfcmplxspecid = 0
    allocate(reaction%eqsurfcmplxstoich(reaction%ncomp,reaction%neqsurfcmplx))
    reaction%eqsurfcmplxstoich = 0.d0
    allocate(reaction%eqsurfcmplxh2oid(reaction%neqsurfcmplx))
    reaction%eqsurfcmplxh2oid = 0
    allocate(reaction%eqsurfcmplxh2ostoich(reaction%neqsurfcmplx))
    reaction%eqsurfcmplxh2ostoich = 0.d0
    allocate(reaction%eqsurfcmplx_free_site_id(reaction%neqsurfcmplx))
    reaction%eqsurfcmplx_free_site_id = 0
    allocate(reaction%eqsurfcmplx_free_site_stoich(reaction%neqsurfcmplx))
    reaction%eqsurfcmplx_free_site_stoich = 0.d0
    allocate(reaction%eqsurfcmplx_mineral_id(reaction%neqsurfcmplx))
    reaction%eqsurfcmplx_mineral_id = 0
    allocate(reaction%eqsurfcmplx_logK(reaction%neqsurfcmplx))
    reaction%eqsurfcmplx_logK = 0.d0
    allocate(reaction%eqsurfcmplx_logKcoef(reaction%num_dbase_temperatures, &
                                           reaction%neqsurfcmplx))
    reaction%eqsurfcmplx_logKcoef = 0.d0
    allocate(reaction%eqsurfcmplx_Z(reaction%neqsurfcmplx))
    reaction%eqsurfcmplx_Z = 0.d0

    isurfcplx = 0
    irxn = 0
    cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_surfcplx_rxn)) exit
      
      irxn = irxn + 1
      reaction%surface_site_names(irxn) = cur_surfcplx_rxn%free_site_name
      reaction%surface_site_print(irxn) = cur_surfcplx_rxn%free_site_print_me .or. &
                                            reaction%print_all_species
      reaction%eqsurfcmplx_rxn_to_mineral(irxn) = &
        GetMineralIDFromName(reaction,cur_surfcplx_rxn%mineral_name)
      reaction%eqsurfcmplx_rxn_site_density(irxn) = cur_surfcplx_rxn%site_density
            
      cur_surfcplx => cur_surfcplx_rxn%complex_list
      do
        if (.not.associated(cur_surfcplx)) exit
        
        isurfcplx = isurfcplx + 1
        
        ! set up integer pointers from site to complexes
        ! increment count for site
        reaction%eqsurfcmplx_rxn_to_complex(0,irxn) = &
          reaction%eqsurfcmplx_rxn_to_complex(0,irxn) + 1
        reaction%eqsurfcmplx_rxn_to_complex( &
          reaction%eqsurfcmplx_rxn_to_complex(0,irxn),irxn) = isurfcplx 
        
        reaction%surface_complex_names(isurfcplx) = cur_surfcplx%name
        reaction%surface_complex_print(isurfcplx) = cur_surfcplx%print_me .or. &
                                            reaction%print_all_species
        reaction%eqsurfcmplx_free_site_id(isurfcplx) = &
          cur_surfcplx_rxn%free_site_id
        reaction%eqsurfcmplx_free_site_stoich(isurfcplx) =  &
          cur_surfcplx%free_site_stoich
          
        if (cur_surfcplx%free_site_stoich > 1.d0) then
          reaction%eqsurfcmplx_rxn_stoich_flag(irxn) = PETSC_TRUE
        endif
 
        ispec = 0
        do i = 1, cur_surfcplx%eqrxn%nspec
          if (cur_surfcplx%eqrxn%spec_ids(i) /= h2o_id) then
            ispec = ispec + 1
            spec_id = cur_surfcplx%eqrxn%spec_ids(i)
            if (spec_id > h2o_id) spec_id = spec_id - 1
            reaction%eqsurfcmplxspecid(ispec,isurfcplx) = spec_id
            reaction%eqsurfcmplxstoich(ispec,isurfcplx) = &
              cur_surfcplx%eqrxn%stoich(i)
            
          else ! fill in h2o id and stoich
            reaction%eqsurfcmplxh2oid(isurfcplx) = h2o_id
            reaction%eqsurfcmplxh2ostoich(isurfcplx) = &
              cur_surfcplx%eqrxn%stoich(i)
          endif
        enddo
        reaction%eqsurfcmplxspecid(0,isurfcplx) = ispec
        reaction%eqsurfcmplx_logKcoef(:,isurfcplx) = &
          cur_surfcplx%eqrxn%logK
        call Interpolate(temp_high,temp_low,option%reference_temperature, &
                         cur_surfcplx%eqrxn%logK(itemp_high), &
                         cur_surfcplx%eqrxn%logK(itemp_low), &
                         reaction%eqsurfcmplx_logK(isurfcplx))
        !reaction%eqsurfcmplx_logK(isurfcplx) = cur_surfcplx%eqrxn%logK(option%itemp_ref)
        reaction%eqsurfcmplx_Z(isurfcplx) = cur_surfcplx%Z

        cur_surfcplx => cur_surfcplx%next
      enddo
      nullify(cur_surfcplx)
      cur_surfcplx_rxn => cur_surfcplx_rxn%next
    enddo
    nullify(cur_surfcplx_rxn)  
  
  endif

  if (reaction%neqionxrxn > 0) then

    ! determine max # cations for a given ionx exchange rxn
    icount = 0
    cur_ionx_rxn => reaction%ion_exchange_rxn_list
    do
      if (.not.associated(cur_ionx_rxn)) exit
      ication = 0
      cur_cation => cur_ionx_rxn%cation_list
      do
        if (.not.associated(cur_cation)) exit
        ication = ication + 1
        cur_cation => cur_cation%next
      enddo
      if (ication > icount) icount = ication
      cur_ionx_rxn => cur_ionx_rxn%next
    enddo
    nullify(cur_ionx_rxn)
    
    allocate(reaction%eqionx_rxn_cationid(0:icount,reaction%neqionxrxn))
    reaction%eqionx_rxn_cationid = 0
    allocate(reaction%eqionx_rxn_Z_flag(reaction%neqionxrxn))
    reaction%eqionx_rxn_Z_flag = PETSC_FALSE
    allocate(reaction%eqionx_rxn_cation_X_offset(reaction%neqionxrxn))
    reaction%eqionx_rxn_cation_X_offset = 0
    allocate(reaction%eqionx_rxn_CEC(reaction%neqionxrxn))
    reaction%eqionx_rxn_CEC = 0.d0
    allocate(reaction%eqionx_rxn_k(icount,reaction%neqionxrxn))
    reaction%eqionx_rxn_k = 0.d0

    irxn = 0
    icount = 0
    cur_ionx_rxn => reaction%ion_exchange_rxn_list
    do
      if (.not.associated(cur_ionx_rxn)) exit
      irxn = irxn + 1
      ication = 0
      reaction%eqionx_rxn_CEC(irxn) = cur_ionx_rxn%CEC
        ! compute the offset to the first cation in rxn
      reaction%eqionx_rxn_cation_X_offset(irxn) = icount
        
      cur_cation => cur_ionx_rxn%cation_list
      do
        if (.not.associated(cur_cation)) exit
        ication = ication + 1
        icount = icount + 1
        reaction%eqionx_rxn_k(ication,irxn) = cur_cation%k

        found = PETSC_FALSE
        do i = 1, reaction%ncomp
          if (StringCompare(cur_cation%name, &
                              reaction%primary_species_names(i), &
                              MAXWORDLENGTH)) then
            reaction%eqionx_rxn_cationid(ication,irxn) = i
            found = PETSC_TRUE        
          endif
        enddo
        if (.not.found) then
          option%io_buffer = 'Cation ' // trim(cur_cation%name) // &
                   ' in ion exchange reaction' // &
                   ' not found in swapped basis.'
          call printErrMsg(option)     
        endif
        cur_cation => cur_cation%next
      enddo
      reaction%eqionx_rxn_cationid(0,irxn) = ication
      ! Find any Zi /= Zj for all species i, j
      found = PETSC_FALSE
      do i = 1, reaction%eqionx_rxn_cationid(0,irxn)
        do j = 1, reaction%eqionx_rxn_cationid(0,irxn)
          if (abs(reaction%primary_spec_Z(reaction%eqionx_rxn_cationid(i,irxn))- &
                  reaction%primary_spec_Z(reaction%eqionx_rxn_cationid(j,irxn))) > &
              0.1d0) then
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (found) exit
      enddo
      reaction%eqionx_rxn_Z_flag(irxn) = found
      cur_ionx_rxn => cur_ionx_rxn%next
    enddo
    nullify(cur_ionx_rxn)

  endif

  call BasisPrint(reaction,'Final Basis',option)
  
  ! locate specific species
  do ispec = 1, reaction%ncomp
    if (reaction%h_ion_id == 0) then
      word = 'H+'
      if (StringCompare(reaction%primary_species_names(ispec), &
                          word,MAXWORDLENGTH)) then
        reaction%h_ion_id = ispec
      endif
    endif
  enddo
  
  do ispec = 1, reaction%neqcmplx
    if (reaction%h_ion_id == 0) then
      word = 'H+'
      if (StringCompare(reaction%secondary_species_names(ispec), &
                          word,MAXWORDLENGTH)) then
        reaction%h_ion_id = -ispec
      endif
    endif
  enddo

  do ispec = 1, reaction%ngas
    if (reaction%o2_gas_id == 0) then
      word = 'O2(g)'
      if (StringCompare(reaction%gas_species_names(ispec), &
                          word,MAXWORDLENGTH)) then
        reaction%o2_gas_id = ispec
      endif
    endif
    if (reaction%co2_gas_id == 0) then
      word = 'CO2(g)'
      if (StringCompare(reaction%gas_species_names(ispec), &
                          word,MAXWORDLENGTH)) then
        reaction%co2_gas_id = ispec
      endif
      word = 'CO2(g)*'
      if (StringCompare(reaction%gas_species_names(ispec), &
                          word,MAXWORDLENGTH)) then
        reaction%co2_gas_id = ispec
      endif

    endif

  enddo
  
90 format(80('-'))
100 format(/,2x,i4,2x,a)
110 format(100(/,14x,3(a20,2x)))

  if (OptionPrintToFile(option)) then
    write(option%fid_out,90)
    write(option%fid_out,100) reaction%ncomp, 'Primary Species'
    write(option%fid_out,110) (reaction%primary_species_names(i),i=1,reaction%ncomp)
    
    write(option%fid_out,100) reaction%neqcmplx, 'Secondary Complex Species'
    write(option%fid_out,110) (reaction%secondary_species_names(i),i=1,reaction%neqcmplx)
    
    write(option%fid_out,100) reaction%ngas, 'Gas Species'
    write(option%fid_out,110) (reaction%gas_species_names(i),i=1,reaction%ngas)
    
    write(option%fid_out,100) reaction%nmnrl, 'Reference Minerals'
    write(option%fid_out,110) (reaction%mineral_names(i),i=1,reaction%nmnrl)
    
    write(option%fid_out,100) reaction%nkinmnrl, 'Kinetic Mineral Reactions'
    write(option%fid_out,110) (reaction%kinmnrl_names(i),i=1,reaction%nkinmnrl)
    
    write(option%fid_out,100) reaction%neqsurfcmplxrxn, 'Surface Complexation Reactions'
    write(option%fid_out,110) (reaction%surface_site_names(i),i=1,reaction%neqsurfcmplxrxn)
    write(option%fid_out,100) reaction%neqsurfcmplx, 'Surface Complexes'
    write(option%fid_out,110) (reaction%surface_complex_names(i),i=1,reaction%neqsurfcmplx)
    
    write(option%fid_out,100) reaction%neqionxrxn, 'Ion Exchange Reactions'
    write(option%fid_out,100) reaction%neqionxcation, 'Ion Exchange Cations'
    write(option%fid_out,90)
  endif
  
  if (allocated(new_basis)) deallocate(new_basis)
  if (allocated(old_basis)) deallocate(old_basis)
  if (allocated(transformation)) deallocate(transformation)
  if (allocated(stoich_prev)) deallocate(stoich_prev)
  if (allocated(stoich_new)) deallocate(stoich_new)
  if (allocated(logKvector)) deallocate(logKvector)
  if (allocated(indices)) deallocate(indices)

  if (allocated(new_basis_names)) deallocate(new_basis_names)
  if (allocated(old_basis_names)) deallocate(old_basis_names)
  
end subroutine BasisInit

! ************************************************************************** !
!
! GetSpeciesBasisID: Reduces redundant coding above
! author: Glenn Hammond
! date: 12/02/08
!
! ************************************************************************** !
function GetSpeciesBasisID(reaction,option,ncomp_h2o,reaction_name, &
                           species_name, &
                           pri_names,sec_names,gas_names)

  use Option_module
  use String_module

  implicit none

  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscInt :: ncomp_h2o
  character(len=MAXWORDLENGTH) :: reaction_name
  character(len=MAXWORDLENGTH) :: species_name
  character(len=MAXWORDLENGTH) :: pri_names(:)
  character(len=MAXWORDLENGTH) :: sec_names(:)
  character(len=MAXWORDLENGTH) :: gas_names(:)

  PetscInt :: GetSpeciesBasisID
  PetscInt :: i

  GetSpeciesBasisID = 0
  do i=1,ncomp_h2o
    if (StringCompare(species_name, &
                        pri_names(i),MAXWORDLENGTH)) then
      GetSpeciesBasisID = i
      return
    endif
  enddo
  ! secondary aqueous and gas species denoted by negative id
  do i=1,reaction%neqcmplx
    if (StringCompare(species_name, &
                        sec_names(i),MAXWORDLENGTH)) then
      GetSpeciesBasisID = -i
      return
    endif
  enddo
  do i=1,reaction%ngas
    if (StringCompare(species_name, &
                        gas_names(i),MAXWORDLENGTH)) then
      GetSpeciesBasisID = -(reaction%neqcmplx+i)
      return
    endif
  enddo
  
  option%io_buffer = 'Species ' // &
           trim(species_name) // &
           ' listed in reaction for ' // &
           trim(reaction_name) // &
           ' not found among primary, secondary, or gas species.'
  call printErrMsg(option)

end function GetSpeciesBasisID

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
                                  rxn_stoich,rxn_species_ids,option)

  use Option_module
  use String_module
  
  implicit none
  
  PetscInt :: num_basis_species
  character(len=MAXWORDLENGTH) :: basis_names(num_basis_species)
  PetscInt :: num_rxn_species
  character(len=MAXWORDLENGTH) :: rxn_species_names(num_rxn_species)
  PetscReal :: rxn_stoich(num_rxn_species)
  PetscInt :: rxn_species_ids(num_rxn_species)
  type(option_type) :: option

  PetscInt :: i_rxn_species
  PetscInt :: i_basis_species
  PetscReal :: stoich_new(num_basis_species)
  PetscTruth :: found
  
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
               ' not found in basis (BasisAlignSpeciesInRxn)'
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
subroutine BasisSubSpeciesInGasOrSecRxn(name1,eqrxn1,eqrxn2)

  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name1
  type(equilibrium_rxn_type) :: eqrxn1
  type(equilibrium_rxn_type) :: eqrxn2
  
  PetscInt :: i, j, tempcount, prevcount
  character(len=MAXWORDLENGTH) :: tempnames(20)
  PetscReal :: tempstoich(20)
  PetscReal :: scale
  PetscTruth :: found

  tempnames = ''
  tempstoich = 0.d0

  ! load species in reaction other than species 1 into new arrays
  scale = 1.d0
  tempcount = 0
  do i=1,eqrxn2%nspec
    if (.not.StringCompare(name1, &
                             eqrxn2%spec_name(i), &
                             MAXWORDLENGTH)) then
      tempcount = tempcount + 1
      tempnames(tempcount) = eqrxn2%spec_name(i)
      tempstoich(tempcount) = eqrxn2%stoich(i)
    else
      scale = eqrxn2%stoich(i)
    endif
  enddo
  
  ! search for duplicate species and add stoichs or add new species
  ! if not duplicated
  do j=1,eqrxn1%nspec
    found = PETSC_FALSE
    do i=1,tempcount
      if (StringCompare(tempnames(i), &
                          eqrxn1%spec_name(j), &
                          MAXWORDLENGTH)) then
        tempstoich(i) = tempstoich(i) + scale*eqrxn1%stoich(j)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      tempcount = tempcount + 1
      tempnames(tempcount) = eqrxn1%spec_name(j)
      tempstoich(tempcount) = scale*eqrxn1%stoich(j)
    endif
  enddo
  
  ! deallocate arrays
  deallocate(eqrxn2%spec_name)
  deallocate(eqrxn2%stoich)
  
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
  allocate(eqrxn2%spec_name(tempcount))
  allocate(eqrxn2%stoich(tempcount))
  
  ! fill arrays in eqrxn
  eqrxn2%nspec = tempcount
  do i=1,tempcount
    eqrxn2%spec_name(i) = tempnames(i)
    eqrxn2%stoich(i) = tempstoich(i)
  enddo
  
  ! adjust the equilibrium coefficient
  eqrxn2%logK = eqrxn2%logK + scale*eqrxn1%logK

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
subroutine BasisSubSpeciesInMineralRxn(name,eqrxn,tstrxn)

  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(equilibrium_rxn_type) :: eqrxn
  type(transition_state_rxn_type) :: tstrxn
  
  PetscInt :: i, j, tempcount, prevcount
  character(len=MAXWORDLENGTH) :: tempnames(20)
  PetscReal :: tempstoich(20)
  PetscReal :: scale
  PetscTruth :: found

  tempnames = ''
  tempstoich = 0.d0
  
  ! load species in reaction other than species 1 into new arrays
  scale = 1.d0
  tempcount = 0
  do i=1,tstrxn%nspec
    if (.not.StringCompare(name, &
                             tstrxn%spec_name(i), &
                             MAXWORDLENGTH)) then
      tempcount = tempcount + 1
      tempnames(tempcount) = tstrxn%spec_name(i)
      tempstoich(tempcount) = tstrxn%stoich(i)
    else
      scale = tstrxn%stoich(i)
    endif
  enddo
  
  ! search for duplicate species and add stoichs or add new species
  ! if not duplicated
  do j=1,eqrxn%nspec
    found = PETSC_FALSE
    do i=1,tstrxn%nspec
      if (StringCompare(tempnames(i), &
                          eqrxn%spec_name(j), &
                          MAXWORDLENGTH)) then
        tempstoich(i) = tempstoich(i) + scale*eqrxn%stoich(j)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      tempcount = tempcount + 1
      tempnames(tempcount) = eqrxn%spec_name(j)
      tempstoich(tempcount) = scale*eqrxn%stoich(j)
    endif
  enddo
  
  ! deallocate arrays
  deallocate(tstrxn%spec_name)
  deallocate(tstrxn%stoich)
  
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
  allocate(tstrxn%spec_name(tempcount))
  allocate(tstrxn%stoich(tempcount))
  
  ! fill arrays in eqrxn
  tstrxn%nspec = tempcount
  do i=1,tempcount
    tstrxn%spec_name(i) = tempnames(i)
    tstrxn%stoich(i) = tempstoich(i)
  enddo
  
  ! adjust the equilibrium coefficient
  tstrxn%logK = tstrxn%logK + scale*eqrxn%logK

end subroutine BasisSubSpeciesInMineralRxn

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
  type(surface_complexation_rxn_type), pointer :: cur_surfcplx_rxn
  type(surface_complex_type), pointer :: cur_surfcplx
  type(ion_exchange_rxn_type), pointer :: cur_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cur_cation

  PetscInt :: ispec, itemp

100 format(a)
110 format(a,f9.4)
120 format(a,f6.2,2x,a)
130 format(a,100f9.4)
140 format(a,f6.2)
150 format(a,es11.4)

  if (OptionPrintToFile(option)) then
    write(option%fid_out,*)
    write(option%fid_out,*) '! *************************************************' // &
                    '************************* !'
    write(option%fid_out,*)
    write(option%fid_out,*) trim(title)
    write(option%fid_out,*)

    write(option%fid_out,*) 'Primary Components:'
    cur_aq_spec => reaction%primary_species_list
    do
      if (.not.associated(cur_aq_spec)) exit
      write(option%fid_out,100) '  ' // trim(cur_aq_spec%name)
      write(option%fid_out,140) '    Charge: ', cur_aq_spec%Z
      write(option%fid_out,110) '    Molar Weight: ', cur_aq_spec%molar_weight
      write(option%fid_out,110) '    Debye-Huckel a0: ', cur_aq_spec%a0
      if (associated(cur_aq_spec%eqrxn)) then
        write(option%fid_out,100) '    Equilibrium Aqueous Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_aq_spec%name
        do ispec = 1, cur_aq_spec%eqrxn%nspec
          write(option%fid_out,120) '      ', cur_aq_spec%eqrxn%stoich(ispec), &
                          cur_aq_spec%eqrxn%spec_name(ispec)
        enddo
        write(option%fid_out,130) '      logK:', (cur_aq_spec%eqrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(option%fid_out,*)
      cur_aq_spec => cur_aq_spec%next
    enddo
    
    cur_aq_spec => reaction%secondary_species_list
    if (associated(cur_aq_spec)) then
      write(option%fid_out,*)
      write(option%fid_out,*) 'Secondary Components:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Secondary Components: None'
    endif
    do
      if (.not.associated(cur_aq_spec)) exit
      write(option%fid_out,100) '  ' // trim(cur_aq_spec%name)
      write(option%fid_out,140) '    Charge: ', cur_aq_spec%Z
      write(option%fid_out,110) '    Molar Weight: ', cur_aq_spec%molar_weight
      write(option%fid_out,110) '    Debye-Huckel a0: ', cur_aq_spec%a0
      if (associated(cur_aq_spec%eqrxn)) then
        write(option%fid_out,100) '    Equilibrium Aqueous Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_aq_spec%name
        do ispec = 1, cur_aq_spec%eqrxn%nspec
          write(option%fid_out,120) '      ', cur_aq_spec%eqrxn%stoich(ispec), &
                          cur_aq_spec%eqrxn%spec_name(ispec)
        enddo
        write(option%fid_out,130) '      logK:', (cur_aq_spec%eqrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(option%fid_out,*)
      cur_aq_spec => cur_aq_spec%next
    enddo
    
    cur_gas_spec => reaction%gas_species_list
    if (associated(cur_gas_spec)) then
      write(option%fid_out,*)
      write(option%fid_out,*) 'Gas Components:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Gas Components: None'
    endif
    do
      if (.not.associated(cur_gas_spec)) exit
      write(option%fid_out,100) '  ' // trim(cur_gas_spec%name)
      write(option%fid_out,110) '    Molar Weight: ', cur_gas_spec%molar_weight
      if (associated(cur_gas_spec%eqrxn)) then
        write(option%fid_out,100) '    Gas Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_gas_spec%name
        do ispec = 1, cur_gas_spec%eqrxn%nspec
          write(option%fid_out,120) '      ', cur_gas_spec%eqrxn%stoich(ispec), &
                          cur_gas_spec%eqrxn%spec_name(ispec)
        enddo
        write(option%fid_out,130) '      logK:', (cur_gas_spec%eqrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(option%fid_out,*)
      cur_gas_spec => cur_gas_spec%next
    enddo
    
    cur_mineral => reaction%mineral_list
    if (associated(cur_mineral)) then    
      write(option%fid_out,*)
      write(option%fid_out,*) 'Minerals:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Minerals: None'    
    endif
    do
      if (.not.associated(cur_mineral)) exit
      write(option%fid_out,100) '  ' // trim(cur_mineral%name)
      write(option%fid_out,110) '    Molar Weight: ', cur_mineral%molar_weight
      write(option%fid_out,110) '    Molar Volume: ', cur_mineral%molar_volume
      if (associated(cur_mineral%tstrxn)) then
        write(option%fid_out,100) '    Mineral Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_mineral%name
        do ispec = 1, cur_mineral%tstrxn%nspec
          write(option%fid_out,120) '      ', cur_mineral%tstrxn%stoich(ispec), &
                          cur_mineral%tstrxn%spec_name(ispec)
        enddo
        write(option%fid_out,130) '      logK:', (cur_mineral%tstrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
      endif
      write(option%fid_out,*)
      cur_mineral => cur_mineral%next
    enddo
    
    cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
    if (associated(cur_surfcplx)) then
      write(option%fid_out,*)
      write(option%fid_out,*) 'Surface Complexation Reactions:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Surface Complexation Reactions: None'
    endif
    do
      if (.not.associated(cur_surfcplx_rxn)) exit
      cur_surfcplx => cur_surfcplx_rxn%complex_list
      if (associated(cur_surfcplx)) then
        write(option%fid_out,*)
        write(option%fid_out,*) '  Surface Complexes:'
      else
        write(option%fid_out,*)
        write(option%fid_out,*) '  Surface Complexes: None'
      endif
      do
        if (.not.associated(cur_surfcplx)) exit
        write(option%fid_out,100) '  ' // trim(cur_surfcplx%name)
        write(option%fid_out,140) '    Charge: ', cur_surfcplx%Z
        write(option%fid_out,100) '    Surface Complex Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_surfcplx%name
        write(option%fid_out,120) '      ', cur_surfcplx%free_site_stoich, &
          cur_surfcplx_rxn%free_site_name
        do ispec = 1, cur_surfcplx%eqrxn%nspec
          write(option%fid_out,120) '      ', cur_surfcplx%eqrxn%stoich(ispec), &
                          cur_surfcplx%eqrxn%spec_name(ispec)
        enddo
        write(option%fid_out,130) '      logK:', (cur_surfcplx%eqrxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
        write(option%fid_out,*)
        cur_surfcplx => cur_surfcplx%next
      enddo
      cur_surfcplx_rxn => cur_surfcplx_rxn%next
    enddo
    
    cur_ionx_rxn => reaction%ion_exchange_rxn_list
    if (associated(cur_ionx_rxn)) then
      write(option%fid_out,*)
      write(option%fid_out,*) 'Ion Exchange Reactions:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Ion Exchange Reactions: None'
    endif
    do
      if (.not.associated(cur_ionx_rxn)) exit
      write(option%fid_out,*) '  Mineral: ', trim(cur_ionx_rxn%mineral_name)
      write(option%fid_out,150) '      CEC: ', cur_ionx_rxn%CEC
      cur_cation => cur_ionx_rxn%cation_list
      if (associated(cur_cation)) then
        write(option%fid_out,*) '  Cations:'
      else
        write(option%fid_out,*) '  Cations: None'
      endif
      do
        if (.not.associated(cur_cation)) exit
        write(option%fid_out,150) '      ' // trim(cur_cation%name), cur_cation%k
        cur_cation => cur_cation%next
      enddo
      write(option%fid_out,*)
      cur_ionx_rxn => cur_ionx_rxn%next
    enddo
    
    write(option%fid_out,*)
    write(option%fid_out,*) '! *************************************************' // &
                    '************************* !'
    write(option%fid_out,*)
  endif

end subroutine BasisPrint

end module Database_module
