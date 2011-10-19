module Reaction_module

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: ReactionCreate, &
            ReactionRead, &
            ReactionReadMineralKinetics, &
            ReactionReadOutput, &
            ReactionReadRedoxSpecies, &
            RTotal, &
            RTotalSorb, &
            CO2AqActCoeff, &
            RActivityCoefficients, &
            RReaction, &
            RReactionDerivative, &
            ReactionProcessConstraint, &
            ReactionEquilibrateConstraint, &
            ReactionPrintConstraint, &
            ReactionFitLogKCoef, &
            ReactionInitializeLogK, &
            ReactionComputeKd, &
            RAccumulationSorb, &
            RAccumulationSorbDerivative, &
            RJumpStartKineticSorption, &
            RAge, &
            RReact, &
            RTAuxVarCompute, &
            RTAccumulation, &
            RTAccumulationDerivative, &
            RTPrintAuxVar, &
            RMineralSaturationIndex, &
            DoubleLayer

contains

! ************************************************************************** !
!
! ReactionRead: Reads chemical species
! author: Glenn Hammond
! date: 05/02/08
!
! ************************************************************************** !
subroutine ReactionRead(reaction,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
  type(aq_species_type), pointer :: species, prev_species
  type(gas_species_type), pointer :: gas, prev_gas
  type(mineral_type), pointer :: mineral, prev_mineral
  type(colloid_type), pointer :: colloid, prev_colloid
  type(surface_complex_type), pointer :: srfcplx, cur_srfcplx, prev_srfcplx
  type(surface_complex_type), pointer :: rate_list, cur_srfcplx_rate, prev_srfcplx_rate
  type(surface_complexation_rxn_type), pointer :: srfcplx_rxn, &
                                                  prev_srfcplx_rxn
  type(ion_exchange_rxn_type), pointer :: ionx_rxn, prev_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cation, prev_cation
  type(general_rxn_type), pointer :: general_rxn, prev_general_rxn
  type(kd_rxn_type), pointer :: kd_rxn, prev_kd_rxn
  PetscInt :: i, tempint
  PetscReal :: tempreal
  PetscInt :: srfcplx_count
  PetscInt :: temp_srfcplx_count
  PetscBool :: found

  nullify(prev_species)
  nullify(prev_gas)
  nullify(prev_mineral)
  nullify(prev_colloid)
  nullify(prev_srfcplx_rxn)
  nullify(cur_srfcplx)
  nullify(prev_srfcplx)
  nullify(rate_list)
  nullify(cur_srfcplx_rate)
  nullify(prev_srfcplx_rate)
  nullify(prev_ionx_rxn)
  nullify(prev_cation)
  nullify(prev_general_rxn)
  nullify(prev_kd_rxn)
  
  srfcplx_count = 0
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY')
    call StringToUpper(word)
    
    select case(trim(word))
    
      case('PRIMARY_SPECIES')
        nullify(prev_species)
        do
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          reaction%naqcomp = reaction%naqcomp + 1
          
          species => AqueousSpeciesCreate()
          call InputReadWord(input,option,species%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,PRIMARY_SPECIES')    
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
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          reaction%neqcplx = reaction%neqcplx + 1
          
          species => AqueousSpeciesCreate()
          call InputReadWord(input,option,species%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,PRIMARY_SPECIES')    
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
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          reaction%ngas = reaction%ngas + 1
          
          gas => GasSpeciesCreate()
          call InputReadWord(input,option,gas%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,GAS_SPECIES')    
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
      case('GENERAL_REACTION')

        reaction%ngeneral_rxn = reaction%ngeneral_rxn + 1
        
        general_rxn => GeneralRxnCreate()
        do 
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,GENERAL_REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('REACTION')
              ! remainder of string should be the reaction equation
              general_rxn%reaction = trim(adjustl(input%buf))
              ! set flag for error message
              if (len_trim(general_rxn%reaction) < 2) input%ierr = 1
              call InputErrorMsg(input,option,'reaction', &
                                 'CHEMISTRY,GENERAL_REACTION,REACTION') 
! For now, the reactants are negative stoich, products positive in reaction equation - geh
#if 0
            case('FORWARD_SPECIES')
              nullify(prev_species)
              do
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                
                species => AqueousSpeciesCreate()
                call InputReadWord(input,option,species%name,PETSC_TRUE)  
                call InputErrorMsg(input,option,'keyword','CHEMISTRY, &
                                   GENERAL_REACTION,FORWARD_SPECIES')    
                if (.not.associated(general_rxn%forward_species_list)) then
                  general_rxn%forward_species_list => species
                  species%id = 1
                endif
                if (associated(prev_species)) then
                  prev_species%next => species
                  species%id = prev_species%id + 1
                endif
                prev_species => species
                nullify(species)
              enddo
              nullify(prev_species)
            case('BACKWARD_SPECIES')
              nullify(prev_species)
              do
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                
                species => AqueousSpeciesCreate()
                call InputReadWord(input,option,species%name,PETSC_TRUE)  
                call InputErrorMsg(input,option,'keyword','CHEMISTRY, &
                                   GENERAL_REACTION,BACKWARD_SPECIES')    
                if (.not.associated(general_rxn%backward_species_list)) then
                  general_rxn%forward_species_list => species
                  species%id = 1
                endif
                if (associated(prev_species)) then
                  prev_species%next => species
                  species%id = prev_species%id + 1
                endif
                prev_species => species
                nullify(species)
              enddo
              nullify(prev_species)
#endif                
            case('FORWARD_RATE')
              call InputReadDouble(input,option,general_rxn%forward_rate)  
              call InputDefaultMsg(input,option, &
                                   'CHEMISTRY,GENERAL_REACTION,FORWARD_RATE') 
            case('BACKWARD_RATE')
              call InputReadDouble(input,option,general_rxn%backward_rate)  
              call InputDefaultMsg(input,option, &
                                   'CHEMISTRY,GENERAL_REACTION,BACKWARD_RATE') 
          end select
        enddo   
        if (.not.associated(reaction%general_rxn_list)) then
          reaction%general_rxn_list => general_rxn
          general_rxn%id = 1
        endif
        if (associated(prev_general_rxn)) then
          prev_general_rxn%next => general_rxn
          general_rxn%id = prev_general_rxn%id + 1
        endif
        prev_general_rxn => general_rxn
        nullify(general_rxn)

      case('MINERALS')
        nullify(prev_mineral)
        do
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          reaction%nmnrl = reaction%nmnrl + 1
          
          mineral => MineralCreate()
          call InputReadWord(input,option,mineral%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,MINERALS')    
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
      case('MINERAL_KINETICS')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,name,PETSC_TRUE)
          call InputErrorMsg(input,option,name,'CHEMISTRY,MINERAL_KINETICS')
          do
            call InputReadFlotranString(input,option)
            call InputReadStringErrorMsg(input,option,card)
            if (InputCheckExit(input,option)) exit
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword', &
                                   'CHEMISTRY,MINERAL_KINETICS')
            call StringToUpper(word)
            select case(word)
              case('PREFACTOR')
                do 
                  call InputReadFlotranString(input,option)
                  call InputReadStringErrorMsg(input,option,card)
                  if (InputCheckExit(input,option)) exit
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'keyword', &
                                     'CHEMISTRY,MINERAL_KINETICS,PREFACTOR')
                  call StringToUpper(word)
                  select case(word)
                    case('PREFACTOR_SPECIES')
                      call InputSkipToEnd(input,option,word)
                  end select
                enddo
            end select
          enddo
        enddo
      case('COLLOIDS')
        nullify(prev_colloid)
        do
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          reaction%ncoll = reaction%ncoll + 1
          
          colloid => ColloidCreate()
          call InputReadWord(input,option,colloid%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,COLLOIDS')    
          call InputReadDouble(input,option,colloid%mobile_fraction)  
          call InputDefaultMsg(input,option,'CHEMISTRY,COLLOIDS,MOBILE_FRACTION')          
          if (.not.associated(reaction%colloid_list)) then
            reaction%colloid_list => colloid
            colloid%id = 1
          endif
          if (associated(prev_colloid)) then
            prev_colloid%next => colloid
            colloid%id = prev_colloid%id + 1
          endif
          prev_colloid => colloid
          nullify(colloid)
        enddo
      case('SORPTION')
        nullify(prev_srfcplx_rxn)
        do
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,SORPTION')
          call StringToUpper(word)   

          select case(trim(word))

            case('ISOTHERM_REACTIONS')
              option%io_buffer = 'Isotherm reactions currently calculated as ' // &
                'a function of free-ion, not totals.  Contact Glenn!'
              call printErrMsg(option)
              do
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                reaction%neqkdrxn = reaction%neqkdrxn + 1

                kd_rxn => KDRxnCreate()
                ! first string is species name
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'species name', &
                                   'CHEMISTRY,ISOTHERM_REACTIONS')
                kd_rxn%species_name = trim(word)
                do 
                  call InputReadFlotranString(input,option)
                  if (InputError(input)) exit
                  if (InputCheckExit(input,option)) exit

                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'keyword', &
                                     'CHEMISTRY,ISOTHERM_REACTIONS')
                  call StringToUpper(word)
                  
                  ! default type is linear
                  kd_rxn%itype = SORPTION_LINEAR
                  select case(trim(word))
                    case('TYPE')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'type', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                      select case(word)
                        case('LINEAR')
                          kd_rxn%itype = SORPTION_LINEAR
                        case('LANGMUIR')
                          kd_rxn%itype = SORPTION_LANGMUIR
                        case('FREUNDLICH')
                          kd_rxn%itype = SORPTION_FREUNDLICH
                      end select
                    case('DISTRIBUTION_COEFFICIENT','KD')
                      call InputReadDouble(input,option,kd_rxn%Kd)
                      call InputErrorMsg(input,option,'DISTRIBUTION_COEFFICIENT', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                    case('LANGMUIR_B')
                      call InputReadDouble(input,option,kd_rxn%Langmuir_B)
                      call InputErrorMsg(input,option,'Langmuir_B', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                      kd_rxn%itype = SORPTION_LANGMUIR
                    case('FREUNDLICH_N')
                      call InputReadDouble(input,option,kd_rxn%Freundlich_N)
                      call InputErrorMsg(input,option,'Freundlich_N', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                      kd_rxn%itype = SORPTION_FREUNDLICH
                  end select
                enddo
                ! add to list
                if (.not.associated(reaction%kd_rxn_list)) then
                  reaction%kd_rxn_list => kd_rxn
                  kd_rxn%id = 1
                endif
                if (associated(prev_kd_rxn)) then
                  prev_kd_rxn%next => kd_rxn
                  kd_rxn%id = prev_kd_rxn%id + 1
                endif
                prev_kd_rxn => kd_rxn
                nullify(kd_rxn)
                
              enddo
            
            case('SURFACE_COMPLEXATION_RXN')

              ! initialization of temporary variables
              temp_srfcplx_count = 0
              srfcplx_rxn => SurfaceComplexationRxnCreate()
              srfcplx_rxn%itype = SRFCMPLX_RXN_EQUILIBRIUM
              do
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword', &
                                   'CHEMISTRY,SURFACE_COMPLEXATION_RXN')
                call StringToUpper(word)
                
                select case(trim(word))
                  case('EQUILIBRIUM')
                    srfcplx_rxn%itype = SRFCMPLX_RXN_EQUILIBRIUM
                  case('MULTIRATE_KINETIC')
                    srfcplx_rxn%itype = SRFCMPLX_RXN_MULTIRATE_KINETIC
                  case('KINETIC')
                    srfcplx_rxn%itype = SRFCMPLX_RXN_KINETIC
                  case('COMPLEX_KINETICS')
                    nullify(prev_srfcplx)
                    do
                      call InputReadFlotranString(input,option)
                      if (InputError(input)) exit
                      if (InputCheckExit(input,option)) exit
                      
                      srfcplx => SurfaceComplexCreate()
                      call InputReadWord(input,option,srfcplx%name,PETSC_TRUE)
                      call InputErrorMsg(input,option,'keyword', &
                        'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE')
                        
                      do
                        call InputReadFlotranString(input,option)
                        call InputReadStringErrorMsg(input,option,card)
                        if (InputCheckExit(input,option)) exit
                        call InputReadWord(input,option,word,PETSC_TRUE)
                        call InputErrorMsg(input,option,'word', &
                               'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE') 
                        select case(trim(word))
                          case('FORWARD_RATE_CONSTANT')
                            call InputReadDouble(input,option,srfcplx%forward_rate)
                            call InputErrorMsg(input,option,'forward_rate', &
                                   'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE')
                          case('BACKWARD_RATE_CONSTANT')
                            call InputReadDouble(input,option,srfcplx%backward_rate)
                            call InputErrorMsg(input,option,'backward_rate', &
                                   'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE')
                          case default
                            option%io_buffer = 'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE keyword: ' // &
                                               trim(word) // ' not recognized'
                            call printErrMsg(option)
                        end select
                      enddo
                                      
                      if (.not.associated(rate_list)) then
                        rate_list => srfcplx
                      endif
                      if (associated(prev_srfcplx)) then
                        prev_srfcplx%next => srfcplx
                      endif
                      prev_srfcplx => srfcplx
                      nullify(srfcplx)
                    enddo
                    nullify(prev_srfcplx)
                  case('RATE','RATES') 
                    srfcplx_rxn%itype = SRFCMPLX_RXN_MULTIRATE_KINETIC
                    string = 'RATES inside SURFACE_COMPLEXATION_RXN'
                    call UtilityReadArray(reaction%kinmr_rate,NEG_ONE_INTEGER,string,input, &
                                          option) 
                  case('SITE_FRACTION') 
                    string = 'SITE_FRACTION inside SURFACE_COMPLEXATION_RXN'
                    call UtilityReadArray(reaction%kinmr_frac,NEG_ONE_INTEGER,string,input, &
                                          option) 
                  case('MULTIRATE_SCALE_FACTOR')
                    call InputReadDouble(input,option,reaction%kinmr_scale_factor)
                    call InputErrorMsg(input,option,'keyword', &
                      'CHEMISTRY,SURFACE_COMPLEXATION_RXN,MULTIRATE_SCALE_FACTOR')
                  case('MINERAL')
                    call InputReadWord(input,option,srfcplx_rxn%mineral_name, &
                      PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword', &
                      'CHEMISTRY,SURFACE_COMPLEXATION_RXN,MINERAL_NAME')
                  case('COLLOID')
                    call InputReadWord(input,option,srfcplx_rxn%colloid_name, &
                      PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword', &
                      'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COLLOID_NAME')
                  case('SITE')
                    call InputReadWord(input,option,srfcplx_rxn%free_site_name, &
                      PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword', &
                      'CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_NAME')
                    ! site density in mol/m^3 bulk
                    call InputReadDouble(input,option,srfcplx_rxn%site_density)
                    call InputErrorMsg(input,option,'keyword', &
                      'CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_DENSITY')                   
                  case('COMPLEXES')
                    nullify(prev_srfcplx)
                    do
                      call InputReadFlotranString(input,option)
                      if (InputError(input)) exit
                      if (InputCheckExit(input,option)) exit
                      
                      temp_srfcplx_count = temp_srfcplx_count + 1
                      srfcplx_count = srfcplx_count + 1
                      srfcplx => SurfaceComplexCreate()
                      srfcplx%id = srfcplx_count
                      call InputReadWord(input,option,srfcplx%name,PETSC_TRUE)
                      call InputErrorMsg(input,option,'keyword', &
                        'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_NAME')
                
                      if (.not.associated(srfcplx_rxn%complex_list)) then
                        srfcplx_rxn%complex_list => srfcplx
                      endif
                      if (associated(prev_srfcplx)) then
                        prev_srfcplx%next => srfcplx
                      endif
                      prev_srfcplx => srfcplx
                      nullify(srfcplx)
                    enddo
                    nullify(prev_srfcplx)
                  case default
                    option%io_buffer = 'CHEMISTRY, SURFACE_COMPLEXATION_RXN keyword: '// &
                                     trim(word)//' not recognized'
                    call printErrMsg(option)
                end select

              enddo
              if (.not.associated(reaction%surface_complexation_rxn_list)) then
                reaction%surface_complexation_rxn_list => srfcplx_rxn
                srfcplx_rxn%id = 1
              endif
              if (associated(prev_srfcplx_rxn)) then
                prev_srfcplx_rxn%next => srfcplx_rxn
                srfcplx_rxn%id = prev_srfcplx_rxn%id + 1
              endif
              prev_srfcplx_rxn => srfcplx_rxn

              select case(srfcplx_rxn%itype)
                ! default (NULL) to EQUILIBRIUM
                case(SRFCMPLX_RXN_NULL,SRFCMPLX_RXN_EQUILIBRIUM, &
                     SRFCMPLX_RXN_MULTIRATE_KINETIC)
                  reaction%neqsrfcplx = reaction%neqsrfcplx + &
                    temp_srfcplx_count
                  reaction%neqsrfcplxrxn = reaction%neqsrfcplxrxn + 1
                case(SRFCMPLX_RXN_KINETIC)
                  ! match up rates with their corresponding surface complex
                  cur_srfcplx => srfcplx_rxn%complex_list
                  do
                    if (.not.associated(cur_srfcplx)) exit
                    found = PETSC_FALSE
                    nullify(prev_srfcplx_rate)
                    cur_srfcplx_rate => rate_list
                    do
                      if (.not.associated(cur_srfcplx_rate)) exit
                      ! check for same name
                      if (StringCompare(cur_srfcplx_rate%name, &
                                        cur_srfcplx%name, &
                                        MAXWORDLENGTH)) then
                        ! set rates
                        cur_srfcplx%forward_rate = cur_srfcplx_rate%forward_rate
                        cur_srfcplx%backward_rate = cur_srfcplx_rate%backward_rate
                        ! remove srfcplx_rate from list of rates
                        if (associated(prev_srfcplx_rate)) then
                          prev_srfcplx_rate%next => cur_srfcplx_rate%next
                        else
                          rate_list => cur_srfcplx_rate%next
                        endif
                        ! destroy the object
                        call SurfaceComplexDestroy(cur_srfcplx_rate)
                        found = PETSC_TRUE
                        exit
                      endif
                      prev_srfcplx_rate => cur_srfcplx_rate
                      cur_srfcplx_rate => cur_srfcplx_rate%next
                    enddo
                    if (.not.found) then
                      option%io_buffer = 'Rates for surface complex ' // &
                        trim(cur_srfcplx%name) // ' not found in kinetic rate list'
                      call printErrMsg(option)
                    endif
                    cur_srfcplx => cur_srfcplx%next
                  enddo
                  ! check to ensure that rates are matched
                  if (associated(rate_list)) then
                    option%io_buffer = '# of rates is greater than # of surface complexes'
                    call printErrMsg(option)
                  endif
                  nullify(cur_srfcplx)
                  nullify(prev_srfcplx)
                  nullify(rate_list)
                  nullify(cur_srfcplx_rate)
                  nullify(prev_srfcplx_rate)                  
                  reaction%nkinsrfcplx = reaction%nkinsrfcplx + temp_srfcplx_count
                  reaction%nkinsrfcplxrxn = reaction%nkinsrfcplxrxn + 1
              end select
              srfcplx_rxn%free_site_id = srfcplx_rxn%id
              
              nullify(srfcplx_rxn)

            case('ION_EXCHANGE_RXN')
            
              ionx_rxn => IonExchangeRxnCreate()
              do
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword','CHEMISTRY,ION_EXCHANGE_RXN')
                call StringToUpper(word)
                
                select case(trim(word))
                  case('MINERAL')
                    call InputReadWord(input,option,ionx_rxn%mineral_name,PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword', &
                      'CHEMISTRY,ION_EXCHANGE_RXN,MINERAL_NAME')
                  case('CEC')
                    call InputReadDouble(input,option,ionx_rxn%CEC)
                    call InputErrorMsg(input,option,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,CEC')                   
                  case('CATIONS')
                    nullify(prev_cation)
                    do
                      call InputReadFlotranString(input,option)
                      if (InputError(input)) exit
                      if (InputCheckExit(input,option)) exit
                      
                      cation => IonExchangeCationCreate()
                      reaction%neqionxcation = reaction%neqionxcation + 1
                      call InputReadWord(input,option,cation%name,PETSC_TRUE)
                      call InputErrorMsg(input,option,'keyword', &
                        'CHEMISTRY,ION_EXCHANGE_RXN,CATION_NAME')
                      call InputReadDouble(input,option,cation%k)
                      call InputErrorMsg(input,option,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,K')                   
    
                      if (.not.associated(ionx_rxn%cation_list)) then
                        ionx_rxn%cation_list => cation
                      endif
                      if (associated(prev_cation)) then
                        prev_cation%next => cation
                      endif
                      prev_cation => cation
                      nullify(cation)
                    enddo
                  case default
                    option%io_buffer = 'CHEMISTRY, ION_EXCHANGE_RXN keyword: '// &
                                     trim(word)//' not recognized'
                    call printErrMsg(option)
                end select
              enddo
              if (.not.associated(reaction%ion_exchange_rxn_list)) then
                reaction%ion_exchange_rxn_list => ionx_rxn
                ionx_rxn%id = 1
              endif
              if (associated(prev_ionx_rxn)) then
                prev_ionx_rxn%next => ionx_rxn
                ionx_rxn%id = prev_ionx_rxn%id + 1
              endif
              prev_ionx_rxn => ionx_rxn

              reaction%neqionxrxn = ionx_rxn%id

              nullify(ionx_rxn)          
            case('JUMPSTART_KINETIC_SORPTION')
              option%jumpstart_kinetic_sorption = PETSC_TRUE
              option%no_restart_kinetic_sorption = PETSC_TRUE
            case('NO_CHECKPOINT_KINETIC_SORPTION')
              option%no_checkpoint_kinetic_sorption = PETSC_TRUE
            case('NO_RESTART_KINETIC_SORPTION')
              option%no_restart_kinetic_sorption = PETSC_TRUE
          end select
        enddo
      case('DATABASE')
        call InputReadNChars(input,option,reaction%database_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)  
        call InputErrorMsg(input,option,'keyword', &
                           'CHEMISTRY,DATABASE FILENAME')  
      case('LOG_FORMULATION')
        reaction%use_log_formulation = PETSC_TRUE        
      case('NO_CHECK_UPDATE')
        reaction%check_update = PETSC_FALSE       
      case('NO_RESTART_MINERAL_VOL_FRAC')
        option%no_restart_mineral_vol_frac = PETSC_TRUE
      case('NO_CHECKPOINT_ACT_COEFS')
        reaction%checkpoint_activity_coefs = PETSC_FALSE
      case('ACTIVITY_COEFFICIENTS')
        reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_LAG        
        reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_TIMESTEP        
        do 
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr /= 0) exit
          select case(trim(word))
            case('OFF')
              reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
            case('LAG')
              reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_LAG    
            case('NEWTON')
              reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_NEWTON       
            case('TIMESTEP')
              reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_TIMESTEP
            case('NEWTON_ITERATION')
              reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_NEWTON_ITER
            case default
              option%io_buffer = 'CHEMISTRY,ACTIVITY_COEFFICIENTS keyword: ' &
                                 //trim(word)//' not recognized'
              call printErrMsg(option)
          end select
        enddo
      case('NO_BDOT')
        reaction%act_coef_use_bdot = PETSC_FALSE
      case('CHUNK_SIZE')
        call InputReadInt(input,option,option%chunk_size)
        call InputErrorMsg(input,option,'chunk_size','CHEMISTRY')
      case('NUM_THREADS')
        call InputReadInt(input,option,option%num_threads)
        call InputErrorMsg(input,option,'num_thread','CHEMISTRY')
      case('UPDATE_POROSITY')
        option%update_porosity = PETSC_TRUE
      case('UPDATE_TORTUOSITY')
        option%update_tortuosity = PETSC_TRUE
      case('UPDATE_PERMEABILITY')
        option%update_permeability = PETSC_TRUE
      case('UPDATE_MINERAL_SURFACE_AREA')
        option%update_mineral_surface_area = PETSC_TRUE
      case('MOLAL','MOLALITY')
        reaction%initialize_with_molality = PETSC_TRUE
      case('ACTIVITY_H2O','ACTIVITY_WATER')
        reaction%use_activity_h2o = PETSC_TRUE
      case('REDOX_SPECIES')
        call InputSkipToEnd(input,option,word)
      case('OUTPUT')
        call InputSkipToEnd(input,option,word)
      case('MAX_DLNC')
        call InputReadDouble(input,option,reaction%max_dlnC)
        call InputErrorMsg(input,option,trim(word),'CHEMISTRY')
      case('OPERATOR_SPLIT','OPERATOR_SPLITTING')
        option%reactive_transport_coupling = OPERATOR_SPLIT    
      case('MAX_RELATIVE_CHANGE_TOLERANCE','REACTION_TOLERANCE')
        call InputReadDouble(input,option,reaction%max_relative_change_tolerance)
        call InputErrorMsg(input,option,'maximum relative change tolerance','CHEMISTRY')
      case('MAX_RESIDUAL_TOLERANCE')
        call InputReadDouble(input,option,reaction%max_residual_tolerance)
        call InputErrorMsg(input,option,'maximum residual  tolerance','CHEMISTRY')
      case default
        option%io_buffer = 'CHEMISTRY keyword: '//trim(word)//' not recognized'
        call printErrMsg(option)
    end select
  enddo
  
  reaction%neqsorb = reaction%neqsrfcplxrxn + reaction%neqionxrxn + &
                     reaction%neqkdrxn

  if (reaction%print_free_conc_type == 0) then
    if (reaction%initialize_with_molality) then
      reaction%print_free_conc_type = PRIMARY_MOLALITY
    else
      reaction%print_free_conc_type = PRIMARY_MOLARITY
    endif
  endif
  if (reaction%print_tot_conc_type == 0) then
    if (reaction%initialize_with_molality) then
      reaction%print_tot_conc_type = TOTAL_MOLALITY
    else
      reaction%print_tot_conc_type = TOTAL_MOLARITY
    endif
  endif
  if (reaction%print_secondary_conc_type == 0) then
    if (reaction%initialize_with_molality) then
      reaction%print_secondary_conc_type = SECONDARY_MOLALITY
    else
      reaction%print_secondary_conc_type = SECONDARY_MOLARITY
    endif
  endif
  if (reaction%neqcplx + reaction%neqsorb + reaction%nmnrl + &
      reaction%ngeneral_rxn > 0) then
    reaction%use_full_geochemistry = PETSC_TRUE
  endif
 
  ! check to ensure that rates for multirate surface complexation are aligned
  ! with surface fractions
  if (associated(reaction%kinmr_rate)) then
    reaction%kinmr_nrate = size(reaction%kinmr_rate)
    if (reaction%kinmr_nrate > 0) then
      if (size(reaction%kinmr_rate) /= size(reaction%kinmr_frac)) then
        write(word,*) size(reaction%kinmr_rate)
        write(string,*) size(reaction%kinmr_frac)
        option%io_buffer = 'Number of kinetic rates (' // &
          trim(adjustl(word)) // &
          ') does not match the number of surface fractions (' // &
          trim(adjustl(string)) // ').'
        call printErrMsg(option)
      endif
      tempreal = 0.d0
      do i = 1, size(reaction%kinmr_frac)
        tempreal = tempreal + reaction%kinmr_frac(i)
        reaction%kinmr_rate(i) = reaction%kinmr_rate(i) * reaction%kinmr_scale_factor
      enddo
    
      if (dabs(1.d0 - tempreal) > 1.d-6) then
        write(string,*) tempreal
        option%io_buffer = 'The sum of the surface fractions for ' // &
          'multirate kinetic sorption does not add up to 1.d0 (' // &
          trim(adjustl(string)) // '.'
        call printErrMsg(option)
      endif
    endif
  endif
  
  if (len_trim(reaction%database_filename) < 2) &
    reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
  
end subroutine ReactionRead

! ************************************************************************** !
!
! ReactionReadMineralKinetics: Reads mineral kinetics
! author: Glenn Hammond
! date: 10/16/08
!
! ************************************************************************** !
subroutine ReactionReadMineralKinetics(reaction,input,option)

  use Input_module
  use String_module  
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
  
  type(mineral_type), pointer :: cur_mineral
  type(transition_state_rxn_type), pointer :: tstrxn, cur_tstrxn
  type(transition_state_prefactor_type), pointer :: prefactor, &
                                                    cur_prefactor
  type(ts_prefactor_species_type), pointer :: prefactor_species, &
                                              cur_prefactor_species
  PetscBool :: found
  PetscInt :: imnrl,icount

  cur_mineral => reaction%mineral_list
  do 
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -1*abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo
  
  input%ierr = 0
  icount = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,name,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MINERAL_KINETICS')
    
    cur_mineral => reaction%mineral_list
    found = PETSC_FALSE
    do 
      if (.not.associated(cur_mineral)) exit
      if (StringCompare(cur_mineral%name,name,MAXWORDLENGTH)) then
        found = PETSC_TRUE
        cur_mineral%itype = MINERAL_KINETIC
        tstrxn => TransitionStateTheoryRxnCreate()
        ! initialize to -999 to ensure that it is set
        tstrxn%rate = -999.d0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          error_string = 'CHEMISTRY,MINERAL_KINETICS'
          call InputErrorMsg(input,option,'word',error_string) 
          select case(trim(word))
            case('RATE_CONSTANT')
!             read rate constant
              call InputReadDouble(input,option,tstrxn%rate)
              call InputErrorMsg(input,option,'rate',error_string)
            case('ACTIVATION_ENERGY')
!             read activation energy for Arrhenius law
              call InputReadDouble(input,option,tstrxn%activation_energy)
              call InputErrorMsg(input,option,'activation',error_string)
            case('AFFINITY_THRESHOLD')
!             read affinity threshold for precipitation
              call InputReadDouble(input,option,tstrxn%affinity_threshold)
              call InputErrorMsg(input,option,'threshold',error_string)
            case('RATE_LIMITER')
!             read rate limiter for precipitation
              call InputReadDouble(input,option,tstrxn%rate_limiter)
              call InputErrorMsg(input,option,'rate_limiter',error_string)
            case('IRREVERSIBLE')
!             read flag for irreversible reaction
              tstrxn%irreversible = 1
              call InputErrorMsg(input,option,'irreversible',error_string)
            case('PREFACTOR')
              error_string = 'CHEMISTRY,MINERAL_KINETICS,PREFACTOR'
              prefactor => TransitionStatePrefactorCreate()
              ! Initialize to -999.d0 to check later whether they were set
              prefactor%rate = -999.d0
              prefactor%activation_energy = -999.d0
              do
                call InputReadFlotranString(input,option)
                call InputReadStringErrorMsg(input,option,card)
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'word',error_string) 
                select case(trim(word))
                  case('RATE_CONSTANT')
    !             read rate constant
                  call InputReadDouble(input,option,prefactor%rate)
                  call InputErrorMsg(input,option,'rate',error_string)
                  case('ACTIVATION_ENERGY')
      !             read activation energy for Arrhenius law
                    call InputReadDouble(input,option,prefactor%activation_energy)
                    call InputErrorMsg(input,option,'activation',error_string)
                  case('PREFACTOR_SPECIES')
                    error_string = 'CHEMISTRY,MINERAL_KINETICS,PREFACTOR,SPECIES'
                    prefactor_species => TSPrefactorSpeciesCreate()
                    call InputReadWord(input,option,prefactor_species%name,PETSC_TRUE)
                    call InputErrorMsg(input,option,'name',error_string)
                    do
                      call InputReadFlotranString(input,option)
                      call InputReadStringErrorMsg(input,option,card)
                      if (InputCheckExit(input,option)) exit
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'keyword',error_string) 
                      select case(trim(word))
                        case('ALPHA')
                          call InputReadDouble(input,option, &
                                               prefactor_species%alpha)
                          call InputErrorMsg(input,option,'alpha',error_string)
                        case('BETA')
                          call InputReadDouble(input,option, &
                                               prefactor_species%beta)
                          call InputErrorMsg(input,option,'beta',error_string)
                        case('ATTENUATION_COEF')
                          call InputReadDouble(input,option, &
                                            prefactor_species%attenuation_coef)
                          call InputErrorMsg(input,option, &
                                             'attenuation coefficient', &
                                             error_string)
                        case default
                          option%io_buffer = 'CHEMISTRY,MINERAL_KINETICS,PREFACTOR, ' // &
                                             'SPECIES keyword: ' // &
                                             trim(word) // ' not recognized'
                          call printErrMsg(option)
                      end select
                    enddo
                    ! add prefactor species
                    if (.not.associated(prefactor%species)) then
                      prefactor%species => prefactor_species
                    else ! append to end of list
                      cur_prefactor_species => prefactor%species
                      do
                        if (.not.associated(cur_prefactor_species%next)) then
                          cur_prefactor_species%next => prefactor_species
                          exit
                        else
                          cur_prefactor_species => cur_prefactor_species%next
                        endif
                      enddo
                    endif                    
                    error_string = 'CHEMISTRY,MINERAL_KINETICS,PREFACTOR'
                  case default
                    option%io_buffer = 'CHEMISTRY,MINERAL_KINETICS,PREFACTOR ' // &
                                 'keyword: ' // trim(word) // ' not recognized'
                    call printErrMsg(option)
                end select
              enddo
              ! add prefactor
              if (.not.associated(tstrxn%prefactor)) then
                tstrxn%prefactor => prefactor
              else ! append to end of list
                cur_prefactor => tstrxn%prefactor
                do
                  if (.not.associated(cur_prefactor%next)) then
                    cur_prefactor%next => prefactor
                    exit
                  else
                    cur_prefactor => cur_prefactor%next
                  endif
                enddo
              endif
              error_string = 'CHEMISTRY,MINERAL_KINETICS'
            case default
              option%io_buffer = 'CHEMISTRY,MINERAL_KINETICS keyword: ' // &
                                 trim(word) // ' not recognized'
              call printErrMsg(option)
          end select
        enddo
        ! Loop over prefactors and set kinetic rates and activation energies
        ! equal to the "outer" values if zero.  
        cur_prefactor => tstrxn%prefactor
        do
          if (.not.associated(cur_prefactor)) exit
          ! if not initialized
          if (dabs(cur_prefactor%rate - (-999.d0)) < 1.d-40) then
            cur_prefactor%rate = tstrxn%rate
            if (dabs(cur_prefactor%rate - (-999.d0)) < 1.d-40) then
              option%io_buffer = 'Both outer and inner prefactor rate ' // &
                'constants uninitialized for kinetic mineral ' // &
                cur_mineral%name // '.'
              call printErrMsg(option)
            endif
          endif
          if (dabs(cur_prefactor%activation_energy - (-999.d0)) < 1.d-40) then
            cur_prefactor%activation_energy = tstrxn%activation_energy
          endif
          cur_prefactor => cur_prefactor%next
        enddo
        ! add tst rxn
        if (.not.associated(cur_mineral%tstrxn)) then
          cur_mineral%tstrxn => tstrxn
        else ! append to end of list
          cur_tstrxn => cur_mineral%tstrxn
          do
            if (.not.associated(cur_tstrxn%next)) then
              cur_tstrxn%next => tstrxn
              exit
            else
              cur_tstrxn => cur_tstrxn%next
            endif
          enddo
        endif
        cur_mineral%id = abs(cur_mineral%id)
        reaction%nkinmnrl = reaction%nkinmnrl + 1
        exit
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (.not.found) then
      option%io_buffer = 'Mineral "' // trim(name) // '" specified under ' // &
        'CHEMISTRY,MINERAL_KINETICS not found in list of available minerals.'
      call printErrMsg(option)
    endif
  enddo
  
  ! allocate kinetic mineral names
  if (reaction%nkinmnrl > 0) then
    if (associated(reaction%kinmnrl_names)) deallocate(reaction%kinmnrl_names)
    allocate(reaction%kinmnrl_names(reaction%nkinmnrl))
    reaction%kinmnrl_names(reaction%nkinmnrl) = ''
  endif
 
  cur_mineral => reaction%mineral_list
  imnrl = 0
  do 
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%id < 0 .and. &
        cur_mineral%itype == MINERAL_KINETIC) then
      option%io_buffer = 'No rate provided in input file for mineral: ' // &
               trim(cur_mineral%name) // '.'
      call printErrMsg(option)
    endif
    if (associated(cur_mineral%tstrxn)) then
      imnrl = imnrl + 1
      reaction%kinmnrl_names(imnrl) = cur_mineral%name
    endif
    cur_mineral => cur_mineral%next
  enddo
  
  cur_mineral => reaction%mineral_list
  do 
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo

end subroutine ReactionReadMineralKinetics

! ************************************************************************** !
!
! ReactionReadRedoxSpecies: Reads names of mineral species and sets flag
! author: Glenn Hammond
! date: 04/01/11
!
! ************************************************************************** !
subroutine ReactionReadRedoxSpecies(reaction,input,option)

  use Input_module
  use String_module  
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: name
  
  type(aq_species_type), pointer :: cur_species

  input%ierr = 0
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,name,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REDOX_SPECIES')
    
    cur_species => reaction%primary_species_list
    do 
      if (.not.associated(cur_species)) exit
      if (StringCompare(cur_species%name,name,MAXWORDLENGTH)) then
        cur_species%is_redox = PETSC_TRUE
        exit
      endif
      cur_species => cur_species%next
    enddo
    
    if (.not.associated(cur_species)) then
      option%io_buffer = 'Redox species "' // trim(name) // &
        '" not found among primary species.'
      call printErrMsg(option)
    endif
  enddo
  
end subroutine ReactionReadRedoxSpecies

! ************************************************************************** !
!
! ReactionProcessConstraint: Initializes constraints based on primary
!                            species in system
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine ReactionProcessConstraint(reaction,constraint_name, &
                                     aq_species_constraint, &
                                     mineral_constraint, &
                                     srfcplx_constraint, &
                                     colloid_constraint, &
                                     option)
  use Option_module
  use Input_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(srfcplx_constraint_type), pointer :: srfcplx_constraint
  type(colloid_constraint_type), pointer :: colloid_constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: icomp, jcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: icoll, jcoll
  PetscInt :: igas
  PetscInt :: isrfcplx, jsrfcplx
  PetscReal :: constraint_conc(reaction%naqcomp)
  PetscInt :: constraint_type(reaction%naqcomp)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(reaction%naqcomp)
  character(len=MAXWORDLENGTH) :: constraint_mnrl_name(reaction%nkinmnrl)
  character(len=MAXWORDLENGTH) :: constraint_srfcplx_name(reaction%nkinsrfcplx)
  character(len=MAXWORDLENGTH) :: constraint_colloid_name(reaction%ncoll)
  PetscInt :: constraint_id(reaction%naqcomp)
  PetscBool :: external_dataset(reaction%naqcomp)
  
  constraint_id = 0
  constraint_aux_string = ''
  constraint_type = 0
  constraint_conc = 0.d0
  
  ! aqueous species
  do icomp = 1, reaction%naqcomp
    found = PETSC_FALSE
    do jcomp = 1, reaction%naqcomp
      if (StringCompare(aq_species_constraint%names(icomp), &
                        reaction%primary_species_names(jcomp), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
               'Species ' // trim(aq_species_constraint%names(icomp)) // &
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among primary species.'
      call printErrMsg(option)
    else
      constraint_type(jcomp) = aq_species_constraint%constraint_type(icomp)
      constraint_aux_string(jcomp) = aq_species_constraint%constraint_aux_string(icomp)
      constraint_conc(jcomp) = aq_species_constraint%constraint_conc(icomp)
      external_dataset(jcomp) = aq_species_constraint%external_dataset(icomp)
      
      ! link constraint species
      select case(constraint_type(jcomp))
        case(CONSTRAINT_MINERAL)
          found = PETSC_FALSE
          do imnrl = 1, reaction%nmnrl
            if (StringCompare(constraint_aux_string(jcomp), &
                                reaction%mineral_names(imnrl), &
                                MAXWORDLENGTH)) then
              constraint_id(jcomp) = imnrl
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            option%io_buffer = 'Constraint mineral: ' // &
                     trim(constraint_aux_string(jcomp)) // &
                     ' for aqueous species: ' // &
                     trim(reaction%primary_species_names(jcomp)) // &
                     ' in constraint: ' // &
                     trim(constraint_name) // ' not found.' 
            call printErrMsg(option)         
          endif
        case(CONSTRAINT_GAS, CONSTRAINT_SUPERCRIT_CO2)
          found = PETSC_FALSE
          do igas = 1, reaction%ngas
            if (StringCompare(constraint_aux_string(jcomp), &
                                reaction%gas_species_names(igas), &
                                MAXWORDLENGTH)) then
              constraint_id(jcomp) = igas
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            option%io_buffer = 'Constraint gas: ' // &
                     trim(constraint_aux_string(jcomp)) // &
                     ' for aqueous species: ' // &
                     trim(reaction%primary_species_names(jcomp)) // &
                     ' in constraint: ' // &
                     trim(constraint_name) // ' not found.' 
            call printErrMsg(option)         
          endif
      end select

    endif
  enddo
  
  ! place ordered constraint parameters back in original arrays
  aq_species_constraint%constraint_type = constraint_type
  aq_species_constraint%constraint_aux_string = constraint_aux_string
  aq_species_constraint%constraint_spec_id = constraint_id
  aq_species_constraint%constraint_conc = constraint_conc
  aq_species_constraint%external_dataset = external_dataset

  ! minerals
  if (reaction%use_full_geochemistry .and. associated(mineral_constraint)) then
    constraint_mnrl_name = ''
    do imnrl = 1, reaction%nkinmnrl
      found = PETSC_FALSE
      do jmnrl = 1, reaction%nkinmnrl
        if (StringCompare(mineral_constraint%names(imnrl), &
                            reaction%kinmnrl_names(jmnrl), &
                            MAXWORDLENGTH)) then
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = &
                 'Mineral ' // trim(mineral_constraint%names(imnrl)) // &
                 'from CONSTRAINT ' // trim(constraint_name) // &
                 ' not found among kinetic minerals.'
        call printErrMsg(option)
      else
        mineral_constraint%basis_vol_frac(jmnrl) = &
          mineral_constraint%constraint_vol_frac(imnrl)
        mineral_constraint%basis_area(jmnrl) = &
          mineral_constraint%constraint_area(imnrl)
        constraint_mnrl_name(jmnrl) = mineral_constraint%names(imnrl)
      endif  
    enddo
    mineral_constraint%names = constraint_mnrl_name
    mineral_constraint%constraint_vol_frac = mineral_constraint%basis_vol_frac
    mineral_constraint%constraint_area = mineral_constraint%basis_area
  endif

  ! surface complexes
  if (reaction%use_full_geochemistry .and. associated(srfcplx_constraint)) then
    constraint_srfcplx_name = ''
    do isrfcplx = 1, reaction%nkinsrfcplx
      found = PETSC_FALSE
      do jsrfcplx = 1, reaction%nkinsrfcplx
        if (StringCompare(srfcplx_constraint%names(isrfcplx), &
                          reaction%kinsrfcplx_names(jsrfcplx), &
                            MAXWORDLENGTH)) then
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = &
                 'Surface complex ' // trim(srfcplx_constraint%names(isrfcplx)) // &
                 'from CONSTRAINT ' // trim(constraint_name) // &
                 ' not found among kinetic surface complexes.'
        call printErrMsg(option)
      else
        srfcplx_constraint%basis_conc(jsrfcplx) = &
          srfcplx_constraint%constraint_conc(isrfcplx)
        constraint_srfcplx_name(jsrfcplx) = srfcplx_constraint%names(isrfcplx)
      endif  
    enddo
    srfcplx_constraint%names = constraint_srfcplx_name
    srfcplx_constraint%constraint_conc = srfcplx_constraint%basis_conc
  endif
  
  ! colloids
  if (reaction%ncoll > 0) then
    if (.not.associated(colloid_constraint)) then
      option%io_buffer = 'Constraint "' // trim(constraint_name) // &
        'missing colloid entries when colloids present in problem.'
      call printErrMsg(option)
    endif
    constraint_colloid_name = ''
    do icoll = 1, reaction%ncoll
      found = PETSC_FALSE
      do jcoll = 1, reaction%ncoll
        if (StringCompare(colloid_constraint%names(icoll), &
                          reaction%colloid_names(jcoll), &
                            MAXWORDLENGTH)) then
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = &
                 'Surface complex ' // trim(colloid_constraint%names(icoll)) // &
                 'from CONSTRAINT ' // trim(constraint_name) // &
                 ' not found among colloids.'
        call printErrMsg(option)
      else
        colloid_constraint%basis_conc_mob(jcoll) = &
          colloid_constraint%constraint_conc_mob(icoll)
        colloid_constraint%basis_conc_imb(jcoll) = &
          colloid_constraint%constraint_conc_imb(icoll)
        constraint_colloid_name(jcoll) = colloid_constraint%names(icoll)
      endif  
    enddo
    colloid_constraint%names = constraint_colloid_name
    colloid_constraint%constraint_conc_mob = colloid_constraint%basis_conc_mob
    colloid_constraint%constraint_conc_imb = colloid_constraint%basis_conc_imb
  endif

end subroutine ReactionProcessConstraint

! ************************************************************************** !
!
! ReactionEquilibrateConstraint: Equilibrates constraint concentrations
!                                with prescribed geochemistry
! author: Glenn Hammond
! date: 10/22/08
!
! ************************************************************************** !
subroutine ReactionEquilibrateConstraint(rt_auxvar,global_auxvar, &
                                         reaction,constraint_name, &
                                         aq_species_constraint, &
                                         srfcplx_constraint, &
                                         colloid_constraint, &
                                         num_iterations, &
                                         initialize_rt_auxvar,option)
  use Option_module
  use Input_module
  use String_module  
  use Utility_module  
#ifdef CHUAN_CO2
  use co2eos_module, only: Henry_duan_sun
  use span_wagner_module, only: co2_span_wagner
  use water_eos_module
#endif  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(srfcplx_constraint_type), pointer :: srfcplx_constraint
  type(colloid_constraint_type), pointer :: colloid_constraint
  PetscInt :: num_iterations
  PetscBool :: initialize_rt_auxvar
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icomp, jcomp, kcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: icplx
  PetscInt :: irxn, isite, ncplx, k
  PetscInt :: igas
  PetscReal :: conc(reaction%naqcomp)
  PetscInt :: constraint_type(reaction%naqcomp)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(reaction%naqcomp)

  PetscReal :: Res(reaction%naqcomp)
  PetscReal :: total_conc(reaction%naqcomp)
  PetscReal :: free_conc(reaction%naqcomp)
  PetscReal :: Jac(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: indices(reaction%naqcomp)
  PetscReal :: norm
  PetscReal :: prev_molal(reaction%naqcomp)
  PetscReal, parameter :: tol = 1.d-12
  PetscReal, parameter :: tol_loose = 1.d-6
  PetscBool :: compute_activity_coefs

  PetscInt :: constraint_id(reaction%naqcomp)
  PetscReal :: lnQK, QK
  PetscReal :: tempreal
  PetscReal :: pres, tc, xphico2, henry, m_na, m_cl, xmass 
  PetscInt :: comp_id
  PetscReal :: convert_molal_to_molar
  PetscReal :: convert_molar_to_molal
  
  PetscBool :: charge_balance_warning_flag = PETSC_FALSE

  PetscReal :: Jac_num(reaction%naqcomp)
  PetscReal :: Res_pert, pert, prev_value

  PetscInt :: iphase
  PetscInt :: idof
  PetscInt :: istartaq, iendaq
  PetscInt :: kinmr_nrate_store, irate

#ifdef CHUAN_CO2  
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
  PetscInt :: iflag
  PetscErrorCode :: ierr
#endif
    
  constraint_type = aq_species_constraint%constraint_type
  constraint_aux_string = aq_species_constraint%constraint_aux_string
  constraint_id = aq_species_constraint%constraint_spec_id
  conc = aq_species_constraint%constraint_conc

  istartaq = reaction%offset_aq
  iendaq = reaction%offset_aq + reaction%naqcomp    

  iphase = 1
  
  xmass = 1.d0  
  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  
  if (reaction%initialize_with_molality) then
    convert_molal_to_molar = global_auxvar%den_kg(iphase)*xmass/1000.d0
    convert_molar_to_molal = 1.d0
  else
    convert_molal_to_molar = 1.d0
    convert_molar_to_molal = 1000.d0/global_auxvar%den_kg(iphase)/xmass
  endif
  
  if (associated(colloid_constraint)) then      
    colloid_constraint%basis_conc_mob = colloid_constraint%constraint_conc_mob        
    colloid_constraint%basis_conc_imb = colloid_constraint%constraint_conc_imb        
    rt_auxvar%colloid%conc_mob = colloid_constraint%basis_conc_mob* &
                                 convert_molar_to_molal
    rt_auxvar%colloid%conc_imb = colloid_constraint%basis_conc_imb* &
                                 convert_molar_to_molal
  endif  
  
  if (.not.reaction%use_full_geochemistry) then
    aq_species_constraint%basis_molarity = conc ! don't need to convert
    rt_auxvar%pri_molal = aq_species_constraint%basis_molarity* &
                          convert_molar_to_molal
    rt_auxvar%total(:,iphase) = aq_species_constraint%basis_molarity
    return
  endif

  ! if using multirate reaction, we need to turn it off to equilibrate the system
  ! then turn it back on
  kinmr_nrate_store = 0
  if (reaction%kinmr_nrate > 0) then
    kinmr_nrate_store = reaction%kinmr_nrate
    reaction%kinmr_nrate = 0
    allocate(rt_auxvar%dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp))
  endif
  
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    if (associated(reaction%eqcplx_logKcoef)) then
      call ReactionInterpolateLogK(reaction%eqcplx_logKcoef,reaction%eqcplx_logK, &
                                   global_auxvar%temp(iphase),reaction%neqcplx)
    endif
    if (associated(reaction%eqgas_logKcoef)) then
      call ReactionInterpolateLogK(reaction%eqgas_logKcoef,reaction%eqgas_logK, &
                                   global_auxvar%temp(iphase),reaction%ngas)
    endif
    if (associated(reaction%eqsrfcplx_logKcoef)) then
      call ReactionInterpolateLogK(reaction%eqsrfcplx_logKcoef,reaction%eqsrfcplx_logK, &
                                   global_auxvar%temp(iphase),reaction%neqsrfcplx)
    endif
    if (associated(reaction%kinmnrl_logKcoef)) then
      call ReactionInterpolateLogK(reaction%kinmnrl_logKcoef,reaction%kinmnrl_logK, &
                                   global_auxvar%temp(iphase),reaction%nkinmnrl)
    endif
    if (associated(reaction%mnrl_logKcoef)) then
      call ReactionInterpolateLogK(reaction%mnrl_logKcoef,reaction%mnrl_logK, &
                                   global_auxvar%temp(iphase),reaction%nmnrl)
    endif
  endif
#endif  
  
  total_conc = 0.d0
  do icomp = 1, reaction%naqcomp
    select case(constraint_type(icomp))
      case(CONSTRAINT_NULL,CONSTRAINT_TOTAL,CONSTRAINT_TOTAL_SORB)
        total_conc(icomp) = conc(icomp)*convert_molal_to_molar
        free_conc(icomp) = 1.d-9
      case(CONSTRAINT_FREE)
        free_conc(icomp) = conc(icomp)*convert_molar_to_molal
      case(CONSTRAINT_LOG)
        free_conc(icomp) = (10.d0**conc(icomp))*convert_molar_to_molal
      case(CONSTRAINT_CHARGE_BAL)
        free_conc(icomp) = conc(icomp)*convert_molar_to_molal ! just a guess
      case(CONSTRAINT_PH)
        ! check if H+ id set
        if (associated(reaction%species_idx)) then
          if (reaction%species_idx%h_ion_id /= 0) then
            ! check if icomp is H+
            if (reaction%species_idx%h_ion_id /= icomp) then
              string = 'OH-'
              if (.not.StringCompare(reaction%primary_species_names(icomp), &
                                     string,MAXWORDLENGTH)) then
                option%io_buffer = &
                         'pH specified as constraint (constraint =' // &
                         trim(constraint_name) // &
                         ') for species other than H+ or OH-: ' // &
                         trim(reaction%primary_species_names(icomp))
                call printErrMsg(option)
              endif
            endif
            free_conc(icomp) = 10.d0**(-conc(icomp))
          else
            option%io_buffer = &
                     'pH specified as constraint (constraint =' // &
                     trim(constraint_name) // &
                     '), but H+ not found in chemical species.'
            call printErrMsg(option)
          endif
        endif        
      case(CONSTRAINT_MINERAL)
        free_conc(icomp) = conc(icomp)*convert_molar_to_molal ! guess
      case(CONSTRAINT_GAS, CONSTRAINT_SUPERCRIT_CO2)
        if (conc(icomp) <= 0.d0) then ! log form
          conc(icomp) = 10.d0**conc(icomp) ! conc log10 partial pressure gas
        endif
        free_conc(icomp) = 1.d-9 ! guess
    end select
  enddo
  
  if (initialize_rt_auxvar) then
    rt_auxvar%pri_molal = free_conc
  endif

  num_iterations = 0
  compute_activity_coefs = PETSC_FALSE
  
  do

    if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF .and. &
        compute_activity_coefs) then
      call RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
      if(option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE)then
            call CO2AqActCoeff(rt_auxvar,global_auxvar,reaction,option)  
       endif
     endif
    call RTotal(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) call RTotalSorb(rt_auxvar,global_auxvar,reaction,option)
    
    ! geh - for debugging
    !call RTPrintAuxVar(rt_auxvar,reaction,option)

    Jac = 0.d0

! for colloids later on    
!    if (reaction%ncoll > 0) then
!      do idof = istartcoll, iendcoll
!        Jac(idof,idof) = 1.d0
!      enddo
!    endif
        
    do icomp = 1, reaction%naqcomp

      select case(constraint_type(icomp))
      
        case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
          ! units = mol/L water
          Res(icomp) = rt_auxvar%total(icomp,1) - total_conc(icomp)
          ! dtotal units = kg water/L water

          ! Jac units = kg water/L water
          Jac(icomp,:) = rt_auxvar%aqueous%dtotal(icomp,:,1)
        case(CONSTRAINT_TOTAL_SORB)
          ! conversion from m^3 bulk -> L water
          tempreal = option%reference_porosity*option%reference_saturation*1000.d0
          ! total = mol/L water  total_sorb = mol/m^3 bulk
          Res(icomp) = rt_auxvar%total(icomp,1) + &
            rt_auxvar%total_sorb_eq(icomp)/tempreal - total_conc(icomp)
          ! dtotal units = kg water/L water
          ! dtotal_sorb units = kg water/m^3 bulk
          ! Jac units = kg water/L water
          Jac(icomp,:) = rt_auxvar%aqueous%dtotal(icomp,:,1) + &
          ! dtotal_sorb units = kg water/m^3 bulk
                         rt_auxvar%dtotal_sorb_eq(icomp,:)/tempreal

        case(CONSTRAINT_FREE,CONSTRAINT_LOG)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
!          Jac(:,icomp) = 0.d0
          Jac(icomp,icomp) = 1.d0
          
        case(CONSTRAINT_CHARGE_BAL)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
          do jcomp = 1, reaction%naqcomp
            Res(icomp) = Res(icomp) + reaction%primary_spec_Z(jcomp) * &
              rt_auxvar%total(jcomp,1)
            do kcomp = 1, reaction%naqcomp
              Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                reaction%primary_spec_Z(kcomp)*rt_auxvar%aqueous%dtotal(kcomp,jcomp,1)
            enddo
          enddo
          if (rt_auxvar%pri_molal(icomp) < 1.d-20 .and. &
              .not.charge_balance_warning_flag) then
            if ((Res(icomp) > 0.d0 .and. &
                 reaction%primary_spec_Z(icomp) > 0.d0) .or. &
                (Res(icomp) < 0.d0 .and. &
                 reaction%primary_spec_Z(icomp) < 0.d0)) then
              option%io_buffer = &
                       'Charge balance species ' // &
                       trim(reaction%primary_species_names(icomp)) // &
                       ' may not satisfy constraint ' // &
                       trim(constraint_name) // &
                       '.  Molality already below 1.e-20.'
              call printMsg(option)
              charge_balance_warning_flag = PETSC_TRUE
              rt_auxvar%pri_molal(icomp) = 1.e-3 ! reset guess
            endif
          endif
          
        case(CONSTRAINT_PH)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
!          Jac(:,icomp) = 0.d0
          Jac(icomp,icomp) = 1.d0
          if (associated(reaction%species_idx)) then
            if (reaction%species_idx%h_ion_id > 0) then ! conc(icomp) = 10**-pH
              rt_auxvar%pri_molal(icomp) = 10.d0**(-conc(icomp)) / &
                                            rt_auxvar%pri_act_coef(icomp)
            else ! H+ is a complex
            
              icplx = abs(reaction%species_idx%h_ion_id)
              
              ! compute secondary species concentration
              ! *note that the sign was flipped below
              lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

              ! activity of water
              if (reaction%eqcplxh2oid(icplx) > 0) then
                lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
              endif

              do jcomp = 1, reaction%eqcplxspecid(0,icplx)
                comp_id = reaction%eqcplxspecid(jcomp,icplx)
                lnQK = lnQK + reaction%eqcplxstoich(jcomp,icplx)* &
                              log(rt_auxvar%pri_molal(comp_id)* &
                              rt_auxvar%pri_act_coef(comp_id))
              enddo
              lnQK = lnQK - log(conc(icomp)) ! this is log activity H+
              QK = exp(lnQK)
              
              Res(icomp) = 1.d0 - QK

              do jcomp = 1,reaction%eqcplxspecid(0,icplx)
                comp_id = reaction%eqcplxspecid(jcomp,icplx)
                Jac(icomp,comp_id) = -exp(lnQK-log(rt_auxvar%pri_molal(comp_id)))* &
                                          reaction%eqcplxstoich(jcomp,icplx)
              enddo
            endif
          endif
                      
        case(CONSTRAINT_MINERAL)

          imnrl = constraint_id(icomp)
          ! compute secondary species concentration
          lnQK = -reaction%mnrl_logK(imnrl)*LOG_TO_LN

          ! activity of water
          if (reaction%mnrlh2oid(imnrl) > 0) then
            lnQK = lnQK + reaction%mnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
          endif

          do jcomp = 1, reaction%mnrlspecid(0,imnrl)
            comp_id = reaction%mnrlspecid(jcomp,imnrl)
            lnQK = lnQK + reaction%mnrlstoich(jcomp,imnrl)* &
                          log(rt_auxvar%pri_molal(comp_id)* &
                          rt_auxvar%pri_act_coef(comp_id))
          enddo
!         QK = exp(lnQK)
          
!         Res(icomp) = 1.d0 - QK
          Res(icomp) = lnQK

          do jcomp = 1,reaction%mnrlspecid(0,imnrl)
            comp_id = reaction%mnrlspecid(jcomp,imnrl)
!           Jac(icomp,comp_id) = -QK/auxvar%primary_spec(comp_id)* &
!                                reaction%mnrlstoich(jcomp,imnrl)
            Jac(icomp,comp_id) = reaction%mnrlstoich(jcomp,imnrl)/ &
              rt_auxvar%pri_molal(comp_id)
                                 
          enddo
  
        case(CONSTRAINT_GAS)

          igas = constraint_id(icomp)
          
          ! compute secondary species concentration
           lnQK = -reaction%eqgas_logK(igas)*LOG_TO_LN
 
          ! divide K by RT
          !lnQK = lnQK - log((auxvar%temp+273.15d0)*IDEAL_GAS_CONST)
          
          ! activity of water
          if (reaction%eqgash2oid(igas) > 0) then
            lnQK = lnQK + reaction%eqgash2ostoich(igas)*rt_auxvar%ln_act_h2o
          endif

          do jcomp = 1, reaction%eqgasspecid(0,igas)
            comp_id = reaction%eqgasspecid(jcomp,igas)
            lnQK = lnQK + reaction%eqgasstoich(jcomp,igas)* &
              log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
          enddo
          
!         QK = exp(lnQK)
          
!         Res(icomp) = QK - conc(icomp)
          Res(icomp) = lnQK - log(conc(icomp)) ! gas pressure
          Jac(icomp,:) = 0.d0
          do jcomp = 1,reaction%eqgasspecid(0,igas)
            comp_id = reaction%eqgasspecid(jcomp,igas)
!           Jac(icomp,comp_id) = QK/auxvar%primary_spec(comp_id)* &
!                                reaction%eqgasstoich(jcomp,igas)
            Jac(icomp,comp_id) = reaction%eqgasstoich(jcomp,igas)/ &
              rt_auxvar%pri_molal(comp_id)

!#ifdef CHUAN_CO2
!            print *,'Gas CO2 constraint Jac,',igas, icomp, comp_id, &
!              reaction%eqgasstoich(jcomp,igas),&
!              Jac(icomp,comp_id), rt_auxvar%pri_molal(comp_id), lnQK
!#endif
          enddo

#ifdef CHUAN_CO2        
        case(CONSTRAINT_SUPERCRIT_CO2)
          
          igas = constraint_id(icomp)
         
          ! compute secondary species concentration
          if(abs(reaction%species_idx%co2_gas_id) == igas) then
          
!           pres = global_auxvar%pres(2)
            pres = conc(icomp)*1.D5
            global_auxvar%pres(2) = pres

!           print *,'reaction-SC: ',icomp,igas,pres,conc(icomp)*1.d5
            
            tc = global_auxvar%temp(1)

            call PSAT(tc, sat_pressure, ierr)
            
            pco2 = conc(icomp)*1.e5
!           pco2 = pres - sat_pressure
            
!            pres = conc(icomp)*1.D5 + sat_pressure
            yco2 = pco2/pres
             
            iflag = 1
            call co2_span_wagner(pres*1D-6,tc+273.15D0,dg,dddt,dddp,fg, &
              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag,option%itable)

!            call co2_span_wagner(pco2*1D-6,tc+273.15D0,dg,dddt,dddp,fg, &
!              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
            
            global_auxvar%den_kg(2) = dg
            
            !compute fugacity coefficient
            fg = fg*1.D6
            xphico2 = fg / pres
!            global_auxvar%fugacoeff(1) = xphico2
            
!           call Henry_duan_sun_0NaCl(pco2*1.d-5, tc, henry)
            m_na = 0.d0
            m_cl = 0.d0
            if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
              m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
              m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
!              call Henry_duan_sun(tc,pco2*1D-5,henry,xphico2,lngamco2, &
!                m_na,m_cl,sat_pressure*1D-5)
              call Henry_duan_sun(tc,pres*1D-5,henry,xphico2,lngamco2, &
                m_na,m_cl,sat_pressure*1D-5)

            else
              call Henry_duan_sun(tc,pres*1D-5,henry,xphico2,lngamco2, &
                option%m_nacl,option%m_nacl,sat_pressure*1D-5)
             !   print *, 'SC: mnacl=', option%m_nacl,'stioh2o=',reaction%eqgash2ostoich(igas)
            endif
            
            lnQk = -log(xphico2*henry)-lngamco2
!            lnQk = -log(xphico2*henry)
!           lnQk = log(fg/henry)

            reaction%eqgas_logK(igas) = -lnQK*LN_TO_LOG
            
            !print *, 'SC CO2 constraint',igas,pres,pco2,tc,xphico2,henry,lnQk,yco2, &
            !   lngamco2,m_na,m_cl,reaction%eqgas_logK(igas),rt_auxvar%ln_act_h2o,&
            !    reaction%eqgash2oid(igas), global_auxvar%fugacoeff(1)
            
            ! activity of water
            if (reaction%eqgash2oid(igas) > 0) then
              lnQK = lnQK + reaction%eqgash2ostoich(igas)*rt_auxvar%ln_act_h2o
            endif
            do jcomp = 1, reaction%eqgasspecid(0,igas)
              comp_id = reaction%eqgasspecid(jcomp,igas)
              lnQK = lnQK + reaction%eqgasstoich(jcomp,igas)* &
!                log(rt_auxvar%pri_molal(comp_id))
               log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
            !    print *,'SC: ',rt_auxvar%pri_molal(comp_id), &
            !      rt_auxvar%pri_act_coef(comp_id),exp(lngamco2)
            enddo
          
!           QK = exp(lnQK)
             
!           Res(icomp) = QK - conc(icomp)
            Res(icomp) = lnQK - log(pco2*1D-5)!log(conc(icomp)) ! gas pressure bars
            Jac(icomp,:) = 0.d0
            do jcomp = 1,reaction%eqgasspecid(0,igas)
              comp_id = reaction%eqgasspecid(jcomp,igas)
!             Jac(icomp,comp_id) = QK/auxvar%primary_spec(comp_id)* &
!                                reaction%eqgasstoich(jcomp,igas)
              Jac(icomp,comp_id) = reaction%eqgasstoich(jcomp,igas)/ &
                rt_auxvar%pri_molal(comp_id)
              
!             print *,'SC CO2 constraint Jac,',igas, icomp, comp_id, &
!               reaction%eqgasstoich(jcomp,igas),&
!               Jac(icomp,comp_id), rt_auxvar%pri_molal(comp_id),conc(icomp) 
            enddo
         endif       
#endif           
      end select
    enddo
    
    ! scale Jacobian
    do icomp = 1, reaction%naqcomp
      norm = max(1.d0,maxval(abs(Jac(icomp,:))))
      norm = 1.d0/norm
      Res(icomp) = Res(icomp)*norm
      Jac(icomp,:) = Jac(icomp,:)*norm
    enddo

    ! for derivatives with respect to ln conc
    do icomp = 1, reaction%naqcomp
      Jac(:,icomp) = Jac(:,icomp)*rt_auxvar%pri_molal(icomp)
    enddo
    
    call ludcmp(Jac,reaction%naqcomp,indices,icomp)
    call lubksb(Jac,reaction%naqcomp,indices,Res)

    prev_molal = rt_auxvar%pri_molal

    Res = dsign(1.d0,Res)*min(dabs(Res),5.d0)
    
    rt_auxvar%pri_molal = rt_auxvar%pri_molal*exp(-Res)

    num_iterations = num_iterations + 1
    
    if (mod(num_iterations,1000) == 0) then
100   format('Constraint iteration count has exceeded: ',i5)
      write(option%io_buffer,100) num_iterations
      call printMsg(option)
      do icomp=1,reaction%naqcomp
        write(option%io_buffer,200) reaction%primary_species_names(icomp), &
        prev_molal(icomp),Res(icomp)
        call printMsg(option)
      enddo
200   format(a12,1x,1p2e12.4)
      if (num_iterations >= 10000) then
        option%io_buffer = 'Stopping due to excessive iteration count!'
        call printErrMsg(option)
      endif
    endif
    
    ! check for convergence
    if (maxval(dabs(res)) < tol) then
      ! need some sort of convergence before we kick in activities
      if (compute_activity_coefs) exit
      compute_activity_coefs = PETSC_TRUE
    endif

  enddo

  ! once equilibrated, compute sorbed concentrations
  if (reaction%neqsorb > 0) call RTotalSorb(rt_auxvar,global_auxvar,reaction,option)

  ! WARNING: below assumes site concentration multiplicative factor
  if (kinmr_nrate_store > 0) then
    reaction%kinmr_nrate = kinmr_nrate_store
    kinmr_nrate_store = 0
    do irate = 1, reaction%kinmr_nrate
      rt_auxvar%kinmr_total_sorb(:,irate) = reaction%kinmr_frac(irate) * rt_auxvar%total_sorb_eq
    enddo
    deallocate(rt_auxvar%dtotal_sorb_eq)
    nullify(rt_auxvar%dtotal_sorb_eq)
  endif

  if (reaction%nkinsrfcplx > 0 .and. associated(srfcplx_constraint)) then
  ! compute surface complex conc. at new time step (5.1-30) 
    rt_auxvar%kinsrfcplx_conc = srfcplx_constraint%constraint_conc
    do irxn = 1, reaction%nkinsrfcplxrxn
      isite = reaction%kinsrfcplx_rxn_to_site(irxn)
      rt_auxvar%kinsrfcplx_free_site_conc(isite) = reaction%kinsrfcplx_rxn_site_density(isite)
      ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
      do k = 1, ncplx ! ncplx in rxn
        icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
        rt_auxvar%kinsrfcplx_free_site_conc(isite) = &
          rt_auxvar%kinsrfcplx_free_site_conc(isite) - &
          rt_auxvar%kinsrfcplx_conc(icplx)
      enddo
    enddo
    do irxn = 1, reaction%nkinsrfcplxrxn
      isite = reaction%kinsrfcplx_rxn_to_site(irxn)
      if (rt_auxvar%kinsrfcplx_free_site_conc(isite) < 0.d0) then
        option%io_buffer = 'Free site concentration for site ' // &
          trim(reaction%kinsrfcplx_site_names(isite)) // &
          ' is less than zero.'
        call printErrMsg(option)
      endif
    enddo
    srfcplx_constraint%constraint_free_site_conc = rt_auxvar%kinsrfcplx_free_site_conc
    srfcplx_constraint%basis_free_site_conc = srfcplx_constraint%constraint_free_site_conc
  endif
  
  ! do not scale by molal_to_molar since it could be 1.d0 if MOLAL flag set
  aq_species_constraint%basis_molarity = rt_auxvar%pri_molal* &
                                         global_auxvar%den_kg(option%liquid_phase)/ &
                                         1000.d0

#if 0
  call RCalculateCompression(global_auxvar,rt_auxvar,reaction,option)
#endif

! this is performed above
!  if (associated(colloid_constraint%colloids)) then                        
!    colloid_constraint%colloids%basis_conc_mob = rt_auxvar%colloid%conc_mob* &
!                           global_auxvar%den_kg(option%liquid_phase)/1000.d0
!    colloid_constraint%colloids%basis_conc_imb = rt_auxvar%colloid%conc_imb* &
!                           global_auxvar%den_kg(option%liquid_phase)/1000.d0
!  endif
    
!  write(option%io_buffer,111) trim(constraint_name),num_iterations
!  call printMsg(option)
!111 format(' Equilibrate Constraint: ',a30,i4)

end subroutine ReactionEquilibrateConstraint

! ************************************************************************** !
!
! ReactionPrintConstraint: Prints a constraint associated with reactive 
!                          transport
! author: Glenn Hammond
! date: 10/28/08
!
! ************************************************************************** !
subroutine ReactionPrintConstraint(constraint_coupler,reaction,option)

  use Option_module
  use Input_module
  use String_module
  use Condition_module

  implicit none
  
  type(option_type) :: option
  type(tran_constraint_coupler_type) :: constraint_coupler
  type(reaction_type), pointer :: reaction
  
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, icomp, irxn, j, jj, ncomp, ncplx
  PetscInt :: icplx, icplx2
  PetscInt :: imnrl,igas
  PetscInt :: eqcplxsort(reaction%neqcplx+1)
  PetscInt :: eqcplxid(reaction%neqcplx+1)
  PetscInt :: eqminsort(reaction%nmnrl)
  PetscInt :: eqsrfcplxsort(reaction%neqsrfcplx+reaction%neqsrfcplxrxn)
  PetscBool :: finished, found
  PetscReal :: conc, conc2
  PetscReal :: lnQK(reaction%nmnrl), QK(reaction%nmnrl)
  PetscReal :: lnQKgas(reaction%ngas), QKgas(reaction%ngas)
  PetscReal :: charge_balance, ionic_strength
  PetscReal :: percent(reaction%neqcplx+1)
  PetscReal :: totj, retardation, kd
  PetscInt :: comp_id, jcomp
  PetscInt :: icount
  PetscInt :: iphase
  PetscReal :: bulk_vol_to_fluid_vol, molar_to_molal, molal_to_molar
  PetscReal :: sum_molality, sum_mass, mole_fraction_h2o, mass_fraction_h2o, &
               mass_fraction_co2, mole_fraction_co2

  aq_species_constraint => constraint_coupler%aqueous_species
  mineral_constraint => constraint_coupler%minerals
  
  iphase = 1

  90 format(2x,76('-'))
  91 format(a)

  write(option%fid_out,'(/,''  Constraint: '',a)') &
    trim(constraint_coupler%constraint_name)

  rt_auxvar => constraint_coupler%rt_auxvar
  global_auxvar => constraint_coupler%global_auxvar

  select case(option%iflowmode)
    case(FLASH2_MODE,MPH_MODE,IMS_MODE)
    case(NULL_MODE)
      global_auxvar%den_kg(iphase) = option%reference_water_density
      global_auxvar%temp(1) = option%reference_temperature
      global_auxvar%sat(iphase) = option%reference_saturation
    case(RICHARDS_MODE)
      global_auxvar%temp(1) = option%reference_temperature
  end select
        
!  global_auxvar%den_kg(iphase) = option%reference_water_density
!  global_auxvar%temp(1) = option%reference_temperature
!  global_auxvar%sat(iphase) = option%reference_saturation
  bulk_vol_to_fluid_vol = option%reference_porosity*global_auxvar%sat(iphase)*1000.d0

! compute mass fraction of H2O
  if (reaction%use_full_geochemistry) then
    sum_molality = 0.d0
    do icomp = 1, reaction%naqcomp
      sum_molality = sum_molality + rt_auxvar%pri_molal(icomp)
    enddo
    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        sum_molality = sum_molality + rt_auxvar%sec_molal(i)
      enddo
    endif
    mole_fraction_h2o = 1.d0/(1.d0+FMWH2O*sum_molality*1.d-3)

    sum_mass = 0.d0
    do icomp = 1, reaction%naqcomp
      sum_mass = sum_mass + reaction%primary_spec_molar_wt(icomp)*rt_auxvar%pri_molal(icomp)
    enddo
    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        sum_mass = sum_mass + reaction%eqcplx_molar_wt(i)*rt_auxvar%sec_molal(i)
      enddo
    endif
    mass_fraction_h2o = 1.d0/(1.d0 + sum_mass*1.d-3)
  endif
  
  molal_to_molar = global_auxvar%den_kg(iphase)/1000.d0
  molar_to_molal = 1.d0/molal_to_molar
    
  if (.not.reaction%use_full_geochemistry) then
    100 format(/,'  species       molality')  
    write(option%fid_out,100)
    101 format(2x,a12,es12.4)
    do icomp = 1, reaction%naqcomp
      write(option%fid_out,101) reaction%primary_species_names(icomp), &
                                rt_auxvar%pri_molal(icomp)
    enddo
  else

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    if (associated(reaction%eqcplx_logKcoef)) then
      call ReactionInterpolateLogK(reaction%eqcplx_logKcoef,reaction%eqcplx_logK, &
                                   global_auxvar%temp(iphase),reaction%neqcplx)
    endif
    if (associated(reaction%eqgas_logKcoef)) then
      call ReactionInterpolateLogK(reaction%eqgas_logKcoef,reaction%eqgas_logK, &
                                   global_auxvar%temp(iphase),reaction%ngas)
    endif
    if (associated(reaction%eqsrfcplx_logKcoef)) then
      call ReactionInterpolateLogK(reaction%eqsrfcplx_logKcoef,reaction%eqsrfcplx_logK, &
                                   global_auxvar%temp(iphase),reaction%neqsrfcplx)
    endif
    if (associated(reaction%kinmnrl_logKcoef)) then
      call ReactionInterpolateLogK(reaction%kinmnrl_logKcoef,reaction%kinmnrl_logK, &
                                   global_auxvar%temp(iphase),reaction%nkinmnrl)
    endif
    if (associated(reaction%mnrl_logKcoef)) then
      call ReactionInterpolateLogK(reaction%mnrl_logKcoef,reaction%mnrl_logK, &
                                   global_auxvar%temp(iphase),reaction%nmnrl)
    endif
  endif
#endif  

    200 format('')
    201 format(a20,i5)
    202 format(a20,f10.2)
    203 format(a20,f8.4)
    204 format(a20,es12.4)
    write(option%fid_out,90)
    write(option%fid_out,201) '      iterations: ', &
      constraint_coupler%num_iterations
    if (associated(reaction%species_idx)) then
      if (reaction%species_idx%h_ion_id > 0) then
        write(option%fid_out,203) '              pH: ', &
          -log10(rt_auxvar%pri_molal(reaction%species_idx%h_ion_id)* &
                 rt_auxvar%pri_act_coef(reaction%species_idx%h_ion_id))
      else if (reaction%species_idx%h_ion_id < 0) then
        write(option%fid_out,203) '              pH: ', &
          -log10(rt_auxvar%sec_molal(abs(reaction%species_idx%h_ion_id))* &
                 rt_auxvar%sec_act_coef(abs(reaction%species_idx%h_ion_id)))
      endif
    endif
    
    ionic_strength = 0.d0
    charge_balance = 0.d0
    do icomp = 1, reaction%naqcomp      
      charge_balance = charge_balance + rt_auxvar%total(icomp,1)* &
                                        reaction%primary_spec_Z(icomp)
      ionic_strength = ionic_strength + rt_auxvar%pri_molal(icomp)* &
        reaction%primary_spec_Z(icomp)*reaction%primary_spec_Z(icomp)
    enddo
    
    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        ionic_strength = ionic_strength + rt_auxvar%sec_molal(i)* &
                                          reaction%eqcplx_Z(i)* &
                                          reaction%eqcplx_Z(i)
      enddo
    endif
    ionic_strength = 0.5d0 * ionic_strength
    
    write(option%fid_out,'(a20,es12.4,a8)') '  ionic strength: ', &
      ionic_strength,' [mol/L]'
    write(option%fid_out,204) '  charge balance: ', charge_balance
    
    write(option%fid_out,'(a20,1pe12.4,a5)') '        pressure: ', &
      global_auxvar%pres(1),' [Pa]'
    write(option%fid_out,'(a20,f8.2,a4)') '     temperature: ', &
      global_auxvar%temp(1),' [C]'
    write(option%fid_out,'(a20,f8.2,a9)') '     density H2O: ', &
      global_auxvar%den_kg(1),' [kg/m^3]'
    write(option%fid_out,'(a20,1p2e12.4,a9)') 'ln / activity H2O: ', &
      rt_auxvar%ln_act_h2o,exp(rt_auxvar%ln_act_h2o),' [---]'
    write(option%fid_out,'(a20,1pe12.4,a9)') 'mole fraction H2O: ', &
      mole_fraction_h2o,' [---]'
    write(option%fid_out,'(a20,1pe12.4,a9)') 'mass fraction H2O: ', &
      mass_fraction_h2o,' [---]'
#ifdef CHUAN_CO2
    if (option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE) then
      if (global_auxvar%den_kg(2) > 0.d0) then
        write(option%fid_out,'(a20,f8.2,a9)') '     density CO2: ', &
          global_auxvar%den_kg(2),' [kg/m^3]'
        write(option%fid_out,'(a20,es12.4,a9)') '            xphi: ', &
          global_auxvar%fugacoeff(1)

        if (reaction%species_idx%co2_aq_id /= 0) then
          icomp = reaction%species_idx%co2_aq_id
          mass_fraction_co2 = reaction%primary_spec_molar_wt(icomp)*rt_auxvar%pri_molal(icomp)* &
            mass_fraction_h2o*1.d-3
          mole_fraction_co2 = rt_auxvar%pri_molal(icomp)*FMWH2O*mole_fraction_h2o*1.e-3
          write(option%fid_out,'(a20,es12.4,a9)') 'mole fraction CO2: ', &
            mole_fraction_co2
          write(option%fid_out,'(a20,es12.4,a9)') 'mass fraction CO2: ', &
            mass_fraction_co2
        endif
      endif
    endif
#endif
    write(option%fid_out,90)

    102 format(/,'                        free        total')  
    103 format('  species               molal       molal       act coef     constraint')  
    write(option%fid_out,102)
    write(option%fid_out,103)
    write(option%fid_out,90)
  
    104 format(2x,a20,es12.4,es12.4,es12.4,4x,a)
    do icomp = 1, reaction%naqcomp
      select case(aq_species_constraint%constraint_type(icomp))
        case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
          string = 'total'
        case(CONSTRAINT_TOTAL_SORB)
          string = 'aq+sorb'
        case(CONSTRAINT_FREE)
          string = 'free'
        case(CONSTRAINT_CHARGE_BAL)
          string = 'chrg'
        case(CONSTRAINT_LOG)
          string = 'log'
        case(CONSTRAINT_PH)
          string = 'pH'
        case(CONSTRAINT_MINERAL,CONSTRAINT_GAS)
          string = aq_species_constraint%constraint_aux_string(icomp)
        case(CONSTRAINT_SUPERCRIT_CO2)
          string = 'SC ' // aq_species_constraint%constraint_aux_string(icomp)
      end select
      write(option%fid_out,104) reaction%primary_species_names(icomp), &
                                rt_auxvar%pri_molal(icomp), &
                                rt_auxvar%total(icomp,1)*molar_to_molal, &
                                rt_auxvar%pri_act_coef(icomp), &
                                trim(string)
    enddo 
  endif 
      
  if (reaction%neqcplx > 0) then    
    ! sort complex concentrations from largest to smallest
    do i = 1, reaction%neqcplx
      eqcplxsort(i) = i
    enddo
    do
      finished = PETSC_TRUE
      do i = 1, reaction%neqcplx-1
        icplx = eqcplxsort(i)
        icplx2 = eqcplxsort(i+1)
        if (rt_auxvar%sec_molal(icplx) < &
            rt_auxvar%sec_molal(icplx2)) then
          eqcplxsort(i) = icplx2
          eqcplxsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo
            
    110 format(/,'  complex               molality    act coef  logK')  
    write(option%fid_out,110)
    write(option%fid_out,90)
    111 format(2x,a20,es12.4,es12.4,2x,es12.4)
    do i = 1, reaction%neqcplx ! for each secondary species
      icplx = eqcplxsort(i)
      write(option%fid_out,111) reaction%secondary_species_names(icplx), &
                                rt_auxvar%sec_molal(icplx), &
                                rt_auxvar%sec_act_coef(icplx), &
                                reaction%eqcplx_logK(icplx)
    enddo 

    !print speciation precentages
    write(option%fid_out,92)
    92 format(/)
    134 format(2x,'complex species       percent   molality')
    135 format(2x,'primary species: ',a20,2x,' total conc: ',1pe12.4)
    136 format(2x,a20,2x,f6.2,2x,1pe12.4,1p2e12.4)
    do icomp = 1, reaction%naqcomp
    
      eqcplxsort = 0
      eqcplxid = 0
      percent = 0.d0
      totj = 0.d0
      
      icount = 0
      do icplx = 1, reaction%neqcplx
        found = PETSC_FALSE
        do i = 1, reaction%naqcomp
          if (reaction%eqcplxspecid(i,icplx) == icomp) then
            icount = icount + 1
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (found) then
          eqcplxid(icount) = icplx
          percent(icount) = dabs(rt_auxvar%sec_molal(icplx)* &
                                 reaction%eqcplxstoich(i,icplx))
          totj = totj + percent(icount)
        endif
      enddo
      icount = icount + 1
      eqcplxid(icount) = -icomp
      percent(icount) = rt_auxvar%pri_molal(icomp)
      totj = totj + percent(icount)
      percent = percent / totj
      
      eqcplxsort = 0
      do i = 1, icount
        eqcplxsort(i) = i
      enddo
      
      do
        finished = PETSC_TRUE
        do i = 1, icount-1
          icplx = eqcplxsort(i)
          icplx2 = eqcplxsort(i+1)
          if (percent(abs(icplx)) < percent(abs(icplx2))) then
            eqcplxsort(i) = icplx2
            eqcplxsort(i+1) = icplx
            finished = PETSC_FALSE
          endif
        enddo
        if (finished) exit
      enddo

      write(option%fid_out,90)
      write(option%fid_out,135) reaction%primary_species_names(icomp), &
                                rt_auxvar%total(icomp,iphase)
      write(option%fid_out,134)
      write(option%fid_out,90)
      do i = 1, icount
        j = eqcplxsort(i)
        if (percent(j) < 0.0001d0) cycle
        icplx = eqcplxid(j)
        if (icplx < 0) then
          icplx = abs(icplx)
          write(option%fid_out,136) reaction%primary_species_names(icplx), &
                                    percent(j)*100.d0, &
                                    rt_auxvar%pri_molal(icplx)
        else
          write(option%fid_out,136) reaction%secondary_species_names(icplx), &
                                    percent(j)*100.d0, &
                                    rt_auxvar%sec_molal(icplx)
        endif
      enddo
    enddo

  endif 
          
  if (reaction%neqsrfcplxrxn > 0) then
    ! sort surface complex concentrations from largest to smallest
    ! note that we include free site concentrations; their ids negated
#if 0
    do i = 1, reaction%neqsrfcplx
      eqsrfcplxsort(i) = i
    enddo
    do i = 1, reaction%neqsrfcplxrxn
      eqsrfcplxsort(reaction%neqsrfcplx+i) = -i
    enddo
    do
      finished = PETSC_TRUE
      do i = 1, reaction%neqsrfcplx+reaction%neqsrfcplxrxn-1
        icplx = eqsrfcplxsort(i)
        icplx2 = eqsrfcplxsort(i+1)
        if (icplx > 0) then
          conc = rt_auxvar%eqsrfcplx_conc(icplx)
        else
          conc = rt_auxvar%eqsrfcplx_free_site_conc(-icplx)
        endif
        if (icplx2 > 0) then
          conc2 = rt_auxvar%eqsrfcplx_conc(icplx2)
        else
          conc2 = rt_auxvar%eqsrfcplx_free_site_conc(-icplx2)
        endif
        if (conc < conc2) then
          eqsrfcplxsort(i) = icplx2
          eqsrfcplxsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo
            
    write(option%fid_out,120)
    write(option%fid_out,90)
    do i = 1, reaction%neqsrfcplxrxn + reaction%neqsrfcplx
      icplx = eqsrfcplxsort(i)
      if (icplx > 0) then
        write(option%fid_out,121) reaction%eqsrfcplx_names(icplx), &
                                  rt_auxvar%eqsrfcplx_conc(icplx), &
                                  reaction%eqsrfcplx_logK(icplx)
      else
        write(option%fid_out,122) reaction%eqsrfcplx_site_names(-icplx), &
                                  rt_auxvar%eqsrfcplx_free_site_conc(-icplx)
      endif
    enddo
#endif

    120 format(/,'  surf complex          mol/m^3 blk logK')  
    121 format(2x,a20,es12.4,es12.4)
    122 format(2x,a20,es12.4,'  free site')

    write(option%fid_out,120)
    write(option%fid_out,90)
    do irxn = 1, reaction%neqsrfcplxrxn
      write(option%fid_out,122) reaction%eqsrfcplx_site_names(irxn), &
                                rt_auxvar%eqsrfcplx_free_site_conc(irxn)
      ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
      do i = 1, ncplx
        icplx = reaction%eqsrfcplx_rxn_to_complex(i,irxn)
        write(option%fid_out,121) reaction%eqsrfcplx_names(icplx), &
                                  rt_auxvar%eqsrfcplx_conc(icplx), &
                                  reaction%eqsrfcplx_logK(icplx)
      enddo
    enddo
  
  endif

! retardation
  if (reaction%neqsrfcplxrxn > 0) then
    write(option%fid_out,123)
    write(option%fid_out,90)
    do j = 1, reaction%naqcomp
      retardation = 1.d0
      do irxn = 1, reaction%neqsrfcplxrxn
        ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
        do i = 1, ncplx
          icplx = reaction%eqsrfcplx_rxn_to_complex(i,irxn)
          ncomp = reaction%eqsrfcplxspecid(0,icplx)
          do jj = 1, ncomp
            jcomp = reaction%eqsrfcplxspecid(jj,icplx)
            if (j == jcomp) then
              if (rt_auxvar%total(j,iphase) /= 0.d0) &
              retardation = retardation + &
                            reaction%eqsrfcplxstoich(jj,icplx)* &
                            rt_auxvar%eqsrfcplx_conc(icplx)/ &
                            bulk_vol_to_fluid_vol/ &
                            rt_auxvar%total(j,iphase)
              exit
            endif
          enddo
        enddo
      enddo
      write(option%fid_out,124) reaction%primary_species_names(j),retardation
    enddo
    123 format(/,'  primary species  retardation')  
    124 format(2x,a12,4x,1pe12.4)

#ifdef DOUBLE_LAYER
    call DoubleLayer (constraint_coupler,reaction,option)
#endif

  endif
  
  ! Ion Exchange
  if (reaction%neqionxrxn > 0) then
    write(option%fid_out,125)
    do irxn = 1, reaction%neqionxrxn
      write(option%fid_out,90)
      write(option%fid_out,126) reaction%eqionx_rxn_CEC(irxn)
      write(option%fid_out,127)
      write(option%fid_out,90)
      ncomp = reaction%eqionx_rxn_cationid(0,irxn)
      do i = 1, ncomp
        icomp = reaction%eqionx_rxn_cationid(i,irxn)
        kd = rt_auxvar%eqionx_conc(i,irxn)/rt_auxvar%total(icomp,iphase) & 
                      /bulk_vol_to_fluid_vol
        write(option%fid_out,128) reaction%primary_species_names(icomp), &
          reaction%eqionx_rxn_k(i,irxn), & 
          rt_auxvar%eqionx_conc(i,irxn), &
          kd
      enddo
    enddo
    125 format(/,2x,'ion-exchange reactions')
    126 format(2x,'CEC = ',1pe12.4)
    127 format(2x,'cation    selectivity coef.   sorbed conc.     Kd',&
               /,30x,'[mol/m^3]')
    128 format(2x,a8,2x,1pe12.4,4x,1pe12.4,4x,1pe12.4,4x,1pe12.4)
  endif
  
! total retardation from ion exchange and surface complexation
  if (reaction%neqsorb > 0) then
    write(option%fid_out,1128)
    write(option%fid_out,90)
    do jcomp = 1, reaction%naqcomp
      if (abs(rt_auxvar%total(jcomp,iphase)) > 0.d0) &
      retardation = 1.d0 + rt_auxvar%total_sorb_eq(jcomp)/bulk_vol_to_fluid_vol &
        /rt_auxvar%total(jcomp,iphase)
      totj = rt_auxvar%total(jcomp,iphase)+rt_auxvar%total_sorb_eq(jcomp)/bulk_vol_to_fluid_vol
      write(option%fid_out,129) reaction%primary_species_names(jcomp), &
        totj,retardation
    enddo
   1128 format(/,2x,'primary species     total(aq+sorbed)    total retardation', &
               /,25x,'[mol/L]',15x,'1+Kd')
    129 format(2x,a12,8x,1pe12.4,8x,1pe12.4)
  endif
  
  if (reaction%nmnrl > 0) then
  
    130 format(/,'  mineral                             log SI    log K')
    131 format(2x,a30,2x,f12.4,2x,1pe12.4)

    do imnrl = 1, reaction%nmnrl
      ! compute saturation
      lnQK(imnrl) = -reaction%mnrl_logK(imnrl)*LOG_TO_LN
      if (reaction%mnrlh2oid(imnrl) > 0) then
        lnQK(imnrl) = lnQK(imnrl) + reaction%mnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
      endif
      do jcomp = 1, reaction%mnrlspecid(0,imnrl)
        comp_id = reaction%mnrlspecid(jcomp,imnrl)
        lnQK(imnrl) = lnQK(imnrl) + reaction%mnrlstoich(jcomp,imnrl)* &
                      log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
      enddo
      QK(imnrl) = exp(lnQK(imnrl))    
    enddo

    ! sort mineral saturation indices from largest to smallest
    do i = 1, reaction%nmnrl
      eqminsort(i) = i
    enddo
    do
      finished = PETSC_TRUE
      do i = 1, reaction%nmnrl-1
        icplx = eqminsort(i)
        icplx2 = eqminsort(i+1)
        if (QK(icplx) < QK(icplx2)) then
          eqminsort(i) = icplx2
          eqminsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo

    write(option%fid_out,130)
    write(option%fid_out,90)
  
    do imnrl = 1, reaction%nmnrl
      i = eqminsort(imnrl)
      write(option%fid_out,131) reaction%mineral_names(i), &
                                lnQK(i)*LN_TO_LOG, &
                                reaction%mnrl_logK(i)
    enddo
  endif
    
  if (reaction%ngas > 0) then
    
    132 format(/,'  gas           log part. press.  part. press. [bars]   log K')
    133 format(2x,a10,2x,1pe12.4,6x,1pe12.4,8x,1pe12.4)
    
    write(option%fid_out,132)
    write(option%fid_out,90)
    
    do igas = 1, reaction%ngas
      
      ! compute gas partial pressure
      lnQKgas(igas) = -reaction%eqgas_logK(igas)*LOG_TO_LN
      
      ! divide K by RT
      !lnQKgas = lnQKgas - log((auxvar%temp+273.15d0)*IDEAL_GAS_CONST)
      
      ! activity of water
      if (reaction%eqgash2oid(igas) > 0) then
        lnQKgas(igas) = lnQKgas(igas) + reaction%eqgash2ostoich(igas)*rt_auxvar%ln_act_h2o
      endif

      do jcomp = 1, reaction%eqgasspecid(0,igas)
        comp_id = reaction%eqgasspecid(jcomp,igas)
        lnQKgas(igas) = lnQKgas(igas) + reaction%eqgasstoich(jcomp,igas)* &
                      log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
      enddo
      
      QKgas(igas) = exp(lnQKgas(igas))
          
      write(option%fid_out,133) reaction%gas_species_names(igas),lnQKgas(igas)*LN_TO_LOG, &
      QKgas(igas),reaction%eqgas_logK(igas)
    enddo
  endif

end subroutine ReactionPrintConstraint

! ************************************************************************** !
!
! DoubleLayer: Calculates double layer potential, surface charge, and
!              sorbed surface complex concentrations
! author: Peter C. Lichtner
! date: 10/28/08
!
! ************************************************************************** !
subroutine DoubleLayer(constraint_coupler,reaction,option)

  use Option_module
  use Input_module
  use String_module
  use Condition_module

  implicit none
  
  type(option_type) :: option
  type(tran_constraint_coupler_type) :: constraint_coupler
  type(reaction_type), pointer :: reaction

  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar

  PetscReal, parameter :: rgas = 8.3144621d0
  PetscReal, parameter :: tk = 273.15d0
  PetscReal, parameter :: epsilon = 78.5d0
  PetscReal, parameter :: epsilon0 = 8.854187817d-12
  PetscReal, parameter :: faraday = 96485.d0
  
  PetscReal :: fac, boltzmann, dbl_charge, surface_charge, ionic_strength, &
               charge_balance, potential, tempk, debye_length, &
               srfchrg_capacitance_model
               
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: srfcplx_conc(reaction%neqsrfcplx)

  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total

  PetscInt :: iphase
  PetscInt :: i, j, icomp, icplx, irxn, ncomp, ncplx

  PetscReal :: site_density(2)
  PetscReal :: mobile_fraction
  PetscInt :: num_types_of_sites
  PetscInt :: isite

  PetscBool :: one_more

    rt_auxvar => constraint_coupler%rt_auxvar
    global_auxvar => constraint_coupler%global_auxvar

    iphase = 1
    global_auxvar%temp(iphase) = option%reference_temperature
    tempk = tk + global_auxvar%temp(iphase)
    
    potential = 0.1d0 ! initial guess
    boltzmann = exp(-faraday*potential/(rgas*tempk))
        
    fac = sqrt(epsilon*epsilon0*rgas*tempk)
    
    ionic_strength = 0.d0
    charge_balance = 0.d0
    dbl_charge = 0.d0
    do icomp = 1, reaction%naqcomp      
      charge_balance = charge_balance + reaction%primary_spec_Z(icomp)* &
                       rt_auxvar%total(icomp,1)
                                        
      ionic_strength = ionic_strength + reaction%primary_spec_Z(icomp)**2* &
                       rt_auxvar%pri_molal(icomp)
      dbl_charge = dbl_charge + rt_auxvar%pri_molal(icomp)* &
                   (boltzmann**reaction%primary_spec_Z(icomp) - 1.d0)
    enddo
    
    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        ionic_strength = ionic_strength + reaction%eqcplx_Z(i)**2* &
                         rt_auxvar%sec_molal(i)
        dbl_charge = dbl_charge + rt_auxvar%sec_molal(i)* &
                     (boltzmann**reaction%eqcplx_Z(i) - 1.d0)
      enddo
    endif
    ionic_strength = 0.5d0*ionic_strength
    if (dbl_charge > 0.d0) then
      dbl_charge = fac*sqrt(2.d0*dbl_charge)
    else
      print *,'neg. dbl_charge: ',dbl_charge
      dbl_charge = fac*sqrt(2.d0*(-dbl_charge))
    endif
    
    srfchrg_capacitance_model = faraday*potential* &
      sqrt(2.d0*epsilon*epsilon0*ionic_strength/(rgas*tempk))
    
    surface_charge = 0.d0
    do irxn = 1, reaction%neqsrfcplxrxn
      ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
      do i = 1, ncplx
        icplx = reaction%eqsrfcplx_rxn_to_complex(i,irxn)
        surface_charge = surface_charge + reaction%eqsrfcplx_Z(icplx)* &
                         rt_auxvar%eqsrfcplx_conc(icplx)
      enddo
    enddo
    surface_charge = faraday*surface_charge
    
    debye_length = sqrt(fac/(2.d0*ionic_strength*1.d3))/faraday
    
    print *,'========================='
    print *,'dbl: debye_length = ',debye_length
    print *,'surface charge = ',dbl_charge,surface_charge, &
      srfchrg_capacitance_model
    print *,'ionic strength = ',ionic_strength
    print *,'chrg bal. = ',charge_balance,' Tk = ',tempk,' Boltz. = ',boltzmann
    print *,'srfcmplx: ',rt_auxvar%eqsrfcplx_conc
    print *,'========================='

!   compute surface complex concentrations  
    ln_conc = log(rt_auxvar%pri_molal)
    ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqsrfcplx_logKcoef, &
      reaction%eqsrfcplx_logK, &
      global_auxvar%temp(iphase),reaction%neqsrfcplx)
  endif
#endif  

  do irxn = 1, reaction%neqsrfcplxrxn
  
    ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
    
    free_site_conc = rt_auxvar%eqsrfcplx_free_site_conc(irxn)

    site_density(1) = reaction%eqsrfcplx_rxn_site_density(irxn)
    num_types_of_sites = 1
    
    do isite = 1, num_types_of_sites
      ! isite == 1 - immobile (colloids, minerals, etc.)
      ! isite == 2 - mobile (colloids)
    
      if (site_density(isite) < 1.d-40) cycle
    
      ! get a pointer to the first complex (there will always be at least 1)
      ! in order to grab free site conc
      one_more = PETSC_FALSE
      do
        total = free_site_conc
        ln_free_site = log(free_site_conc)
        
!       call srfcmplx(irxn,icplx,lnQK,reaction%eqsrfcplx_logK,reaction%eqsrfcplx_Z,potential, &
!               tempk,ln_act,rt_auxvar%ln_act_h2o,ln_free_site,srfcplx_conc)

#if 0
        do j = 1, ncplx
          icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
          
          ! compute ion activity product
          lnQK = -reaction%eqsrfcplx_logK(icplx)*LOG_TO_LN &
                 + reaction%eqsrfcplx_Z(icplx)*faraday*potential &
                 /(rgas*tempk)/LOG_TO_LN

          ! activity of water
          if (reaction%eqsrfcplxh2oid(icplx) > 0) then
            lnQK = lnQK + reaction%eqsrfcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
          endif

          lnQK = lnQK + reaction%eqsrfcplx_free_site_stoich(icplx)* &
                        ln_free_site
        
          ncomp = reaction%eqsrfcplxspecid(0,icplx)
          do i = 1, ncomp
            icomp = reaction%eqsrfcplxspecid(i,icplx)
            lnQK = lnQK + reaction%eqsrfcplxstoich(i,icplx)*ln_act(icomp)
          enddo
          srfcplx_conc(icplx) = exp(lnQK)
          
          total = total + reaction%eqsrfcplx_free_site_stoich(icplx)*srfcplx_conc(icplx) 
          
        enddo
#endif
        if (one_more) exit
        
        total = total / free_site_conc
        free_site_conc = site_density(isite) / total  
          
        one_more = PETSC_TRUE 

      enddo ! generic do
    enddo
  enddo
  
  print *,'srfcmplx1: ',srfcplx_conc

end subroutine DoubleLayer

#if 0
subroutine srfcmplx(irxn,icplx,lnQK,logK,Z,potential,tempk, &
                ln_act,ln_act_h2o,ln_free_site,srfcplx_conc)

implicit none

  PetscReal, parameter :: rgas = 8.3144621d0
  PetscReal, parameter :: tk = 273.15d0
  PetscReal, parameter :: faraday = 96485.d0
  
  PetscReal :: fac, boltzmann, dbl_charge, surface_charge, ionic_strength, &
               charge_balance, potential, tempk, debye_length, &
               srfchrg_capacitance_model
               
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: srfcplx_conc(reaction%neqsrfcplx)

  PetscReal :: free_site_conc
  PetscReal :: ln_free_site, ln_act_h2o
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total

  PetscInt :: i, j, icomp, icplx, irxn, ncomp, ncplx

        do j = 1, ncplx
          icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
          ! compute secondary species concentration
          lnQK = -logK(icplx)*LOG_TO_LN &
                 + Z(icplx)*faraday*potential &
                 /(rgas*tempk)/LOG_TO_LN

          ! activity of water
          if (reaction%eqsrfcplxh2oid(icplx) > 0) then
            lnQK = lnQK + reaction%eqsrfcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
          endif

          lnQK = lnQK + reaction%eqsrfcplx_free_site_stoich(icplx)* &
                        ln_free_site
        
          ncomp = reaction%eqsrfcplxspecid(0,icplx)
          do i = 1, ncomp
            icomp = reaction%eqsrfcplxspecid(i,icplx)
            lnQK = lnQK + reaction%eqsrfcplxstoich(i,icplx)*ln_act(icomp)
          enddo
          srfcplx_conc(icplx) = exp(lnQK)
        enddo
end subroutine srfcmplx
#endif

! ************************************************************************** !
!
! ReactionReadOutput: Reads species to be printed in output
! author: Glenn Hammond
! date: 01/24/09
!
! ************************************************************************** !
subroutine ReactionReadOutput(reaction,input,option)

  use Input_module
  use String_module  
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  PetscBool :: found
  PetscInt :: temp_int

  type(aq_species_type), pointer :: cur_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_type), pointer :: cur_mineral
  type(surface_complex_type), pointer :: cur_srfcplx
  type(surface_complexation_rxn_type), pointer :: cur_srfcplx_rxn
  
  nullify(cur_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_srfcplx)
  nullify(cur_srfcplx_rxn)
  
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,OUTPUT,SPECIES_NAME')
    
    word = name
    call StringToUpper(word)
    select case(word)
      case('OFF')
        reaction%print_all_species = PETSC_FALSE
        reaction%print_all_primary_species = PETSC_FALSE
        reaction%print_all_secondary_species = PETSC_FALSE
        reaction%print_all_gas_species = PETSC_FALSE
        reaction%print_all_mineral_species = PETSC_FALSE
        reaction%print_pH = PETSC_FALSE
        reaction%print_kd = PETSC_FALSE
        reaction%print_total_sorb = PETSC_FALSE
        reaction%print_total_sorb_mobile = PETSC_FALSE
        reaction%print_colloid = PETSC_FALSE
        reaction%print_act_coefs = PETSC_FALSE
        reaction%print_total_component = PETSC_TRUE
        reaction%print_free_ion = PETSC_FALSE
      case('ALL')
        reaction%print_all_species = PETSC_TRUE
        reaction%print_all_primary_species = PETSC_TRUE
 !       reaction%print_all_secondary_species = PETSC_TRUE
 !       reaction%print_all_gas_species = PETSC_TRUE
        reaction%print_all_mineral_species = PETSC_TRUE
        reaction%print_pH = PETSC_TRUE
      case('PRIMARY_SPECIES')
        reaction%print_all_primary_species = PETSC_TRUE
        reaction%print_pH = PETSC_TRUE
      case('SECONDARY_SPECIES')
        reaction%print_all_secondary_species = PETSC_TRUE
      case('GASES')
        reaction%print_all_gas_species = PETSC_TRUE
      case('MINERALS')
        reaction%print_all_mineral_species = PETSC_TRUE
      case('PH')
        reaction%print_pH = PETSC_TRUE
      case('KD')
        reaction%print_kd = PETSC_TRUE
      case('COLLOIDS')
        reaction%print_colloid = PETSC_TRUE
      case('TOTAL_SORBED')
        reaction%print_total_sorb = PETSC_TRUE
      case('TOTAL_SORBED_MOBILE')
        reaction%print_total_sorb_mobile = PETSC_TRUE
      case('FREE_ION')
        reaction%print_total_component = PETSC_FALSE
        reaction%print_free_ion = PETSC_TRUE
      case('ACTIVITY_COEFFICIENTS')
        reaction%print_act_coefs = PETSC_TRUE
      case('MOLARITY')
        reaction%print_free_conc_type = PRIMARY_MOLARITY
        reaction%print_tot_conc_type = TOTAL_MOLARITY
      case('MOLALITY')
        reaction%print_free_conc_type = PRIMARY_MOLALITY
        reaction%print_tot_conc_type = TOTAL_MOLALITY
      case('AGE')
        reaction%print_age = PETSC_TRUE
      case default        
        found = PETSC_FALSE
        if (.not.found) then
          cur_aq_spec => reaction%primary_species_list
          do
            if (.not.associated(cur_aq_spec)) exit
            if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
              cur_aq_spec%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_aq_spec => cur_aq_spec%next
          enddo
        endif
        if (.not.found) then
          cur_aq_spec => reaction%secondary_species_list
          do
            if (.not.associated(cur_aq_spec)) exit
            if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
              cur_aq_spec%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_aq_spec => cur_aq_spec%next
          enddo  
        endif
        if (.not.found) then
          cur_gas_spec => reaction%gas_species_list
          do
            if (.not.associated(cur_gas_spec)) exit
            if (StringCompare(name,cur_gas_spec%name,MAXWORDLENGTH)) then
              cur_gas_spec%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_gas_spec => cur_gas_spec%next
          enddo  
        endif
        if (.not.found) then
          cur_mineral => reaction%mineral_list
          do
            if (.not.associated(cur_mineral)) exit
            if (StringCompare(name,cur_mineral%name,MAXWORDLENGTH)) then
              cur_mineral%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_mineral => cur_mineral%next
          enddo
        endif
        if (.not.found) then
          cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
          do
            if (.not.associated(cur_srfcplx_rxn)) exit
            if (StringCompare(name,cur_srfcplx_rxn%free_site_name,MAXWORDLENGTH)) then
              cur_srfcplx_rxn%free_site_print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            if (.not.found) then
              cur_srfcplx => cur_srfcplx_rxn%complex_list
              do  
                if (.not.associated(cur_srfcplx)) exit
              if (StringCompare(name,cur_srfcplx%name,MAXWORDLENGTH)) then
                cur_srfcplx%print_me = PETSC_TRUE
                found = PETSC_TRUE
                exit
              endif
                cur_srfcplx => cur_srfcplx%next
              enddo
            endif
            cur_srfcplx_rxn => cur_srfcplx_rxn%next
          enddo  
        endif

        if (.not.found) then
          option%io_buffer = 'CHEMISTRY,OUTPUT species name: '//trim(name)// &
                             ' not found among chemical species'
          call printErrMsg(option)
        endif
    end select

  enddo

end subroutine ReactionReadOutput

! ************************************************************************** !
!
! RJumpStartKineticSorption: Calculates the concentrations of species sorbing
!                            through kinetic sorption processes based
!                            on equilibrium with the aqueous phase.
! author: Glenn Hammond
! date: 08/05/09
!
! ************************************************************************** !
subroutine RJumpStartKineticSorption(rt_auxvar,global_auxvar, &
                                     reaction,option)

  use Option_module
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option

  PetscInt :: irate
  
  ! WARNING: below assumes site concentration multiplicative factor
  allocate(rt_auxvar%dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp))
  rt_auxvar%dtotal_sorb_eq = 0.d0
  call RTotalSorbEqSurfCplx(rt_auxvar,global_auxvar,reaction,option)
  do irate = 1, reaction%kinmr_nrate
    rt_auxvar%kinmr_total_sorb(:,irate) = reaction%kinmr_frac(irate) * &
                                          rt_auxvar%total_sorb_eq
  enddo
  deallocate(rt_auxvar%dtotal_sorb_eq)
  nullify(rt_auxvar%dtotal_sorb_eq)

end subroutine RJumpStartKineticSorption

! ************************************************************************** !
!
! RReact: Solves reaction portion of operator splitting using Newton-Raphson
! author: Glenn Hammond
! date: 05/04/10
!
! ************************************************************************** !
subroutine RReact(rt_auxvar,global_auxvar,total,volume,porosity, &
                  num_iterations_,reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: total(reaction%ncomp)
  type(option_type) :: option
  PetscReal :: volume
  PetscReal :: porosity
  PetscInt :: num_iterations_
  PetscReal :: sign_(reaction%ncomp)
  
  PetscReal :: residual(reaction%ncomp)
  PetscReal :: res(reaction%ncomp)
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  PetscReal :: one_over_dt
  PetscReal :: prev_molal(reaction%ncomp)
  PetscReal :: update(reaction%ncomp)
  PetscReal :: maximum_relative_change
  PetscReal :: accumulation_coef
  PetscReal :: fixed_accum(reaction%ncomp)
  PetscInt :: num_iterations
  PetscInt :: icomp
  PetscReal :: ratio, min_ratio
  
  PetscInt, parameter :: iphase = 1

  one_over_dt = 1.d0/option%tran_dt
  num_iterations = 0

  ! calculate fixed portion of accumulation term
  ! fixed_accum is overwritten in RTAccumulation
  ! Since RTAccumulation uses rt_auxvar%total, we must overwrite the 
  ! rt_auxvar total variables
  ! aqueous
  rt_auxvar%total(:,iphase) = total(1:reaction%naqcomp)

! skip chemistry if species nonreacting 
#if 1  
  if (.not.reaction%use_full_geochemistry) then
    rt_auxvar%pri_molal(:) = total(:)/global_auxvar%den_kg(iphase)*1.d3
    return
  endif
#endif  
  
  ! still need code to overwrite other phases
  call RTAccumulation(rt_auxvar,global_auxvar,porosity,volume,reaction, &
                      option,fixed_accum)
  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    call RAccumulationSorb(rt_auxvar,global_auxvar,volume,reaction, &
                           option,fixed_accum)  
  endif

  ! now update activity coefficients
  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    call RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
  endif
  
  do
  
    num_iterations = num_iterations + 1

    if (reaction%act_coef_update_frequency == ACT_COEF_FREQUENCY_NEWTON_ITER) then
      call RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
    endif
    call RTAuxVarCompute(rt_auxvar,global_auxvar,reaction,option)
    
    ! Accumulation
    ! residual is overwritten in RTAccumulation()
    call RTAccumulation(rt_auxvar,global_auxvar,porosity,volume,reaction, &
                        option,residual)
    residual = residual-fixed_accum

    ! J is overwritten in RTAccumulationDerivative()
    call RTAccumulationDerivative(rt_auxvar,global_auxvar,porosity,volume, &
                                  reaction,option,J)

    if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
      call RAccumulationSorb(rt_auxvar,global_auxvar,volume,reaction, &
                             option,residual)
      call RAccumulationSorbDerivative(rt_auxvar,global_auxvar,volume, &
                                       reaction,option,J)
    endif

                         ! derivative
    call RReaction(residual,J,PETSC_TRUE,rt_auxvar,global_auxvar, porosity, volume, &
                   reaction,option)
    
    if (maxval(abs(residual)) < reaction%max_residual_tolerance) exit


    call RSolve(residual,J,rt_auxvar%pri_molal,update,reaction%ncomp, &
                reaction%use_log_formulation)
    
    prev_molal = rt_auxvar%pri_molal

    if (reaction%use_log_formulation) then
      update = dsign(1.d0,update)*min(dabs(update),reaction%max_dlnC)
      rt_auxvar%pri_molal = rt_auxvar%pri_molal*exp(-update)    
    else ! linear upage
      ! ensure non-negative concentration
      min_ratio = 1.d20 ! large number
      do icomp = 1, reaction%ncomp
        if (prev_molal(icomp) <= update(icomp)) then
          ratio = abs(prev_molal(icomp)/update(icomp))
          if (ratio < min_ratio) min_ratio = ratio
        endif
      enddo
      if (min_ratio < 1.d0) then
        ! scale by 0.99 to make the update slightly smaller than the min_ratio
        update = update*min_ratio*0.99d0
      endif
      rt_auxvar%pri_molal = prev_molal - update
    endif
    

    maximum_relative_change = maxval(abs((rt_auxvar%pri_molal-prev_molal)/ &
                                         prev_molal))
    
    if (maximum_relative_change < reaction%max_relative_change_tolerance) exit

    if (num_iterations > 50) then
      if (num_iterations > 50) then
        rt_auxvar%pri_molal(:) = (rt_auxvar%pri_molal(:)-prev_molal(:))* &
                                 0.1d0+prev_molal(:)
      else if (num_iterations > 100) then
        rt_auxvar%pri_molal(:) = (rt_auxvar%pri_molal(:)-prev_molal(:))* &
                                 0.01d0+prev_molal(:)
      else if (num_iterations > 150) then
        rt_auxvar%pri_molal(:) = (rt_auxvar%pri_molal(:)-prev_molal(:))* &
                                 0.001d0+prev_molal(:)
      else if (num_iterations > 500) then
        print *, 'Maximum iterations in RReact: stop: ',num_iterations
        print *, 'Maximum iterations in RReact: residual: ',residual
        print *, 'Maximum iterations in RReact: primary species: ',rt_auxvar%pri_molal
        stop
      endif

    endif
  
  enddo

  ! one last update
  call RTAuxVarCompute(rt_auxvar,global_auxvar,reaction,option)

  num_iterations_ = num_iterations
  
end subroutine RReact
      
! ************************************************************************** !
!
! RReaction: Computes reactions
! author: Glenn Hammond
! date: 09/30/08
!
! ************************************************************************** !
subroutine RReaction(Res,Jac,derivative,rt_auxvar,global_auxvar,porosity, &
                     volume,reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscBool :: derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume

  if (reaction%nkinmnrl > 0) then
    call RKineticMineral(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                         volume,reaction,option)
  endif
  
  if (reaction%kinmr_nrate > 0) then
    call RMultiRateSorption(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                            volume,reaction,option)
  endif
  
  if (reaction%nkinsrfcplxrxn > 0) then
    call RKineticSurfCplx(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                          volume,reaction,option)
  endif
  
  if (reaction%ngeneral_rxn > 0) then
    call RGeneral(Res,Jac,derivative,rt_auxvar,global_auxvar,porosity, &
                  volume,reaction,option)
  endif
  
  ! add new reactions here

end subroutine RReaction

! ************************************************************************** !
!
! RReaction: Computes reactions
! author: Glenn Hammond
! date: 09/30/08
!
! ************************************************************************** !
subroutine RReactionDerivative(Res,Jac,rt_auxvar,global_auxvar,porosity, &
                               volume,reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(reactive_transport_auxvar_type) :: rt_auxvar_pert
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
   
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp
  PetscReal :: Jac_dummy(reaction%ncomp,reaction%ncomp)
  PetscReal :: pert
  PetscBool :: compute_derivative

  ! add new reactions in the 3 locations below

  if (.not.option%numerical_derivatives) then ! analytical derivative
  !if (PETSC_FALSE) then
    compute_derivative = PETSC_TRUE
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res,Jac,compute_derivative,rt_auxvar, &
                           global_auxvar,volume,reaction,option)
    endif
    if (reaction%kinmr_nrate > 0) then
      call RMultiRateSorption(Res,Jac,compute_derivative,rt_auxvar, &
                              global_auxvar,volume,reaction,option)
    endif
    if (reaction%nkinsrfcplxrxn > 0) then
      call RKineticSurfCplx(Res,Jac,compute_derivative,rt_auxvar, &
                            global_auxvar,volume,reaction,option)
    endif
    if (reaction%ngeneral_rxn > 0) then
      call RGeneral(Res,Jac,compute_derivative,rt_auxvar, &
                    global_auxvar,porosity,volume,reaction,option)
    endif    

    ! #1: add new reactions here

  else ! numerical derivative
    compute_derivative = PETSC_FALSE
    Res_orig = 0.d0
    option%iflag = 0 ! be sure not to allocate mass_balance array
    call RTAuxVarInit(rt_auxvar_pert,reaction,option)
    call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                           global_auxvar,volume,reaction,option)
    endif
    if (reaction%kinmr_nrate > 0) then
      call RMultiRateSorption(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                              global_auxvar,volume,reaction,option)
    endif     
    if (reaction%nkinsrfcplxrxn > 0) then
      call RKineticSurfCplx(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                            global_auxvar,volume,reaction,option)
    endif
    if (reaction%ngeneral_rxn > 0) then
      call RGeneral(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                    global_auxvar,porosity,volume,reaction,option)
    endif    

    ! #2: add new reactions here

    do jcomp = 1, reaction%ncomp
      Res_pert = 0.d0
      call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
      pert = rt_auxvar_pert%pri_molal(jcomp)*perturbation_tolerance
      rt_auxvar_pert%pri_molal(jcomp) = rt_auxvar_pert%pri_molal(jcomp) + pert
      
      call RTotal(rt_auxvar_pert,global_auxvar,reaction,option)
      if (reaction%neqsorb > 0) call RTotalSorb(rt_auxvar_pert,global_auxvar, &
                                              reaction,option)

      if (reaction%nkinmnrl > 0) then
        call RKineticMineral(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                             global_auxvar,volume,reaction,option)
      endif
      if (reaction%kinmr_nrate > 0) then
        call RMultiRateSorption(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                                global_auxvar,volume,reaction,option)
      endif      
      if (reaction%nkinsrfcplxrxn > 0) then
        call RKineticSurfCplx(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                              global_auxvar,volume,reaction,option)
      endif
      if (reaction%ngeneral_rxn > 0) then
        call RGeneral(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                      global_auxvar,porosity,volume,reaction,option)
      endif  
      
      ! #3: add new reactions here

      do icomp = 1, reaction%ncomp
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + (Res_pert(icomp)-Res_orig(icomp))/pert
      enddo
    enddo
    do icomp = 1, reaction%ncomp
      do jcomp = 1, reaction%ncomp
        if (dabs(Jac(icomp,jcomp)) < 1.d-40)  Jac(icomp,jcomp) = 0.d0
      enddo
    enddo
    call RTAuxVarDestroy(rt_auxvar_pert)
  endif

end subroutine RReactionDerivative
                               
! ************************************************************************** !
!
! CO2AqActCoeff: Computes activity coefficients of aqueous CO2
! author: Chuan Lu
! date: 07/13/09
!
! ************************************************************************** !
subroutine CO2AqActCoeff(rt_auxvar,global_auxvar,reaction,option)
    
  use Option_module
  use co2eos_module
 
  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
   
  PetscReal :: m_na, m_cl, tc, co2aqact, lngamco2, henry, xphico2, pco2
  PetscReal :: sat_pressure
  PetscErrorCode :: ierr 

  tc = global_auxvar%temp(1)
  pco2 = global_auxvar%pres(2)
  sat_pressure =0D0

  m_na = option%m_nacl; m_cl = m_na
  if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
     m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
     m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
  endif

  call Henry_duan_sun(tc,pco2*1D-5,henry, 1.D0,lngamco2, &
         m_na,m_cl,sat_pressure*1D-5, co2aqact)
         
  rt_auxvar%pri_act_coef(reaction%species_idx%co2_aq_id) = co2aqact 
 ! print *, 'CO2AqActCoeff', tc, pco2, m_na,m_cl, sat_pressure,co2aqact
end subroutine CO2AqActCoeff

! ************************************************************************** !
!
! RActivityCoefficients: Computes the ionic strength and activity coefficients
! author: Glenn Hammond
! date: 09/30/08
!
! ************************************************************************** !
subroutine RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)

  use Option_module
  
  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: icplx, icomp, it, j, jcomp, ncomp
  PetscReal :: I, sqrt_I, II, sqrt_II, f, fpri, didi, dcdi, den, dgamdi, &
    lnQK, sum, sum_pri_molal, sum_sec_molal
  PetscReal :: sum_molality
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)

  if (reaction%use_activity_h2o) then
    sum_pri_molal = 0.d0
    do j = 1, reaction%naqcomp
      sum_pri_molal = sum_pri_molal + rt_auxvar%pri_molal(j)
    enddo
  endif

  if (reaction%act_coef_update_algorithm == ACT_COEF_ALGORITHM_NEWTON) then

    ln_conc = log(rt_auxvar%pri_molal)
    ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  
#ifdef TEMP_DEPENDENT_LOGK
    if (.not.option%use_isothermal) then
      call ReactionInterpolateLogK(reaction%eqcplx_logKcoef,reaction%eqcplx_logK, &
                               global_auxvar%temp(1),reaction%neqcplx)
    endif
#endif  
  
  ! compute primary species contribution to ionic strength
    fpri = 0.d0
    sum_molality = 0.d0
    do j = 1, reaction%naqcomp
      fpri = fpri + rt_auxvar%pri_molal(j)*reaction%primary_spec_Z(j)* &
                                         reaction%primary_spec_Z(j)
    enddo
  
    it = 0
    II = 0
    do
      it = it + 1
      
      if (it > 50) then
        print *,' too many iterations in computing activity coefficients-stop',it,f,I
        stop
      endif
    
  ! add secondary species contribution to ionic strength
      I = fpri
      do icplx = 1, reaction%neqcplx ! for each secondary species
        I = I + rt_auxvar%sec_molal(icplx)*reaction%eqcplx_Z(icplx)* &
                                         reaction%eqcplx_Z(icplx)
      enddo
      I = 0.5d0*I
      f = I
    
      if (abs(I-II) < 1.d-6*I) exit
    
      if (reaction%neqcplx > 0) then
        didi = 0.d0
        sqrt_I = sqrt(I)
        do icplx = 1, reaction%neqcplx
          if (abs(reaction%eqcplx_Z(icplx)) > 0.d0) then
            sum = 0.5d0*reaction%debyeA*reaction%eqcplx_Z(icplx)* &
            reaction%eqcplx_Z(icplx) &
            /(sqrt_I*(1.d0+reaction%debyeB*reaction%eqcplx_a0(icplx)*sqrt_I)**2) &
            -reaction%debyeBdot
            ncomp = reaction%eqcplxspecid(0,icplx)
            do jcomp = 1, ncomp
              j = reaction%eqcplxspecid(jcomp,icplx)
              if(abs(reaction%primary_spec_Z(j)) > 0.d0) then
                dgamdi = -0.5d0*reaction%debyeA*reaction%primary_spec_Z(j)**2/(sqrt_I* &
                (1.d0+reaction%debyeB*reaction%primary_spec_a0(j)*sqrt_I)**2)+ &
                reaction%debyeBdot 
                sum = sum + reaction%eqcplxstoich(jcomp,icplx)*dgamdi
              endif
            enddo
            dcdi = rt_auxvar%sec_molal(icplx)*LOG_TO_LN*sum
            didi = didi+0.5d0*reaction%eqcplx_Z(icplx)*reaction%eqcplx_Z(icplx)*dcdi
          endif
        enddo
        den = 1.d0-didi
        if (abs(den) > 0.d0) then
          II = (f-I*didi)/den
        else
          II = f
        endif
      else
        II = f
      endif
    
      if (II < 0.d0) then
        write(option%io_buffer,*) 'ionic strength negative! it =',it, &
          ' I= ',I,II,den,didi,dcdi,sum
        call printErrMsg(option)        
      endif
    
  ! compute activity coefficients
  ! primary species
      I = II
      sqrt_I = sqrt(I)
      do icomp = 1, reaction%naqcomp
        if (abs(reaction%primary_spec_Z(icomp)) > 0.d0) then
          rt_auxvar%pri_act_coef(icomp) = exp((-reaction%primary_spec_Z(icomp)* &
                                        reaction%primary_spec_Z(icomp)* &
                                        sqrt_I*reaction%debyeA/ &
                                        (1.d0+reaction%primary_spec_a0(icomp)* &
                                        reaction%debyeB*sqrt_I)+ &
                                        reaction%debyeBdot*I)* &
                                        LOG_TO_LN)
        else
          rt_auxvar%pri_act_coef(icomp) = 1.d0
        endif
      enddo
                
  ! secondary species
      sum_sec_molal = 0.d0
      do icplx = 1, reaction%neqcplx
        if (abs(reaction%eqcplx_Z(icplx)) > 0.d0) then
          rt_auxvar%sec_act_coef(icplx) = exp((-reaction%eqcplx_Z(icplx)* &
                                        reaction%eqcplx_Z(icplx)* &
                                        sqrt_I*reaction%debyeA/ &
                                        (1.d0+reaction%eqcplx_a0(icplx)* &
                                        reaction%debyeB*sqrt_I)+ &
                                        reaction%debyeBdot*I)* &
                                        LOG_TO_LN)
        else
          rt_auxvar%sec_act_coef(icplx) = 1.d0
        endif
    
    ! compute secondary species concentration
        lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

    ! activity of water
        if (reaction%eqcplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
        endif

        ncomp = reaction%eqcplxspecid(0,icplx)
        do jcomp = 1, ncomp
          icomp = reaction%eqcplxspecid(jcomp,icplx)
          lnQK = lnQK + reaction%eqcplxstoich(jcomp,icplx)*ln_act(icomp)
        enddo
        rt_auxvar%sec_molal(icplx) = exp(lnQK)/rt_auxvar%sec_act_coef(icplx)
        sum_sec_molal = sum_sec_molal + rt_auxvar%sec_molal(icplx)
      
      enddo
      
      if (reaction%use_activity_h2o) then
        rt_auxvar%ln_act_h2o = 1.d0-0.017d0*(sum_pri_molal+sum_sec_molal)
        if (rt_auxvar%ln_act_h2o > 0.d0) then
          rt_auxvar%ln_act_h2o = log(rt_auxvar%ln_act_h2o)
        else
          rt_auxvar%ln_act_h2o = 0.d0
          write(option%io_buffer,*) 'activity of H2O negative! ln act H2O =', &
            rt_auxvar%ln_act_h2o
          call printMsg(option)
        endif
      endif

    enddo
  
  else
  
  ! compute ionic strength
  ! primary species
    I = 0.d0
    do icomp = 1, reaction%naqcomp
      I = I + rt_auxvar%pri_molal(icomp)*reaction%primary_spec_Z(icomp)* &
                                       reaction%primary_spec_Z(icomp)
    enddo
  
  ! secondary species
    do icplx = 1, reaction%neqcplx ! for each secondary species
      I = I + rt_auxvar%sec_molal(icplx)*reaction%eqcplx_Z(icplx)* &
                                       reaction%eqcplx_Z(icplx)
    enddo
    I = 0.5d0*I
    sqrt_I = sqrt(I)
  
  ! compute activity coefficients
  ! primary species
    do icomp = 1, reaction%naqcomp
      if (abs(reaction%primary_spec_Z(icomp)) > 1.d-10) then
        rt_auxvar%pri_act_coef(icomp) = exp((-reaction%primary_spec_Z(icomp)* &
                                      reaction%primary_spec_Z(icomp)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%primary_spec_a0(icomp)* &
                                      reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                      LOG_TO_LN)
      else
        rt_auxvar%pri_act_coef(icomp) = 1.d0
      endif
    enddo
                
  ! secondary species
    sum_sec_molal = 0.d0
    do icplx = 1, reaction%neqcplx
      if (dabs(reaction%eqcplx_Z(icplx)) > 1.d-10) then
        rt_auxvar%sec_act_coef(icplx) = exp((-reaction%eqcplx_Z(icplx)* &
                                      reaction%eqcplx_Z(icplx)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%eqcplx_a0(icplx)* &
                                      reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                      LOG_TO_LN)
      else
        rt_auxvar%sec_act_coef(icplx) = 1.d0
      endif
      sum_sec_molal = sum_sec_molal + rt_auxvar%sec_molal(icplx)
    enddo
    
    if (reaction%use_activity_h2o) then
      rt_auxvar%ln_act_h2o = 1.d0-0.017d0*(sum_pri_molal+sum_sec_molal)
      if (rt_auxvar%ln_act_h2o > 0.d0) then
        rt_auxvar%ln_act_h2o = log(rt_auxvar%ln_act_h2o)
      else
        rt_auxvar%ln_act_h2o = 0.d0
      endif
    endif
  endif
  
end subroutine RActivityCoefficients

! ************************************************************************** !
!
! RTotal: Computes the total component concentrations and derivative with
!         respect to free-ion
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine RTotal(rt_auxvar,global_auxvar,reaction,option)

  use Option_module
  use co2eos_module, only: Henry_duan_sun
  use water_eos_module
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, icplx, icomp, jcomp, iphase, ncomp, ieqgas
  PetscErrorCode :: ierr
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal
  PetscReal :: den_kg_per_L, xmass
  PetscReal :: pressure, temperature, xphico2, muco2, den, m_na, m_cl
  
#ifdef CHUAN_CO2  
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
  rt_auxvar%total = 0.d0 !debugging 
#endif
  
  iphase = 1           
!  den_kg_per_L = global_auxvar%den_kg(iphase)*1.d-3              
  xmass = 1.d0
  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  den_kg_per_L = global_auxvar%den_kg(iphase)*xmass*1.d-3

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  rt_auxvar%total(:,iphase) = rt_auxvar%pri_molal(:)
  
  ! initialize derivatives
  rt_auxvar%aqueous%dtotal = 0.d0
  do icomp = 1, reaction%naqcomp
    rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) = 1.d0
  enddo
  
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal .and. reaction%neqcplx > 0) then
    call ReactionInterpolateLogK(reaction%eqcplx_logKcoef,reaction%eqcplx_logK, &
                                 global_auxvar%temp(iphase),reaction%neqcplx)
  endif
#endif  
  
  do icplx = 1, reaction%neqcplx ! for each secondary species
    ! compute secondary species concentration
    lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

    ! activity of water
    if (reaction%eqcplxh2oid(icplx) > 0) then
      lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
    endif

    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      lnQK = lnQK + reaction%eqcplxstoich(i,icplx)*ln_act(icomp)
    enddo
    rt_auxvar%sec_molal(icplx) = exp(lnQK)/rt_auxvar%sec_act_coef(icplx)
  
    ! add contribution to primary totals
    ! units of total = mol/L
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
                                      reaction%eqcplxstoich(i,icplx)* &
                                      rt_auxvar%sec_molal(icplx)
    enddo
    
    ! add contribution to derivatives of total with respect to free
    ! bear in mind that the water density portion is scaled below
    do j = 1, ncomp
      jcomp = reaction%eqcplxspecid(j,icplx)
      tempreal = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 rt_auxvar%sec_act_coef(icplx)
      do i = 1, ncomp
        icomp = reaction%eqcplxspecid(i,icplx)
        rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) = rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) + &
                                               reaction%eqcplxstoich(i,icplx)*tempreal
      enddo
    enddo
  enddo

  ! convert molality -> molarity
  ! unit of total = mol/L water
  rt_auxvar%total(:,iphase) = rt_auxvar%total(:,iphase)*den_kg_per_L
  
  ! units of dtotal = kg water/L water
  rt_auxvar%aqueous%dtotal = rt_auxvar%aqueous%dtotal*den_kg_per_L
 !*********** Add SC phase contribution ***************************  
#ifdef CHUAN_CO2

  iphase = 2           
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal .and. reaction%ngas > 0) then
    call ReactionInterpolateLogK(reaction%eqgas_logKcoef,reaction%eqgas_logK, &
                                 global_auxvar%temp(1),reaction%ngas)
  endif
#endif  

  if(iphase > option%nphase) return 
  rt_auxvar%total(:,iphase) = 0D0
  rt_auxvar%aqueous%dtotal(:,:,iphase)=0D0
!  do icomp = 1, reaction%naqcomp
!    rt_auxvar%dtotal(icomp,icomp,iphase) = 1.d0
!  enddo
    
!  den_kg_per_L = global_auxvar%den_kg(iphase)*1.d-3     
  if(global_auxvar%sat(iphase)>1D-20)then
    do ieqgas = 1, reaction%ngas ! all gas phase species are secondary
   
      pressure = global_auxvar%pres(2)
      temperature = global_auxvar%temp(1)
      xphico2 = global_auxvar%fugacoeff(1)
      den = global_auxvar%den(2)
 
      call PSAT(temperature, sat_pressure, ierr)
      pco2 = pressure - sat_pressure
!     call co2_span_wagner(pressure*1.D-6,temperature+273.15D0,dg,dddt,dddp,fg, &
!              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
!
!            fg = fg*1D6
!            xphico2 = fg / pco2
!            global_auxvar%fugacoeff(1) = xphico2


      if(abs(reaction%species_idx%co2_gas_id) == ieqgas )then
!          call Henry_duan_sun_0NaCl(pco2*1D-5, temperature, henry)
        if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
          m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
          m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
          call Henry_duan_sun(temperature,pressure*1D-5,muco2,xphico2, &
                lngamco2,m_na,m_cl,sat_pressure*1D-5)
        else
          call Henry_duan_sun(temperature,pressure*1D-5,muco2,xphico2, &
                lngamco2,option%m_nacl,option%m_nacl,sat_pressure*1D-5)
        endif
        !lnQk = - log(muco2) 
        lnQk = - log(muco2)-lngamco2
           
      else   
        lnQK = -reaction%eqgas_logK(ieqgas)*LOG_TO_LN
      endif 
          
      if (reaction%eqgash2oid(ieqgas) > 0) then
        lnQK = lnQK + reaction%eqgash2ostoich(ieqgas)*rt_auxvar%ln_act_h2o
!       print *,'Ttotal', reaction%eqgash2ostoich(ieqgas), rt_auxvar%ln_act_h2o
      endif
   
   ! contribute to %total          
   !     do i = 1, ncomp
   ! removed loop over species, suppose only one primary species is related
      icomp = reaction%eqgasspecid(1,ieqgas)
      pressure =pressure *1D-5
        
      rt_auxvar%gas_molal(ieqgas) = &
          exp(lnQK+lngamco2)*rt_auxvar%pri_molal(icomp)&
!          rt_auxvar%pri_act_coef(icomp)*exp(lnQK)*rt_auxvar%pri_molal(icomp)&
          /pressure /xphico2* den
      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
                                        reaction%eqgasstoich(1,ieqgas)* &
                                        rt_auxvar%gas_molal(ieqgas)
!       print *,'Ttotal',pressure, temperature, xphico2, den, lnQk,rt_auxvar%pri_molal(icomp),&
!        global_auxvar%sat(2),rt_auxvar%gas_molal(ieqgas)
   !     if(rt_auxvar%total(icomp,iphase) > den)rt_auxvar%total(icomp,iphase) = den* .99D0
   !     enddo

   ! contribute to %dtotal
   !      tempreal = exp(lnQK+lngamco2)/pressure /xphico2* den 
      tempreal = rt_auxvar%pri_act_coef(icomp)*exp(lnQK)/pressure /xphico2* den 
      rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) = rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) + &
                                               reaction%eqgasstoich(1,ieqgas)*tempreal
    
    enddo
   ! rt_auxvar%total(:,iphase) = rt_auxvar%total(:,iphase)!*den_kg_per_L
    ! units of dtotal = kg water/L water
   ! rt_auxvar%dtotal(:, :,iphase) = rt_auxvar%dtotal(:,:,iphase)!*den_kg_per_L
   endif   
  
#endif  
end subroutine RTotal

! ************************************************************************** !
!
! RTotalSorb: Computes the total sorbed component concentrations and 
!             derivative with respect to free-ion
! author: Glenn Hammond
! date: 10/22/08
!
! ************************************************************************** !
subroutine RTotalSorb(rt_auxvar,global_auxvar,reaction,option)

  use Option_module
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  ! initialize total sorbed concentrations and derivatives
  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    rt_auxvar%total_sorb_eq = 0.d0
    rt_auxvar%dtotal_sorb_eq = 0.d0  
  endif

  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    call RTotalSorbEqSurfCplx(rt_auxvar,global_auxvar,reaction,option)
  endif
  
  if (reaction%neqionxrxn > 0) then
    call RTotalSorbEqIonx(rt_auxvar,global_auxvar,reaction,option)
  endif
  
  if (reaction%neqkdrxn > 0) then
    call RTotalSorbKD(rt_auxvar,global_auxvar,reaction,option)
  endif
  
end subroutine RTotalSorb

! ************************************************************************** !
!
! RTotalSorbKD: Computes the total sorbed component concentrations and 
!               derivative with respect to free-ion for the linear 
!               K_D model
! author: Glenn Hammond
! date: 09/30/2010
!
! ************************************************************************** !
subroutine RTotalSorbKD(rt_auxvar,global_auxvar,reaction,option)

  use Option_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: irxn
  PetscInt :: icomp
  PetscReal :: res
  PetscReal :: dres_dc
  PetscReal :: activity
  PetscReal :: molality
  PetscReal :: tempreal
  PetscReal :: one_over_n
  PetscReal :: activity_one_over_n

  ! Surface Complexation
  do irxn = 1, reaction%neqkdrxn
    icomp = reaction%eqkdspecid(irxn)
    molality = rt_auxvar%pri_molal(icomp)
    activity = molality*rt_auxvar%pri_act_coef(icomp)
    select case(reaction%eqkdtype(irxn))
      case(SORPTION_LINEAR)
        ! Csorb = Kd*Caq
       res = reaction%eqkddistcoef(irxn)*activity
        dres_dc = res/molality
      case(SORPTION_LANGMUIR)
        ! Csorb = K*Caq*b/(1+K*Caq)
        tempreal = reaction%eqkddistcoef(irxn)*activity
        res = tempreal*reaction%eqkdlangmuirb(irxn) / (1.d0 + tempreal)
        dres_dc = res/molality - &
                  res / (1.d0 + tempreal) * tempreal / molality
      case(SORPTION_FREUNDLICH)
        ! Csorb = Kd*Caq**(1/n)
        one_over_n = 1.d0/reaction%eqkdfreundlichn(irxn)
        activity_one_over_n = activity**one_over_n
        res = reaction%eqkddistcoef(irxn)* &
                activity**one_over_n
        dres_dc = res/molality*one_over_n
      case default
        res = 0.d0
        dres_dc = 0.d0
    end select
    rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + res
    rt_auxvar%dtotal_sorb_eq(icomp,icomp) = &
      rt_auxvar%dtotal_sorb_eq(icomp,icomp) + dres_dc 
  enddo

end subroutine RTotalSorbKD

! ************************************************************************** !
!
! RTotalSorbEqSurfCplx: Computes the total sorbed component concentrations and 
!                       derivative with respect to free-ion for equilibrium 
!                       surface complexation
! author: Glenn Hammond
! date: 10/22/08; 05/26/09
!
! ************************************************************************** !
subroutine RTotalSorbEqSurfCplx(rt_auxvar,global_auxvar,reaction,option)

  use Option_module
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, ncomp, ncplx
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: srfcplx_conc(reaction%neqsrfcplx)
  PetscReal :: dSx_dmi(reaction%naqcomp)
  PetscReal :: nui_Si_over_Sx
  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: tol = 1.d-12
  PetscBool :: one_more
  PetscReal :: res, dres_dfree_site, dfree_site_conc
  PetscReal :: site_density(2)
  PetscReal :: mobile_fraction
  PetscInt :: num_types_of_sites
  PetscInt :: isite
  
  if (reaction%ncollcomp > 0) then  
    rt_auxvar%colloid%total_eq_mob = 0.d0
    rt_auxvar%colloid%dRj_dCj%dtotal = 0.d0
  endif
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqsrfcplx_logKcoef,reaction%eqsrfcplx_logK, &
                               global_auxvar%temp(iphase),reaction%neqsrfcplx)
  endif
#endif  

  ! Surface Complexation
  do irxn = 1, reaction%neqsrfcplxrxn
  
    ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
    
    free_site_conc = rt_auxvar%eqsrfcplx_free_site_conc(irxn)

    select case(reaction%eqsrfcplx_rxn_surf_type(irxn))
      case(MINERAL_SURFACE)
        site_density(1) = reaction%eqsrfcplx_rxn_site_density(irxn)
!        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)* &
!                       rt_auxvar%mnrl_volfrac(reaction%eqsrfcplx_rxn_to_surf(irxn))
        num_types_of_sites = 1
      case(COLLOID_SURFACE)
        mobile_fraction = reaction%colloid_mobile_fraction(reaction%eqsrfcplx_rxn_to_surf(irxn))
        site_density(1) = (1.d0-mobile_fraction)*reaction%eqsrfcplx_rxn_site_density(irxn)
        site_density(2) = mobile_fraction*reaction%eqsrfcplx_rxn_site_density(irxn)
!        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)* &
!                       rt_auxvar%colloid%total_colloid_conc(reaction%eqsrfcplx_rxn_to_surf(irxn))
        num_types_of_sites = 2 ! two types of sites (mobile and immobile) with separate
                               ! site densities
      case(NULL_SURFACE)
        site_density(1) = reaction%eqsrfcplx_rxn_site_density(irxn)
        num_types_of_sites = 1
    end select
    
    do isite=1, num_types_of_sites
      ! isite == 1 - immobile (colloids, minerals, etc.)
      ! isite == 2 - mobile (colloids)
    
      if (site_density(isite) < 1.d-40) cycle
    
      ! get a pointer to the first complex (there will always be at least 1)
      ! in order to grab free site conc
      one_more = PETSC_FALSE
      do

        total = free_site_conc
        ln_free_site = log(free_site_conc)
        do j = 1, ncplx
          icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
          ! compute secondary species concentration
          lnQK = -reaction%eqsrfcplx_logK(icplx)*LOG_TO_LN

          ! activity of water
          if (reaction%eqsrfcplxh2oid(icplx) > 0) then
            lnQK = lnQK + reaction%eqsrfcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
          endif

          lnQK = lnQK + reaction%eqsrfcplx_free_site_stoich(icplx)* &
                        ln_free_site
        
          ncomp = reaction%eqsrfcplxspecid(0,icplx)
          do i = 1, ncomp
            icomp = reaction%eqsrfcplxspecid(i,icplx)
            lnQK = lnQK + reaction%eqsrfcplxstoich(i,icplx)*ln_act(icomp)
          enddo
          srfcplx_conc(icplx) = exp(lnQK)
          total = total + reaction%eqsrfcplx_free_site_stoich(icplx)*srfcplx_conc(icplx) 
          
        enddo
        
        if (one_more) exit
        
        if (reaction%eqsrfcplx_rxn_stoich_flag(irxn)) then 
          ! stoichiometry for free sites in one of reactions is not 1, thus must
          ! use nonlinear iteration to solve
          res = site_density(isite)-total
          
          dres_dfree_site = 1.d0

          do j = 1, ncplx
            icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
            dres_dfree_site = dres_dfree_site + &
              reaction%eqsrfcplx_free_site_stoich(icplx)* &
              srfcplx_conc(icplx)/free_site_conc
          enddo

          dfree_site_conc = res / dres_dfree_site
          free_site_conc = free_site_conc + dfree_site_conc
        
          if (dabs(dfree_site_conc/free_site_conc) < tol) then
            one_more = PETSC_TRUE
          endif
        
        else
        
          total = total / free_site_conc
          free_site_conc = site_density(isite) / total  
          
          one_more = PETSC_TRUE 
        
        endif

      enddo ! generic do
      
      rt_auxvar%eqsrfcplx_free_site_conc(irxn) = free_site_conc
   
  !!!!!!!!!!!!
      ! 2.3-46

      ! Sx = free site
      ! mi = molality of component i
      dSx_dmi = 0.d0
      tempreal = 0.d0
      do j = 1, ncplx
        icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
        ncomp = reaction%eqsrfcplxspecid(0,icplx)
        do i = 1, ncomp
          icomp = reaction%eqsrfcplxspecid(i,icplx)
          ! sum of nu_li * nu_i * S_i
          dSx_dmi(icomp) = dSx_dmi(icomp) + reaction%eqsrfcplxstoich(i,icplx)* &
                                            reaction%eqsrfcplx_free_site_stoich(icplx)* &
                                            srfcplx_conc(icplx)
        enddo
        ! sum of nu_i^2 * S_i
        tempreal = tempreal + reaction%eqsrfcplx_free_site_stoich(icplx)* & 
                              reaction%eqsrfcplx_free_site_stoich(icplx)* &
                              srfcplx_conc(icplx)
      enddo 
      ! divide denominator by Sx
      tempreal = tempreal / free_site_conc
      ! add 1.d0 to denominator
      tempreal = tempreal + 1.d0
      ! divide numerator by denominator
      dSx_dmi = -dSx_dmi / tempreal
      ! convert from dlogm to dm
      dSx_dmi = dSx_dmi / rt_auxvar%pri_molal
  !!!!!!!!!!!!
   
      do k = 1, ncplx
        icplx = reaction%eqsrfcplx_rxn_to_complex(k,irxn)

        rt_auxvar%eqsrfcplx_conc(icplx) = srfcplx_conc(icplx)

        ncomp = reaction%eqsrfcplxspecid(0,icplx)
        if (isite == 1) then ! immobile sites  
          do i = 1, ncomp
            icomp = reaction%eqsrfcplxspecid(i,icplx)
            rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + &
              reaction%eqsrfcplxstoich(i,icplx)*srfcplx_conc(icplx)
          enddo
        else ! mobile sites
          do i = 1, ncomp
            icomp = reaction%pri_spec_to_coll_spec(reaction%eqsrfcplxspecid(i,icplx))
            rt_auxvar%colloid%total_eq_mob(icomp) = rt_auxvar%colloid%total_eq_mob(icomp) + &
              reaction%eqsrfcplxstoich(i,icplx)*srfcplx_conc(icplx)
          enddo
        endif
        
        ! for 2.3-47 which feeds into 2.3-50
        nui_Si_over_Sx = reaction%eqsrfcplx_free_site_stoich(icplx)* &
                         srfcplx_conc(icplx)/ &
                         free_site_conc

        do j = 1, ncomp
          jcomp = reaction%eqsrfcplxspecid(j,icplx)
          tempreal = reaction%eqsrfcplxstoich(j,icplx)*srfcplx_conc(icplx) / &
                     rt_auxvar%pri_molal(jcomp)+ &
                     nui_Si_over_Sx*dSx_dmi(jcomp)
          if (isite == 1) then ! immobile sites                  
            do i = 1, ncomp
              icomp = reaction%eqsrfcplxspecid(i,icplx)
              rt_auxvar%dtotal_sorb_eq(icomp,jcomp) = &
                rt_auxvar%dtotal_sorb_eq(icomp,jcomp) + &
                                                   reaction%eqsrfcplxstoich(i,icplx)* &
                                                   tempreal
            enddo ! i
          else ! mobile sites
            do i = 1, ncomp
              icomp = reaction%eqsrfcplxspecid(i,icplx)
              rt_auxvar%colloid%dRj_dCj%dtotal(icomp,jcomp,1) = &
                rt_auxvar%colloid%dRj_dCj%dtotal(icomp,jcomp,1) + &
                                       reaction%eqsrfcplxstoich(i,icplx)* &
                                       tempreal
            enddo ! i
          endif
        enddo ! j
      enddo ! k
    enddo ! isite
  enddo ! irxn
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
end subroutine RTotalSorbEqSurfCplx

! ************************************************************************** !
!
! RTotalSorbEqIonx: Computes the total sorbed component concentrations and 
!                   derivative with respect to free-ion for equilibrium ion
!                   exchange
! author: Glenn Hammond
! date: 10/22/08; 05/26/09
!
! ************************************************************************** !
subroutine RTotalSorbEqIonx(rt_auxvar,global_auxvar,reaction,option)

  use Option_module
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, iphase, ncomp, ncplx
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscReal, parameter :: tol = 1.d-12
  PetscBool :: one_more
  PetscReal :: res
    
  PetscReal :: omega
  PetscReal :: ref_cation_X, ref_cation_conc, ref_cation_Z, ref_cation_k, &
               ref_cation_quotient
  PetscReal :: cation_X(reaction%naqcomp)
  PetscReal :: dres_dref_cation_X, dref_cation_X
  PetscReal :: sumZX, sumkm
    
  PetscReal :: total_pert, ref_cation_X_pert, pert
  PetscReal :: ref_cation_quotient_pert, dres_dref_cation_X_pert

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
    
  ! Ion Exchange
  if (associated(rt_auxvar%eqionx_conc)) rt_auxvar%eqionx_conc = 0.d0
  do irxn = 1, reaction%neqionxrxn

    ncomp = reaction%eqionx_rxn_cationid(0,irxn)

    ! for now we assume that omega is equal to CEC.
    omega = reaction%eqionx_rxn_CEC(irxn)

    if (reaction%eqionx_rxn_Z_flag(irxn)) then ! Zi /= Zj for any i,j

      icomp = reaction%eqionx_rxn_cationid(1,irxn)
      ref_cation_conc = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)
      ref_cation_Z = reaction%primary_spec_Z(icomp)
      ref_cation_k = reaction%eqionx_rxn_k(1,irxn)
      ref_cation_X = ref_cation_Z*rt_auxvar%eqionx_ref_cation_sorbed_conc(irxn)/omega

      one_more = PETSC_FALSE
      cation_X = 0.d0
      do

        if (ref_cation_X <= 0.d0) ref_cation_X = 0.99d0
        cation_X(1) = ref_cation_X
        ref_cation_quotient = ref_cation_X*ref_cation_k/ref_cation_conc
        total = ref_cation_X

        do j = 2, ncomp
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          cation_X(j) = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)/ &
                        reaction%eqionx_rxn_k(j,irxn)* &
                        ref_cation_quotient** &
                        (reaction%primary_spec_Z(icomp)/ref_cation_Z)
          total = total + cation_X(j)
        enddo
        
        if (one_more) exit
        
        res = 1.d0-total
          
        dres_dref_cation_X = 1.d0

#if 0
  ! test derivative
        pert = 1.d-6 * ref_cation_X
        ref_cation_X_pert = ref_cation_X + pert
        ref_cation_quotient_pert = ref_cation_X_pert*ref_cation_k/ref_cation_conc
        total_pert = ref_cation_X_pert

          do j = 2, ncomp
            icomp = reaction%eqionx_rxn_cationid(j,irxn)
            total_pert = total_pert + &
                         rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)/ &
                         reaction%eqionx_rxn_k(j,irxn)* &
                         ref_cation_quotient_pert** &
                         (reaction%primary_spec_Z(icomp)/ref_cation_Z)
          enddo
        dres_dref_cation_X_pert = (1.d0-total_pert-res)/pert
  ! test
#endif

        do j = 2, ncomp
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          dres_dref_cation_X = dres_dref_cation_X + &
            (reaction%primary_spec_Z(icomp)/ref_cation_Z)* &
            cation_X(j)/ref_cation_X
        enddo

        dref_cation_X = res / (-dres_dref_cation_X)
!        dref_cation_X = res / dres_dref_cation_X_pert
        ref_cation_X = ref_cation_X - dref_cation_X
      
        if (dabs(dref_cation_X/ref_cation_X) < tol) then
          one_more = PETSC_TRUE
        endif
      
      enddo

      rt_auxvar%eqionx_ref_cation_sorbed_conc(irxn) = ref_cation_X*omega/ref_cation_Z

    else ! Zi == Zj for all i,j
        
      sumkm = 0.d0
      cation_X = 0.d0
      
      do j = 1, ncomp  
        icomp = reaction%eqionx_rxn_cationid(j,irxn)
        cation_X(j) = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)* &
                      reaction%eqionx_rxn_k(j,irxn)
        sumkm = sumkm + cation_X(j)
      enddo
          
      cation_X = cation_X / sumkm

    endif
                
    ! sum up charges
    sumZX = 0.d0
    do i = 1, ncomp
      icomp = reaction%eqionx_rxn_cationid(i,irxn)
      sumZX = sumZX + reaction%primary_spec_Z(icomp)*cation_X(i)
    enddo

    ! compute totals based on sorbed ions
    do i = 1, ncomp
      icomp = reaction%eqionx_rxn_cationid(i,irxn)
      tempreal1 = cation_X(i)*omega/reaction%primary_spec_Z(icomp)
      ! residual function entry
      
      rt_auxvar%eqionx_conc(i,irxn) = rt_auxvar%eqionx_conc(i,irxn) + tempreal1

      rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + tempreal1

      tempreal2 = reaction%primary_spec_Z(icomp)/sumZX
      do j = 1, ncomp
        jcomp = reaction%eqionx_rxn_cationid(j,irxn)
        if (i == j) then
          rt_auxvar%dtotal_sorb_eq(icomp,jcomp) = rt_auxvar%dtotal_sorb_eq(icomp,jcomp) + &
                                               tempreal1*(1.d0-(tempreal2*cation_X(j)))/ &
                                               rt_auxvar%pri_molal(jcomp)
        else
          rt_auxvar%dtotal_sorb_eq(icomp,jcomp) = rt_auxvar%dtotal_sorb_eq(icomp,jcomp) + &
                                               (-tempreal1)*tempreal2*cation_X(j)/ &
                                               rt_auxvar%pri_molal(jcomp)
        endif
      enddo
    enddo    

  enddo
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
end subroutine RTotalSorbEqIonx

! ************************************************************************** !
!
! RMultiRateSorption: Computes contribution to the accumualtion term due
!                     due to multirate sorption
! author: Glenn Hammond
! date: 05/20/09
!
! ************************************************************************** !
subroutine RMultiRateSorption(Res,Jac,compute_derivative,rt_auxvar, &
                              global_auxvar,volume,reaction,option)

  use Option_module

  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: volume
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(option_type) :: option
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, ncomp, ncplx
  PetscInt, parameter :: iphase = 1
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: srfcplx_conc(reaction%neqsrfcplx)
  PetscReal :: dSx_dmi(reaction%naqcomp)
  PetscReal :: nui_Si_over_Sx
  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscReal, parameter :: tol = 1.d-12
  PetscBool :: one_more
  PetscReal :: residual, dres_dfree_site, dfree_site_conc
  PetscReal :: site_density
  
  PetscInt :: irate
  PetscReal :: kdt, one_plus_kdt, k_over_one_plus_kdt
  PetscReal :: total_sorb_eq(reaction%naqcomp)
  PetscReal :: dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp)

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
    
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqsrfcplx_logKcoef,reaction%eqsrfcplx_logK, &
                               global_auxvar%temp(iphase),reaction%neqsrfcplx)
  endif
#endif  

  rt_auxvar%total_sorb_eq = 0.d0
  rt_auxvar%eqsrfcplx_conc = 0.d0

  ! Surface Complexation
  do irxn = 1, reaction%neqsrfcplxrxn
  
    !WARNING! the below assumes site density multiplicative factor
    select case(reaction%eqsrfcplx_rxn_surf_type(irxn))
      case(MINERAL_SURFACE)
        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)
!        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)* &
!                       rt_auxvar%mnrl_volfrac(reaction%eqsrfcplx_rxn_to_surf(irxn))
      case(COLLOID_SURFACE)
        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)
!        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)* &
!                       rt_auxvar%colloid%total(reaction%eqsrfcplx_rxn_to_surf(irxn))
      case(NULL_SURFACE)
        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)
    end select  

    ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
    free_site_conc = rt_auxvar%eqsrfcplx_free_site_conc(irxn)

    ! get a pointer to the first complex (there will always be at least 1)
    ! in order to grab free site conc
    one_more = PETSC_FALSE
    do
      total = free_site_conc
      ln_free_site = log(free_site_conc)
      do j = 1, ncplx
        icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
        ! compute secondary species concentration
        lnQK = -reaction%eqsrfcplx_logK(icplx)*LOG_TO_LN

        ! activity of water
        if (reaction%eqsrfcplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqsrfcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
        endif

        lnQK = lnQK + reaction%eqsrfcplx_free_site_stoich(icplx)* &
                      ln_free_site
      
        ncomp = reaction%eqsrfcplxspecid(0,icplx)
        do i = 1, ncomp
          icomp = reaction%eqsrfcplxspecid(i,icplx)
          lnQK = lnQK + reaction%eqsrfcplxstoich(i,icplx)*ln_act(icomp)
        enddo
        srfcplx_conc(icplx) = exp(lnQK)
        total = total + reaction%eqsrfcplx_free_site_stoich(icplx)*srfcplx_conc(icplx) 
        
      enddo
      
      if (one_more) exit
      
      if (reaction%eqsrfcplx_rxn_stoich_flag(irxn)) then 
        ! stoichiometry for free sites in one of reactions is not 1, thus must
        ! use nonlinear iteration to solve
        residual = site_density-total
        
        dres_dfree_site = 1.d0

        do j = 1, ncplx
          icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
          dres_dfree_site = dres_dfree_site + &
            reaction%eqsrfcplx_free_site_stoich(icplx)* &
            srfcplx_conc(icplx)/free_site_conc
        enddo

        dfree_site_conc = residual / dres_dfree_site
        free_site_conc = free_site_conc + dfree_site_conc
      
        if (dabs(dfree_site_conc/free_site_conc) < tol) then
          one_more = PETSC_TRUE
        endif
      
      else
      
        total = total / free_site_conc
        free_site_conc = site_density / total  
        
        one_more = PETSC_TRUE
      endif
    enddo

    rt_auxvar%eqsrfcplx_free_site_conc(irxn) = free_site_conc
   
    dSx_dmi = 0.d0
    tempreal = 0.d0
    do j = 1, ncplx
      icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
      ncomp = reaction%eqsrfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsrfcplxspecid(i,icplx)
        ! numerator of 4.39
        dSx_dmi(icomp) = dSx_dmi(icomp) + reaction%eqsrfcplxstoich(i,icplx)* &
          reaction%eqsrfcplx_free_site_stoich(icplx)*srfcplx_conc(icplx)
      enddo
      ! denominator of 4.39
      tempreal = tempreal + reaction%eqsrfcplx_free_site_stoich(icplx)* & 
        reaction%eqsrfcplx_free_site_stoich(icplx)*srfcplx_conc(icplx)
    enddo 
    ! divide denominator by Sx
    tempreal = tempreal / free_site_conc
    ! add 1.d0 to denominator
    tempreal = tempreal + 1.d0
    ! divide numerator by denominator
    dSx_dmi = -dSx_dmi / tempreal
    ! convert from dlogm to dm
    dSx_dmi = dSx_dmi / rt_auxvar%pri_molal

    ! initialize total sorbed concentrations and derivatives
    total_sorb_eq = 0.d0
    dtotal_sorb_eq = 0.d0
      
    do k = 1, ncplx
      icplx = reaction%eqsrfcplx_rxn_to_complex(k,irxn)

      rt_auxvar%eqsrfcplx_conc(k) = &
        rt_auxvar%eqsrfcplx_conc(k) + srfcplx_conc(icplx)

      ncomp = reaction%eqsrfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsrfcplxspecid(i,icplx)
        total_sorb_eq(icomp) = total_sorb_eq(icomp) + &
          reaction%eqsrfcplxstoich(i,icplx)*srfcplx_conc(icplx)
      enddo
      
      if (compute_derivative) then
        nui_Si_over_Sx = reaction%eqsrfcplx_free_site_stoich(icplx)* &
                         srfcplx_conc(icplx)/free_site_conc

        do j = 1, ncomp
          jcomp = reaction%eqsrfcplxspecid(j,icplx)
          tempreal = reaction%eqsrfcplxstoich(j,icplx)*srfcplx_conc(icplx) / &
            rt_auxvar%pri_molal(jcomp) + nui_Si_over_Sx*dSx_dmi(jcomp)
                      
          do i = 1, ncomp
            icomp = reaction%eqsrfcplxspecid(i,icplx)
            dtotal_sorb_eq(icomp,jcomp) = dtotal_sorb_eq(icomp,jcomp) + &
              reaction%eqsrfcplxstoich(i,icplx)*tempreal
          enddo
        enddo
      endif
    enddo
  enddo
      
  ! WARNING: this assumes site fraction multiplicative factor 
  do irate = 1, reaction%kinmr_nrate
    kdt = reaction%kinmr_rate(irate) * option%tran_dt
    one_plus_kdt = 1.d0 + kdt
    k_over_one_plus_kdt = reaction%kinmr_rate(irate)/one_plus_kdt
        
    Res(:) = Res(:) + volume * k_over_one_plus_kdt * &
      (reaction%kinmr_frac(irate)*total_sorb_eq(:) - rt_auxvar%kinmr_total_sorb(:,irate))
      
    if (compute_derivative) then
      Jac = Jac + volume * k_over_one_plus_kdt * reaction%kinmr_frac(irate) * dtotal_sorb_eq
    endif

  enddo
  
  ! store the target equilibrium concentration to update the sorbed 
  ! concentration at the end of the time step.
  rt_auxvar%total_sorb_eq = total_sorb_eq
  
end subroutine RMultiRateSorption

! ************************************************************************** !
!
! RKineticSurfCplx: Computes contribution to residual and jacobian for 
!                   kinetic surface complexation reactions
! author: Glenn Hammond
! date: 12/07/09
!
! ************************************************************************** !
subroutine RKineticSurfCplx(Res,Jac,compute_derivative,rt_auxvar, &
                            global_auxvar,volume,reaction,option)

  use Option_module
  
  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: volume
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(option_type) :: option

  PetscInt :: i, j, k, l, icplx, icomp, jcomp, lcomp, ncomp, ncplx
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscInt :: irxn, isite
  PetscReal :: dt

  PetscReal :: numerator_sum(reaction%nkinsrfcplxrxn)
  PetscReal :: denominator_sum(reaction%nkinsrfcplxrxn)

  PetscReal :: denominator
  PetscReal :: fac
  PetscReal :: fac_sum(reaction%naqcomp)
  PetscReal :: lnQ(reaction%nkinsrfcplx)
  PetscReal :: Q(reaction%nkinsrfcplx)
  PetscReal :: srfcplx_conc_k(reaction%nkinsrfcplx)
  PetscReal :: srfcplx_conc_kp1(reaction%nkinsrfcplx)
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

! Members of the rt aux var object: mol/m^3
! PetscReal, pointer :: kinsrfcplx_conc(:)          ! S_{i\alpha}^k
! PetscReal, pointer :: kinsrfcplx_conc_kp1(:)      ! S_{i\alpha}^k+1
! PetscReal, pointer :: kinsrfcplx_free_site_conc(:) ! S_\alpha
  
! units
! k_f: dm^3/mol/sec
! k_b: 1/sec
! Res: mol/sec

  dt = option%tran_dt
  
! compute ion activity product and store: units mol/L
  lnQ = 0.d0
  do irxn = 1, reaction%nkinsrfcplxrxn
    ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
      if (reaction%kinsrfcplxh2oid(icplx) > 0) then
        lnQ(icplx) = lnQ(icplx) + reaction%kinsrfcplxh2ostoich(icplx)* &
          rt_auxvar%ln_act_h2o
      endif
    
      ncomp = reaction%kinsrfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%kinsrfcplxspecid(i,icplx)
        lnQ(icplx) = lnQ(icplx) + reaction%kinsrfcplxstoich(i,icplx)* &
          ln_act(icomp)
      enddo
      Q(icplx) = exp(lnQ(icplx))
    enddo
  enddo
    
  ! compute summation in numerator of 5.1-29: units mol/m^3
  numerator_sum = 0.d0
  do irxn = 1, reaction%nkinsrfcplxrxn
    isite = reaction%kinsrfcplx_rxn_to_site(irxn)
    ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
      numerator_sum(isite) = numerator_sum(isite) + &
                  rt_auxvar%kinsrfcplx_conc(icplx)/ &
                  (1.d0+reaction%kinsrfcplx_backward_rate(icplx)*dt)
    enddo
  enddo

  do irxn = 1, reaction%nkinsrfcplxrxn
    isite = reaction%kinsrfcplx_rxn_to_site(irxn)
    numerator_sum(isite) = reaction%kinsrfcplx_rxn_site_density(isite) - &
                           numerator_sum(isite)
  enddo
  
  ! compute summation in denominator of 5.1-29
  denominator_sum = 1.d0
  do irxn = 1, reaction%nkinsrfcplxrxn
    isite = reaction%kinsrfcplx_rxn_to_site(irxn)
    ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
      denominator_sum(isite) = denominator_sum(isite) + &
                               (reaction%kinsrfcplx_forward_rate(icplx)*dt)/ &
                               (1.d0+reaction%kinsrfcplx_backward_rate(icplx)*dt)* &
                               Q(icplx)
    enddo
  enddo

! compute surface complex conc. at new time step (5.1-30)  
  do irxn = 1, reaction%nkinsrfcplxrxn
    isite = reaction%kinsrfcplx_rxn_to_site(irxn)
    ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
      srfcplx_conc_k(icplx) = rt_auxvar%kinsrfcplx_conc(icplx)
      denominator = 1.d0 + reaction%kinsrfcplx_backward_rate(icplx)*dt
      srfcplx_conc_kp1(icplx) = (srfcplx_conc_k(icplx) + &
                                reaction%kinsrfcplx_forward_rate(icplx)*dt * &
                                numerator_sum(isite)/denominator_sum(isite)* &
                                Q(icplx))/denominator
      rt_auxvar%kinsrfcplx_conc_kp1(icplx) = srfcplx_conc_kp1(icplx)
    enddo
    rt_auxvar%kinsrfcplx_free_site_conc(isite) = numerator_sum(isite)/ &
                                               denominator_sum(isite)
  enddo

! compute residual (5.1-34)
  
  do irxn = 1, reaction%nkinsrfcplxrxn
    ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
      ncomp = reaction%kinsrfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%kinsrfcplxspecid(i,icplx)
        Res(icomp) = Res(icomp) + reaction%kinsrfcplxstoich(i,icplx)* &
                     (srfcplx_conc_kp1(icplx)-rt_auxvar%kinsrfcplx_conc(icplx))/ &
                     dt * volume
      enddo
    enddo
  enddo

  if (compute_derivative) then
  ! compute jacobian (5.1-39)
    fac_sum = 0.d0
    do irxn = 1, reaction%nkinsrfcplxrxn
      ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
      do k = 1, ncplx ! ncplx in rxn
        icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
        denominator = 1.d0 + reaction%kinsrfcplx_backward_rate(icplx)*dt
        fac = reaction%kinsrfcplx_forward_rate(icplx)/denominator
        ncomp = reaction%kinsrfcplxspecid(0,icplx)
        do j = 1, ncomp
          jcomp = reaction%kinsrfcplxspecid(j,icplx)
          fac_sum(jcomp) = fac_sum(jcomp) + reaction%kinsrfcplxstoich(j,icplx)* &
            fac * Q(icplx)
        enddo
      enddo
    enddo

    do irxn = 1, reaction%nkinsrfcplxrxn
      isite = reaction%kinsrfcplx_rxn_to_site(irxn)
      ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
      do k = 1, ncplx ! ncplx in rxn
        icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
        denominator = 1.d0 + reaction%kinsrfcplx_backward_rate(icplx)*dt
        fac = reaction%kinsrfcplx_forward_rate(icplx)/denominator
        ncomp = reaction%kinsrfcplxspecid(0,icplx)
        do j = 1, ncomp
          jcomp = reaction%kinsrfcplxspecid(j,icplx)
          do l = 1, ncomp
            lcomp = reaction%kinsrfcplxspecid(l,icplx)
            Jac(jcomp,lcomp) = Jac(jcomp,lcomp) + &
              (reaction%kinsrfcplxstoich(j,icplx) * fac * numerator_sum(isite) * &
              Q(icplx) * (reaction%kinsrfcplxstoich(l,icplx) - &
              dt * fac_sum(lcomp)/denominator_sum(isite)))/denominator_sum(isite) * &
              exp(-ln_conc(lcomp)) * volume
          enddo
        enddo
      enddo
    enddo
  endif
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
end subroutine RKineticSurfCplx

! ************************************************************************** !
!
! RKineticMineral: Computes the kinetic mineral precipitation/dissolution
!                  rates
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RKineticMineral(Res,Jac,compute_derivative,rt_auxvar, &
                           global_auxvar,volume,reaction,option)

  use Option_module
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscInt :: i, j, k, imnrl, icomp, jcomp, kcplx, iphase, ncomp
  PetscInt :: ipref, ipref_species
  ! I am assuming a maximum of 10 prefactors and 5 species per prefactor
  PetscReal :: tempreal, tempreal2
  PetscReal :: affinity_factor, sign_
  PetscReal :: Im, Im_const, dIm_dQK
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_sec(reaction%neqcplx) 
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: ln_sec_act(reaction%neqcplx)
  PetscReal :: QK, lnQK, dQK_dCj, dQK_dmj, den

  PetscReal :: ln_spec_act, spec_act_coef, ln_spec_conc
  PetscReal :: ln_prefactor, ln_numerator, ln_denominator
  PetscReal :: prefactor(10), ln_prefactor_spec(5,10)
  PetscReal :: sum_prefactor_rate
  PetscReal :: dIm_dsum_prefactor_rate, dIm_dspec
  PetscReal :: dprefactor_dprefactor_spec, dprefactor_spec_dspec
  PetscReal :: dprefactor_spec_dspec_numerator
  PetscReal :: dprefactor_spec_dspec_denominator
  PetscReal :: denominator
  PetscInt ::  icplx
  PetscReal :: ln_gam_m_beta

  PetscInt, parameter :: needs_to_be_fixed = 1
  
  PetscReal :: arrhenius_factor, rgas = 8.3144621d-3

  iphase = 1                         

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  if (reaction%neqcplx > 0) then
    ln_sec = log(rt_auxvar%sec_molal)
    ln_sec_act = ln_sec+log(rt_auxvar%sec_act_coef)
  endif

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%kinmnrl_logKcoef,reaction%kinmnrl_logK, &
                                 global_auxvar%temp(iphase),reaction%nkinmnrl)
  endif
#endif  

  do imnrl = 1, reaction%nkinmnrl ! for each mineral
    ! compute ion activity product
    lnQK = -reaction%kinmnrl_logK(imnrl)*LOG_TO_LN

    ! activity of water
    if (reaction%kinmnrlh2oid(imnrl) > 0) then
      lnQK = lnQK + reaction%kinmnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
    endif

    ncomp = reaction%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = reaction%kinmnrlspecid(i,imnrl)
      lnQK = lnQK + reaction%kinmnrlstoich(i,imnrl)*ln_act(icomp)
    enddo
    
    if (lnQK <= 6.90776d0) then
      QK = exp(lnQK)
    else
      QK = 1.d3
    endif
    
    if (associated(reaction%kinmnrl_Tempkin_const)) then
      affinity_factor = 1.d0-QK**(1.d0/reaction%kinmnrl_Tempkin_const(imnrl))
    else
      affinity_factor = 1.d0-QK
    endif
    
    sign_ = sign(1.d0,affinity_factor)

    if (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then

!   if ((reaction%kinmnrl_irreversible(imnrl) == 0 &
!     .and. (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0)) &
!     .or. (reaction%kinmnrl_irreversible(imnrl) == 1 .and. sign_ < 0.d0)) then
    
!     check for supersaturation threshold for precipitation
!     if (associated(reaction%kinmnrl_affinity_threshold)) then
      if (reaction%kinmnrl_affinity_threshold(imnrl) > 0.d0) then
        if (sign_ < 0.d0 .and. QK < reaction%kinmnrl_affinity_threshold(imnrl)) cycle
      endif
    
!     check for rate limiter for precipitation
      if (reaction%kinmnrl_rate_limiter(imnrl) > 0.d0) then
        affinity_factor = affinity_factor/(1.d0+(1.d0-affinity_factor) &
          /reaction%kinmnrl_rate_limiter(imnrl))
      endif

      ! compute prefactor
      if (reaction%kinmnrl_num_prefactors(imnrl) > 0) then
        sum_prefactor_rate = 0
        prefactor = 0.d0
        ln_prefactor_spec = 0.d0
        ! sum over parallel prefactors
        do ipref = 1, reaction%kinmnrl_num_prefactors(imnrl)
          ln_prefactor = 0.d0
          ! product of "monod" equations
          do ipref_species = 1, reaction%kinmnrl_prefactor_id(0,ipref,imnrl)
            icomp = reaction%kinmnrl_prefactor_id(ipref_species,ipref,imnrl)
            if (icomp > 0) then ! primary species
              ln_spec_act = ln_act(icomp)
            else ! secondary species (given a negative id to differentiate)
              ln_spec_act = ln_sec_act(-icomp)
            endif
            ln_numerator = &
              reaction%kinmnrl_pref_alpha(ipref_species,ipref,imnrl)* &
              ln_spec_act
            ln_denominator = log(1.d0 + &
              exp(log(reaction%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl)) + &
                  reaction%kinmnrl_pref_beta(ipref_species,ipref,imnrl)* &
                  ln_spec_act))
            ln_prefactor = ln_prefactor + ln_numerator
            ln_prefactor = ln_prefactor - ln_denominator
            ln_prefactor_spec(ipref_species,ipref) = ln_numerator - ln_denominator
          enddo
          prefactor(ipref) = exp(ln_prefactor)
        ! Arrhenius factor
          arrhenius_factor = 1.d0
          if (reaction%kinmnrl_pref_activation_energy(ipref,imnrl) > 0.d0) then
            arrhenius_factor = &
              exp(reaction%kinmnrl_pref_activation_energy(ipref,imnrl)/rgas &
                  *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+ &
                                                273.15d0)))
          endif
          sum_prefactor_rate = sum_prefactor_rate + prefactor(ipref)* &
                               reaction%kinmnrl_pref_rate(ipref,imnrl)* &
                               arrhenius_factor
        enddo
      else
        ! Arrhenius factor
        arrhenius_factor = 1.d0
        if (reaction%kinmnrl_activation_energy(imnrl) > 0.d0) then
          arrhenius_factor = exp(reaction%kinmnrl_activation_energy(imnrl)/rgas &
            *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+273.15d0)))
        endif
        sum_prefactor_rate = reaction%kinmnrl_rate(imnrl)*arrhenius_factor
      endif

      ! compute rate
      ! rate = mol/cm^2 mnrl/sec
      ! area = cm^2 mnrl/cm^3 bulk
      ! volume = m^3 bulk
      ! units = cm^2 mnrl/m^3 bulk
      
      Im_const = -rt_auxvar%mnrl_area(imnrl)*1.d6 ! convert cm^3->m^3
      ! units = mol/sec/m^3 bulk

      if (associated(reaction%kinmnrl_affinity_power)) then
        Im = Im_const*sign_*abs(affinity_factor)**reaction%kinmnrl_affinity_power(imnrl)*sum_prefactor_rate
      else
        Im = Im_const*sign_*abs(affinity_factor)*sum_prefactor_rate
      endif
      rt_auxvar%mnrl_rate(imnrl) = Im ! mol/sec/m^3
    else
      rt_auxvar%mnrl_rate(imnrl) = 0.d0
      cycle
    endif
    
    ! units = cm^2 mnrl
    Im_const = Im_const*volume
    ! units = mol/sec
    Im = Im*volume
    
!   print *,'RKineticMineral: ',imnrl,Im,QK

    ncomp = reaction%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = reaction%kinmnrlspecid(i,imnrl)
      Res(icomp) = Res(icomp) + reaction%kinmnrlstoich(i,imnrl)*Im
    enddo 
    
    if (.not. compute_derivative) cycle   

    ! calculate derivatives of rate with respect to free
    ! units = mol/sec
    if (associated(reaction%kinmnrl_affinity_power)) then
      dIm_dQK = -Im*reaction%kinmnrl_affinity_power(imnrl)/abs(affinity_factor)
    else
      dIm_dQK = -Im_const*sum_prefactor_rate
    endif
    
    if (associated(reaction%kinmnrl_Tempkin_const)) then
      dIm_dQK = dIm_dQK*(1.d0/reaction%kinmnrl_Tempkin_const(imnrl))/QK
    endif
    
    ! derivatives with respect to primary species in reaction quotient
    if (reaction%kinmnrl_rate_limiter(imnrl) <= 0.d0) then
      do j = 1, ncomp
        jcomp = reaction%kinmnrlspecid(j,imnrl)
        ! unit = L water/mol
        dQK_dCj = reaction%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
        ! units = (L water/mol)*(kg water/m^3 water)*(m^3 water/1000 L water) = kg water/mol
        dQK_dmj = dQK_dCj*global_auxvar%den_kg(iphase)*1.d-3 ! the multiplication by density could be moved
                                     ! outside the loop
        do i = 1, ncomp
          icomp = reaction%kinmnrlspecid(i,imnrl)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                             reaction%kinmnrlstoich(i,imnrl)*dIm_dQK*dQK_dmj
        enddo
      enddo
      
    else

      den = 1.d0+(1.d0-affinity_factor)/reaction%kinmnrl_rate_limiter(imnrl)
      do j = 1, ncomp
        jcomp = reaction%kinmnrlspecid(j,imnrl)
        ! unit = L water/mol
        dQK_dCj = reaction%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
        ! units = (L water/mol)*(kg water/m^3 water)*(m^3 water/1000 L water) = kg water/mol
        dQK_dmj = dQK_dCj*global_auxvar%den_kg(iphase)*1.d-3 ! the multiplication by density could be moved
                                     ! outside the loop
        do i = 1, ncomp
          icomp = reaction%kinmnrlspecid(i,imnrl)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
            reaction%kinmnrlstoich(i,imnrl)*dIm_dQK  &
            *(1.d0 + QK/reaction%kinmnrl_rate_limiter(imnrl)/den)*dQK_dmj/den
        enddo
      enddo
    endif
    
    if (reaction%kinmnrl_num_prefactors(imnrl) > 0) then ! add contribution of derivative in prefactor - messy
#if 1      
      dIm_dsum_prefactor_rate = Im/sum_prefactor_rate
      ! summation over parallel reactions (prefactors)
      do ipref = 1, reaction%kinmnrl_num_prefactors(imnrl)
        arrhenius_factor = 1.d0
        if (reaction%kinmnrl_pref_activation_energy(ipref,imnrl) > 0.d0) then
          arrhenius_factor = &
            exp(reaction%kinmnrl_pref_activation_energy(ipref,imnrl)/rgas &
                *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+ &
                                              273.15d0)))
        endif
        ! prefactor() saved in residual calc above
        ln_prefactor = log(prefactor(ipref))
        ! product of "monod" equations
        do ipref_species = 1, reaction%kinmnrl_prefactor_id(0,ipref,imnrl)
          ! derivative of 54 with respect to a single "monod" equation
          ! ln_prefactor_spec(,) saved in residual calc above
          dprefactor_dprefactor_spec = ln_prefactor-ln_prefactor_spec(ipref_species,ipref)
          icomp = reaction%kinmnrl_prefactor_id(ipref_species,ipref,imnrl)
          if (icomp > 0) then ! primary species
            ln_spec_conc = ln_act(icomp)
            spec_act_coef = rt_auxvar%pri_act_coef(icomp)
          else ! secondary species
            ln_spec_conc = ln_sec_act(-icomp)
            spec_act_coef = rt_auxvar%sec_act_coef(-icomp)
          endif
          ! derivative of numerator in eq. 54 with respect to species activity
          dprefactor_spec_dspec_numerator = &
            reaction%kinmnrl_pref_alpha(ipref_species,ipref,imnrl) * &
            exp(ln_prefactor_spec(ipref_species,ipref) - ln_spec_act)
          ln_gam_m_beta = reaction%kinmnrl_pref_beta(ipref_species,ipref,imnrl)* &
                          ln_spec_act
          ! denominator
          denominator = 1.d0 + &
              exp(log(reaction%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl)) + &
                  ln_gam_m_beta)
          ! derivative of denominator in eq. 54 with respect to species activity
          dprefactor_spec_dspec_denominator = -1.d0 * &
            exp(ln_prefactor_spec(ipref_species,ipref)) / denominator * &
            reaction%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl) * &
            reaction%kinmnrl_pref_beta(ipref_species,ipref,imnrl) * &
            exp(ln_gam_m_beta - ln_spec_act)

          ! chain rule for derivative of "monod" equation
          dprefactor_spec_dspec = dprefactor_spec_dspec_numerator + &
            dprefactor_spec_dspec_denominator

          ! thus far the derivative is with respect to the activity, convert to with
          ! respect to molality
          dprefactor_spec_dspec = dprefactor_spec_dspec * spec_act_coef

          dIm_dspec = dIm_dsum_prefactor_rate * dprefactor_dprefactor_spec * &
                      dprefactor_spec_dspec * &
                      reaction%kinmnrl_pref_rate(ipref,imnrl)* &
                      arrhenius_factor

           
          if (icomp > 0) then 
            ! add derivative for primary species
            Jac(icomp,icomp) = Jac(icomp,icomp) + dIm_dspec
          else ! secondary species -- have to calculate the derivative
            ! have to recalculate the reaction quotient (QK) for secondary species
            icplx = -icomp

            ! compute secondary species concentration
            lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

            ! activity of water
            if (reaction%eqcplxh2oid(icplx) > 0) then
              lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
            endif

            ncomp = reaction%eqcplxspecid(0,icplx)
            do i = 1, ncomp
              icomp = reaction%eqcplxspecid(i,icplx)
              lnQK = lnQK + reaction%eqcplxstoich(i,icplx)*ln_act(icomp)
            enddo
            ! add contribution to derivatives secondary prefactor with respect to free
            do j = 1, ncomp
              jcomp = reaction%eqcplxspecid(j,icplx)
              tempreal = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                rt_auxvar%sec_act_coef(icplx)
              do i = 1, ncomp
                icomp = reaction%eqcplxspecid(i,icplx)
                Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                  reaction%eqcplxstoich(i,icplx)*tempreal*dIm_dspec
              enddo
            enddo
          endif
        enddo
      enddo  ! loop over prefactors
#endif
    endif
  enddo  ! loop over minerals
    
end subroutine RKineticMineral

! ************************************************************************** !
!
! RMineralSaturationIndex: Calculates the mineral saturation index
! author: Glenn Hammond
! date: 08/29/11
!
! ************************************************************************** !
function RMineralSaturationIndex(imnrl,rt_auxvar,global_auxvar,reaction,option)

  use Option_module
  
  type(option_type) :: option
  PetscInt :: imnrl
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscReal :: RMineralSaturationIndex
  PetscInt :: i, icomp
  PetscReal :: lnQK
  PetscInt, parameter :: iphase = 1

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%mnrl_logKcoef,reaction%mnrl_logK, &
                                 global_auxvar%temp(iphase),reaction%nmnrl)
  endif
#endif  

  ! compute saturation
  lnQK = -reaction%mnrl_logK(imnrl)*LOG_TO_LN
  if (reaction%mnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + reaction%mnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
  endif
  do i = 1, reaction%mnrlspecid(0,imnrl)
    icomp = reaction%mnrlspecid(i,imnrl)
    lnQK = lnQK + reaction%mnrlstoich(i,imnrl)* &
           log(rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp))
  enddo
  RMineralSaturationIndex = exp(lnQK)    

end function RMineralSaturationIndex

! ************************************************************************** !
!
! RAccumulationSorb: Computes non-aqueous portion of the accumulation term in 
!                    residual function
! author: Glenn Hammond
! date: 05/26/09
!
! ************************************************************************** !
subroutine RAccumulationSorb(rt_aux_var,global_aux_var,vol,reaction, &
                             option,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscReal :: vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  
  PetscReal :: v_t
  
  ! units = (mol solute/m^3 bulk)*(m^3 bulk)/(sec) = mol/sec
  ! all residual entries should be in mol/sec
  v_t = vol/option%tran_dt
  Res(1:reaction%naqcomp) = Res(1:reaction%naqcomp) + &
    rt_aux_var%total_sorb_eq(:)*v_t

end subroutine RAccumulationSorb


! ************************************************************************** !
!
! RAccumulationSorbDerivative: Computes derivative of non-aqueous portion of 
!                              the accumulation term in residual function 
! author: Glenn Hammond
! date: 05/26/09
!
! ************************************************************************** !
subroutine RAccumulationSorbDerivative(rt_aux_var,global_aux_var, &
                                       vol,reaction,option,J)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var  
  PetscReal :: vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp
  PetscReal :: v_t
  
  ! units = (kg water/m^3 bulk)*(m^3 bulk)/(sec) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  v_t = vol/option%tran_dt
  J(1:reaction%naqcomp,1:reaction%naqcomp) = &
    J(1:reaction%naqcomp,1:reaction%naqcomp) + &
    rt_aux_var%dtotal_sorb_eq(:,:)*v_t

end subroutine RAccumulationSorbDerivative

! ************************************************************************** !
!
! RGeneral: Computes the general reaction rates
! author: Glenn Hammond
! date: 09/08/10
!
! ************************************************************************** !
subroutine RGeneral(Res,Jac,compute_derivative,rt_auxvar,global_auxvar, &
                    porosity,volume,reaction,option)

  use Option_module
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  
  PetscInt :: i, j, icomp, jcomp, irxn, ncomp
  PetscReal :: tempreal
  PetscReal :: kf, kr
  PetscReal :: Qkf, lnQkf
  PetscReal :: Qkr, lnQkr
  PetscReal :: por_den_sat_vol

  PetscInt, parameter :: iphase = 1

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  do irxn = 1, reaction%ngeneral_rxn ! for each mineral
    
    ! units
    ! for nth-order reaction
    ! kf/kr = kg^(n-1)/mol^(n-1)-sec
    ! thus for a 1st-order reaction, kf units = 1/sec
    
    kf = reaction%general_kf(irxn)
    kr = reaction%general_kr(irxn)

    if (kf > 0.d0) then
      ! compute ion activity product
      lnQkf = log(kf)

  ! currently not accommodating activity of water
      ! activity of water
  !    if (reaction%kinmnrlh2oid(irxn) > 0) then
  !      lnQkf = lnQkf + reaction%generalh2ostoich(irxn)*rt_auxvar%ln_act_h2o
  !    endif

      ncomp = reaction%generalforwardspecid(0,irxn)
      do i = 1, ncomp
        icomp = reaction%generalforwardspecid(i,irxn)
        lnQkf = lnQkf + reaction%generalforwardstoich(i,irxn)*ln_act(icomp)
      enddo
      Qkf = exp(lnQkf)
    else
      Qkf = 0.d0
    endif
    
    if (kr > 0.d0) then
      lnQkr = log(kr)

  ! currently not accommodating activity of water
      ! activity of water
  !    if (reaction%kinmnrlh2oid(irxn) > 0) then
  !      lnQkr = lnQkr + reaction%generalh2ostoich(irxn)*rt_auxvar%ln_act_h2o
  !    endif

      ncomp = reaction%generalbackwardspecid(0,irxn)
      do i = 1, ncomp
        icomp = reaction%generalbackwardspecid(i,irxn)
        lnQkr = lnQkr + reaction%generalbackwardstoich(i,irxn)*ln_act(icomp)
      enddo
      Qkr = exp(lnQkr)
    else
      Qkr = 0.d0
    endif

    ! Qkf/Qkr units are now mol/kg(water)-sec

    por_den_sat_vol = porosity*global_auxvar%den_kg(iphase)* &
                      global_auxvar%sat(iphase)*volume

    ncomp = reaction%generalspecid(0,irxn)
    do i = 1, ncomp
      icomp = reaction%generalspecid(i,irxn)
      ! units = mol/sec
      Res(icomp) = Res(icomp) - reaction%generalstoich(i,irxn)*(Qkf-Qkr)* &
                                por_den_sat_vol
    enddo 

    if (.not. compute_derivative) cycle   

    ! calculate derivatives of rate with respect to free
    ! units = mol/sec

    if (kf > 0.d0) then
      ! derivatives with respect to primary species in forward reaction
      do j = 1, reaction%generalforwardspecid(0,irxn)
        jcomp = reaction%generalforwardspecid(j,irxn)
        tempreal = -1.d0*reaction%generalforwardstoich(j,irxn)*exp(lnQkf-ln_conc(jcomp))* &
                   por_den_sat_vol
        do i = 1, reaction%generalspecid(0,irxn)
          icomp = reaction%generalspecid(i,irxn)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                             reaction%generalstoich(i,irxn)*tempreal
        enddo
      enddo
    endif
    
    if (kr > 0.d0) then
      ! derivatives with respect to primary species in forward reaction
      do j = 1, reaction%generalbackwardspecid(0,irxn)
        jcomp = reaction%generalbackwardspecid(j,irxn)
        tempreal = reaction%generalbackwardstoich(j,irxn)*exp(lnQkr-ln_conc(jcomp))* &
                   por_den_sat_vol
        do i = 1, reaction%generalspecid(0,irxn)
          icomp = reaction%generalspecid(i,irxn)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                             reaction%generalstoich(i,irxn)*tempreal
        enddo
      enddo
    endif

  enddo  ! loop over reactions
    
end subroutine RGeneral

! ************************************************************************** !
!
! RSolve: Computes the kinetic mineral precipitation/dissolution
!                  rates
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RSolve(Res,Jac,conc,update,ncomp,use_log_formulation)

  use Utility_module
  
  implicit none

  PetscInt :: ncomp
  PetscReal :: Res(ncomp)
  PetscReal :: Jac(ncomp,ncomp)
  PetscReal :: update(ncomp)
  PetscReal :: conc(ncomp)
  PetscBool :: use_log_formulation
  
  PetscInt :: indices(ncomp)
  PetscReal :: rhs(ncomp)
  PetscInt :: icomp
  PetscReal :: norm

  ! scale Jacobian
  do icomp = 1, ncomp
    norm = max(1.d0,maxval(abs(Jac(icomp,:))))
    norm = 1.d0/norm
    rhs(icomp) = Res(icomp)*norm
    Jac(icomp,:) = Jac(icomp,:)*norm
  enddo
    
  if (use_log_formulation) then
    ! for derivatives with respect to ln conc
    do icomp = 1, ncomp
      Jac(:,icomp) = Jac(:,icomp)*conc(icomp)
    enddo
  endif

  call ludcmp(Jac,ncomp,indices,icomp)
  call lubksb(Jac,ncomp,indices,rhs)
  
  update = rhs

end subroutine RSolve

! ************************************************************************** !
!
! ReactionFitLogKCoef: Least squares fit to log K over database temperature range
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionFitLogKCoef(coefs,logK,name,option,reaction)

  use Option_module
  use Utility_module

  implicit none
  
  type(reaction_type) :: reaction
  PetscReal :: coefs(FIVE_INTEGER)
  character(len=MAXWORDLENGTH) :: name 
  PetscReal :: logK(reaction%num_dbase_temperatures)
  type(option_type) :: option

  PetscInt :: temp_int(reaction%num_dbase_temperatures), &
              indx(reaction%num_dbase_temperatures)
  PetscReal :: a(FIVE_INTEGER,FIVE_INTEGER), &
               vec(FIVE_INTEGER,reaction%num_dbase_temperatures), temperature_kelvin

  PetscInt :: i, j, k, iflag
  
  ! need to fill in vec with equations for temperatures vs coefs.
  
  do i = 1, reaction%num_dbase_temperatures
    temperature_kelvin = reaction%dbase_temperatures(i) + 273.15d0
    vec(1,i) = log(temperature_kelvin)
    vec(2,i) = 1.d0
    vec(3,i) = temperature_kelvin
    vec(4,i) = 1.d0/temperature_kelvin
    vec(5,i) = 1.d0/(temperature_kelvin*temperature_kelvin)
  enddo
  
  iflag = 0
  do j = 1, FIVE_INTEGER
    coefs(j) = 0.d0
    do i = 1, reaction%num_dbase_temperatures
      if (dabs(logK(i) - 500.) < 1.d-10) then
        iflag = 1
        temp_int(i) = ZERO_INTEGER
        option%io_buffer = 'In ReactionFitLogKCoef: log K .gt. 500 for ' // &
                           trim(name)
        call printWrnMsg(option)
      else
        coefs(j) = coefs(j) + vec(j,i)*logK(i)
        temp_int(i) = ONE_INTEGER
      endif
    enddo
  enddo

  do j = 1, FIVE_INTEGER
    do k = j, FIVE_INTEGER
      a(j,k) = 0.d0
      do i = 1, reaction%num_dbase_temperatures
        if (temp_int(i) .eq. 1) then
          a(j,k) = a(j,k) + vec(j,i)*vec(k,i)
        endif
      enddo
      if (j .ne. k) a(k,j) = a(j,k)
    enddo
  enddo

  call ludcmp(a,FIVE_INTEGER,indx,i)
  call lubksb(a,FIVE_INTEGER,indx,coefs)

end subroutine ReactionFitLogKCoef

! ************************************************************************** !
!
! ReactionInitializeLogK: Least squares fit to log K over database temperature range
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionInitializeLogK(logKcoef,logKs,logK,option,reaction)

  use Option_module

  implicit none
  
  type(reaction_type) :: reaction
  PetscReal :: logKcoef(FIVE_INTEGER)
  PetscReal :: logKs(reaction%num_dbase_temperatures)
  PetscReal :: logK, logK_1D_Array(ONE_INTEGER)
  type(option_type) :: option
  
  PetscReal :: coefs(FIVE_INTEGER,ONE_INTEGER)
  PetscReal :: temperature
  PetscInt :: itemperature
  PetscInt :: i
  
  ! we always initialize on reference temperature
  temperature = option%reference_temperature
  
  itemperature = 0
  if (option%use_isothermal) then ! find database temperature if relevant
    do i = 1, reaction%num_dbase_temperatures
      if (dabs(option%reference_temperature - &
               reaction%dbase_temperatures(i)) < 1.d-10) then
        itemperature = i
        exit
      endif
    enddo
  endif
  
  if (itemperature > 0) then ! use database temperature
    logK = logKs(itemperature)
  else                       ! interpolate
    coefs(:,ONE_INTEGER) = logKcoef(:)
    call ReactionInterpolateLogK(coefs,logK_1D_Array,temperature,ONE_INTEGER)
    logK = logK_1D_Array(ONE_INTEGER)
  endif

end subroutine ReactionInitializeLogK

! ************************************************************************** !
!
! ReactionInterpolateLogK: Interpolation log K function: temp - temperature [C]
!                             b - fit coefficients determined from fit(...)
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionInterpolateLogK(coefs,logKs,temp,n)

  PetscInt :: n
  PetscReal :: coefs(5,n), logKs(n), temp

  PetscInt :: i
  PetscReal :: temp_kelvin
  
  temp_kelvin = temp + 273.15d0
  
  do i = 1, n
    logKs(i) = coefs(1,i)*log(temp_kelvin) &
             + coefs(2,i)           &
             + coefs(3,i)*temp_kelvin      &
             + coefs(4,i)/temp_kelvin      &
             + coefs(5,i)/(temp_kelvin*temp_kelvin)
  enddo
  
end subroutine ReactionInterpolateLogK

! ************************************************************************** !
!
! RComputeKd: Computes the Kd for a given chemical component
! author: Glenn Hammond
! date: 05/14/09
!
! ************************************************************************** !
subroutine ReactionComputeKd(icomp,retardation,rt_auxvar,global_auxvar, &
                             porosity,reaction,option)

  use Option_module
  
  PetscInt :: icomp
  PetscReal :: retardation
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar  
  PetscReal :: porosity
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscReal :: bulk_vol_to_fluid_vol
  PetscInt :: i, j, jcomp, irxn, icplx, irate
  PetscInt, parameter :: iphase = 1

  retardation = 0.d0
  if (reaction%neqsorb == 0) return
  
  bulk_vol_to_fluid_vol = porosity*global_auxvar%sat(iphase)*1000.d0

  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
#if 0
    ! we should be able to use total_sorb instead of summing complexes
    do irxn = 1, reaction%neqsrfcplxrxn
      do i = 1, reaction%eqsrfcplx_rxn_to_complex(0,irxn)
        icplx = reaction%eqsrfcplx_rxn_to_complex(i,irxn)
        do j = 1, reaction%eqsrfcplxspecid(0,icplx)
          jcomp = reaction%eqsrfcplxspecid(j,icplx)
          if (icomp == jcomp) then
            retardation = retardation + &
              reaction%eqsrfcplxstoich(j,icplx) * &
              rt_auxvar%eqsrfcplx_conc(icplx)
            exit
          endif
        enddo
      enddo
    enddo
#else
    retardation = rt_auxvar%total_sorb_eq(icomp)
#endif
  else
    do irate = 1, reaction%kinmr_nrate
      retardation = retardation + rt_auxvar%kinmr_total_sorb(icomp,irate)
    enddo
  endif

  if (dabs(rt_auxvar%total(icomp,iphase)) > 1.d-40) &
    retardation = retardation/bulk_vol_to_fluid_vol/ &
      rt_auxvar%total(icomp,iphase)

end subroutine ReactionComputeKd

! ************************************************************************** !
!
! RAge: Computes the ages of the groundwater
! author: Glenn Hammond
! date: 02/22/10
!
! ************************************************************************** !
subroutine RAge(rt_aux_var,global_aux_var,por,vol,option,reaction,Res)

  use Option_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var  
  PetscReal :: por,vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscInt, parameter :: iphase = 1
  
  Res(:) = 0.d0
  if (reaction%calculate_water_age) then
    Res(reaction%species_idx%water_age_id) = por*global_aux_var%sat(iphase)* &
      1000.d0 * vol
  endif
  if (reaction%calculate_tracer_age) then
    Res(reaction%species_idx%tracer_age_id) = &
      -rt_aux_var%total(reaction%species_idx%tracer_aq_id,iphase)* &
      por*global_aux_var%sat(iphase)*1000.d0*vol
  endif
end subroutine RAge

! ************************************************************************** !
!
! RTAuxVarCompute: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine RTAuxVarCompute(rt_aux_var,global_aux_var,reaction,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var
  
#if 0  
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp
  PetscReal :: dtotal(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dtotalsorb(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: pert
  type(reactive_transport_auxvar_type) :: rt_auxvar_pert
#endif

  ! any changes to the below must also be updated in 
  ! Reaction.F90:RReactionDerivative()
  
!already set  rt_aux_var%pri_molal = x

  call RTotal(rt_aux_var,global_aux_var,reaction,option)
  if (reaction%neqsorb > 0) then
    call RTotalSorb(rt_aux_var,global_aux_var,reaction,option)
  endif

#if 0
! numerical check
  Res_orig = 0.d0
  dtotal = 0.d0
  dtotalsorb = 0.d0
  option%iflag = 0 ! be sure not to allocate mass_balance array
  call RTAuxVarInit(rt_auxvar_pert,reaction,option)
  do jcomp = 1, reaction%naqcomp
    Res_pert = 0.d0
    call RTAuxVarCopy(rt_auxvar_pert,rt_aux_var,option)
    if (reaction%neqcplx > 0) then
      rt_aux_var%sec_molal = 0.d0
    endif
    if (reaction%ngas > 0) then
      rt_aux_var%gas_molal = 0.d0
    endif
    if (reaction%neqsrfcplxrxn > 0) then
      rt_auxvar_pert%eqsrfcplx_free_site_conc = 1.d-9
      rt_auxvar_pert%eqsrfcplx_conc = 0.d0
    endif
    if (reaction%neqionxrxn > 0) then
      rt_aux_var%eqionx_ref_cation_sorbed_conc = 1.d-9
    endif
    pert = rt_auxvar_pert%pri_molal(jcomp)*perturbation_tolerance
    rt_auxvar_pert%pri_molal(jcomp) = rt_auxvar_pert%pri_molal(jcomp) + pert
    
    call RTotal(rt_auxvar_pert,global_aux_var,reaction,option)
    dtotal(:,jcomp) = (rt_auxvar_pert%total(:,1) - rt_aux_var%total(:,1))/pert
    if (reaction%neqsorb > 0) then
      call RTotalSorb(rt_auxvar_pert,global_aux_var,reaction,option)
      if (reaction%kinmr_nrate <= 0) &
        dtotalsorb(:,jcomp) = (rt_auxvar_pert%total_sorb_eq(:) - &
                               rt_aux_var%total_sorb_eq(:))/pert
    endif
  enddo
  do icomp = 1, reaction%naqcomp
    do jcomp = 1, reaction%naqcomp
      if (dabs(dtotal(icomp,jcomp)) < 1.d-16) dtotal(icomp,jcomp) = 0.d0
      if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
        if (dabs(dtotalsorb(icomp,jcomp)) < 1.d-16) dtotalsorb(icomp,jcomp) = 0.d0
      endif
    enddo
  enddo
  rt_aux_var%aqueous%dtotal(:,:,1) = dtotal
  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) &
    rt_aux_var%dtotal_sorb_eq = dtotalsorb
  call RTAuxVarDestroy(rt_auxvar_pert)
#endif
  
end subroutine RTAuxVarCompute

! ************************************************************************** !
!
! RTAccumulation: Computes aqueous portion of the accumulation term in 
!                 residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTAccumulation(rt_aux_var,global_aux_var,por,vol,reaction,option,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscReal :: por, vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  
  PetscInt :: iphase
  PetscInt :: istart, iend
  PetscInt :: idof
  PetscInt :: icoll
  PetscInt :: icollcomp
  PetscInt :: iaqcomp
  PetscReal :: psv_t
  PetscReal :: v_t
  
  iphase = 1
  ! units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  !         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  ! 1000.d0 converts vol from m^3 -> L
  ! all residual entries should be in mol/sec
  psv_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt  
  istart = 1
  iend = reaction%naqcomp
  Res(istart:iend) = psv_t*rt_aux_var%total(:,iphase) 

  if (reaction%ncoll > 0) then
    do icoll = 1, reaction%ncoll
      idof = reaction%offset_coll + icoll
      Res(idof) = psv_t*rt_aux_var%colloid%conc_mob(icoll)
    enddo
  endif
  if (reaction%ncollcomp > 0) then
    do icollcomp = 1, reaction%ncollcomp
      iaqcomp = reaction%coll_spec_to_pri_spec(icollcomp)
      Res(iaqcomp) = Res(iaqcomp) + &
        psv_t*rt_aux_var%colloid%total_eq_mob(icollcomp)
    enddo
  endif

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2
  do 
    iphase = iphase + 1
    if (iphase > option%nphase) exit

! super critical CO2 phase
    if (iphase == 2) then
      psv_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt 
      Res(istart:iend) = Res(istart:iend) + psv_t*rt_aux_var%total(:,iphase) 
      ! should sum over gas component only need more implementations
    endif 
! add code for other phases here
  enddo
#endif
  
end subroutine RTAccumulation

! ************************************************************************** !
!
! RTAccumulationDerivative: Computes derivative of aqueous portion of the 
!                           accumulation term in residual function 
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTAccumulationDerivative(rt_aux_var,global_aux_var, &
                                    por,vol,reaction,option,J)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var  
  PetscReal :: por, vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp, iphase
  PetscInt :: istart, iendaq
  PetscInt :: idof
  PetscInt :: icoll
  PetscReal :: psvd_t, v_t

  iphase = 1
  istart = 1
  iendaq = reaction%naqcomp 
  ! units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  !         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  J = 0.d0
  if (associated(rt_aux_var%aqueous%dtotal)) then ! units of dtotal = kg water/L water
    psvd_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt
    J(istart:iendaq,istart:iendaq) = rt_aux_var%aqueous%dtotal(:,:,iphase)*psvd_t
  else
    psvd_t = por*global_aux_var%sat(iphase)* &
             global_aux_var%den_kg(iphase)*vol/option%tran_dt ! units of den = kg water/m^3 water
    do icomp=istart,iendaq
      J(icomp,icomp) = psvd_t
    enddo
  endif

  if (reaction%ncoll > 0) then
    do icoll = 1, reaction%ncoll
      idof = reaction%offset_coll + icoll
      ! shouldn't have to sum a this point
      J(idof,idof) = psvd_t
    enddo
  endif
  if (reaction%ncollcomp > 0) then
    ! dRj_dCj - mobile
    J(istart:iendaq,istart:iendaq) = J(istart:iendaq,istart:iendaq) + &
      rt_aux_var%colloid%dRj_dCj%dtotal(:,:,1)*psvd_t
    ! need the below
    ! dRj_dSic
    ! dRic_dCj                                 
  endif

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2
  do
    iphase = iphase +1 
    if (iphase > option%nphase) exit
! super critical CO2 phase
    if (iphase == 2) then
      if (associated(rt_aux_var%aqueous%dtotal)) then
        psvd_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt  
        J(istart:iendaq,istart:iendaq) = J(istart:iendaq,istart:iendaq) + &
          rt_aux_var%aqueous%dtotal(:,:,iphase)*psvd_t
      else
        psvd_t = por*global_aux_var%sat(iphase)* &
          global_aux_var%den_kg(iphase)*vol/option%tran_dt ! units of den = kg water/m^3 water
        do icomp=istart,iendaq
          J(icomp,icomp) = J(icomp,icomp) + psvd_t
        enddo
      endif   
    endif
  enddo
#endif

end subroutine RTAccumulationDerivative

! ************************************************************************** !
!
! RCalculateCompression: Calculates the compression for the Jacobian block 
! author: Glenn Hammond
! date: 07/12/10
!
! ************************************************************************** !
subroutine RCalculateCompression(global_auxvar,rt_auxvar,reaction,option)

  use Option_module

  implicit none
  
  type(reaction_type), pointer :: reaction
  type(option_type) :: option

  PetscInt :: dfill(reaction%ncomp,reaction%ncomp)
  PetscInt :: ofill(reaction%ncomp,reaction%ncomp)
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  PetscReal :: residual(reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscInt :: i, jj
  PetscReal :: vol = 1.d0
  PetscReal :: por = 0.25d0
  PetscReal :: sum
  
  dfill = 0
  ofill = 0
  J = 0.d0
  residual = 0.d0

  call RTAuxVarCompute(rt_auxvar,global_auxvar,reaction,option)
  call RTAccumulationDerivative(rt_auxvar,global_auxvar, &
                                por,vol,reaction,option,J)
    
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (dabs(J(i,jj)) > 1.d-20) ofill(i,jj) = 1
    enddo
  enddo

  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    call RAccumulationSorbDerivative(rt_auxvar,global_auxvar,vol, &
                                     reaction,option,J)
  endif

  call RReaction(residual,J,PETSC_TRUE,rt_auxvar,global_auxvar,por,vol, &
                 reaction,option)
 
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (dabs(J(i,jj)) > 1.d-20) dfill(i,jj) = 1
    enddo
  enddo

  sum = 0.d0
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (dfill(i,jj) == 1) sum = sum + 1.d0
    enddo
  enddo
  write(option%io_buffer,'(''Diagonal Fill (%): '',f6.2)') &
    sum / (reaction%ncomp*reaction%ncomp) * 100.d0
  call printMsg(option)
 
 
  sum = 0.d0
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (ofill(i,jj) == 1) sum = sum + 1.d0
    enddo
  enddo
  write(option%io_buffer,'(''Off-Diagonal Fill (%): '',f6.2)') &
    sum / (reaction%ncomp*reaction%ncomp) * 100.d0
  call printMsg(option)

end subroutine RCalculateCompression

! ************************************************************************** !
!
! PrintRTAuxVar: Prints data from RTAuxVar object
! author: Glenn Hammond
! date: 05/18/2011
!
! ************************************************************************** !
subroutine RTPrintAuxVar(rt_auxvar,reaction,option)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  
  10 format(a20,':',10es19.11)
  20 format(a20,':',a20)
  30 format(/)

  if (OptionPrintToScreen(option)) write(*,30)
  if (OptionPrintToFile(option)) write(option%fid_out,30)

  if (OptionPrintToScreen(option)) &
    write(*,20) 'Primary', 'free molal., total molar., act. coef.'
  if (OptionPrintToFile(option)) &
    write(option%fid_out,20) 'Primary', 'free molal., total molar., act. coef.'
  do i = 1, reaction%naqcomp  
    if (OptionPrintToScreen(option)) &
      write(*,10) reaction%primary_species_names(i), &
        rt_auxvar%pri_molal(i), &
        rt_auxvar%total(i,1), &
        rt_auxvar%pri_act_coef(i)
    if (OptionPrintToFile(option)) &
      write(option%fid_out,10) reaction%primary_species_names(i), &
        rt_auxvar%pri_molal(i), &
        rt_auxvar%total(i,1), &
        rt_auxvar%pri_act_coef(i)
  enddo
  if (OptionPrintToScreen(option)) write(*,30)
  if (OptionPrintToFile(option)) write(option%fid_out,30)

!  aux_var%aqueous => MatrixBlockAuxVarCreate(option)
!  call MatrixBlockAuxVarInit(aux_var%aqueous,reaction%naqcomp, &
!                             reaction%naqcomp,option%nphase,option)
  
  if (reaction%neqcplx > 0) then
    if (OptionPrintToScreen(option)) &
      write(*,20) 'Secondary Complex', 'molal., act. coef.'
    if (OptionPrintToFile(option)) &
      write(option%fid_out,20) 'Secondary Complex', 'molal., act. coef.'
    do i = 1, reaction%neqcplx  
      if (OptionPrintToScreen(option)) &
        write(*,10) reaction%secondary_species_names(i), &
          rt_auxvar%sec_molal(i), &
          rt_auxvar%sec_act_coef(i)
      if (OptionPrintToFile(option)) &
        write(option%fid_out,10) reaction%secondary_species_names(i), &
          rt_auxvar%sec_molal(i), &
          rt_auxvar%sec_act_coef(i)
    enddo
    if (OptionPrintToScreen(option)) write(*,30)
    if (OptionPrintToFile(option)) write(option%fid_out,30)
  endif

  if (reaction%neqsorb > 0) then  
    if (OptionPrintToScreen(option)) &
      write(*,20) 'Total Sorbed', 'mol/m^3'
    if (OptionPrintToFile(option)) &
      write(option%fid_out,20) 'Total Sorbed', 'mol/m^3'
    do i = 1, reaction%naqcomp  
      if (OptionPrintToScreen(option)) &
        write(*,10) reaction%primary_species_names(i), rt_auxvar%total_sorb_eq(i)
      if (OptionPrintToFile(option)) &
        write(option%fid_out,10) reaction%primary_species_names(i), &
          rt_auxvar%total_sorb_eq(i)
    enddo
    if (OptionPrintToScreen(option)) write(*,30)
    if (OptionPrintToFile(option)) write(option%fid_out,30)
  endif    
  
  if (reaction%neqsrfcplx > 0) then
    if (OptionPrintToScreen(option)) &
      write(*,20) 'Surface Complex Conc.', 'mol/m^3'
    if (OptionPrintToFile(option)) &
      write(option%fid_out,20) 'Surface Complex Conc.', 'mol/m^3'
    do i = 1, reaction%neqsrfcplx
      if (OptionPrintToScreen(option)) &
        write(*,10) reaction%eqsrfcplx_names(i), rt_auxvar%eqsrfcplx_conc(i)
      if (OptionPrintToFile(option)) &
        write(option%fid_out,10) reaction%eqsrfcplx_names(i), &
          rt_auxvar%eqsrfcplx_conc(i)
    enddo
    if (OptionPrintToScreen(option)) write(*,30)
    if (OptionPrintToFile(option)) write(option%fid_out,30)
  endif

  if (reaction%nkinsrfcplxrxn > 0) then
  endif
  
  if (reaction%neqionxrxn > 0) then
  endif
  
  if (reaction%nkinmnrl > 0) then
    if (OptionPrintToScreen(option)) &
      write(*,20) 'Kinetic Minerals', 'vol frac, area, rate'
    if (OptionPrintToFile(option)) &
      write(option%fid_out,20) 'Kinetic Minerals', 'vol frac, area, rate'
    do i = 1, reaction%nkinmnrl
      if (OptionPrintToScreen(option)) &
        write(*,10) reaction%kinmnrl_names(i), &
          rt_auxvar%mnrl_volfrac(i), &
          rt_auxvar%mnrl_area(i), &
          rt_auxvar%mnrl_rate(i)

      if (OptionPrintToFile(option)) &
        write(option%fid_out,10) reaction%kinmnrl_names(i), &
          rt_auxvar%mnrl_volfrac(i), &
          rt_auxvar%mnrl_area(i), &
          rt_auxvar%mnrl_rate(i)
    enddo
    if (OptionPrintToScreen(option)) write(*,30)
    if (OptionPrintToFile(option)) write(option%fid_out,30)
  endif

end subroutine RTPrintAuxVar

end module Reaction_module
