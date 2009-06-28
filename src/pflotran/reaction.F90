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
            RTotal, &
            RTotalSorb, &
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
            RAccumulationSorbDerivative

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
  type(aq_species_type), pointer :: species, prev_species
  type(gas_species_type), pointer :: gas, prev_gas
  type(mineral_type), pointer :: mineral, prev_mineral
  type(surface_complex_type), pointer :: srfcmplx, prev_srfcmplx
  type(surface_complexation_rxn_type), pointer :: srfcmplx_rxn, &
                                                  prev_srfcmplx_rxn
  type(ion_exchange_rxn_type), pointer :: ionx_rxn, prev_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cation, prev_cation
  PetscInt :: srfcmplx_count

  nullify(prev_species)
  nullify(prev_gas)
  nullify(prev_mineral)
  nullify(prev_srfcmplx_rxn)
  nullify(prev_srfcmplx)
  nullify(prev_ionx_rxn)
  nullify(prev_cation)

  srfcmplx_count = 0
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
          
          reaction%ncomp = reaction%ncomp + 1
          
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
          
          reaction%neqcmplx = reaction%neqcmplx + 1
          
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
        call InputSkipToEnd(input,option,word)
      case('SORPTION')
        nullify(prev_srfcmplx_rxn)
        do
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,SORPTION')
          call StringToUpper(word)   

          select case(trim(word))
          
            case('SURFACE_COMPLEXATION_RXN')
          
              srfcmplx_rxn => SurfaceComplexationRXNCreate()
              do
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadWord(input,option,word,.true.)
                call InputErrorMsg(input,option,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN')
                call StringToUpper(word)
                
                select case(trim(word))
                  case('RATE','RATES') 
                    string = 'RATES inside SURFACE_COMPLEXATION_RXN'
                    call UtilityReadArray(reaction%kinmr_rate,-1,string,input,option) 
                    reaction%kinmr_nrate = size(reaction%kinmr_rate)
                  case('MINERAL')
                    call InputReadWord(input,option,srfcmplx_rxn%mineral_name,PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,MINERAL_NAME')
                  case('SITE')
                    call InputReadWord(input,option,srfcmplx_rxn%free_site_name,PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_NAME')
                    call InputReadDouble(input,option,srfcmplx_rxn%site_density)
                    call InputErrorMsg(input,option,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_DENSITY')                   
                  case('COMPLEXES')
                    nullify(prev_srfcmplx)
                    do
                      call InputReadFlotranString(input,option)
                      if (InputError(input)) exit
                      if (InputCheckExit(input,option)) exit
                      
                      srfcmplx_count = srfcmplx_count + 1
                      reaction%neqsurfcmplx = srfcmplx_count
                      srfcmplx => SurfaceComplexCreate()
                      srfcmplx%id = srfcmplx_count
                      call InputReadWord(input,option,srfcmplx%name,PETSC_TRUE)
                      call InputErrorMsg(input,option,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_NAME')
                
                      if (.not.associated(srfcmplx_rxn%complex_list)) then
                        srfcmplx_rxn%complex_list => srfcmplx
                      endif
                      if (associated(prev_srfcmplx)) then
                        prev_srfcmplx%next => srfcmplx
                      endif
                      prev_srfcmplx => srfcmplx
                      nullify(srfcmplx)
                
                    enddo
                  case default
                    option%io_buffer = 'CHEMISTRY, SURFACE_COMPLEXATION_RXN keyword: '// &
                                     trim(word)//' not recognized'
                    call printErrMsg(option)
                end select

              enddo
              if (.not.associated(reaction%surface_complexation_rxn_list)) then
                reaction%surface_complexation_rxn_list => srfcmplx_rxn
                srfcmplx_rxn%id = 1
              endif
              if (associated(prev_srfcmplx_rxn)) then
                prev_srfcmplx_rxn%next => srfcmplx_rxn
                srfcmplx_rxn%id = prev_srfcmplx_rxn%id + 1
              endif
              prev_srfcmplx_rxn => srfcmplx_rxn

              srfcmplx_rxn%free_site_id = srfcmplx_rxn%id
              reaction%neqsurfcmplxrxn = srfcmplx_rxn%id

              nullify(srfcmplx_rxn)

            case('ION_EXCHANGE_RXN')
            
              ionx_rxn => IonExchangeRxnCreate()
              do
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadWord(input,option,word,.true.)
                call InputErrorMsg(input,option,'keyword','CHEMISTRY,ION_EXCHANGE_RXN')
                call StringToUpper(word)
                
                select case(trim(word))
                  case('MINERAL')
                    call InputReadWord(input,option,ionx_rxn%mineral_name,PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,MINERAL_NAME')
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
                      call InputErrorMsg(input,option,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,CATION_NAME')
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
            case('DISTRIBUTION_COEF')
          !   call DistributionCoefRead
          end select
        enddo
      case('DATABASE')
        call InputReadNChars(input,option,reaction%database_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)  
        call InputErrorMsg(input,option,'keyword', &
                           'CHEMISTRY,DATABASE FILENAME')  
      case('LOG_FORMULATION')
        reaction%use_log_formulation = PETSC_TRUE        
      case('NO_CHECKPOINT_ACT_COEFS')
        reaction%checkpoint_activity_coefs = PETSC_FALSE
      case('ACTIVITY_COEFFICIENTS')
        reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_LAG        
        reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_TIMESTEP        
        do 
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr /= 0) exit
          select case(trim(word))
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
      case('MOLAL','MOLARITY')
        option%initialize_with_molality = PETSC_TRUE
      case('ACTIVITY_H2O','ACTIVITY_WATER')
        reaction%use_activity_h2o = PETSC_TRUE
      case('OUTPUT')
        call InputSkipToEnd(input,option,word)
      case('MAX_DLNC')
        call InputReadDouble(input,option,reaction%max_dlnC)
        call InputErrorMsg(input,option,trim(word),'CHEMISTRY')
      case default
        option%io_buffer = 'CHEMISTRY keyword: '//trim(word)//' not recognized'
        call printErrMsg(option)
    end select
  enddo
  
  reaction%neqsorb = reaction%neqsurfcmplxrxn + reaction%neqionxrxn

  if (reaction%neqcmplx + reaction%neqsorb + reaction%nmnrl > 0) then
    reaction%use_full_geochemistry = PETSC_TRUE
  endif

  
  if (len_trim(reaction%database_filename) < 2) &
    reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
 
end subroutine ReactionRead

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
                                        mineral_constraint,option)
  use Option_module
  use Input_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(option_type) :: option
  
  PetscTruth :: found
  PetscInt :: icomp, jcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: igas
  PetscReal :: constraint_conc(reaction%ncomp)
  PetscInt :: constraint_type(reaction%ncomp)
  character(len=MAXWORDLENGTH) :: constraint_spec_name(reaction%ncomp)
  character(len=MAXWORDLENGTH) :: constraint_mnrl_name(reaction%nkinmnrl)
  PetscInt :: constraint_id(reaction%ncomp)
    
  constraint_id = 0
  constraint_spec_name = ''
  constraint_type = 0
  constraint_conc = 0.d0
  
  ! aqueous species
  do icomp = 1, reaction%ncomp
    found = PETSC_FALSE
    do jcomp = 1, reaction%ncomp
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
      constraint_spec_name(jcomp) = aq_species_constraint%constraint_spec_name(icomp)
      constraint_conc(jcomp) = aq_species_constraint%constraint_conc(icomp)
      
      ! link constraint species
      select case(constraint_type(jcomp))
        case(CONSTRAINT_MINERAL)
          found = PETSC_FALSE
          do imnrl = 1, reaction%nmnrl
            if (StringCompare(constraint_spec_name(jcomp), &
                                reaction%mineral_names(imnrl), &
                                MAXWORDLENGTH)) then
              constraint_id(jcomp) = imnrl
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            option%io_buffer = 'Constraint mineral: ' // &
                     trim(constraint_spec_name(jcomp)) // &
                     ' for aqueous species: ' // &
                     trim(reaction%primary_species_names(jcomp)) // &
                     ' in constraint: ' // &
                     trim(constraint_name) // ' not found.' 
            call printErrMsg(option)         
          endif
        case(CONSTRAINT_GAS, CONSTRAINT_SUPERCRIT_CO2)
          found = PETSC_FALSE
          do igas = 1, reaction%ngas
            if (StringCompare(constraint_spec_name(jcomp), &
                                reaction%gas_species_names(igas), &
                                MAXWORDLENGTH)) then
              constraint_id(jcomp) = igas
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            option%io_buffer = 'Constraint gas: ' // &
                     trim(constraint_spec_name(jcomp)) // &
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
  aq_species_constraint%constraint_spec_name = constraint_spec_name
  aq_species_constraint%constraint_spec_id = constraint_id
  aq_species_constraint%constraint_conc = constraint_conc

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
  PetscInt :: num_iterations
  PetscTruth :: initialize_rt_auxvar
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icomp, jcomp, kcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: icplx
  PetscInt :: igas
  PetscReal :: conc(reaction%ncomp)
  PetscInt :: constraint_type(reaction%ncomp)
  character(len=MAXWORDLENGTH) :: constraint_spec_name(reaction%ncomp)

  PetscReal :: Res(reaction%ncomp)
  PetscReal :: total_conc(reaction%ncomp)
  PetscReal :: free_conc(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscInt :: indices(reaction%ncomp)
  PetscReal :: norm
  PetscReal :: prev_molal(reaction%ncomp)
  PetscReal, parameter :: tol = 1.d-12
  PetscReal, parameter :: tol_loose = 1.d-6
  PetscTruth :: compute_activity_coefs

  PetscInt :: constraint_id(reaction%ncomp)
  PetscReal :: lnQK, QK
  PetscReal :: tempreal
  PetscReal :: pres, tc, xphico2, henry, m_na, m_cl 
  PetscInt :: comp_id
  PetscReal :: convert_molal_to_molar
  PetscReal :: convert_molar_to_molal
  
  PetscTruth :: charge_balance_warning_flag = PETSC_FALSE

  PetscReal :: Jac_num(reaction%ncomp)
  PetscReal :: Res_pert, pert, prev_value

  PetscInt :: iphase
  PetscInt :: kinmr_nrate_store, irate

#ifdef CHUAN_CO2  
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
  PetscInt :: ierr
#endif
    
  constraint_type = aq_species_constraint%constraint_type
  constraint_spec_name = aq_species_constraint%constraint_spec_name
  constraint_id = aq_species_constraint%constraint_spec_id
  conc = aq_species_constraint%constraint_conc

  iphase = 1
  if (option%initialize_with_molality) then
    convert_molal_to_molar = global_auxvar%den_kg(iphase)/1000.d0
    convert_molar_to_molal = 1.d0
  else
    convert_molal_to_molar = 1.d0
    convert_molar_to_molal = 1000.d0/global_auxvar%den_kg(iphase)
  endif
  
  if (.not.reaction%use_full_geochemistry) then
!    aq_species_constraint%basis_molarity = conc*convert_molar_to_molal
    aq_species_constraint%basis_molarity = conc ! don't need to convert
    rt_auxvar%pri_molal = aq_species_constraint%basis_molarity
    rt_auxvar%total(:,iphase) = aq_species_constraint%basis_molarity
    return
  endif

  ! if using multirate reaction, we need to turn it off to equilibrate the system
  ! then turn it back on
  kinmr_nrate_store = 0
  if (reaction%kinmr_nrate > 0) then
    kinmr_nrate_store = reaction%kinmr_nrate
    reaction%kinmr_nrate = 0
    allocate(rt_auxvar%dtotal_sorb_eq(reaction%ncomp,reaction%ncomp))
  endif
  
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
  call ReactionInterpolateLogK(reaction%eqcmplx_logKcoef,reaction%eqcmplx_logK, &
                               global_auxvar%temp(iphase),reaction%neqcmplx)
  call ReactionInterpolateLogK(reaction%eqgas_logKcoef,reaction%eqgas_logK, &
                               global_auxvar%temp(iphase),reaction%ngas)
  call ReactionInterpolateLogK(reaction%eqsurfcmplx_logKcoef,reaction%eqsurfcmplx_logK, &
                               global_auxvar%temp(iphase),reaction%neqsurfcmplx)
  call ReactionInterpolateLogK(reaction%kinmnrl_logKcoef,reaction%kinmnrl_logK, &
                               global_auxvar%temp(iphase),reaction%nkinmnrl)
  call ReactionInterpolateLogK(reaction%mnrl_logKcoef,reaction%mnrl_logK, &
                               global_auxvar%temp(iphase),reaction%nmnrl)
  endif
#endif  
  
  total_conc = 0.d0
  do icomp = 1, reaction%ncomp
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
        ! check if h+ id set
        if (reaction%h_ion_id /= 0) then
          ! check if icomp is h+
          if (reaction%h_ion_id /= icomp) then
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
    endif
    call RTotal(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) call RTotalSorb(rt_auxvar,global_auxvar,reaction,option)
    
    Jac = 0.d0
        
    do icomp = 1, reaction%ncomp
      select case(constraint_type(icomp))
      
        case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
          ! units = mol/L water
          Res(icomp) = rt_auxvar%total(icomp,1) - total_conc(icomp)
          ! dtotal units = kg water/L water

          ! Jac units = kg water/L water
          Jac(icomp,:) = rt_auxvar%dtotal(icomp,:,1)
      
        case(CONSTRAINT_TOTAL_SORB)
          ! conversion from m^3 bulk -> L water
          tempreal = option%reference_porosity*option%reference_saturation*1000.d0
          ! total = mol/L water  total_sorb = mol/m^3 bulk
          Res(icomp) = rt_auxvar%total(icomp,1) + &
            rt_auxvar%total_sorb_eq(icomp)/tempreal - total_conc(icomp)
          ! dtotal units = kg water/L water
          ! dtotal_sorb units = kg water/m^3 bulk
          ! Jac units = kg water/L water
          Jac(icomp,:) = rt_auxvar%dtotal(icomp,:,1) + &
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
          do jcomp = 1, reaction%ncomp
            Res(icomp) = Res(icomp) + reaction%primary_spec_Z(jcomp) * &
              rt_auxvar%total(jcomp,1)
            do kcomp = 1, reaction%ncomp
              Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                reaction%primary_spec_Z(kcomp)*rt_auxvar%dtotal(kcomp,jcomp,1)
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
          if (reaction%h_ion_id > 0) then ! conc(icomp) = 10**-pH
            rt_auxvar%pri_molal(icomp) = 10.d0**(-conc(icomp)) / &
                                          rt_auxvar%pri_act_coef(icomp)
          else ! H+ is a complex
          
            icplx = abs(reaction%h_ion_id)
            
            ! compute secondary species concentration
            ! *note that the sign was flipped below
            lnQK = -reaction%eqcmplx_logK(icplx)*LOG_TO_LN

            ! activity of water
            if (reaction%eqcmplxh2oid(icplx) > 0) then
              lnQK = lnQK + reaction%eqcmplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
            endif

            do jcomp = 1, reaction%eqcmplxspecid(0,icplx)
              comp_id = reaction%eqcmplxspecid(jcomp,icplx)
              lnQK = lnQK + reaction%eqcmplxstoich(jcomp,icplx)* &
                            log(rt_auxvar%pri_molal(comp_id)* &
                            rt_auxvar%pri_act_coef(comp_id))
            enddo
            lnQK = lnQK - log(conc(icomp)) ! this is log activity H+
            QK = exp(lnQK)
            
            Res(icomp) = 1.d0 - QK

            do jcomp = 1,reaction%eqcmplxspecid(0,icplx)
              comp_id = reaction%eqcmplxspecid(jcomp,icplx)
              Jac(icomp,comp_id) = -exp(lnQK-log(rt_auxvar%pri_molal(comp_id)))* &
                                        reaction%eqcmplxstoich(jcomp,icplx)
            enddo
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

#ifdef CHUAN_CO2
            print *,'Gas CO2 constraint Jac,',igas, icomp, comp_id, &
              reaction%eqgasstoich(jcomp,igas),&
              Jac(icomp,comp_id), rt_auxvar%pri_molal(comp_id), lnQK
#endif
          enddo

#ifdef CHUAN_CO2        
        case(CONSTRAINT_SUPERCRIT_CO2)
          
          igas = constraint_id(icomp)
         
          ! compute secondary species concentration
          if(abs(reaction%co2_gas_id) == igas) then
!           pres = global_auxvar%pres(2)
!           pres = conc(icomp)*1.D5
            tc = global_auxvar%temp(1)

            call PSAT(tc, sat_pressure, ierr)
            
            pco2 = conc(icomp)*1.e5
!           pco2 = pres - sat_pressure
            
            pres = conc(icomp)*1.D5 + sat_pressure
            yco2 = pco2/pres
                        
!           call co2_span_wagner(pres*1D-6,tc+273.15D0,dg,dddt,dddp,fg, &
!             dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)

            call co2_span_wagner(pco2*1D-6,tc+273.15D0,dg,dddt,dddp,fg, &
              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)

            
            global_auxvar%den_kg(2) = dg
            
            !compute fugacity coefficient
            fg = fg*1.D6
            xphico2 = fg / pco2
            global_auxvar%fugacoeff(1) = xphico2
            
!           call Henry_duan_sun_0NaCl(pco2*1.d-5, tc, henry)
            if (reaction%na_ion_id /= 0 .and. reaction%cl_ion_id /= 0) then
              m_na = rt_auxvar%pri_molal(reaction%na_ion_id)
              m_cl = rt_auxvar%pri_molal(reaction%cl_ion_id)
              call Henry_duan_sun(tc,pco2*1D-5,henry,xphico2,lngamco2, &
                m_na,m_cl,sat_pressure*1D-5)
            else
              call Henry_duan_sun(tc,pco2*1D-5,henry,xphico2,lngamco2, &
                option%m_nacl,option%m_nacl,sat_pressure*1D-5)
            endif
            
            lnQk = -log(xphico2*henry)
!           lnQk = log(fg/henry)

            reaction%eqgas_logK(igas) = -lnQK*LN_TO_LOG
            
            print *, 'SC CO2 constraint',igas,pres,pco2,tc,xphico2,henry,lnQk,yco2, &
              lngamco2,m_na,m_cl,reaction%eqgas_logK(igas)
            
            ! activity of water
            if (reaction%eqgash2oid(igas) > 0) then
              lnQK = lnQK + reaction%eqgash2ostoich(igas)*rt_auxvar%ln_act_h2o
            endif
            do jcomp = 1, reaction%eqgasspecid(0,igas)
              comp_id = reaction%eqgasspecid(jcomp,igas)
              lnQK = lnQK + reaction%eqgasstoich(jcomp,igas)* &
                log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
                print *,'SC: ',rt_auxvar%pri_molal(comp_id), &
                  rt_auxvar%pri_act_coef(comp_id)
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
    do icomp = 1, reaction%ncomp
      norm = max(1.d0,maxval(abs(Jac(icomp,:))))
      norm = 1.d0/norm
      Res(icomp) = Res(icomp)*norm
      Jac(icomp,:) = Jac(icomp,:)*norm
    enddo

    ! for derivatives with respect to ln conc
    do icomp = 1, reaction%ncomp
      Jac(:,icomp) = Jac(:,icomp)*rt_auxvar%pri_molal(icomp)
    enddo

    call ludcmp(Jac,reaction%ncomp,indices,icomp)
    call lubksb(Jac,reaction%ncomp,indices,Res)

    prev_molal = rt_auxvar%pri_molal

    Res = dsign(1.d0,Res)*min(dabs(Res),5.d0)
      
    rt_auxvar%pri_molal = rt_auxvar%pri_molal*exp(-Res)

    num_iterations = num_iterations + 1
    
    if (mod(num_iterations,1000) == 0) then
100   format('Constraint iteration count has exceeded: ',i5)
      write(option%io_buffer,100) num_iterations
      call printMsg(option)
      do icomp=1,reaction%ncomp
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

  if (kinmr_nrate_store > 0) then
    reaction%kinmr_nrate = kinmr_nrate_store
    kinmr_nrate_store = 0
    do irate = 1, reaction%kinmr_nrate
      rt_auxvar%kinmr_total_sorb(:,irate) = rt_auxvar%total_sorb_eq/dble(reaction%kinmr_nrate)
    enddo
    deallocate(rt_auxvar%dtotal_sorb_eq)
    nullify(rt_auxvar%dtotal_sorb_eq)
  endif
  
  ! remember that a density of 1 kg/L was assumed, thus molal and molarity are equal
  ! do not scale by molal_to_molar since it could be 1.d0 if MOLAL flag set
  aq_species_constraint%basis_molarity = rt_auxvar%pri_molal* &
                                         global_auxvar%den_kg(option%liquid_phase)/ &
                                         1000.d0
  
  write(option%io_buffer,111) trim(constraint_name),num_iterations
  call printMsg(option)
111 format(' Equilibrate Constraint: ',a30,i4)

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
  PetscInt :: eqcmplxsort(reaction%neqcmplx+1)
  PetscInt :: eqcmplxid(reaction%neqcmplx+1)
  PetscInt :: eqminsort(reaction%nmnrl)
  PetscInt :: eqsurfcmplxsort(reaction%neqsurfcmplx+reaction%neqsurfcmplxrxn)
  PetscTruth :: finished, found
  PetscReal :: conc, conc2
  PetscReal :: lnQK(reaction%nmnrl), QK(reaction%nmnrl)
  PetscReal :: charge_balance, ionic_strength
  PetscReal :: percent(reaction%neqcmplx+1)
  PetscReal :: totj, retardation, kd
  PetscInt :: comp_id, jcomp
  PetscInt :: icount
  PetscInt :: iphase
  PetscReal :: bulk_vol_to_fluid_vol, molar_to_molal, molal_to_molar

  aq_species_constraint => constraint_coupler%aqueous_species
  mineral_constraint => constraint_coupler%minerals
  
  iphase = 1

  90 format(2x,76('-'))
  91 format(a)

  write(option%fid_out,'(/,''  Constraint: '',a)') &
    trim(constraint_coupler%constraint_name)

  rt_auxvar => constraint_coupler%rt_auxvar
  global_auxvar => constraint_coupler%global_auxvar
  
  global_auxvar%den_kg(iphase) = option%reference_water_density
  global_auxvar%temp(1) = option%reference_temperature
  global_auxvar%sat(iphase) = option%reference_saturation
  bulk_vol_to_fluid_vol = option%reference_porosity*option%reference_saturation*1000.d0

  molal_to_molar = global_auxvar%den_kg(iphase)/1000.d0
  molar_to_molal = 1.d0/molal_to_molar
    
  if (.not.reaction%use_full_geochemistry) then
    100 format(/,'  species       molality')  
    write(option%fid_out,100)
    101 format(2x,a12,es12.4)
    do icomp = 1, reaction%ncomp
      write(option%fid_out,101) reaction%primary_species_names(icomp), &
                                rt_auxvar%pri_molal(icomp)
    enddo
  else

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqcmplx_logKcoef,reaction%eqcmplx_logK, &
                                 global_auxvar%temp(iphase),reaction%neqcmplx)
    call ReactionInterpolateLogK(reaction%eqgas_logKcoef,reaction%eqgas_logK, &
                                 global_auxvar%temp(iphase),reaction%ngas)
    call ReactionInterpolateLogK(reaction%eqsurfcmplx_logKcoef,reaction%eqsurfcmplx_logK, &
                                 global_auxvar%temp(iphase),reaction%neqsurfcmplx)
    call ReactionInterpolateLogK(reaction%kinmnrl_logKcoef,reaction%kinmnrl_logK, &
                                 global_auxvar%temp(iphase),reaction%nkinmnrl)
    call ReactionInterpolateLogK(reaction%mnrl_logKcoef,reaction%mnrl_logK, &
                                 global_auxvar%temp(iphase),reaction%nmnrl)
  endif
#endif  

    200 format('')
    201 format(a20,i5)
    202 format(a20,f10.2)
    203 format(a20,f8.2)
    204 format(a20,es12.4)
    write(option%fid_out,90)
    write(option%fid_out,201) '      iterations: ', &
      constraint_coupler%num_iterations
    if (reaction%h_ion_id > 0) then
      write(option%fid_out,203) '              pH: ', &
        -log10(rt_auxvar%pri_molal(reaction%h_ion_id)* &
               rt_auxvar%pri_act_coef(reaction%h_ion_id))
    else if (reaction%h_ion_id < 0) then
      write(option%fid_out,203) '              pH: ', &
        -log10(rt_auxvar%sec_molal(abs(reaction%h_ion_id))* &
               rt_auxvar%sec_act_coef(abs(reaction%h_ion_id)))
    endif
    
    ionic_strength = 0.d0
    charge_balance = 0.d0
    do icomp = 1, reaction%ncomp
      charge_balance = charge_balance + rt_auxvar%total(icomp,1)* &
                                        reaction%primary_spec_Z(icomp)
      ionic_strength = ionic_strength + rt_auxvar%pri_molal(icomp)* &
        reaction%primary_spec_Z(icomp)*reaction%primary_spec_Z(icomp)
    enddo
    
    if (reaction%neqcmplx > 0) then    
      do i = 1, reaction%neqcmplx
        ionic_strength = ionic_strength + rt_auxvar%sec_molal(i)* &
                                          reaction%eqcmplx_Z(i)* &
                                          reaction%eqcmplx_Z(i)
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
#ifdef CHUAN_CO2
    if (global_auxvar%den_kg(2) > 0.d0) then
      write(option%fid_out,'(a20,f8.2,a9)') '     density CO2: ', &
        global_auxvar%den_kg(2),' [kg/m^3]'
      write(option%fid_out,'(a20,es12.4,a9)') '            xphi: ', &
        global_auxvar%fugacoeff(1)
    endif
#endif
    write(option%fid_out,90)

    102 format(/,'  species               molality    total       act coef  constraint')  
    write(option%fid_out,102)
    write(option%fid_out,90)
  
    103 format(2x,a20,es12.4,es12.4,es12.4,4x,a)
    do icomp = 1, reaction%ncomp
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
          string = aq_species_constraint%constraint_spec_name(icomp)
        case(CONSTRAINT_SUPERCRIT_CO2)
          string = 'SC ' // aq_species_constraint%constraint_spec_name(icomp)
      end select
      write(option%fid_out,103) reaction%primary_species_names(icomp), &
                                rt_auxvar%pri_molal(icomp), &
                                rt_auxvar%total(icomp,1)*molar_to_molal, &
                                rt_auxvar%pri_act_coef(icomp), &
                                trim(string)
    enddo 
  endif 
      
  if (reaction%neqcmplx > 0) then    
    ! sort complex concentrations from largest to smallest
    do i = 1, reaction%neqcmplx
      eqcmplxsort(i) = i
    enddo
    do
      finished = PETSC_TRUE
      do i = 1, reaction%neqcmplx-1
        icplx = eqcmplxsort(i)
        icplx2 = eqcmplxsort(i+1)
        if (rt_auxvar%sec_molal(icplx) < &
            rt_auxvar%sec_molal(icplx2)) then
          eqcmplxsort(i) = icplx2
          eqcmplxsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo
            
    110 format(/,'  complex               molality    act coef  logK')  
    write(option%fid_out,110)
    write(option%fid_out,90)
    111 format(2x,a20,es12.4,es12.4,2x,es12.4)
    do i = 1, reaction%neqcmplx ! for each secondary species
      icplx = eqcmplxsort(i)
      write(option%fid_out,111) reaction%secondary_species_names(icplx), &
                                rt_auxvar%sec_molal(icplx), &
                                rt_auxvar%sec_act_coef(icplx), &
                                reaction%eqcmplx_logK(icplx)
    enddo 

    !print speciation precentages
    write(option%fid_out,92)
    92 format(/)
    134 format(2x,'complex species       percent   molality')
    135 format(2x,'primary species: ',a20,2x,' total conc: ',1pe12.4)
    136 format(2x,a20,2x,f6.2,2x,1pe12.4,1p2e12.4)
    do icomp = 1, reaction%ncomp
    
      eqcmplxsort = 0
      eqcmplxid = 0
      percent = 0.d0
      totj = 0.d0
      
      icount = 0
      do icplx = 1, reaction%neqcmplx
        found = PETSC_FALSE
        do i = 1, reaction%ncomp
          if (reaction%eqcmplxspecid(i,icplx) == icomp) then
            icount = icount + 1
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (found) then
          eqcmplxid(icount) = icplx
          percent(icount) = dabs(rt_auxvar%sec_molal(icplx)* &
                                 reaction%eqcmplxstoich(i,icplx))
          totj = totj + percent(icount)
        endif
      enddo
      icount = icount + 1
      eqcmplxid(icount) = -icomp
      percent(icount) = rt_auxvar%pri_molal(icomp)
      totj = totj + percent(icount)
      percent = percent / totj
      
      eqcmplxsort = 0
      do i = 1, icount
        eqcmplxsort(i) = i
      enddo
      
      do
        finished = PETSC_TRUE
        do i = 1, icount-1
          icplx = eqcmplxsort(i)
          icplx2 = eqcmplxsort(i+1)
          if (percent(abs(icplx)) < percent(abs(icplx2))) then
            eqcmplxsort(i) = icplx2
            eqcmplxsort(i+1) = icplx
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
        j = eqcmplxsort(i)
        if (percent(j) < 0.0001d0) cycle
        icplx = eqcmplxid(j)
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
          
  if (reaction%neqsurfcmplxrxn > 0) then
    ! sort surface complex concentrations from largest to smallest
    ! note that we include free site concentrations; their ids negated
    do i = 1, reaction%neqsurfcmplx
      eqsurfcmplxsort(i) = i
    enddo
    do i = 1, reaction%neqsurfcmplxrxn
      eqsurfcmplxsort(reaction%neqsurfcmplx+i) = -i
    enddo
    do
      finished = PETSC_TRUE
      do i = 1, reaction%neqsurfcmplx+reaction%neqsurfcmplxrxn-1
        icplx = eqsurfcmplxsort(i)
        icplx2 = eqsurfcmplxsort(i+1)
        if (icplx > 0) then
          conc = rt_auxvar%eqsurfcmplx_conc(icplx)
        else
          conc = rt_auxvar%eqsurfcmplx_freesite_conc(-icplx)
        endif
        if (icplx2 > 0) then
          conc2 = rt_auxvar%eqsurfcmplx_conc(icplx2)
        else
          conc2 = rt_auxvar%eqsurfcmplx_freesite_conc(-icplx2)
        endif
        if (conc < conc2) then
          eqsurfcmplxsort(i) = icplx2
          eqsurfcmplxsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo
            
    120 format(/,'  surf complex          mol/m^3 blk logK')  
    write(option%fid_out,120)
    write(option%fid_out,90)
    121 format(2x,a20,es12.4,es12.4)
    122 format(2x,a20,es12.4,'  free site')
    do i = 1, reaction%neqsurfcmplx+reaction%neqsurfcmplxrxn
      icplx = eqsurfcmplxsort(i)
      if (icplx > 0) then
        write(option%fid_out,121) reaction%surface_complex_names(icplx), &
                                  rt_auxvar%eqsurfcmplx_conc(icplx), &
                                  reaction%eqsurfcmplx_logK(icplx)
      else
        write(option%fid_out,122) reaction%surface_site_names(-icplx), &
                                  rt_auxvar%eqsurfcmplx_freesite_conc(-icplx)
      endif
    enddo 
  endif

! retardation
  if (reaction%neqsurfcmplxrxn > 0) then
    write(option%fid_out,123)
    write(option%fid_out,90)
    do j = 1, reaction%ncomp
      retardation = 1.d0
      do irxn = 1, reaction%neqsurfcmplxrxn
        ncplx = reaction%eqsurfcmplx_rxn_to_complex(0,irxn)
        do i = 1, ncplx
          icplx = reaction%eqsurfcmplx_rxn_to_complex(i,irxn)
          ncomp = reaction%eqsurfcmplxspecid(0,icplx)
          do jj = 1, ncomp
            jcomp = reaction%eqsurfcmplxspecid(jj,icplx)
            if (j == jcomp) then
              if (rt_auxvar%total(j,iphase) /= 0.d0) &
              retardation = retardation + &
                            reaction%eqsurfcmplxstoich(jj,icplx)* &
                            rt_auxvar%eqsurfcmplx_conc(icplx)/ &
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
      do jcomp = 1, ncomp
        icomp = reaction%eqionx_rxn_cationid(jcomp,irxn)
        kd = rt_auxvar%eqionx_conc(icomp,irxn)/rt_auxvar%total(icomp,iphase) & 
                      /bulk_vol_to_fluid_vol
        write(option%fid_out,128) reaction%primary_species_names(icomp), &
          reaction%eqionx_rxn_k(jcomp,irxn), & 
          rt_auxvar%eqionx_conc(icomp,irxn), &
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
    do jcomp = 1, reaction%ncomp
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
    131 format(2x,a30,2x,f10.4,2x,1pe12.4)

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
      lnQK(igas) = -reaction%eqgas_logK(igas)*LOG_TO_LN
      
      ! divide K by RT
      !lnQK = lnQK - log((auxvar%temp+273.15d0)*IDEAL_GAS_CONST)
      
      ! activity of water
      if (reaction%eqgash2oid(igas) > 0) then
        lnQK(igas) = lnQK(igas) + reaction%eqgash2ostoich(igas)*rt_auxvar%ln_act_h2o
      endif

      do jcomp = 1, reaction%eqgasspecid(0,igas)
        comp_id = reaction%eqgasspecid(jcomp,igas)
        lnQK(igas) = lnQK(igas) + reaction%eqgasstoich(jcomp,igas)* &
                      log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
      enddo
      
      QK(igas) = exp(lnQK(igas))
          
      write(option%fid_out,133) reaction%gas_species_names(igas),lnQK(igas)*LN_TO_LOG, &
      QK(igas),reaction%eqgas_logK(igas)
    enddo
  endif

end subroutine ReactionPrintConstraint

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
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  
  type(mineral_type), pointer :: cur_mineral
  PetscInt :: imnrl

  cur_mineral => reaction%mineral_list
  do 
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -1*abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,name,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MINERAL_KINETICS')
    
    cur_mineral => reaction%mineral_list
    do 
      if (.not.associated(cur_mineral)) exit
      if (StringCompare(cur_mineral%name,name,MAXWORDLENGTH)) then
        cur_mineral%itype = MINERAL_KINETIC
        if (.not.associated(cur_mineral%tstrxn)) then
          cur_mineral%tstrxn => TransitionStateTheoryRxnCreate()
        endif
        ! read rate
        call InputReadDouble(input,option,cur_mineral%tstrxn%rate)
        call InputErrorMsg(input,option,'rate','CHEMISTRY,MINERAL_KINETICS')
        cur_mineral%id = abs(cur_mineral%id)
        reaction%nkinmnrl = reaction%nkinmnrl + 1
        exit
      endif
      cur_mineral => cur_mineral%next
    enddo
    
  enddo
  
  ! allocated kinetic mineral names
  if (reaction%nkinmnrl > 0) then
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
  PetscTruth :: found

  type(aq_species_type), pointer :: cur_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_type), pointer :: cur_mineral
  type(surface_complex_type), pointer :: cur_surfcplx
  type(surface_complexation_rxn_type), pointer :: cur_surfcplx_rxn
  
  nullify(cur_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_surfcplx)
  nullify(cur_surfcplx_rxn)
  
  reaction%print_all_species = PETSC_FALSE

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,OUTPUT,SPECIES_NAME')
    
    found = PETSC_FALSE
    
    word = name
    call StringToLower(word)
    if (StringCompare(word,'all',THREE_INTEGER)) then
      reaction%print_all_species = PETSC_TRUE
      reaction%print_pH = PETSC_TRUE
      found = PETSC_TRUE
    endif

    if (StringCompare(name,'pH',TWO_INTEGER)) then
      reaction%print_pH = PETSC_TRUE
      found = PETSC_TRUE
    endif

    if (StringCompare(word,'kd',TWO_INTEGER)) then
      reaction%print_kd = PETSC_TRUE
      found = PETSC_TRUE
    endif

    if (StringCompare(word,'total_sorbed',TWELVE_INTEGER)) then
      reaction%print_total_sorb = PETSC_TRUE
      found = PETSC_TRUE
    endif

    if (StringCompare(word,'free_ion',EIGHT_INTEGER)) then
      reaction%print_total_component = PETSC_FALSE
      reaction%print_free_ion = PETSC_TRUE
      found = PETSC_TRUE
    endif

    if (StringCompare(word,'activity_coefficients',TWO_INTEGER)) then
      reaction%print_act_coefs = PETSC_TRUE
    endif    

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
      cur_surfcplx_rxn => reaction%surface_complexation_rxn_list
      do
        if (.not.associated(cur_surfcplx_rxn)) exit
        if (StringCompare(name,cur_surfcplx_rxn%free_site_name,MAXWORDLENGTH)) then
          cur_surfcplx_rxn%free_site_print_me = PETSC_TRUE
          found = PETSC_TRUE
          exit
        endif
        if (.not.found) then
          cur_surfcplx => cur_surfcplx_rxn%complex_list
          do  
            if (.not.associated(cur_surfcplx)) exit
          if (StringCompare(name,cur_surfcplx%name,MAXWORDLENGTH)) then
            cur_surfcplx%print_me = PETSC_TRUE
            found = PETSC_TRUE
            exit
          endif
            cur_surfcplx => cur_surfcplx%next
          enddo
        endif
        cur_surfcplx_rxn => cur_surfcplx_rxn%next
      enddo  
    endif

    if (.not.found) then
      option%io_buffer = 'CHEMISTRY,OUTPUT species name: '//trim(name)// &
                         ' not found among chemical species'
      call printErrMsg(option)
    endif

  enddo

end subroutine ReactionReadOutput
      
! ************************************************************************** !
!
! RReaction: Computes reactions
! author: Glenn Hammond
! date: 09/30/08
!
! ************************************************************************** !
subroutine RReaction(Res,Jac,derivative,rt_auxvar,global_auxvar,volume, &
                     reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscTruth :: derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: volume

  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    call RAccumulationSorb(rt_auxvar,global_auxvar, &
                           volume,reaction,option,Res) 
  endif
   
  if (reaction%nkinmnrl > 0) then
    call RKineticMineral(Res,Jac,derivative,rt_auxvar,global_auxvar,volume, &
                         reaction,option)
  endif
  
  if (reaction%kinmr_nrate > 0) then
    call RMultiRateSorption(Res,Jac,derivative,rt_auxvar,global_auxvar, &
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
subroutine RReactionDerivative(Res,Jac,rt_auxvar,global_auxvar, &
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
  PetscReal :: volume
   
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp
  PetscReal :: Jac_dummy(reaction%ncomp,reaction%ncomp)
  PetscReal :: pert
  PetscTruth :: compute_derivative

  ! add new reactions in the 3 locations below

  if (.not.option%numerical_derivatives) then ! analytical derivative
  !if (PETSC_FALSE) then
    compute_derivative = PETSC_TRUE
    if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
      call RAccumulationSorbDerivative(rt_auxvar,global_auxvar, &
                                       volume,reaction,option,Jac)
    endif
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res,Jac,compute_derivative,rt_auxvar, &
                           global_auxvar,volume,reaction,option)
    endif
    if (reaction%kinmr_nrate > 0) then
      call RMultiRateSorption(Res,Jac,compute_derivative,rt_auxvar, &
                              global_auxvar,volume,reaction,option)
    endif    

    ! #1: add new reactions here

  else ! numerical derivative
    compute_derivative = PETSC_FALSE
    Res_orig = 0.d0
    option%iflag = 0 ! be sure not to allocate mass_balance array
    call RTAuxVarInit(rt_auxvar_pert,reaction,option)
    call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
    if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
      call RAccumulationSorb(rt_auxvar,global_auxvar, &
                             volume,reaction,option,Res_orig) 
    endif    
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                           global_auxvar,volume,reaction,option)
    endif
    if (reaction%kinmr_nrate > 0) then
      call RMultiRateSorption(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                              global_auxvar,volume,reaction,option)
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

      if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
        call RAccumulationSorb(rt_auxvar_pert,global_auxvar, &
                               volume,reaction,option,Res_pert) 
      endif       
      if (reaction%nkinmnrl > 0) then
        call RKineticMineral(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                             global_auxvar,volume,reaction,option)
      endif
      if (reaction%kinmr_nrate > 0) then
        call RMultiRateSorption(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                                global_auxvar,volume,reaction,option)
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
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)

  if (reaction%use_activity_h2o) then
    sum_pri_molal = 0.d0
    do j = 1, reaction%ncomp
      sum_pri_molal = sum_pri_molal + rt_auxvar%pri_molal(j)
    enddo
  endif

  if (reaction%act_coef_update_algorithm == ACT_COEF_ALGORITHM_NEWTON) then

    ln_conc = log(rt_auxvar%pri_molal)
    ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  
#ifdef TEMP_DEPENDENT_LOGK
    if (.not.option%use_isothermal) then
      call ReactionInterpolateLogK(reaction%eqcmplx_logKcoef,reaction%eqcmplx_logK, &
                               global_auxvar%temp(1),reaction%neqcmplx)
    endif
#endif  
  
  ! compute primary species contribution to ionic strength
    fpri = 0.d0
    sum_molality = 0.d0
    do j = 1, reaction%ncomp
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
      do icplx = 1, reaction%neqcmplx ! for each secondary species
        I = I + rt_auxvar%sec_molal(icplx)*reaction%eqcmplx_Z(icplx)* &
                                         reaction%eqcmplx_Z(icplx)
      enddo
      I = 0.5d0*I
      f = I
    
      if (abs(I-II) < 1.d-6*I) exit
    
      if (reaction%neqcmplx > 0) then
        didi = 0.d0
        sqrt_I = sqrt(I)
        do icplx = 1, reaction%neqcmplx
          if (abs(reaction%eqcmplx_Z(icplx)) > 0.d0) then
            sum = 0.5d0*reaction%debyeA*reaction%eqcmplx_Z(icplx)* &
            reaction%eqcmplx_Z(icplx) &
            /(sqrt_I*(1.d0+reaction%debyeB*reaction%eqcmplx_a0(icplx)*sqrt_I)**2) &
            -reaction%debyeBdot
            ncomp = reaction%eqcmplxspecid(0,icplx)
            do jcomp = 1, ncomp
              j = reaction%eqcmplxspecid(jcomp,icplx)
              if(abs(reaction%primary_spec_Z(j)) > 0.d0) then
                dgamdi = -0.5d0*reaction%debyeA*reaction%primary_spec_Z(j)**2/(sqrt_I* &
                (1.d0+reaction%debyeB*reaction%primary_spec_a0(j)*sqrt_I)**2)+ &
                reaction%debyeBdot 
                sum = sum + reaction%eqcmplxstoich(jcomp,icplx)*dgamdi
              endif
            enddo
            dcdi = rt_auxvar%sec_molal(icplx)*LOG_TO_LN*sum
            didi = didi+0.5d0*reaction%eqcmplx_Z(icplx)*reaction%eqcmplx_Z(icplx)*dcdi
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
        write(option%io_buffer,*) 'ionic strength negative! it =',it,' I= ',I,II,den,didi,dcdi,sum
        call printErrMsg(option)        
      endif
    
  ! compute activity coefficients
  ! primary species
      I = II
      sqrt_I = sqrt(I)
      do icomp = 1, reaction%ncomp
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
      do icplx = 1, reaction%neqcmplx
        if (abs(reaction%eqcmplx_Z(icplx)) > 0.d0) then
          rt_auxvar%sec_act_coef(icplx) = exp((-reaction%eqcmplx_Z(icplx)* &
                                        reaction%eqcmplx_Z(icplx)* &
                                        sqrt_I*reaction%debyeA/ &
                                        (1.d0+reaction%eqcmplx_a0(icplx)* &
                                        reaction%debyeB*sqrt_I)+ &
                                        reaction%debyeBdot*I)* &
                                        LOG_TO_LN)
        else
          rt_auxvar%sec_act_coef(icplx) = 1.d0
        endif
    
    ! compute secondary species concentration
        lnQK = -reaction%eqcmplx_logK(icplx)*LOG_TO_LN

    ! activity of water
        if (reaction%eqcmplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqcmplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
        endif

        ncomp = reaction%eqcmplxspecid(0,icplx)
        do jcomp = 1, ncomp
          icomp = reaction%eqcmplxspecid(jcomp,icplx)
          lnQK = lnQK + reaction%eqcmplxstoich(jcomp,icplx)*ln_act(icomp)
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
          write(option%io_buffer,*) 'activity of H2O negative! ln act H2O =',rt_auxvar%ln_act_h2o
          call printMsg(option)
        endif
      endif

    enddo
  
  else
  
  ! compute ionic strength
  ! primary species
    I = 0.d0
    do icomp = 1, reaction%ncomp
      I = I + rt_auxvar%pri_molal(icomp)*reaction%primary_spec_Z(icomp)* &
                                       reaction%primary_spec_Z(icomp)
    enddo
  
  ! secondary species
    do icplx = 1, reaction%neqcmplx ! for each secondary species
      I = I + rt_auxvar%sec_molal(icplx)*reaction%eqcmplx_Z(icplx)* &
                                       reaction%eqcmplx_Z(icplx)
    enddo
    I = 0.5d0*I
    sqrt_I = sqrt(I)
  
  ! compute activity coefficients
  ! primary species
    do icomp = 1, reaction%ncomp
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
    do icplx = 1, reaction%neqcmplx
      if (dabs(reaction%eqcmplx_Z(icplx)) > 1.d-10) then
        rt_auxvar%sec_act_coef(icplx) = exp((-reaction%eqcmplx_Z(icplx)* &
                                      reaction%eqcmplx_Z(icplx)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%eqcmplx_a0(icplx)* &
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
  
  PetscInt :: i, j, icplx, icomp, jcomp, iphase, ncomp, ieqgas,ierr
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: lnQK, tempreal
  PetscReal :: den_kg_per_L
  PetscReal :: pressure, temperature, xphico2, muco2, den, m_na, m_cl
  
#ifdef CHUAN_CO2  
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
#endif
  

  rt_auxvar%total =0D0 !debugging 
  iphase = 1           
  
  den_kg_per_L = global_auxvar%den_kg(iphase)*1.d-3              

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  rt_auxvar%total(:,iphase) = rt_auxvar%pri_molal(:)
  ! initialize derivatives
  rt_auxvar%dtotal = 0.d0
  do icomp = 1, reaction%ncomp
    rt_auxvar%dtotal(icomp,icomp,iphase) = 1.d0
  enddo
  
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqcmplx_logKcoef,reaction%eqcmplx_logK, &
                               global_auxvar%temp(iphase),reaction%neqcmplx)
  endif
#endif  
  
  do icplx = 1, reaction%neqcmplx ! for each secondary species
    ! compute secondary species concentration
    lnQK = -reaction%eqcmplx_logK(icplx)*LOG_TO_LN

    ! activity of water
    if (reaction%eqcmplxh2oid(icplx) > 0) then
      lnQK = lnQK + reaction%eqcmplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
    endif

    ncomp = reaction%eqcmplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = reaction%eqcmplxspecid(i,icplx)
      lnQK = lnQK + reaction%eqcmplxstoich(i,icplx)*ln_act(icomp)
    enddo
    rt_auxvar%sec_molal(icplx) = exp(lnQK)/rt_auxvar%sec_act_coef(icplx)
  
    ! add contribution to primary totals
    ! units of total = mol/L
    do i = 1, ncomp
      icomp = reaction%eqcmplxspecid(i,icplx)
      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
                                      reaction%eqcmplxstoich(i,icplx)* &
                                      rt_auxvar%sec_molal(icplx)
    enddo
    
    ! add contribution to derivatives of total with respect to free
    ! bear in mind that the water density portion is scaled below
    do j = 1, ncomp
      jcomp = reaction%eqcmplxspecid(j,icplx)
      tempreal = reaction%eqcmplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 rt_auxvar%sec_act_coef(icplx)
      do i = 1, ncomp
        icomp = reaction%eqcmplxspecid(i,icplx)
        rt_auxvar%dtotal(icomp,jcomp,iphase) = rt_auxvar%dtotal(icomp,jcomp,iphase) + &
                                               reaction%eqcmplxstoich(i,icplx)*tempreal
      enddo
    enddo
  enddo

  ! convert molality -> molarity
  ! unit of total = mol/L water
  rt_auxvar%total(:,iphase) = rt_auxvar%total(:,iphase)*den_kg_per_L
  
  ! units of dtotal = kg water/L water
  rt_auxvar%dtotal = rt_auxvar%dtotal*den_kg_per_L
  
 !*********** Add SC phase contribution ***************************  
#ifdef CHUAN_CO2

  iphase = 2           
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqgas_logKcoef,reaction%eqgas_logK, &
                               global_auxvar%temp(1),reaction%ngas)
  endif
#endif  

  if(iphase > option%nphase) return 
  rt_auxvar%total(:,iphase) = 0D0
  rt_auxvar%dtotal(:,:,iphase)=0D0
!  do icomp = 1, reaction%ncomp
!    rt_auxvar%dtotal(icomp,icomp,iphase) = 1.d0
!  enddo
    
  den_kg_per_L = global_auxvar%den_kg(iphase)*1.d-3     
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


      if(abs(reaction%co2_gas_id) == ieqgas )then
!          call Henry_duan_sun_0NaCl(pco2*1D-5, temperature, henry)
        if (reaction%na_ion_id /= 0 .and. reaction%cl_ion_id /= 0) then
          m_na = rt_auxvar%pri_molal(reaction%na_ion_id)
          m_cl = rt_auxvar%pri_molal(reaction%cl_ion_id)
          call Henry_duan_sun(temperature,pressure*1D-5,muco2,xphico2, &
                lngamco2,m_na,m_cl,sat_pressure*1D-5)
        else
          call Henry_duan_sun(temperature,pressure*1D-5,muco2,xphico2, &
                lngamco2,option%m_nacl,option%m_nacl,sat_pressure*1D-5)
        endif
        lnQk = - log(muco2)
           
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
          rt_auxvar%pri_act_coef(icomp)*exp(lnQK)*rt_auxvar%pri_molal(icomp)&
          /pressure /xphico2* den
   
      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
                                        reaction%eqgasstoich(1,ieqgas)* &
                                        rt_auxvar%gas_molal(ieqgas)
!       print *,'Ttotal',pressure, temperature, xphico2, den, lnQk,rt_auxvar%pri_molal(icomp),&
!        global_auxvar%sat(2),rt_auxvar%gas_molal(ieqgas)
   !     if(rt_auxvar%total(icomp,iphase) > den)rt_auxvar%total(icomp,iphase) = den* .99D0
   !     enddo

   ! contribute to %dtotal 
      tempreal = rt_auxvar%pri_act_coef(icomp)*exp(lnQK)/pressure /xphico2* den 
         rt_auxvar%dtotal(icomp,icomp,iphase) = rt_auxvar%dtotal(icomp,icomp,iphase) + &
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
  
end subroutine RTotalSorb

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
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: surfcmplx_conc(reaction%neqsurfcmplx)
  PetscReal :: dSx_dmi(reaction%ncomp)
  PetscReal :: dSi_dSx
  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: tol = 1.d-12
  PetscTruth :: one_more
  PetscReal :: res, dres_dfree_site, dfree_site_conc
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqsurfcmplx_logKcoef,reaction%eqsurfcmplx_logK, &
                               global_auxvar%temp(iphase),reaction%neqsurfcmplx)
  endif
#endif  

  ! Surface Complexation
  do irxn = 1, reaction%neqsurfcmplxrxn
  
    ncplx = reaction%eqsurfcmplx_rxn_to_complex(0,irxn)
    
    free_site_conc = rt_auxvar%eqsurfcmplx_freesite_conc(irxn)

    ! get a pointer to the first complex (there will always be at least 1)
    ! in order to grab free site conc
    one_more = PETSC_FALSE
    do

      total = free_site_conc
      ln_free_site = log(free_site_conc)
      do j = 1, ncplx
        icplx = reaction%eqsurfcmplx_rxn_to_complex(j,irxn)
        ! compute secondary species concentration
        lnQK = -reaction%eqsurfcmplx_logK(icplx)*LOG_TO_LN

        ! activity of water
        if (reaction%eqsurfcmplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqsurfcmplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
        endif

        lnQK = lnQK + reaction%eqsurfcmplx_free_site_stoich(icplx)* &
                      ln_free_site
      
        ncomp = reaction%eqsurfcmplxspecid(0,icplx)
        do i = 1, ncomp
          icomp = reaction%eqsurfcmplxspecid(i,icplx)
          lnQK = lnQK + reaction%eqsurfcmplxstoich(i,icplx)*ln_act(icomp)
        enddo
        surfcmplx_conc(icplx) = exp(lnQK)
        total = total + reaction%eqsurfcmplx_free_site_stoich(icplx)*surfcmplx_conc(icplx) 
        
      enddo
      
      if (one_more) exit
      
      if (reaction%eqsurfcmplx_rxn_stoich_flag(irxn)) then 
        ! stoichiometry for free sites in one of reactions is not 1, thus must
        ! use nonlinear iteration to solve
        res = reaction%eqsurfcmplx_rxn_site_density(irxn)-total
        
        dres_dfree_site = 1.d0

        do j = 1, ncplx

          icplx = reaction%eqsurfcmplx_rxn_to_complex(j,irxn)
          dres_dfree_site = dres_dfree_site + &
            reaction%eqsurfcmplx_free_site_stoich(icplx)* &
            surfcmplx_conc(icplx)/free_site_conc
        enddo

        dfree_site_conc = res / dres_dfree_site
        free_site_conc = free_site_conc - dfree_site_conc
      
        if (dabs(dfree_site_conc/free_site_conc) < tol) then
          one_more = PETSC_TRUE
        endif
      
      else
      
        total = total / free_site_conc
        free_site_conc = reaction%eqsurfcmplx_rxn_site_density(irxn) / total  
        
        one_more = PETSC_TRUE 
      
      endif

    enddo
    
    rt_auxvar%eqsurfcmplx_freesite_conc(irxn) = free_site_conc
 
    dSx_dmi = 0.d0
    tempreal = 0.d0
    do j = 1, ncplx
      icplx = reaction%eqsurfcmplx_rxn_to_complex(j,irxn)
      ncomp = reaction%eqsurfcmplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsurfcmplxspecid(i,icplx)
        ! numerator of 4.39
        dSx_dmi(icomp) = dSx_dmi(icomp) + reaction%eqsurfcmplxstoich(i,icplx)* &
                                          reaction%eqsurfcmplx_free_site_stoich(icplx)* &
                                          surfcmplx_conc(icplx)
      enddo
      ! denominator of 4.39
      tempreal = tempreal + reaction%eqsurfcmplx_free_site_stoich(icplx)* & 
                            reaction%eqsurfcmplx_free_site_stoich(icplx)* &
                            surfcmplx_conc(icplx)
    enddo 
    ! divide denominator by Sx
    tempreal = tempreal / free_site_conc
    ! add 1.d0 to denominator
    tempreal = tempreal + 1.d0
    ! divide numerator by denominator
    dSx_dmi = -dSx_dmi / tempreal
    ! convert from dlogm to dm
    dSx_dmi = dSx_dmi / rt_auxvar%pri_molal
 
    do k = 1, ncplx
      icplx = reaction%eqsurfcmplx_rxn_to_complex(k,irxn)

!     rt_auxvar%eqsurfcmplx_conc(icplx) = surfcmplx_conc(icplx)
      rt_auxvar%eqsurfcmplx_conc(k) = surfcmplx_conc(icplx)

      ncomp = reaction%eqsurfcmplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsurfcmplxspecid(i,icplx)
        rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + &
          reaction%eqsurfcmplxstoich(i,icplx)*surfcmplx_conc(icplx)
      enddo
      
      dSi_dSx = reaction%eqsurfcmplx_free_site_stoich(icplx)* &
                surfcmplx_conc(icplx)/ &
                free_site_conc

      do j = 1, ncomp
        jcomp = reaction%eqsurfcmplxspecid(j,icplx)
        tempreal = reaction%eqsurfcmplxstoich(j,icplx)*surfcmplx_conc(icplx) / &
                   rt_auxvar%pri_molal(jcomp)+ &
                   dSi_dSx*dSx_dmi(jcomp)
                  
        do i = 1, ncomp
          icomp = reaction%eqsurfcmplxspecid(i,icplx)
          rt_auxvar%dtotal_sorb_eq(icomp,jcomp) = rt_auxvar%dtotal_sorb_eq(icomp,jcomp) + &
                                               reaction%eqsurfcmplxstoich(i,icplx)* &
                                               tempreal
        enddo
      enddo
    enddo
  enddo
  
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
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscReal, parameter :: tol = 1.d-12
  PetscTruth :: one_more
  PetscReal :: res
    
  PetscReal :: omega
  PetscReal :: ref_cation_X, ref_cation_conc, ref_cation_Z, ref_cation_k, &
               ref_cation_quotient
  PetscReal :: cation_X(reaction%ncomp)
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
      
      rt_auxvar%eqionx_conc(icomp,irxn) = rt_auxvar%eqionx_conc(icomp,irxn) + tempreal1

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

  PetscTruth :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: volume
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(option_type) :: option
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, ncomp, ncplx
  PetscInt, parameter :: iphase = 1
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: surfcmplx_conc(reaction%neqsurfcmplx)
  PetscReal :: dSx_dmi(reaction%ncomp)
  PetscReal :: dSi_dSx
  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscReal, parameter :: tol = 1.d-12
  PetscTruth :: one_more
  PetscReal :: residual, dres_dfree_site, dfree_site_conc
  PetscReal :: site_density
  
  PetscInt :: irate
  PetscReal :: kdt, one_plus_kdt, k_over_one_plus_kdt
  PetscReal :: total_sorb_eq(reaction%ncomp)
  PetscReal :: dtotal_sorb_eq(reaction%ncomp,reaction%ncomp)

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
    
#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqsurfcmplx_logKcoef,reaction%eqsurfcmplx_logK, &
                               global_auxvar%temp(iphase),reaction%neqsurfcmplx)
  endif
#endif  

  rt_auxvar%total_sorb_eq = 0.d0
  rt_auxvar%eqsurfcmplx_conc = 0.d0

  ! Surface Complexation
  do irxn = 1, reaction%neqsurfcmplxrxn
  
    ! the below assumes equal site density for each multi-rate reaction
    site_density = reaction%eqsurfcmplx_rxn_site_density(irxn)/dble(reaction%kinmr_nrate)
  
    ncplx = reaction%eqsurfcmplx_rxn_to_complex(0,irxn)
    free_site_conc = rt_auxvar%eqsurfcmplx_freesite_conc(irxn)

    ! get a pointer to the first complex (there will always be at least 1)
    ! in order to grab free site conc
    one_more = PETSC_FALSE
    do
      total = free_site_conc
      ln_free_site = log(free_site_conc)
      do j = 1, ncplx
        icplx = reaction%eqsurfcmplx_rxn_to_complex(j,irxn)
        ! compute secondary species concentration
        lnQK = -reaction%eqsurfcmplx_logK(icplx)*LOG_TO_LN

        ! activity of water
        if (reaction%eqsurfcmplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqsurfcmplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
        endif

        lnQK = lnQK + reaction%eqsurfcmplx_free_site_stoich(icplx)* &
                      ln_free_site
      
        ncomp = reaction%eqsurfcmplxspecid(0,icplx)
        do i = 1, ncomp
          icomp = reaction%eqsurfcmplxspecid(i,icplx)
          lnQK = lnQK + reaction%eqsurfcmplxstoich(i,icplx)*ln_act(icomp)
        enddo
        surfcmplx_conc(icplx) = exp(lnQK)
        total = total + reaction%eqsurfcmplx_free_site_stoich(icplx)*surfcmplx_conc(icplx) 
        
      enddo
      
      if (one_more) exit
      
      if (reaction%eqsurfcmplx_rxn_stoich_flag(irxn)) then 
        ! stoichiometry for free sites in one of reactions is not 1, thus must
        ! use nonlinear iteration to solve
        residual = site_density-total
        
        dres_dfree_site = 1.d0

        do j = 1, ncplx
          icplx = reaction%eqsurfcmplx_rxn_to_complex(j,irxn)
          dres_dfree_site = dres_dfree_site + &
            reaction%eqsurfcmplx_free_site_stoich(icplx)* &
            surfcmplx_conc(icplx)/free_site_conc
        enddo

        dfree_site_conc = residual / dres_dfree_site
        free_site_conc = free_site_conc - dfree_site_conc
      
        if (dabs(dfree_site_conc/free_site_conc) < tol) then
          one_more = PETSC_TRUE
        endif
      
      else
      
        total = total / free_site_conc
        free_site_conc = site_density / total  
        
        one_more = PETSC_TRUE
      endif
    enddo

    rt_auxvar%eqsurfcmplx_freesite_conc(irxn) = free_site_conc
   
    dSx_dmi = 0.d0
    tempreal = 0.d0
    do j = 1, ncplx
      icplx = reaction%eqsurfcmplx_rxn_to_complex(j,irxn)
      ncomp = reaction%eqsurfcmplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsurfcmplxspecid(i,icplx)
        ! numerator of 4.39
        dSx_dmi(icomp) = dSx_dmi(icomp) + reaction%eqsurfcmplxstoich(i,icplx)* &
          reaction%eqsurfcmplx_free_site_stoich(icplx)*surfcmplx_conc(icplx)
      enddo
      ! denominator of 4.39
      tempreal = tempreal + reaction%eqsurfcmplx_free_site_stoich(icplx)* & 
        reaction%eqsurfcmplx_free_site_stoich(icplx)*surfcmplx_conc(icplx)
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
      icplx = reaction%eqsurfcmplx_rxn_to_complex(k,irxn)

      rt_auxvar%eqsurfcmplx_conc(k) = &
        rt_auxvar%eqsurfcmplx_conc(k) + surfcmplx_conc(icplx)

      ncomp = reaction%eqsurfcmplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsurfcmplxspecid(i,icplx)
        total_sorb_eq(icomp) = total_sorb_eq(icomp) + &
          reaction%eqsurfcmplxstoich(i,icplx)*surfcmplx_conc(icplx)
      enddo
      
      if (compute_derivative) then
        dSi_dSx = reaction%eqsurfcmplx_free_site_stoich(icplx)* &
          surfcmplx_conc(icplx)/free_site_conc

        do j = 1, ncomp
          jcomp = reaction%eqsurfcmplxspecid(j,icplx)
          tempreal = reaction%eqsurfcmplxstoich(j,icplx)*surfcmplx_conc(icplx) / &
            rt_auxvar%pri_molal(jcomp)+dSi_dSx*dSx_dmi(jcomp)
                      
          do i = 1, ncomp
            icomp = reaction%eqsurfcmplxspecid(i,icplx)
            dtotal_sorb_eq(icomp,jcomp) = dtotal_sorb_eq(icomp,jcomp) + &
              reaction%eqsurfcmplxstoich(i,icplx)*tempreal
          enddo
        enddo
      endif
    enddo
  enddo
      
  ! WARNING: this assumes equal site distribution 
  do irate = 1, reaction%kinmr_nrate
    kdt = reaction%kinmr_rate(irate) * option%tran_dt
    one_plus_kdt = 1.d0 + kdt
    k_over_one_plus_kdt = reaction%kinmr_rate(irate)/one_plus_kdt
        
    Res(:) = Res(:) + volume * k_over_one_plus_kdt * &
      (total_sorb_eq(:) - rt_auxvar%kinmr_total_sorb(:,irate))
      
    if (compute_derivative) then
      Jac = Jac + volume * k_over_one_plus_kdt * dtotal_sorb_eq
    endif

  enddo
  
  ! store the target equilibrium concentration to update the sorbed 
  ! concentration at the end of the time step.
  rt_auxvar%total_sorb_eq = total_sorb_eq
  
end subroutine RMultiRateSorption

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
  PetscTruth :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscInt :: i, j, k, imnrl, icomp, jcomp, kcplx, iphase, ncomp, ipref
  PetscReal :: prefactor(10), sum_prefactor_rate
  PetscReal :: dIm_dsum_prefactor_rate, dIm_dprefactor_rate
  PetscReal :: dprefactor_dcomp_numerator, dprefactor_dcomp_denominator
  PetscReal :: tempreal, tempreal2
  PetscReal :: affinity_factor, sign_
  PetscReal :: Im, Im_const, dIm_dQK
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_sec(reaction%neqcmplx)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: ln_sec_act(reaction%neqcmplx)
  PetscReal :: QK, lnQK, dQK_dCj, dQK_dmj
  PetscTruth :: prefactor_exists

  iphase = 1                         

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  if (reaction%neqcmplx > 0) then
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
    ! compute secondary species concentration
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
    QK = exp(lnQK)
    
    if (associated(reaction%kinmnrl_Tempkin_const)) then
      affinity_factor = 1.d0-QK**(1/reaction%kinmnrl_Tempkin_const(imnrl))
    else
      affinity_factor = 1.d0-QK
    endif
    
    sign_ = sign(1.d0,affinity_factor)

    if (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then
      ! compute prefactor
      if (reaction%kinmnrl_num_prefactors(imnrl) > 0) then
        print *, 'Kinetic mineral reaction prefactor calculations have not been verified.  Ask Glenn.'
        stop
        sum_prefactor_rate = 0
        do ipref = 1, reaction%kinmnrl_num_prefactors(imnrl)
          prefactor(ipref) = 1.d0
          do i = 1, reaction%kinmnrl_pri_prefactor_id(0,ipref,imnrl) ! primary contribution
            icomp = reaction%kinmnrl_pri_prefactor_id(i,ipref,imnrl)
            prefactor(ipref) = prefactor(ipref) * &
                               exp(reaction%kinmnrl_pri_pref_alpha_stoich(i,ipref,imnrl)* &
                                   ln_act(icomp))/ &
                               ((1.d0+reaction%kinmnrl_pri_pref_atten_coef(i,ipref,imnrl))* &
                                 exp(reaction%kinmnrl_pri_pref_beta_stoich(i,ipref,imnrl)* &
                                     ln_act(icomp)))
          enddo
          if (reaction%neqcmplx > 0) then
            do k = 1, reaction%kinmnrl_sec_prefactor_id(0,ipref,imnrl) ! secondary contribution
              kcplx = reaction%kinmnrl_sec_prefactor_id(k,ipref,imnrl)
              prefactor(ipref) = prefactor(ipref) * &
                                 exp(reaction%kinmnrl_sec_pref_alpha_stoich(k,ipref,imnrl)* &
                                     ln_sec_act(kcplx))/ &
                                 ((1.d0+reaction%kinmnrl_sec_pref_atten_coef(i,ipref,imnrl))* &
                                   exp(reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                       ln_sec_act(kcplx)))
            enddo
          endif
          sum_prefactor_rate = sum_prefactor_rate + prefactor(ipref)*reaction%kinmnrl_rate(ipref,imnrl)
        enddo
      else
        sum_prefactor_rate = reaction%kinmnrl_rate(1,imnrl)
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
      dIm_dQK = dIm_dQK*(1/reaction%kinmnrl_Tempkin_const(imnrl))/QK
    endif
    
    ! derivatives with respect to primary species in reaction quotient
    do j = 1, ncomp
      jcomp = reaction%kinmnrlspecid(j,imnrl)
      ! unit = L water/mol
      dQK_dCj = reaction%kinmnrlstoich(j,imnrl)*exp(lnQK-ln_conc(jcomp))
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

    if (reaction%kinmnrl_num_prefactors(imnrl) > 0) then ! add contribution of derivative in prefactor - messy
      print *, 'Kinetic mineral reaction prefactor calculations have not been verified.  Ask Glenn.'
      stop
      
      dIm_dsum_prefactor_rate = Im/sum_prefactor_rate
      do ipref = 1, reaction%kinmnrl_num_prefactors(imnrl)
        dIm_dprefactor_rate = dIm_dsum_prefactor_rate*reaction%kinmnrl_rate(ipref,imnrl)
        do j = 1, reaction%kinmnrl_pri_prefactor_id(0,ipref,imnrl) ! primary contribution
          jcomp = reaction%kinmnrl_pri_prefactor_id(j,ipref,imnrl)
          ! numerator
          dprefactor_dcomp_numerator = reaction%kinmnrl_pri_pref_alpha_stoich(j,ipref,imnrl)* &
                                       prefactor(ipref)/rt_auxvar%pri_molal(jcomp) ! dR_dm
          ! denominator
          dprefactor_dcomp_denominator = -prefactor(ipref)/ &
                                         ((1.d0+reaction%kinmnrl_pri_pref_atten_coef(j,ipref,imnrl))* &
                                           exp(reaction%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)* &
                                               ln_act(jcomp)))* & 
                                         reaction%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)* &
                                         reaction%kinmnrl_pri_pref_atten_coef(j,ipref,imnrl)* &
                                         exp((reaction%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)-1.d0)* &
                                             ln_act(jcomp))* & ! dR_da
                                         rt_auxvar%pri_act_coef(jcomp) ! da_dc
          tempreal = dIm_dprefactor_rate*(dprefactor_dcomp_numerator+ &
                     dprefactor_dcomp_denominator)*global_auxvar%den_kg(iphase)
          do i = 1, ncomp
            icomp = reaction%kinmnrlspecid(i,imnrl)
            Jac(icomp,jcomp) = Jac(icomp,jcomp) + reaction%kinmnrlstoich(i,imnrl)*tempreal
          enddo  ! loop over col
        enddo !loop over row
        if (reaction%neqcmplx > 0) then
          do k = 1, reaction%kinmnrl_sec_prefactor_id(0,ipref,imnrl) ! secondary contribution
            kcplx = reaction%kinmnrl_sec_prefactor_id(k,ipref,imnrl)
            ! numerator
            dprefactor_dcomp_numerator = reaction%kinmnrl_sec_pref_alpha_stoich(k,ipref,imnrl)* &
                                         prefactor(ipref)/(rt_auxvar%sec_molal(kcplx)* &
                                                           rt_auxvar%sec_act_coef(kcplx)) ! dR_dax
            ! denominator
            dprefactor_dcomp_denominator = -prefactor(ipref)/ &
                                           (1.d0+reaction%kinmnrl_sec_pref_atten_coef(k,ipref,imnrl)* &
                                            exp(reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                                ln_sec_act(kcplx)))* &
                                           reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                           reaction%kinmnrl_sec_pref_atten_coef(k,ipref,imnrl)* &
                                           exp((reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)-1.d0)* &
                                                ln_sec_act(kcplx)) ! dR_dax
            tempreal = dIm_dprefactor_rate*(dprefactor_dcomp_numerator+ &
                       dprefactor_dcomp_denominator)*global_auxvar%den_kg(iphase)
            do j = 1, reaction%eqcmplxspecid(0,kcplx)
              jcomp = reaction%eqcmplxspecid(j,kcplx)
              tempreal2 = reaction%eqcmplxstoich(j,kcplx)*exp(ln_sec_act(kcplx)-ln_conc(jcomp)) !dax_dc
              do i = 1, ncomp
                icomp = reaction%kinmnrlspecid(i,imnrl)
                Jac(icomp,jcomp) = Jac(icomp,jcomp) + reaction%kinmnrlstoich(i,imnrl)*tempreal* &
                                                      tempreal2
              enddo  ! loop over col
            enddo  ! loop over row
          enddo  ! loop over complexes
        endif
      enddo  ! loop over prefactors
    endif
  enddo  ! loop over minerals
    
end subroutine RKineticMineral

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
  Res(:) = Res(:) + rt_aux_var%total_sorb_eq(:)*v_t

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
  J = J + rt_aux_var%dtotal_sorb_eq(:,:)*v_t

end subroutine RAccumulationSorbDerivative

! ************************************************************************** !
!
! RSolve: Computes the kinetic mineral precipitation/dissolution
!                  rates
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RSolve(Res,Jac,conc,update,ncomp)

  use Utility_module
  
  implicit none

  PetscInt :: ncomp
  PetscReal :: Res(ncomp)
  PetscReal :: Jac(ncomp,ncomp)
  PetscReal :: update(ncomp)
  PetscReal :: conc(ncomp)
  
  PetscInt :: indices(ncomp)
  PetscInt :: icomp
  PetscReal :: norm

  ! scale Jacobian
  do icomp = 1, ncomp
    norm = max(1.d0,maxval(abs(Jac(icomp,:))))
    norm = 1.d0/norm
    res(icomp) = res(icomp)*norm
    Jac(icomp,:) = Jac(icomp,:)*norm
  enddo
    
  ! for derivatives with respect to ln conc
  do icomp = 1, ncomp
    Jac(:,icomp) = Jac(:,icomp)*conc(icomp)
  enddo
  call ludcmp(Jac,ncomp,indices,icomp)
  call lubksb(Jac,ncomp,indices,res)

  update = dsign(1.d0,res)*min(dabs(res),5.d0)
  
end subroutine RSolve

! ************************************************************************** !
!
! ReactionFitLogKCoef: Least squares fit to log K over database temperature range
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
subroutine ReactionFitLogKCoef(coefs,logK,option,reaction)

  use Option_module
  use Utility_module

  implicit none
  
  type(reaction_type) :: reaction
  PetscReal :: coefs(FIVE_INTEGER)
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
      else if (logK(i) .gt. 500.) then
        option%io_buffer = 'In ReactionFitLogKCoef: log K .gt. 500---stop!'
        call printErrMsg(option)
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
    do irxn = 1, reaction%neqsurfcmplxrxn
      do i = 1, reaction%eqsurfcmplx_rxn_to_complex(0,irxn)
        icplx = reaction%eqsurfcmplx_rxn_to_complex(i,irxn)
        do j = 1, reaction%eqsurfcmplxspecid(0,icplx)
          jcomp = reaction%eqsurfcmplxspecid(j,icplx)
          if (icomp == jcomp) then
            retardation = retardation + &
              reaction%eqsurfcmplxstoich(j,icplx) * &
              rt_auxvar%eqsurfcmplx_conc(icplx)
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

end module Reaction_module
