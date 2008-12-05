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
            RTotal, &
            RTotalSorb, &
            RActivityCoefficients, &
            RReaction, &
            RReactionDerivative, &
            ReactionProcessConstraint, &
            ReactionEquilibrateConstraint, &
            ReactionPrintConstraint

contains

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
  type(surface_complex_type), pointer :: srfcmplx, prev_srfcmplx
  type(surface_complexation_rxn_type), pointer :: srfcmplx_rxn, prev_srfcmplx_rxn
  type(ion_exchange_rxn_type), pointer :: ionx_rxn, prev_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cation, prev_cation
  PetscInt :: length
  PetscInt :: srfcmplx_count = 0
  PetscErrorCode :: ierr

  nullify(prev_srfcmplx_rxn)
  nullify(prev_ionx_rxn)

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)
    if (ierr /= 0) exit
    if (fiCheckExit(string)) exit
    
    call fiReadWord(string,word,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY', ierr)
    call fiWordToUpper(word)   

    select case(trim(word))
    
      case('PRIMARY_SPECIES')
        nullify(prev_species)
        do
          call fiReadFlotranString(fid,string,ierr)
          if (ierr /= 0) exit
          if (fiCheckExit(string)) exit
          
          reaction%ncomp = reaction%ncomp + 1
          
          species => AqueousSpeciesCreate()
          call fiReadWord(string,species%name,PETSC_TRUE,ierr)  
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
          if (fiCheckExit(string)) exit
          
          reaction%neqcmplx = reaction%neqcmplx + 1
          
          species => AqueousSpeciesCreate()
          call fiReadWord(string,species%name,PETSC_TRUE,ierr)  
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
          if (fiCheckExit(string)) exit

          reaction%ngas = reaction%ngas + 1
          
          gas => GasSpeciesCreate()
          call fiReadWord(string,gas%name,PETSC_TRUE,ierr)  
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
          if (fiCheckExit(string)) exit
          
          reaction%nmnrl = reaction%nmnrl + 1
          
          mineral => MineralCreate()
          call fiReadWord(string,mineral%name,PETSC_TRUE,ierr)  
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,MINERALS', ierr)    
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
        call fiSkipToEND(fid,option%myrank,word)
      case('SORPTION')
        nullify(prev_srfcmplx_rxn)
        do
          call fiReadFlotranString(fid,string,ierr)
          if (ierr /= 0) exit
          if (fiCheckExit(string)) exit

          call fiReadWord(string,word,PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SORPTION', ierr)
          call fiWordToUpper(word)   

          select case(trim(word))
          
            case('SURFACE_COMPLEXATION_RXN')
          
              srfcmplx_rxn => SurfaceComplexationRXNCreate()
              do
                call fiReadFlotranString(fid,string,ierr)
                if (ierr /= 0) exit
                if (fiCheckExit(string)) exit

                call fiReadWord(string,word,.true.,ierr)
                call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN', ierr)
                call fiWordToUpper(word)
                
                select case(trim(word))
                  case('MINERAL')
                    call fiReadWord(string,srfcmplx_rxn%mineral_name,PETSC_TRUE,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,MINERAL_NAME', ierr)
                  case('SITE')
                    call fiReadWord(string,srfcmplx_rxn%free_site_name,PETSC_TRUE,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_NAME', ierr)
                    call fiReadDouble(string,srfcmplx_rxn%site_density,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_DENSITY',ierr)                   
                  case('COMPLEXES')
                    nullify(prev_srfcmplx)
                    do
                      call fiReadFlotranString(fid,string,ierr)
                      if (ierr /= 0) exit
                      if (fiCheckExit(string)) exit
                      
                      srfcmplx_count = srfcmplx_count + 1
                      reaction%neqsurfcmplx = srfcmplx_count
                      srfcmplx => SurfaceComplexCreate()
                      srfcmplx%id = srfcmplx_count
                      call fiReadWord(string,srfcmplx%name,PETSC_TRUE,ierr)
                      call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_NAME', ierr)
                
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
                    call printErrMsg(option,'CHEMISTRY, SURFACE_COMPLEXATION_RXN keyword: '// &
                                     trim(word)//' not recognized')
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
                call fiReadFlotranString(fid,string,ierr)
                if (ierr /= 0) exit
                if (fiCheckExit(string)) exit

                call fiReadWord(string,word,.true.,ierr)
                call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,ION_EXCHANGE_RXN', ierr)
                call fiWordToUpper(word)
                
                select case(trim(word))
                  case('MINERAL')
                    call fiReadWord(string,ionx_rxn%mineral_name,PETSC_TRUE,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,MINERAL_NAME', ierr)
                  case('CEC')
                    call fiReadDouble(string,ionx_rxn%CEC,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,CEC',ierr)                   
                  case('CATIONS')
                    nullify(prev_cation)
                    do
                      call fiReadFlotranString(fid,string,ierr)
                      if (ierr /= 0) exit
                      if (fiCheckExit(string)) exit
                      
                      cation => IonExchangeCationCreate()
                      reaction%neqionxcation = reaction%neqionxcation + 1
                      call fiReadWord(string,cation%name,PETSC_TRUE,ierr)
                      call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,CATION_NAME', ierr)
                      call fiReadDouble(string,cation%k,ierr)
                      call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,ION_EXCHANGE_RXN,K',ierr)                   
    
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
                    call printErrMsg(option,'CHEMISTRY, ION_EXCHANGE_RXN keyword: '// &
                                     trim(word)//' not recognized')
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
        call fiReadNChars(string,reaction%database_filename,MAXSTRINGLENGTH,PETSC_TRUE,ierr)  
        call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,DATABASE FILENAME', ierr)  
      case('LOG_FORMULATION')
        reaction%use_log_formulation = PETSC_TRUE        
      case('ACTIVITY_COEFFICIENTS')
        call fiReadWord(string,word,PETSC_TRUE,ierr)
        call fiDefaultMsg(option%myrank,'CHEMISTRY,ACTIVITY COEFFICIENTS',ierr)        
        select case(trim(word))
          case('ITERATION')
            reaction%compute_activity_coefs = ACTIVITY_COEFFICIENTS_ITERATION    
          case('NEWTON')
            reaction%compute_activity_coefs = ACTIVITY_COEFFICIENTS_NEWTON       
          case default
            reaction%compute_activity_coefs = ACTIVITY_COEFFICIENTS_TIMESTEP        
        end select
      case default
        call printErrMsg(option,'CHEMISTRY keyword: '//trim(word)//' not recognized')
    end select
  enddo
  
  reaction%nsorb = reaction%neqsurfcmplxrxn + reaction%neqionxrxn

  if (reaction%neqcmplx + reaction%nsorb + reaction%nmnrl > 0) then
    reaction%use_full_geochemistry = PETSC_TRUE
  endif

  
  if (len_trim(reaction%database_filename) < 2) &
    reaction%compute_activity_coefs = ACTIVITY_COEFFICIENTS_OFF
 
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
  use Fileio_module
  use Utility_module  
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: found
  PetscInt :: icomp, jcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: igas, dummy_int
  PetscReal :: value
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
      if (fiStringCompare(aq_species_constraint%names(icomp), &
                          reaction%primary_species_names(jcomp), &
                          MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      string = 'Species ' // trim(aq_species_constraint%names(icomp)) // &
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among primary species.'
      call printErrMsg(option,string)
    else
      constraint_type(jcomp) = aq_species_constraint%constraint_type(icomp)
      constraint_spec_name(jcomp) = aq_species_constraint%constraint_spec_name(icomp)
      constraint_conc(jcomp) = aq_species_constraint%constraint_conc(icomp)
      
      ! link constraint species
      select case(constraint_type(jcomp))
        case(CONSTRAINT_MINERAL)
          found = PETSC_FALSE
          do imnrl = 1, reaction%nmnrl
            if (fiStringCompare(constraint_spec_name(jcomp), &
                                reaction%mineral_names(imnrl), &
                                MAXWORDLENGTH)) then
              constraint_id(jcomp) = imnrl
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            string = 'Constraint mineral: ' // &
                     trim(constraint_spec_name(jcomp)) // &
                     ' for aqueous species: ' // &
                     trim(reaction%primary_species_names(jcomp)) // &
                     ' in constraint: ' // &
                     trim(constraint_name) // ' not found.' 
            call printErrMsg(option,string)         
          endif
        case(CONSTRAINT_GAS)
          found = PETSC_FALSE
          do igas = 1, reaction%ngas
            if (fiStringCompare(constraint_spec_name(jcomp), &
                                reaction%gas_species_names(igas), &
                                MAXWORDLENGTH)) then
              constraint_id(jcomp) = igas
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            string = 'Constraint gas: ' // &
                     trim(constraint_spec_name(jcomp)) // &
                     ' for aqueous species: ' // &
                     trim(reaction%primary_species_names(jcomp)) // &
                     ' in constraint: ' // &
                     trim(constraint_name) // ' not found.' 
            call printErrMsg(option,string)         
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
        if (fiStringCompare(mineral_constraint%names(imnrl), &
                            reaction%kinmnrl_names(jmnrl), &
                            MAXWORDLENGTH)) then
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        string = 'Mineral ' // trim(mineral_constraint%names(imnrl)) // &
                 'from CONSTRAINT ' // trim(constraint_name) // &
                 ' not found among kinetic minerals.'
        call printErrMsg(option,string)
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
                                         num_iterations,option)
  use Option_module
  use Fileio_module
  use Utility_module  
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  PetscInt :: num_iterations
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
  PetscReal :: ln_act_h2o
  PetscReal :: lnQK, QK
  PetscReal :: tempreal
  PetscInt :: comp_id
  PetscReal :: convert_molal_to_molar
  PetscReal :: convert_molar_to_molal
  
  PetscTruth :: charge_balance_warning_flag = PETSC_FALSE

  PetscReal :: Jac_num(reaction%ncomp)
  PetscReal :: Res_pert, pert, prev_value

  PetscInt :: iphase
    
  constraint_type = aq_species_constraint%constraint_type
  constraint_spec_name = aq_species_constraint%constraint_spec_name
  constraint_id = aq_species_constraint%constraint_spec_id
  conc = aq_species_constraint%constraint_conc

  iphase = 1
  
  if (.not.reaction%use_full_geochemistry) then
    aq_species_constraint%basis_molarity = conc
    return
  endif
  
  if (option%initialize_with_molality) then
    convert_molal_to_molar = global_auxvar%den_kg(iphase)/1000.d0
    convert_molar_to_molal = 1.d0
  else
    convert_molal_to_molar = 1.d0
    convert_molar_to_molal = 1000.d0/global_auxvar%den_kg(iphase)
  endif
  
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
            if (.not.fiStringCompare(reaction%primary_species_names(icomp), &
                                     string,MAXWORDLENGTH)) then
              string = 'pH specified as constraint (constraint =' // &
                       trim(constraint_name) // &
                       ') for species other than H+ or OH-: ' // &
                       trim(reaction%primary_species_names(icomp))
              call printErrMsg(option,string)
            endif
          endif
          free_conc(icomp) = 10.d0**(-conc(icomp))
        else
          string = 'pH specified as constraint (constraint =' // &
                   trim(constraint_name) // &
                   '), but H+ not found in chemical species.'
          call printErrMsg(option,string)
        endif        
      case(CONSTRAINT_MINERAL)
        free_conc(icomp) = conc(icomp)*convert_molar_to_molal ! guess
      case(CONSTRAINT_GAS)
        if (conc(icomp) <= 0.d0) then ! log form
          conc(icomp) = 10.d0**conc(icomp) ! conc log10 partial pressure gas
        endif
        free_conc(icomp) = 1.d-9 ! guess
    end select
  enddo
  
  rt_auxvar%pri_molal = free_conc

  num_iterations = 0
  compute_activity_coefs = PETSC_FALSE
  
  do

    if (reaction%compute_activity_coefs /= ACTIVITY_COEFFICIENTS_OFF .and. &
        compute_activity_coefs) then
      call RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
    endif
    call RTotal(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%nsorb > 0) call RTotalSorb(rt_auxvar,global_auxvar,reaction,option)
    
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
          Res(icomp) = rt_auxvar%total(icomp,1) + rt_auxvar%total_sorb(icomp)/tempreal - &
                       total_conc(icomp)
          ! dtotal units = kg water/L water
          ! dtotal_sorb units = kg water/m^3 bulk
          ! Jac units = kg water/L water
          Jac(icomp,:) = rt_auxvar%dtotal(icomp,:,1) + &
          ! dtotal_sorb units = kg water/m^3 bulk
                         rt_auxvar%dtotal_sorb(icomp,:)/tempreal

        case(CONSTRAINT_FREE,CONSTRAINT_LOG)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
!          Jac(:,icomp) = 0.d0
          Jac(icomp,icomp) = 1.d0
          
        case(CONSTRAINT_CHARGE_BAL)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
          do jcomp = 1, reaction%ncomp
            Res(icomp) = Res(icomp) + reaction%primary_spec_Z(jcomp)*rt_auxvar%total(jcomp,1)
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
              string = 'Charge balance species ' // &
                       trim(reaction%primary_species_names(icomp)) // &
                       ' may not satify constraint ' // &
                       trim(constraint_name) // &
                       '.  Molality already below 1.e-20.'
              call printMsg(option,string)
              charge_balance_warning_flag = PETSC_TRUE
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
          
            ln_act_h2o = 0.d0
          
            icplx = abs(reaction%h_ion_id)
            
            ! compute secondary species concentration
            ! *note that the sign was flipped below
            lnQK = -reaction%eqcmplx_logK(icplx)*LOG_TO_LN

            ! activity of water
            if (reaction%eqcmplxh2oid(icplx) > 0) then
              lnQK = lnQK + reaction%eqcmplxh2ostoich(icplx)*ln_act_h2o
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

          ln_act_h2o = 0.d0
  
          imnrl = constraint_id(icomp)
          ! compute secondary species concentration
          lnQK = -reaction%mnrl_logK(imnrl)*LOG_TO_LN

          ! activity of water
          if (reaction%mnrlh2oid(imnrl) > 0) then
            lnQK = lnQK + reaction%mnrlh2ostoich(imnrl)*ln_act_h2o
          endif

          do jcomp = 1, reaction%mnrlspecid(0,imnrl)
            comp_id = reaction%mnrlspecid(jcomp,imnrl)
            lnQK = lnQK + reaction%mnrlstoich(jcomp,imnrl)* &
                          log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
          enddo
!         QK = exp(lnQK)
          
!         Res(icomp) = 1.d0 - QK
          Res(icomp) = lnQK

          do jcomp = 1,reaction%mnrlspecid(0,imnrl)
            comp_id = reaction%mnrlspecid(jcomp,imnrl)
!           Jac(icomp,comp_id) = -QK/auxvar%primary_spec(comp_id)* &
!                                reaction%mnrlstoich(jcomp,imnrl)
            Jac(icomp,comp_id) = reaction%mnrlstoich(jcomp,imnrl)/rt_auxvar%pri_molal(comp_id)
                                 
          enddo
  
        case(CONSTRAINT_GAS)

          ln_act_h2o = 0.d0
  
          igas = constraint_id(icomp)
          
          ! compute secondary species concentration
          lnQK = -reaction%eqgas_logK(igas)*LOG_TO_LN
          
          ! divide K by RT
          !lnQK = lnQK - log((auxvar%temp+273.15d0)*IDEAL_GAS_CONST)
          
          ! activity of water
          if (reaction%eqgash2oid(igas) > 0) then
            lnQK = lnQK + reaction%eqgash2ostoich(igas)*ln_act_h2o
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
            Jac(icomp,comp_id) = reaction%eqgasstoich(jcomp,igas)/rt_auxvar%pri_molal(comp_id)
          enddo
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
100 format('Constraint iteration count has exceeded: ',i5)
      write(string,100) num_iterations
      call printMsg(option,string)
      if (num_iterations >= 10000) then
        string = 'Stopping due to excessive iteration count!'
        call printErrMsg(option,string)
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
  if (reaction%nsorb > 0) call RTotalSorb(rt_auxvar,global_auxvar,reaction,option)
  
  ! remember that a density of 1 kg/L was assumed, thus molal and molarity are equal
  ! do not scale by molal_to_molar since it could be 1.d0 if MOLAL flag set
  aq_species_constraint%basis_molarity = rt_auxvar%pri_molal* &
                                         global_auxvar%den_kg(option%liquid_phase)/ &
                                         1000.d0
  
  if (option%myrank == 0) &
    print *,'ReactionEquilibrateConstraint: ' // trim(constraint_name) // &
            '  iterations: ',num_iterations

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
  use Fileio_module
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
  PetscReal :: ln_act_h2o
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
  
  global_auxvar%den_kg(iphase) = option%reference_density
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
    
    write(option%fid_out,204) '  ionic strength: ', ionic_strength
    write(option%fid_out,204) '  charge balance: ', charge_balance
    
    write(option%fid_out,202) '        pressure: ', global_auxvar%pres(1)
    write(option%fid_out,203) '     temperature: ', global_auxvar%temp(1)
    write(option%fid_out,90)

    102 format(/,'  species               molality    total       act coef  constraint')  
    write(option%fid_out,102)
    write(option%fid_out,90)
  
    103 format(2x,a20,es12.4,es12.4,f8.4,4x,a)
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
      end select
      write(option%fid_out,103) reaction%primary_species_names(icomp), &
                                rt_auxvar%pri_molal(icomp), &
                                rt_auxvar%total(icomp,1), &
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
    111 format(2x,a20,es12.4,f8.4,2x,es12.4)
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
  if (reaction%nsorb > 0) then
    write(option%fid_out,1128)
    write(option%fid_out,90)
    do jcomp = 1, reaction%ncomp
      if (abs(rt_auxvar%total(jcomp,iphase)) > 0.d0) &
      retardation = 1.d0 + rt_auxvar%total_sorb(jcomp)/bulk_vol_to_fluid_vol &
        /rt_auxvar%total(jcomp,iphase)
      totj = rt_auxvar%total(jcomp,iphase)+rt_auxvar%total_sorb(jcomp)/bulk_vol_to_fluid_vol
      write(option%fid_out,129) reaction%primary_species_names(jcomp), &
        totj,retardation
    enddo
   1128 format(/,2x,'primary species     total(aq+sorbed)    total retardation', &
               /,25x,'[mol/L]',15x,'1+Kd')
    129 format(2x,a12,8x,1pe12.4,8x,1pe12.4)
  endif
  
  if (reaction%nmnrl > 0) then
  
    130 format(/,'  mineral                   log SI    log K')
    131 format(2x,a20,2x,f10.4,2x,1pe12.4)

    do imnrl = 1, reaction%nmnrl
      ! compute saturation
      ln_act_h2o = 0.d0
      lnQK(imnrl) = -reaction%mnrl_logK(imnrl)*LOG_TO_LN
      if (reaction%mnrlh2oid(imnrl) > 0) then
        lnQK(imnrl) = lnQK(imnrl) + reaction%mnrlh2ostoich(imnrl)*ln_act_h2o
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

      ln_act_h2o = 0.d0
      
      ! compute gas partial pressure
      lnQK(igas) = -reaction%eqgas_logK(igas)*LOG_TO_LN
      
      ! divide K by RT
      !lnQK = lnQK - log((auxvar%temp+273.15d0)*IDEAL_GAS_CONST)
      
      ! activity of water
      if (reaction%eqgash2oid(igas) > 0) then
        lnQK(igas) = lnQK(igas) + reaction%eqgash2ostoich(igas)*ln_act_h2o
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
subroutine ReactionReadMineralKinetics(reaction,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  
  type(mineral_type), pointer :: cur_mineral
  PetscInt :: imnrl
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

    if (fiCheckExit(string)) exit  

    call fiReadWord(string,name,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,MINERAL_KINETICS', ierr)
    
    cur_mineral => reaction%mineral_list
    do 
      if (.not.associated(cur_mineral)) exit
      if (fiStringCompare(cur_mineral%name,name,MAXWORDLENGTH)) then
        cur_mineral%itype = MINERAL_KINETIC
        if (.not.associated(cur_mineral%tstrxn)) then
          cur_mineral%tstrxn => TransitionStateTheoryRxnCreate()
        endif
        ! read rate
        call fiReadDouble(string,cur_mineral%tstrxn%rate,ierr)
        call fiErrorMsg(option%myrank,'rate','CHEMISTRY,MINERAL_KINETICS', ierr)
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
      string = 'No rate provided in input file for mineral: ' // &
               trim(cur_mineral%name) // '.'
      call printErrMsg(option,string)
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
   
  if (reaction%nkinmnrl > 0) then
    call RKineticMineral(Res,Jac,derivative,rt_auxvar,global_auxvar,volume, &
                         reaction,option)
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
  type(reactive_transport_auxvar_type) :: rt_auxvar_pert
  PetscTruth :: compute_derivative

  ! add new reactions in the 3 locations below

  if (.not.option%numerical_derivatives) then ! analytical derivative
  !if (PETSC_FALSE) then
    compute_derivative = PETSC_TRUE
    ! #1: add new reactions here
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res,Jac,compute_derivative,rt_auxvar, &
                           global_auxvar,volume,reaction,option)
    endif
  else ! numerical derivative
    compute_derivative = PETSC_FALSE
    Res_orig = 0.d0
    call RTAuxVarInit(rt_auxvar_pert,reaction,option)
    call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
    ! #2: add new reactions here
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                           global_auxvar,volume,reaction,option)
    endif
    do jcomp = 1, reaction%ncomp
      Res_pert = 0.d0
      call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
      pert = rt_auxvar_pert%pri_molal(jcomp)*perturbation_tolerance
      rt_auxvar_pert%pri_molal(jcomp) = rt_auxvar_pert%pri_molal(jcomp) + pert
      
      call RTotal(rt_auxvar_pert,global_auxvar,reaction,option)
      if (reaction%nsorb > 0) call RTotalSorb(rt_auxvar_pert,global_auxvar, &
                                              reaction,option)

      ! #3: add new reactions here
      if (reaction%nkinmnrl > 0) then
        call RKineticMineral(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                             global_auxvar,volume,reaction,option)
      endif
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
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: icplx, icomp, it, j, jcomp, ncomp
  PetscReal :: I, sqrt_I, II, sqrt_II, f, fpri, didi, dcdi, den, dgamdi, lnQK, sum
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: ln_act_h2o


  if (reaction%compute_activity_coefs == ACTIVITY_COEFFICIENTS_NEWTON) then

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  ln_act_h2o = 0.d0  ! assume act h2o = 1 for now
  
  ! compute primary species contribution to ionic strength
  fpri = 0.d0
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
          sum = 0.5d0*reaction%debyeA*reaction%eqcmplx_Z(icplx)*reaction%eqcmplx_Z(icplx) &
            /(sqrt_I*(1.d0+reaction%debyeB*reaction%eqcmplx_a0(icplx)*sqrt_I)**2) &
            -reaction%debyeBdot
          ncomp = reaction%eqcmplxspecid(0,icplx)
          do jcomp = 1, ncomp
            j = reaction%eqcmplxspecid(jcomp,icplx)
            if(abs(reaction%primary_spec_Z(j)) > 0.d0) then
               dgamdi = -0.5d0*reaction%debyeA*reaction%primary_spec_Z(j)**2/(sqrt_I* &
                 (1.d0+reaction%debyeB*reaction%primary_spec_a0(j)*sqrt_I)**2)+reaction%debyeBdot 
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
      print *,'ionic strength negative! it =',it,' I= ',I,II,den,didi,dcdi,sum
      stop
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
        lnQK = lnQK + reaction%eqcmplxh2ostoich(icplx)*ln_act_h2o
      endif

      ncomp = reaction%eqcmplxspecid(0,icplx)
      do jcomp = 1, ncomp
        icomp = reaction%eqcmplxspecid(jcomp,icplx)
        lnQK = lnQK + reaction%eqcmplxstoich(jcomp,icplx)*ln_act(icomp)
      enddo
      rt_auxvar%sec_molal(icplx) = exp(lnQK)/rt_auxvar%sec_act_coef(icplx)
    enddo
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
  enddo
  
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
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, icplx, icomp, jcomp, iphase, ncomp
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: ln_act_h2o
  PetscReal :: lnQK, tempreal
  PetscReal :: den_kg_per_L

  iphase = 1           
  
  den_kg_per_L = global_auxvar%den_kg(iphase)*1.d-3              

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  ln_act_h2o = 0.d0  ! assume act h2o = 1 for now
  rt_auxvar%total(:,iphase) = rt_auxvar%pri_molal(:)
  ! initialize derivatives
  rt_auxvar%dtotal = 0.d0
  do icomp = 1, reaction%ncomp
    rt_auxvar%dtotal(icomp,icomp,iphase) = 1.d0
  enddo
  
  do icplx = 1, reaction%neqcmplx ! for each secondary species
    ! compute secondary species concentration
    lnQK = -reaction%eqcmplx_logK(icplx)*LOG_TO_LN

    ! activity of water
    if (reaction%eqcmplxh2oid(icplx) > 0) then
      lnQK = lnQK + reaction%eqcmplxh2ostoich(icplx)*ln_act_h2o
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
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, iphase, ncomp, ncplx
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: surfcmplx_conc(reaction%neqsurfcmplx)
  PetscReal :: dSx_dmi(reaction%ncomp)
  PetscReal :: dSi_dSx
  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: ln_act_h2o
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscReal, parameter :: tol = 1.d-12
  PetscTruth :: one_more
  PetscReal :: res, dres_dfree_site, dfree_site_conc
  
  PetscReal :: omega
  PetscReal :: ref_cation_X, ref_cation_conc, ref_cation_Z, ref_cation_k, &
               ref_cation_quotient
  PetscReal :: cation_X(reaction%ncomp)
  PetscReal :: dres_dref_cation_X, dref_cation_X
  PetscReal :: sumZX, sumkm
  
  PetscReal :: total_pert, ref_cation_X_pert, pert
  PetscReal :: ref_cation_quotient_pert, dres_dref_cation_X_pert

  iphase = 1                         

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  ln_act_h2o = 0.d0  ! assume act h2o = 1 for now
    
  ! initialize total sorbed concentrations and derivatives
  rt_auxvar%total_sorb = 0.d0
  rt_auxvar%dtotal_sorb = 0.d0

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
          lnQK = lnQK + reaction%eqsurfcmplxh2ostoich(icplx)*ln_act_h2o
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
        rt_auxvar%total_sorb(icomp) = rt_auxvar%total_sorb(icomp) + &
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
          rt_auxvar%dtotal_sorb(icomp,jcomp) = rt_auxvar%dtotal_sorb(icomp,jcomp) + &
                                               reaction%eqsurfcmplxstoich(i,icplx)* &
                                               tempreal
        enddo
      enddo
    enddo
  enddo
  
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

      rt_auxvar%total_sorb(icomp) = rt_auxvar%total_sorb(icomp) + tempreal1

      tempreal2 = reaction%primary_spec_Z(icomp)/sumZX
      do j = 1, ncomp
        jcomp = reaction%eqionx_rxn_cationid(j,irxn)
        if (i == j) then
          rt_auxvar%dtotal_sorb(icomp,jcomp) = rt_auxvar%dtotal_sorb(icomp,jcomp) + &
                                               tempreal1*(1.d0-(tempreal2*cation_X(j)))/ &
                                               rt_auxvar%pri_molal(jcomp)
        else
          rt_auxvar%dtotal_sorb(icomp,jcomp) = rt_auxvar%dtotal_sorb(icomp,jcomp) + &
                                               (-tempreal1)*tempreal2*cation_X(j)/ &
                                               rt_auxvar%pri_molal(jcomp)
        endif
      enddo
    enddo    

  enddo
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
end subroutine RTotalSorb

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
  PetscReal :: ln_act_h2o
  PetscReal :: QK, lnQK, dQK_dCj, dQK_dmj
  PetscTruth :: prefactor_exists

  iphase = 1                         

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  if (reaction%neqcmplx > 0) then
    ln_sec = log(rt_auxvar%sec_molal)
    ln_sec_act = ln_sec+log(rt_auxvar%sec_act_coef)
  endif
  
  ln_act_h2o = 0.d0
  
  do imnrl = 1, reaction%nkinmnrl ! for each mineral
    ! compute secondary species concentration
    lnQK = -reaction%kinmnrl_logK(imnrl)*LOG_TO_LN

    ! activity of water
    if (reaction%kinmnrlh2oid(imnrl) > 0) then
      lnQK = lnQK + reaction%kinmnrlh2ostoich(imnrl)*ln_act_h2o
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

end module Reaction_module
