module Reaction_module

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: ReactionCreate, &
            ReactionRead, &
            ReactionReadMineralKinetics, &
!           ReactionReadSurfaceComplexes, &
            ReactionInitializeConstraint, &
            RTotal, &
            RTotalSorb, &
            RActivity, &
            RReaction, &
            RReactionDerivative, &
            RPrintConstraint

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
          !   call SurfaceComplexRXNRead
          
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
          !   call IonExchangeRXNRead
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
        option%use_log_formulation = PETSC_TRUE        
      case('ACTIVITY')
        reaction%compute_activity = PETSC_TRUE        
      case default
        call printErrMsg(option,'CHEMISTRY keyword: '//trim(word)//' not recognized')
    end select
  enddo
  
  reaction%nsorb = reaction%neqsurfcmplxrxn + reaction%neqionxrxn
  
  if (len_trim(reaction%database_filename) < 2) &
    reaction%compute_activity = PETSC_FALSE
 
end subroutine ReactionRead

! ************************************************************************** !
!
! ReactionInitializeConstraint: Initializes constraints based on primary
!                               species in system
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine ReactionInitializeConstraint(reaction,constraint_name, &
                                        aq_species_constraint, &
                                        mineral_constraint,option)
  use Option_module
  use Fileio_module
  use Utility_module  
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  character(len=MAXNAMELENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(option_type) :: option
  
  type(reactive_transport_auxvar_type) :: auxvar
  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: found
  PetscInt :: icomp, jcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: igas
  PetscReal :: value
  PetscReal :: constraint_conc(reaction%ncomp)
  PetscInt :: constraint_type(reaction%ncomp)
  character(len=MAXNAMELENGTH) :: constraint_spec_name(reaction%ncomp)
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
                          MAXNAMELENGTH)) then
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
                                MAXNAMELENGTH)) then
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
                                MAXNAMELENGTH)) then
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
  if (associated(reaction)) then
    if (associated(mineral_constraint)) then
      do imnrl = 1, reaction%nmnrl
        found = PETSC_FALSE
        do jmnrl = 1, reaction%nmnrl
          if (fiStringCompare(mineral_constraint%names(imnrl), &
                              reaction%mineral_names(jmnrl), &
                              MAXNAMELENGTH)) then
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (.not.found) then
          string = 'Mineral ' // trim(mineral_constraint%names(imnrl)) // &
                   'from CONSTRAINT ' // trim(constraint_name) // &
                   ' not found among primary species.'
          call printErrMsg(option,string)
        else
          mineral_constraint%basis_mol_frac(jmnrl) = &
            mineral_constraint%constraint_mol_frac(imnrl)
        endif  
      enddo
    endif
  endif
  
  call RTAuxVarInit(auxvar,option)
  call ReactionEquilibrateConstraint(auxvar,reaction,constraint_name, &
                                     aq_species_constraint, &
                                     option)
  call RTAuxVarDestroy(auxvar)

end subroutine ReactionInitializeConstraint

! ************************************************************************** !
!
! ReactionEquilibrateConstraint: Equilibrates constraint concentrations
!                                with prescribed geochemistry
! author: Glenn Hammond
! date: 10/22/08
!
! ************************************************************************** !
subroutine ReactionEquilibrateConstraint(auxvar,reaction,constraint_name, &
                                         aq_species_constraint, &
                                         option)
  use Option_module
  use Fileio_module
  use Utility_module  
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type), pointer :: reaction
  character(len=MAXNAMELENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icomp, jcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: igas
  PetscReal :: conc(reaction%ncomp)
  PetscInt :: constraint_type(reaction%ncomp)
  character(len=MAXNAMELENGTH) :: constraint_spec_name(reaction%ncomp)

  PetscReal :: Res(reaction%ncomp)
  PetscReal :: total_conc(reaction%ncomp)
  PetscReal :: free_conc(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscInt :: indices(reaction%ncomp)
  PetscInt :: num_it
  PetscReal :: norm
  PetscReal :: prev_molal(reaction%ncomp)
  PetscReal, parameter :: tol = 1.d-12
  PetscReal, parameter :: tol_loose = 1.d-6
  PetscTruth :: compute_activity

  PetscInt :: constraint_id(reaction%ncomp)
  PetscReal :: ln_act_h2o
  PetscReal :: lnQK, QK
  PetscInt :: comp_id
  PetscReal, parameter :: log_to_ln = 2.30258509299d0
    
  if (.not.associated(reaction)) return
  
  constraint_type = aq_species_constraint%constraint_type
  constraint_spec_name = aq_species_constraint%constraint_spec_name
  constraint_id = aq_species_constraint%constraint_spec_id
  conc = aq_species_constraint%constraint_conc
  
  total_conc = 0.d0
  do icomp = 1, reaction%ncomp
    select case(constraint_type(icomp))
      case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
        total_conc(icomp) = conc(icomp)
        free_conc(icomp) = 1.d-9
      case(CONSTRAINT_FREE)
        free_conc(icomp) = conc(icomp)
      case(CONSTRAINT_LOG)
        free_conc(icomp) = 10**conc(icomp)
      case(CONSTRAINT_MINERAL)
        free_conc(icomp) = conc(icomp) ! guess
      case(CONSTRAINT_GAS)
        free_conc(icomp) = conc(icomp) ! guess
    end select
  enddo
  
  auxvar%den(1) = 1000.d0 ! assume a density of 1 kg/L (1000 kg/m^3)
  auxvar%primary_molal = free_conc

  num_it = 0
  compute_activity = PETSC_FALSE
  
  do

    auxvar%primary_spec = auxvar%primary_molal ! assume a density of 1 kg/L
    if (reaction%compute_activity .and. compute_activity) then
      call RActivity(auxvar,reaction,option)
    endif
    call RTotal(auxvar,reaction,option)
    if (reaction%nsorb > 0) call RTotalSorb(auxvar,reaction,option)
    
    Res = auxvar%total(:,1)
    Jac = auxvar%dtotal(:,:,1)
    if (reaction%neqsurfcmplxrxn > 0) Jac = Jac + auxvar%dtotal_sorb(:,:)
    ! dtotal must be scaled by 1.d-3 to scale density in RTotal from kg/m^3 -> kg/L
    Jac = Jac * 1.d-3
        
    do icomp = 1, reaction%ncomp
      select case(constraint_type(icomp))
        case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
        case(CONSTRAINT_FREE,CONSTRAINT_LOG)
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
          Jac(:,icomp) = 0.d0
          Jac(icomp,icomp) = 1.d0
        case(CONSTRAINT_MINERAL)

          ln_act_h2o = 0.d0
  
          imnrl = constraint_id(icomp)
          ! compute secondary species concentration
          lnQK = -1.d0*reaction%kinmnrl_logK(imnrl)*log_to_ln

          ! activity of water
          if (reaction%kinmnrlh2oid(imnrl) > 0) then
            lnQK = lnQK + reaction%kinmnrlh2ostoich(imnrl)*ln_act_h2o
          endif

          do jcomp = 1, reaction%kinmnrlspecid(0,imnrl)
            comp_id = reaction%kinmnrlspecid(jcomp,imnrl)
            lnQK = lnQK + reaction%kinmnrlstoich(jcomp,imnrl)* &
                          log(auxvar%primary_spec(comp_id)*auxvar%pri_act_coef(comp_id))
          enddo
          QK = exp(lnQK)
          
          Res(icomp) = 1.d0 - QK
          Jac(icomp,:) = 0.d0
          do jcomp = 1,reaction%kinmnrlspecid(0,imnrl)
            comp_id = reaction%kinmnrlspecid(jcomp,imnrl)
            Jac(icomp,comp_id) = -exp(lnQK-log(auxvar%primary_molal(comp_id)))* &
                                 reaction%kinmnrlstoich(jcomp,imnrl)
          enddo
        case(CONSTRAINT_GAS)

          ln_act_h2o = 0.d0
  
          igas = constraint_id(icomp)
          ! compute secondary species concentration
          lnQK = -1.d0*reaction%eqgas_logK(igas)*log_to_ln

          ! activity of water
          if (reaction%eqgash2oid(igas) > 0) then
            lnQK = lnQK + reaction%eqgash2ostoich(igas)*ln_act_h2o
          endif

          do jcomp = 1, reaction%eqgasspecid(0,igas)
            comp_id = reaction%eqgasspecid(jcomp,igas)
            lnQK = lnQK + reaction%eqgasstoich(jcomp,igas)* &
                          log(auxvar%primary_spec(comp_id)*auxvar%pri_act_coef(comp_id))
          enddo
          QK = exp(lnQK)
          
          Res(icomp) = 1.d0 - QK
          Jac(icomp,:) = 0.d0
          do jcomp = 1,reaction%eqgasspecid(0,igas)
            comp_id = reaction%eqgasspecid(jcomp,igas)
            Jac(icomp,comp_id) = -exp(lnQK-log(auxvar%primary_molal(comp_id)))* &
                                 reaction%eqgasstoich(jcomp,igas)
          enddo
      end select
    enddo
    
    Res = total_conc - Res

    ! scale Jacobian
    do icomp = 1, reaction%ncomp
      norm = max(1.d0,maxval(abs(Jac(icomp,:))))
      norm = 1.d0/norm
      Res(icomp) = Res(icomp)*norm
      Jac(icomp,:) = Jac(icomp,:)*norm
    enddo
      
    ! for derivatives with respect to ln conc
    do icomp = 1, reaction%ncomp
      Jac(:,icomp) = Jac(:,icomp)*auxvar%primary_spec(icomp)
    enddo

    call ludcmp(Jac,reaction%ncomp,indices,icomp)
    call lubksb(Jac,reaction%ncomp,indices,Res)

    prev_molal = auxvar%primary_molal

    Res = dsign(1.d0,Res)*min(dabs(Res),5.d0)
      
    auxvar%primary_molal = auxvar%primary_molal*exp(Res)

    num_it = num_it + 1
    
    ! need some sort of convergence before we kick in activities
    if (maxval(dabs(auxvar%primary_molal-prev_molal)/ &
               auxvar%primary_molal) < tol_loose) then
      compute_activity = PETSC_TRUE
    endif

    ! check for convergence
    if (maxval(dabs(auxvar%primary_molal-prev_molal)/ &
               auxvar%primary_molal) < tol) exit
                     
  enddo
  
  ! remember that a density of 1 kg/L was assumed, thus molal and molarity are equal
  aq_species_constraint%basis_molarity = auxvar%primary_molal

end subroutine ReactionEquilibrateConstraint

! ************************************************************************** !
!
! RPrintConstraint: Prints a constraint associated with reactive transport
! author: Glenn Hammond
! date: 10/28/08
!
! ************************************************************************** !
subroutine RPrintConstraint(constraint_coupler,reaction,option)

  use Option_module
  use Condition_module

  implicit none
  
  type(option_type) :: option
  type(tran_constraint_coupler_type) :: constraint_coupler
  type(reaction_type), pointer :: reaction
  
  type(reactive_transport_auxvar_type) :: auxvar
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, icomp
  PetscInt :: icplx, icplx2
  PetscInt :: eqcmplxsort(reaction%neqcmplx)
  PetscInt :: eqsurfcmplxsort(reaction%neqsurfcmplx+reaction%neqsurfcmplxrxn)
  PetscTruth :: finished
  PetscReal :: conc, conc2

  aq_species_constraint => constraint_coupler%aqueous_species

99 format(80('-'))
100 format(a)

  write(option%fid_out,100) '  Constraint: ' // trim(constraint_coupler%constraint_name)
101 format(/,'  species       molality    total       act coef  constraint')  
  write(option%fid_out,101)
  write(option%fid_out,99)
  
 call RTAuxVarInit(auxvar,option)
  
102 format(2x,a12,es12.4,es12.4,f8.4,4x,a)

  call ReactionEquilibrateConstraint(auxvar,reaction, &
                                     constraint_coupler%constraint_name, &
                                     aq_species_constraint, &
                                     option)
                                     
  do icomp = 1, reaction%ncomp
    select case(aq_species_constraint%constraint_type(icomp))
      case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
        string = 'total'
      case(CONSTRAINT_FREE)
        string = 'free'
      case(CONSTRAINT_LOG)
        string = 'log'
      case(CONSTRAINT_MINERAL,CONSTRAINT_GAS)
        string = aq_species_constraint%constraint_spec_name(icomp)
    end select
    write(option%fid_out,102) reaction%primary_species_names(icomp), &
                              auxvar%primary_molal(icomp), &
                              auxvar%total(icomp,1)/auxvar%den(1)*1000.d0, &
                              auxvar%pri_act_coef(icomp), &
                              trim(string)
  enddo  
      
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
        if (auxvar%secondary_spec(icplx) < &
            auxvar%secondary_spec(icplx2)) then
          eqcmplxsort(i) = icplx2
          eqcmplxsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo
            
  103 format(/,'  complex       molality    act coef  logK')  
    write(option%fid_out,103)
    write(option%fid_out,99)
  104 format(2x,a12,es12.4,f8.4,2x,es12.4)
    do i = 1, reaction%neqcmplx ! for each secondary species
      icplx = eqcmplxsort(i)
      write(option%fid_out,104) reaction%secondary_species_names(icplx), &
                                auxvar%secondary_spec(icplx)/ &
                                auxvar%den(1)*1000.d0, &
                                auxvar%sec_act_coef(icplx), &
                                reaction%eqcmplx_logK(icplx)
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
          conc = auxvar%eqsurfcmplx_spec(icplx)
        else
          conc = auxvar%eqsurfcmplx_freesite_conc(-icplx)
        endif
        if (icplx2 > 0) then
          conc2 = auxvar%eqsurfcmplx_spec(icplx2)
        else
          conc2 = auxvar%eqsurfcmplx_freesite_conc(-icplx2)
        endif
        if (conc < conc2) then
          eqsurfcmplxsort(i) = icplx2
          eqsurfcmplxsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo
            
  105 format(/,'  surf complex  molality    logK')  
    write(option%fid_out,105)
    write(option%fid_out,99)
  106 format(2x,a12,es12.4,es12.4)
  107 format(2x,a12,es12.4,'  free site')
    do i = 1, reaction%neqsurfcmplx+reaction%neqsurfcmplxrxn
      icplx = eqsurfcmplxsort(i)
      if (icplx > 0) then
        write(option%fid_out,106) reaction%surface_complex_names(icplx), &
                                  auxvar%eqsurfcmplx_spec(icplx)/ &
                                  auxvar%den(1)*1000.d0, &
                                  reaction%eqsurfcmplx_logK(icplx)
      else
        write(option%fid_out,107) reaction%surface_site_names(-icplx), &
                                  auxvar%eqsurfcmplx_freesite_conc(-icplx)/ &
                                  auxvar%den(1)*1000.d0
      endif
    enddo 
  endif
          
  call RTAuxVarDestroy(auxvar)
            
end subroutine RPrintConstraint

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
  character(len=MAXNAMELENGTH) :: name
  
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

    call fiReadNChars(string,name,MAXNAMELENGTH,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY', ierr)
    
    cur_mineral => reaction%mineral_list
    do 
      if (.not.associated(cur_mineral)) exit
      if (fiStringCompare(cur_mineral%name,name,MAXNAMELENGTH)) then
        if (.not.associated(cur_mineral%tstrxn)) then
          cur_mineral%tstrxn => TransitionStateTheoryRxnCreate()
        endif
        ! read rate
        call fiReadDouble(string,cur_mineral%tstrxn%rate,ierr)
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
    if (cur_mineral%id < 0) then
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
  
end subroutine ReactionReadMineralKinetics

! ************************************************************************** !
!
! ReactionReadSurfaceComplexes: Reads surface complexation reactions
! author: Peter Lichtner
! date: 10/16/08
!
! ************************************************************************** !
#if 0
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

    if (fiCheckExit(string)) exit  

    call fiReadNChars(string,name,MAXNAMELENGTH,PETSC_TRUE,ierr)
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
#endif

! ************************************************************************** !
!
! RReaction: Computes reactions
! author: Glenn Hammond
! date: 09/30/08
!
! ************************************************************************** !
subroutine RReaction(Res,Jac,derivative,auxvar,volume,reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: auxvar 
  type(option_type) :: option
  PetscTruth :: derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: volume
   
  if (reaction%nkinmnrl > 0) then
    call RKineticMineral(Res,Jac,derivative,auxvar,volume,reaction,option)
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
subroutine RReactionDerivative(Res,Jac,auxvar,volume,reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: auxvar 
  type(option_type) :: option
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: volume
   
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp
  PetscReal :: Jac_dummy(reaction%ncomp,reaction%ncomp)
  PetscReal :: pert
  type(reactive_transport_auxvar_type) :: auxvar_pert
  PetscTruth :: compute_derivative

  ! add new reactions in the 3 locations below

  if (.not.option%numerical_derivatives) then ! analytical derivative
!  if (PETSC_FALSE) then
    compute_derivative = PETSC_TRUE
    ! #1: add new reactions here
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res,Jac,compute_derivative,auxvar,volume,reaction,option)
    endif
  else ! numerical derivative
    compute_derivative = PETSC_FALSE
    Res_orig = 0.d0
    call RTAuxVarInit(auxvar_pert,option)
    call RTAuxVarCopy(auxvar_pert,auxvar,option)
    ! #2: add new reactions here
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res_orig,Jac_dummy,compute_derivative,auxvar,volume,reaction,option)
    endif
    do jcomp = 1, reaction%ncomp
      Res_pert = 0.d0
      call RTAuxVarCopy(auxvar_pert,auxvar,option)
      pert = auxvar_pert%primary_molal(jcomp)*perturbation_tolerance
      auxvar_pert%primary_molal(jcomp) = auxvar_pert%primary_molal(jcomp) + pert
      
      ! this is essentially what RTAuxVarCompute() performs
!      call RTAuxVarCompute(auxvar_pert%primary_molal,auxvar_pert,option)      
      auxvar_pert%primary_spec = auxvar_pert%primary_molal* &
                                 auxvar_pert%den(1)*1.d-3
      call RTotal(auxvar_pert,reaction,option)
      if (reaction%nsorb > 0) call RTotalSorb(auxvar_pert,reaction,option)

      ! #3: add new reactions here
      if (reaction%nkinmnrl > 0) then
        call RKineticMineral(Res_pert,Jac_dummy,compute_derivative,auxvar_pert, &
                             volume,reaction,option)
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
    call RTAuxVarDestroy(auxvar_pert)
  endif

end subroutine RReactionDerivative
                               
! ************************************************************************** !
!
! RActivity: Computes the ionic strength and activity coefficients
! author: Glenn Hammond
! date: 09/30/08
!
! ************************************************************************** !
subroutine RActivity(auxvar,reaction,option)

  use Option_module
  
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: icplx, icomp
  PetscReal :: I, sqrt_I
  PetscReal, parameter :: log_to_ln = 2.30258509299d0
  PetscReal, parameter :: ln_to_log = 0.434294481904d0
  ! compute ionic strength
  ! primary species
  I = 0.d0
  do icomp = 1, reaction%ncomp
    I = I + auxvar%primary_spec(icomp)*reaction%primary_spec_Z(icomp)* &
                                       reaction%primary_spec_Z(icomp)
  enddo
  
  ! secondary species
  do icplx = 1, reaction%neqcmplx ! for each secondary species
    I = I + auxvar%secondary_spec(icplx)*reaction%eqcmplx_Z(icplx)* &
                                         reaction%eqcmplx_Z(icplx)
  enddo
  I = I/auxvar%den(1)*1000.d0 ! molarity -> molality
  I = 0.5d0*I
  sqrt_I = sqrt(I)
  
  ! compute activity coefficients
  ! primary species
  do icomp = 1, reaction%ncomp
    auxvar%pri_act_coef(icomp) = exp((-1.d0* &
                                      reaction%primary_spec_Z(icomp)* &
                                      reaction%primary_spec_Z(icomp)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%primary_spec_a0(icomp)* &
                                            reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                     log_to_ln)
  enddo
                
  ! secondary species
  do icplx = 1, reaction%neqcmplx
    auxvar%sec_act_coef(icplx) = exp((-1.d0* &
                                      reaction%eqcmplx_Z(icplx)* &
                                      reaction%eqcmplx_Z(icplx)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%eqcmplx_a0(icplx)* &
                                            reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                     log_to_ln)
  enddo
  
end subroutine RActivity

! ************************************************************************** !
!
! RTotal: Computes the total component concentrations and derivative with
!         respect to free-ion
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine RTotal(auxvar,reaction,option)

  use Option_module
  
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, icplx, icomp, jcomp, iphase, ncomp
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: ln_act_h2o
  PetscReal :: lnQK, tempreal
  PetscReal, parameter :: log_to_ln = 2.30258509299d0

  iphase = 1                         

  ln_conc = log(auxvar%primary_spec)
  ln_act = ln_conc+log(auxvar%pri_act_coef)
  ln_act_h2o = 0.d0  ! assume act h2o = 1 for now
  auxvar%total(:,iphase) = auxvar%primary_spec
  ! initialize derivatives
  auxvar%dtotal = 0.d0
  do icomp = 1, reaction%ncomp
    auxvar%dtotal(icomp,icomp,iphase) = 1.d0
  enddo
  
  do icplx = 1, reaction%neqcmplx ! for each secondary species
    ! compute secondary species concentration
    lnQK = -1.d0*reaction%eqcmplx_logK(icplx)*log_to_ln

    ! activity of water
    if (reaction%eqcmplxh2oid(icplx) > 0) then
      lnQK = lnQK + reaction%eqcmplxh2ostoich(icplx)*ln_act_h2o
    endif

    ncomp = reaction%eqcmplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = reaction%eqcmplxspecid(i,icplx)
      lnQK = lnQK + reaction%eqcmplxstoich(i,icplx)*ln_act(icomp)
    enddo
    auxvar%secondary_spec(icplx) = exp(lnQK)/auxvar%sec_act_coef(icplx)
  
    ! add contribution to primary totals
    ! units of total = mol/L
    do i = 1, ncomp
      icomp = reaction%eqcmplxspecid(i,icplx)
      auxvar%total(icomp,iphase) = auxvar%total(icomp,iphase) + &
                                   reaction%eqcmplxstoich(i,icplx)* &
                                   auxvar%secondary_spec(icplx)
    enddo
    
    ! add contribution to derivatives of total with respect to free
    do j = 1, ncomp
      jcomp = reaction%eqcmplxspecid(j,icplx)
      tempreal = reaction%eqcmplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 auxvar%sec_act_coef(icplx)
      do i = 1, ncomp
        icomp = reaction%eqcmplxspecid(i,icplx)
        auxvar%dtotal(icomp,jcomp,iphase) = auxvar%dtotal(icomp,jcomp,iphase) + &
                                      reaction%eqcmplxstoich(i,icplx)*tempreal
      enddo
    enddo
  enddo
  
  ! convert from dpsi/dc to dpsi/dm where c = rho*m
  ! units of dtotal = kg water/m^3 water
  auxvar%dtotal = auxvar%dtotal*auxvar%den(iphase)
  
end subroutine RTotal

! ************************************************************************** !
!
! RTotalSorb: Computes the total sorbed component concentrations and 
!             derivative with respect to free-ion
! author: Glenn Hammond
! date: 10/22/08
!
! ************************************************************************** !
subroutine RTotalSorb(auxvar,reaction,option)

  use Option_module
  
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, iphase, ncomp, ncplx
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: surfcmplx_conc(reaction%neqsurfcmplx)
  PetscReal :: free_site_conc
  PetscReal :: ln_free_site(reaction%neqsurfcmplxrxn)
  PetscReal :: ln_act_h2o
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscReal, parameter :: log_to_ln = 2.30258509299d0
  PetscReal, parameter :: tol = 1.d-12
  PetscTruth :: one_more
  PetscReal :: res, dres_dfree_site, dfree_site_conc
  
  PetscReal :: omega
  PetscReal :: ref_cation_X, ref_cation_conc, ref_cation_Z, ref_cation_k, &
               ref_cation_quotient
  PetscReal :: cation_X(reaction%ncomp)
  PetscReal :: dres_dref_cation_X, dref_cation_X
  PetscReal :: sumZX
  
  PetscReal :: total_pert, ref_cation_X_pert, pert
  PetscReal :: ref_cation_quotient_pert, dres_dref_cation_X_pert

  iphase = 1                         

  ln_conc = log(auxvar%primary_spec)
  ln_act = ln_conc+log(auxvar%pri_act_coef)
  ln_act_h2o = 0.d0  ! assume act h2o = 1 for now
  ln_free_site = log(auxvar%eqsurfcmplx_freesite_conc)
    
  auxvar%total_sorb(:) = 0.d0
  ! initialize derivatives
  auxvar%dtotal_sorb = 0.d0

  ! Surface Complexation
  do irxn = 1, reaction%neqsurfcmplxrxn
  
    ncplx = reaction%eqsurfcmplx_rxn_to_complex(0,irxn)
    
    free_site_conc = auxvar%eqsurfcmplx_freesite_conc(irxn)

    ! get a pointer to the first complex (there will always be at least 1)
    ! in order to grab free site conc
    one_more = PETSC_FALSE
    do

      total = free_site_conc
      do j = 1, ncplx
        icplx = reaction%eqsurfcmplx_rxn_to_complex(j,irxn)
        ! compute secondary species concentration
        lnQK = -1.d0*reaction%eqsurfcmplx_logK(icplx)*log_to_ln

        ! activity of water
        if (reaction%eqsurfcmplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqsurfcmplxh2ostoich(icplx)*ln_act_h2o
        endif

        lnQK = lnQK + reaction%eqsurfcmplx_free_site_stoich(icplx)* &
                      ln_free_site(irxn)
      
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
    
    auxvar%eqsurfcmplx_freesite_conc(irxn) = free_site_conc
 
    do k = 1, ncplx
      icplx = reaction%eqsurfcmplx_rxn_to_complex(k,irxn)

      auxvar%eqsurfcmplx_spec(icplx) = surfcmplx_conc(icplx)

      ncomp = reaction%eqsurfcmplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsurfcmplxspecid(i,icplx)
        auxvar%total_sorb(icomp) = auxvar%total_sorb(icomp) + &
          reaction%eqsurfcmplxstoich(i,icplx)*surfcmplx_conc(icplx)
      enddo
      
      do j = 1, ncomp
        jcomp = reaction%eqsurfcmplxspecid(j,icplx)
        tempreal = reaction%eqsurfcmplxstoich(j,icplx)*surfcmplx_conc(icplx) / &
          auxvar%primary_spec(jcomp)
        do i = 1, ncomp
          icomp = reaction%eqsurfcmplxspecid(i,icplx)
          auxvar%dtotal_sorb(icomp,jcomp) = auxvar%dtotal_sorb(icomp,jcomp) + &
                                        reaction%eqsurfcmplxstoich(i,icplx)*tempreal
        enddo
      enddo
    enddo
  enddo
  
  ! Ion Exchange
  do irxn = 1, reaction%neqionxrxn

    ncomp = reaction%eqionx_rxn_cationid(0,irxn)

    ! for now we assume that omega is equal to CEC.
    omega = reaction%eqionx_rxn_CEC(irxn)

    icomp = reaction%eqionx_rxn_cationid(1,irxn)
    ref_cation_conc = auxvar%primary_spec(icomp)
    ref_cation_Z = reaction%primary_spec_Z(icomp)
    ref_cation_k = reaction%eqionx_rxn_k(1,irxn)
    ref_cation_X = ref_cation_Z*auxvar%eqionx_ref_cation_sorbed_conc(irxn)/omega

    one_more = PETSC_FALSE
    cation_X = 0.d0
    do

      if (ref_cation_X <= 0.d0) ref_cation_X = 0.99d0
      cation_X(1) = ref_cation_X
      ref_cation_quotient = ref_cation_X*ref_cation_k/ref_cation_conc
      total = ref_cation_X

      if (reaction%eqionx_rxn_Z_flag(irxn)) then ! Zi /= Zj for any i,j

        do j = 2, ncomp
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          cation_X(j) = auxvar%primary_spec(icomp)/reaction%eqionx_rxn_k(j,irxn)* &
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
      total_pert = ref_cation_X

        do j = 2, ncomp
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          total_pert = total_pert + &
                       auxvar%primary_spec(icomp)/reaction%eqionx_rxn_k(j,irxn)* &
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

        dref_cation_X = res / -dres_dref_cation_X
!        dref_cation_X = res / dres_dref_cation_X_pert
        ref_cation_X = ref_cation_X - dref_cation_X
      
        if (dabs(dref_cation_X/ref_cation_X) < tol) then
          one_more = PETSC_TRUE
        endif
    
      else
      
        do j = 2, ncomp  ! Zi == Zj for all i,j
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          cation_X(j) = auxvar%primary_spec(icomp)/reaction%eqionx_rxn_k(j,irxn)* &
                        ref_cation_quotient
          total = total + cation_X(j)
        enddo
        
        if (one_more) exit
        
        total = total / ref_cation_X
        ref_cation_X = omega / total  
        
        one_more = PETSC_TRUE 
      
      endif
      
    enddo
    auxvar%eqionx_ref_cation_sorbed_conc(irxn) = ref_cation_X*omega/ref_cation_Z

    ! sum up charges
    do i = 1, ncomp
    icomp = reaction%eqionx_rxn_cationid(i,irxn)
      sumZX = sumZX + reaction%primary_spec_Z(icomp)*cation_X(i)
    enddo

    ! compute totals based on sorbed ions
    do j = 1, ncomp
      jcomp = reaction%eqionx_rxn_cationid(j,irxn)
      tempreal1 = cation_X(j)*omega/reaction%primary_spec_Z(jcomp)
      ! residual function entry
      auxvar%total_sorb(jcomp) = auxvar%total_sorb(jcomp) + tempreal1

      tempreal2 = reaction%primary_spec_Z(jcomp)/sumZX
      do i = 1, ncomp
        icomp = reaction%eqionx_rxn_cationid(i,irxn)
        if (i == j) then
          auxvar%dtotal_sorb(icomp,jcomp) = auxvar%dtotal_sorb(icomp,jcomp) + &
                                            tempreal1*(1.d0-(tempreal2*cation_X(i)))/ &
                                            auxvar%primary_spec(icomp)
        else
          auxvar%dtotal_sorb(icomp,jcomp) = auxvar%dtotal_sorb(icomp,jcomp) + &
                                            -tempreal1*tempreal2*cation_X(i)/ &
                                            auxvar%primary_spec(icomp)
        endif
      enddo
    enddo    

  enddo

  ! convert from dpsi/dc to dpsi/dm where c = rho*m
  ! units of dtotal = kg water/m^3 water
  auxvar%dtotal_sorb = auxvar%dtotal_sorb*auxvar%den(iphase)
  
end subroutine RTotalSorb

! ************************************************************************** !
!
! RKineticMineral: Computes the kinetic mineral precipitation/dissolution
!                  rates
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RKineticMineral(Res,Jac,compute_derivative,auxvar,volume, &
                           reaction,option)

  use Option_module
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscTruth :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: auxvar
  
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
  PetscReal, parameter :: log_to_ln = 2.30258509299d0
  PetscTruth :: prefactor_exists

  iphase = 1                         

  ln_conc = log(auxvar%primary_spec)
  ln_sec = log(auxvar%secondary_spec)
  
  ln_act = ln_conc+log(auxvar%pri_act_coef)
  ln_sec_act = ln_sec+log(auxvar%sec_act_coef)
  ln_act_h2o = 0.d0
  
  do imnrl = 1, reaction%nkinmnrl ! for each mineral
    ! compute secondary species concentration
    lnQK = -1.d0*reaction%kinmnrl_logK(imnrl)*log_to_ln

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

    if (auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then
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
          do k = 1, reaction%kinmnrl_sec_prefactor_id(0,ipref,imnrl) ! secondary contribution
            kcplx = reaction%kinmnrl_sec_prefactor_id(k,ipref,imnrl)
            prefactor(ipref) = prefactor(ipref) * &
                               exp(reaction%kinmnrl_sec_pref_alpha_stoich(k,ipref,imnrl)* &
                                   ln_sec_act(kcplx))/ &
                               ((1.d0+reaction%kinmnrl_sec_pref_atten_coef(i,ipref,imnrl))* &
                                 exp(reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                     ln_sec_act(kcplx)))
          enddo
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
      Im_const = -1.d0*auxvar%mnrl_area0(imnrl)*1.d6 ! convert cm^3->m^3
      ! units = mol/sec/m^3 bulk
      if (associated(reaction%kinmnrl_affinity_power)) then
        Im = Im_const*sign_*abs(affinity_factor)**reaction%kinmnrl_affinity_power(imnrl)*sum_prefactor_rate
      else
        Im = Im_const*sign_*abs(affinity_factor)*sum_prefactor_rate
      endif
      auxvar%mnrl_rate(imnrl) = Im ! mol/sec/m^3
    else
      auxvar%mnrl_rate(imnrl) = 0.d0
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
      dIm_dQK = -1.d0*Im*reaction%kinmnrl_affinity_power(imnrl)/abs(affinity_factor)
    else
      dIm_dQK = -1.d0*Im_const*sum_prefactor_rate
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
      dQK_dmj = dQK_dCj*auxvar%den(iphase)*1.d-3 ! the multiplication by density could be moved
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
                                       prefactor(ipref)/auxvar%primary_spec(jcomp) ! dR_dc
          ! denominator
          dprefactor_dcomp_denominator = -1.d0*prefactor(ipref)/ &
                                         ((1.d0+reaction%kinmnrl_pri_pref_atten_coef(j,ipref,imnrl))* &
                                           exp(reaction%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)* &
                                               ln_act(jcomp)))* & 
                                         reaction%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)* &
                                         reaction%kinmnrl_pri_pref_atten_coef(j,ipref,imnrl)* &
                                         exp((reaction%kinmnrl_pri_pref_beta_stoich(j,ipref,imnrl)-1.d0)* &
                                             ln_act(jcomp))* & ! dR_da
                                         auxvar%pri_act_coef(jcomp) ! da_dc
          tempreal = dIm_dprefactor_rate*(dprefactor_dcomp_numerator+dprefactor_dcomp_denominator)*auxvar%den(iphase)
          do i = 1, ncomp
            icomp = reaction%kinmnrlspecid(i,imnrl)
            Jac(icomp,jcomp) = Jac(icomp,jcomp) + reaction%kinmnrlstoich(i,imnrl)*tempreal
          enddo  ! loop over col
        enddo !loop over row
        do k = 1, reaction%kinmnrl_sec_prefactor_id(0,ipref,imnrl) ! secondary contribution
          kcplx = reaction%kinmnrl_sec_prefactor_id(k,ipref,imnrl)
          ! numerator
          dprefactor_dcomp_numerator = reaction%kinmnrl_sec_pref_alpha_stoich(k,ipref,imnrl)* &
                                       prefactor(ipref)/(auxvar%secondary_spec(kcplx)* &
                                                         auxvar%sec_act_coef(kcplx)) ! dR_dax
          ! denominator
          dprefactor_dcomp_denominator = -1.d0*prefactor(ipref)/ &
                                         (1.d0+reaction%kinmnrl_sec_pref_atten_coef(k,ipref,imnrl)* &
                                          exp(reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                              ln_sec_act(kcplx)))* &
                                         reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)* &
                                         reaction%kinmnrl_sec_pref_atten_coef(k,ipref,imnrl)* &
                                         exp((reaction%kinmnrl_sec_pref_beta_stoich(k,ipref,imnrl)-1.d0)* &
                                              ln_sec_act(kcplx)) ! dR_dax
          tempreal = dIm_dprefactor_rate*(dprefactor_dcomp_numerator+dprefactor_dcomp_denominator)*auxvar%den(iphase)
          do j = 1, reaction%eqcmplxstoich(0,kcplx)
            jcomp = reaction%eqcmplxstoich(j,kcplx)
            tempreal2 = reaction%eqcmplxstoich(j,kcplx)*exp(ln_sec_act(kcplx)-ln_conc(jcomp)) !dax_dc
            do i = 1, ncomp
              icomp = reaction%kinmnrlspecid(i,imnrl)
              Jac(icomp,jcomp) = Jac(icomp,jcomp) + reaction%kinmnrlstoich(i,imnrl)*tempreal* &
                                                    tempreal2
            enddo  ! loop over col
          enddo  ! loop over row
        enddo  ! loop over complexes
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
