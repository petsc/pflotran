module Reaction_module

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: ReactionCreate, &
            ReactionRead, &
            CarbonateTestProblemCreate, &
            ReactionReadMineralKinetics, &
!           ReactionReadSurfaceComplexes, &
            ReactionInitializeConstraint, &
            RTotal, &
            RActivity, &
            RReaction, &
            RReactionDerivative

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
  type(surface_complex_type), pointer :: srfcmplx, prev_srfcmplx
  type(surface_complexation_rxn_type), pointer :: srfcmplx_rxn, prev_srfcmplx_rxn
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
              nullify(prev_srfcmplx)
              do
                call fiReadFlotranString(fid,string,ierr)
                if (ierr /= 0) exit
                if (fiCheckExit(string)) exit

                call fiReadNChars(string,word,MAXNAMELENGTH,.true.,ierr)
                call fiErrorMsg(option%myrank,'keyword','CHEMISTRY', ierr)
                call fiWordToUpper(word)
                
                select case(trim(word))
                  case('MINERAL')
                    call fiReadNChars(string,srfcmplx_rxn%mineral_name,MAXNAMELENGTH,PETSC_TRUE,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,MINERAL_NAME', ierr)
                  case('SITE')
                    call fiReadNChars(string,srfcmplx_rxn%free_site_name,MAXNAMELENGTH,PETSC_TRUE,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_NAME', ierr)
                    call fiReadDouble(string,srfcmplx_rxn%site_density,ierr)
                    call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_DENSITY',ierr)                   
                  case('COMPLEXES')
                    nullify(srfcmplx)
                    do
                      call fiReadFlotranString(fid,string,ierr)
                      if (ierr /= 0) exit
                      if (fiCheckExit(string)) exit
                      
                      srfcmplx_count = srfcmplx_count + 1
                      srfcmplx => SurfaceComplexCreate()
                      srfcmplx%id = srfcmplx_count
                      call fiReadNChars(string,srfcmplx%name,MAXNAMELENGTH,PETSC_TRUE,ierr)
                      call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_NAME', ierr)
                
                      if (.not.associated(srfcmplx_rxn%complex_list)) then
                        srfcmplx_rxn%complex_list => srfcmplx
                      endif
                      if (associated(prev_srfcmplx_rxn)) then
                        prev_srfcmplx%next => srfcmplx
                      endif
                      prev_srfcmplx => srfcmplx
                      nullify(srfcmplx)
                
                    enddo
                end select

                if (.not.associated(reaction%surface_complexation_rxn_list)) then
                  reaction%surface_complexation_rxn_list => srfcmplx_rxn
                  srfcmplx%id = 1
                endif
                if (associated(prev_srfcmplx_rxn)) then
                  prev_srfcmplx_rxn%next => srfcmplx_rxn
                  srfcmplx%id = prev_srfcmplx%id + 1
                endif
                prev_srfcmplx_rxn => srfcmplx_rxn
                nullify(srfcmplx_rxn)
              enddo
            case('ION_EXCHANGE_RXN')
          !   call IonExchangeRXNRead
            case('DISTRIBUTION_COEF')
          !   call DistributionCoefRead
          end select
        enddo
      case('DATABASE')
        call fiReadNChars(string,reaction%database_filename,MAXSTRINGLENGTH,PETSC_TRUE,ierr)  
        call fiErrorMsg(option%myrank,'keyword','CHEMISTRY,DATABASE FILENAME', ierr)          
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
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: found
  PetscInt :: icomp, jcomp
  PetscInt :: imnrl, jmnrl
  PetscReal :: value
  PetscReal :: conc(option%ncomp)
  PetscInt :: constraint_type(option%ncomp)
  character(len=MAXNAMELENGTH) :: constraint_spec_name(option%ncomp)


  type(reactive_transport_auxvar_type) :: auxvar
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: total_conc(reaction%ncomp)
  PetscReal :: free_conc(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscInt :: indices(reaction%ncomp)
  PetscInt :: num_it
  PetscReal :: norm
  PetscReal :: prev_molal(reaction%ncomp)
  PetscReal, parameter :: tol = 1.d-6
    
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
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among primary species.'
      call printErrMsg(option,string)
    else
      constraint_type(jcomp) = aq_species_constraint%constraint_type(icomp)
      constraint_spec_name(jcomp) = aq_species_constraint%constraint_spec_name(icomp)
      conc(jcomp) = aq_species_constraint%conc(icomp)
    endif
  enddo
  
  if (.not.associated(reaction)) then ! simply tracer transport
    aq_species_constraint%basis_conc = conc
  else
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
          mineral_constraint%basis_conc(jmnrl) = mineral_constraint%conc(imnrl)
        endif  
      enddo
    endif

    total_conc = 0.d0
    do icomp = 1, reaction%ncomp
      select case(constraint_type(icomp))
        case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
          total_conc(icomp) = conc(icomp)
          free_conc(icomp) = 1.d-9
        case(CONSTRAINT_FREE)
          free_conc(icomp) = conc(icomp)
        case(CONSTRAINT_P)
          free_conc(icomp) = 10**(-1.d0*conc(icomp))
        case(CONSTRAINT_MINERAL,CONSTRAINT_GAS)
          free_conc(icomp) = conc(icomp) ! guess
      end select
    enddo
    
    call RTAuxVarInit(auxvar,option)
    auxvar%den(1) = 997.d0
    auxvar%primary_molal = free_conc

    num_it = 0
    
    do

      auxvar%primary_spec = auxvar%primary_molal* &
                            auxvar%den(1)*1.d-3
      call RTotal(auxvar,reaction,option)
      
      Res = auxvar%total(:,1)
      Jac = auxvar%dtotal(:,:,1)
          
      do icomp = 1, reaction%ncomp
        select case(constraint_type(icomp))
          case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
          case(CONSTRAINT_FREE,CONSTRAINT_P)
            Res(icomp) = 0.d0
            Jac(icomp,:) = 0.d0
            Jac(:,icomp) = 0.d0
            Jac(icomp,icomp) = 1.d0
          case(CONSTRAINT_MINERAL)
          case(CONSTRAINT_GAS)
        end select
      enddo
      
      Res = 1.d0*(total_conc - Res)

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
      print *, num_it, Res
      
      ! check for convergence
      if (maxval(dabs(auxvar%primary_molal-prev_molal)/ &
                 auxvar%primary_molal) < tol) exit
                       
    enddo
    
    aq_species_constraint%basis_conc = auxvar%primary_molal

    call RTAuxVarDestroy(auxvar)

  endif

end subroutine ReactionInitializeConstraint

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
    compute_derivative = PETSC_TRUE
    ! #1: add new reactions here
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res,Jac,compute_derivative,auxvar,volume,reaction,option)
    endif
  else ! numerical derivative
    compute_derivative = PETSC_FALSE
    call RTAuxVarInit(auxvar_pert,option)
    call RTAuxVarCopy(auxvar_pert,auxvar,option)
    ! #2: add new reactions here
    if (reaction%nkinmnrl > 0) then
      call RKineticMineral(Res_orig,Jac_dummy,compute_derivative,auxvar,volume,reaction,option)
    endif
    do jcomp = 1, reaction%ncomp
      call RTAuxVarCopy(auxvar_pert,auxvar,option)
      pert = auxvar_pert%primary_molal(jcomp)*perturbation_tolerance
      auxvar_pert%primary_molal(jcomp) = auxvar_pert%primary_molal(jcomp) + pert
      
      ! this is essentially what RTAuxVarCompute() performs
!      call RTAuxVarCompute(auxvar_pert%primary_molal,auxvar_pert,option)      
      auxvar_pert%primary_spec = auxvar_pert%primary_molal* &
                                 auxvar_pert%den(1)*1.d-3
      call RTotal(auxvar_pert,reaction,option)
      !

      ! #3: add new reactions here
      if (reaction%nkinmnrl > 0) then
        call RKineticMineral(Res_pert,Jac_dummy,compute_derivative,auxvar_pert, &
                             volume,reaction,option)
      endif
      do icomp = 1, reaction%ncomp
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + (Res_pert(icomp)-Res_orig(icomp))/pert
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
  I = 0.5d0*I
  sqrt_I = sqrt(I)
  
  ! compute activity coefficients
  ! primary species
  do icomp = 1, reaction%ncomp
    auxvar%pri_act_coef(icomp) = exp((reaction%primary_spec_Z(icomp)* &
                                      reaction%primary_spec_Z(icomp)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%primary_spec_a0(icomp)* &
                                            reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                     log_to_ln)
  enddo
                
  ! secondary species
  do icplx = 1, reaction%neqcmplx
    auxvar%sec_act_coef(icplx) = exp((reaction%eqcmplx_Z(icplx)* &
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
! RInitConcentration: Initializaes concentrations based on constraints
! author: Glenn Hammond
! date: 10/20/08
!
! ************************************************************************** !
subroutine RInitConcentration(auxvar,reaction,option)

  use Option_module
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(reactive_transport_auxvar_type) :: auxvar_tmp
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: update(reaction%ncomp)
  PetscInt :: num_it
  PetscInt :: iphase
  PetscInt :: icomp
  PetscReal :: norm
  PetscReal :: prev_molal(reaction%ncomp)
  PetscReal, parameter :: tol = 1.d-12

  call RTAuxVarInit(auxvar_tmp,option)
  call RTAuxVarCopy(auxvar_tmp,auxvar,option)

  num_it = 0
  do
  
    call RTotal(auxvar_tmp,reaction,option)

    Res = 0.d0
    Jac = 0.d0
    do iphase = 1, option%nphase
      Res =  auxvar_tmp%total(:,iphase) - auxvar%total(:,iphase)
      Jac = Jac + auxvar_tmp%dtotal(:,:,iphase)
    enddo
    
    call RSolve(Res,Jac,update,auxvar_tmp%primary_spec,reaction%ncomp)
 
    auxvar_tmp%primary_molal = auxvar_tmp%primary_molal*exp(update)
    
    num_it = num_it + 1
    
    ! check for convergence
    if (maxval(dabs(auxvar_tmp%primary_molal-prev_molal)/ &
               auxvar_tmp%primary_molal) < tol) exit
    
  enddo

  auxvar%primary_molal = auxvar_tmp%primary_molal
  auxvar%primary_spec = auxvar_tmp%primary_spec

  call RTAuxVarDestroy(auxvar_tmp)

end subroutine RInitConcentration

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
