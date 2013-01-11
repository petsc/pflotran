module Microbial_module

  use Microbial_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  public :: MicrobialRead, &
            RMicrobial, &
            MicrobialProcessConstraint

contains

! ************************************************************************** !
!
! MicrobialRead: Reads chemical species
! author: Glenn Hammond
! date: 08/16/12
!
! ************************************************************************** !
subroutine MicrobialRead(microbial,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(microbial_type) :: microbial
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  type(microbial_rxn_type), pointer :: microbial_rxn, cur_microbial_rxn
  type(monod_type), pointer :: monod, prev_monod
  type(biomass_type), pointer :: biomass
  type(inhibition_type), pointer :: inhibition, prev_inhibition
  
  microbial%nrxn = microbial%nrxn + 1
        
  microbial_rxn => MicrobialRxnCreate()
  nullify(prev_monod)
  nullify(prev_inhibition)
  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MICROBIAL_REACTION')
    call StringToUpper(word)   

    select case(trim(word))
      case('REACTION')
        ! remainder of string should be the reaction equation
        microbial_rxn%reaction = trim(adjustl(input%buf))
        ! set flag for error message
        if (len_trim(microbial_rxn%reaction) < 2) input%ierr = 1
        call InputErrorMsg(input,option,'reaction string', &
                            'CHEMISTRY,MICROBIAL_REACTION,REACTION')     
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,microbial_rxn%rate_constant)  
        call InputDefaultMsg(input,option, &
                              'CHEMISTRY,MICROBIAL_REACTION,RATE_CONSTANT') 
      case('MONOD')
        monod => MicrobialMonodCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,MICROBIAL_REACTION,MONOD')
        monod%species_name = word
        call InputReadDouble(input,option,monod%half_saturation_constant)  
        call InputErrorMsg(input,option,'half saturation constant', &
                           'CHEMISTRY,MICROBIAL_REACTION,MONOD')
        ! append to list
        if (.not.associated(microbial_rxn%monod)) then
          microbial_rxn%monod => monod
        else
          prev_monod%next => monod
        endif
        prev_monod => monod
        nullify(monod)
      case('INHIBITION')
        inhibition => MicrobialInhibitionCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
        inhibition%species_name = word
        call InputReadDouble(input,option,inhibition%inhibition_constant)  
        call InputErrorMsg(input,option,'inhibition constant', &
                           'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
        ! append to list
        if (.not.associated(microbial_rxn%inhibition)) then
          microbial_rxn%inhibition => inhibition
        else
          prev_inhibition%next => inhibition
        endif
        prev_inhibition => inhibition
        nullify(inhibition)
      case('BIOMASS')
        biomass => MicrobialBiomassCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
        biomass%species_name = word
        call InputReadDouble(input,option,biomass%yield)  
        call InputErrorMsg(input,option,'yield', &
                           'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
      case default
        option%io_buffer = 'CHEMISTRY,MICROBIAL_REACTION keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo   
  
  ! add to list
  if (.not.associated(microbial%microbial_rxn_list)) then
    microbial%microbial_rxn_list => microbial_rxn
    microbial_rxn%id = 1
  else
    cur_microbial_rxn => microbial%microbial_rxn_list
    do
      if (.not.associated(cur_microbial_rxn%next)) then
        cur_microbial_rxn%next => microbial_rxn
        microbial_rxn%id = cur_microbial_rxn%id + 1
        exit
      endif
      cur_microbial_rxn => cur_microbial_rxn%next
    enddo
  endif
  nullify(microbial_rxn)

end subroutine MicrobialRead

! ************************************************************************** !
!
! MicrobialBiomassRead: Reads biomass species
! author: Glenn Hammond
! date: 01/02/13
!
! ************************************************************************** !
subroutine MicrobialBiomassRead(microbial,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(microbial_type) :: microbial
  type(input_type) :: input
  type(option_type) :: option
  
  type(biomass_species_type), pointer :: biomass, prev_biomass
           
  nullify(prev_biomass)
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
          
    microbial%nbiomass = microbial%nbiomass + 1
          
    biomass => MicrobialBiomassSpeciesCreate()
    call InputReadWord(input,option,biomass%name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MINERALS')    
    if (.not.associated(microbial%biomass_list)) then
      microbial%biomass_list => biomass
      biomass%id = 1
    endif
    if (associated(prev_biomass)) then
      prev_biomass%next => biomass
      biomass%id = prev_biomass%id + 1
    endif
    prev_biomass => biomass
    nullify(biomass)
  enddo

end subroutine MicrobialBiomassRead

! ************************************************************************** !
!
! MicrobialProcessConstraint: Initializes constraints based on biomass
!                             species in system
! author: Glenn Hammond
! date: 01/07/13
!
! ************************************************************************** !
subroutine MicrobialProcessConstraint(microbial,constraint_name, &
                                      constraint,option)
  use Option_module
  use Input_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(microbial_type), pointer :: microbial
  character(len=MAXWORDLENGTH) :: constraint_name
  type(biomass_constraint_type), pointer :: constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: ibiomass, jbiomass
  
  character(len=MAXWORDLENGTH) :: biomass_name(microbial%nbiomass)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(microbial%nbiomass)
  PetscReal :: constraint_conc(microbial%nbiomass)
  PetscBool :: external_dataset(microbial%nbiomass)
  
  if (.not.associated(constraint)) return
  
  biomass_name = ''
  constraint_aux_string = ''
  external_dataset = PETSC_FALSE
  do ibiomass = 1, microbial%nbiomass
    found = PETSC_FALSE
    do jbiomass = 1, microbial%nbiomass
      if (StringCompare(constraint%names(ibiomass), &
                        microbial%biomass_names(jbiomass), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
                'Biomass species "' // trim(constraint%names(ibiomass)) // &
                '" from CONSTRAINT "' // trim(constraint_name) // &
                '" not found among biomass species.'
      call printErrMsg(option)
    else
      biomass_name(ibiomass) = constraint%names(ibiomass)
      constraint_conc(ibiomass) = &
        constraint%constraint_conc(ibiomass)
      constraint_aux_string(ibiomass) = &
        constraint%constraint_aux_string(ibiomass)
      external_dataset(ibiomass) = constraint%external_dataset(ibiomass)
    endif  
  enddo
  constraint%names = biomass_name
  constraint%constraint_conc = constraint_conc
  constraint%constraint_aux_string = constraint_aux_string
  constraint%external_dataset = external_dataset

end subroutine MicrobialProcessConstraint

! ************************************************************************** !
!
! RMicrobial: Computes the microbial reaction
! author: Glenn Hammond
! date: 10/31/12
!
! ************************************************************************** !
subroutine RMicrobial(Res,Jac,compute_derivative,rt_auxvar, &
                      global_auxvar,porosity,volume,reaction,option)

  use Option_module, only : option_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type
  use Reaction_Aux_module, only : reaction_type
  
  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscInt, parameter :: iphase = 1
  PetscReal :: por_sat_vol
  PetscInt :: irxn, i, ii, icomp, jcomp, ncomp
  PetscInt :: imonod, iinhibition, ibiomass
  PetscReal :: Im
  PetscReal :: rate_constant
  PetscReal :: activity
  PetscReal :: act_coef
  PetscReal :: monod(10)
  PetscReal :: inhibition(10)
  PetscReal :: biomass
  PetscReal :: denominator, dR_dX, dX_dc, dR_dc
  type(microbial_type), pointer :: microbial
  
  microbial => reaction%microbial
  
  ! units:
  ! Residual: mol/sec
  ! Jacobian: (mol/sec)*(kg water/mol) = kg water/sec
  
  do irxn = 1, microbial%nrxn
  
    ! units:
    !   without biomass: mol/L-sec
    !   with biomass: mol/L-sec * (m^3 bulk / mol biomass)
    rate_constant = microbial%rate_constant(irxn)
    Im = rate_constant

    ! monod expressions
    do ii = 1, microbial%monodid(0,irxn)
      imonod = microbial%monodid(ii,irxn)
      icomp = microbial%monod_specid(imonod)
      activity = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)
      monod(ii) = activity / (microbial%monod_K(imonod) + activity)
      Im = Im*monod(ii)
    enddo

    ! inhibition expressions
    do ii = 1, microbial%inhibitionid(0,irxn)
      iinhibition = microbial%inhibitionid(ii,irxn)
      icomp = microbial%inhibition_specid(iinhibition)
      activity = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)
      inhibition(ii) = microbial%inhibition_C(iinhibition) / &
                      (microbial%inhibition_C(iinhibition) + activity)
      Im = Im*inhibition(ii)
    enddo
    
    ! biomass term
    ibiomass = microbial%biomassid(irxn)
    if (ibiomass > 0) then
      biomass = rt_auxvar%immobile(microbial%biomassid(ibiomass))
      Im = Im*biomass
    endif
    
    ! por_sat_vol units: m^3 water
    por_sat_vol = porosity*global_auxvar%sat(iphase)*volume

    ! Im units (before): mol/L-sec
    Im = Im * 1.d3*por_sat_vol
    ! Im units (after): mol/sec
    
    ncomp = microbial%specid(0,irxn)
    do i = 1, ncomp
      icomp = microbial%specid(i,irxn)
      Res(icomp) = Res(icomp) - microbial%stoich(i,irxn)*Im
    enddo
    
    if (.not. compute_derivative) cycle
    
    ! monod expressions
    do ii = 1, microbial%monodid(0,irxn)
      imonod = microbial%monodid(ii,irxn)
      jcomp = microbial%monod_specid(imonod)
      act_coef = rt_auxvar%pri_act_coef(jcomp)
      activity = rt_auxvar%pri_molal(jcomp)*act_coef
        
      dR_dX = Im / monod(ii)
        
      denominator = microbial%monod_K(imonod) + activity
        
      dX_dc = act_coef / denominator - &
              act_coef * activity / (denominator*denominator)
        
      dR_dc = -1.d0*dR_dX*dX_dc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                            microbial%stoich(i,irxn)*dR_dc
      enddo
    enddo

    ! inhibition expressions
    do ii = 1, microbial%inhibitionid(0,irxn)
      iinhibition = microbial%inhibitionid(ii,irxn)
      jcomp = microbial%inhibition_specid(iinhibition)
      act_coef = rt_auxvar%pri_act_coef(jcomp)
      activity = rt_auxvar%pri_molal(jcomp)*act_coef

      dR_dX = Im / inhibition(ii)
        
      denominator = microbial%inhibition_C(iinhibition) + activity
        
      dX_dc = -1.d0 * act_coef *microbial%inhibition_C(iinhibition) / &
              (denominator*denominator)
        
      dR_dc = -1.d0*dR_dX*dX_dc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                            microbial%stoich(i,irxn)*dR_dc
      enddo
    enddo

  enddo
    
end subroutine RMicrobial

end module Microbial_module
