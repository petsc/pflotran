module Microbial_module

  use Microbial_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  public :: MicrobialRead, &
            RMicrobial

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
  type(microbial_biomass_type), pointer :: microbial_biomass
  type(inhibition_type), pointer :: inhibition, prev_inhibition
  
  microbial%nrxn = microbial%nrxn + 1
        
  microbial_rxn => MicrobialRxnCreate()
  nullify(prev_monod)
  nullify(prev_inhibition)
  nullify(microbial_biomass)
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
        if (associated(microbial_biomass)) then
          call MicrobialBiomassDestroy(microbial_biomass)
        endif
        microbial_biomass => MicrobialBiomassCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
        microbial_biomass%species_name = word
        call InputReadDouble(input,option,microbial_biomass%yield)  
        call InputErrorMsg(input,option,'yield', &
                           'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
      case default
        option%io_buffer = 'CHEMISTRY,MICROBIAL_REACTION keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
  ! add linkage to biomass if exists
  microbial_rxn%biomass => microbial_biomass
  
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
! RMicrobial: Computes the microbial reaction
! author: Glenn Hammond
! date: 10/31/12
!
! ************************************************************************** !
subroutine RMicrobial(Res,Jac,compute_derivative,rt_auxvar, &
                      global_auxvar,porosity,volume,reaction,option)

  use Option_module, only : option_type, printErrMsg
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type
  use Reaction_Aux_module, only : reaction_type
  use Immobile_Aux_module, only : immobile_type
  
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
  PetscReal :: biomass_conc
  PetscInt :: immobile_id
  PetscReal :: denominator, dR_dX, dX_dc, dR_dc, dR_dbiomass
  type(microbial_type), pointer :: microbial
  type(immobile_type), pointer :: immobile
  
  microbial => reaction%microbial
  immobile => reaction%immobile
  
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
      immobile_id = reaction%offset_immobile + ibiomass
      biomass_conc = rt_auxvar%immobile(ibiomass)
      Im = Im*biomass_conc
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

    if (ibiomass > 0) then
      Res(immobile_id) = Res(immobile_id) - microbial%biomass_yield(irxn)*Im
    endif
    
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
      if (ibiomass > 0) then
        Jac(immobile_id,jcomp) = Jac(immobile_id,jcomp) + &
          microbial%biomass_yield(irxn)*dR_dc      
      endif      
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
      if (ibiomass > 0) then
        Jac(immobile_id,jcomp) = Jac(immobile_id,jcomp) + &
          microbial%biomass_yield(irxn)*dR_dc      
      endif
    enddo

    ! biomass expression
    if (ibiomass > 0) then
      dR_dbiomass = Im / biomass_conc
      option%io_buffer = "Shouldn't biomass contribution be negative"
      call printErrMsg(option)
      do i = 1, ncomp
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,immobile_id) = Jac(icomp,immobile_id) + &
                            microbial%stoich(i,irxn)*dR_dbiomass
      enddo
      Jac(immobile_id,immobile_id) = Jac(immobile_id,immobile_id) + &
        microbial%biomass_yield(irxn)*dR_dbiomass
    endif
    
  enddo
    
end subroutine RMicrobial

end module Microbial_module
