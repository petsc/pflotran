module Mineral_module

  use Mineral_Aux_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: MineralRead, &
            MineralReadKinetics, &
            MineralReadFromDatabase, &
            MineralProcessConstraint, &
            RKineticMineral, &
            RMineralSaturationIndex
            
contains

! ************************************************************************** !
!
! MineralRead: Reads chemical species
! author: Glenn Hammond
! date: 08/16/12
!
! ************************************************************************** !
subroutine MineralRead(mineral,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(mineral_type) :: mineral
  type(input_type) :: input
  type(option_type) :: option
  
  type(mineral_rxn_type), pointer :: cur_mineral, prev_mineral
           
  nullify(prev_mineral)
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
          
    mineral%nmnrl = mineral%nmnrl + 1
          
    cur_mineral => MineralRxnCreate()
    call InputReadWord(input,option,cur_mineral%name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MINERALS')    
    if (.not.associated(mineral%mineral_list)) then
      mineral%mineral_list => cur_mineral
      cur_mineral%id = 1
    endif
    if (associated(prev_mineral)) then
      prev_mineral%next => cur_mineral
      cur_mineral%id = prev_mineral%id + 1
    endif
    prev_mineral => cur_mineral
    nullify(cur_mineral)
  enddo

end subroutine MineralRead

! ************************************************************************** !
!
! MineralReadKinetics: Reads mineral kinetics
! author: Glenn Hammond
! date: 10/16/08
!
! ************************************************************************** !
subroutine MineralReadKinetics(mineral,input,option)

  use Input_module
  use String_module  
  use Option_module
  use Units_module
  
  implicit none
  
  type(mineral_type) :: mineral
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
  
  type(mineral_rxn_type), pointer :: cur_mineral
  type(transition_state_rxn_type), pointer :: tstrxn, cur_tstrxn
  type(transition_state_prefactor_type), pointer :: prefactor, &
                                                    cur_prefactor
  type(ts_prefactor_species_type), pointer :: prefactor_species, &
                                              cur_prefactor_species
  PetscBool :: found
  PetscInt :: imnrl,icount
  PetscReal :: temp_real

  cur_mineral => mineral%mineral_list
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
    
    cur_mineral => mineral%mineral_list
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
              if (tstrxn%rate < 0.d0) then
                tstrxn%rate = 10.d0**tstrxn%rate
              endif
              call InputErrorMsg(input,option,'rate',error_string)
              ! read units if they exist
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = trim(cur_mineral%name) // 'RATE UNITS'
                call InputDefaultMsg(input,option)
              else
                tstrxn%rate = tstrxn%rate * UnitsConvertToInternal(word,option)
              endif
            case('ACTIVATION_ENERGY')
!             read activation energy for Arrhenius law
              call InputReadDouble(input,option,tstrxn%activation_energy)
              call InputErrorMsg(input,option,'activation',error_string)
            case('AFFINITY_THRESHOLD')
!             read affinity threshold for precipitation
              call InputReadDouble(input,option,tstrxn%affinity_threshold)
              call InputErrorMsg(input,option,'affinity threshold', &
                                 error_string)
            case('AFFINITY_POWER')
!             reads exponent on affinity term
              call InputReadDouble(input,option,tstrxn%affinity_factor_sigma)
              call InputErrorMsg(input,option,'affinity power',error_string)
            case('TEMPKINS_CONSTANT')
!             reads exponent on affinity term
              call InputReadDouble(input,option,tstrxn%affinity_factor_beta)
              call InputErrorMsg(input,option,"Tempkin's constant", &
                                 error_string)
            case('SURFACE_AREA_POROSITY_POWER')
              call InputReadDouble(input,option,tstrxn%surf_area_porosity_pwr)
              call InputErrorMsg(input,option,'surface area porosity power', &
                                 error_string)
            case('SURFACE_AREA_VOL_FRAC_POWER')
              call InputReadDouble(input,option,tstrxn%surf_area_vol_frac_pwr)
              call InputErrorMsg(input,option, &
                                 'surface area voluem fraction power', &
                                 error_string)
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
                    ! read rate constant
                    call InputReadDouble(input,option,prefactor%rate)
                    call InputErrorMsg(input,option,'rate',error_string)
                    if (prefactor%rate < 0.d0) then
                      prefactor%rate = 10.d0**prefactor%rate
                    endif
                    ! read units if they exist
                    call InputReadWord(input,option,word,PETSC_TRUE)
                    if (InputError(input)) then
                      input%err_buf = trim(cur_mineral%name) // &
                                      'PREFACTOR RATE UNITS'
                      call InputDefaultMsg(input,option)
                    else
                      prefactor%rate = prefactor%rate * &
                                       UnitsConvertToInternal(word,option)
                    endif
                  case('ACTIVATION_ENERGY')
                    ! read activation energy for Arrhenius law
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
 
  cur_mineral => mineral%mineral_list
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
!geh  reaction%kinmnrl_names(imnrl) = cur_mineral%name
    endif
    cur_mineral => cur_mineral%next
  enddo
  
  cur_mineral => mineral%mineral_list
  do 
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo

end subroutine MineralReadKinetics

! ************************************************************************** !
!
! MineralReadFromDatabase: Reads mineral from database
! author: Glenn Hammond
! date: 10/16/08
!
! ************************************************************************** !
subroutine MineralReadFromDatabase(mineral,num_dbase_temperatures,input, &
                                   option)

  use Input_module
  use String_module  
  use Option_module
  use Database_Aux_module
  
  implicit none
  
  type(mineral_rxn_type) :: mineral
  PetscInt :: num_dbase_temperatures
  type(input_type) :: input
  type(option_type) :: option
  
  PetscInt :: ispec
  PetscInt :: itemp

  ! read the molar volume
  call InputReadDouble(input,option,mineral%molar_volume)
  call InputErrorMsg(input,option,'MINERAL molar volume','DATABASE')            
  ! convert from cm^3/mol to m^3/mol
  mineral%molar_volume = mineral%molar_volume*1.d-6
  ! create mineral reaction
  if (.not.associated(mineral%tstrxn)) then
    mineral%tstrxn => TransitionStateTheoryRxnCreate()
  endif
  ! read the number of aqueous species in mineral rxn
  mineral%dbaserxn => DatabaseRxnCreate()
  call InputReadInt(input,option,mineral%dbaserxn%nspec)
  call InputErrorMsg(input,option,'Number of species in mineral reaction', &
                  'DATABASE')  
  ! allocate arrays for rxn
  allocate(mineral%dbaserxn%spec_name(mineral%dbaserxn%nspec))
  mineral%dbaserxn%spec_name = ''
  allocate(mineral%dbaserxn%stoich(mineral%dbaserxn%nspec))
  mineral%dbaserxn%stoich = 0.d0
  allocate(mineral%dbaserxn%logK(num_dbase_temperatures))
  mineral%dbaserxn%logK = 0.d0
  ! read in species and stoichiometries
  do ispec = 1, mineral%dbaserxn%nspec
    call InputReadDouble(input,option,mineral%dbaserxn%stoich(ispec))
    call InputErrorMsg(input,option,'MINERAL species stoichiometry','DATABASE')            
    call InputReadQuotedWord(input,option,mineral%dbaserxn% &
                              spec_name(ispec),PETSC_TRUE)
    call InputErrorMsg(input,option,'MINERAL species name','DATABASE')            
  enddo
  do itemp = 1, num_dbase_temperatures
    call InputReadDouble(input,option,mineral%dbaserxn%logK(itemp))
    call InputErrorMsg(input,option,'MINERAL logKs','DATABASE')            
  enddo
  ! read the molar weight
  call InputReadDouble(input,option,mineral%molar_weight)
  call InputErrorMsg(input,option,'MINERAL molar weight','DATABASE')
        
end subroutine MineralReadFromDatabase

! ************************************************************************** !
!
! MineralProcessConstraint: Initializes constraints based on mineral
!                           species in system
! author: Glenn Hammond
! date: 01/07/13
!
! ************************************************************************** !
subroutine MineralProcessConstraint(mineral,constraint_name,constraint,option)

  use Option_module
  use Input_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(mineral_type), pointer :: mineral
  character(len=MAXWORDLENGTH) :: constraint_name
  type(mineral_constraint_type), pointer :: constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: imnrl, jmnrl
  
  character(len=MAXWORDLENGTH) :: mineral_name(mineral%nkinmnrl)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(mineral%nkinmnrl)
  PetscReal :: constraint_vol_frac(mineral%nkinmnrl)
  PetscReal :: constraint_area(mineral%nkinmnrl)
  PetscBool :: external_dataset(mineral%nkinmnrl)
  
  if (.not.associated(constraint)) return

  mineral_name = ''
  constraint_aux_string = ''
  external_dataset = PETSC_FALSE
  do imnrl = 1, mineral%nkinmnrl
    found = PETSC_FALSE
    do jmnrl = 1, mineral%nkinmnrl
      if (StringCompare(constraint%names(imnrl), &
                        mineral%kinmnrl_names(jmnrl), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
                'Mineral ' // trim(constraint%names(imnrl)) // &
                'from CONSTRAINT ' // trim(constraint_name) // &
                ' not found among kinetic minerals.'
      call printErrMsg(option)
    else
      constraint_vol_frac(jmnrl) = &
        constraint%constraint_vol_frac(imnrl)
      constraint_area(jmnrl) = &
        constraint%constraint_area(imnrl)
      mineral_name(jmnrl) = constraint%names(imnrl)
      constraint_aux_string(jmnrl) = constraint%constraint_aux_string(imnrl)
      external_dataset(jmnrl) = constraint%external_dataset(imnrl)
    endif  
  enddo
  constraint%names = mineral_name
  constraint%constraint_vol_frac = constraint_vol_frac
  constraint%constraint_area = constraint_area
  constraint%constraint_aux_string = constraint_aux_string
  constraint%external_dataset = external_dataset
  
end subroutine MineralProcessConstraint

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
#ifdef SOLID_SOLUTION
  use Solid_Solution_Aux_module
#endif
  
  implicit none
  
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

  PetscReal :: ln_spec_act, spec_act_coef
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
  
#ifdef SOLID_SOLUTION
  PetscBool :: cycle_
  PetscReal :: rate_scale(reaction%mineral%nkinmnrl)
  type(solid_solution_type), pointer :: cur_solid_soln
  PetscInt :: istoich_solid
#endif  

  type(mineral_type), pointer :: mineral

  PetscInt, parameter :: needs_to_be_fixed = 1
  
  PetscReal :: arrhenius_factor, rgas = 8.3144621d-3

  iphase = 1  
  mineral => reaction%mineral

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  if (reaction%neqcplx > 0) then
    ln_sec = log(rt_auxvar%sec_molal)
    ln_sec_act = ln_sec+log(rt_auxvar%sec_act_coef)
  endif

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    if (.not.reaction%use_geothermal_hpt)then
      call ReactionInterpolateLogK(mineral%kinmnrl_logKcoef, &
                                   mineral%kinmnrl_logK, &
                                   global_auxvar%temp(iphase), &
                                   mineral%nkinmnrl)
    else
     ! print *,'RKineticMineral:: ', global_auxvar%pres(iphase)
      call ReactionInterpolateLogK_hpt(mineral%kinmnrl_logKcoef, &
                                       mineral%kinmnrl_logK, &
                                       global_auxvar%temp(iphase), &
                                       global_auxvar%pres(iphase),&
                                       mineral%nkinmnrl)
    endif
  endif
#endif

#ifdef SOLID_SOLUTION
  rate_scale = 1.d0
  if (associated(reaction%solid_solution_list)) then
    do imnrl = 1, mineral%nkinmnrl ! for each mineral
      call RMineralRate(imnrl,ln_act,ln_sec_act,rt_auxvar,global_auxvar, &
                        QK,Im,Im_const,sum_prefactor_rate,affinity_factor, &
                        prefactor,ln_prefactor_spec,cycle_, &
                        reaction,mineral,option)
    enddo
  
    cur_solid_soln => reaction%solid_solution_list
    do
      if (.not.associated(cur_solid_soln)) exit
      tempreal = 0.d0
      do istoich_solid = 1, cur_solid_soln%num_stoich_solid
        imnrl = cur_solid_soln%stoich_solid_ids(istoich_solid)
        ! do something with mineral ikinmnrl rate, e.g. 
        tempreal = tempreal + rt_auxvar%mnrl_rate(imnrl)
        !tempreal = max(tempreal,dabs(rt_auxvar%mnrl_rate(ikinmnrl))
      enddo
      tempreal = tempreal / dble(cur_solid_soln%num_stoich_solid)
      do istoich_solid = 1, cur_solid_soln%num_stoich_solid
        imnrl = cur_solid_soln%stoich_solid_ids(istoich_solid)
        rate_scale(imnrl) = tempreal
      enddo
      cur_solid_soln => cur_solid_soln%next
    enddo
  endif
#endif

  ! Zero all rates as default 
  rt_auxvar%mnrl_rate(:) = 0.d0

  do imnrl = 1, mineral%nkinmnrl ! for each mineral

#ifdef SOLID_SOLUTION
    call RMineralRate(imnrl,ln_act,ln_sec_act,rt_auxvar,global_auxvar, &
                      QK,Im,Im_const,sum_prefactor_rate,affinity_factor, &
                      prefactor,ln_prefactor_spec,cycle_, &
                      reaction,mineral,option)
    if (cycle_) cycle
    
    Im = rate_scale(imnrl)*Im
    Im_const = rate_scale(imnrl)*Im_const
#else
    ! compute ion activity product
    lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN

    ! activity of water
    if (mineral%kinmnrlh2oid(imnrl) > 0) then
      lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)* &
                    rt_auxvar%ln_act_h2o
    endif

    ncomp = mineral%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = mineral%kinmnrlspecid(i,imnrl)
      lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
    enddo
    
    if (lnQK <= 6.90776d0) then
      QK = exp(lnQK)
    else
      QK = 1.d3
    endif
    
    if (associated(mineral%kinmnrl_Tempkin_const)) then
      affinity_factor = 1.d0-QK**(1.d0/ &
                                 mineral%kinmnrl_Tempkin_const(imnrl))
    else
      affinity_factor = 1.d0-QK
    endif
    
    sign_ = sign(1.d0,affinity_factor)

    if (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then

!   if ((mineral%kinmnrl_irreversible(imnrl) == 0 &
!     .and. (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0)) &
!     .or. (mineral%kinmnrl_irreversible(imnrl) == 1 .and. sign_ < 0.d0)) then
    
!     check for supersaturation threshold for precipitation
!     if (associated(mineral%kinmnrl_affinity_threshold)) then
      if (mineral%kinmnrl_affinity_threshold(imnrl) > 0.d0) then
        if (sign_ < 0.d0 .and. &
            QK < mineral%kinmnrl_affinity_threshold(imnrl)) cycle
      endif
    
!     check for rate limiter for precipitation
      if (mineral%kinmnrl_rate_limiter(imnrl) > 0.d0) then
        affinity_factor = affinity_factor/(1.d0+(1.d0-affinity_factor) &
          /mineral%kinmnrl_rate_limiter(imnrl))
      endif

      ! compute prefactor
      if (mineral%kinmnrl_num_prefactors(imnrl) > 0) then
        sum_prefactor_rate = 0
        prefactor = 0.d0
        ln_prefactor_spec = 0.d0
        ! sum over parallel prefactors
        do ipref = 1, mineral%kinmnrl_num_prefactors(imnrl)
          ln_prefactor = 0.d0
          ! product of "monod" equations
          do ipref_species = 1, mineral%kinmnrl_prefactor_id(0,ipref,imnrl)
            icomp = mineral%kinmnrl_prefactor_id(ipref_species,ipref,imnrl)
            if (icomp > 0) then ! primary species
              ln_spec_act = ln_act(icomp)
            else ! secondary species (given a negative id to differentiate)
              ln_spec_act = ln_sec_act(-icomp)
            endif
            ln_numerator = &
              mineral%kinmnrl_pref_alpha(ipref_species,ipref,imnrl)* &
              ln_spec_act
            ln_denominator = log(1.d0 + &
              exp(log(mineral%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl)) + &
                  mineral%kinmnrl_pref_beta(ipref_species,ipref,imnrl)* &
                  ln_spec_act))
            ln_prefactor = ln_prefactor + ln_numerator
            ln_prefactor = ln_prefactor - ln_denominator
            ln_prefactor_spec(ipref_species,ipref) = ln_numerator - ln_denominator
          enddo
          prefactor(ipref) = exp(ln_prefactor)
        ! Arrhenius factor
          arrhenius_factor = 1.d0
          if (mineral%kinmnrl_pref_activation_energy(ipref,imnrl) > 0.d0) then
            arrhenius_factor = &
              exp(mineral%kinmnrl_pref_activation_energy(ipref,imnrl)/rgas &
                  *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+ &
                                                273.15d0)))
          endif
          sum_prefactor_rate = sum_prefactor_rate + prefactor(ipref)* &
                               mineral%kinmnrl_pref_rate(ipref,imnrl)* &
                               arrhenius_factor
        enddo
      else
        ! Arrhenius factor
        arrhenius_factor = 1.d0
        if (mineral%kinmnrl_activation_energy(imnrl) > 0.d0) then
          arrhenius_factor = exp(mineral%kinmnrl_activation_energy(imnrl)/rgas &
            *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+273.15d0)))
        endif
        sum_prefactor_rate = mineral%kinmnrl_rate(imnrl)*arrhenius_factor
      endif

      ! compute rate
      ! rate: mol/m^2 mnrl/sec
      ! area: m^2 mnrl/m^3 bulk
      ! volume: m^3 bulk
      Im_const = -rt_auxvar%mnrl_area(imnrl)

      ! units: mol/sec/m^3 bulk
      if (associated(mineral%kinmnrl_affinity_power)) then
        ! Im_const: m^2 mnrl/m^3 bulk
        ! sum_prefactor_rate: mol/m^2 mnrl/sec
        Im = Im_const*sign_* &
             abs(affinity_factor)**mineral%kinmnrl_affinity_power(imnrl)* &
             sum_prefactor_rate
      else
        Im = Im_const*sign_*abs(affinity_factor)*sum_prefactor_rate
      endif
      ! store volumetric rate to be used in updating mineral volume fractions
      ! at end of time step
      rt_auxvar%mnrl_rate(imnrl) = Im ! mol/sec/m^3 bulk
    else ! rate is already zero by default; move on to next mineral
      cycle
    endif
#endif

    ! scale Im_const by volume for calculating derivatives below
    ! units: m^2 mnrl
    Im_const = Im_const*volume

    ! convert rate from volumetric (mol/sec/m^3 bulk) to mol/sec
    ! units: mol/sec
    Im = Im*volume
    
    ncomp = mineral%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = mineral%kinmnrlspecid(i,imnrl)
      Res(icomp) = Res(icomp) + mineral%kinmnrlstoich(i,imnrl)*Im
    enddo 
    
    if (.not. compute_derivative) cycle   

    ! calculate derivatives of rate with respect to free
    ! units = mol/sec
    if (associated(mineral%kinmnrl_affinity_power)) then
      dIm_dQK = -Im*mineral%kinmnrl_affinity_power(imnrl)/abs(affinity_factor)
    else
      dIm_dQK = -Im_const*sum_prefactor_rate
    endif
    
    if (associated(mineral%kinmnrl_Tempkin_const)) then
      dIm_dQK = dIm_dQK*(1.d0/mineral%kinmnrl_Tempkin_const(imnrl))/QK
    endif
    
    ! derivatives with respect to primary species in reaction quotient
    if (mineral%kinmnrl_rate_limiter(imnrl) <= 0.d0) then
      do j = 1, ncomp
        jcomp = mineral%kinmnrlspecid(j,imnrl)
        ! unit = L water/mol
        dQK_dCj = mineral%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
        ! units = (L water/mol)*(kg water/m^3 water)*(m^3 water/1000 L water) = kg water/mol
        dQK_dmj = dQK_dCj*global_auxvar%den_kg(iphase)*1.d-3 ! the multiplication by density could be moved
                                     ! outside the loop
        do i = 1, ncomp
          icomp = mineral%kinmnrlspecid(i,imnrl)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                             mineral%kinmnrlstoich(i,imnrl)*dIm_dQK*dQK_dmj
        enddo
      enddo
      
    else

      den = 1.d0+(1.d0-affinity_factor)/mineral%kinmnrl_rate_limiter(imnrl)
      do j = 1, ncomp
        jcomp = mineral%kinmnrlspecid(j,imnrl)
        ! unit = L water/mol
        dQK_dCj = mineral%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
        ! units = (L water/mol)*(kg water/m^3 water)*(m^3 water/1000 L water) = kg water/mol
        dQK_dmj = dQK_dCj*global_auxvar%den_kg(iphase)*1.d-3 ! the multiplication by density could be moved
                                     ! outside the loop
        do i = 1, ncomp
          icomp = mineral%kinmnrlspecid(i,imnrl)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
            mineral%kinmnrlstoich(i,imnrl)*dIm_dQK  &
            *(1.d0 + QK/mineral%kinmnrl_rate_limiter(imnrl)/den)*dQK_dmj/den
        enddo
      enddo
    endif
    
    if (mineral%kinmnrl_num_prefactors(imnrl) > 0) then ! add contribution of derivative in prefactor - messy
#if 1      
      dIm_dsum_prefactor_rate = Im/sum_prefactor_rate
      ! summation over parallel reactions (prefactors)
      do ipref = 1, mineral%kinmnrl_num_prefactors(imnrl)
        arrhenius_factor = 1.d0
        if (mineral%kinmnrl_pref_activation_energy(ipref,imnrl) > 0.d0) then
          arrhenius_factor = &
            exp(mineral%kinmnrl_pref_activation_energy(ipref,imnrl)/rgas &
                *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+ &
                                              273.15d0)))
        endif
        ! prefactor() saved in residual calc above
        ln_prefactor = log(prefactor(ipref))
        ! product of "monod" equations
        do ipref_species = 1, mineral%kinmnrl_prefactor_id(0,ipref,imnrl)
          ! derivative of 54 with respect to a single "monod" equation
          ! ln_prefactor_spec(,) saved in residual calc above
          dprefactor_dprefactor_spec = ln_prefactor-ln_prefactor_spec(ipref_species,ipref)
          icomp = mineral%kinmnrl_prefactor_id(ipref_species,ipref,imnrl)
          if (icomp > 0) then ! primary species
            ln_spec_act = ln_act(icomp)
            spec_act_coef = rt_auxvar%pri_act_coef(icomp)
          else ! secondary species
            ln_spec_act = ln_sec_act(-icomp)
            spec_act_coef = rt_auxvar%sec_act_coef(-icomp)
          endif
          ! derivative of numerator in eq. 54 with respect to species activity
          dprefactor_spec_dspec_numerator = &
            mineral%kinmnrl_pref_alpha(ipref_species,ipref,imnrl) * &
            exp(ln_prefactor_spec(ipref_species,ipref) - ln_spec_act)
          ln_gam_m_beta = mineral%kinmnrl_pref_beta(ipref_species,ipref,imnrl)* &
                          ln_spec_act
          ! denominator
          denominator = 1.d0 + &
              exp(log(mineral%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl)) + &
                  ln_gam_m_beta)
          ! derivative of denominator in eq. 54 with respect to species activity
          dprefactor_spec_dspec_denominator = -1.d0 * &
            exp(ln_prefactor_spec(ipref_species,ipref)) / denominator * &
            mineral%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl) * &
            mineral%kinmnrl_pref_beta(ipref_species,ipref,imnrl) * &
            exp(ln_gam_m_beta - ln_spec_act)

          ! chain rule for derivative of "monod" equation
          dprefactor_spec_dspec = dprefactor_spec_dspec_numerator + &
            dprefactor_spec_dspec_denominator

          ! thus far the derivative is with respect to the activity, convert to with
          ! respect to molality
          dprefactor_spec_dspec = dprefactor_spec_dspec * spec_act_coef

          dIm_dspec = dIm_dsum_prefactor_rate * dprefactor_dprefactor_spec * &
                      dprefactor_spec_dspec * &
                      mineral%kinmnrl_pref_rate(ipref,imnrl)* &
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
! RMineralRate: Calculates the mineral saturation index
! author: Glenn Hammond
! date: 08/29/11
!
! ************************************************************************** !
subroutine RMineralRate(imnrl,ln_act,ln_sec_act,rt_auxvar,global_auxvar, &
                        QK,Im,Im_const,sum_prefactor_rate,affinity_factor, &
                        prefactor,ln_prefactor_spec,cycle_, &
                        reaction,mineral,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(mineral_type) :: mineral
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: ln_sec_act(reaction%neqcplx)
  PetscReal :: QK
  PetscReal :: Im, Im_const
  PetscReal :: sum_prefactor_rate
  PetscReal :: affinity_factor
  PetscReal :: prefactor(10), ln_prefactor_spec(5,10)
  PetscBool :: cycle_
  
  PetscReal :: lnQK
  PetscInt :: i, imnrl, icomp, ncomp, ipref, ipref_species
  PetscReal :: sign_

  PetscReal :: ln_spec_act
  PetscReal :: ln_prefactor, ln_numerator, ln_denominator
  
  PetscReal :: arrhenius_factor
  PetscReal, parameter :: rgas = 8.3144621d-3
  PetscInt, parameter :: iphase = 1
  
  cycle_ = PETSC_FALSE
  
  ! compute ion activity product
  lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN

  ! activity of water
  if (mineral%kinmnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)* &
                  rt_auxvar%ln_act_h2o
  endif

  ncomp = mineral%kinmnrlspecid(0,imnrl)
  do i = 1, ncomp
    icomp = mineral%kinmnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
  enddo
    
  if (lnQK <= 6.90776d0) then
    QK = exp(lnQK)
  else
    QK = 1.d3
  endif
    
  if (associated(mineral%kinmnrl_Tempkin_const)) then
    affinity_factor = 1.d0-QK**(1.d0/ &
                                mineral%kinmnrl_Tempkin_const(imnrl))
  else
    affinity_factor = 1.d0-QK
  endif
    
  sign_ = sign(1.d0,affinity_factor)

  if (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then

!   if ((mineral%kinmnrl_irreversible(imnrl) == 0 &
!     .and. (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0)) &
!     .or. (mineral%kinmnrl_irreversible(imnrl) == 1 .and. sign_ < 0.d0)) then
    
!     check for supersaturation threshold for precipitation
!     if (associated(mineral%kinmnrl_affinity_threshold)) then
    if (mineral%kinmnrl_affinity_threshold(imnrl) > 0.d0) then
      if (sign_ < 0.d0 .and. &
          QK < mineral%kinmnrl_affinity_threshold(imnrl)) then
        cycle_ = PETSC_TRUE
        return
      endif
    endif
    
!     check for rate limiter for precipitation
    if (mineral%kinmnrl_rate_limiter(imnrl) > 0.d0) then
      affinity_factor = affinity_factor/(1.d0+(1.d0-affinity_factor) &
        /mineral%kinmnrl_rate_limiter(imnrl))
    endif

    ! compute prefactor
    if (mineral%kinmnrl_num_prefactors(imnrl) > 0) then
      sum_prefactor_rate = 0
      prefactor = 0.d0
      ln_prefactor_spec = 0.d0
      ! sum over parallel prefactors
      do ipref = 1, mineral%kinmnrl_num_prefactors(imnrl)
        ln_prefactor = 0.d0
        ! product of "monod" equations
        do ipref_species = 1, mineral%kinmnrl_prefactor_id(0,ipref,imnrl)
          icomp = mineral%kinmnrl_prefactor_id(ipref_species,ipref,imnrl)
          if (icomp > 0) then ! primary species
            ln_spec_act = ln_act(icomp)
          else ! secondary species (given a negative id to differentiate)
            ln_spec_act = ln_sec_act(-icomp)
          endif
          ln_numerator = &
            mineral%kinmnrl_pref_alpha(ipref_species,ipref,imnrl)* &
            ln_spec_act
          ln_denominator = log(1.d0 + &
            exp(log(mineral%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl)) + &
                mineral%kinmnrl_pref_beta(ipref_species,ipref,imnrl)* &
                ln_spec_act))
          ln_prefactor = ln_prefactor + ln_numerator
          ln_prefactor = ln_prefactor - ln_denominator
          ln_prefactor_spec(ipref_species,ipref) = ln_numerator - ln_denominator
        enddo
        prefactor(ipref) = exp(ln_prefactor)
      ! Arrhenius factor
        arrhenius_factor = 1.d0
        if (mineral%kinmnrl_pref_activation_energy(ipref,imnrl) > 0.d0) then
          arrhenius_factor = &
            exp(mineral%kinmnrl_pref_activation_energy(ipref,imnrl)/rgas &
                *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+ &
                                              273.15d0)))
        endif
        sum_prefactor_rate = sum_prefactor_rate + prefactor(ipref)* &
                              mineral%kinmnrl_pref_rate(ipref,imnrl)* &
                              arrhenius_factor
      enddo
    else
      ! Arrhenius factor
      arrhenius_factor = 1.d0
      if (mineral%kinmnrl_activation_energy(imnrl) > 0.d0) then
        arrhenius_factor = exp(mineral%kinmnrl_activation_energy(imnrl)/rgas &
          *(1.d0/(25.d0+273.15d0)-1.d0/(global_auxvar%temp(iphase)+273.15d0)))
      endif
      sum_prefactor_rate = mineral%kinmnrl_rate(imnrl)*arrhenius_factor
    endif

    ! compute rate
    ! rate: mol/m^2 mnrl/sec
    ! area: m^2 mnrl/m^3 bulk
    ! volume: m^3 bulk
    Im_const = -rt_auxvar%mnrl_area(imnrl)

    ! units: mol/sec/m^3 bulk
    if (associated(mineral%kinmnrl_affinity_power)) then
      ! Im_const: m^2 mnrl/m^3 bulk
      ! sum_prefactor_rate: mol/m^2 mnrl/sec
      Im = Im_const*sign_* &
            abs(affinity_factor)**mineral%kinmnrl_affinity_power(imnrl)* &
            sum_prefactor_rate
    else
      Im = Im_const*sign_*abs(affinity_factor)*sum_prefactor_rate
    endif
    ! store volumetric rate to be used in updating mineral volume fractions
    ! at end of time step
    rt_auxvar%mnrl_rate(imnrl) = Im ! mol/sec/m^3 bulk
  else ! rate is already zero by default; move on to next mineral
    cycle_ = PETSC_TRUE
  endif

end subroutine RMineralRate

! ************************************************************************** !
!
! RMineralSaturationIndex: Calculates the mineral saturation index
! author: Glenn Hammond
! date: 08/29/11
!
! ************************************************************************** !
function RMineralSaturationIndex(imnrl,rt_auxvar,global_auxvar,reaction,option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: imnrl
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscReal :: RMineralSaturationIndex
  PetscInt :: i, icomp
  PetscReal :: lnQK
  PetscInt, parameter :: iphase = 1
  type(mineral_type), pointer :: mineral
  
  mineral => reaction%mineral

#ifdef TEMP_DEPENDENT_LOGK
  if (.not.option%use_isothermal) then
    if (.not.reaction%use_geothermal_hpt)then
      call ReactionInterpolateLogK(mineral%mnrl_logKcoef, &
                                   mineral%mnrl_logK, &
                                   global_auxvar%temp(iphase), &
                                   mineral%nmnrl)
    else
    !  print *,'RMineralSaturationIndex:: ',global_auxvar%pres(iphase)
      call ReactionInterpolateLogK_hpt(mineral%mnrl_logKcoef, &
                                       mineral%mnrl_logK, &
                                       global_auxvar%temp(iphase), &
                                       global_auxvar%pres(iphase),&
                                       mineral%nmnrl)
    endif
  endif
#endif  

  ! compute saturation
  lnQK = -mineral%mnrl_logK(imnrl)*LOG_TO_LN
  if (mineral%mnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%mnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
  endif
  do i = 1, mineral%mnrlspecid(0,imnrl)
    icomp = mineral%mnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%mnrlstoich(i,imnrl)* &
           log(rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp))
  enddo
  RMineralSaturationIndex = exp(lnQK)    

end function RMineralSaturationIndex

end module Mineral_module
