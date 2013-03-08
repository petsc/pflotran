module Saturation_Function_module
 

  implicit none

  private

#include "definitions.h"
 
  type, public :: saturation_function_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: saturation_function_ctype
    PetscInt :: saturation_function_itype
    character(len=MAXWORDLENGTH) :: permeability_function_ctype
    PetscInt :: permeability_function_itype
    PetscReal, pointer :: Sr(:)
    PetscReal :: m
    PetscReal :: lambda
    PetscReal :: alpha
    PetscReal :: pcwmax
    PetscReal :: betac
    PetscReal :: power
    PetscInt :: hysteresis_id
    PetscInt :: hysteresis_params(6)
    PetscReal :: spline_low
    PetscReal :: spline_high
    PetscReal :: spline_coefficients(4)
    PetscReal :: ani_A       ! parameters for anisotropic relative permeability model
    PetscReal :: ani_B       ! parameters for anisotropic relative permeability model
    PetscReal :: ani_C       ! parameters for anisotropic relative permeability model

    type(saturation_function_type), pointer :: next
  end type saturation_function_type
  
  type, public :: saturation_function_ptr_type
    type(saturation_function_type), pointer :: ptr
  end type saturation_function_ptr_type

  interface SaturationFunctionCompute
    module procedure SaturationFunctionCompute1
    module procedure SaturationFunctionCompute2
    module procedure SaturationFunctionCompute3
  end interface
  
  public :: SaturationFunctionCreate, &
            SaturationFunctionDestroy, &
            SaturationFunctionAddToList, &
            SaturationFunctionCompute, &
            SaturatFuncConvertListToArray, &
            SaturationFunctionComputeSpline, &
            PermFunctionComputeSpline, &
            SaturationFunctionRead, &
            SatFuncGetRelPermFromSat, &
            SatFuncGetCapillaryPressure, &
            SaturationFunctionGetID, &
            SaturationFunctionComputeIce, &
            CapillaryPressureThreshold

  ! Saturation function 
  PetscInt, parameter :: VAN_GENUCHTEN = 1
  PetscInt, parameter :: BROOKS_COREY = 2
  PetscInt, parameter :: THOMEER_COREY = 3
  PetscInt, parameter :: NMT_EXP = 4
  PetscInt, parameter :: PRUESS_1 = 5

  ! Permeability function
  PetscInt, parameter :: DEFAULT = 0
  PetscInt, parameter :: BURDINE = 1
  PetscInt, parameter :: MUALEM = 2
  
contains

! ************************************************************************** !
!
! SaturationFunctionCreate: Creates a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function SaturationFunctionCreate(option)
  
  use Option_module
  
  implicit none

  type(saturation_function_type), pointer :: SaturationFunctionCreate
  type(option_type) :: option
  
  type(saturation_function_type), pointer :: saturation_function
  
  allocate(saturation_function)
  saturation_function%id = 0
  saturation_function%name = ''
  saturation_function%saturation_function_ctype = 'VAN_GENUCHTEN'
  saturation_function%saturation_function_itype = VAN_GENUCHTEN
  saturation_function%permeability_function_ctype = 'MUALEM'
  saturation_function%permeability_function_itype = MUALEM
  allocate(saturation_function%Sr(option%nphase))
  saturation_function%Sr = 0.d0
  saturation_function%m = 0.d0
  saturation_function%lambda = 0.d0
  saturation_function%alpha = 0.d0
  saturation_function%pcwmax = 0.d0
  saturation_function%betac = 0.d0
  saturation_function%power = 0.d0
  saturation_function%hysteresis_id = 0
  saturation_function%hysteresis_params = 0
  saturation_function%spline_low = 0.d0
  saturation_function%spline_high = 0.d0
  saturation_function%spline_coefficients = 0.d0
  saturation_function%ani_A = 0.d0
  saturation_function%ani_B = 0.d0
  saturation_function%ani_C = 0.d0
  nullify(saturation_function%next)
  SaturationFunctionCreate => saturation_function

end function SaturationFunctionCreate

! ************************************************************************** !
!
! SaturationFunctionRead: Reads in contents of a saturation_function card
! author: Glenn Hammond
! date: 01/21/09
! 
! ************************************************************************** !
subroutine SaturationFunctionRead(saturation_function,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  type(saturation_function_type) :: saturation_function
  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: iphase
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('PERMEABILITY_FUNCTION_TYPE') 
        call InputReadWord(input,option, &
                           saturation_function%permeability_function_ctype, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'permeability function type', &
                           'SATURATION_FUNCTION')
      case('SATURATION_FUNCTION_TYPE') 
        call InputReadWord(input,option, &
                           saturation_function%saturation_function_ctype, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'saturation function type', &
                           'SATURATION_FUNCTION')
      case('RESIDUAL_SATURATION') 
        select case(option%iflowmode)
          case(FLASH2_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
            end select
            call InputReadDouble(input,option,saturation_function%Sr(iphase))
            word = trim(keyword) // ' residual saturation'
            call InputErrorMsg(input,option,word,'SATURATION_FUNCTION')
          case(MPH_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
            end select
            call InputReadDouble(input,option,saturation_function%Sr(iphase))
            word = trim(keyword) // ' residual saturation'
            call InputErrorMsg(input,option,word,'SATURATION_FUNCTION')
          case(IMS_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
              case('OIL','OIL_PHASE','NAPL','NAPL_PHASE')
                iphase = 3
            end select
            call InputReadDouble(input,option,saturation_function%Sr(iphase))
            word = trim(keyword) // ' residual saturation'
            call InputErrorMsg(input,option,word,'SATURATION_FUNCTION')
          case(RICHARDS_MODE,TH_MODE,THC_MODE,THMC_MODE,G_MODE)
            call InputReadDouble(input,option,saturation_function%Sr(1))
            call InputErrorMsg(input,option,'residual saturation','SATURATION_FUNCTION')
        end select
      case('LAMBDA') 
        call InputReadDouble(input,option,saturation_function%lambda)
        call InputErrorMsg(input,option,'residual saturation','SATURATION_FUNCTION')
        saturation_function%m = saturation_function%lambda
      case('ALPHA') 
        call InputReadDouble(input,option,saturation_function%alpha)
        call InputErrorMsg(input,option,'alpha','SATURATION_FUNCTION')
      case('MAX_CAPILLARY_PRESSURE') 
        call InputReadDouble(input,option,saturation_function%pcwmax)
        call InputErrorMsg(input,option,'maximum capillary pressure','SATURATION_FUNCTION')
      case('BETAC') 
        call InputReadDouble(input,option,saturation_function%betac)
        call InputErrorMsg(input,option,'betac','SATURATION_FUNCTION')
      case('POWER') 
        call InputReadDouble(input,option,saturation_function%power)
        call InputErrorMsg(input,option,'power','SATURATION_FUNCTION')
      case('ANI_A') 
        call InputReadDouble(input,option,saturation_function%ani_A)
        call InputErrorMsg(input,option,'ani_A','SATURATION_FUNCTION')
      case('ANI_B') 
        call InputReadDouble(input,option,saturation_function%ani_B)
        call InputErrorMsg(input,option,'ani_B','SATURATION_FUNCTION')
      case('ANI_C') 
        call InputReadDouble(input,option,saturation_function%ani_C)
        call InputErrorMsg(input,option,'ani_C','SATURATION_FUNCTION')
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in saturation_function'    
        call printErrMsg(option)
    end select 
  
  enddo 
  
  if (saturation_function%m < 1.d-40 .and. &
      .not.StringCompare(saturation_function%name,'default',SEVEN_INTEGER)) then
    option%io_buffer = 'Saturation function parameter "m" not set ' // &
                       'properly in saturation function "' // &
                       trim(saturation_function%name) // '".'
    call printErrMsg(option)
  endif 

end subroutine SaturationFunctionRead

! ************************************************************************** !
!
! SaturationFunctionComputeSpline: Computes a spline spanning the 
!                                  discontinuity in Brooks Corey
! author: Glenn Hammond
! date: 02/27/08
!
! ************************************************************************** !
subroutine SaturationFunctionComputeSpline(option,saturation_function)
  
  use Option_module
  use Utility_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  
  PetscReal :: b(4)
  PetscReal :: pressure_high, pressure_low
  
  PetscReal :: n

  select case(saturation_function%saturation_function_itype) 
    case(BROOKS_COREY)

      ! fill matix with values
      pressure_high = 1.d0/saturation_function%alpha*2.d0
      pressure_low = 1.d0/saturation_function%alpha*0.5d0
  
      saturation_function%spline_low = pressure_low
      saturation_function%spline_high = pressure_high
  
      b(1) = (pressure_high*saturation_function%alpha)** &
               (-saturation_function%lambda)
      b(2) = 1.d0
      b(3) = -saturation_function%lambda/pressure_high* &
               (pressure_high*saturation_function%alpha)** &
                 (-saturation_function%lambda)
      b(4) = 0.d0

      call CubicPolynomialSetup(pressure_high,pressure_low,b)
      
      saturation_function%spline_coefficients(1:4) = b(1:4)
      
  case(VAN_GENUCHTEN)
 
      ! return for now
      return

      !geh: keep for now
#if 0
      ! fill matix with values
      ! these are capillary pressures
      pressure_low = 0  ! saturated
      pressure_high = 0.01d0*option%reference_pressure  ! just below saturated
  
      saturation_function%spline_low = pressure_low
      saturation_function%spline_high = pressure_high
    
      n = 1.d0/(1.d0 - saturation_function%m)
      b(1) = (1.d0+(pressure_high*saturation_function%alpha)**n)** &
               (-saturation_function%m)
      b(2) = 1.d0
      b(3) = -saturation_function%m*n*saturation_function%alpha* &
             (saturation_function%alpha*pressure_high)**(n-1.d0)* &
             (1.d0+(saturation_function%alpha*pressure_high)**n)** &
               (-saturation_function%m-1.d0)
      b(4) = 0.d0
  
      call CubicPolynomialSetup(pressure_high,pressure_low,b)

      saturation_function%spline_coefficients(1:4) = b(1:4)
#endif

  end select
  
end subroutine SaturationFunctionComputeSpline

! ************************************************************************** !
!
! PermFunctionComputeSpline: Computes a spline spanning the 
!                                  discontinuity in Brooks Corey
! author: Glenn Hammond
! date: 02/27/12
!
! ************************************************************************** !
subroutine PermFunctionComputeSpline(option,saturation_function)
  
  use Option_module
  use Utility_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  
  PetscReal :: b(4)
  PetscReal :: Se_high, Se_low, one_over_m, Se_one_over_m, m
  
  select case(saturation_function%saturation_function_itype) 

    case(BROOKS_COREY)

    case(VAN_GENUCHTEN)
 
#ifdef MUALEM_SPLINE
      ! fill matix with values
      Se_low = 0.99d0  ! saturated
      Se_high = 1.d0  ! just below saturated
  
      saturation_function%spline_low = Se_low
      saturation_function%spline_high = Se_high
    
      m = saturation_function%m
      one_over_m = 1.d0/m
      Se_one_over_m = Se_low**one_over_m
      b(1) = 1.d0
      b(2) = sqrt(Se_low)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
      b(3) = 0.d0
      b(4) = 0.5d0*b(2)/Se_low+ &
             2.d0*Se_low**(one_over_m-0.5d0)* &
             (1.d0-Se_one_over_m)**(m-1.d0)* &
             (1.d0-(1.d0-Se_one_over_m)**m)
  
      call CubicPolynomialSetup(Se_high,Se_low,b)
  
      saturation_function%spline_coefficients(1:4) = b(1:4)
#endif

  end select
  
end subroutine PermFunctionComputeSpline

! ************************************************************************** !
!
! SaturationFunctionAddToList: Adds a saturation function to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine SaturationFunctionAddToList(saturation_function,list)

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  type(saturation_function_type), pointer :: list

  type(saturation_function_type), pointer :: cur_saturation_function
  
  if (associated(list)) then
    cur_saturation_function => list
    ! loop to end of list
    do
      if (.not.associated(cur_saturation_function%next)) exit
      cur_saturation_function => cur_saturation_function%next
    enddo
    cur_saturation_function%next => saturation_function
  else
    list => saturation_function
  endif
  
end subroutine SaturationFunctionAddToList

! ************************************************************************** !
!
! SaturatFuncConvertListToArray: Creates an array of pointers to the 
!                                saturation functions in the list
! author: Glenn Hammond
! date: 12/11/07
!
! ************************************************************************** !
subroutine SaturatFuncConvertListToArray(list,array,option)

  use String_module
  use Option_module
  
  implicit none
  
  type(saturation_function_type), pointer :: list
  type(saturation_function_ptr_type), pointer :: array(:)
  type(option_type) :: option
    
  type(saturation_function_type), pointer :: cur_saturation_function
  PetscInt :: count

  count = 0
  cur_saturation_function => list
  do 
    if (.not.associated(cur_saturation_function)) exit
    count = count + 1
    
    ! set permeability function integer type
    call StringToUpper(cur_saturation_function%permeability_function_ctype)
    select case(trim(cur_saturation_function%permeability_function_ctype))
      case('DEFAULT')
        cur_saturation_function%permeability_function_itype = DEFAULT
      case('BURDINE')
        cur_saturation_function%permeability_function_itype = BURDINE
      case('MUALEM')
        cur_saturation_function%permeability_function_itype = MUALEM
      case('NMT_EXP')
        cur_saturation_function%permeability_function_itype = NMT_EXP
      case('PRUESS_1')
        cur_saturation_function%permeability_function_itype = PRUESS_1
      case default
        option%io_buffer = 'Permeability function type "' // &
                           trim(cur_saturation_function%permeability_function_ctype) // &
                           '" not recognized ' // &
                           ' in saturation function ' // &
                           trim(cur_saturation_function%name)
        call printErrMsg(option)
    end select
    
    ! set saturation function integer type
    call StringToUpper(cur_saturation_function%saturation_function_ctype)
    select case(trim(cur_saturation_function%saturation_function_ctype))
      case('VAN_GENUCHTEN')
        cur_saturation_function%saturation_function_itype = VAN_GENUCHTEN
      case('BROOKS_COREY')
        cur_saturation_function%saturation_function_itype = BROOKS_COREY
      case('THOMEER_COREY')
        cur_saturation_function%saturation_function_itype = THOMEER_COREY
      case('NMT_EXP')
        cur_saturation_function%saturation_function_itype = NMT_EXP
      case('PRUESS_1')
        cur_saturation_function%saturation_function_itype = PRUESS_1
       
      case default
        option%io_buffer = 'Saturation function type "' // &
                           trim(cur_saturation_function%saturation_function_ctype) // &
                           '" not recognized ' // &
                           ' in saturation function ' // &
                           trim(cur_saturation_function%name)
        call printErrMsg(option)
    end select
    
    cur_saturation_function => cur_saturation_function%next
  enddo
  
  allocate(array(count))
  
  count = 0
  cur_saturation_function => list
  do 
    if (.not.associated(cur_saturation_function)) exit
    count = count + 1
    cur_saturation_function%id = count
    array(count)%ptr => cur_saturation_function
    cur_saturation_function => cur_saturation_function%next
  enddo

end subroutine SaturatFuncConvertListToArray

! ************************************************************************** !
!
! SaturationFunctionCompute1: Computes the saturation and relative permeability
!                             (and associated derivatives) as a function of 
!                             capillary pressure
! author: Glenn Hammond
! date: 2/9/12
!
! ************************************************************************** !
subroutine SaturationFunctionCompute1(pressure,saturation,relative_perm, &
                                     dsat_pres,dkr_pres, &
                                     saturation_function, &
                                     auxvar1,auxvar2, &
                                     option)
  use Option_module
  
  implicit none

  PetscReal :: pressure, saturation, relative_perm, dsat_pres, dkr_pres
  type(saturation_function_type) :: saturation_function
  PetscReal :: auxvar1,auxvar2
  type(option_type) :: option

  PetscBool :: switch_to_saturated
  
  call SaturationFunctionCompute2(pressure,saturation,relative_perm, &
                                  dsat_pres,dkr_pres, &
                                  saturation_function, &
                                  auxvar1,auxvar2, &
                                  switch_to_saturated,option)

end subroutine SaturationFunctionCompute1

! ************************************************************************** !
!
! SaturationFunctionCompute2: Computes the saturation and relative permeability
!                             (and associated derivatives) as a function of 
!                             capillary pressure
! author: Glenn Hammond
! date: 12/11/07
!
! ************************************************************************** !
subroutine SaturationFunctionCompute2(pressure,saturation,relative_perm, &
                                      dsat_pres,dkr_pres, &
                                      saturation_function, &
                                      auxvar1,auxvar2, &
                                      switch_to_saturated,option)
  use Option_module
  use Utility_module, only:CubicPolynomialEvaluate
  
  implicit none

  PetscReal :: pressure, saturation, relative_perm, dsat_pres, dkr_pres
  type(saturation_function_type) :: saturation_function
  PetscReal :: auxvar1,auxvar2
  PetscBool :: switch_to_saturated
  type(option_type) :: option

  PetscInt :: iphase
  PetscReal :: alpha, lambda, m, n, Sr, one_over_alpha
  PetscReal :: pc, Se, one_over_m, Se_one_over_m, dSe_pc, dsat_pc, dkr_pc
  PetscReal :: dkr_Se, power
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: por, perm
  PetscReal :: Fg, a, Pd, PHg

  PetscReal, parameter :: pc_alpha_n_epsilon = 1.d-15
  
  iphase = 1
  dsat_pres = 0.d0
  dkr_pres = 0.d0
  
  ! compute saturation
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
#if 0
      if (pc < saturation_function%spline_low) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
      else if (pc < saturation_function%spline_high) then
        Sr = saturation_function%Sr(iphase)
        call CubicPolynomialEvaluate(saturation_function%spline_coefficients, &
                                     pc,Se,dSe_pc)
        saturation = Sr + (1.d0-Sr)*Se
        dsat_pc = (1.d0-Sr)*dSe_pc
#else
      if (pressure >= option%reference_pressure) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
#endif        
      else
        alpha = saturation_function%alpha
        pc = option%reference_pressure-pressure
        m = saturation_function%m
        n = 1.d0/(1.d0-m)
        pc_alpha = pc*alpha
        pc_alpha_n = pc_alpha**n
        !geh:  This conditional does not catch potential cancelation in 
        !      the dkr_Se deriviative calculation.  Therefore, I am setting
        !      an epsilon here
!        if (1.d0 + pc_alpha_n == 1.d0) then ! check for zero perturbation
        if (pc_alpha_n < pc_alpha_n_epsilon) then 
          saturation = 1.d0
          relative_perm = 1.d0
          switch_to_saturated = PETSC_TRUE
          return
        endif
        one_plus_pc_alpha_n = 1.d0+pc_alpha_n
        Sr = saturation_function%Sr(iphase)
        Se = one_plus_pc_alpha_n**(-m)
        dSe_pc = -m*n*alpha*pc_alpha_n/(pc_alpha*one_plus_pc_alpha_n**(m+1))
        saturation = Sr + (1.d0-Sr)*Se
        dsat_pc = (1.d0-Sr)*dSe_pc
      endif
      if (saturation > 1.d0) then
        print *, option%myrank, 'vG Saturation > 1:', saturation
      else if (saturation > 1.d0 .or. saturation < Sr) then
        print *, option%myrank, 'vG Saturation < Sr:', saturation, Sr
      endif
      ! compute relative permeability
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = Se*Se*(1.d0-(1.d0-Se_one_over_m)**m)
          dkr_Se = 2.d0*relative_perm/Se + &
                   Se*Se_one_over_m*(1.d0-Se_one_over_m)**(m-1.d0)
          dkr_pc = dkr_Se*dSe_pc
        case(MUALEM)
#ifdef MUALEM_SPLINE
          if (Se > saturation_function%spline_low) then
            call CubicPolynomialEvaluate( &
              saturation_function%spline_coefficients, &
              Se,relative_perm,dkr_Se)
          else
#endif          
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
          dkr_Se = 0.5d0*relative_perm/Se+ &
                   2.d0*Se**(one_over_m-0.5d0)* &
                        (1.d0-Se_one_over_m)**(m-1.d0)* &
                        (1.d0-(1.d0-Se_one_over_m)**m)
#ifdef MUALEM_SPLINE
          endif
#endif          
          dkr_pc = dkr_Se*dSe_pc
        case default
          option%io_buffer = 'Unknown relative permeabilty function' 
          call printErrMsg(option)
      end select
    case(BROOKS_COREY)
      alpha = saturation_function%alpha
      one_over_alpha = 1.d0/alpha
      pc = option%reference_pressure-pressure
#if 0
      if (pc < saturation_function%spline_low) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
      else if (pc < saturation_function%spline_high) then
        Sr = saturation_function%Sr(iphase)
        call CublicPolynomialEvaluate(saturation_function%spline_coefficients, &
                                      pc,Se,dSe_pc)
        saturation = Sr + (1.d0-Sr)*Se
        dsat_pc = (1.d0-Sr)*dSe_pc
#else
      if (pc < one_over_alpha) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
#endif        
      else
        lambda = saturation_function%lambda
        Sr = saturation_function%Sr(iphase)
        pc_alpha_neg_lambda = (pc*alpha)**(-lambda)
        Se = pc_alpha_neg_lambda
        dSe_pc = -lambda/pc*pc_alpha_neg_lambda
        saturation = Sr + (1.d0-Sr)*Se
!        dsat_pc = -lambda*(1.d0-Sr)/pc*pc_alpha_neg_lambda
        dsat_pc = (1.d0-Sr)*dSe_pc
      endif
      if (saturation > 1.d0) then
        print *, option%myrank, 'BC Saturation > 1:', saturation
      else if (saturation > 1.d0 .or. saturation < Sr) then
        print *, option%myrank, 'BC Saturation < Sr:', saturation, Sr
      endif
      ! compute relative permeability
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          power = 3.d0+2.d0/lambda
          relative_perm = Se**power
          dkr_Se = power*relative_perm/Se
          dkr_pc = dkr_Se*dSe_pc
        case(MUALEM)
          power = 2.5d0+2.d0/lambda
          relative_perm = Se**power
          dkr_Se = power*relative_perm/Se
          dkr_pc = dkr_Se*dSe_pc
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(THOMEER_COREY)
      pc = option%reference_pressure-pressure
      por = auxvar1
      perm = auxvar2*1.013202d15 ! convert from m^2 to mD
      Fg = saturation_function%alpha
      a = saturation_function%m
      Pd = 100.d0*por/sqrt(perm/(3.8068d0*(Fg**(-1.334d0)))) ! psi
      PHg = 9.63051d-4*pc
      if (PHg > Pd) then
        saturation = 1.d0-exp(-Fg/log10(PHg/Pd))
#if 0
        alpha = pc*(1.d0+1.d-8)
        m = 9.63051d-4*alpha
        n = 1.d0-exp(-Fg/log10(m/Pd))
        n = (n-saturation)/(alpha-pc)
#endif        
        dsat_pc = (saturation-1.d0)*Fg/(log10(PHg/Pd)**2.d0)/(pc*2.30258509d0)
        ! Sr assumed to be zero
        relative_perm = saturation**a
        dkr_pc = a*saturation**(a-1.d0)*dsat_pc
      else
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
      endif
    case default
      option%io_buffer = 'Unknown saturation function'
      call printErrMsg(option)
  end select

  dsat_pres = -dsat_pc 
  dkr_pres = -dkr_pc

end subroutine SaturationFunctionCompute2

! ************************************************************************** !
!
! SaturationFunctionCompute3: Computes just saturation as a function of 
!                             capillary pressure
! author: Glenn Hammond
! date: 2/9/12
!
! ************************************************************************** !
subroutine SaturationFunctionCompute3(pressure,saturation, &
                                      saturation_function, &
                                      option)
  use Option_module
  
  implicit none

  PetscReal :: pressure, saturation
  PetscReal :: relative_perm_dummy, dsat_pres_dummy, dkr_pres_dummy
  type(saturation_function_type) :: saturation_function
  PetscReal :: auxvar1_dummy, auxvar2_dummy
  type(option_type) :: option

  PetscBool :: switch_to_saturated_dummy
  
  call SaturationFunctionCompute2(pressure,saturation,relative_perm_dummy, &
                                  dsat_pres_dummy,dkr_pres_dummy, &
                                  saturation_function, &
                                  auxvar1_dummy,auxvar2_dummy, &
                                  switch_to_saturated_dummy,option)

end subroutine SaturationFunctionCompute3

! ************************************************************************** !
!
! SaturationFunctionComputeIce:Computes the saturation of ice, water vapor 
!                              and liquid water given the saturation function
!                              temperature, water vapor pressure and liquid
!                              water pressure 
! author: Satish Karra
! date: 11/14/11
!
! ************************************************************************** !
subroutine SaturationFunctionComputeIce(liquid_pressure, temperature, &
                                        ice_saturation, &
                                        liquid_saturation, gas_saturation, &
                                        liquid_relative_perm, dsl_pl, & 
                                        dsl_temp, dsg_pl, dsg_temp, dsi_pl, &
                                        dsi_temp, dkr_pl, dkr_temp, &
                                        saturation_function, pth, option)

  use Option_module
 
implicit none

  PetscReal :: liquid_pressure, temperature
  PetscReal :: ice_saturation, liquid_saturation, gas_saturation
  PetscReal :: liquid_relative_perm
  PetscReal :: dsl_pl, dsl_temp
  PetscReal :: dsg_pl, dsg_temp
  PetscReal :: dsi_pl, dsi_temp
  PetscReal :: dkr_pl
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option

  PetscReal :: alpha, lambda, m, n
  PetscReal :: pc, Se, one_over_m, Se_one_over_m, dSe_pc, dkr_pc
  PetscReal :: dkr_Se, power
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: function_A, function_B
  PetscReal :: pc_il, gamma, pc_il_alpha, pc_il_alpha_n, Se_temp
  PetscReal :: one_plus_pc_il_alpha_n
  PetscReal :: dfunc_A_temp
  PetscReal :: dfunc_B_pl
  PetscReal :: liq_sat_one_over_m, dkr_ds_liq, dkr_temp
  PetscReal :: pth, dSe_pc_at_pth
  PetscReal, parameter :: den_ice = 9.167d2 !in kg/m3 at 273.15K
  PetscReal, parameter :: heat_of_fusion = 3.34d5 !in J/kg at 273.15K
  PetscReal, parameter :: interfacial_tensions_ratio = 2.33
  PetscReal, parameter :: T_0 = 273.15d0 !in K
  
  dsl_pl = 0.d0
  dsl_temp = 0.d0
  dsg_pl = 0.d0
  dsg_temp = 0.d0
  dsi_pl = 0.d0
  dsi_temp = 0.d0
  dkr_pl = 0.d0
  dkr_temp = 0.d0
  dkr_ds_liq = 0.d0
  
  ! compute saturation
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      if (liquid_pressure >= option%reference_pressure) then
        function_B = 1.d0
        dfunc_B_pl = 0.d0
      else
        alpha = saturation_function%alpha
        pc = option%reference_pressure - liquid_pressure
        m = saturation_function%m
        n = 1.d0/(1.d0 - m)
        pc_alpha = pc*alpha
        pc_alpha_n = pc_alpha**n
        one_plus_pc_alpha_n = 1.d0 + pc_alpha_n
        Se = one_plus_pc_alpha_n**(-m)
        dSe_pc = -m*n*alpha*pc_alpha_n/(pc_alpha*one_plus_pc_alpha_n**(m+1))
        if (pc >= pth) then
          dSe_pc_at_pth = -m*n*(1.d0 + (alpha*pth)**n)**(-1.d0-m)*(alpha**n*pth**(n-1.d0))
          Se = (pc - 1.d8)*dSe_pc_at_pth
          dSe_pc = dSe_pc_at_pth
        ! write (*,*) option%myrank, 'pc:', pc, 'Se:', Se, 'dSe_pc', dSe_pc 
        endif 
        function_B = 1.d0/Se
        dfunc_B_pl = 1.d0/(Se**(2.d0))*dSe_pc        
      endif
      if (temperature >= 0.d0) then
        function_A = 1.d0
        dfunc_A_temp = 0.d0
      else
        gamma = den_ice*heat_of_fusion*interfacial_tensions_ratio
        pc_il = gamma*(-(temperature))/T_0
        alpha = saturation_function%alpha
        m = saturation_function%m
        n = 1.d0/(1.d0 - m)
        pc_il_alpha = pc_il*alpha
        pc_il_alpha_n = pc_il_alpha**n
        one_plus_pc_il_alpha_n = 1.d0 + pc_il_alpha_n
        Se_temp = one_plus_pc_il_alpha_n**(-m)
        function_A = 1.d0/Se_temp
        dfunc_A_temp = (gamma/T_0)*1.d0/(Se_temp**(2.d0))*(-m)* &
                       ((one_plus_pc_il_alpha_n)**(-m - 1.d0))*n* &
                       (pc_il**(n - 1.d0))*(alpha**n)
      endif           
    case default
      option%io_buffer = 'Ice module only supports Van Genuchten'
      call printErrMsg(option)
  end select
  
  liquid_saturation = 1.d0/(function_A + function_B - 1.d0)
  gas_saturation = liquid_saturation*(function_B - 1.d0)
  ice_saturation = liquid_saturation*(function_A - 1.d0)

  dsl_pl = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_B_pl)
  dsl_temp = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_A_temp)
  
  dsg_pl = dsl_pl*(function_B - 1.d0) + liquid_saturation*dfunc_B_pl
  dsg_temp = dsl_temp*(function_B - 1.d0)
  
  dsi_pl = dsl_pl*(function_A - 1.d0)
  dsi_temp = dsl_temp*(function_A - 1.d0) + liquid_saturation*dfunc_A_temp
  
  if (liquid_saturation > 1.d0) then
    print *, option%myrank, 'vG Liquid Saturation > 1:', liquid_saturation
  else if (liquid_saturation < 0.d0) then
    print *, option%myrank, 'vG Liquid Saturation < 0:', liquid_saturation
  endif

  if (gas_saturation > 1.d0) then
    print *, option%myrank, 'vG Gas Saturation > 1:', gas_saturation
  else if (liquid_saturation < 0.d0) then
    print *, option%myrank, 'vG Gas Saturation < 0:', gas_saturation
  endif
 
  if (ice_saturation > 1.d0) then
    print *, option%myrank, 'vG Ice Saturation > 1:', ice_saturation
  else if (liquid_saturation < 0.d0) then
    print *, option%myrank, 'vG Ice Saturation < 0:', ice_saturation
  endif

  select case(saturation_function%permeability_function_itype)
    case(MUALEM)
      one_over_m = 1.d0/m
      liq_sat_one_over_m = liquid_saturation**one_over_m
      liquid_relative_perm = sqrt(liquid_saturation)* &
                             (1.d0 - (1.d0 - liq_sat_one_over_m)**m)**2.d0
      if (liquid_saturation == 1.d0) then
        dkr_ds_liq = 0.d0
      else
        dkr_ds_liq = 0.5d0*liquid_relative_perm/liquid_saturation + &
                     2.d0*liquid_saturation**(one_over_m - 0.5d0)* &
                     (1.d0 - liq_sat_one_over_m)**(m - 1.d0)* &
                     (1.d0 - (1.d0 - liq_sat_one_over_m)**m)
      endif
        dkr_pl = dkr_ds_liq*dsl_pl
        dkr_temp = dkr_ds_liq*dsl_temp
    case default
      option%io_buffer = 'Ice module only supports Mualem' 
      call printErrMsg(option)
  end select
  
!  write(*,*) 'rank:', option%myrank, 'sl:', liquid_saturation, &
!  'sg:', gas_saturation, 'si:', ice_saturation, 'dsl_pl:', dsl_pl, &
! 'dsl_temp:', dsl_temp, 'dsg_pl:', dsg_pl, 'dsg_temp:', dsg_temp, &
!  'dsi_pl:', dsi_pl, 'dsi_temp:', dsi_temp, 'kr:', liquid_relative_perm, &
!  'dkr_pl:', dkr_pl, 'dkr_temp:', dkr_temp

end subroutine SaturationFunctionComputeIce

! ************************************************************************** !
!
! SatFuncGetRelPermFromSat: Calculates relative permeability from
!                           phase saturation
! author: Glenn Hammond
! date: 03/05/11
!
! ************************************************************************** !
subroutine SatFuncGetRelPermFromSat(saturation,relative_perm,dkr_Se, &
                                    saturation_function,iphase, &
                                    derivative,option)

  use Option_module
  
  implicit none

  PetscReal :: saturation, relative_perm, dkr_Se
  PetscInt :: iphase
  type(saturation_function_type) :: saturation_function
  PetscBool :: derivative
  type(option_type) :: option

  PetscReal :: m, Sr
  PetscReal :: Se, one_over_m, Se_one_over_m
  
  Sr = saturation_function%Sr(iphase)
  Se = (saturation-Sr)/(1.d0-Sr)
  
  ! if saturation is below residual saturation (this is possible with 
  ! two-phase with dry-out), need to bail out.
  if (Se <= 0.d0) then
    relative_perm = 0.d0
    dkr_Se = 0.d0
    return
  else if (Se > 1.d0) then
    Se = 1.d0
  endif
    
  ! compute relative permeability
  select case(saturation_function%permeability_function_itype)
    case(BURDINE)
      m = saturation_function%m
      one_over_m = 1.d0/m
      Se_one_over_m = Se**one_over_m
      relative_perm = Se*Se*(1.d0-(1.d0-Se_one_over_m)**m)
      if (derivative) then
        dkr_Se = 2.d0*relative_perm/Se + &
                 Se*Se_one_over_m*(1.d0-Se_one_over_m)**(m-1.d0)
      endif
    case(MUALEM)
      m = saturation_function%m
      one_over_m = 1.d0/m
      Se_one_over_m = Se**one_over_m
      relative_perm = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
      if (derivative) then
        dkr_Se = 0.5d0*relative_perm/Se+ &
                 2.d0*Se**(one_over_m-0.5d0)* &
                      (1.d0-Se_one_over_m)**(m-1.d0)* &
                      (1.d0-(1.d0-Se_one_over_m)**m)
      endif
    case default
      option%io_buffer = 'Unknown relative permeabilty function' 
      call printErrMsg(option)
  end select

end subroutine SatFuncGetRelPermFromSat

! ************************************************************************** !
!
! CapillaryPressureThreshold: Computes the capillary pressure threshold 
! after which instead of van Genuchten a linear function is used upto 100 Mpa
! capillary pressure. The saturation at 100 Mpa is set to zero
! This threshold value depends only on van Genuchten parameters alpha and lambda
! This is used mainly for ice problem, so that the pressure doesnt go to large
! negative values
! author: Satish Karra
! date: 09/12/12
!
! ************************************************************************** !

subroutine CapillaryPressureThreshold(saturation_function,cap_threshold,option)

  use Option_module
  
  implicit none
  
  PetscReal :: alpha, lambda, cap_threshold
  type(option_type) :: option
  type(saturation_function_type) :: saturation_function

  
  PetscReal :: gamma, p_new, res_value, jac_value, p_old
  PetscReal, parameter :: eps = 1.d-8
  PetscInt, parameter :: maxit = 100
  PetscInt :: i 
  
  alpha = saturation_function%alpha
  lambda = saturation_function%m
  alpha = alpha*1.d6
  gamma = 1.d0/(1.d0 - lambda)
  
  p_old = 99.d0
  
  
  do i = 1, maxit
    call ResidualCapPressThre(p_old,alpha,lambda,gamma,res_value)
    call JacobianCapPressThre(p_old,alpha,lambda,gamma,jac_value)
    p_new = p_old - res_value/jac_value
    !write (*,*) 'rank:', option%myrank, 'iter:', i, 'p_new:', p_new, 'p_old:', p_old, &
    !  'residual:', res_value, 'jacobian:', jac_value
    if ((abs(p_new - p_old) < eps)) exit
    p_old = p_new
  enddo
  
  cap_threshold = p_old*1.d6 ! convert to Pa
  
end subroutine CapillaryPressureThreshold

! ************************************************************************** !
!
! ResidualCapPressThre: Computes the residual to calculate capillary pressure
! thresold in the subroutine CapillaryPressureThreshold
! author: Satish Karra
! date: 09/12/12
!
! ************************************************************************** !

subroutine ResidualCapPressThre(p,alpha,lambda,gamma,res)

  implicit none
  
  PetscReal :: p, alpha, lambda, gamma, res
  
  res = lambda*gamma*alpha**gamma*p**(gamma-1.0)*1.d2 - &
        (alpha*p)**gamma*(1.d0 + gamma*lambda) - 1.d0

end subroutine ResidualCapPressThre


! ************************************************************************** !
!
! JacobianCapPressThre: Computes the jacobian to calculate capillary pressure
! thresold in the subroutine CapillaryPressureThreshold
! author: Satish Karra
! date: 09/12/12
!
! ************************************************************************** !

subroutine JacobianCapPressThre(p,alpha,lambda,gamma,jac)

  implicit none
  
  PetscReal :: p, alpha, lambda, gamma, jac
  
  jac = lambda*gamma*alpha**gamma*(gamma - 1.d0)*p**(gamma - 2.d0)*1.d2 - &
        alpha**gamma*gamma*p**(gamma - 1.d0)*(1.d0 + gamma*lambda)
  
 
end subroutine JacobianCapPressThre


! ************************************************************************** !
!
! SatFuncGetCapillaryPressure: Computes the capillary pressure as a function of 
!                           pressure
! author: Glenn Hammond
! date: 06/03/09
!
! ************************************************************************** !
subroutine SatFuncGetCapillaryPressure(capillary_pressure,saturation, &
                                       saturation_function,option)

  use Option_module
  
  implicit none

  PetscReal :: capillary_pressure, saturation
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option

  PetscInt :: iphase
  PetscReal :: alpha, lambda, m, n, Sr, one_over_alpha
  PetscReal :: pc, Se
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda
  
  iphase = 1

  Sr = saturation_function%Sr(iphase)
  if (saturation <= Sr) then
    capillary_pressure = saturation_function%pcwmax
    return
  endif
    
  ! compute saturation
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      if (saturation >= 1.d0) then
        capillary_pressure = 0.d0
        return
      else
        alpha = saturation_function%alpha
        m = saturation_function%m
        n = 1.d0/(1.d0-m)
        Se = (saturation-Sr)/(1.d0-Sr)
        one_plus_pc_alpha_n = Se**(-1.d0/m)
        pc_alpha_n = one_plus_pc_alpha_n - 1.d0
        pc_alpha = pc_alpha_n**(1.d0/n)
        capillary_pressure = pc_alpha/alpha
      endif
      ! compute relative permeability
    case(BROOKS_COREY)
      alpha = saturation_function%alpha
      one_over_alpha = 1.d0/alpha
!      pc = option%reference_pressure-pressure
      if (saturation >= 1.d0) then
        capillary_pressure = one_over_alpha
        return
      else
        lambda = saturation_function%lambda
        Sr = saturation_function%Sr(iphase)
        Se = (saturation-Sr)/(1.d0-Sr)
        pc_alpha_neg_lambda = Se
        capillary_pressure = (pc_alpha_neg_lambda**(-1.d0/lambda))/alpha
      endif
#if 0
    case(THOMEER_COREY)
      pc = option%reference_pressure-pressure
      por = auxvar1
      perm = auxvar2*1.013202d15 ! convert from m^2 to mD
      Fg = saturation_function%alpha
      a = saturation_function%m
      Pd = 100.d0*por/sqrt(perm/(3.8068d0*(Fg**(-1.334d0)))) ! psi
      PHg = 9.63051d-4*pc
      if (PHg > Pd) then
        saturation = 1.d0-exp(-Fg/log10(PHg/Pd))
#if 0
        alpha = pc*(1.d0+1.d-8)
        m = 9.63051d-4*alpha
        n = 1.d0-exp(-Fg/log10(m/Pd))
        n = (n-saturation)/(alpha-pc)
#endif        
        dsat_pc = (saturation-1.d0)*Fg/(log10(PHg/Pd)**2.d0)/(pc*2.30258509d0)
        ! Sr assumed to be zero
        relative_perm = saturation**a
        dkr_pc = a*saturation**(a-1.d0)*dsat_pc
      else
        saturation = 1.d0
        relative_perm = 1.d0
        return
      endif
#endif
    case default
      option%io_buffer = 'Unknown saturation function'
      call printErrMsg(option)
  end select

end subroutine SatFuncGetCapillaryPressure

! ************************************************************************** !
!
! SaturationFunctionGetID: Returns the ID of the saturation function named
!                          "saturation_function_name"
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
function SaturationFunctionGetID(saturation_function_list, &
                                 saturation_function_name, &
                                 material_property_name, option)

  use Option_module
  use String_module
  
  type(saturation_function_type), pointer :: saturation_function_list
  character(len=MAXWORDLENGTH) :: saturation_function_name
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: SaturationFunctionGetID
  PetscBool :: found
  type(saturation_function_type), pointer :: cur_saturation_function

  found = PETSC_FALSE
  cur_saturation_function => saturation_function_list
  do 
    if (.not.associated(cur_saturation_function)) exit
    if (StringCompare(saturation_function_name, &
                      cur_saturation_function%name,MAXWORDLENGTH)) then
      found = PETSC_TRUE
      SaturationFunctionGetID = cur_saturation_function%id
      return
    endif
    cur_saturation_function => cur_saturation_function%next
  enddo
  if (.not.found) then
    option%io_buffer = 'Saturation function "' // &
             trim(saturation_function_name) // &
             '" in material property "' // &
             trim(material_property_name) // &
             '" not found among available saturation functions.'
    call printErrMsg(option)    
  endif

end function SaturationFunctionGetID

! ************************************************************************** !
!
! SaturationFunctionDestroy: Destroys a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine SaturationFunctionDestroy(saturation_function)

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  
  if (.not.associated(saturation_function)) return
  
  call SaturationFunctionDestroy(saturation_function%next)
    
  if (associated(saturation_function%Sr)) deallocate(saturation_function%Sr)
  nullify(saturation_function%Sr)
    
  deallocate(saturation_function)
  nullify(saturation_function)
  
end subroutine SaturationFunctionDestroy

end module Saturation_Function_module
