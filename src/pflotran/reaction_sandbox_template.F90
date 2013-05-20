module Reaction_Sandbox_Template_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_template_type
  contains
    procedure, public :: Init => RSandboxInit
    procedure, public :: ReadInput => RSandboxRead
    procedure, public :: Evaluate => RSandbox
    procedure, public :: Destroy => RSandboxDestroy
  end type reaction_sandbox_template_type

contains

! ************************************************************************** !
!
! RSandboxInit: Initializes reaction sandbox at beginning of simulation
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxInit(reaction_sandbox)

  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: reaction_sandbox
      
end subroutine RSandboxInit

! ************************************************************************** !
!
! RSandboxRead: Reads input deck for reaction sandbox parameters
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxRead(reaction_sandbox,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: reaction_sandbox
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
!  PetscReal :: example_real
!  PetscInt :: example_int
!  PetscBool :: example_bool
  
  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REACTION_SANDBOX')
    call StringToUpper(word)   

    select case(trim(word))

      ! Example Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     RATE_CONSTANT 1.d-5
      !   : end user defined input
      !   END
      !   ...
      ! END

      ! How to implement read in code:
      ! 1. Add module variable 'rate_constant' to top of module file.  The 
      !    following as data types are recommended (though the user is free 
      !    to use any Fortran data type as long as accept the risk of compiler 
      !    dependencies).
      !   PetscReal - real*8
      !   PetscInt - integer*4
      !   PetscBool - logical
      !   character(len=128) - String of length 128, best for reaction strings
      !   character(len=32) - String of length 32, species names, etc.

      !   Example (see also above at top of module):
      !     PetscReal :: rate_constant

      ! 2. Add case statement to this select case statement

      !   Example:
      !     case('RATE_CONSTANT')
      !       call InputReadDouble(input,option,rate_constant)  
      !       call InputDefaultMsg(input,option, &
      !                             'CHEMISTRY,REACTION_SANDBOX,RATE_CONSTANT') 

      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine RSandboxRead

! ************************************************************************** !
!
! RMicrobial: Evaluates reaction storing residual and/or Jacobian
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandbox(reaction_sandbox,Residual,Jacobian,compute_derivative, &
                    rt_auxvar,global_auxvar,porosity,volume,reaction,option)

  use Option_module
  use Reaction_Aux_module
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: reaction_sandbox  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: saturation
  
  ! Description of subroutine arguments:

  ! Residual - 1D array storing residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entires in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analtical derivative should
  !   be calculated.  The user must provide either the analytical derivatives 
  !   or a numerical approximation unless always running with 
  !   NUMERICAL_JACOBIAN_RXN defined in input deck.  If the use of 
  !   NUMERICAL_JACOBIAN_RXN is assumed, the user should provide an error 
  !   message when compute_derivative is true.  E.g.
  !
  !   option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
  !                      'due to assumptions in RSandbox'
  !   call printErrMsg(option)
  !
  ! rt_auxvar - Object holding chemistry information (e.g. concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.  
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations 
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g. saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - fluid density [mol/m^3] 
  !     global_auxvar%den_kg(iphase) - fluid density [kg/m^3] 
  !     global_auxvar%sat(iphase) - saturation 
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]                     
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.
  
  saturation = global_auxvar%sat(iphase)

  ! add code for Residual evaluation
  
  if (compute_derivative) then
    ! remove this error statement if analytical derivatives are provided.
    option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
                       'due to assumptions in RSandbox'
    call printErrMsg(option)
    ! add code for Jacobian evaluation
    !geh: If you study other reactions built into PFLOTRAN you will see that I
    !     first calculate the residual and then modify it to get the Jacobian
    !     entries.  This is my preference that you do not need to follow.
  endif
  
end subroutine RSandbox

! ************************************************************************** !
!
! RSandboxDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxDestroy(reaction_sandbox)

  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: reaction_sandbox  

end subroutine RSandboxDestroy

end module Reaction_Sandbox_Template_class
