module Reaction_Sandbox_Template_class

! 1. Change all references to "Template" as desired to rename the module and
!    and subroutines within the module. 

#include "finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
! 2. Add module variables here.  Note that one must use the PETSc data types 
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0
  
  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_template_type
! 3. Add variables/arrays associated with new reaction.  E.g.
!   PetscInt :: example_integer
!   PetscInt, pointer :: example_integer_array(:)
!    Character strings must be sized to either MAXWORDLENGTH or 
!    MAXSTRINGLENGTH where max string length is a very long string.
  contains
    procedure, public :: ReadInput => TemplateRead
    procedure, public :: Setup => TemplateSetup
    procedure, public :: Evaluate => TemplateReact
    procedure, public :: Destroy => TemplateDestroy
  end type reaction_sandbox_template_type

  public :: TemplateCreate

contains

! ************************************************************************** !

function TemplateCreate()
  ! 
  ! Allocates template reaction object.
  ! 
  ! Author: John Doe (replace in all subroutine headers with name of developer)
  ! Date: 00/00/00 (replace in all subroutine headers with current date)
  ! 

  implicit none
  
  class(reaction_sandbox_template_type), pointer :: TemplateCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(TemplateCreate)
! TemplateCreate%example_integer = 0
! nullify(TemplateCreate%example_integer_array)
  nullify(TemplateCreate%next)
  
end function TemplateCreate

! ************************************************************************** !

subroutine TemplateRead(this,input,option)
  ! 
  ! Reads input deck for template reaction parameters (if any)
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  
  implicit none
  
  class(reaction_sandbox_template_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,TEMPLATE')
    call StringToUpper(word)   

    select case(trim(word))

      ! Example Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     TEMPLATE
      !       EXAMPLE_INTEGER 1
      !       EXAMPLE_INTEGER_ARRAY 2 3 4
      !     END
      !   : end user defined input
      !   END
      !   ...
      ! END

! 5. Add case statement for reading variables.  E.g.
!     case('EXAMPLE_INTEGER')
! 6. Read the variable
!       call InputReadInt(input,option,this%example_integer)  
! 7. Inform the user of any errors if not read correctly.
!       call InputErrorMsg(input,option,'example_integer', & 
!                          'CHEMISTRY,REACTION_SANDBOX,TEMPLATE') 
! 8. Repeat for other variables
!     case('EXAMPLE_INTEGER_Array')
!       allocate(this%example_integer_array(3))
!       this%example_integer_array = 0
!       do i = 1, 3
!         call InputReadInt(input,option,this%example_integer_array(i))  
!         call InputErrorMsg(input,option,'example_integer_array', & 
!                            'CHEMISTRY,REACTION_SANDBOX,TEMPLATE') 
!       
!       enddo  
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,TEMPLATE keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine TemplateRead

! ************************************************************************** !

subroutine TemplateSetup(this,reaction,option)
  ! 
  ! Sets up the template reaction either with parameters either
  ! read from the input deck or hardwired.
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Reaction_Aux_module, only : reaction_type
  use Option_module

  implicit none
  
  class(reaction_sandbox_template_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize 
      
end subroutine TemplateSetup

! ************************************************************************** !

subroutine TemplateReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,porosity,volume,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Option_module
  use Reaction_Aux_module
  
  implicit none
  
  class(reaction_sandbox_template_type) :: this  
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
  !                      'due to assumptions in Template'
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

! 10. Add code for Residual evaluation
  
  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

    ! remove this error statement if analytical derivatives are provided.
    option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
                       'due to assumptions in Template'
    call printErrMsg(option)
    ! add code for Jacobian evaluation
    !geh: If you study other reactions built into PFLOTRAN you will see that I
    !     first calculate the residual and then modify it to get the Jacobian
    !     entries.  This is my preference that you do not need to follow.
  endif
  
end subroutine TemplateReact

! ************************************************************************** !

subroutine TemplateDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  implicit none
  
  class(reaction_sandbox_template_type) :: this  

! 12. Add code to deallocate contents of the template object
! deallocate(this%example_integer_array)
! nullify(this%example_integer_array) 

end subroutine TemplateDestroy

end module Reaction_Sandbox_Template_class
