module Reaction_Sandbox_Example_class

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
!    float/real*8, and logical respectively.

! None

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_example_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
  contains
    procedure, public :: ReadInput => ExampleRead
    procedure, public :: Setup => ExampleSetup
    procedure, public :: Evaluate => ExampleReact
    procedure, public :: Destroy => ExampleDestroy
  end type reaction_sandbox_example_type

  public :: ExampleCreate

contains

! ************************************************************************** !

function ExampleCreate()
  ! 
  ! Allocates example reaction object.
  ! 
  ! Author: John Doe (replace in all subroutine headers with name of developer)
  ! Date: 00/00/00 (replace in all subroutine headers with current date)
  ! 

  implicit none
  
  class(reaction_sandbox_example_type), pointer :: ExampleCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers.
  allocate(ExampleCreate)
  ExampleCreate%species_name = ''
  ExampleCreate%species_id = 0
  ExampleCreate%rate_constant = 0.d0
  nullify(ExampleCreate%next)
  
end function ExampleCreate

! ************************************************************************** !

subroutine ExampleRead(this,input,option)
  ! 
  ! Reads input deck for example reaction parameters (if any)
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  
  implicit none
  
  class(reaction_sandbox_example_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,EXAMPLE')
    call StringToUpper(word)   

    select case(trim(word))

! Nothing to add given the reaction is hardwired.  But note that the empty
! block for the sandbox reaction EXAMPLE must still exist in the input file, 
! i.e.

! CHEMISTRY
!   ...
!   REACTION_SANDBOX
!     EXAMPLE
!     /
!   /
!   ...
! END
 
! 5. Add case statement for reading variables.  E.g.
!     case('EXAMPLE_INTEGER')
! 6. Read the variable
!       call InputReadInt(input,option,this%example_integer)  
! 7. Inform the user of any errors if not read correctly.
!       call InputErrorMsg(input,option,'example_integer', & 
!                          'CHEMISTRY,REACTION_SANDBOX,EXAMPLE') 
! 8. Repeat for other variables
!     case('EXAMPLE_INTEGER_Array')
!       allocate(this%example_integer_array(3))
!       this%example_integer_array = 0
!       do i = 1, 3
!         call InputReadInt(input,option,this%example_integer_array(i))  
!         call InputErrorMsg(input,option,'example_integer_array', & 
!                            'CHEMISTRY,REACTION_SANDBOX,EXAMPLE') 
!       
!       enddo

      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,EXAMPLE keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine ExampleRead

! ************************************************************************** !

subroutine ExampleSetup(this,reaction,option)
  ! 
  ! Sets up the example reaction either with parameters either
  ! read from the input deck or hardwired.
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_example_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize   
  this%species_name = 'A(aq)'
  this%species_id = &
    GetPrimarySpeciesIDFromName(this%species_name,reaction,option)
  this%rate_constant = 4.3959d-9 ! 5 yr half-life
  
end subroutine ExampleSetup

! ************************************************************************** !

subroutine ExampleReact(this,Residual,Jacobian,compute_derivative, &
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
  
  class(reaction_sandbox_example_type) :: this  
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
  PetscReal :: L_water
  
! 10. Add code for Residual evaluation
  L_water = porosity*global_auxvar%sat(iphase)*volume*1.d3
  ! alway subtract contribution from residual
  ! units = (mol/sec)
  Residual(this%species_id) = Residual(this%species_id) - &
    this%rate_constant * &  ! 1/sec
    L_water * & ! L water
    rt_auxvar%total(this%species_id,iphase) ! mol/L water
    
  
  
  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

    ! always add contribution to Jacobian
    ! units = (mol/sec)*(kg water/mol) = kg water/sec
    Jacobian(this%species_id,this%species_id) = &
    Jacobian(this%species_id,this%species_id) + &
      this%rate_constant * & ! 1/sec
      L_water * & ! L water
                  ! kg water/L water
      rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) 

  endif
  
end subroutine ExampleReact

! ************************************************************************** !

subroutine ExampleDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  implicit none
  
  class(reaction_sandbox_example_type) :: this  

! 12. Add code to deallocate contents of the example object

end subroutine ExampleDestroy

end module Reaction_Sandbox_Example_class
