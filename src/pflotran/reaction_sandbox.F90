module Reaction_Sandbox_module

  use Reaction_Sandbox_Base_class
  use Reaction_Sandbox_CLM_CN_class
  
  implicit none
  
  private
  
#include "definitions.h"

  class(reaction_sandbox_base_type), pointer, public :: sandbox_list

  interface RSandboxRead
    module procedure :: RSandboxRead1
    module procedure :: RSandboxRead2
  end interface
  
  interface RSandboxDestroy
    module procedure :: RSandboxDestroy1
    module procedure :: RSandboxDestroy2
  end interface
  
  public :: RSandboxInit, &
            RSandboxRead, &
            RSandboxSkipInput, &
            RSandboxSetup, &
            RSandbox, &
            RSandboxDestroy

contains

! ************************************************************************** !
!
! RSandboxInit: Initializes the sandbox list
! author: Glenn Hammond
! date: 01/28/13
!
! ************************************************************************** !
subroutine RSandboxInit(option)
  use Option_module
  implicit none
  type(option_type) :: option
  if (associated(sandbox_list)) then
    call RSandboxDestroy()
  endif
  nullify(sandbox_list)

end subroutine RSandboxInit

! ************************************************************************** !
!
! RSandboxSetup: Calls all the initialization routines for all reactions in
!                the sandbox list
! author: Glenn Hammond
! date: 01/28/13
!
! ************************************************************************** !
subroutine RSandboxSetup(reaction,option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type 
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  class(reaction_sandbox_base_type), pointer :: cur_sandbox  

  ! sandbox reactions
  cur_sandbox => sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%init(reaction,option)
    cur_sandbox => cur_sandbox%next
  enddo 

end subroutine RSandboxSetup

! ************************************************************************** !
!
! RSandboxRead1: Reads input deck for reaction sandbox parameters
! author: Glenn Hammond
! date: 05/16/13
!
! ************************************************************************** !
subroutine RSandboxRead1(input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(input_type) :: input
  type(option_type) :: option

  call RSandboxRead(sandbox_list,input,option)

end subroutine RSandboxRead1

! ************************************************************************** !
!
! RSandboxRead: Reads input deck for reaction sandbox parameters
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxRead2(local_sandbox_list,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_base_type), pointer :: local_sandbox_list  
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  class(reaction_sandbox_base_type), pointer :: new_sandbox, cur_sandbox
  
  call RSandboxInit(option)
  
  nullify(new_sandbox)
  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REACTION_SANDBOX')
    call StringToUpper(word)   

    select case(trim(word))
      case('CLM-CN')
        new_sandbox => CLM_CN_Create()
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
    
    call new_sandbox%ReadInput(input,option)
    
    if (.not.associated(local_sandbox_list)) then
      local_sandbox_list => new_sandbox
    else
      cur_sandbox => local_sandbox_list
      do
        if (.not.associated(cur_sandbox%next)) exit
        cur_sandbox => cur_sandbox%next
      enddo
      cur_sandbox%next => new_sandbox
    endif
  enddo
  
end subroutine RSandboxRead2

! ************************************************************************** !
!
! RSandboxSkipInput: Intelligently skips over REACTION_SANDBOX block
! author: Glenn Hammond
! date: 02/04/13
!
! ************************************************************************** !
subroutine RSandboxSkipInput(input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  
  class(reaction_sandbox_base_type), pointer :: dummy_list
  
  nullify(dummy_list)
  call RSandboxRead(dummy_list,input,option)
  call RSandboxDestroy(dummy_list)
  
end subroutine RSandboxSkipInput

! ************************************************************************** !
!
! RSandbox: Evaluates reaction storing residual and/or Jacobian
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandbox(Residual,Jacobian,compute_derivative,rt_auxvar, &
                    global_auxvar,porosity,volume,reaction,option)

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  implicit none

  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  
  class(reaction_sandbox_base_type), pointer :: cur_reaction
  
  cur_reaction => sandbox_list
  do
    if (.not.associated(cur_reaction)) exit
!    select type(cur_reaction)
!      class is(reaction_sandbox_clm_cn_type)
        call cur_reaction%Evaluate(Residual,Jacobian,compute_derivative, &
                                   rt_auxvar,global_auxvar,porosity,volume, &
                                   reaction,option)
!    end select
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandbox

! ************************************************************************** !
!
! RSandboxDestroy1: Destroys master sandbox list
! author: Glenn Hammond
! date: 05/16/13
!
! ************************************************************************** !
subroutine RSandboxDestroy1()

  implicit none

  call RSandboxDestroy(sandbox_list)
  
end subroutine RSandboxDestroy1

! ************************************************************************** !
!
! RSandboxDestroy2: Destroys arbitrary sandbox list
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxDestroy2(local_sandbox_list)

  implicit none

  class(reaction_sandbox_base_type), pointer :: local_sandbox_list

  class(reaction_sandbox_base_type), pointer :: cur_sandbox, prev_sandbox
  
  ! sandbox reactions
  cur_sandbox => local_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    prev_sandbox => cur_sandbox%next
    call cur_sandbox%Destroy()
    cur_sandbox => prev_sandbox
  enddo  

end subroutine RSandboxDestroy2

end module Reaction_Sandbox_module
