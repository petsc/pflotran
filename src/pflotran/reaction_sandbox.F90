module Reaction_Sandbox_module

  use Reaction_Sandbox_Base_class
  use Reaction_Sandbox_CLM_CN_class
  
  implicit none
  
  private
  
#include "definitions.h"

  public :: RSandboxRead, &
            RSandbox, &
            RSandboxDestroy

contains

! ************************************************************************** !
!
! RSandboxRead: Reads input deck for reaction sandbox parameters
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxRead(sandbox_list,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_base_type), pointer :: sandbox_list
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  class(reaction_sandbox_base_type), pointer :: new_sandbox, cur_sandbox
  
  class(reaction_sandbox_clm_cn_type), pointer :: clm_cn
  
  nullify(new_sandbox)
  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REACTION_SANDBOX')
    call StringToUpper(word)   

    select case(trim(word))
      case('CLM_CN')
        clm_cn => CLM_CN_Create()
        call clm_cn%Init()
        new_sandbox => clm_cn
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
    
    if (.not.associated(sandbox_list)) then
      sandbox_list => new_sandbox
    else
      cur_sandbox => sandbox_list
      do
        if (.not.associated(cur_sandbox%next)) exit
        cur_sandbox => cur_sandbox%next
      enddo
      cur_sandbox%next => new_sandbox
    endif
  enddo
  
end subroutine RSandboxRead

! ************************************************************************** !
!
! RSandbox: Evaluates reaction storing residual and/or Jacobian
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandbox(sandbox_list, &
                    Residual,Jacobian,compute_derivative,rt_auxvar, &
                    global_auxvar,porosity,volume,reaction,option)

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  implicit none

  class(reaction_sandbox_base_type), pointer :: sandbox_list
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
    select type(cur_reaction)
      class is(reaction_sandbox_clm_cn_type)
        call cur_reaction%Evaluate(Residual,Jacobian,compute_derivative, &
                                   rt_auxvar,global_auxvar,porosity,volume, &
                                   reaction,option)
    end select
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandbox

! ************************************************************************** !
!
! RSandboxDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxDestroy()

  implicit none

end subroutine RSandboxDestroy

end module Reaction_Sandbox_module
