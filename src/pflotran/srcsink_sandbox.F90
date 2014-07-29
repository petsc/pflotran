module SrcSink_Sandbox_module

  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_WIPP_Gas_class
  use SrcSink_Sandbox_Mass_Rate_class
  use SrcSink_Sandbox_Downreg_class
  use SrcSink_Sandbox_WIPP_Well_class
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  class(srcsink_sandbox_base_type), pointer, public :: sandbox_list

  interface SSSandboxRead
    module procedure SSSandboxRead1
    module procedure SSSandboxRead2
  end interface
  
  interface SSSandboxDestroy
    module procedure SSSandboxDestroy1
    module procedure SSSandboxDestroy2
  end interface
  
  public :: SSSandboxInit, &
            SSSandboxRead, &
            SSSandboxSetup, &
            SSSandbox, &
            SSSandboxDestroy

contains

! ************************************************************************** !

subroutine SSSandboxInit(option)
  ! 
  ! Initializes the sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 
  use Option_module
  implicit none
  type(option_type) :: option

  if (associated(sandbox_list)) then
    call SSSandboxDestroy()
  endif
  nullify(sandbox_list)

end subroutine SSSandboxInit

! ************************************************************************** !

subroutine SSSandboxSetup(region_list,option)
  ! 
  ! Calls all the initialization routines for all source/sinks in
  ! the sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Region_module
  
  implicit none
  
  type(option_type) :: option
  type(region_list_type) :: region_list
  
  class(srcsink_sandbox_base_type), pointer :: cur_sandbox  

  ! sandbox source/sinks
  cur_sandbox => sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%Setup(region_list,option)
    cur_sandbox => cur_sandbox%next
  enddo 

end subroutine SSSandboxSetup

! ************************************************************************** !

subroutine SSSandboxRead1(input,option)
  ! 
  ! Reads input deck for source/sink sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(input_type) :: input
  type(option_type) :: option

  call SSSandboxRead(sandbox_list,input,option)

end subroutine SSSandboxRead1

! ************************************************************************** !

subroutine SSSandboxRead2(local_sandbox_list,input,option)
  ! 
  ! Reads input deck for src/sink sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  class(srcsink_sandbox_base_type), pointer :: local_sandbox_list  
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  class(srcsink_sandbox_base_type), pointer :: new_sandbox, cur_sandbox
  
  nullify(new_sandbox)
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SOURCE_SINK_SANDBOX')
    call StringToUpper(word)   

    select case(trim(word))
      case('WIPP-WELL')
        new_sandbox => WIPPWellCreate()
      case('WIPP-GAS_GENERATION')
        new_sandbox => WIPPGasGenerationCreate()
      case('MASS_RATE')
        new_sandbox => MassRateCreate()
      case('MASS_RATE_DOWNREGULATED')
        new_sandbox => DownregCreate()
      case default
        option%io_buffer = 'SRCSINK_SANDBOX keyword: ' // &
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
  
end subroutine SSSandboxRead2

! ************************************************************************** !

subroutine SSSandbox(residual,Jacobian,compute_derivative, &
                     grid,material_auxvars,option)
  ! 
  ! Evaluates source/sink term storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Grid_module
  use Material_Aux_class, only: material_auxvar_type
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  PetscBool :: compute_derivative
  Vec :: residual
  Mat :: Jacobian
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscReal, pointer :: r_p(:)
  PetscReal :: res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  PetscInt :: i, local_id, ghosted_id, istart, iend
  PetscReal :: aux_real(0)
  PetscErrorCode :: ierr
  
  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif
  
  cur_srcsink => sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
      do i = 1, size(cur_srcsink%region%cell_ids)
        local_id = cur_srcsink%region%cell_ids(i)
        ghosted_id = grid%nL2G(local_id)
        res = 0.d0
        Jac = 0.d0
        call cur_srcsink%Evaluate(res,Jac,compute_derivative, &
                                  material_auxvars(ghosted_id), &
                                  aux_real,option)
        if (compute_derivative) then
          call MatSetValuesBlockedLocal(Jacobian,1,ghosted_id-1,1, &
                                        ghosted_id-1,Jac,ADD_VALUES, &
                                        ierr);CHKERRQ(ierr)
        else
          iend = local_id*option%nflowdof
          istart = iend - option%nflowdof + 1
          r_p(istart:iend) = r_p(istart:iend) + res
        endif
      enddo
    cur_srcsink => cur_srcsink%next
  enddo
  
  if (.not.compute_derivative) then
    call VecRestoreArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine SSSandbox

! ************************************************************************** !

subroutine SSSandboxDestroy1()
  ! 
  ! Destroys master sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none

  call SSSandboxDestroy(sandbox_list)
  
end subroutine SSSandboxDestroy1

! ************************************************************************** !

subroutine SSSandboxDestroy2(local_sandbox_list)
  ! 
  ! Destroys arbitrary sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none

  class(srcsink_sandbox_base_type), pointer :: local_sandbox_list

  class(srcsink_sandbox_base_type), pointer :: cur_sandbox, prev_sandbox
  
  ! sandbox source/sinks
  cur_sandbox => local_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    prev_sandbox => cur_sandbox%next
    call cur_sandbox%Destroy()
    deallocate(cur_sandbox)
    cur_sandbox => prev_sandbox
  enddo  

end subroutine SSSandboxDestroy2

end module SrcSink_Sandbox_module
