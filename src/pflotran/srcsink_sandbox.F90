module SrcSink_Sandbox_module

  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_WIPP_Gas_class
  use SrcSink_Sandbox_Mass_Rate_class
  use SrcSink_Sandbox_Downreg_class
  use SrcSink_Sandbox_WIPP_Well_class
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  class(srcsink_sandbox_base_type), pointer, public :: ss_sandbox_list

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
            SSSandboxUpdate, &
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

  if (associated(ss_sandbox_list)) then
    call SSSandboxDestroy()
  endif
  nullify(ss_sandbox_list)

end subroutine SSSandboxInit

! ************************************************************************** !

subroutine SSSandboxSetup(region_list,grid,option)
  ! 
  ! Calls all the initialization routines for all source/sinks in
  ! the sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Region_module
  use Grid_module
  
  implicit none
  
  type(region_list_type) :: region_list
  type(grid_type) :: grid
  type(option_type) :: option
  
  class(srcsink_sandbox_base_type), pointer :: cur_sandbox  
  class(srcsink_sandbox_base_type), pointer :: prev_sandbox  
  class(srcsink_sandbox_base_type), pointer :: next_sandbox  

  ! sandbox source/sinks
  cur_sandbox => ss_sandbox_list
  nullify(prev_sandbox)
  do
    if (.not.associated(cur_sandbox)) exit
    next_sandbox => cur_sandbox%next
    call cur_sandbox%Setup(region_list,grid,option)
    ! destory if not on process
    if (.not.associated(cur_sandbox%region%cell_ids)) then
      if (associated(prev_sandbox)) then
        prev_sandbox%next => next_sandbox
      else
        ss_sandbox_list => next_sandbox
      endif
      nullify(cur_sandbox%next)
      call SSSandboxDestroy(cur_sandbox)
    endif
    if (associated(cur_sandbox)) prev_sandbox => cur_sandbox
    cur_sandbox => next_sandbox
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
  
  type(input_type), pointer :: input
  type(option_type) :: option

  call SSSandboxRead(ss_sandbox_list,input,option)

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
  type(input_type), pointer :: input
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
        call InputKeywordUnrecognized(word,'SRCSINK_SANDBOX',option)
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
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"

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
  
  cur_srcsink => ss_sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
      do i = 1, cur_srcsink%region%num_cells
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

subroutine SSSandboxUpdate(sandbox_list,time,option,output_option)
  ! 
  ! Updates datasets associated with a sandbox, if they exist
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/22/15
  ! 
  use Option_module
  use Output_Aux_module

  implicit none

  class(srcsink_sandbox_base_type), pointer :: sandbox_list
  PetscReal :: time
  type(option_type) :: option
  type(output_option_type) :: output_option

  class(srcsink_sandbox_base_type), pointer :: cur_sandbox
  
  cur_sandbox => sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%Update(time,option)
    cur_sandbox => cur_sandbox%next
  enddo  

end subroutine SSSandboxUpdate

! ************************************************************************** !

function SSSandboxOutputFilename(option)
  ! 
  ! Generates filename for source/sink sandbox output
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/16

  use Option_module

  implicit none
  
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: SSSandboxOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  SSSandboxOutputFilename = trim(option%global_prefix) // &
                            trim(option%group_prefix) // &
                            '-ss_mass-' // trim(adjustl(word)) // '.dat'
  
end function SSSandboxOutputFilename  

! ************************************************************************** !

subroutine SSSandboxOutputHeader(sandbox_list,grid,option,output_option)
  ! 
  ! Writes header for source/sink sandbox output
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/16

  use Option_module
  use Output_Aux_module
  use Grid_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  class(srcsink_sandbox_base_type), pointer :: sandbox_list
  type(grid_type) :: grid
  type(option_type) :: option
  type(output_option_type) :: output_option

  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  character(len=MAXSTRINGLENGTH) :: cell_string
!  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: icolumn, i
  PetscInt :: local_id
  
  filename = SSSandboxOutputFilename(option)
  open(unit=IUNIT_TEMP,file=filename,action="write",status="replace")  
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif 
  
  write(IUNIT_TEMP,'(a)',advance="no") ' "Time [' // &
    trim(output_option%tunit) // ']"'

  cur_srcsink => sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
      do i = 1, cur_srcsink%region%num_cells
        local_id = cur_srcsink%region%cell_ids(i)

        ! cell natural id
        write(cell_string,*) grid%nG2A(grid%nL2G(local_id))
        cell_string = ' (' // trim(cur_srcsink%region%name) // ' ' // &
                      trim(adjustl(cell_string)) // ')'
        ! coordinate of cell
!        x_string = BestFloat(cur_waste_form%coordinate%x,1.d4,1.d-2)
!        y_string = BestFloat(cur_waste_form%coordinate%y,1.d4,1.d-2)
!        z_string = BestFloat(cur_waste_form%coordinate%z,1.d4,1.d-2)
!        cell_string = trim(cell_string) // &
!                 ' (' // trim(adjustl(x_string)) // &
!                 ' ' // trim(adjustl(y_string)) // &
!                 ' ' // trim(adjustl(z_string)) // ')'
        variable_string = ' Mass Flux'
        ! cumulative
        units_string = 'mol'
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                  cell_string,icolumn)
        ! instantaneous
        units_string = 'mol/' // trim(adjustl(output_option%tunit))
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                  cell_string,icolumn)
                                  
      enddo
    cur_srcsink => cur_srcsink%next
  enddo
  
  close(IUNIT_TEMP)
  
end subroutine SSSandboxOutputHeader

! ************************************************************************** !

subroutine SSSandboxOutput(sandbox_list,option,output_option)
  ! 
  ! Writes output for for source/sink sandbox
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/16

  use Option_module
  use Output_Aux_module

  implicit none
  
  class(srcsink_sandbox_base_type), pointer :: sandbox_list
  type(option_type) :: option
  type(output_option_type) :: output_option
  
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: i
  
  if (.not.associated(sandbox_list)) return
  
100 format(100es16.8)

  filename = SSSandboxOutputFilename(option)
  open(unit=IUNIT_TEMP,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the reactive transport step
  write(IUNIT_TEMP,100,advance="no") option%time / output_option%tconv
  
  cur_srcsink => sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    do i = 1, cur_srcsink%region%num_cells
      write(IUNIT_TEMP,100,advance="no") cur_srcsink%cumulative_mass(i), &
                                  cur_srcsink%instantaneous_mass_rate(i) * &
                                  output_option%tconv
    enddo
    cur_srcsink => cur_srcsink%next
  enddo
  close(IUNIT_TEMP)
  
end subroutine SSSandboxOutput

! ************************************************************************** !

subroutine SSSandboxDestroy1()
  ! 
  ! Destroys master sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none

  call SSSandboxDestroy(ss_sandbox_list)
  
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
