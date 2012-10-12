module Regression_module
 
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
 
  type, public :: regression_type
    type(regression_variable_type), pointer :: variable_list
    PetscInt, pointer :: natural_cell_ids(:)
    PetscInt :: num_cells_per_process
    PetscInt, pointer :: cells_per_process_natural_ids(:)
    Vec :: natural_cell_id_vec
    Vec :: cells_per_process_vec
    VecScatter :: scatter_natural_cell_id_gtos
    VecScatter :: scatter_cells_per_process_gtos
    type(regression_type), pointer :: next
  end type regression_type

  type, public :: regression_variable_type
    character(len=MAXSTRINGLENGTH) :: name
    type(regression_variable_type), pointer :: next
  end type regression_variable_type
  
  public :: RegressionRead, &
            RegressionCreateMapping, &
            RegressionOutput, &
            RegressionDestroy
  
contains

! ************************************************************************** !
!
! RegressionCreate: Creates a regression object
! author: Glenn Hammond
! date: 10/11/12
!
! ************************************************************************** !
function RegressionCreate()
  
  implicit none

  type(regression_type), pointer :: RegressionCreate
  
  type(regression_type), pointer :: regression
  
  allocate(regression)
  nullify(regression%variable_list)
  nullify(regression%natural_cell_ids)
  regression%num_cells_per_process = 0
  nullify(regression%cells_per_process_natural_ids)
  regression%natural_cell_id_vec = 0
  regression%cells_per_process_vec = 0
  regression%scatter_natural_cell_id_gtos = 0
  regression%scatter_cells_per_process_gtos = 0
  nullify(regression%next)
  RegressionCreate => regression

end function RegressionCreate

! ************************************************************************** !
!
! RegressionVariableCreate: Creates a regression variable object
! author: Glenn Hammond
! date: 10/11/12
!
! ************************************************************************** !
function RegressionVariableCreate()
  
  implicit none

  type(regression_variable_type), pointer :: RegressionVariableCreate
  
  type(regression_variable_type), pointer :: regression_variable
  
  allocate(regression_variable)
  regression_variable%name = ''
  nullify(regression_variable%next)
  RegressionVariableCreate => regression_variable

end function RegressionVariableCreate

! ************************************************************************** !
!
! RegressionRead: Reads in contents of a regression card
! author: Glenn Hammond
! date: 10/11/12
! 
! ************************************************************************** !
subroutine RegressionRead(regression,input,option)

  use Option_module
  use Input_module
  use String_module
  use Utility_module

  implicit none
  
  type(regression_type), pointer :: regression
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  type(regression_variable_type), pointer :: cur_variable, new_variable
  PetscInt :: count, max_cells
  PetscInt, pointer :: int_array(:)
  PetscErrorCode :: ierr

  regression => RegressionCreate()
  
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','REGRESSION')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('VARIABLES') 
        count = 0
        do 
          call InputReadFlotranString(input,option)
          if (InputCheckExit(input,option)) exit  

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'variable','REGRESSION,VARIABLES')
          call StringToUpper(word)
          new_variable => RegressionVariableCreate()
          new_variable%name = word
          if (.not.associated(regression%variable_list)) then
            regression%variable_list => new_variable
          else
            cur_variable%next => new_variable
          endif
          cur_variable => new_variable
        enddo
      case('CELLS')
        max_cells = 100
        allocate(int_array(max_cells))
        count = 0
        do 
          call InputReadFlotranString(input,option)
          if (InputCheckExit(input,option)) exit  

          count = count + 1
          if (count > max_cells) then
            call reallocateIntArray(int_array,max_cells)
          endif
          call InputReadInt(input,option,int_array(count))
          call InputErrorMsg(input,option,'natural cell id','REGRESSION,CELLS')
        enddo
        allocate(regression%natural_cell_ids(count))
        regression%natural_cell_ids = int_array(1:count)
        call PetscSortInt(count,regression%natural_cell_ids,ierr)        
        deallocate(int_array)
      case('CELLS_PER_PROCESS')
        call InputReadInt(input,option,regression%num_cells_per_process)
        call InputErrorMsg(input,option,'num cells per process','REGRESSION')
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in regression'    
        call printErrMsg(option)
    end select
    
  enddo
  
end subroutine RegressionRead

! ************************************************************************** !
!
! RegressionCreateMapping: Creates mapping between a natural mpi vec and a 
!                          sequential vec on io_rank
! author: Glenn Hammond
! date: 10/12/12
!
! ************************************************************************** !
subroutine RegressionCreateMapping(regression,realization)

  use Option_module
  use Realization_module
  use Grid_module
  use Discretization_module
  
  implicit none
  
#include "finclude/petscis.h"
#include "finclude/petscis.h90"

  type(regression_type), pointer :: regression
  type(realization_type) :: realization
  
  IS :: is_global
  PetscInt, allocatable :: int_array(:), int_array2(:)
  PetscInt :: i, upper_bound, lower_bound, count, temp_int
  PetscInt :: local_id
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  if (.not.associated(regression)) return
  
  grid => realization%patch%grid
  option => realization%option
  
  ! natural cell ids
  if (size(regression%natural_cell_ids) > 0) then
    call VecCreate(option%mycomm,regression%natural_cell_id_vec,ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(regression%natural_cell_id_vec, &
                       size(regression%natural_cell_ids), &
                       PETSC_DECIDE,ierr)  
    else
      call VecSetSizes(regression%natural_cell_id_vec,0, &
                       PETSC_DECIDE,ierr)  
    endif
    call VecSetFromOptions(regression%natural_cell_id_vec,ierr)
  
    ! determine how many of the natural cell ids are local
    allocate(int_array(size(regression%natural_cell_ids)))
    int_array = regression%natural_cell_ids
    ! convert to zero based
    int_array = int_array - 1
    call DiscretAOApplicationToPetsc(realization%discretization,int_array)
    ! convert back to one based
    int_array = int_array + 1
    allocate(int_array2(size(int_array)))
    int_array2 = -999
    lower_bound = grid%global_offset
    upper_bound = lower_bound + grid%nlmax
    count = 0
    do i = 1, size(int_array)
      if (int_array(i) > lower_bound .and. int_array(i) <= upper_bound) then
        count = count + 1
        int_array2(count) = int_array(i)
      endif
    enddo
    deallocate(int_array)
  
    ! create IS for global petsc cell ids
    allocate(int_array(count))
    int_array = int_array2(1:count) - 1 ! zero-based
    deallocate(int_array2)
    call ISCreateGeneral(option%mycomm,count,int_array,PETSC_COPY_VALUES, &
                         is_global,ierr)
    deallocate(int_array)
  
    ! create scatter context
    call VecScatterCreate(realization%field%work,is_global, &
                          regression%natural_cell_id_vec,PETSC_NULL_OBJECT, &
                          regression%scatter_natural_cell_id_gtos,ierr)
    call ISDestroy(is_global,ierr)
  endif

  if (regression%num_cells_per_process > 0) then
    ! determine minimum number of cells per process
    i = grid%nlmax
    call MPI_Allreduce(i,count,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (count < regression%num_cells_per_process) then
      option%io_buffer = 'Number of cells per process for Regression ' // &
        'exceeds minimum number of cells per process.  Truncating.'
      call printMsg(option)
      regression%num_cells_per_process = count
    endif
  
    ! cells ids per processor
    call VecCreate(option%mycomm,regression%cells_per_process_vec,ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(regression%cells_per_process_vec, &
                       regression%num_cells_per_process*option%mycommsize, &
                       PETSC_DECIDE,ierr)  
    else
      call VecSetSizes(regression%cells_per_process_vec,0, &
                       PETSC_DECIDE,ierr)  
    endif
    call VecSetFromOptions(regression%cells_per_process_vec,ierr)
  
    ! calculate interval
    temp_int = grid%nlmax / regression%num_cells_per_process
  
    ! calculate local cells ids
    allocate(int_array(regression%num_cells_per_process))
    do i = 1, regression%num_cells_per_process
      int_array(i) = temp_int*(i-1) + 1
    enddo
  
    ! create IS for global petsc cell ids
    ! convert local ids to global petsc ids
    int_array = int_array + lower_bound ! the global offset for this process
    int_array = int_array - 1 ! zero-based
    call ISCreateGeneral(option%mycomm,regression%num_cells_per_process, &
                         int_array,PETSC_COPY_VALUES,is_global,ierr)
    deallocate(int_array)
  
    ! create scatter context
    call VecScatterCreate(realization%field%work,is_global, &
                          regression%cells_per_process_vec, &
                          PETSC_NULL_OBJECT, &
                          regression%scatter_cells_per_process_gtos,ierr)
    call ISDestroy(is_global,ierr)
  
    ! fill in natural ids of these cells on the io_rank
    if (option%myrank == option%io_rank) then
      allocate(regression%cells_per_process_natural_ids( &
               regression%num_cells_per_process*option%mycommsize))
    endif
    call VecGetArrayF90(realization%field%work,vec_ptr,ierr)
    do local_id = 1, grid%nlmax
      vec_ptr(local_id) = grid%nG2A(grid%nL2G(local_id))
    enddo
    call VecRestoreArrayF90(realization%field%work,vec_ptr,ierr)
    call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                          realization%field%work, &
                          regression%cells_per_process_vec, &
                          INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                       realization%field%work, &
                       regression%cells_per_process_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
    if (option%myrank == option%io_rank) then
      call VecGetArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
      regression%cells_per_process_natural_ids(:) = int(vec_ptr(:)+0.1)
      call VecRestoreArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
    endif
  endif
  
end subroutine RegressionCreateMapping

! ************************************************************************** !
!
! RegressionOutput: Prints regression output through the io_rank
! author: Glenn Hammond
! date: 10/12/12
!
! ************************************************************************** !
subroutine RegressionOutput(regression,realization,flow_stepper, &
                            tran_stepper)

  use Realization_module
  use Timestepper_module
  use Option_module
  use Discretization_module
  use Output_module
  
  implicit none
  
  type(regression_type), pointer :: regression
  type(realization_type) :: realization
  ! these must be pointers as they can be null
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper  
  
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  PetscInt :: ivar, isubvar
  type(option_type), pointer :: option
  type(regression_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr
  
  if (.not.associated(regression)) return
  
  option => realization%option
  
  if (option%myrank == option%io_rank) then
    string = trim(option%global_prefix) // &
             trim(option%group_prefix) // &  
             '.regression'
    option%io_buffer = '--> write regression output file: ' // trim(string)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=string,action="write")
  endif
  
  call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)  
  
  cur_variable => regression%variable_list
  do 
    if (.not.associated(cur_variable)) exit
    ! need to move this elsewhere
    ivar = 0
    isubvar = 0
    select case(cur_variable%name)
      case('PRESSURE')
        ivar = LIQUID_PRESSURE
      case default
        option%io_buffer = 'Variable "' // & trim(cur_variable%name) // &
          '" not recognized in regression suite.'
        call printErrMsg(option)
    end select
  
    call OutputGetVarFromArray(realization,global_vec,ivar,isubvar)
    
    ! list of natural ids
    if (size(regression%natural_cell_ids) > 0) then
      call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                           global_vec, &
                           regression%natural_cell_id_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                         global_vec, &
                         regression%natural_cell_id_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
    endif
    if (regression%num_cells_per_process > 0) then
      ! cells per process
      call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                           global_vec, &
                           regression%cells_per_process_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                         global_vec, &
                         regression%cells_per_process_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
    endif

100 format(i9,1x,es20.13)    
    
    if (option%myrank == option%io_rank) then
      write(OUTPUT_UNIT,'(''-- '',a,'' --'')') trim(cur_variable%name)
      ! natural cell ids
      call VecGetArrayF90(regression%natural_cell_id_vec,vec_ptr,ierr)
      do i = 1, size(regression%natural_cell_ids)
        write(OUTPUT_UNIT,100) &
          regression%natural_cell_ids(i),vec_ptr(i)
      enddo
      call VecRestoreArrayF90(regression%natural_cell_id_vec,vec_ptr,ierr)

      ! cell ids per process
      call VecGetArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
      do i = 1, regression%num_cells_per_process*option%mycommsize
        write(OUTPUT_UNIT,100) &
          regression%cells_per_process_natural_ids(i),vec_ptr(i)
      enddo
      call VecRestoreArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
    endif
  
    cur_variable => cur_variable%next
  enddo
  
  call VecDestroy(global_vec,ierr)
  close(OUTPUT_UNIT)
  
end subroutine RegressionOutput

! ************************************************************************** !
!
! RegressionVariableDestroy: Destroys a regression variable object
! author: Glenn Hammond
! date: 10/11/12
!
! ************************************************************************** !
recursive subroutine RegressionVariableDestroy(regression_variable)

  implicit none
  
  type(regression_variable_type), pointer :: regression_variable
  
  if (.not.associated(regression_variable)) return
  
  call RegressionVariableDestroy(regression_variable%next)

  deallocate(regression_variable)
  nullify(regression_variable)
  
end subroutine RegressionVariableDestroy

! ************************************************************************** !
!
! RegressionDestroy: Destroys a regression object
! author: Glenn Hammond
! date: 10/11/12
!
! ************************************************************************** !
subroutine RegressionDestroy(regression)

  use Utility_module
  
  implicit none
  
  type(regression_type), pointer :: regression
  
  PetscErrorCode :: ierr
  
  if (.not.associated(regression)) return
  
  call RegressionVariableDestroy(regression%variable_list)
  call DeallocateArray(regression%natural_cell_ids)
  regression%num_cells_per_process = 0
  call DeallocateArray(regression%cells_per_process_natural_ids)
  if (regression%natural_cell_id_vec /= 0) &
    call VecDestroy(regression%natural_cell_id_vec,ierr)
  if (regression%cells_per_process_vec /= 0) &
    call VecDestroy(regression%cells_per_process_vec,ierr)
  if (regression%scatter_natural_cell_id_gtos /= 0) &
    call VecScatterDestroy(regression%scatter_natural_cell_id_gtos,ierr)
  if (regression%scatter_cells_per_process_gtos /= 0) &
    call VecScatterDestroy(regression%scatter_cells_per_process_gtos,ierr)

  deallocate(regression)
  nullify(regression)
  
end subroutine RegressionDestroy

end module Regression_module
