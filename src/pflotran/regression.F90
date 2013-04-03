module Regression_module
 
  use Output_Aux_module
  
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
  use Realization_class
  use Grid_module
  use Discretization_module
  
  implicit none
  
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(regression_type), pointer :: regression
  type(realization_type) :: realization
  
  IS :: is_petsc
  PetscInt, allocatable :: int_array(:)
  PetscInt :: i, upper_bound, lower_bound, count, temp_int
  PetscInt :: local_id
  PetscReal, pointer :: vec_ptr(:)
  Vec :: temp_vec
  VecScatter :: temp_scatter
  IS :: temp_is
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  if (.not.associated(regression)) return
  
  grid => realization%patch%grid
  option => realization%option
  
  ! natural cell ids
  if (associated(regression%natural_cell_ids)) then
    ! ensure that natural ids are within problem domain
    if (maxval(regression%natural_cell_ids) > grid%nmax) then
      option%io_buffer = 'Natural IDs outside problem domain requested ' // &
        'for regression output.  Removing non-existent IDs.'
      call printWrnMsg(option)
      count = 0
      allocate(int_array(size(regression%natural_cell_ids)))
      do i = 1, size(regression%natural_cell_ids)
        if (regression%natural_cell_ids(i) <= grid%nmax) then
          count = count + 1
          int_array(count) = regression%natural_cell_ids(i)
        endif
      enddo
      ! reallocate array
      deallocate(regression%natural_cell_ids)
      allocate(regression%natural_cell_ids(count))
      regression%natural_cell_ids = int_array
      deallocate(int_array)
    endif
    call VecCreate(PETSC_COMM_SELF,regression%natural_cell_id_vec,ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(regression%natural_cell_id_vec, &
                       size(regression%natural_cell_ids), &
                       PETSC_DECIDE,ierr)  
    else
      call VecSetSizes(regression%natural_cell_id_vec,0, &
                       PETSC_DECIDE,ierr)  
    endif
    call VecSetFromOptions(regression%natural_cell_id_vec,ierr)
  
    if (option%myrank == option%io_rank) then
      count = size(regression%natural_cell_ids)
      ! determine how many of the natural cell ids are local
      allocate(int_array(count))
      int_array = regression%natural_cell_ids
      ! convert to zero based
      int_array = int_array - 1
    else
      count = 0
      allocate(int_array(count))
    endif
    call DiscretAOApplicationToPetsc(realization%discretization,int_array)
  
  ! create IS for global petsc cell ids
    call ISCreateGeneral(option%mycomm,count,int_array,PETSC_COPY_VALUES, &
                         is_petsc,ierr)
    deallocate(int_array)
  
#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'is_petsc_natural_cell_id.out', &
                              viewer,ierr)
    call ISView(is_petsc,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! create scatter context
    call VecScatterCreate(realization%field%work,is_petsc, &
                          regression%natural_cell_id_vec,PETSC_NULL_OBJECT, &
                          regression%scatter_natural_cell_id_gtos,ierr)

    call ISDestroy(is_petsc,ierr)

#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'regression_scatter_nat_cell_ids.out',viewer,ierr)
    call VecScatterView(regression%scatter_natural_cell_id_gtos,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

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
    call VecCreate(PETSC_COMM_SELF,regression%cells_per_process_vec,ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(regression%cells_per_process_vec, &
                       regression%num_cells_per_process*option%mycommsize, &
                       PETSC_DECIDE,ierr)  
    else
      call VecSetSizes(regression%cells_per_process_vec,ZERO_INTEGER, &
                       PETSC_DECIDE,ierr)  
    endif
    call VecSetFromOptions(regression%cells_per_process_vec,ierr)

    ! create temporary vec to transfer down ids of cells
    call VecCreate(option%mycomm,temp_vec,ierr)
    call VecSetSizes(temp_vec, &
                     regression%num_cells_per_process, &
                     PETSC_DECIDE,ierr)  
    call VecSetFromOptions(temp_vec,ierr)
  
    ! calculate interval
    call VecGetArrayF90(temp_vec,vec_ptr,ierr)
    temp_int = grid%nlmax / regression%num_cells_per_process
    do i = 1, regression%num_cells_per_process
      vec_ptr(i) = temp_int*(i-1) + 1 + grid%global_offset
    enddo
    call VecRestoreArrayF90(temp_vec,vec_ptr,ierr)

    ! create temporary scatter to transfer values to io_rank
    if (option%myrank == option%io_rank) then
      count = option%mycommsize*regression%num_cells_per_process
      ! determine how many of the natural cell ids are local
      allocate(int_array(count))
      do i = 1, count
        int_array(i) = i
      enddo
      ! convert to zero based
      int_array = int_array - 1
    else
      count = 0
      allocate(int_array(count))
    endif
    call ISCreateGeneral(option%mycomm,count, &
                         int_array,PETSC_COPY_VALUES,temp_is,ierr)

    call VecScatterCreate(temp_vec,temp_is, &
                          regression%cells_per_process_vec,PETSC_NULL_OBJECT, &
                          temp_scatter,ierr)
    call ISDestroy(temp_is,ierr)
 
    ! scatter ids to io_rank
    call VecScatterBegin(temp_scatter,temp_vec, &
                         regression%cells_per_process_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(temp_scatter,temp_vec, &
                       regression%cells_per_process_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterDestroy(temp_scatter,ierr) 
    call VecDestroy(temp_vec,ierr)
   
    ! transfer cell ids into array for creating new scatter
    if (option%myrank == option%io_rank) then
      count = option%mycommsize*regression%num_cells_per_process
      call VecGetArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
      do i = 1, count
        int_array(i) = int(vec_ptr(i)+0.1d0) ! tolerance to ensure int value
      enddo
      call VecRestoreArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
      ! convert to zero based
      int_array = int_array - 1
    endif

    call ISCreateGeneral(option%mycomm,count, &
                         int_array,PETSC_COPY_VALUES,is_petsc,ierr)
    deallocate(int_array)

#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'is_petsc_cells_per_process.out', &
                              viewer,ierr)
    call ISView(is_petsc,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    call VecScatterCreate(realization%field%work,is_petsc, &
                          regression%cells_per_process_vec, &
                          PETSC_NULL_OBJECT, &
                          regression%scatter_cells_per_process_gtos,ierr)
    call ISDestroy(is_petsc,ierr)

#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'regression_scatter_cells_per_process.out', &
                              viewer,ierr)
    call VecScatterView(regression%scatter_cells_per_process_gtos,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif
  
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

  use Realization_class
  use Timestepper_module
  use Option_module
  use Discretization_module
  use Output_module
  use Output_Aux_module
  use Output_Common_module, only : OutputGetCellCenteredVelocities, &
                                   OutputGetVarFromArray
  
  implicit none
  
  type(regression_type), pointer :: regression
  type(realization_type) :: realization
  ! these must be pointers as they can be null
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper  
  
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  Vec :: x_vel_natural, y_vel_natural, z_vel_natural
  Vec :: x_vel_process, y_vel_process, z_vel_process
  PetscInt :: ivar, isubvar
  type(option_type), pointer :: option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:), y_ptr(:), z_ptr(:)
  PetscInt :: i
  PetscInt :: iphase
  PetscReal :: r_norm, x_norm
  PetscReal :: max, min, mean
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
  
  cur_variable => realization%output_option%output_variable_list%first
  do 
    if (.not.associated(cur_variable)) exit
    
    ivar = cur_variable%ivar
    isubvar = cur_variable%isubvar
  
    call OutputGetVarFromArray(realization,global_vec,ivar,isubvar)
    
    call VecMax(global_vec,PETSC_NULL_INTEGER,max,ierr)
    call VecMin(global_vec,PETSC_NULL_INTEGER,min,ierr)
    call VecSum(global_vec,mean,ierr)
    mean = mean / realization%patch%grid%nmax
    
    ! list of natural ids
    if (associated(regression%natural_cell_ids)) then
      call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                           global_vec, &
                           regression%natural_cell_id_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
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

100 format(i9,': ',es21.13)    
101 format(i9,': ',i9)    
    
    if (option%myrank == option%io_rank) then
      string = OutputVariableToCategoryString(cur_variable%icategory)
      write(OUTPUT_UNIT,'(''-- '',a,'': '',a,'' --'')') &
        trim(string), trim(cur_variable%name)
      
      ! max, min, mean
      if (cur_variable%iformat == 0) then
        write(OUTPUT_UNIT,'(6x,''Max: '',es21.13)') max
        write(OUTPUT_UNIT,'(6x,''Min: '',es21.13)') min
      else
        write(OUTPUT_UNIT,'(6x,''Max: '',i9)') int(max)
        write(OUTPUT_UNIT,'(6x,''Min: '',i9)') int(min)
      endif
      write(OUTPUT_UNIT,'(5x,''Mean: '',es21.13)') mean
      
      ! natural cell ids
      if (associated(regression%natural_cell_ids)) then
        if (size(regression%natural_cell_ids) > 0) then
          call VecGetArrayF90(regression%natural_cell_id_vec,vec_ptr,ierr)
          if (cur_variable%iformat == 0) then
            do i = 1, size(regression%natural_cell_ids)
              write(OUTPUT_UNIT,100) &
                regression%natural_cell_ids(i),vec_ptr(i)
            enddo
          else
            do i = 1, size(regression%natural_cell_ids)
              write(OUTPUT_UNIT,101) &
                regression%natural_cell_ids(i),int(vec_ptr(i))
            enddo
          endif
          call VecRestoreArrayF90(regression%natural_cell_id_vec,vec_ptr,ierr)
        endif
      endif
      
      ! cell ids per process
      if (regression%num_cells_per_process > 0) then
        call VecGetArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
        if (cur_variable%iformat == 0) then
          do i = 1, regression%num_cells_per_process*option%mycommsize
            write(OUTPUT_UNIT,100) &
              regression%cells_per_process_natural_ids(i),vec_ptr(i)
          enddo
        else
          do i = 1, regression%num_cells_per_process*option%mycommsize
            write(OUTPUT_UNIT,101) &
              regression%cells_per_process_natural_ids(i),int(vec_ptr(i))
          enddo
        endif
        call VecRestoreArrayF90(regression%cells_per_process_vec,vec_ptr,ierr)
      endif
    endif
  
    cur_variable => cur_variable%next
  enddo
  
  ! velocities
  if ((realization%output_option%print_tecplot_velocities .or. &
       realization%output_option%print_hdf5_velocities) .and. &
      option%nflowdof > 0) then
    if (associated(regression%natural_cell_ids)) then
      call VecDuplicate(regression%natural_cell_id_vec,x_vel_natural,ierr)
      call VecDuplicate(regression%natural_cell_id_vec,y_vel_natural,ierr)
      call VecDuplicate(regression%natural_cell_id_vec,z_vel_natural,ierr)
      call VecZeroEntries(x_vel_natural,ierr)
      call VecZeroEntries(y_vel_natural,ierr)
      call VecZeroEntries(z_vel_natural,ierr)
    endif
    if (regression%num_cells_per_process > 0) then
      call VecDuplicate(regression%cells_per_process_vec,x_vel_process,ierr)
      call VecDuplicate(regression%cells_per_process_vec,y_vel_process,ierr)
      call VecDuplicate(regression%cells_per_process_vec,z_vel_process,ierr)
      call VecZeroEntries(x_vel_process,ierr)
      call VecZeroEntries(y_vel_process,ierr)
      call VecZeroEntries(z_vel_process,ierr)
    endif

    do iphase = 1, option%nphase
      if (associated(regression%natural_cell_ids) .or. &
          regression%num_cells_per_process > 0) then
    
        if (iphase == 1) then
          string = 'LIQUID'
        else
          string = 'GAS'
        endif
        if (option%myrank == option%io_rank) then
          write(OUTPUT_UNIT,'(''-- GENERIC: '',a,'' VELOCITY ['',a, &
                              &''] --'')') &
            trim(string), 'm/' // trim(realization%output_option%tunit)
        endif
    
        ! X
        call OutputGetCellCenteredVelocities(realization,global_vec,iphase, &
                                             X_DIRECTION)
        if (associated(regression%natural_cell_ids)) then
          call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                               global_vec,x_vel_natural,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr)
          call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
                             global_vec,x_vel_natural,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr)
        endif
        if (regression%num_cells_per_process > 0) then
          call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                               global_vec,x_vel_process,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr)
          call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                             global_vec,x_vel_process,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr)
        endif
        ! Y
        call OutputGetCellCenteredVelocities(realization,global_vec,iphase, &
                                             Y_DIRECTION)
        if (associated(regression%natural_cell_ids)) then
          call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                               global_vec,y_vel_natural,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr)
          call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
                             global_vec,y_vel_natural,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr)
        endif
        if (regression%num_cells_per_process > 0) then
          call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                               global_vec,y_vel_process,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr)
          call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                             global_vec,y_vel_process,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr)
        endif
        ! Z
        call OutputGetCellCenteredVelocities(realization,global_vec,iphase, &
                                             Z_DIRECTION)
        if (associated(regression%natural_cell_ids)) then
          call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                               global_vec,z_vel_natural,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr)
          call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
                             global_vec,z_vel_natural,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr)
        endif
        if (regression%num_cells_per_process > 0) then
          call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                               global_vec,z_vel_process,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr)
          call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                             global_vec,z_vel_process,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr)
        endif
      
104 format(i9,': ',3es21.13) 

        ! natural cell ids
        if (option%myrank == option%io_rank) then
          if (associated(regression%natural_cell_ids)) then
            if (size(regression%natural_cell_ids) > 0) then
              call VecGetArrayF90(x_vel_natural,vec_ptr,ierr)
              call VecGetArrayF90(y_vel_natural,y_ptr,ierr)
              call VecGetArrayF90(z_vel_natural,z_ptr,ierr)
              do i = 1, size(regression%natural_cell_ids)
                write(OUTPUT_UNIT,104) &
                  regression%natural_cell_ids(i),vec_ptr(i),y_ptr(i),z_ptr(i)
              enddo
              call VecRestoreArrayF90(x_vel_natural,vec_ptr,ierr)
              call VecRestoreArrayF90(y_vel_natural,y_ptr,ierr)
              call VecRestoreArrayF90(z_vel_natural,z_ptr,ierr)
            endif
          endif
      
          ! cell ids per process
          if (regression%num_cells_per_process > 0) then
            call VecGetArrayF90(x_vel_process,vec_ptr,ierr)
            call VecGetArrayF90(y_vel_process,y_ptr,ierr)
            call VecGetArrayF90(z_vel_process,z_ptr,ierr)
            do i = 1, regression%num_cells_per_process*option%mycommsize
              write(OUTPUT_UNIT,104) &
                regression%cells_per_process_natural_ids(i),vec_ptr(i), &
                  y_ptr(i),z_ptr(i)
            enddo
            call VecRestoreArrayF90(x_vel_process,vec_ptr,ierr)
            call VecRestoreArrayF90(y_vel_process,y_ptr,ierr)
            call VecRestoreArrayF90(z_vel_process,z_ptr,ierr)
          endif
        endif
      endif
    enddo

    if (associated(regression%natural_cell_ids)) then
      call VecDestroy(x_vel_natural,ierr)
      call VecDestroy(y_vel_natural,ierr)
      call VecDestroy(z_vel_natural,ierr)
    endif
    if (regression%num_cells_per_process > 0) then
      call VecDestroy(x_vel_process,ierr)
      call VecDestroy(y_vel_process,ierr)
      call VecDestroy(z_vel_process,ierr)
    endif
  endif ! option%nflowdof > 0
  
  call VecDestroy(global_vec,ierr)
  
102 format(i12)    
103 format(es21.13)

  ! timestep, newton iteration, solver iteration output
  if (associated(flow_stepper)) then
    call VecNorm(realization%field%flow_xx,NORM_2,x_norm,ierr)
    call VecNorm(realization%field%flow_r,NORM_2,r_norm,ierr)
    if (option%myrank == option%io_rank) then
      write(OUTPUT_UNIT,'(''-- SOLUTION: Flow --'')')
      write(OUTPUT_UNIT,'(''   Time (seconds): '',es21.13)') &
        flow_stepper%cumulative_solver_time
      write(OUTPUT_UNIT,'(''   Time Steps: '',i12)') flow_stepper%steps
      write(OUTPUT_UNIT,'(''   Newton Iterations: '',i12)') &
        flow_stepper%cumulative_newton_iterations
      write(OUTPUT_UNIT,'(''   Solver Iterations: '',i12)') &
        flow_stepper%cumulative_linear_iterations
      write(OUTPUT_UNIT,'(''   Time Step Cuts: '',i12)') &
        flow_stepper%cumulative_time_step_cuts
      write(OUTPUT_UNIT,'(''   Solution 2-Norm: '',es21.13)') x_norm
      write(OUTPUT_UNIT,'(''   Residual 2-Norm: '',es21.13)') r_norm
    endif
  endif
  if (associated(tran_stepper)) then
    call VecNorm(realization%field%tran_xx,NORM_2,x_norm,ierr)
    if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
      call VecNorm(realization%field%tran_r,NORM_2,r_norm,ierr)
    endif
    if (option%myrank == option%io_rank) then
      write(OUTPUT_UNIT,'(''-- SOLUTION: Transport --'')')
      write(OUTPUT_UNIT,'(''   Time (seconds): '',es21.13)') &
        tran_stepper%cumulative_solver_time
      write(OUTPUT_UNIT,'(''   Time Steps: '',i12)') tran_stepper%steps
      write(OUTPUT_UNIT,'(''   Newton Iterations: '',i12)') &
        tran_stepper%cumulative_newton_iterations
      write(OUTPUT_UNIT,'(''   Solver Iterations: '',i12)') &
        tran_stepper%cumulative_linear_iterations
      write(OUTPUT_UNIT,'(''   Time Step Cuts: '',i12)') &
        tran_stepper%cumulative_time_step_cuts
      write(OUTPUT_UNIT,'(''   Solution 2-Norm: '',es21.13)') x_norm
      if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
        write(OUTPUT_UNIT,'(''   Residual 2-Norm: '',es21.13)') r_norm
      endif
    endif
  endif
  
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
