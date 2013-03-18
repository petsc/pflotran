module Output_HDF5_module

  use Logging_module 
  use Output_Aux_module
  use Output_Common_module
  
  implicit none

  private

#include "definitions.h"

#if defined(SCORPIO_WRITE)
  include "scorpiof.h"
#endif

  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: hdf5_first
  
  public :: OutputHDF5Init, &
            OutputHDF5, &
            OutputHDF5UGridXDMF

contains

! ************************************************************************** !
!
! OutputHDF5Init: Initializes module variables for HDF5 output
! author: Glenn Hammond
! date: 01/16/13
!
! ************************************************************************** !
subroutine OutputHDF5Init(realization_base,num_steps)

  use Realization_Base_class, only : realization_base_type
  use Option_module

  implicit none
  
  class(realization_base_type) :: realization_base
  PetscInt :: num_steps
  
  if (num_steps == 0) then
    hdf5_first = PETSC_TRUE
  else
    hdf5_first = PETSC_FALSE
  endif

end subroutine OutputHDF5Init

! ************************************************************************** !
!
! OutputHDF5: Print to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputHDF5(realization_base,var_list_type)

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use String_module
  
#if !defined(PETSC_HAVE_HDF5)
  implicit none
  
  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type  

  call printMsg(realization_base%option,'')
  write(realization_base%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to &
        &write HDF5 formatted structured grids Darn.")')
  call printErrMsg(realization_base%option)
#else

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use hdf5
  use HDF5_module
  use HDF5_Aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type

#if defined(SCORPIO_WRITE)
  integer:: file_id
  integer:: grp_id
  integer:: file_space_id
  integer:: realization_set_id
  integer:: prop_id
  PetscMPIInt :: rank
  integer:: dims(3)
  integer:: pio_dataset_groupid
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  integer(HSIZE_T) :: dims(3)
#endif
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  
  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string,string2,string3
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_err
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      string2=''
      write(string3,'(i4)') output_option%plot_number
    case (AVERAGED_VARS)
      string2='-aveg'
      write(string3,'(i4)') int(option%time/output_option%periodic_output_time_incr)
  end select

  if (output_option%print_single_h5_file) then
    first = hdf5_first
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               trim(string2) // '.h5'
  else
    string = OutputHDF5FilenameID(output_option,option,var_list_type)
    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        if (mod(output_option%plot_number,output_option%times_per_h5_file)==0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
      case (AVERAGED_VARS)
        if (mod((option%time-output_option%periodic_output_time_incr)/ &
                output_option%periodic_output_time_incr, &
                real(output_option%times_per_h5_file))==0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
    end select

    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                '-' // trim(string) // trim(string2) // '.h5'
  endif

  grid => patch%grid
#if defined(SCORPIO_WRITE)
  if (.not.first) then
    filename = trim(filename) // CHAR(0)
    call scorpio_open_file(filename, option%iowrite_group_id, &
                              SCORPIO_FILE_READWRITE, file_id, ierr)
    if (file_id == -1) first = PETSC_TRUE
  endif
  if (first) then
    filename = trim(filename) // CHAR(0)
    call scorpio_open_file(filename, option%iowrite_group_id, &
                              SCORPIO_FILE_CREATE, file_id, ierr)
  endif

#else
! SCORPIO_WRITE is not defined

    ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  if (.not.first) then
    call h5eset_auto_f(OFF,hdf5_err)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
    if (hdf5_err /= 0) first = PETSC_TRUE
    call h5eset_auto_f(ON,hdf5_err)
  endif
  if (first) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                      H5P_DEFAULT_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)
#endif
! SCORPIO_WRITE

  if (first) then
    option%io_buffer = '--> creating hdf5 output file: ' // trim(filename)
  else
    option%io_buffer = '--> appending to hdf5 output file: ' // trim(filename)
  endif
  call printMsg(option)

  if (first) then

    ! create a group for the coordinates data set
#if defined(SCORPIO_WRITE)
    string = "Coordinates" // CHAR(0)
    call scorpio_create_dataset_group(pio_dataset_groupid, string, file_id, &
                                        option%iowrite_group_id, ierr)
        ! set grp_id here
        ! As we already created the group, we will use file_id as group_id
    grp_id = file_id
#else
    string = "Coordinates"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
#endif

    !GEH - Structured Grid Dependence - Begin
    ! write out coordinates in x, y, and z directions
    string = "X [m]"
    allocate(array(grid%structured_grid%nx+1))
    array(1) = grid%structured_grid%origin(X_DIRECTION)
    do i=2,grid%structured_grid%nx+1
      array(i) = array(i-1) + grid%structured_grid%dx_global(i-1)
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%nx+1,array,grp_id)
    deallocate(array)

    string = "Y [m]"
    allocate(array(grid%structured_grid%ny+1))
    array(1) = grid%structured_grid%origin(Y_DIRECTION)
    do i=2,grid%structured_grid%ny+1
      array(i) = array(i-1) + grid%structured_grid%dy_global(i-1)
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%ny+1,array,grp_id)
    deallocate(array)

    string = "Z [m]"
    allocate(array(grid%structured_grid%nz+1))
    array(1) = grid%structured_grid%origin(Z_DIRECTION)
    do i=2,grid%structured_grid%nz+1
      array(i) = array(i-1) + grid%structured_grid%dz_global(i-1)
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%nz+1,array,grp_id)
    deallocate(array)
    !GEH - Structured Grid Dependence - End

#if defined(SCORPIO_WRITE)
    call scorpio_close_dataset_group(pio_dataset_groupid, file_id, &
                                        option%iowrite_group_id, ierr)
#else
    call h5gclose_f(grp_id,hdf5_err)
#endif

  endif
        
  ! create a group for the data set
  write(string,'(''Time:'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  !string = trim(string3) // ' ' // trim(string)
#if defined(SCORPIO_WRITE)
  string = trim(string) //CHAR(0)
    ! This opens existing dataset and creates it if needed
  call scorpio_create_dataset_group(pio_dataset_groupid, string, file_id, &
                                        option%iowrite_group_id, ierr)
  grp_id = file_id
#else
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)
#endif
! SCORPIO_WRITE
  
  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)


  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over variables and write to file
      cur_variable => output_option%output_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVarFromArray(realization_base,global_vec,cur_variable%ivar, &
                                   cur_variable%isubvar)
        string = cur_variable%name
        call StringSwapChar(string," ","_")
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                            global_vec,grp_id,H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                            global_vec,grp_id,H5T_NATIVE_INTEGER)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      cur_variable => output_option%aveg_output_variable_list%first
      do ivar = 1,output_option%aveg_output_variable_list%nvars
        string = 'Aveg. ' // cur_variable%name
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif

        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           field%avg_vars_vec(ivar),grp_id, &
                                           H5T_NATIVE_DOUBLE)

        cur_variable => cur_variable%next
      enddo

  end select

  if (output_option%print_hdf5_velocities.and.(var_list_type==INSTANTANEOUS_VARS)) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)
    call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    if (option%nphase > 1) then
        call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,X_DIRECTION)
        string = "Gas X-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,Y_DIRECTION)
        string = "Gas Y-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,Z_DIRECTION)
        string = "Gas Z-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)
    endif
  endif

  if (output_option%print_hdf5_flux_velocities.and.(var_list_type==INSTANTANEOUS_VARS)) then

    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
        string = "Liquid X-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE,X_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas X-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE,X_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%ny > 1) then
        string = "Liquid Y-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE,Y_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Y-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE,Y_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%nz > 1) then
        string = "Liquid Z-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE,Z_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Z-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE,Z_DIRECTION,grp_id)
        endif
    endif
   
  endif

  call VecDestroy(global_vec,ierr)

#if defined(SCORPIO_WRITE)
    call scorpio_close_dataset_group(pio_dataset_groupid, file_id, &
            option%iowrite_group_id, ierr)
    call scorpio_close_file(file_id, option%iowrite_group_id, ierr)
#else
    call h5gclose_f(grp_id,hdf5_err)
    call h5fclose_f(file_id,hdf5_err)
     call h5close_f(hdf5_err)
#endif
!SCORPIO_WRITE
#endif
!PETSC_HAVE_HDF5

  hdf5_first = PETSC_FALSE

end subroutine OutputHDF5

! ************************************************************************** !
!> This subroutine prints a HDF5 file.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/31/12
! ************************************************************************** !
subroutine OutputHDF5UGrid(realization_base)

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use Variables_module
  use String_module

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  
  class(realization_base_type) :: realization_base

  call printMsg(realization_base%option,'')
  write(realization_base%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to &
        &write HDF5 formatted structured grids Darn.")')
  call printErrMsg(realization_base%option)
#else

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use hdf5
  use HDF5_module
  use HDF5_Aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"

  class(realization_base_type) :: realization_base

#if defined(SCORPIO_WRITE)
  integer:: file_id
  integer:: data_type
  integer:: grp_id
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: realization_set_id
  integer:: prop_id
  PetscMPIInt :: rank
  integer :: rank_mpi,file_space_rank_mpi
  integer:: dims(3)
  integer :: start(3), length(3), stride(3),istart
  integer :: pio_dataset_groupid
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3),istart
#endif
  
  PetscMPIInt :: hdf5_flag

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  
  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscInt :: local_size
  PetscReal, pointer :: double_array(:)
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscMPIInt :: hdf5_err  
  PetscFortranAddr :: app_ptr
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr  

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option


  if (output_option%print_single_h5_file) then
    first = hdf5_first
    filename = trim(option%global_prefix) // trim(option%group_prefix) // '.h5'
  else
    string = OutputFilenameID(output_option,option)
    first = PETSC_TRUE
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                '-' // trim(string) // '.h5'
  endif

  grid => patch%grid

#if defined(SCORPIO_WRITE)

  if (.not.first) then
    filename = trim(filename) // CHAR(0)
    call scorpio_open_file(filename, option%iowrite_group_id, &
                              SCORPIO_FILE_READWRITE, file_id, ierr)
    if (file_id == -1) first = PETSC_TRUE
  endif
  if (first) then
    filename = trim(filename) // CHAR(0)
    call scorpio_open_file(filename, option%iowrite_group_id, &
                              SCORPIO_FILE_CREATE, file_id, ierr)
  endif

  if (first) then
    option%io_buffer = '--> creating hdf5 output file: ' // trim(filename)
  else
    option%io_buffer = '--> appending to hdf5 output file: ' // trim(filename)
  endif
  call printMsg(option)

  if (first) then
    ! create a group for the coordinates data set
    string = "Domain" // CHAR(0)
    call scorpio_create_dataset_group(pio_dataset_groupid, string, file_id, &
                                         option%iowrite_group_id, ierr)
    ! set grp_id here
    ! As we already created the group, we will use file_id as group_id
    grp_id = file_id
    call WriteHDF5CoordinatesUGrid(grid,option,grp_id)
    call scorpio_close_dataset_group(pio_dataset_groupid, file_id, &
            option%iowrite_group_id, ierr)
  endif

#else

  !
  !        not(SCORPIO_WRITE)
  !

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  if (.not.first) then
    call h5eset_auto_f(OFF,hdf5_err)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
    if (hdf5_err /= 0) first = PETSC_TRUE
    call h5eset_auto_f(ON,hdf5_err)
  endif
  if (first) then
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                      H5P_DEFAULT_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (first) then
    option%io_buffer = '--> creating hdf5 output file: ' // trim(filename)
  else
    option%io_buffer = '--> appending to hdf5 output file: ' // trim(filename)
  endif
  call printMsg(option)

  if (first) then
    ! create a group for the coordinates data set
    string = "Domain"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
    call WriteHDF5CoordinatesUGrid(grid,option,grp_id)
    call h5gclose_f(grp_id,hdf5_err)
  endif
#endif
! SCORPIO_WRITE

    ! create a group for the data set
    write(string,'(''Time:'',es13.5,x,a1)') &
          option%time/output_option%tconv,output_option%tunit
    if (len_trim(output_option%plot_name) > 2) then
      string = trim(string) // ' ' // output_option%plot_name
    endif
#if defined(SCORPIO_WRITE)
    string = trim(string) //CHAR(0)
      ! This opens existing dataset and creates it if needed
    call scorpio_create_dataset_group(pio_dataset_groupid, string, file_id, &
                                          option%iowrite_group_id, ierr)
    grp_id = file_id
#else
    call h5eset_auto_f(OFF,hdf5_err)
    call h5gopen_f(file_id,string,grp_id,hdf5_err)
    if (hdf5_err /= 0) then
      call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
    endif
    call h5eset_auto_f(ON,hdf5_err)
#endif
! SCORPIO_WRITE

  ! write out data sets
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)

  ! loop over variables and write to file
  cur_variable => output_option%output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputGetVarFromArray(realization_base,global_vec,cur_variable%ivar, &
                                cur_variable%isubvar)
    string = cur_variable%name
    call StringSwapChar(string," ","_")
    if (len_trim(cur_variable%units) > 0) then
      word = cur_variable%units
      call HDF5MakeStringCompatible(word)
      string = trim(string) // ' [' // trim(word) // ']'
    endif
    if (cur_variable%iformat == 0) then
      call HDF5WriteUnstructuredDataSetFromVec(string,option, &
                                         global_vec,grp_id,H5T_NATIVE_DOUBLE)
    else
      call HDF5WriteUnstructuredDataSetFromVec(string,option, &
                                         global_vec,grp_id,H5T_NATIVE_INTEGER)
    endif
    cur_variable => cur_variable%next
  enddo

  if (output_option%print_hdf5_velocities) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call HDF5WriteUnstructuredDataSetFromVec(string,option,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)
    call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call HDF5WriteUnstructuredDataSetFromVec(string,option,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call HDF5WriteUnstructuredDataSetFromVec(string,option,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    if (option%nphase > 1) then
        call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,X_DIRECTION)
        string = "Gas X-Velocity"
        call HDF5WriteUnstructuredDataSetFromVec(string,option,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,Y_DIRECTION)
        string = "Gas Y-Velocity"
        call HDF5WriteUnstructuredDataSetFromVec(string,option,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,Z_DIRECTION)
        string = "Gas Z-Velocity"
        call HDF5WriteUnstructuredDataSetFromVec(string,option,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)
    endif
  endif

  if (output_option%print_hdf5_flux_velocities) then

    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
        string = "Liquid X-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE,X_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas X-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE,X_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%ny > 1) then
        string = "Liquid Y-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE,Y_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Y-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE,Y_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%nz > 1) then
        string = "Liquid Z-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE,Z_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Z-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE,Z_DIRECTION,grp_id)
        endif
    endif
   
  endif

  call VecDestroy(global_vec,ierr)

#if defined(SCORPIO_WRITE)
!    call scorpio_close_dataset_group(pio_dataset_groupid, file_id, &
!            option%iowrite_group_id, ierr)
!    call scorpio_close_file(file_id, option%iowrite_group_id, ierr)
#else
  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif
!SCORPIO_WRITE

  hdf5_first = PETSC_FALSE

#endif
! !defined(PETSC_HAVE_HDF5)

end subroutine OutputHDF5UGrid

! ************************************************************************** !
!> This routine writes unstructured grid data in HDF5 XDMF format.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/29/2012
! ************************************************************************** !
subroutine OutputHDF5UGridXDMF(realization_base,var_list_type)

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module

#if  !defined(PETSC_HAVE_HDF5)

  implicit none
  
  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type

  call printMsg(realization_base%option,'')
  write(realization_base%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to &
        &write HDF5 formatted structured grids Darn.")')
  call printErrMsg(realization_base%option)

#else

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use hdf5
  use HDF5_module
  use HDF5_Aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type

#if defined(SCORPIO_WRITE)
  integer:: file_id
  integer:: data_type
  integer:: grp_id
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: realization_set_id
  integer:: prop_id
  PetscMPIInt :: rank
  integer :: rank_mpi,file_space_rank_mpi
  integer:: dims(3)
  integer :: start(3), length(3), stride(3),istart
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3),istart
#endif

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: xmf_filename, att_datasetname, group_name
  character(len=MAXSTRINGLENGTH) :: string, string2,string3
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_err  
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type
  PetscInt :: vert_count
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      string2=''
      write(string3,'(i4)') output_option%plot_number
      xmf_filename = OutputFilename(output_option,option,'xmf','')
    case (AVERAGED_VARS)
      string2='-aveg'
      write(string3,'(i4)') int(option%time/output_option%periodic_output_time_incr)
      xmf_filename = OutputFilename(output_option,option,'xmf','aveg')
  end select

  if (output_option%print_single_h5_file) then
    first = hdf5_first
    filename = trim(option%global_prefix) // trim(string2) // trim(option%group_prefix) // '.h5'
  else
    string = OutputHDF5FilenameID(output_option,option,var_list_type)
    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        if (mod(output_option%plot_number,output_option%times_per_h5_file)==0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
      case (AVERAGED_VARS)
        if (mod((option%time-output_option%periodic_output_time_incr)/ &
                output_option%periodic_output_time_incr, &
                real(output_option%times_per_h5_file))==0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
    end select

    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               trim(string2) // '-' // trim(string) // '.h5'
  endif

  grid => patch%grid

#ifdef SCORPIO_WRITE
   option%io_buffer='OutputHDF5UGridXDMF not supported with SCORPIO_WRITE'
   call printErrMsg(option)
#endif

    ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif

  if (.not.first) then
    call h5eset_auto_f(OFF,hdf5_err)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
    if (hdf5_err /= 0) first = PETSC_TRUE
    call h5eset_auto_f(ON,hdf5_err)
  endif
  if (first) then
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                     H5P_DEFAULT_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (first) then
    option%io_buffer = '--> creating hdf5 output file: ' // trim(filename)
  else
    option%io_buffer = '--> appending to hdf5 output file: ' // trim(filename)
  endif
  call printMsg(option)

  if (first) then
    ! create a group for the coordinates data set
    string = "Domain"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
    call WriteHDF5CoordinatesUGridXDMF(realization_base,option,grp_id)
    call h5gclose_f(grp_id,hdf5_err)
  endif

  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write xmf output file: ' // trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=xmf_filename,action="write")
    !call OutputXMFHeader(OUTPUT_UNIT,realization_base,filename)
    call OutputXMFHeader(OUTPUT_UNIT, &
                         option%time/output_option%tconv, &
                         grid%nmax, &
                         realization_base%output_option%xmf_vert_len, &
                         grid%unstructured_grid%num_vertices_global,filename)
  endif

  ! create a group for the data set
  write(string,'(''Time'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  string = trim(string3) // ' ' // trim(string)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  group_name=string
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over variables and write to file
      cur_variable => output_option%output_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVarFromArray(realization_base,global_vec,cur_variable%ivar, &
                                   cur_variable%isubvar)
        call DiscretizationGlobalToNatural(discretization,global_vec, &
                                           natural_vec,ONEDOF)
        string = cur_variable%name
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteUnstructuredDataSetFromVec(string,option, &
                                              natural_vec,grp_id,H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteUnstructuredDataSetFromVec(string,option, &
                                              natural_vec,grp_id,H5T_NATIVE_INTEGER)
        endif
        att_datasetname = trim(filename) // ":/" // trim(group_name) // "/" // trim(string)
        if (option%myrank == option%io_rank) then
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string,att_datasetname)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      cur_variable => output_option%aveg_output_variable_list%first
      do ivar = 1,output_option%aveg_output_variable_list%nvars
        string = 'Aveg. ' // cur_variable%name
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif

        call DiscretizationGlobalToNatural(discretization,field%avg_vars_vec(ivar), &
                                           natural_vec,ONEDOF)
        call HDF5WriteUnstructuredDataSetFromVec(string,option, &
                                           natural_vec,grp_id,H5T_NATIVE_DOUBLE)
        att_datasetname = trim(filename) // ":/" // trim(group_name) // "/" // trim(string)
        if (option%myrank == option%io_rank) then
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string,att_datasetname)
        endif
        cur_variable => cur_variable%next
      enddo

  end select

  call VecDestroy(global_vec,ierr)
  call VecDestroy(natural_vec,ierr)
  call h5gclose_f(grp_id,hdf5_err)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

  if (option%myrank == option%io_rank) then
    call OutputXMFFooter(OUTPUT_UNIT)
    close(OUTPUT_UNIT)
  endif

  hdf5_first = PETSC_FALSE

#endif
! !defined(PETSC_HAVE_HDF5)

end subroutine OutputHDF5UGridXDMF

#if defined(PETSC_HAVE_HDF5)
! ************************************************************************** !
!> This subroutine creates an ID for HDF5 filename for:
!! - Instantaneous, or
!! - Temporally averaged variables.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/10/13
! ************************************************************************** !  
function OutputHDF5FilenameID(output_option,option,var_list_type)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option
  PetscInt :: var_list_type

  character(len=MAXWORDLENGTH) :: OutputHDF5FilenameID
  PetscInt :: file_number

  select case(var_list_type)
    case (INSTANTANEOUS_VARS)
      file_number = floor(real(output_option%plot_number)/ &
                               output_option%times_per_h5_file)
    case (AVERAGED_VARS)
      file_number = floor((option%time - &
                           output_option%periodic_output_time_incr)/ &
                          output_option%periodic_output_time_incr/ &
                          output_option%times_per_h5_file)
  end select

  if (file_number < 10) then
    write(OutputHDF5FilenameID,'("00",i1)') file_number
  else if (output_option%plot_number < 100) then
    write(OutputHDF5FilenameID,'("0",i2)') file_number  
  else if (output_option%plot_number < 1000) then
    write(OutputHDF5FilenameID,'(i3)') file_number  
  else if (output_option%plot_number < 10000) then
    write(OutputHDF5FilenameID,'(i4)') file_number
  endif 
  
  OutputHDF5FilenameID = adjustl(OutputHDF5FilenameID)

end function OutputHDF5FilenameID

! ************************************************************************** !
!
! WriteHDF5FluxVelocities: Print flux velocities to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5FluxVelocities(name,realization_base,iphase,direction,file_id)

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  use hdf5
  use HDF5_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"

  character(len=32) :: name
  class(realization_base_type) :: realization_base
  PetscInt :: iphase
  PetscInt :: direction
  integer(HID_T) :: file_id

  PetscInt :: i, j, k
  PetscInt :: count, iconn
  PetscInt :: local_id, ghosted_id
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
    
  PetscReal, allocatable :: array(:)
  PetscReal, pointer :: vec_ptr(:)

  PetscBool, save :: trick_flux_vel_x = PETSC_FALSE
  PetscBool, save :: trick_flux_vel_y = PETSC_FALSE
  PetscBool, save :: trick_flux_vel_z = PETSC_FALSE

  Vec :: global_vec

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
    
  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option  

  ! in a few cases (i.e. for small test problems), some processors may
  ! have no velocities to print.  This results in zero-length arrays
  ! in collective H5Dwrite().  To avoid, we switch to independent
  ! H5Dwrite() and don't write from the zero-length procs. 
!GEH - Structured Grid Dependence - Begin
  if (hdf5_first) then
    trick_flux_vel_x = PETSC_FALSE
    trick_flux_vel_y = PETSC_FALSE
    trick_flux_vel_z = PETSC_FALSE
    
    nx_local = grid%structured_grid%nlx
    ny_local = grid%structured_grid%nly
    nz_local = grid%structured_grid%nlz
    if (grid%structured_grid%gxe-grid%structured_grid%lxe == 0) then
      nx_local = grid%structured_grid%nlx-1
    endif
    call MPI_Allreduce(nx_local,i,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (i == 0) trick_flux_vel_x = PETSC_TRUE
    if (grid%structured_grid%gye-grid%structured_grid%lye == 0) then
      ny_local = grid%structured_grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (j == 0) trick_flux_vel_y = PETSC_TRUE
    if (grid%structured_grid%gze-grid%structured_grid%lze == 0) then
      nz_local = grid%structured_grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (k == 0) trick_flux_vel_z = PETSC_TRUE
  endif

  nx_local = grid%structured_grid%nlx
  ny_local = grid%structured_grid%nly
  nz_local = grid%structured_grid%nlz
  nx_global = grid%structured_grid%nx
  ny_global = grid%structured_grid%ny
  nz_global = grid%structured_grid%nz

  select case(direction)
    case(X_DIRECTION)
      nx_global = grid%structured_grid%nx-1
      if (grid%structured_grid%gxe-grid%structured_grid%lxe == 0) then
        nx_local = grid%structured_grid%nlx-1
      endif
      if (trick_flux_vel_x) trick_hdf5 = PETSC_TRUE
    case(Y_DIRECTION)
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%gye-grid%structured_grid%lye == 0) then
        ny_local = grid%structured_grid%nly-1
      endif
      if (trick_flux_vel_y) trick_hdf5 = PETSC_TRUE
    case(Z_DIRECTION)
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%gze-grid%structured_grid%lze == 0) then
        nz_local = grid%structured_grid%nlz-1
      endif
      if (trick_flux_vel_z) trick_hdf5 = PETSC_TRUE
  end select  
  allocate(array(nx_local*ny_local*nz_local))


  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option) 
  call VecZeroEntries(global_vec,ierr)
  call VecGetArrayF90(global_vec,vec_ptr,ierr)
  
  ! place interior velocities in a vector
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      ghosted_id = cur_connection_set%id_up(iconn)
      local_id = grid%nG2L(ghosted_id) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      if (local_id <= 0 .or. &
          dabs(cur_connection_set%dist(direction,iconn)) < 0.99d0) cycle
      vec_ptr(local_id) = patch%internal_velocities(iphase,iconn)
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  count = 0
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        array(count) = vec_ptr(local_id) 
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
  
  call VecDestroy(global_vec,ierr)
  
  array(1:nx_local*ny_local*nz_local) = &  ! convert time units
    array(1:nx_local*ny_local*nz_local) * output_option%tconv

  call HDF5WriteStructuredDataSet(name,array,file_id,H5T_NATIVE_DOUBLE,option, &
                        nx_global,ny_global,nz_global, &
                        nx_local,ny_local,nz_local, &
                        grid%structured_grid%lxs,grid%structured_grid%lys,grid%structured_grid%lzs)
!GEH - Structured Grid Dependence - End

  deallocate(array)
  trick_hdf5 = PETSC_FALSE

end subroutine WriteHDF5FluxVelocities

! ************************************************************************** !
!
! WriteHDF5Coordinates: Writes structured coordinates to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5Coordinates(name,option,length,array,file_id)

  use hdf5
  use Option_module
  
  implicit none
  
#if defined(SCORPIO_WRITE)
  character(len=32) :: name
  type(option_type) :: option
  PetscInt :: length
  PetscReal :: array(:)
  integer:: file_id

  integer:: file_space_id
  integer:: data_set_id
  integer:: prop_id
  integer:: dims(3)
  PetscMPIInt :: rank
  integer:: globaldims(3)
#else
  character(len=32) :: name
  type(option_type) :: option
  PetscInt :: length
  PetscReal :: array(:)
  integer(HID_T) :: file_id
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  PetscMPIInt :: rank
#endif

  PetscMPIInt :: hdf5_flag
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  
  call PetscLogEventBegin(logging%event_output_coordinates_hdf5,ierr) 
#if defined(SCORPIO_WRITE)

  name = trim(name) // CHAR(0)
  ! write out grid structure
  rank = 1
  dims = 0
  globaldims = 0
  ! x-direction

  ! Only process 0 writes coordinates
  if (option%myrank == 0 ) then
     dims(1) = length
     globaldims(1) = length
  else
     dims(1) = 0
     globaldims(1) = length
  endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call scorpio_write_dataset(array, SCORPIO_DOUBLE, rank, globaldims, dims, &
       file_id, name, option%iowrite_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE, &
       ierr)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

#else
!SCORPIO_WRITE is not defined

  ! write out grid structure
  rank = 1
  dims = 0
  ! x-direction
  dims(1) = length
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims)
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err) ! must be independent and only from p0
#endif
  if (option%myrank == option%io_rank) then
     call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)     
     call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                    hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
     call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

#endif
! SCORPIO_WRITE

  call PetscLogEventEnd(logging%event_output_coordinates_hdf5,ierr) 

end subroutine WriteHDF5Coordinates

! ************************************************************************** !
!> This subroutine writes unstructured coordinates to HDF5 file
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/31/12
! ************************************************************************** !
subroutine WriteHDF5CoordinatesUGrid(grid,option,file_id)

  use hdf5
  use HDF5_module
  use Grid_module
  use Option_module
  use Unstructured_Grid_Aux_module
  use Variables_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

#if defined(SCORPIO_WRITE)
  integer:: file_id
  integer:: data_type
  integer:: grp_id
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: realization_set_id
  integer:: prop_id
  integer:: dims(3)
  integer :: start(3), length(3), stride(3),istart
  integer :: rank_mpi,file_space_rank_mpi
  integer :: hdf5_flag
  integer, parameter :: ON=1, OFF=0
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3),istart
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0
#endif

  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err  

  PetscInt :: local_size
  PetscInt :: i,j
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  PetscInt, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element
  PetscErrorCode :: ierr

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_x_vertex_vec,ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_y_vertex_vec,ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_z_vertex_vec,ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr)

  call GetVertexCoordinates(grid, global_x_vertex_vec,X_COORDINATE,option)
  call GetVertexCoordinates(grid, global_y_vertex_vec,Y_COORDINATE,option)
  call GetVertexCoordinates(grid, global_z_vertex_vec,Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr)

#if defined(SCORPIO_WRITE)
  write(*,*),'SCORPIO_WRITE'
  option%io_buffer = 'WriteHDF5CoordinatesUGrid not supported for SCORPIO_WRITE'
  call printErrMsg(option)
#else

  !
  !        not(SCORPIO_WRITE)
  !
   
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 3D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%num_vertices_global
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)
  
  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = 3
  
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  allocate(double_array(local_size*3))
  
  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo
  
  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

  deallocate(double_array)
  call h5pclose_f(prop_id,hdf5_err)

  
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)


#endif

  call VecRestoreArrayF90(global_x_vertex_vec,vec_x_ptr,ierr)
  call VecRestoreArrayF90(global_y_vertex_vec,vec_y_ptr,ierr)
  call VecRestoreArrayF90(global_z_vertex_vec,vec_z_ptr,ierr)


  call VecDestroy(global_x_vertex_vec,ierr)
  call VecDestroy(global_y_vertex_vec,ierr)
  call VecDestroy(global_z_vertex_vec,ierr)


  !
  !  Write elements
  !
  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call GetCellConnections(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)

  local_size = grid%unstructured_grid%nlmax
#if defined(SCORPIO_WRITE)
  write(*,*),'SCORPIO_WRITE'
  option%io_buffer = 'WriteHDF5CoordinatesUGrid not supported for SCORPIO_WRITE'
  call printErrMsg(option)
#else

  !
  !        not(SCORPIO_WRITE)
  !
   
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size*NINE_INTEGER
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 3D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%nmax
  dims(1) = NINE_INTEGER
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then 
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)
  
  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = NINE_INTEGER
  
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  allocate(int_array(local_size*NINE_INTEGER))
  
  do i=1,local_size
    int_array((i-1)*9 + 1) = 0
    int_array((i-1)*9 + 2) = INT(vec_ptr((i-1)*8+1))
    int_array((i-1)*9 + 3) = INT(vec_ptr((i-1)*8+2))
    int_array((i-1)*9 + 4) = INT(vec_ptr((i-1)*8+3))
    int_array((i-1)*9 + 5) = INT(vec_ptr((i-1)*8+4))
    int_array((i-1)*9 + 6) = INT(vec_ptr((i-1)*8+5))
    int_array((i-1)*9 + 7) = INT(vec_ptr((i-1)*8+6))
    int_array((i-1)*9 + 8) = INT(vec_ptr((i-1)*8+7))
    int_array((i-1)*9 + 9) = INT(vec_ptr((i-1)*8+8))
    do j=2,9
      if(int_array((i-1)*9 + j)>0) int_array((i-1)*9 + 1)= int_array((i-1)*9 + 1) +1
    enddo
  enddo
  
  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

  deallocate(int_array)
  call h5pclose_f(prop_id,hdf5_err)

  
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

#endif

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)
  call VecDestroy(global_vec,ierr)
  call VecDestroy(natural_vec,ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine WriteHDF5CoordinatesUGrid

! ************************************************************************** !
!> This routine writes unstructured coordinates to HDF5 file in XDMF format
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/29/2012
! ************************************************************************** !
subroutine WriteHDF5CoordinatesUGridXDMF(realization_base,option,file_id)

  use hdf5
  use HDF5_module
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Unstructured_Grid_Aux_module
  use Variables_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petsclog.h"

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

#if defined(SCORPIO_WRITE)
  integer:: file_id
  integer:: data_type
  integer:: grp_id
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: realization_set_id
  integer:: prop_id
  integer:: dims(3)
  integer :: start(3), length(3), stride(3),istart
  integer :: rank_mpi,file_space_rank_mpi
  integer :: hdf5_flag
  integer, parameter :: ON=1, OFF=0
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3),istart
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0
#endif

  type(grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err  

  PetscInt :: local_size,vert_count,nverts
  PetscInt :: i,j
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec
  Vec :: global_x_cell_vec,global_y_cell_vec,global_z_cell_vec
  Vec :: natural_x_cell_vec,natural_y_cell_vec,natural_z_cell_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  PetscInt, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element, ugdm_cell
  PetscErrorCode :: ierr

  PetscInt :: TET_ID_XDMF = 6
  PetscInt :: PYR_ID_XDMF = 7
  PetscInt :: WED_ID_XDMF = 8
  PetscInt :: HEX_ID_XDMF = 9

  grid => realization_base%patch%grid

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_x_vertex_vec,ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_y_vertex_vec,ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_z_vertex_vec,ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr)

  call GetVertexCoordinates(grid, global_x_vertex_vec,X_COORDINATE,option)
  call GetVertexCoordinates(grid, global_y_vertex_vec,Y_COORDINATE,option)
  call GetVertexCoordinates(grid, global_z_vertex_vec,Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr)

#if defined(SCORPIO_WRITE)
  write(*,*),'SCORPIO_WRITE'
  option%io_buffer = 'WriteHDF5CoordinatesUGrid not supported for SCORPIO_WRITE'
  call printErrMsg(option)
#else

  !
  !        not(SCORPIO_WRITE)
  !

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%num_vertices_global
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = 3

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(double_array(local_size*3))
  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

  deallocate(double_array)
  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(global_x_vertex_vec,vec_x_ptr,ierr)
  call VecRestoreArrayF90(global_y_vertex_vec,vec_y_ptr,ierr)
  call VecRestoreArrayF90(global_z_vertex_vec,vec_z_ptr,ierr)


  call VecDestroy(global_x_vertex_vec,ierr)
  call VecDestroy(global_y_vertex_vec,ierr)
  call VecDestroy(global_z_vertex_vec,ierr)

  !
  !  Write elements
  !
  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call GetCellConnections(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)

  local_size = grid%unstructured_grid%nlmax

  vert_count=0
  do i=1,local_size*EIGHT_INTEGER
    if(int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%nlmax

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = vert_count
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  call MPI_Allreduce(vert_count,dims(1),ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  realization_base%output_option%xmf_vert_len=int(dims(1))

  ! file space which is a 2D block
  rank_mpi = 1
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(vert_count, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(1) = istart
  length(1) = vert_count
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(int_array(vert_count))

  vert_count=0
  do i=1,local_size
    nverts=0
    do j=1,8
      if(vec_ptr((i-1)*8+j)>0) nverts=nverts+1
    enddo
    vert_count=vert_count+1
    select case (nverts)
      case (4) ! Tetrahedron
        int_array(vert_count) = TET_ID_XDMF
      case (5) ! Pyramid
        int_array(vert_count) = PYR_ID_XDMF
      case (6) ! Wedge
        int_array(vert_count) = WED_ID_XDMF
      case (8) ! Hexahedron
        int_array(vert_count) = HEX_ID_XDMF
    end select

    do j=1,8
      if(vec_ptr((i-1)*8+j)>0) then
        vert_count=vert_count+1
        int_array(vert_count) = INT(vec_ptr((i-1)*8+j))-1
      endif
    enddo
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

  deallocate(int_array)
  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)
  call VecDestroy(global_vec,ierr)
  call VecDestroy(natural_vec,ierr)
  call UGridDMDestroy(ugdm_element)

  ! Cell center X/Y/Z
  call VecCreateMPI(option%mycomm,grid%nlmax, &
                    PETSC_DETERMINE, &
                    global_x_cell_vec,ierr)
  call VecCreateMPI(option%mycomm,grid%nlmax, &
                    PETSC_DETERMINE, &
                    global_y_cell_vec,ierr)
  call VecCreateMPI(option%mycomm,grid%nlmax, &
                    PETSC_DETERMINE, &
                    global_z_cell_vec,ierr)

  call GetCellCoordinates(grid, global_x_cell_vec,X_COORDINATE)
  call GetCellCoordinates(grid, global_y_cell_vec,Y_COORDINATE)
  call GetCellCoordinates(grid, global_z_cell_vec,Z_COORDINATE)


  call UGridCreateUGDM(grid%unstructured_grid,ugdm_cell,ONE_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell,natural_x_cell_vec, &
                           NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell,natural_y_cell_vec, &
                           NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell,natural_z_cell_vec, &
                           NATURAL,option)
                           
  call VecScatterBegin(ugdm_cell%scatter_gton,global_x_cell_vec,natural_x_cell_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_x_cell_vec,natural_x_cell_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecScatterBegin(ugdm_cell%scatter_gton,global_y_cell_vec,natural_y_cell_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_y_cell_vec,natural_y_cell_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecScatterBegin(ugdm_cell%scatter_gton,global_z_cell_vec,natural_z_cell_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_z_cell_vec,natural_z_cell_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(natural_x_cell_vec,vec_x_ptr,ierr)
  call VecGetArrayF90(natural_y_cell_vec,vec_y_ptr,ierr)
  call VecGetArrayF90(natural_z_cell_vec,vec_z_ptr,ierr)
  local_size = grid%unstructured_grid%nlmax

  ! XC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "XC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_x_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  ! YC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "YC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_y_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  ! ZC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "ZC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_z_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)


  call VecRestoreArrayF90(natural_x_cell_vec,vec_x_ptr,ierr)
  call VecRestoreArrayF90(natural_y_cell_vec,vec_y_ptr,ierr)
  call VecRestoreArrayF90(natural_z_cell_vec,vec_z_ptr,ierr)

  call VecDestroy(global_x_cell_vec,ierr)
  call VecDestroy(global_y_cell_vec,ierr)
  call VecDestroy(global_z_cell_vec,ierr)

  call VecDestroy(natural_x_cell_vec,ierr)
  call VecDestroy(natural_y_cell_vec,ierr)
  call VecDestroy(natural_z_cell_vec,ierr)

#endif
!if defined(SCORPIO_WRITE)

end subroutine WriteHDF5CoordinatesUGridXDMF
#endif

end module Output_HDF5_module
