module Output_module

  use Logging_module 
  use Output_Aux_module

 ! use Output_Surface_module
  use Output_HDF5_module
  use Output_Tecplot_module
  use Output_VTK_module
  use Output_Observation_module
  
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petsclog.h"

#if defined(SCORPIO_WRITE)
  include "scorpiof.h"
#endif

  PetscInt, parameter :: TECPLOT_INTEGER = 0
  PetscInt, parameter :: TECPLOT_REAL = 1

  PetscInt, parameter :: VTK_INTEGER = 0
  PetscInt, parameter :: VTK_REAL = 1

  PetscInt, parameter :: TECPLOT_FILE = 0
  PetscInt, parameter ::  HDF5_FILE = 1

  
  PetscBool :: observation_first
  PetscBool :: hdf5_first
  PetscBool :: mass_balance_first

  public :: OutputInit, &
            Output, &
            OutputPrintCouplers

contains

! ************************************************************************** !
!
! OutputInit: Initializes variables
! author: Glenn Hammond
! date: 01/22/09
!
! ************************************************************************** !
subroutine OutputInit(realization_base,num_steps)

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Output_Common_module

  implicit none
  
  class(realization_base_type) :: realization_base
  PetscInt :: num_steps
  
  call OutputCommonInit(realization_base,num_steps)
  call OutputObservationInit(realization_base,num_steps)
  call OutputHDF5Init(realization_base,num_steps)

end subroutine OutputInit

! ************************************************************************** !
!
! Output: Main driver for all output subroutines
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine Output(realization_base,plot_flag,transient_plot_flag)

  use Realization_Base_class, only : realization_base_type
  use Option_module, only : OptionCheckTouch, option_type, printMsg
  use Grid_module, only : UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC
  
  implicit none
  
  class(realization_base_type) :: realization_base
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend
  type(option_type), pointer :: option

  option => realization_base%option
  
  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr)

  ! check for plot request from active directory
  if (.not.plot_flag) then

    if (option%use_touch_options) then
      string = 'plot'
      if (OptionCheckTouch(option,string)) then
        realization_base%output_option%plot_name = 'plot'
        plot_flag = PETSC_TRUE
      endif
    endif

  endif

  if (plot_flag) then
  
    if (realization_base%output_option%print_hdf5) then
      call PetscTime(tstart,ierr)
      call PetscLogEventBegin(logging%event_output_hdf5,ierr)    
      if (realization_base%discretization%itype == UNSTRUCTURED_GRID .or. &
          realization_base%discretization%itype == UNSTRUCTURED_GRID_MIMETIC) then
        call OutputHDF5UGridXDMF(realization_base,INSTANTANEOUS_VARS)
      else
        call OutputHDF5(realization_base,INSTANTANEOUS_VARS)
      endif      
      call PetscLogEventEnd(logging%event_output_hdf5,ierr)    
      call PetscTime(tend,ierr)
#ifdef SCORPIO_WRITE
      if (option%myrank == 0) write (*,'(" Parallel IO Write method is used in & 
                               &writing the output, HDF5_WRITE_GROUP_SIZE = ",i5)') &
                               option%hdf5_write_group_size
#endif
      write(option%io_buffer,'(f10.2," Seconds to write HDF5 file.")') tend-tstart
      call printMsg(option)
    endif
   
    if (realization_base%output_option%print_tecplot) then
      call PetscTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_tecplot,ierr) 
      select case(realization_base%output_option%tecplot_format)
        case (TECPLOT_POINT_FORMAT)
          call OutputTecplotPoint(realization_base)
        case (TECPLOT_BLOCK_FORMAT,TECPLOT_FEBRICK_FORMAT)
          call OutputTecplotBlock(realization_base)
      end select
      call PetscLogEventEnd(logging%event_output_tecplot,ierr)    
      call PetscTime(tend,ierr) 
      write(option%io_buffer,'(f10.2," Seconds to write to Tecplot file(s)")') &
            tend-tstart
      call printMsg(option)        
    endif
    
    if (realization_base%output_option%print_explicit_flowrate) then
      call PetscTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_tecplot,ierr) 
      call OutputPrintExplicitFlowrates(realization_base)
      call PetscLogEventEnd(logging%event_output_tecplot,ierr)    
      call PetscTime(tend,ierr) 
      write(option%io_buffer,'(f10.2," Seconds to write to Rates file.")') &
            tend-tstart
      call printMsg(option)        
    endif

    if (realization_base%output_option%print_vtk) then
      call PetscTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_vtk,ierr) 
      call OutputVTK(realization_base)

      call PetscLogEventEnd(logging%event_output_vtk,ierr)    
      call PetscTime(tend,ierr) 
      write(option%io_buffer,'(f10.2," Seconds to write to VTK file(s)")') &
            tend-tstart
      call printMsg(option) 
    endif
      
    if (realization_base%output_option%print_mad) then
      call PetscTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_mad,ierr) 
      call OutputMAD(realization_base)

      call PetscLogEventEnd(logging%event_output_mad,ierr)    
      call PetscTime(tend,ierr) 
      write(option%io_buffer,'(f10.2," Seconds to write to MAD HDF5 file(s)")') &
            tend-tstart
      call printMsg(option) 
    endif
    
    ! Print secondary continuum variables vs sec. continuum dist.
    if (option%use_mc) then
      if (realization_base%output_option%print_tecplot) then
        call PetscTime(tstart,ierr) 
        call PetscLogEventBegin(logging%event_output_secondary_tecplot,ierr) 
        call OutputSecondaryContinuumTecplot(realization_base)
        call PetscLogEventEnd(logging%event_output_secondary_tecplot,ierr)    
        call PetscTime(tend,ierr) 
        write(option%io_buffer,'(f10.2," Seconds to write to secondary' // &
              ' continuum Tecplot file(s)")') &
              tend-tstart
        call printMsg(option) 
      endif
    endif
      
    if (option%compute_statistics) then
      call ComputeFlowCellVelocityStats(realization_base)
      call ComputeFlowFluxVelocityStats(realization_base)
    endif

  endif
  
  if (transient_plot_flag) then
    if (option%compute_mass_balance_new) then
      call OutputMassBalance(realization_base)
    endif
    call OutputObservation(realization_base)
  endif

  ! Output temporally average variables 
  call OutputAvegVars(realization_base)

  if(plot_flag) then
    realization_base%output_option%plot_number = realization_base%output_option%plot_number + 1
  endif

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  realization_base%output_option%plot_name = ''

  call PetscLogStagePop(ierr)


end subroutine Output

! ************************************************************************** !
!
! OutputMAD: Print to HDF5 file for MAD final output
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputMAD(realization_base)

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use Variables_module
  use Output_Common_module, only : OutputGetVarFromArray
 
#if !defined(PETSC_HAVE_HDF5)
  implicit none
  
  class(realization_base_type) :: realization_base

  call printMsg(realization_base%option,'')
  write(realization_base%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to ", &
        &"write HDF5 formatted structured grids.")')
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
  
  implicit none

  class(realization_base_type) :: realization_base

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer(HSIZE_T) :: dims(3)
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  
  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_flag 
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  reaction => realization_base%reaction
  output_option => realization_base%output_option

#define ALL
#ifdef ALL
  write(string,'(i6)') option%mygroup_id
  filename = trim(option%global_prefix) // '-MAD-G' // trim(adjustl(string)) // '.h5'
!  filename = trim(option%global_prefix) // '-MAD.h5'

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  ! turn off error reporting
  call h5eset_auto_f(OFF,hdf5_err)
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                     prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)
#else
  filename = trim(option%global_prefix) // trim(option%group_prefix) // '.h5'

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                   prop_id)
  call h5pclose_f(prop_id,hdf5_err)
#endif

  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)   

  ! pressure
  call OutputGetVarFromArray(realization_base,global_vec,LIQUID_PRESSURE,ZERO_INTEGER)
#ifdef ALL
  string = 'Pressure' // trim(option%group_prefix)
#else
  string = 'Pressure'
#endif
  call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec,file_id,H5T_NATIVE_DOUBLE)

  call VecDestroy(global_vec,ierr)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif
end subroutine OutputMAD

! ************************************************************************** !
!
! ComputeFlowCellVelocityStats: 
! author: Glenn Hammond
! date: 03/11/08
!
! ************************************************************************** !
subroutine ComputeFlowCellVelocityStats(realization_base)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Patch_module
  use Discretization_module

  implicit none
  
  class(realization_base_type) :: realization_base
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(discretization_type), pointer :: discretization
  type(output_option_type), pointer :: output_option
  PetscInt :: iconn, i, direction, iphase, sum_connection
  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscReal :: flux
  Vec :: global_vec, global_vec2

  PetscReal :: average, sum, max, min, std_dev
  PetscInt :: max_loc, min_loc
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:), den_loc_p(:)
  PetscReal, allocatable :: sum_area(:)
  PetscErrorCode :: ierr
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  discretization => realization_base%discretization
    
  allocate(sum_area(grid%nlmax))
  call DiscretizationDuplicateVector(discretization,field%porosity0,global_vec)
  call DiscretizationDuplicateVector(discretization,field%porosity0,global_vec2)

  do iphase = 1,option%nphase

    do direction = 1,3
    
      sum_area(1:grid%nlmax) = 0.d0
      call VecSet(global_vec,0.d0,ierr)
      call VecGetArrayF90(global_vec,vec_ptr,ierr)

      ! interior velocities  
      connection_set_list => grid%internal_connection_set_list
      cur_connection_set => connection_set_list%first
      sum_connection = 0
      do 
        if (.not.associated(cur_connection_set)) exit
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          ghosted_id_up = cur_connection_set%id_up(iconn)
          ghosted_id_dn = cur_connection_set%id_dn(iconn)
          local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
          local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
          ! velocities are stored as the downwind face of the upwind cell
          flux = patch%internal_velocities(iphase,sum_connection)* &
                   cur_connection_set%area(iconn)* &
                   cur_connection_set%dist(direction,iconn)
          if (local_id_up > 0) then
            vec_ptr(local_id_up) = vec_ptr(local_id_up) - flux
          endif
          if (local_id_dn > 0) then
            vec_ptr(local_id_dn) = vec_ptr(local_id_dn) + flux
          endif
        enddo
        cur_connection_set => cur_connection_set%next
      enddo

      ! boundary velocities
      boundary_condition => patch%boundary_conditions%first
      sum_connection = 0
      do
        if (.not.associated(boundary_condition)) exit
        cur_connection_set => boundary_condition%connection_set
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          local_id = cur_connection_set%id_dn(iconn)
          vec_ptr(local_id) = vec_ptr(local_id)+ &
                              cur_connection_set%dist(direction,iconn)* &
                              patch%boundary_velocities(iphase,sum_connection)* &
                              cur_connection_set%area(iconn)
        enddo
        boundary_condition => boundary_condition%next
      enddo

      call VecRestoreArrayF90(global_vec,vec_ptr,ierr)

      call VecSum(global_vec,sum,ierr)
      average = sum/real(grid%nmax)
      call VecSet(global_vec2,average,ierr)
      call VecMax(global_vec,max_loc,max,ierr)
      call VecMin(global_vec,min_loc,min,ierr)
      call VecAYPX(global_vec2,-1.d0,global_vec,ierr)
      call VecNorm(global_vec2,NORM_2,std_dev,ierr)
      select case(direction)
        case(X_DIRECTION)
          string = 'X-Direction,'
        case(Y_DIRECTION)
          string = 'Y-Direction,'
        case(Z_DIRECTION)
          string = 'Z-Direction,'
      end select
      select case(iphase)
        case(LIQUID_PHASE)
          string = trim(string) // ' Liquid Phase'
        case(GAS_PHASE)
          string = trim(string) // ' Gas Phase'
      end select
      string = trim(string) // ' Velocity Statistics [m/' // &
               trim(output_option%tunit) // ']:'

      if (option%myrank == option%io_rank) then
        write(*,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
        write(option%fid_out,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
      endif

    enddo
  enddo
  
  if (allocated(sum_area)) deallocate(sum_area)
  call VecDestroy(global_vec,ierr)
  call VecDestroy(global_vec2,ierr)

end subroutine ComputeFlowCellVelocityStats

! ************************************************************************** !
!
! ComputeFlowFluxVelocityStats: Print flux statistics
! author: Glenn Hammond
! date: 03/11/08
!
! ************************************************************************** !
subroutine ComputeFlowFluxVelocityStats(realization_base)
!geh - specifically, the flow velocities at the interfaces between cells
 
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  
  implicit none

  class(realization_base_type) :: realization_base
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: discretization  
  type(output_option_type), pointer :: output_option
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscInt :: iphase
  PetscInt :: direction
  PetscInt :: local_id, ghosted_id
  PetscInt :: iconn, sum_connection
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, global_vec2
  PetscReal :: sum, average, max, min , std_dev
  PetscInt :: max_loc, min_loc
  PetscErrorCode :: ierr

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
    
  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  
  call DiscretizationDuplicateVector(discretization,field%porosity0,global_vec) 
  call DiscretizationDuplicateVector(discretization,field%porosity0,global_vec2) 

  do iphase = 1,option%nphase
    do direction = 1,3
    
      call VecZeroEntries(global_vec,ierr)
      call VecGetArrayF90(global_vec,vec_ptr,ierr)
      
      ! place interior velocities in a vector
      connection_set_list => grid%internal_connection_set_list
      cur_connection_set => connection_set_list%first
      sum_connection = 0
      do 
        if (.not.associated(cur_connection_set)) exit
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          ghosted_id = cur_connection_set%id_up(iconn)
          local_id = grid%nG2L(ghosted_id) ! = zero for ghost nodes
          ! velocities are stored as the downwind face of the upwind cell
          if (local_id <= 0 .or. &
              dabs(cur_connection_set%dist(direction,iconn)) < 0.99d0) cycle
          vec_ptr(local_id) = patch%internal_velocities(iphase,sum_connection)
        enddo
        cur_connection_set => cur_connection_set%next
      enddo

      call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
      
      ! compute stats
      call VecSum(global_vec,sum,ierr)
      average = sum/real(grid%nmax)
      call VecSet(global_vec2,average,ierr)
      call VecMax(global_vec,max_loc,max,ierr)
      call VecMin(global_vec,min_loc,min,ierr)
      call VecAYPX(global_vec2,-1.d0,global_vec,ierr)
      call VecNorm(global_vec2,NORM_2,std_dev,ierr)
      select case(direction)
        case(X_DIRECTION)
          string = 'X-Direction,'
        case(Y_DIRECTION)
          string = 'Y-Direction,'
        case(Z_DIRECTION)
          string = 'Z-Direction,'
      end select
      select case(iphase)
        case(LIQUID_PHASE)
          string = trim(string) // ' Liquid Phase'
        case(GAS_PHASE)
          string = trim(string) // ' Gas Phase'
      end select
      string = trim(string) // ' Flux Velocity Statistics [m/' // &
               trim(output_option%tunit) // ']:'
      if (option%myrank == option%io_rank) then
        write(*,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
        write(option%fid_out,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
      endif
    enddo
  enddo
  
  call VecDestroy(global_vec,ierr)
  call VecDestroy(global_vec2,ierr)
  
end subroutine ComputeFlowFluxVelocityStats

! ************************************************************************** !
!
! OutputPrintCouplers: Prints values of auxiliary variables associated with
!                      couplers (boundary and initial conditions, source
!                      sinks).  Note that since multiple connections for
!                      couplers can exist for a single cell, the latter will
!                      overwrite the former.
! author: Glenn Hammond
! date: 11/02/11
!
! ************************************************************************** !
subroutine OutputPrintCouplers(realization_base,istep)

  use Realization_Base_class, only : realization_base_type
  use Coupler_module
  use Connection_module
  use Option_module
  use Debug_module
  use Field_module
  use Patch_module
  use Level_module
  use Grid_module
  use Input_module

  class(realization_base_type) :: realization_base
  PetscInt :: istep
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: cur_patch
  type(level_type), pointer :: cur_level
  type(field_type), pointer :: field
  type(coupler_type), pointer :: coupler
  type(flow_debug_type), pointer :: flow_debug
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string, coupler_string
  type(connection_set_type), pointer :: cur_connection_set
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id, iconn, iauxvar
  PetscErrorCode :: ierr
  
  
  option => realization_base%option
  flow_debug => realization_base%debug
  field => realization_base%field

  if (len_trim(flow_debug%coupler_string) == 0) then
    option%io_buffer = &
      'Coupler debugging requested, but no string of coupler names was included.'
    call printErrMsg(option)
  endif
  
  coupler_string = flow_debug%coupler_string
  ierr = 0
  do
    call InputReadWord(coupler_string,word,PETSC_TRUE,ierr)
    if (ierr /= 0) exit
    
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        iauxvar = RICHARDS_PRESSURE_DOF
      case default
        option%io_buffer = &
          'OutputPrintCouplers() not yet supported for this flow mode'
        call printErrMsg(option)
    end select
    
    cur_level => realization_base%level_list%first
    do 
      if (.not.associated(cur_level)) exit
      cur_patch => cur_level%patch_list%first
      do
        if (.not.associated(cur_patch)) exit
        grid => cur_patch%grid
        coupler => CouplerGetPtrFromList(word,cur_patch%boundary_conditions)
        call VecZeroEntries(field%work,ierr)
        call GridVecGetArrayF90(grid,field%work,vec_ptr,ierr)
        if (associated(coupler)) then
          cur_connection_set => coupler%connection_set
          do iconn = 1, cur_connection_set%num_connections
            local_id = cur_connection_set%id_dn(iconn)
            if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
            vec_ptr(local_id) = coupler%flow_aux_real_var(iauxvar,iconn)
          enddo
        endif
        call GridVecRestoreArrayF90(grid,field%work,vec_ptr,ierr)
        cur_patch => cur_patch%next
      enddo
      cur_level => cur_level%next
    enddo

    if (istep > 0) then
      write(string,*) istep
      string = adjustl(string)
    else 
      string = ''
    endif
    string = trim(word) // trim(string)
    if (len_trim(option%group_prefix) > 1) then
      string = trim(string) // trim(option%group_prefix)
    endif
    string = trim(string) // '.tec'
    call OutputVectorTecplot(string,word,realization_base,field%work)
      
  enddo

end subroutine OutputPrintCouplers

! ************************************************************************** !
!> This routine temporally averages variables and outputs thems
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/10/13
! ************************************************************************** !
subroutine OutputAvegVars(realization_base)

  use Realization_Base_class, only : realization_base_type
  use Option_module, only : OptionCheckTouch, option_type, printMsg
  use Output_Aux_module
  use Output_Common_module, only : OutputGetVarFromArray  
  use Field_module
  use Grid_module, only : UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC

  implicit none
  
  class(realization_base_type) :: realization_base

  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  type(field_type), pointer :: field  

  PetscReal :: dtime
  PetscBool :: aveg_plot_flag
  PetscInt :: ivar
  PetscReal,pointer :: aval_p(:),ival_p(:)
  PetscErrorCode :: ierr  
  PetscLogDouble :: tstart, tend

  option => realization_base%option
  output_option => realization_base%output_option
  field => realization_base%field

  ! 
  if(option%time<1.d-10) return
  
  dtime = option%time-output_option%aveg_var_time
  output_option%aveg_var_dtime = output_option%aveg_var_dtime + dtime
  output_option%aveg_var_time = output_option%aveg_var_time + dtime
  
  if(abs(output_option%aveg_var_dtime-output_option%periodic_output_time_incr)<1.d0) then
    aveg_plot_flag=PETSC_TRUE
  else
    aveg_plot_flag=PETSC_FALSE
  endif

  if(.not.associated(output_option%aveg_output_variable_list%first)) then
    if(output_option%print_hdf5_aveg_mass_flowrate.or. &
       output_option%print_hdf5_aveg_energy_flowrate) then
      ! There is a possibility to output average-flowrates, thus
      ! call output subroutine depending on mesh type
      if (realization_base%discretization%itype == UNSTRUCTURED_GRID.or. &
          realization_base%discretization%itype == UNSTRUCTURED_GRID_MIMETIC) then
        call OutputHDF5UGridXDMF(realization_base,AVERAGED_VARS)
      else
      !  call OutputHDF5(realization_base,AVERAGED_VARS)
      endif
    endif
    return
  endif
  
  ivar = 0
  cur_variable => output_option%aveg_output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit

    ! Get the variable
    call OutputGetVarFromArray(realization_base,field%work, &
                               cur_variable%ivar, &
                               cur_variable%isubvar)

    ! Cumulatively add the variable*dtime
    ivar = ivar + 1
    call VecGetArrayF90(field%work,ival_p,ierr)
    call VecGetArrayF90(field%avg_vars_vec(ivar),aval_p,ierr)
    aval_p = aval_p + ival_p*dtime
    call VecRestoreArrayF90(field%work,ival_p,ierr)
    call VecRestoreArrayF90(field%avg_vars_vec(ivar),aval_p,ierr)

    ! Check if it is time to output the temporally average variable
    if(aveg_plot_flag) then

      ! Divide vector values by 'time'
      call VecGetArrayF90(field%avg_vars_vec(ivar),aval_p,ierr)
      aval_p = aval_p/output_option%periodic_output_time_incr
      call VecRestoreArrayF90(field%avg_vars_vec(ivar),aval_p,ierr)

    endif
    
    cur_variable => cur_variable%next
  enddo

  if(aveg_plot_flag) then

    if (realization_base%output_option%print_hdf5) then
      call PetscTime(tstart,ierr)
      call PetscLogEventBegin(logging%event_output_hdf5,ierr)    
      if (realization_base%discretization%itype == UNSTRUCTURED_GRID.or. &
          realization_base%discretization%itype == UNSTRUCTURED_GRID_MIMETIC) then
        call OutputHDF5UGridXDMF(realization_base,AVERAGED_VARS)
      else
        call OutputHDF5(realization_base,AVERAGED_VARS)
      endif      
      call PetscLogEventEnd(logging%event_output_hdf5,ierr)    
      call PetscTime(tend,ierr)
      write(option%io_buffer,'(f10.2," Seconds to write HDF5 file.")') tend-tstart
      call printMsg(option)
    endif

    ! Reset the vectors to zero
    do ivar=1,output_option%aveg_output_variable_list%nvars
      call VecSet(field%avg_vars_vec(ivar),0.d0,ierr)
    enddo

    output_option%aveg_var_dtime=0.d0

  endif


end subroutine OutputAvegVars

end module Output_module
