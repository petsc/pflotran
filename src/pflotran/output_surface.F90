#ifdef SURFACE_FLOW
module Output_Surface_module

  use Logging_module 
  use Output_Aux_module
  use Output_Common_module
  use Output_HDF5_module
  use Output_Tecplot_module
  
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

  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: hydrograph_first
  PetscBool :: surf_hdf5_first
  
  public :: OutputSurfaceInit, &
            OutputSurface

contains

! ************************************************************************** !
!
! OutputSurfaceInit: Initializes module variables for surface variables
! author: Glenn Hammond
! date: 01/16/13
!
! ************************************************************************** !
subroutine OutputSurfaceInit(realization_base,num_steps)

  use Realization_Base_class, only : realization_base_type
  use Option_module

  implicit none
  
  class(realization_base_type) :: realization_base
  PetscInt :: num_steps

  if (num_steps == 0) then
    hydrograph_first = PETSC_TRUE
    surf_hdf5_first = PETSC_TRUE
  else
    hydrograph_first = PETSC_FALSE
    surf_hdf5_first = PETSC_FALSE
  endif

end subroutine OutputSurfaceInit

! ************************************************************************** !
!> This subroutine is main driver for all output subroutines related to
!! surface flows.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/29/12
! ************************************************************************** !
subroutine OutputSurface(surf_realization,realization,plot_flag, &
                         transient_plot_flag)

  use Surface_Realization_class, only : surface_realization_type
  use Realization_class, only : realization_type
  use Option_module, only : OptionCheckTouch, option_type, &
                            printMsg, printErrMsg

  implicit none

  class(surface_realization_type) :: surf_realization
  class(realization_type) :: realization
  PetscBool                      :: plot_flag
  PetscBool                      :: transient_plot_flag

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode                 :: ierr
  PetscLogDouble                 :: tstart, tend
  type(option_type), pointer     :: option

  option => surf_realization%option

  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr)

  ! check for plot request from active directory
  if (.not.plot_flag) then

    if (option%use_touch_options) then
      string = 'plot'
      if (OptionCheckTouch(option,string)) then
        surf_realization%output_option%plot_name = 'plot'
        plot_flag = PETSC_TRUE
      endif
    endif
  endif

  if (plot_flag) then
    if (surf_realization%output_option%print_hdf5) then
      call OutputSurfaceHDF5UGridXDMF(surf_realization,realization)
    endif
  
    if (surf_realization%output_option%print_tecplot) then
      call PetscGetTime(tstart,ierr)
      call PetscLogEventBegin(logging%event_output_tecplot,ierr) 
      select case(surf_realization%output_option%tecplot_format)
        case (TECPLOT_FEQUADRILATERAL_FORMAT)
          call OutputTecplotFEQUAD(surf_realization,realization)
      end select
      call PetscGetTime(tend,ierr)
      call PetscLogEventEnd(logging%event_output_tecplot,ierr)
    endif

    if (surf_realization%output_option%print_hydrograph) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_hydrograph,ierr)
      call OutputHydrograph(surf_realization)
      call PetscGetTime(tend,ierr)
      call PetscLogEventEnd(logging%event_output_hydrograph,ierr)
    endif
    surf_realization%output_option%plot_number = surf_realization%output_option%plot_number + 1
  endif

  call PetscLogStagePop(ierr)

end subroutine OutputSurface

! ************************************************************************** !
!> This subroutine print to Tecplot file in FEQUADRILATERAL format for surface
!! flows.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/29/12
! ************************************************************************** !
subroutine OutputTecplotFEQUAD(surf_realization,realization)

  use Surface_Realization_class, only : surface_realization_type
  use Realization_class, only : realization_type
  use Discretization_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Surface_Field_module
  use Patch_module
  
  implicit none

  type(surface_realization_type) :: surf_realization
  class(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  PetscInt, parameter :: icolumn = -1
  character(len=MAXSTRINGLENGTH) :: filename, string, string2
  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXSTRINGLENGTH) :: tmp_global_prefix
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(patch_type), pointer :: patch 
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  Vec :: global_cconn_vec
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr  
  
  type(ugdm_type), pointer :: ugdm_element
  
  discretization => surf_realization%discretization
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  output_option => surf_realization%output_option

  tmp_global_prefix = option%global_prefix 
  option%global_prefix = trim(tmp_global_prefix) // '-surf'
  filename = OutputFilename(output_option,option,'tec','')
  option%global_prefix = tmp_global_prefix
    
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=filename,action="write")
    call OutputTecplotHeader(OUTPUT_UNIT,surf_realization,icolumn)
  endif

  ! write blocks
  ! write out data sets
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)

  ! write out coordinates
  !TODO(gautam): check this.  It was flipped (realization was uncommented) before.
  call WriteTecplotUGridVertices(OUTPUT_UNIT,surf_realization)
  !call WriteTecplotUGridVertices(OUTPUT_UNIT,realization)

  ! loop over variables and write to file
  cur_variable => output_option%output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputGetVarFromArray(realization,global_vec,cur_variable%ivar, &
                                cur_variable%isubvar)
    call DiscretizationGlobalToNatural(discretization,global_vec, &
                                        natural_vec,ONEDOF)
    if (cur_variable%iformat == 0) then
      call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec, &
                                      TECPLOT_REAL)
    else
      call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec, &
                                      TECPLOT_INTEGER)
    endif
    cur_variable => cur_variable%next
  enddo

  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)
  
  ! write vertices
  call WriteTecplotUGridElements(OUTPUT_UNIT,surf_realization)

  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)

end subroutine OutputTecplotFEQUAD

! ************************************************************************** !
!> This subroutine prints header to Tecplot file for surface grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine OutputTecplotHeader(fid,surf_realization,icolumn)

  use Surface_Realization_class
  use Grid_module
  use Option_module
  use Patch_module

  implicit none

  PetscInt :: fid
  type(surface_realization_type) :: surf_realization
  PetscInt :: icolumn
  
  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  PetscInt :: comma_count, quote_count, variable_count
  PetscInt :: i
  
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  output_option => surf_realization%output_option

  ! write header
  ! write title
  write(fid,'(''TITLE = "'',1es13.5," [",a1,'']"'')') &
                option%surf_flow_time/output_option%tconv,output_option%tunit

  ! initial portion of header
  header = 'VARIABLES=' // &
            '"X [m]",' // &
            '"Y [m]",' // &
            '"Z [m]"'
  header2=''
  header2 = OutputVariableListToHeader(output_option%output_variable_list,'', &
                                      icolumn,PETSC_TRUE)

  header = trim(header) // trim(header2)
  write(fid,'(a)') trim(header)

  ! count vars in header
  quote_count = 0
  comma_count = 0
  do i=1,len_trim(header)
    ! 34 = '"'
    if (iachar(header(i:i)) == 34) then
      quote_count = quote_count + 1
    ! 44 = ','
    else if (iachar(header(i:i)) == 44 .and. mod(quote_count,2) == 0) then
      comma_count = comma_count + 1
    endif
  enddo
  
  variable_count = comma_count + 1

  !geh: due to pgi bug, cannot embed functions with calls to write() within
  !     write statement
  string = OutputTecplotZoneHeader(surf_realization,variable_count, &
                                   output_option%tecplot_format)
  write(fid,'(a)') trim(string)

end subroutine OutputTecplotHeader

! ************************************************************************** !
!> This subroutine prints zone header to Tecplot file.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
function OutputTecplotZoneHeader(surf_realization,variable_count,tecplot_format)

  use Surface_Realization_class
  use Grid_module
  use Option_module
  
  implicit none

  type(surface_realization_type) :: surf_realization
  PetscInt :: variable_count
  PetscInt :: tecplot_format
  
  character(len=MAXSTRINGLENGTH) :: OutputTecplotZoneHeader

  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  grid => surf_realization%patch%grid
  option => surf_realization%option
  output_option => surf_realization%output_option
  
  
  string = 'ZONE T="' // &
           trim(OutputFormatDouble(option%time/output_option%tconv)) // &
           '"'
  string2 = ''
  select case(tecplot_format)
    case (TECPLOT_POINT_FORMAT)
      if ((surf_realization%discretization%itype == STRUCTURED_GRID).or. &
          (surf_realization%discretization%itype == STRUCTURED_GRID_MIMETIC)) then
        string2 = ', I=' // &
                  trim(OutputFormatInt(grid%structured_grid%nx)) // &
                  ', J=' // &
                  trim(OutputFormatInt(grid%structured_grid%ny)) // &
                  ', K=' // &
                  trim(OutputFormatInt(grid%structured_grid%nz))
      else
        string2 = 'POINT format currently not supported for unstructured'
      endif
      string2 = trim(string2) // &
              ', DATAPACKING=POINT'
    case default !(TECPLOT_BLOCK_FORMAT,TECPLOT_FEBRICK_FORMAT)
      if ((surf_realization%discretization%itype == STRUCTURED_GRID).or. &
          (surf_realization%discretization%itype == STRUCTURED_GRID_MIMETIC)) then
        string2 = ', I=' // &
                  trim(OutputFormatInt(grid%structured_grid%nx+1)) // &
                  ', J=' // &
                  trim(OutputFormatInt(grid%structured_grid%ny+1)) // &
                  ', K=' // &
                  trim(OutputFormatInt(grid%structured_grid%nz+1))
      else
        string2 = ', N=' // &
                  trim(OutputFormatInt(grid%unstructured_grid%num_vertices_global)) // &
                  ', ELEMENTS=' // &
                  trim(OutputFormatInt(grid%unstructured_grid%nmax))
        string2 = trim(string2) // ', ZONETYPE=FEQUADRILATERAL'
      endif
  
      if (variable_count > 4) then
        string3 = ', VARLOCATION=([4-' // &
                  trim(OutputFormatInt(variable_count)) // &
                  ']=CELLCENTERED)'
      else
        string3 = ', VARLOCATION=([4]=CELLCENTERED)'
      endif
      string2 = trim(string2) // trim(string3) // ', DATAPACKING=BLOCK'
  end select
  
  OutputTecplotZoneHeader = trim(string) // string2

end function OutputTecplotZoneHeader

! ************************************************************************** !
!> This subroutine writes unstructured grid elements for surface grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine WriteTecplotUGridElements(fid, &
                                      surf_realization)

  use Surface_Realization_class
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module
  
  implicit none

  PetscInt :: fid
  type(surface_realization_type) :: surf_realization

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  Vec :: global_cconn_vec
  type(ugdm_type), pointer :: ugdm_element
  PetscReal, pointer :: vec_ptr(:),vec_ptr2(:)
  PetscInt :: ii
  PetscErrorCode :: ierr
  
  Vec :: global_vec
  Vec :: natural_vec
  Vec :: surface_natural_vec
  PetscInt :: natural_vec_local_size, surface_natural_vec_local_size

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  
  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call GetCellConnectionsTecplot(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  surface_natural_vec_local_size = 0
  call VecGetLocalSize(natural_vec,natural_vec_local_size,ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)
  do ii = 1,natural_vec_local_size
    if(vec_ptr(ii) /= -999) & 
      surface_natural_vec_local_size = surface_natural_vec_local_size + 1
  enddo
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)
  
  call VecCreateMPI(option%mycomm,surface_natural_vec_local_size,PETSC_DETERMINE,surface_natural_vec,ierr)
  call VecGetArrayF90(surface_natural_vec,vec_ptr2,ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)
  do ii = 1,surface_natural_vec_local_size
    vec_ptr2(ii) = vec_ptr(ii)
  enddo
  call VecRestoreArrayF90(surface_natural_vec,vec_ptr2,ierr)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)
  
  call VecGetArrayF90(surface_natural_vec,vec_ptr2,ierr)
  call WriteTecplotDataSetNumPerLine(fid,surf_realization,vec_ptr2, &
                                     TECPLOT_INTEGER, &
                                     surface_natural_vec_local_size, &
                                     FOUR_INTEGER)
  call VecRestoreArrayF90(surface_natural_vec,vec_ptr2,ierr)
  call VecDestroy(global_vec,ierr)
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(surface_natural_vec,ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine WriteTecplotUGridElements

! ************************************************************************** !
!> This subroutine writes unstructure grid vertices for surface grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine WriteTecplotUGridVertices(fid,surf_realization)

  use Surface_Realization_class
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module
  use Variables_module

  implicit none

  PetscInt :: fid
  type(surface_realization_type) :: surf_realization 
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch 
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  PetscInt :: local_size
  PetscErrorCode :: ierr
  
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_vertex_vec,ierr)
  call VecGetLocalSize(global_vertex_vec,local_size,ierr)
  call GetVertexCoordinates(grid, global_vertex_vec,X_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,surf_realization,vec_ptr,TECPLOT_REAL, &
                           local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

  call GetVertexCoordinates(grid,global_vertex_vec,Y_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,surf_realization,vec_ptr,TECPLOT_REAL, &
                           local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

  call GetVertexCoordinates(grid,global_vertex_vec, Z_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,surf_realization,vec_ptr,TECPLOT_REAL, &
                           local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

  call VecDestroy(global_vertex_vec, ierr)

end subroutine WriteTecplotUGridVertices

! ************************************************************************** !
!> This routine outputs hydrograph fluxes.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date:
! ************************************************************************** !

subroutine OutputHydrograph(surf_realization)

  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Connection_module
  use Utility_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(output_option_type), pointer :: output_option
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscInt :: icol
  PetscReal :: sum_flux, sum_flux_global

  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word, units
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  PetscInt :: fid = 86

  patch => surf_realization%patch
  option => surf_realization%option
  output_option => surf_realization%output_option
  
  if (output_option%print_column_ids) then
   icol = 1
  else
    icol = -1
  endif


  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    filename = trim(boundary_condition%name) // '_hydrograph.dat'
    
    if (option%myrank == option%io_rank) then
      if (hydrograph_first .or. .not.FileExists(filename)) then
        open(unit=fid,file=filename,action="write",status="replace")

        ! write header
        write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
      
        header = ''
        if (option%iflowmode > 0) then
          call OutputAppendToHeader(header,'dt_flow',output_option%tunit,'',icol)
        endif
        
        write(fid,'(a)',advance="no") trim(header)

        header = ''
        !string = 'Outflow '
        !call OutputAppendToHeader(header,string,'[m^2/s]','',icol)
        string = 'Outflow'
        call OutputAppendToHeader(header,string,'[m^3/s]','',icol)
        write(fid,'(a)',advance="no") trim(header)

        write(fid,'(a)') '' 
      else
        open(unit=fid,file=filename,action="write",status="old",position="append")
      endif
    endif

100 format(100es16.8)
110 format(100es16.8)

    ! write time
    if (option%myrank == option%io_rank) then
      write(fid,100,advance="no") option%surf_flow_time/output_option%tconv
    endif
  
    if (option%nflowdof > 0) then
      if (option%myrank == option%io_rank) &
        write(fid,100,advance="no") option%flow_dt/output_option%tconv
    endif

    sum_flux = 0.d0
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      !patch%boundary_velocities(1,sum_connection)
      sum_flux = sum_flux + patch%surf_boundary_fluxes(sum_connection)
    enddo
    
    call MPI_Reduce(sum_flux,sum_flux_global, &
                    ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)

    if (option%myrank == option%io_rank) then
      ! change sign for positive in / negative out
      write(fid,110,advance="no") -sum_flux_global
      write(fid,'(a)') ''
      close(fid)
    endif
    
    boundary_condition => boundary_condition%next
  enddo

  hydrograph_first = PETSC_FALSE

end subroutine OutputHydrograph

! ************************************************************************** !
!> This routine writes unstructured grid data in HDF5 XDMF format for 
!! surface flow.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/29/2012
! ************************************************************************** !
subroutine OutputSurfaceHDF5UGridXDMF(surf_realization,realization)

  use Surface_Realization_class
  use Realization_class
  use Discretization_module
  use Option_module
  use Grid_module
  use Surface_Field_module
  use Patch_module
  use Reaction_Aux_module

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  
  class(surface_realization_type) :: surf_realization
  class(realization_type) :: realization

  call printMsg(surf_realization%option,'')
  write(surf_realization%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to &
        &write HDF5 formatted structured grids Darn.")')
  call printErrMsg(surf_realization%option)
#endif

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

  class(surface_realization_type) :: surf_realization
  class(realization_type) :: realization

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

  type(grid_type), pointer :: subsurf_grid
  type(grid_type), pointer :: surf_grid
  type(discretization_type), pointer :: surf_discretization
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: xmf_filename, att_datasetname, group_name
  character(len=MAXSTRINGLENGTH) :: string
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

  surf_discretization => surf_realization%discretization
  patch => surf_realization%patch
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  output_option => surf_realization%output_option

  xmf_filename = OutputFilename(output_option,option,'xmf','surf')

  if (output_option%print_single_h5_file) then
    first = surf_hdf5_first
    filename = trim(option%global_prefix) // trim(option%group_prefix) // '-surf.h5'
  else
    string = OutputFilenameID(output_option,option)
    first = PETSC_TRUE
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                '-' // trim(string) // '-surf.h5'
  endif

  subsurf_grid => realization%patch%grid
  surf_grid    => surf_realization%patch%grid

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
    call WriteHDF5CoordinatesUGridXDMF(surf_realization,realization,option,grp_id)
    call h5gclose_f(grp_id,hdf5_err)
  endif

  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write xmf output file: ' // trim(xmf_filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=xmf_filename,action="write")
    call OutputXMFHeader(OUTPUT_UNIT, &
                         option%time/output_option%tconv, &
                         surf_grid%nmax, &
                         surf_realization%output_option%surf_xmf_vert_len, &
                         subsurf_grid%unstructured_grid%num_vertices_global,filename)

  endif

  ! create a group for the data set
  write(string,'(''Time'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif

  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  group_name=string
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  ! write out data sets 
  call DiscretizationCreateVector(surf_discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(surf_discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)

  ! loop over variables and write to file
  cur_variable => output_option%output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputGetVarFromArray(surf_realization,global_vec,cur_variable%ivar, &
                                cur_variable%isubvar)
    call DiscretizationGlobalToNatural(surf_discretization,global_vec, &
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
      call OutputXMFAttribute(OUTPUT_UNIT,surf_grid%nmax,string,att_datasetname)
    endif
    cur_variable => cur_variable%next
  enddo

  call VecDestroy(global_vec,ierr)
  call VecDestroy(natural_vec,ierr)
  call h5gclose_f(grp_id,hdf5_err)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

  if (option%myrank == option%io_rank) then
    call OutputXMFFooter(OUTPUT_UNIT)
    close(OUTPUT_UNIT)
  endif

  surf_hdf5_first = PETSC_FALSE

end subroutine OutputSurfaceHDF5UGridXDMF

! ************************************************************************** !
!> This routine writes unstructured coordinates to HDF5 file in XDMF format
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/29/2012
! ************************************************************************** !
subroutine WriteHDF5CoordinatesUGridXDMF(surf_realization,realization, &
                                          option,file_id)

  use hdf5
  use HDF5_module
  use Surface_Realization_class
  use Realization_class
  use Grid_module
  use Option_module
  use Unstructured_Grid_Aux_module
  use Variables_module
  
  implicit none
  
  class(surface_realization_type) :: surf_realization
  class(realization_type) :: realization
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

  PetscMPIInt :: hdf5_err
  type(grid_type), pointer :: surf_grid
  type(grid_type), pointer :: subsurf_grid
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: local_size,vert_count,nverts
  PetscInt :: i,j
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec
  PetscErrorCode :: ierr

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  PetscInt, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element

  PetscInt :: TRI_ID_XDMF = 4
  PetscInt :: QUA_ID_XDMF = 5

  surf_grid => surf_realization%patch%grid
  subsurf_grid => realization%patch%grid

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    subsurf_grid%unstructured_grid%num_vertices_global, &
                    global_x_vertex_vec,ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    subsurf_grid%unstructured_grid%num_vertices_global, &
                    global_y_vertex_vec,ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    subsurf_grid%unstructured_grid%num_vertices_global, &
                    global_z_vertex_vec,ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr)

  call GetVertexCoordinates(subsurf_grid, global_x_vertex_vec,X_COORDINATE,option)
  call GetVertexCoordinates(subsurf_grid, global_y_vertex_vec,Y_COORDINATE,option)
  call GetVertexCoordinates(subsurf_grid, global_z_vertex_vec,Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr)

#if defined(SCORPIO_WRITE)
  write(*,*),'SCORPIO_WRITE'
  option%io_buffer = 'WriteHDF5CoordinatesUGrid not supported for SCORPIO_WRITE'
  call printErrMsg(option)
#endif

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 2
  dims = 0
  dims(2) = subsurf_grid%unstructured_grid%num_vertices_global
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
  call UGridCreateUGDM(surf_grid%unstructured_grid,ugdm_element,FOUR_INTEGER,option)
  call UGridDMCreateVector(surf_grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(surf_grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call GetCellConnections(surf_grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)

  local_size = surf_grid%unstructured_grid%nlmax

  vert_count=0
  do i=1,local_size*FOUR_INTEGER
    if(int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+surf_grid%nlmax

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = vert_count
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  call MPI_Allreduce(vert_count,dims(1),ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  surf_realization%output_option%surf_xmf_vert_len=int(dims(1))

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
    do j=1,FOUR_INTEGER
      if(vec_ptr((i-1)*FOUR_INTEGER+j)>0) nverts=nverts+1
    enddo
    vert_count=vert_count+1
    select case (nverts)
      case (3) ! Triangle
        int_array(vert_count) = TRI_ID_XDMF
      case (4) ! Quadrilateral
        int_array(vert_count) = QUA_ID_XDMF
    end select

    do j=1,FOUR_INTEGER
      if(vec_ptr((i-1)*FOUR_INTEGER+j)>0) then
        vert_count=vert_count+1
        int_array(vert_count) = INT(vec_ptr((i-1)*FOUR_INTEGER+j))-1
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

end subroutine WriteHDF5CoordinatesUGridXDMF

end module Output_Surface_module
#endif
