#ifdef GEOMECH
module Output_Geomechanics_module

  use Geomechanics_Logging_module
  use Output_Aux_module  
  use Output_Tecplot_module
  use Output_Common_module
  
  implicit none
  
  private
  
#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petsclog.h"

  public :: OutputGeomechanics
      
contains

! ************************************************************************** !
!
! OutputGeomechanics: Main driver for all geomechanics output
! author: Satish Karra, LANL
! date: 07/2/13
!
! ************************************************************************** !
subroutine OutputGeomechanics(geomech_realization,plot_flag, &
                              transient_plot_flag)

  use Geomechanics_Realization_module
  use Option_module, only : OptionCheckTouch, option_type, &
                            printMsg, printErrMsg

  implicit none

  class(geomech_realization_type) :: geomech_realization
  PetscBool                       :: plot_flag
  PetscBool                       :: transient_plot_flag

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode                 :: ierr
  PetscLogDouble                 :: tstart, tend
  type(option_type), pointer     :: option

  option => geomech_realization%option

  call PetscLogStagePush(geomech_logging%stage(GEOMECH_OUTPUT_STAGE),ierr)

  ! check for plot request from active directory
  if (.not.plot_flag) then

    if (option%use_touch_options) then
      string = 'plot'
      if (OptionCheckTouch(option,string)) then
        geomech_realization%output_option%plot_name = 'plot'
        plot_flag = PETSC_TRUE
      endif
    endif
  endif

  if (plot_flag) then
   
    if (geomech_realization%output_option%print_tecplot) then
      call PetscTime(tstart,ierr)
      call PetscLogEventBegin(geomech_logging%event_output_tecplot,ierr) 
      call OutputTecplotGeomechanics(geomech_realization)
      call PetscTime(tend,ierr)
      call PetscLogEventEnd(geomech_logging%event_output_tecplot,ierr)
    endif

  endif

  ! Increment the plot number
  if(plot_flag) then
    geomech_realization%output_option%plot_number = &
      geomech_realization%output_option%plot_number + 1
  endif

  call PetscLogStagePop(ierr)

end subroutine OutputGeomechanics

! ************************************************************************** !
!
! OutputTecplotGeomechanics: Tecplot output for geomechanics
! author: Satish Karra, LANL
! date: 07/2/13
!
! ************************************************************************** !
subroutine OutputTecplotGeomechanics(geomech_realization)

  use Geomechanics_Realization_module
  use Geomechanics_Discretization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module
  use Option_module
  
  implicit none

  type(geomech_realization_type) :: geomech_realization
  
  PetscInt :: i, comma_count, quote_count
  PetscInt, parameter :: icolumn = -1
  character(len=MAXSTRINGLENGTH) :: filename, string, string2
  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXSTRINGLENGTH) :: tmp_global_prefix
  character(len=MAXWORDLENGTH) :: word
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_field_type), pointer :: geomech_field
  type(geomech_patch_type), pointer :: patch 
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  Vec :: global_cconn_vec
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr  
  
  type(gmdm_type), pointer :: gmdm_element
  
  discretization => geomech_realization%discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_field => geomech_realization%geomech_field
  output_option => geomech_realization%output_option

  tmp_global_prefix = option%global_prefix 
  option%global_prefix = trim(tmp_global_prefix) // '-geomech'
  filename = OutputFilename(output_option,option,'tec','')
  option%global_prefix = tmp_global_prefix
    
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=filename,action="write")
    call OutputTecplotHeader(OUTPUT_UNIT,geomech_realization,icolumn)
  endif

  ! write blocks
  ! write out data sets
  call GeomechDiscretizationCreateVector(discretization,ONEDOF, &
                                         global_vec,GLOBAL,option)
  call GeomechDiscretizationCreateVector(discretization,ONEDOF, &
                                         natural_vec,NATURAL,option)

  ! write out coordinates
  call WriteTecplotGeomechGridVertices(OUTPUT_UNIT,geomech_realization)

#if 0
  ! loop over variables and write to file
  cur_variable => output_option%output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputSurfaceGetVarFromArray(surf_realization,global_vec,cur_variable%ivar, &
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
#endif  

  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)

end subroutine OutputTecplotGeomechanics

! ************************************************************************** !
!
! OutputTecplotHeader: Prints Tecplot header for geomechanics
! author: Satish Karra, LANL
! date: 07/2/13
!
! ************************************************************************** !
subroutine OutputTecplotHeader(fid,geomech_realization,icolumn)

  use Geomechanics_Realization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Patch_module

  implicit none

  PetscInt :: fid
  type(geomech_realization_type) :: geomech_realization
  PetscInt :: icolumn
  
  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  PetscInt :: comma_count, quote_count, variable_count
  PetscInt :: i
  
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  output_option => geomech_realization%output_option

  ! write header
  ! write title
  write(fid,'(''TITLE = "'',1es13.5," [",a1,'']"'')') &
                option%geomech_time/output_option%tconv,output_option%tunit

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
  string = OutputTecplotZoneHeader(geomech_realization,variable_count, &
                                   output_option%tecplot_format)
  write(fid,'(a)') trim(string)

end subroutine OutputTecplotHeader

! ************************************************************************** !
!
! OutputTecplotZoneHeader: Prints zone header to a tecplot file
! author: Satish Karra, LANL
! date: 07/2/13
!
! ************************************************************************** !
function OutputTecplotZoneHeader(geomech_realization,variable_count, &
                                 tecplot_format)

  use Geomechanics_Realization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  
  implicit none

  type(geomech_realization_type) :: geomech_realization
  PetscInt :: variable_count
  PetscInt :: tecplot_format
  
  character(len=MAXSTRINGLENGTH) :: OutputTecplotZoneHeader

  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  grid => geomech_realization%geomech_patch%geomech_grid
  option => geomech_realization%option
  output_option => geomech_realization%output_option
  
  
  string = 'ZONE T="' // &
           trim(OutputFormatDouble(option%time/output_option%tconv)) // &
           '"'
  string2 = ''
  select case(tecplot_format)
    case (TECPLOT_POINT_FORMAT)
      string2 = 'POINT format not supported for geomechanics ' // &
                'unstructured'
      string2 = trim(string2) // &
              ', DATAPACKING=POINT'
    case default !(TECPLOT_BLOCK_FORMAT,TECPLOT_FEBRICK_FORMAT)
      string2 = ', N=' // &
                trim(OutputFormatInt(grid%nmax_node)) // &
                ', ELEMENTS=' // &
                trim(OutputFormatInt(grid%nmax_elem))
      string2 = trim(string2) // ', ZONETYPE=FEBRICK' 
      
      string3 = ', VARLOCATION=(NODAL)'

      string2 = trim(string2) // trim(string3) // ', DATAPACKING=BLOCK'
  end select
  
  OutputTecplotZoneHeader = trim(string) // string2

end function OutputTecplotZoneHeader

! ************************************************************************** !
!
! WriteTecplotGeomechGridVertices: Prints zone header to a tecplot file
! author: Satish Karra, LANL
! date: 07/2/13
!
! ************************************************************************** !
subroutine WriteTecplotGeomechGridVertices(fid,geomech_realization)

  use Geomechanics_Realization_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Option_module
  use Geomechanics_Patch_module
  use Variables_module

  implicit none

  PetscInt :: fid
  type(geomech_realization_type) :: geomech_realization 
  
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch 
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  PetscInt :: local_size
  PetscErrorCode :: ierr
  
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%nmax_node, &
                    global_vertex_vec,ierr)
  call VecGetLocalSize(global_vertex_vec,local_size,ierr)
  call GetVertexCoordinatesGeomech(grid,global_vertex_vec,X_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
  call WriteTecplotDataSetGeomech(fid,geomech_realization,vec_ptr, &
                                  TECPLOT_REAL,local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

  call GetVertexCoordinatesGeomech(grid,global_vertex_vec,Y_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
  call WriteTecplotDataSetGeomech(fid,geomech_realization,vec_ptr, &
                                  TECPLOT_REAL,local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

  call GetVertexCoordinatesGeomech(grid,global_vertex_vec,Z_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
  call WriteTecplotDataSetGeomech(fid,geomech_realization,vec_ptr, &
                                  TECPLOT_REAL,local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

  call VecDestroy(global_vertex_vec,ierr)

end subroutine WriteTecplotGeomechGridVertices

! ************************************************************************** !
!
! GetVertexCoordinatesGeomech: Extracts vertex coordinates of cells into 
!                               a PetscVec
! author: Satish Karra, LANL
! date: 07/02/2013
!
! ************************************************************************** !
subroutine GetVertexCoordinatesGeomech(grid,vec,direction,option)

  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Variables_module, only : X_COORDINATE, Y_COORDINATE, Z_COORDINATE
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(geomech_grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  type(option_type) :: option
  
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: indices(:)
  PetscReal, allocatable :: values(:)
  PetscErrorCode :: ierr
  
  if (option%mycommsize == 1) then
    call VecGetArrayF90(vec,vec_ptr,ierr)
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          vec_ptr(ivertex) = grid%nodes(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          vec_ptr(ivertex) = grid%nodes(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          vec_ptr(ivertex) = grid%nodes(ivertex)%z
        enddo
    end select
    call VecRestoreArrayF90(vec,vec_ptr,ierr)
  else
    ! initialize to -999 to catch bugs
    call VecSet(vec,-999.d0,ierr)
    allocate(values(grid%nlmax_node))
    allocate(indices(grid%nlmax_node))
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          values(ivertex) = grid%nodes(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          values(ivertex) = grid%nodes(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          values(ivertex) = grid%nodes(ivertex)%z
        enddo
    end select
    indices(:) = grid%node_ids_local_natural(:)-1
    call VecSetValues(vec,grid%nlmax_node, &
                      indices,values,INSERT_VALUES,ierr)
    call VecAssemblyBegin(vec,ierr)
    deallocate(values)
    deallocate(indices)
    call VecAssemblyEnd(vec,ierr)
  endif
  
end subroutine GetVertexCoordinatesGeomech

! ************************************************************************** !
!
! WriteTecplotDataSetGeomech: Writes data from an array within a block
!                      of a Tecplot file
! author: Satish Karra
! date: 07/02//13
!
! ************************************************************************** !
subroutine WriteTecplotDataSetGeomech(fid,geomech_realization,array,datatype, &
                                      size_flag)

  use Geomechanics_Realization_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Patch_module

  implicit none

  PetscInt :: fid
  type(geomech_realization_type) :: geomech_realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size

  PetscInt, parameter :: num_per_line = 10

  call WriteTecplotDataSetNumPerLineGeomech(fid,geomech_realization,array, &
                                            datatype,size_flag,num_per_line) 
  
end subroutine WriteTecplotDataSetGeomech

! ************************************************************************** !
!
! WriteTecplotDataSetNumPerLine: Writes data from an array within a block
!                                of a Tecplot file with a specified number
!                                of values per line
! author: Glenn Hammond
! date: 10/25/07, 12/02/11, Satish Karra 07/02/13
!
! ************************************************************************** !
subroutine WriteTecplotDataSetNumPerLineGeomech(fid,geomech_realization, &
                                            array,datatype, &
                                            size_flag,num_per_line)

  use Geomechanics_Realization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Patch_module

  implicit none
  
  PetscInt :: fid
  type(geomech_realization_type) :: geomech_realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  PetscInt :: num_per_line
  
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch  
  PetscInt :: i
  PetscInt :: max_proc, max_proc_prefetch
  PetscMPIInt :: iproc_mpi, recv_size_mpi
  PetscInt :: max_local_size
  PetscMPIInt :: local_size_mpi
  PetscInt :: istart, iend, num_in_array
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)
  PetscErrorCode :: ierr  
  
1000 format(100(i2,1x))
1001 format(100(i4,1x))
1002 format(100(i6,1x))
1003 format(100(i8,1x))
1004 format(100(i10,1x))
1010 format(100(es13.6,1x))

  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option

  call PetscLogEventBegin(geomech_logging%event_output_write_tecplot,ierr)    

  ! if num_per_line exceeds 100, need to change the format statement below
  if (num_per_line > 100) then
    option%io_buffer = 'Number of values to be written to line in ' // &
      'WriteTecplotDataSetNumPerLine() exceeds 100.  ' // &
      'Must fix format statements.'
    call printErrMsg(option)
  endif

  ! maximum number of initial messages  
#define HANDSHAKE  
  max_proc = option%io_handshake_buffer_size
  max_proc_prefetch = option%io_handshake_buffer_size / 10

  if (size_flag /= 0) then
    call MPI_Allreduce(size_flag,max_local_size,ONE_INTEGER_MPI,MPIU_INTEGER, &
                       MPI_MAX,option%mycomm,ierr)
    local_size_mpi = size_flag
  else 
  ! if first time, determine the maximum size of any local array across 
  ! all procs
    if (max_local_size_saved < 0) then
      call MPI_Allreduce(grid%nlmax_node,max_local_size,ONE_INTEGER_MPI, &
                         MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
      max_local_size_saved = max_local_size
      write(option%io_buffer,'("max_local_size_saved: ",i9)') max_local_size
      call printMsg(option)
    endif
    max_local_size = max_local_size_saved
    local_size_mpi = grid%nlmax_node
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == TECPLOT_INTEGER) then
    allocate(integer_data(max_local_size+10))
    allocate(integer_data_recv(max_local_size))
    do i=1,local_size_mpi
      integer_data(i) = int(array(i))
    enddo
  else
    allocate(real_data(max_local_size+10))
    allocate(real_data_recv(max_local_size))
    do i=1,local_size_mpi
      real_data(i) = array(i)
    enddo
  endif
  
  ! communicate data to processor 0, round robin style
  if (option%myrank == option%io_rank) then
    if (datatype == TECPLOT_INTEGER) then
      ! This approach makes output files identical, regardless of processor
      ! distribution.  It is necessary when diffing files.
      iend = 0
      do
        istart = iend+1
        if (iend+num_per_line > local_size_mpi) exit
        iend = istart+(num_per_line-1)
        i = abs(maxval(integer_data(istart:iend)))
        if (i < 10) then
          write(fid,1000) integer_data(istart:iend)
        else if (i < 1000) then
          write(fid,1001) integer_data(istart:iend)
        else if (i < 100000) then
          write(fid,1002) integer_data(istart:iend)
        else if (i < 10000000) then
          write(fid,1003) integer_data(istart:iend)
        else
          write(fid,1004) integer_data(istart:iend)
        endif
      enddo
      ! shift remaining data to front of array
      integer_data(1:local_size_mpi-iend) = integer_data(iend+1:local_size_mpi)
      num_in_array = local_size_mpi-iend
    else
      iend = 0
      do
        istart = iend+1
        if (iend+num_per_line > local_size_mpi) exit
        iend = istart+(num_per_line-1)
        ! if num_per_line exceeds 100, need to change the format statement below
        write(fid,1010) real_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      real_data(1:local_size_mpi-iend) = real_data(iend+1:local_size_mpi)
      num_in_array = local_size_mpi-iend
    endif
    do iproc_mpi=1,option%mycommsize-1
#ifdef HANDSHAKE    
      if (option%io_handshake_buffer_size > 0 .and. &
          iproc_mpi+max_proc_prefetch >= max_proc) then
        max_proc = max_proc + option%io_handshake_buffer_size
        call MPI_Bcast(max_proc,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                       option%mycomm,ierr)
      endif
#endif      
      call MPI_Probe(iproc_mpi,MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
      recv_size_mpi = status_mpi(MPI_TAG)
      if (datatype == TECPLOT_INTEGER) then
        call MPI_Recv(integer_data_recv,recv_size_mpi,MPIU_INTEGER,iproc_mpi, &
                      MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
        if (recv_size_mpi > 0) then
          integer_data(num_in_array+1:num_in_array+recv_size_mpi) = &
                                             integer_data_recv(1:recv_size_mpi)
          num_in_array = num_in_array+recv_size_mpi
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+num_per_line > num_in_array) exit
          iend = istart+(num_per_line-1)
          i = abs(maxval(integer_data(istart:iend)))
          if (i < 10) then
            write(fid,1000) integer_data(istart:iend)
          else if (i < 1000) then
            write(fid,1001) integer_data(istart:iend)
          else if (i < 100000) then
            write(fid,1002) integer_data(istart:iend)
          else if (i < 10000000) then
            write(fid,1003) integer_data(istart:iend)
          else
            write(fid,1004) integer_data(istart:iend)
          endif
        enddo
        if (iend > 0) then
          integer_data(1:num_in_array-iend) = integer_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      else
        call MPI_Recv(real_data_recv,recv_size_mpi,MPI_DOUBLE_PRECISION,iproc_mpi, &
                      MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
        if (recv_size_mpi > 0) then
          real_data(num_in_array+1:num_in_array+recv_size_mpi) = &
                                             real_data_recv(1:recv_size_mpi)
          num_in_array = num_in_array+recv_size_mpi
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+num_per_line > num_in_array) exit
          iend = istart+(num_per_line-1)
          ! if num_per_line exceeds 100, need to change the format statement below
          write(fid,1010) real_data(istart:iend)
        enddo
        if (iend > 0) then
          real_data(1:num_in_array-iend) = real_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      endif
    enddo
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      max_proc = -1
      call MPI_Bcast(max_proc,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    endif
#endif      
    ! Print the remaining values, if they exist
    if (datatype == TECPLOT_INTEGER) then
      if (num_in_array > 0) then
        i = abs(maxval(integer_data(1:num_in_array)))
        if (i < 10) then
          write(fid,1000) integer_data(1:num_in_array)
        else if (i < 1000) then
          write(fid,1001) integer_data(1:num_in_array)
        else if (i < 100000) then
          write(fid,1002) integer_data(1:num_in_array)
        else if (i < 10000000) then
          write(fid,1003) integer_data(1:num_in_array)
        else
          write(fid,1004) integer_data(1:num_in_array)
        endif
      endif
    else
      if (num_in_array > 0) &
        write(fid,1010) real_data(1:num_in_array)
    endif
  else
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        if (option%myrank < max_proc) exit
        call MPI_Bcast(max_proc,1,MPIU_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
      enddo
    endif
#endif    
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,local_size_mpi,MPIU_INTEGER,option%io_rank, &
                    local_size_mpi,option%mycomm,ierr)
    else
      call MPI_Send(real_data,local_size_mpi,MPI_DOUBLE_PRECISION,option%io_rank, &
                    local_size_mpi,option%mycomm,ierr)
    endif
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        call MPI_Bcast(max_proc,1,MPIU_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
        if (max_proc < 0) exit
      enddo
    endif
#endif
#undef HANDSHAKE
  endif
      
  if (datatype == TECPLOT_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

  call PetscLogEventEnd(geomech_logging%event_output_write_tecplot,ierr)    

end subroutine WriteTecplotDataSetNumPerLineGeomech


end module Output_Geomechanics_module

#endif
