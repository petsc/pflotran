module Output_module

  use Logging_module 
  use Output_Aux_module
  
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petsclog.h"

#if defined(PARALLELIO_LIB_WRITE)
  include "piof.h"
#endif

  PetscInt, parameter :: TECPLOT_INTEGER = 0
  PetscInt, parameter :: TECPLOT_REAL = 1

  PetscInt, parameter :: VTK_INTEGER = 0
  PetscInt, parameter :: VTK_REAL = 1

  PetscInt, parameter :: TECPLOT_FILE = 0
  PetscInt, parameter ::  HDF5_FILE = 1

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  PetscInt, save :: max_local_size_saved = -1
  
  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: observation_first
  PetscBool :: hdf5_first
  PetscBool :: mass_balance_first
  PetscBool :: hydrograph_first

  interface Output
    module procedure Output1
#ifdef SURFACE_FLOW
    module procedure Output2
#endif
  end interface

  interface OutputTecplotHeader
    module procedure OutputTecplotHeader1
#ifdef SURFACE_FLOW
    module procedure OutputTecplotHeader2
#endif
  end interface

  interface OutputTecplotZoneHeader
    module procedure OutputTecplotZoneHeader1
#ifdef SURFACE_FLOW
    module procedure OutputTecplotZoneHeader2
#endif
  end interface

  interface OutputGetVarFromArray
    module procedure OutputGetVarFromArray1
#ifdef SURFACE_FLOW
    module procedure OutputGetVarFromArray2
#endif
  end interface

  interface WriteTecplotDataSet
    module procedure WriteTecplotDataSet1
#ifdef SURFACE_FLOW
    module procedure WriteTecplotDataSet2
#endif
  end interface

  interface WriteTecplotDataSetNumPerLine
    module procedure WriteTecplotDataSetNumPerLine1
#ifdef SURFACE_FLOW
    module procedure WriteTecplotDataSetNumPerLine2
#endif
  end interface

  interface WriteTecplotDataSetFromVec
    module procedure WriteTecplotDataSetFromVec1
#ifdef SURFACE_FLOW
    module procedure WriteTecplotDataSetFromVec2
#endif
  end interface

  interface WriteTecplotUGridElements
    module procedure WriteTecplotUGridElements1
#ifdef SURFACE_FLOW
    module procedure WriteTecplotUGridElements2
#endif
  end interface

  interface WriteTecplotUGridVertices
    module procedure WriteTecplotUGridVertices1
#ifdef SURFACE_FLOW
    module procedure WriteTecplotUGridVertices2
#endif
  end interface

  public :: OutputInit, Output, OutputVectorTecplot, &
            OutputObservation, OutputGetVarFromArray, &
            OutputPermeability, OutputPrintCouplers, &
            OutputGetCellCenteredVelocities

contains

! ************************************************************************** !
!
! OutputInit: Initializes variables
! author: Glenn Hammond
! date: 01/22/09
!
! ************************************************************************** !
subroutine OutputInit(realization,num_steps)

  use Realization_module
  use Option_module

  implicit none
  
  type(realization_type) :: realization
  PetscInt :: num_steps
  
  ! set size to -1 in order to re-initialize parallel communication blocks
  max_local_size_saved = -1

  if (num_steps == 0) then
    observation_first = PETSC_TRUE
    hdf5_first = PETSC_TRUE
    mass_balance_first = PETSC_TRUE
    hydrograph_first = PETSC_TRUE
  else
    observation_first = PETSC_FALSE
    hdf5_first = PETSC_FALSE
    mass_balance_first = PETSC_FALSE
    hydrograph_first = PETSC_FALSE
  endif

end subroutine OutputInit

! ************************************************************************** !
!
! Output: Main driver for all output subroutines
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine Output1(realization,plot_flag,transient_plot_flag)

  use Realization_module, only : realization_type
  use Option_module, only : OptionCheckTouch, option_type, printMsg
  
  implicit none
  
  type(realization_type) :: realization
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend
  type(option_type), pointer :: option

  option => realization%option

#ifdef VAMSI_STAGE_BARRIER
  ! barrier to calculate the accurate timing of Output Stage
  call MPI_Barrier(option%mycomm,ierr)
#endif 

  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr)

  ! check for plot request from active directory
  if (.not.plot_flag) then

    if (option%use_touch_options) then
      string = 'plot'
      if (OptionCheckTouch(option,string)) then
        realization%output_option%plot_name = 'plot'
        plot_flag = PETSC_TRUE
      endif
    endif

  endif

  if (plot_flag) then
  
    if (realization%output_option%print_hdf5) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_hdf5,ierr)    
      call OutputHDF5(realization)
      call PetscLogEventEnd(logging%event_output_hdf5,ierr)    
      call PetscGetTime(tend,ierr)
#ifdef PARALLELIO_LIB_WRITE
      if (option%myrank == 0) write (*,'(" Parallel IO Write method is used in & 
                                          writing the output, HDF5_WRITE_GROUP_SIZE = ",i5)') option%hdf5_write_group_size
#endif
#ifdef VAMSI_HDF5_WRITE
      if (option%myrank == 0) write (*,'(" Vamsi''s HDF5 method is used in & 
                                          writing the output, HDF5_WRITE_GROUP_SIZE = ",i5)') option%hdf5_write_group_size
#endif      
      write(option%io_buffer,'(f10.2," Seconds to write HDF5 file.")') tend-tstart
      call printMsg(option)
    endif
   
    if (realization%output_option%print_tecplot) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_tecplot,ierr) 
      select case(realization%output_option%tecplot_format)
        case (TECPLOT_POINT_FORMAT)
          call OutputTecplotPoint(realization)
        case (TECPLOT_BLOCK_FORMAT,TECPLOT_FEBRICK_FORMAT)
          call OutputTecplotBlock(realization)
      end select
      call PetscLogEventEnd(logging%event_output_tecplot,ierr)    
      call PetscGetTime(tend,ierr) 
      write(option%io_buffer,'(f10.2," Seconds to write to Tecplot file(s)")') &
            tend-tstart
      call printMsg(option)        
    endif

    if (realization%output_option%print_vtk) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_vtk,ierr) 
      call OutputVTK(realization)

      call PetscLogEventEnd(logging%event_output_vtk,ierr)    
      call PetscGetTime(tend,ierr) 
      write(option%io_buffer,'(f10.2," Seconds to write to VTK file(s)")') &
            tend-tstart
      call printMsg(option) 
    endif
      
    if (realization%output_option%print_mad) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_mad,ierr) 
      call OutputMAD(realization)

      call PetscLogEventEnd(logging%event_output_mad,ierr)    
      call PetscGetTime(tend,ierr) 
      write(option%io_buffer,'(f10.2," Seconds to write to MAD HDF5 file(s)")') &
            tend-tstart
      call printMsg(option) 
    endif
      
    if (option%compute_statistics) then
      call ComputeFlowCellVelocityStats(realization)
      call ComputeFlowFluxVelocityStats(realization)
    endif
  
    realization%output_option%plot_number = realization%output_option%plot_number + 1

  endif
  
  if (transient_plot_flag) then
    if (option%compute_mass_balance_new) then
      call OutputMassBalance(realization)
    endif
    call OutputObservation(realization)
  endif
  
  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  realization%output_option%plot_name = ''

#ifdef VAMSI_STAGE_BARRIER
  call MPI_Barrier(option%mycomm,ierr)
  ! barrier to calculate the accurate timing of Output Stage
#endif 

  call PetscLogStagePop(ierr)
  
end subroutine Output1

! ************************************************************************** !
!
! OutputObservation: Main driver for all observation output subroutines
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine OutputObservation(realization)
                           ! for some flakey reason, current Intel 10.1 reports
                           ! error if 'only' statement not used.
  use Realization_module, only : realization_type 
  use Option_Module
  
  implicit none
  
  type(realization_type) :: realization

!  if (realization%output_option%print_hdf5) then
!    call OutputObservationHDF5(realization)
!    call OutputObservationTecplot(realization)
!  endif
 
!  if (realization%output_option%print_tecplot .or. &
!      realization%output_option%print_hdf5) then
  if (realization%output_option%print_observation) then
    call OutputObservationTecplot(realization)
  endif

end subroutine OutputObservation

! ************************************************************************** !
!
! OutputFilenameID: Creates an ID for filename
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFilenameID(output_option,option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option

  character(len=MAXWORDLENGTH) :: OutputFilenameID
  
  if (output_option%plot_number < 10) then
    write(OutputFilenameID,'("00",i1)') output_option%plot_number  
  else if (output_option%plot_number < 100) then
    write(OutputFilenameID,'("0",i2)') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(OutputFilenameID,'(i3)') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(OutputFilenameID,'(i4)') output_option%plot_number  
  endif 
  
  OutputFilenameID = adjustl(OutputFilenameID)

end function OutputFilenameID

! ************************************************************************** !
!
! OutputFilename: Creates a filename for a Tecplot file
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFilename(output_option,option,suffix,optional_string)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option
  character(len=*) :: suffix
  character(len=*) :: optional_string
  
  character(len=MAXSTRINGLENGTH) :: OutputFilename

  character(len=MAXWORDLENGTH) :: final_suffix
  character(len=MAXSTRINGLENGTH) :: final_optional_string


  if (len_trim(optional_string) > 0) then
    final_optional_string = '-' // optional_string
  else
    final_optional_string = ''
  endif
  final_suffix = '.' // suffix
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    OutputFilename = trim(output_option%plot_name) // &
            trim(final_optional_string) // &
            final_suffix
  else  
    OutputFilename = trim(option%global_prefix) // &
            trim(option%group_prefix) // &
            trim(final_optional_string) // &
            '-' // &
            trim(OutputFilenameID(output_option,option)) // &
            final_suffix
  endif
  
end function OutputFilename

! ************************************************************************** !
!
! OutputTecplotHeader: Print header to Tecplot file
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
subroutine OutputTecplotHeader1(fid,realization,icolumn)

  use Realization_module
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module

  use Mphase_module
  use Immis_module
  use THC_module
  use THMC_module
  use Richards_module
  use Flash2_module
  use Miscible_module
  use General_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
  
  implicit none

  PetscInt :: fid
  type(realization_type) :: realization
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
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  output_option => realization%output_option

  ! write header
  ! write title
  write(fid,'(''TITLE = "'',1es13.5," [",a1,'']"'')') &
                option%time/output_option%tconv,output_option%tunit

  ! initial portion of header
  header = 'VARIABLES=' // &
            '"X [m]",' // &
            '"Y [m]",' // &
            '"Z [m]"'

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
  string = OutputTecplotZoneHeader(realization,variable_count, &
                                   output_option%tecplot_format)
  write(fid,'(a)') trim(string)

end subroutine OutputTecplotHeader1

! ************************************************************************** !
!
! OutputTecplotZoneHeader: Print zone header to Tecplot file
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputTecplotZoneHeader1(realization,variable_count,tecplot_format)

  use Realization_module
  use Grid_module
  use Option_module
  
  implicit none

  type(realization_type) :: realization
  PetscInt :: variable_count
  PetscInt :: tecplot_format
  
  character(len=MAXSTRINGLENGTH) :: OutputTecplotZoneHeader1

  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  grid => realization%patch%grid
  option => realization%option
  output_option => realization%output_option  

  string = 'ZONE T="' // &
           trim(OutputFormatDouble(option%time/output_option%tconv)) // &
           '"'
  string2 = ''
  select case(tecplot_format)
    case (TECPLOT_POINT_FORMAT)
      if ((realization%discretization%itype == STRUCTURED_GRID).or. &
          (realization%discretization%itype == STRUCTURED_GRID_MIMETIC)) then
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
      if ((realization%discretization%itype == STRUCTURED_GRID).or. &
          (realization%discretization%itype == STRUCTURED_GRID_MIMETIC)) then
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
        string2 = trim(string2) // ', ZONETYPE=FEBRICK'
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
  
  OutputTecplotZoneHeader1 = trim(string) // string2  

end function OutputTecplotZoneHeader1

! ************************************************************************** !
!
! OutputTecplotBlock: Print to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !  
subroutine OutputTecplotBlock(realization)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Field_module
  use Patch_module
  
  use Mphase_module
  use Immis_module
  use THC_module
  use THMC_module
  use Richards_module
  use Flash2_module
  use Miscible_module
  use General_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  PetscInt, parameter :: icolumn = -1
  character(len=MAXSTRINGLENGTH) :: filename, string, string2
  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch 
  type(reaction_type), pointer :: reaction 
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option
  
  filename = OutputFilename(output_option,option,'tec','')
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=filename,action="write")
    call OutputTecplotHeader(OUTPUT_UNIT,realization,icolumn)
  endif
    
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)  

  ! write out coordinates
  if (realization%discretization%itype == STRUCTURED_GRID .or. &
      realization%discretization%itype == STRUCTURED_GRID_MIMETIC ) then
    call WriteTecplotStructuredGrid(OUTPUT_UNIT,realization)
  else
    call WriteTecplotUGridVertices(OUTPUT_UNIT,realization)
  endif

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

  if (realization%discretization%itype == UNSTRUCTURED_GRID .and. &
      realization%discretization%grid%itype == &
      IMPLICIT_UNSTRUCTURED_GRID)  then
    call WriteTecplotUGridElements(OUTPUT_UNIT,realization)
  endif
  
  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)
  
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesTecplotBlock(realization)
  endif
  
  if (output_option%print_tecplot_flux_velocities) then
    if (grid%structured_grid%nx > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization,LIQUID_PHASE, &
                                          X_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization,GAS_PHASE, &
                                              X_DIRECTION)
      end select
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization,LIQUID_PHASE, &
                                          Y_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE, IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization,GAS_PHASE, &
                                              Y_DIRECTION)
      end select
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization,LIQUID_PHASE, &
                                          Z_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE, IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization,GAS_PHASE, &
                                              Z_DIRECTION)
      end select
    endif
  endif
      
end subroutine OutputTecplotBlock

! ************************************************************************** !
!
! OutputVelocitiesTecplotBlock: Print velocities to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVelocitiesTecplotBlock(realization)
 
  use Realization_module
  use Discretization_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Field_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  Vec :: natural_vec

  PetscReal, pointer :: vec_ptr(:)
  
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  option => realization%option
  output_option => realization%output_option
  discretization => realization%discretization

  filename = OutputFilename(output_option,option,'tec','vel')
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot velocity output file: ' // &
                       trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=filename,action="write")
  
    ! write header
    ! write title
    write(OUTPUT_UNIT,'(''TITLE = "'',1es13.5," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",' // &
             '"vlx [m/' // trim(output_option%tunit) // ']",' // &
             '"vly [m/' // trim(output_option%tunit) // ']",' // &
             '"vlz [m/' // trim(output_option%tunit) // ']"'
    if (option%nphase > 1) then
      string = trim(string) // &
               ',"vgx [m/' // trim(output_option%tunit) // ']",' // &
               '"vgy [m/' // trim(output_option%tunit) // ']",' // &
               '"vgz [m/' // trim(output_option%tunit) // ']"'
    endif

    string = trim(string) // ',"Material_ID"'
    write(OUTPUT_UNIT,'(a)') trim(string)
  
    if (option%nphase > 1) then
      string = OutputTecplotZoneHeader(realization,TEN_INTEGER, &
                                       TECPLOT_BLOCK_FORMAT)
    else
      string = OutputTecplotZoneHeader(realization,SEVEN_INTEGER, &
                                       TECPLOT_BLOCK_FORMAT)
    endif
    write(OUTPUT_UNIT,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)    

  ! write out coorindates
  if (realization%discretization%itype == STRUCTURED_GRID .or. &
      realization%discretization%itype == STRUCTURED_GRID_MIMETIC)  then
    call WriteTecplotStructuredGrid(OUTPUT_UNIT,realization)
  else
    call WriteTecplotUGridVertices(OUTPUT_UNIT,realization)
  endif
  
  call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec,TECPLOT_REAL)

  call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec,TECPLOT_REAL)

  call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec,TECPLOT_REAL)

  if (option%nphase > 1) then
    call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec,TECPLOT_REAL)

    call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec,TECPLOT_REAL)

    call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec,TECPLOT_REAL)
  endif

  ! material id
  call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization,natural_vec,TECPLOT_INTEGER)
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if (realization%discretization%itype == UNSTRUCTURED_GRID .and. &
      realization%discretization%grid%itype == &
      IMPLICIT_UNSTRUCTURED_GRID)  then
    call WriteTecplotUGridElements(OUTPUT_UNIT,realization)
  endif

  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)
  
end subroutine OutputVelocitiesTecplotBlock

! ************************************************************************** !
!
! OutputFluxVelocitiesTecplotBlk: Print intercellular fluxes to Tecplot file 
!                                 in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputFluxVelocitiesTecplotBlk(realization,iphase, &
                                            direction)
!geh - specifically, the flow velocities at the interfaces between cells
 
  use Realization_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  PetscInt :: iphase
  PetscInt :: direction
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: discretization  
  type(output_option_type), pointer :: output_option
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscInt :: local_size, global_size
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscInt :: i, j, k
  PetscInt :: local_id, ghosted_id
  PetscInt :: adjusted_size
  PetscInt :: count, iconn, sum_connection
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: array(:)
  PetscInt, allocatable :: indices(:)
  Vec :: global_vec, global_vec2
  PetscReal :: sum, average, max, min , std_dev
  PetscInt :: max_loc, min_loc

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
    
  nullify(array)

  call PetscLogEventBegin(logging%event_output_write_flux_tecplot,ierr) 
                          
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-'
  else  
    filename = trim(option%global_prefix) // trim(option%group_prefix) // '-'
  endif
  
  select case(iphase)
    case(LIQUID_PHASE)
      filename = trim(filename) // 'vl'
    case(GAS_PHASE)
      filename = trim(filename) // 'vg'
  end select
  
  select case(direction)
    case(X_DIRECTION)
      filename = trim(filename) // 'x'
    case(Y_DIRECTION)
      filename = trim(filename) // 'y'
    case(Z_DIRECTION)
      filename = trim(filename) // 'z'
  end select 
  
  string = trim(OutputFilenameID(output_option,option))
  
  filename = trim(filename) // '-' // trim(string) // '.tec'
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot velocity flux output file: ' // &
                       trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=filename,action="write")
  
    ! write header
    ! write title
    write(OUTPUT_UNIT,'(''TITLE = "'',1es13.5," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",'
    select case(iphase)
      case(LIQUID_PHASE)
        string = trim(string) // '"Liquid'
      case(GAS_PHASE)
        string = trim(string) // '"Gas'
    end select
  
    select case(direction)
      case(X_DIRECTION)
        string = trim(string) // ' vlx [m/' // trim(output_option%tunit) // ']"'
      case(Y_DIRECTION)
        string = trim(string) // ' vly [m/' // trim(output_option%tunit) // ']"'
      case(Z_DIRECTION)
        string = trim(string) // ' vlz [m/' // trim(output_option%tunit) // ']"'
    end select 
    
    write(OUTPUT_UNIT,'(a)') trim(string)
  
    ! write zone header
    select case(direction)
      case(X_DIRECTION)
        write(string,'(''ZONE T= "'',1es13.5,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4)') &
                     option%time/output_option%tconv,grid%structured_grid%nx-1,grid%structured_grid%ny,grid%structured_grid%nz 
      case(Y_DIRECTION)
        write(string,'(''ZONE T= "'',1es13.5,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4)') &
                     option%time/output_option%tconv,grid%structured_grid%nx,grid%structured_grid%ny-1,grid%structured_grid%nz 
      case(Z_DIRECTION)
        write(string,'(''ZONE T= "'',1es13.5,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4)') &
                     option%time/output_option%tconv,grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz-1
    end select 
    string = trim(string) // ', DATAPACKING=BLOCK'
    write(OUTPUT_UNIT,'(a)') trim(string)

  endif
  
  ! write blocks'
  
  ! face coordinates
  local_size = grid%nlmax
  global_size = grid%nmax
!GEH - Structured Grid Dependence - Begin
  nx_local = grid%structured_grid%nlx
  ny_local = grid%structured_grid%nly
  nz_local = grid%structured_grid%nlz
  nx_global = grid%structured_grid%nx
  ny_global = grid%structured_grid%ny
  nz_global = grid%structured_grid%nz
  select case(direction)
    case(X_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%ny*grid%structured_grid%nz
      nx_global = grid%structured_grid%nx-1
      if (grid%structured_grid%gxe-grid%structured_grid%lxe == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlyz
        nx_local = grid%structured_grid%nlx-1
      endif
    case(Y_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%nx*grid%structured_grid%nz
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%gye-grid%structured_grid%lye == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlxz
        ny_local = grid%structured_grid%nly-1
      endif
    case(Z_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%nxy
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%gze-grid%structured_grid%lze == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlxy
        nz_local = grid%structured_grid%nlz-1
      endif
  end select  
  allocate(indices(local_size))

  ! fill indices array with natural ids in newly sized array
  count = 0
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        indices(count) = i+grid%structured_grid%lxs+(j-1+grid%structured_grid%lys)*nx_global+ &
                         (k-1+grid%structured_grid%lzs)*nx_global*ny_global
      enddo
    enddo
  enddo
  
  ! X-coordinates
  count = 0
  allocate(array(local_size))
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        ghosted_id = grid%nL2G(local_id)
        array(count) = grid%x(ghosted_id)
        if (direction == X_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dx(ghosted_id)
      enddo
    enddo
  enddo
  ! warning: adjusted size will be changed in ConvertArrayToNatural
  ! thus, you cannot pass in local_size, since it is needed later
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size,option)
  call WriteTecplotDataSet(OUTPUT_UNIT,realization,array,TECPLOT_REAL,adjusted_size)
  ! since the array has potentially been resized, must reallocate
  deallocate(array)
  nullify(array)

  ! Y-coordinates
  count = 0
  allocate(array(local_size))
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        ghosted_id = grid%nL2G(local_id)        
        array(count) = grid%y(ghosted_id)
        if (direction == Y_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dy(ghosted_id)
      enddo
    enddo
  enddo
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size,option)
  call WriteTecplotDataSet(OUTPUT_UNIT,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)

  ! Z-coordinates
  count = 0
  allocate(array(local_size))
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        ghosted_id = grid%nL2G(local_id)        
        array(count) = grid%z(ghosted_id)
        if (direction == Z_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dz(ghosted_id)
      enddo
    enddo
  enddo
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size,option)
  call WriteTecplotDataSet(OUTPUT_UNIT,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)

  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option) 
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

  ! write out data set 
  count = 0 
  allocate(array(local_size)) 
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

!GEH - Structured Grid Dependence - End
  
  array(1:local_size) = array(1:local_size)*output_option%tconv ! convert time units
  
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size,option)
  call WriteTecplotDataSet(OUTPUT_UNIT,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  deallocate(indices)

  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)

  call PetscLogEventEnd(logging%event_output_write_flux_tecplot,ierr) 
  
end subroutine OutputFluxVelocitiesTecplotBlk

! ************************************************************************** !
!
! OutputTecplotPoint: Print to Tecplot file in POINT format
! author: Glenn Hammond
! date: 11/03/08
!
! ************************************************************************** !  
subroutine OutputTecplotPoint(realization)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Structured_Grid_module
  use Option_module
  use Field_module
  use Patch_module
  
  use Mphase_module
  use Immis_module
  use THC_module
  use THMC_module
  use Richards_module
  use Flash2_module
  use Miscible_module
  use General_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  PetscInt :: icolumn
  character(len=MAXSTRINGLENGTH) :: filename, string
  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch 
  type(reaction_type), pointer :: reaction 
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: value
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option

  filename = OutputFilename(output_option,option,'tec','')
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // &
                       trim(filename)
    call printMsg(option)                       
    open(unit=OUTPUT_UNIT,file=filename,action="write")
  
    if (output_option%print_column_ids) then
      icolumn = 3
    else
      icolumn = -1
    endif
    call OutputTecplotHeader(OUTPUT_UNIT,realization,icolumn)
  endif
  
1000 format(es13.6,1x)
1001 format(i4,1x)
1009 format('')

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    write(OUTPUT_UNIT,1000,advance='no') grid%x(ghosted_id)
    write(OUTPUT_UNIT,1000,advance='no') grid%y(ghosted_id)
    write(OUTPUT_UNIT,1000,advance='no') grid%z(ghosted_id)

    ! loop over variables and write to file
    cur_variable => output_option%output_variable_list%first
    do
      if (.not.associated(cur_variable)) exit
      value = RealizGetDatasetValueAtCell(realization,cur_variable%ivar, &
                                          cur_variable%isubvar,ghosted_id)
      if (cur_variable%iformat == 0) then
        write(OUTPUT_UNIT,1000,advance='no') value
      else
        write(OUTPUT_UNIT,1001,advance='no') int(value)
      endif
      cur_variable => cur_variable%next
    enddo

    write(OUTPUT_UNIT,1009) 

  enddo
  
  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)
  
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesTecplotPoint(realization)
  endif
  
end subroutine OutputTecplotPoint

! ************************************************************************** !
!
! OutputVelocitiesTecplotPoint: Print velocities to Tecplot file in POINT format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVelocitiesTecplotPoint(realization)
 
  use Realization_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: value  
  Vec :: global_vec_vx, global_vec_vy, global_vec_vz

  PetscReal, pointer :: vec_ptr_vx(:), vec_ptr_vy(:), vec_ptr_vz(:)
  
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  option => realization%option
  output_option => realization%output_option
  discretization => realization%discretization
  
  filename = OutputFilename(output_option,option,'tec','vel')
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot velocity output file: ' // &
                       trim(filename)
    call printMsg(option)                       
    open(unit=OUTPUT_UNIT,file=filename,action="write")
  
    ! write header
    ! write title
    write(OUTPUT_UNIT,'(''TITLE = "'',1es13.4," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",' // &
             '"vlx [m/' // trim(output_option%tunit) // ']",' // &
             '"vly [m/' // trim(output_option%tunit) // ']",' // &
             '"vlz [m/' // trim(output_option%tunit) // ']"'
    if (option%nphase > 1) then
      string = trim(string) // &
               ',"vgx [m/' // trim(output_option%tunit) // ']",' // &
               '"vgy [m/' // trim(output_option%tunit) // ']",' // &
               '"vgz [m/' // trim(output_option%tunit) // ']"'
    endif
    
    string = trim(string) // ',"Material_ID"'
    write(OUTPUT_UNIT,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es13.5,''",'','' I='',i5,'', J='',i5, &
                 &'', K='',i5)') &
                 option%time/output_option%tconv, &
                 grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
    string = trim(string) // ', DATAPACKING=POINT'
    write(OUTPUT_UNIT,'(a)') trim(string)

  endif
  
  ! currently supported for only liquid phase'
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec_vx,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec_vy,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec_vz,GLOBAL, &
                                  option)  
  
  call OutputGetCellCenteredVelocities(realization,global_vec_vx,LIQUID_PHASE,X_DIRECTION)
  call OutputGetCellCenteredVelocities(realization,global_vec_vy,LIQUID_PHASE,Y_DIRECTION)
  call OutputGetCellCenteredVelocities(realization,global_vec_vz,LIQUID_PHASE,Z_DIRECTION)

  call GridVecGetArrayF90(grid,global_vec_vx,vec_ptr_vx,ierr)
  call GridVecGetArrayF90(grid,global_vec_vy,vec_ptr_vy,ierr)
  call GridVecGetArrayF90(grid,global_vec_vz,vec_ptr_vz,ierr)

  ! write points
1000 format(es13.6,1x)
1001 format(i4,1x)
1002 format(3(es13.6,1x))
1009 format('')

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)  ! local and ghosted are same for non-parallel
    write(OUTPUT_UNIT,1000,advance='no') grid%x(ghosted_id)
    write(OUTPUT_UNIT,1000,advance='no') grid%y(ghosted_id)
    write(OUTPUT_UNIT,1000,advance='no') grid%z(ghosted_id)
    
    write(OUTPUT_UNIT,1000,advance='no') vec_ptr_vx(ghosted_id)
    write(OUTPUT_UNIT,1000,advance='no') vec_ptr_vy(ghosted_id)
    write(OUTPUT_UNIT,1000,advance='no') vec_ptr_vz(ghosted_id)

    ! material id
    value = RealizGetDatasetValueAtCell(realization,MATERIAL_ID, &
                                        ZERO_INTEGER,ghosted_id)
    write(OUTPUT_UNIT,1001,advance='no') int(value)
  
    write(OUTPUT_UNIT,1009)
    
  enddo
  
  call GridVecRestoreArrayF90(grid,global_vec_vx,vec_ptr_vx,ierr)
  call GridVecRestoreArrayF90(grid,global_vec_vy,vec_ptr_vy,ierr)
  call GridVecRestoreArrayF90(grid,global_vec_vz,vec_ptr_vz,ierr)
  
  call VecDestroy(global_vec_vx,ierr)
  call VecDestroy(global_vec_vy,ierr)
  call VecDestroy(global_vec_vz,ierr)

  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)
  
end subroutine OutputVelocitiesTecplotPoint

! ************************************************************************** !
!
! OutputVectorTecplot: Print a vector to a Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVectorTecplot(filename,dataset_name,realization,vector)
 
  use Realization_module
  use Discretization_module
  use Option_module
  use Field_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Patch_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: dataset_name
  type(realization_type) :: realization
  Vec :: vector

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: patch  
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt, parameter :: fid=86

  call PetscLogEventBegin(logging%event_output_vec_tecplot,ierr) 

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  discretization => realization%discretization
  
  ! open file
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // trim(filename)
    call printMsg(option)
    open(unit=fid,file=filename,action="write")
  
    ! write header
    ! write title
    write(fid,'(''TITLE = "PFLOTRAN Vector"'')')
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",'
    string = trim(string) // '"' // trim(dataset_name) // '"'
    string = trim(string) // ',"Material_ID"'
    write(fid,'(a)') trim(string)
  
    !geh: due to pgi bug, cannot embed functions with calls to write() within
    !     write statement
    string = OutputTecplotZoneHeader(realization,FIVE_INTEGER, &
                                     TECPLOT_BLOCK_FORMAT)
    write(fid,'(a)') trim(string)
  endif
  
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)  
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)    

  ! write out coorindates

  if (realization%discretization%itype == STRUCTURED_GRID .or. &
      realization%discretization%itype == STRUCTURED_GRID_MIMETIC)  then
    call WriteTecplotStructuredGrid(fid,realization)
  else  
    call WriteTecplotUGridVertices(fid,realization)
  endif    

  call DiscretizationGlobalToNatural(discretization,vector,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_INTEGER)
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if (realization%discretization%itype == UNSTRUCTURED_GRID .and. &
      realization%discretization%grid%itype == &
      IMPLICIT_UNSTRUCTURED_GRID)  then
    call WriteTecplotUGridElements(fid,realization)
  endif    

  close(fid)

  call PetscLogEventEnd(logging%event_output_vec_tecplot,ierr) 
                            
end subroutine OutputVectorTecplot

! ************************************************************************** !
!
! WriteTecplotDataSetFromVec: Writes data from a Petsc Vec within a block
!                             of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSetFromVec1(fid,realization,vec,datatype)

  use Realization_module
  
  implicit none

  PetscInt :: fid
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: datatype
  
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,realization,vec_ptr,datatype,ZERO_INTEGER) ! 0 implies grid%nlmax
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteTecplotDataSetFromVec1

! ************************************************************************** !
!
! WriteTecplotUGridVertices: Writes unstructured grid vertices
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine WriteTecplotUGridVertices1(fid,realization)

  use Realization_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module

  implicit none

  PetscInt :: fid
  type(realization_type) :: realization 
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch 
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  PetscInt :: local_size
  PetscErrorCode :: ierr
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
    call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                      grid%unstructured_grid%num_vertices_global, &
                      global_vertex_vec,ierr)
    call VecGetLocalSize(global_vertex_vec,local_size,ierr)
    call GetVertexCoordinates(grid, global_vertex_vec,X_COORDINATE,option)
    call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
    call WriteTecplotDataSet(fid,realization,vec_ptr,TECPLOT_REAL, &
                             local_size)
    call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

    call GetVertexCoordinates(grid,global_vertex_vec,Y_COORDINATE,option)
    call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
    call WriteTecplotDataSet(fid,realization,vec_ptr,TECPLOT_REAL, &
                             local_size)
    call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

    call GetVertexCoordinates(grid,global_vertex_vec, Z_COORDINATE,option)
    call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
    call WriteTecplotDataSet(fid,realization,vec_ptr,TECPLOT_REAL, &
                             local_size)
    call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

    call VecDestroy(global_vertex_vec, ierr)
  else
    if (option%myrank == option%io_rank) then
      write(fid,'(">",/,"Add explicit mesh vertex information here",/,">")')
    endif
  endif

end subroutine WriteTecplotUGridVertices1

! ************************************************************************** !
!
! WriteTecplotUGridVertices: Writes unstructured grid elements
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine WriteTecplotUGridElements1(fid,realization)

  use Realization_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module
  
  implicit none

  PetscInt :: fid
  type(realization_type) :: realization

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch 
  Vec :: global_cconn_vec
  type(ugdm_type), pointer :: ugdm_element
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  Vec :: global_vec
  Vec :: natural_vec

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
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
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)
  call WriteTecplotDataSetNumPerLine(fid,realization,vec_ptr, &
                                     TECPLOT_INTEGER, &
                                     grid%unstructured_grid%nlmax*8, &
                                     EIGHT_INTEGER)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)
  call VecDestroy(global_vec,ierr)
  call VecDestroy(natural_vec,ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine WriteTecplotUGridElements1

! ************************************************************************** !
!
! WriteTecplotStructuredGrid: Writes structured grid face coordinates 
! author: Glenn Hammond
! date: 02/26/08
!
! ************************************************************************** !
subroutine WriteTecplotStructuredGrid(fid,realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscInt :: i, j, k, count, nx, ny, nz
  PetscReal :: temp_real

1000 format(es13.6,1x)
1001 format(10(es13.6,1x))
  
  call PetscLogEventBegin(logging%event_output_str_grid_tecplot,ierr) 
                              
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  nx = grid%structured_grid%nx
  ny = grid%structured_grid%ny
  nz = grid%structured_grid%nz
  
  if (option%myrank == option%io_rank) then
    ! x-dir
    count = 0
    do k=1,nz+1
      do j=1,ny+1
        temp_real = realization%discretization%origin(X_DIRECTION)
        write(fid,1000,advance='no') temp_real
        count = count + 1
        if (mod(count,10) == 0) then
          write(fid,'(a)') ""
          count = 0
        endif
        do i=1,nx
          temp_real = temp_real + grid%structured_grid%dx_global(i)
          write(fid,1000,advance='no') temp_real
          count = count + 1
          if (mod(count,10) == 0) then
            write(fid,'(a)') ""
            count = 0
          endif
        enddo
      enddo
    enddo
    if (count /= 0) write(fid,'(a)') ""
    ! y-dir
    count = 0
    do k=1,nz+1
      temp_real = realization%discretization%origin(Y_DIRECTION)
      do i=1,nx+1
        write(fid,1000,advance='no') temp_real
        count = count + 1
        if (mod(count,10) == 0) then
          write(fid,'(a)') ""
          count = 0
        endif
      enddo
      do j=1,ny
        temp_real = temp_real + grid%structured_grid%dy_global(j)
        do i=1,nx+1
          write(fid,1000,advance='no') temp_real
          count = count + 1
          if (mod(count,10) == 0) then
            write(fid,'(a)') ""
            count = 0
          endif
        enddo
      enddo
    enddo
    if (count /= 0) write(fid,'(a)') ""
    ! z-dir
    count = 0
    temp_real = realization%discretization%origin(Z_DIRECTION)
    do i=1,(nx+1)*(ny+1)
      write(fid,1000,advance='no') temp_real
      count = count + 1
      if (mod(count,10) == 0) then
        write(fid,'(a)') ""
        count = 0
      endif
    enddo
    do k=1,nz
      temp_real = temp_real + grid%structured_grid%dz_global(k)
      do j=1,ny+1
        do i=1,nx+1
          write(fid,1000,advance='no') temp_real
          count = count + 1
          if (mod(count,10) == 0) then
            write(fid,'(a)') ""
            count = 0
          endif
        enddo
      enddo
    enddo
    if (count /= 0) write(fid,'(a)') ""

  endif

  call PetscLogEventEnd(logging%event_output_str_grid_tecplot,ierr) 
                            
end subroutine WriteTecplotStructuredGrid

! ************************************************************************** !
!
! WriteTecplotDataSet: Writes data from an array within a block
!                      of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSet1(fid,realization,array,datatype,size_flag)

  use Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size

  PetscInt, parameter :: num_per_line = 10

  call WriteTecplotDataSetNumPerLine(fid,realization,array,datatype, &
                                     size_flag,num_per_line) 
  
end subroutine WriteTecplotDataSet1

! ************************************************************************** !
!
! WriteTecplotDataSetNumPerLine: Writes data from an array within a block
!                                of a Tecplot file with a specified number
!                                of values per line
! author: Glenn Hammond
! date: 10/25/07, 12/02/11
!
! ************************************************************************** !
subroutine WriteTecplotDataSetNumPerLine1(fid,realization,array,datatype, &
                                         size_flag,num_per_line)

  use Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  PetscInt :: num_per_line
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscInt :: i
  PetscInt :: max_proc, max_proc_prefetch
  PetscMPIInt :: iproc_mpi, recv_size_mpi
  PetscInt :: max_local_size
  PetscMPIInt :: local_size_mpi
  PetscInt :: istart, iend, num_in_array
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)

  
1000 format(100(i2,1x))
1001 format(100(i4,1x))
1002 format(100(i6,1x))
1003 format(100(i8,1x))
1004 format(100(i10,1x))
1010 format(100(es13.6,1x))

  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call PetscLogEventBegin(logging%event_output_write_tecplot,ierr)    

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
      call MPI_Allreduce(grid%nlmax,max_local_size,ONE_INTEGER_MPI, &
                         MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
      max_local_size_saved = max_local_size
      write(option%io_buffer,'("max_local_size_saved: ",i9)') max_local_size
      call printMsg(option)
    endif
    max_local_size = max_local_size_saved
    local_size_mpi = grid%nlmax
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

  call PetscLogEventEnd(logging%event_output_write_tecplot,ierr)    

end subroutine WriteTecplotDataSetNumPerLine1

! ************************************************************************** !
!
! OutputFormatInt: Writes a integer to a string
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFormatInt(int_value)

  implicit none
  
  PetscInt :: int_value
  
  character(len=MAXWORDLENGTH) :: OutputFormatInt

  write(OutputFormatInt,'(1i12)') int_value
  
  OutputFormatInt = adjustl(OutputFormatInt)
  
end function OutputFormatInt

! ************************************************************************** !
!
! OutputFormatDouble: Writes a double or real to a string
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputFormatDouble(real_value)

  implicit none
  
  PetscReal :: real_value
  
  character(len=MAXWORDLENGTH) :: OutputFormatDouble

  write(OutputFormatDouble,'(1es13.5)') real_value
  
  OutputFormatDouble = adjustl(OutputFormatDouble)
  
end function OutputFormatDouble

! ************************************************************************** !
!
! OutputObservationTecplot: Print to observation data to TECPLOT file
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine OutputObservationTecplot(realization)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  use Observation_module
  use Utility_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: fid, icell
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(observation_type), pointer :: observation
  PetscBool, save :: check_for_observation_points = PETSC_TRUE
  PetscBool, save :: open_file = PETSC_FALSE
  PetscInt :: local_id
  PetscInt :: icolumn

  call PetscLogEventBegin(logging%event_output_observation,ierr)    
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  if (check_for_observation_points) then
    open_file = PETSC_FALSE
    observation => patch%observation%first
    do
      if (.not.associated(observation)) exit
      if (observation%itype == OBSERVATION_SCALAR .or. &
          (observation%itype == OBSERVATION_FLUX .and. &
           option%myrank == option%io_rank)) then
        open_file = PETSC_TRUE
        exit
      endif
      observation => observation%next
    enddo
    check_for_observation_points = PETSC_FALSE
  endif
  
  
  if (open_file) then

    if (option%myrank < 10) then
      write(string,'(i1)') option%myrank  
    else if (option%myrank < 100) then
      write(string,'(i2)') option%myrank  
    else if (option%myrank < 1000) then
      write(string,'(i3)') option%myrank  
    else if (option%myrank < 10000) then
      write(string,'(i4)') option%myrank  
    else if (option%myrank < 100000) then
      write(string,'(i5)') option%myrank  
    endif
    filename = 'observation' // trim(option%group_prefix) // '-' // &
               trim(string) // '.tec'
  
    ! open file
    fid = 86
    if (observation_first .or. .not.FileExists(filename)) then
      open(unit=fid,file=filename,action="write",status="replace")
      ! write header
      ! write title
      write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
      observation => patch%observation%first

      ! must initialize icolumn here so that icolumn does not restart with
      ! each observation point
      if (output_option%print_column_ids) then
        icolumn = 1
      else
        icolumn = -1
      endif

      do 
        if (.not.associated(observation)) exit
        
        select case(observation%itype)
          case(OBSERVATION_SCALAR)
            if (associated(observation%region%coordinates) .and. &
                .not.observation%at_cell_center) then
              option%io_buffer = 'Writing of data at coordinates not ' // &
                'functioning properly for minerals.  Perhaps due to ' // &
                'non-ghosting of vol frac....>? - geh'
              call printErrMsg(option)
              call WriteObservationHeaderForCoord(fid,realization, &
                                                  observation%region, &
                                                  observation%print_velocities, &
                                                  observation% &
                                                  print_secondary_data, &
                                                  icolumn)
            else
              do icell=1,observation%region%num_cells
                call WriteObservationHeaderForCell(fid,realization, &
                                                   observation%region,icell, &
                                                   observation%print_velocities, &
                                                   observation% &
                                                   print_secondary_data, &
                                                   icolumn)
              enddo
            endif
          case(OBSERVATION_FLUX)
            if (option%myrank == option%io_rank) then
              call WriteObservationHeaderForBC(fid,realization, &
                                                observation%linkage_name)
            endif
        end select
        observation => observation%next
      enddo
      write(fid,'(a)',advance="yes") ""
    else
      open(unit=fid,file=filename,action="write",status="old", &
           position="append")
    endif
  
    observation => patch%observation%first
    write(fid,'(1es14.6)',advance="no") option%time/output_option%tconv
    do 
      if (.not.associated(observation)) exit
        select case(observation%itype)
          case(OBSERVATION_SCALAR)
            if (associated(observation%region%coordinates) .and. &
                .not.observation%at_cell_center) then
              call WriteObservationDataForCoord(fid,realization, &
                                                 observation%region)
              if (observation%print_velocities) then
                call WriteVelocityAtCoord(fid,realization, &
                                          observation%region)
              endif
            else
              do icell=1,observation%region%num_cells
                local_id = observation%region%cell_ids(icell)
                call WriteObservationDataForCell(fid,realization,local_id)
                if (observation%print_velocities) then
                  call WriteVelocityAtCell(fid,realization,local_id)
                endif
                if (observation%print_secondary_data(1)) then
                  call WriteObservationSecondaryDataAtCell(fid,realization, &
                                                           local_id, &
                                                           PRINT_SEC_TEMP)
                endif
                if (observation%print_secondary_data(2)) then
                  call WriteObservationSecondaryDataAtCell(fid,realization, &
                                                           local_id, &
                                                           PRINT_SEC_CONC)
                endif
                if (observation%print_secondary_data(3)) then
                  call WriteObservationSecondaryDataAtCell(fid,realization, &
                                                        local_id, &
                                                        PRINT_SEC_MIN_VOLFRAC)
                endif
              enddo
            endif
          case(OBSERVATION_FLUX)
            call WriteObservationDataForBC(fid,realization, &
                                            patch, &
                                            observation%connection_set)
      end select
      observation => observation%next
    enddo
    write(fid,'(a)',advance="yes") ""
    close(fid)

  endif

  observation_first = PETSC_FALSE
  
  call PetscLogEventEnd(logging%event_output_observation,ierr)    
      
end subroutine OutputObservationTecplot

! ************************************************************************** !
!
! WriteObservationHeaderForCell: Print a header for data at a cell
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine WriteObservationHeaderForCell(fid,realization,region,icell, &
                                         print_velocities, &
                                         print_secondary_data, &
                                         icolumn)

  use Realization_module
  use Grid_module
  use Option_module
  use Output_Aux_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(region_type) :: region
  PetscInt :: icell
  PetscBool :: print_velocities
  PetscBool :: print_secondary_data(3)
  PetscInt :: icolumn
  
  PetscInt :: local_id
  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  type(grid_type), pointer :: grid

  grid => realization%patch%grid
  
  local_id = region%cell_ids(icell)
  write(cell_string,*) grid%nG2A(grid%nL2G(region%cell_ids(icell)))
  cell_string = trim(region%name) // ' (' // trim(adjustl(cell_string)) // ')'

  ! add coordinate of cell center
  x_string = BestFloat(grid%x(grid%nL2G(local_id)),1.d4,1.d-2)
  y_string = BestFloat(grid%y(grid%nL2G(local_id)),1.d4,1.d-2)
  z_string = BestFloat(grid%z(grid%nL2G(local_id)),1.d4,1.d-2)
  cell_string = trim(cell_string) // ' (' // trim(adjustl(x_string)) // &
                ' ' // trim(adjustl(y_string)) // &
                ' ' // trim(adjustl(z_string)) // ')'
  
  call WriteObservationHeader(fid,realization,cell_string,print_velocities, &
                              print_secondary_data,icolumn)

end subroutine WriteObservationHeaderForCell

! ************************************************************************** !
!
! WriteObservationHeaderForCoord: Print a header for data at a coordinate
! author: Glenn Hammond
! date: 04/11/08
!
! ************************************************************************** !  
subroutine WriteObservationHeaderForCoord(fid,realization,region, &
                                         print_velocities, &
                                         print_secondary_data, &
                                         icolumn)

  use Realization_module
  use Option_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(region_type) :: region
  PetscBool :: print_velocities
  PetscBool :: print_secondary_data(3)
  PetscInt :: icolumn
  
  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  
  cell_string = trim(region%name)
  
  x_string = BestFloat(region%coordinates(ONE_INTEGER)%x,1.d4,1.d-2)
  y_string = BestFloat(region%coordinates(ONE_INTEGER)%y,1.d4,1.d-2)
  z_string = BestFloat(region%coordinates(ONE_INTEGER)%z,1.d4,1.d-2)
  cell_string = trim(cell_string) // ' (' // trim(adjustl(x_string)) // ' ' // &
                trim(adjustl(y_string)) // ' ' // &
                trim(adjustl(z_string)) // ')'

  call WriteObservationHeader(fid,realization,cell_string,print_velocities, &
                              print_secondary_data,icolumn)

end subroutine WriteObservationHeaderForCoord

! ************************************************************************** !
!
! WriteObservationHeader: Print a header for data
! author: Glenn Hammond
! date: 10/27/11
!
! ************************************************************************** !  
subroutine WriteObservationHeader(fid,realization,cell_string, &
                                         print_velocities, &
                                         print_secondary_data, &
                                         icolumn)
  use Realization_module
  use Option_module
  use Reactive_Transport_module
  use Secondary_Continuum_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscBool :: print_velocities
  PetscBool :: print_secondary_data(3)
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscInt :: icolumn
  
  PetscInt :: i
  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option  
  
  option => realization%option
  output_option => realization%output_option
  
  header = OutputVariableListToHeader(output_option%output_variable_list, &
                                      cell_string,icolumn,PETSC_FALSE)
  write(fid,'(a)',advance="no") trim(header)

  if (print_velocities) then
    header = ''
    write(string,'(''[m/'',a,'']'')') trim(realization%output_option%tunit)
    call OutputAppendToHeader(header,'vlx',string,cell_string,icolumn)
    call OutputAppendToHeader(header,'vly',string,cell_string,icolumn)
    call OutputAppendToHeader(header,'vlz',string,cell_string,icolumn)
    write(fid,'(a)',advance="no") trim(header)
  endif
  
  ! add secondary temperature to header
  if (print_secondary_data(1)) then
    select case (option%iflowmode) 
      case (THC_MODE, MPH_MODE) 
        header = ''
        do i = 1, option%nsec_cells
          write(string,'(i2)') i
          string = 'T_sec(' // trim(adjustl(string)) // ')'
          call OutputAppendToHeader(header,string,'',cell_string,icolumn)
        enddo
      case default
        header = ''
    end select
    write(fid,'(a)',advance="no") trim(header)
  endif
  
  ! add secondary concentrations to header
  if (print_secondary_data(2)) then
        header = ''
        do i = 1, option%nsec_cells
          write(string,'(i2)') i
          string = 'C_sec(' // trim(adjustl(string)) // ')'
          call OutputAppendToHeader(header,string,'',cell_string,icolumn)
        enddo
    write(fid,'(a)',advance="no") trim(header)
  endif
  
  ! add secondary mineral volume fractions to header
  if (print_secondary_data(3)) then
        header = ''
        do i = 1, option%nsec_cells
          write(string,'(i2)') i
          string = 'min_vol_frac_sec(' // trim(adjustl(string)) // ')'
          call OutputAppendToHeader(header,string,'',cell_string,icolumn)
        enddo
    write(fid,'(a)',advance="no") trim(header)
  endif  
  
end subroutine WriteObservationHeader

! ************************************************************************** !
!
! WriteObservationHeaderForBC: Print a header for data over a region
! author: Glenn Hammond
! date: 12/18/08
!
! ************************************************************************** !  
subroutine WriteObservationHeaderForBC(fid,realization,coupler_name)

  use Realization_module
  use Option_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: coupler_name
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction 
  
  option => realization%option
  reaction => realization%reaction
  
  select case(option%iflowmode)
    case(FLASH2_MODE)
    case(MPH_MODE)
    case(IMS_MODE)
    case(THC_MODE)
    case(THMC_MODE)
    case(MIS_MODE)
    case(RICHARDS_MODE)
      string = ',"Darcy flux ' // trim(coupler_name) // &
               ' [m^3/' // trim(realization%output_option%tunit) // ']"'
    case default
  end select
  write(fid,'(a)',advance="no") trim(string)

  if (associated(reaction)) then
    do i=1, reaction%naqcomp 
      ! may need to modify for molality vs molarity, but I believe molarity is correct
      write(fid,'(a)',advance="no") ',"' // &
        trim(reaction%primary_species_names(i)) // ' ' // &
        trim(coupler_name) // &
        ' [mol/' // trim(realization%output_option%tunit) // ']"'
    enddo
  endif

end subroutine WriteObservationHeaderForBC

! ************************************************************************** !
!
! WriteObservationDataForCell: Print data for data at a cell
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine WriteObservationDataForCell(fid,realization,local_id)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid, i
  type(realization_type) :: realization
  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option    
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  output_option => realization%output_option

100 format(es14.6)
101 format(i2)
110 format(es14.6)
111 format(i2)

  ghosted_id = grid%nL2G(local_id)
  ! write out coorindates
  !write(fid,110,advance="no") grid%x(ghosted_id)
  !write(fid,110,advance="no") grid%y(ghosted_id)
  !write(fid,110,advance="no") grid%z(ghosted_id)

  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,THMC_MODE,RICHARDS_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
    
      ! temperature
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,TEMPERATURE,ZERO_INTEGER,ghosted_id)
      end select

      ! liquid pressure
      write(fid,110,advance="no") &
        RealizGetDatasetValueAtCell(realization,LIQUID_PRESSURE,ZERO_INTEGER,ghosted_id)

      ! gas pressure
      select case(option%iflowmode)
        case(MPH_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,GAS_PRESSURE,ZERO_INTEGER,ghosted_id)
      end select

      ! state
      select case(option%iflowmode)
        case(G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,STATE,ZERO_INTEGER,ghosted_id)
      end select

      ! liquid saturation
      write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,LIQUID_SATURATION,ZERO_INTEGER,ghosted_id)

     ! gas saturation
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE,THC_MODE,THMC_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,GAS_SATURATION,ZERO_INTEGER,ghosted_id)
      end select

#ifdef ICE
     ! ice saturation
      select case(option%iflowmode)
        case(THC_MODE,THMC_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,ICE_SATURATION,ZERO_INTEGER,ghosted_id)
      end select

      ! ice density
      select case(option%iflowmode)
        case(THC_MODE,THMC_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,ICE_DENSITY,ZERO_INTEGER,ghosted_id)
      end select
#endif

      ! liquid density
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,LIQUID_DENSITY,ZERO_INTEGER,ghosted_id)
      end select

      ! gas density
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,GAS_DENSITY,ZERO_INTEGER,ghosted_id)
      end select

      ! liquid energy
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,LIQUID_ENERGY,ZERO_INTEGER,ghosted_id)
      end select

      ! gas energy
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,GAS_ENERGY,ZERO_INTEGER,ghosted_id)
      end select

      ! liquid viscosity
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,LIQUID_VISCOSITY,ZERO_INTEGER,ghosted_id)
      end select

      ! gas viscosity
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,GAS_VISCOSITY,ZERO_INTEGER,ghosted_id)
      end select

      ! liquid mobility
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,LIQUID_MOBILITY,ZERO_INTEGER,ghosted_id)
      end select

      ! gas mobility
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,GAS_MOBILITY,ZERO_INTEGER,ghosted_id)
      end select

      ! liquid mole fractions
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,FLASH2_MODE,MIS_MODE,G_MODE)
          do i=1,option%nflowspec
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,LIQUID_MOLE_FRACTION,i,ghosted_id)
          enddo
      end select

      ! gas mole fractions
      select case(option%iflowmode)
        case(MPH_MODE,FLASH2_MODE,G_MODE)
          do i=1,option%nflowspec
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,GAS_MOLE_FRACTION,i,ghosted_id)
          enddo
      end select 

      ! phase
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,111,advance="no") &
            int(RealizGetDatasetValueAtCell(realization,PHASE,ZERO_INTEGER,ghosted_id))
      end select
  end select

  if (option%ntrandof > 0) then
    reaction => realization%reaction
    ghosted_id = grid%nL2G(local_id)
    if (associated(reaction)) then
      if (reaction%print_pH .and. associated(reaction%species_idx)) then
        if (reaction%species_idx%h_ion_id > 0) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,PH,reaction%species_idx%h_ion_id,ghosted_id)
        endif
      endif
      if (reaction%print_total_component) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,reaction%print_tot_conc_type,i,ghosted_id)
          endif
        enddo
      endif
      if (reaction%print_free_ion) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,reaction%print_free_conc_type,i,ghosted_id)
          endif
        enddo
      endif      
      if (reaction%print_total_bulk) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,TOTAL_BULK,i,ghosted_id)
          endif
        enddo
      endif
      if (reaction%print_act_coefs) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,PRIMARY_ACTIVITY_COEF,i,ghosted_id)
          endif
        enddo
      endif
      do i=1,reaction%neqcplx
        if (reaction%secondary_species_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,reaction%print_secondary_conc_type,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,MINERAL_VOLUME_FRACTION,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
           write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,MINERAL_RATE,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%mineral%nmnrl
        if (reaction%mineral%mnrl_print(i)) then
           write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,MINERAL_SATURATION_INDEX,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_density_print(i)) then
           write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,SURFACE_SITE_DENSITY,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
           write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,SURFACE_CMPLX_FREE,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,SURFACE_CMPLX,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%surface_complexation%nkinsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
        !TODO(geh): fix
           write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,KIN_SURFACE_CMPLX_FREE,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%surface_complexation%nkinsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,KIN_SURFACE_CMPLX,i,ghosted_id)
        endif
      enddo
      if (associated(reaction%kd_print)) then
        do i=1,reaction%naqcomp
          if (reaction%kd_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,PRIMARY_KD,i,ghosted_id)
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_print)) then
        do i=1,reaction%naqcomp
          if (reaction%total_sorb_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,TOTAL_SORBED,i,ghosted_id)
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_mobile_print)) then
        do i=1,reaction%ncollcomp
          if (reaction%total_sorb_mobile_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,TOTAL_SORBED_MOBILE,i,ghosted_id)
          endif
        enddo
      endif
      if (reaction%print_colloid) then
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,COLLOID_MOBILE,i,ghosted_id)
          endif
        enddo
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,COLLOID_IMMOBILE,i,ghosted_id)
          endif
        enddo      
      endif
      
      if (reaction%print_age) then
        if (reaction%species_idx%tracer_age_id > 0) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,AGE, &
            reaction%species_idx%tracer_age_id,ghosted_id, &
            reaction%species_idx%tracer_aq_id)
        endif
      endif
    endif
  endif
          
  ! porosity
  if (output_option%print_porosity) then
    write(fid,110,advance="no") &
      RealizGetDatasetValueAtCell(realization,POROSITY,ZERO_INTEGER,ghosted_id)
  endif  

end subroutine WriteObservationDataForCell



! ************************************************************************** !
!
! WriteObservationDataForCoord: Print data for data at a coordinate
! author: Glenn Hammond
! date: 04/11/08
!
! ************************************************************************** !  
subroutine WriteObservationDataForCoord(fid,realization,region)

  use Realization_module
  use Option_module
  use Region_module  
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  
  use Structured_Grid_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(region_type) :: region

  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option  
    
  PetscInt :: ghosted_ids(8)
  PetscInt :: count
  PetscInt :: i, j, k
  PetscInt :: istart, iend, jstart, jend, kstart, kend
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  output_option => realization%output_option

100 format(es14.6)
101 format(i2)
110 format(es14.6)
111 format(i2)


  ! write out coorindates
!geh - do not print coordinates for now
  !write(fid,110,advance="no") region%coordinate(X_DIRECTION)
  !write(fid,110,advance="no") region%coordinate(Y_DIRECTION)
  !write(fid,110,advance="no") region%coordinate(Z_DIRECTION)
  
  count = 0
  local_id = region%cell_ids(1)
  ghosted_id = grid%nL2G(local_id)
  call StructGridGetIJKFromGhostedID(grid%structured_grid,ghosted_id,i,j,k)
  istart = i
  iend = i
  jstart = j
  jend = j
  kstart = k
  kend = k
  ! find the neighboring cells, between which to interpolate
  if (grid%x(ghosted_id) > region%coordinates(ONE_INTEGER)%x) then
    if (i > 1) then
      istart = i-1
    endif
  else
    if (i < grid%structured_grid%ngx) then
      iend = i+1
    endif
  endif
  if (grid%y(ghosted_id) > region%coordinates(ONE_INTEGER)%y) then
    if (j > 1) then
      jstart = j-1
    endif
  else
    if (j < grid%structured_grid%ngy) then
      jend = j+1
    endif
  endif
  if (grid%z(ghosted_id) > region%coordinates(ONE_INTEGER)%z) then
    if (k > 1) then
      kstart = k-1
    endif
  else
    if (k < grid%structured_grid%ngz) then
      kend = k+1
    endif
  endif
  count = 0
  do k=kstart,kend
    do j=jstart,jend
      do i=istart,iend
        count = count + 1
        ghosted_ids(count) = i + (j-1)*grid%structured_grid%ngx + &
                             (k-1)*grid%structured_grid%ngxy
      enddo
    enddo
  enddo

  select case(option%iflowmode)
    case(RICHARDS_MODE,MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)

    ! temperature
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,TEMPERATURE,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! liquid pressure
      write(fid,110,advance="no") &
        OutputGetVarFromArrayAtCoord(realization,LIQUID_PRESSURE,ZERO_INTEGER, &
                                     region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     count,ghosted_ids)

      ! gas pressure
      select case(option%iflowmode)
        case(MPH_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,GAS_PRESSURE,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! state
      select case(option%iflowmode)
        case(G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,STATE,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,RICHARDS_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,LIQUID_SATURATION,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! gas saturation
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE,THC_MODE,THMC_MODE)
          ! gas saturation
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,GAS_SATURATION,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

#ifdef ICE
    ! ice saturation
      select case(option%iflowmode)
        case(THC_MODE,THMC_MODE)
          ! ice saturation
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,ICE_SATURATION,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

    ! ice density
      select case(option%iflowmode)
        case(THC_MODE,THMC_MODE)
          ! ice saturation
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,ICE_DENSITY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select
#endif

      ! liquid density
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,LIQUID_DENSITY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! gas density
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,GAS_DENSITY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! liquid energy
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,LIQUID_ENERGY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! gas energy
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,GAS_ENERGY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! liquid viscosity
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,LIQUID_VISCOSITY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! gas viscosity
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,GAS_VISCOSITY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! liquid mobility
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,LIQUID_MOBILITY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! gas mobility
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,GAS_MOBILITY,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
      end select

      ! liquid mole fraction
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          do i=1,option%nflowspec
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,LIQUID_MOLE_FRACTION,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          enddo
      end select

     ! gas mole fractions
      select case(option%iflowmode)
        case(MPH_MODE,FLASH2_MODE,G_MODE)
          do i=1,option%nflowspec
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,GAS_MOLE_FRACTION,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          enddo
      end select 

      ! phase
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
         write(fid,111,advance="no") &
           int(OutputGetVarFromArrayAtCoord(realization,PHASE,ZERO_INTEGER, &
                                            region%coordinates(ONE_INTEGER)%x, &
                                            region%coordinates(ONE_INTEGER)%y, &
                                            region%coordinates(ONE_INTEGER)%z, &
                                            count,ghosted_ids))
      end select
  
  end select

  if (option%ntrandof > 0) then
    reaction => realization%reaction
    if (associated(reaction)) then
      if (reaction%print_pH .and. associated(reaction%species_idx)) then
        if (reaction%species_idx%h_ion_id > 0) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,PH,reaction%species_idx%h_ion_id, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      endif
      if (reaction%print_total_component) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,reaction%print_tot_conc_type,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
      endif
      if (reaction%print_free_ion) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,reaction%print_free_conc_type,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
      endif      
      if (reaction%print_total_bulk) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,TOTAL_BULK,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
      endif
      if (reaction%print_act_coefs) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,PRIMARY_ACTIVITY_COEF,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
          endif
        enddo
      endif
      do i=1,reaction%neqcplx
        if (reaction%secondary_species_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,reaction%print_secondary_conc_type,i, &
                                          region%coordinates(ONE_INTEGER)%x, &
                                          region%coordinates(ONE_INTEGER)%y, &
                                          region%coordinates(ONE_INTEGER)%z, &
                                          count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,MINERAL_VOLUME_FRACTION,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,MINERAL_RATE,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%mineral%nmnrl
        if (reaction%mineral%mnrl_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,MINERAL_SATURATION_INDEX,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_density_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,SURFACE_SITE_DENSITY,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,SURFACE_CMPLX_FREE,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,SURFACE_CMPLX,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
        !TODO(geh): fix
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,KIN_SURFACE_CMPLX_FREE,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,KIN_SURFACE_CMPLX,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      if (associated(reaction%kd_print)) then
        do i=1,reaction%naqcomp
          if (reaction%kd_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,PRIMARY_KD,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_print)) then
        do i=1,reaction%naqcomp
          if (reaction%total_sorb_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,TOTAL_SORBED,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_mobile_print)) then
        do i=1,reaction%ncollcomp
          if (reaction%total_sorb_mobile_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,TOTAL_SORBED_MOBILE,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
      endif
      if (reaction%print_colloid) then
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,COLLOID_MOBILE,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCoord(realization,COLLOID_IMMOBILE,i, &
                                           region%coordinates(ONE_INTEGER)%x, &
                                           region%coordinates(ONE_INTEGER)%y, &
                                           region%coordinates(ONE_INTEGER)%z, &
                                           count,ghosted_ids)
          endif
        enddo
      endif
      if (reaction%print_age) then
        if (reaction%species_idx%tracer_age_id > 0) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,AGE, &
              reaction%species_idx%tracer_age_id, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids, &
                                         reaction%species_idx%tracer_aq_id)
        endif
      endif
    endif
  endif
  
  ! porosity
  if (output_option%print_porosity) then
    write(fid,110,advance="no") &
      OutputGetVarFromArrayAtCoord(realization,POROSITY,ZERO_INTEGER, &
                                   region%coordinates(ONE_INTEGER)%x, &
                                   region%coordinates(ONE_INTEGER)%y, &
                                   region%coordinates(ONE_INTEGER)%z, &
                                   count,ghosted_ids)
  endif

end subroutine WriteObservationDataForCoord

! ************************************************************************** !
!
! WriteObservationDataForBC: Print flux data for a boundary condition
! author: Glenn Hammond
! date: 12/18/08
!
! ************************************************************************** !  
subroutine WriteObservationDataForBC(fid,realization,patch,connection_set)

  use Realization_module
  use Option_module
  use Connection_module  
  use Patch_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(patch_type), pointer :: patch
  type(connection_set_type), pointer :: connection_set

  PetscInt :: i
  PetscInt :: iconn
  PetscInt :: offset
  PetscInt :: iphase
  PetscMPIInt :: int_mpi
  PetscReal :: sum_volumetric_flux(realization%option%nphase)
  PetscReal :: sum_volumetric_flux_global(realization%option%nphase)
  PetscReal :: sum_solute_flux(realization%option%ntrandof)
  PetscReal :: sum_solute_flux_global(realization%option%ntrandof)
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  
  option => realization%option
  reaction => realization%reaction

100 format(es14.6)
!100 format(es16.9)
101 format(i2)
110 format(es14.6)
!110 format(',',es16.9)
111 format(i2)
 
  iphase = 1

  ! sum up fluxes across region
  if (associated(connection_set)) then
    offset = connection_set%offset
    select case(option%iflowmode)
      case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
      case(MIS_MODE)
      case(RICHARDS_MODE)
        sum_volumetric_flux = 0.d0
        if (associated(connection_set)) then
          do iconn = 1, connection_set%num_connections
            sum_volumetric_flux(:) = sum_volumetric_flux(:) + &
                                  patch%boundary_velocities(iphase,offset+iconn)* &
                                  connection_set%area(iconn)
          enddo
        endif
        int_mpi = option%nphase
        call MPI_Reduce(sum_volumetric_flux,sum_volumetric_flux_global, &
                        int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        option%io_rank,option%mycomm,ierr)
        if (option%myrank == option%io_rank) then
          do i = 1, option%nphase
            write(fid,110,advance="no") sum_volumetric_flux_global(i)
          enddo
        endif
    end select

    if (associated(reaction)) then
      sum_solute_flux = 0.d0
      if (associated(connection_set)) then
        do iconn = 1, connection_set%num_connections
          sum_solute_flux(:) = sum_solute_flux(:) + &
                               patch%boundary_fluxes(iphase,:,offset+iconn)* &
                               connection_set%area(iconn)
        enddo
      endif
      int_mpi = option%ntrandof
      call MPI_Reduce(sum_solute_flux,sum_solute_flux_global, &
                      int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
      if (option%myrank == option%io_rank) then
        !we currently only print the aqueous components
        do i = 1, reaction%naqcomp
          write(fid,110,advance="no") sum_solute_flux_global(i)
        enddo
      endif
    endif

  endif

end subroutine WriteObservationDataForBC

! ************************************************************************** !
!
! WriteVelocityAtCell: Computes velocities at a grid cell
! author: Glenn Hammond
! note: limited to structured grids
! date: 03/20/08
!
! ************************************************************************** !  
subroutine WriteVelocityAtCell(fid,realization,local_id)

  use Realization_module
  use Option_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscInt :: local_id

  PetscReal :: velocity(1:3)
  
200 format(3(es14.6))
  
  velocity = GetVelocityAtCell(fid,realization,local_id)
  
  write(fid,200,advance="no") velocity(1:3)*realization%output_option%tconv   

end subroutine WriteVelocityAtCell

! ************************************************************************** !
!
! GetVelocityAtCell: Computes velocities at a grid cell
! author: Glenn Hammond
! note: limited to structured grids
! date: 03/20/08
!
! ************************************************************************** !  
function GetVelocityAtCell(fid,realization,local_id)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Connection_module
  use Coupler_module

  implicit none
  
  PetscReal :: GetVelocityAtCell(3)
  PetscInt :: fid
  type(realization_type) :: realization
  PetscInt :: local_id

  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(connection_set_list_type), pointer :: connection_set_list
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: direction, iphase
  PetscReal :: area
  PetscReal :: sum_velocity(1:3), sum_area(1:3), velocity(1:3)
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  sum_velocity = 0.d0
  sum_area = 0.d0
  iphase = 1

  ! interior velocities  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_up = grid%nG2L(cur_connection_set%id_up(iconn)) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(cur_connection_set%id_dn(iconn)) ! = zero for ghost nodes
      if (local_id_up == local_id .or. local_id_dn == local_id) then
        do direction=1,3        
          area = cur_connection_set%area(iconn)* &
                 !geh: no dabs() here
                 cur_connection_set%dist(direction,iconn)
          sum_velocity(direction) = sum_velocity(direction) + &
                                    patch%internal_velocities(iphase,sum_connection)* &
                                    area
          sum_area(direction) = sum_area(direction) + dabs(area)
        enddo
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
      if (cur_connection_set%id_dn(iconn) == local_id) then
        do direction=1,3        
          area = cur_connection_set%area(iconn)* &
                 !geh: no dabs() here
                 cur_connection_set%dist(direction,iconn)
          sum_velocity(direction) = sum_velocity(direction) + &
                                    patch%boundary_velocities(iphase,sum_connection)* &
                                    area
          sum_area(direction) = sum_area(direction) + dabs(area)
        enddo
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  velocity = 0.d0
  do direction = 1,3
    if (abs(sum_area(direction)) > 1.d-40) &
      velocity(direction) = sum_velocity(direction)/sum_area(direction)
  enddo

  GetVelocityAtCell = velocity  

end function GetVelocityAtCell

! ************************************************************************** !
!
! WriteVelocityAtCoord: Computes velocities at a coordinate
! author: Glenn Hammond
! note: limited to structured grids
! date: 03/20/08
!
! ************************************************************************** !  
subroutine WriteVelocityAtCoord(fid,realization,region)

  use Realization_module
  use Region_module
  use Option_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(region_type) :: region
  PetscInt :: local_id
  PetscReal :: coordinate(3)

  PetscReal :: velocity(1:3)
  
200 format(3(es14.6))
  
  velocity = GetVelocityAtCoord(fid,realization,region%cell_ids(1), &
                                region%coordinates(ONE_INTEGER)%x, &
                                region%coordinates(ONE_INTEGER)%y, &
                                region%coordinates(ONE_INTEGER)%z)
  
  write(fid,200,advance="no") velocity(1:3)*realization%output_option%tconv   

end subroutine WriteVelocityAtCoord

! ************************************************************************** !
!
! GetVelocityAtCoord: Computes velocities at a coordinate
! author: Glenn Hammond
! note: limited to structured grids
! date: 03/20/08
!
! ************************************************************************** !  
function GetVelocityAtCoord(fid,realization,local_id,x,y,z)
  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Connection_module
  use Coupler_module

  implicit none
  
  PetscReal :: GetVelocityAtCoord(3)
  PetscInt :: fid
  type(realization_type) :: realization
  PetscInt :: local_id
  PetscReal :: x, y, z
  
  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(connection_set_list_type), pointer :: connection_set_list
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn
  PetscReal :: cell_coord(3), face_coord
  PetscReal :: coordinate(3)
  PetscInt :: direction, iphase
  PetscReal :: area, weight, distance
  PetscReal :: sum_velocity(1:3), velocity(1:3)
  PetscReal :: sum_weight(1:3)
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  sum_velocity = 0.d0
  sum_weight = 0.d0
  iphase = 1

  ghosted_id = grid%nL2G(local_id)
  
  coordinate(X_DIRECTION) = x
  coordinate(Y_DIRECTION) = y
  coordinate(Z_DIRECTION) = z

  cell_coord(X_DIRECTION) = grid%x(ghosted_id)
  cell_coord(Y_DIRECTION) = grid%y(ghosted_id)
  cell_coord(Z_DIRECTION) = grid%z(ghosted_id)

  ! interior velocities  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_up = grid%nG2L(cur_connection_set%id_up(iconn)) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(cur_connection_set%id_dn(iconn)) ! = zero for ghost nodes
      if (local_id_up == local_id .or. local_id_dn == local_id) then
        do direction=1,3
          if (local_id_up == local_id) then
            face_coord = cell_coord(direction) + &
                         cur_connection_set%dist(-1,iconn)* &
                         cur_connection_set%dist(0,iconn)* &
                         cur_connection_set%dist(direction,iconn)
          else
            face_coord = cell_coord(direction) - &
                         (1.d0-cur_connection_set%dist(-1,iconn))* &
                         cur_connection_set%dist(0,iconn)* &
                         cur_connection_set%dist(direction,iconn)
          endif
          distance = dabs(face_coord-coordinate(direction))
          if (distance < 1.d-40) distance = 1.d-40
          weight = cur_connection_set%area(iconn)* &
                 dabs(cur_connection_set%dist(direction,iconn))/ &
                 distance
 
          sum_velocity(direction) = sum_velocity(direction) + &
                                    cur_connection_set%dist(direction,iconn)* &
                                    patch%internal_velocities(iphase,sum_connection)* &
                                    weight
          sum_weight(direction) = sum_weight(direction) + weight
       enddo
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
      if (cur_connection_set%id_dn(iconn) == local_id) then
        do direction=1,3        
          face_coord = cell_coord(direction) - &
                    !   (1.d0-cur_connection_set%dist(-1,iconn))* & ! fraction upwind is always 0.d0
                       cur_connection_set%dist(0,iconn)* &
                       cur_connection_set%dist(direction,iconn)
          distance = dabs(face_coord-coordinate(direction))
          if (distance < 1.d-40) distance = 1.d-40
          weight = cur_connection_set%area(iconn)* &
                   dabs(cur_connection_set%dist(direction,iconn))/ &
                   distance
          sum_velocity(direction) = sum_velocity(direction) + &
                                    cur_connection_set%dist(direction,iconn)* &
                                    patch%boundary_velocities(iphase,sum_connection)* &
                                    weight
          sum_weight(direction) = sum_weight(direction) + weight
        enddo
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  velocity = 0.d0
  do direction = 1,3
    if (abs(sum_weight(direction)) > 1.d-40) &
      velocity(direction) = sum_velocity(direction)/sum_weight(direction)
  enddo

  GetVelocityAtCoord = velocity  

end function GetVelocityAtCoord


! ************************************************************************** !
!
! WriteObservationSecondaryDataAtCell: Print data for data at a cell
! author: Satish Karra
! date: 10/4/12
!
! ************************************************************************** !  
subroutine WriteObservationSecondaryDataAtCell(fid,realization,local_id,ivar)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module

  implicit none
  
  PetscInt :: fid, i
  type(realization_type) :: realization
  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option    
  PetscInt :: ivar
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  output_option => realization%output_option

100 format(es14.6)
101 format(i2)
110 format(es14.6)
111 format(i2)

  ghosted_id = grid%nL2G(local_id)

  if (option%nsec_cells > 0) then
    if (ivar == PRINT_SEC_TEMP) then
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE)
          do i = 1, option%nsec_cells 
            write(fid,110,advance="no") &
              RealizGetDatasetValueAtCell(realization,SECONDARY_TEMPERATURE,i, &
                                          ghosted_id)
          enddo
        end select
     endif
     if (ivar == PRINT_SEC_CONC) then
       do i = 1, option%nsec_cells 
         write(fid,110,advance="no") &
           RealizGetDatasetValueAtCell(realization,SECONDARY_CONCENTRATION,i, &
                                       ghosted_id)
       enddo
     endif
     if (ivar == PRINT_SEC_MIN_VOLFRAC) then
       do i = 1, option%nsec_cells 
         write(fid,110,advance="no") &
           RealizGetDatasetValueAtCell(realization,SEC_MIN_VOLFRAC,i, &
                                       ghosted_id)
       enddo
     endif
   endif 
   
  
end subroutine WriteObservationSecondaryDataAtCell

! ************************************************************************** !
!
! OutputVTK: Print to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !  
subroutine OutputVTK(realization)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Structured_Grid_module
  use Option_module
  use Field_module
  use Patch_module
  
  use Flash2_module
  use Mphase_module
  use Immis_module
  use Miscible_module
  use THC_module
  use THMC_module
  use Richards_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch 
  type(reaction_type), pointer :: reaction 
  type(output_option_type), pointer :: output_option
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec
  Vec :: natural_vec
  
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option
  
  if (reaction%print_free_conc_type == PRIMARY_MOLALITY) then
    free_mol_char = 'm'
  else
    free_mol_char = 'M'
  endif
  
  if (reaction%print_tot_conc_type == TOTAL_MOLALITY) then
    tot_mol_char = 'm'
  else
    tot_mol_char = 'M'
  endif
  
  if (reaction%print_secondary_conc_type == SECONDARY_MOLALITY) then
    sec_mol_char = 'm'
  else
    sec_mol_char = 'M'
  endif
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.vtk'
  else
    string = OutputFilenameID(output_option,option)
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-' // trim(string) // '.vtk'    
  endif
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write vtk output file: ' // trim(filename)
    call printMsg(option)    
    open(unit=OUTPUT_UNIT,file=filename,action="write")
  
    ! write header
    write(OUTPUT_UNIT,'(''# vtk DataFile Version 2.0'')')
    ! write title
    write(OUTPUT_UNIT,'(''PFLOTRAN output'')')
    write(OUTPUT_UNIT,'(''ASCII'')')
    write(OUTPUT_UNIT,'(''DATASET UNSTRUCTURED_GRID'')')
  endif

  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)  

  ! write out coordinates
  call WriteVTKGrid(OUTPUT_UNIT,realization)

  write(OUTPUT_UNIT,'(''CELL_DATA'',i8)') grid%nmax

  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,THMC_MODE,RICHARDS_MODE,IMS_MODE,MIS_MODE, &
         FLASH2_MODE,G_MODE)

      ! temperature
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          word = 'Temperature'
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
      end select

      ! pressure
      word = 'Liquid Pressure'
      call OutputGetVarFromArray(realization,global_vec,LIQUID_PRESSURE,ZERO_INTEGER)
      call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
      call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)

      ! phase
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          word = 'Phase'
          call OutputGetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_INTEGER)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,RICHARDS_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          word = 'Liquid_Saturation'
          call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
      end select

      ! gas saturation
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE,THC_MODE,THMC_MODE)
          word = 'Gas_Saturation'
          call OutputGetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
      end select

#ifdef ICE
      ! ice saturation
      select case(option%iflowmode)
        case(THC_MODE,THMC_MODE)
          word = 'Ice_Saturation'
          call OutputGetVarFromArray(realization,global_vec,ICE_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
      end select

      ! ice density
      select case(option%iflowmode)
        case(THC_MODE,THMC_MODE)
          word = 'Ice_Density'
          call OutputGetVarFromArray(realization,global_vec,ICE_DENSITY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
      end select
#endif
    
      ! liquid energy
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          word = 'Liquid_Energy'
          call OutputGetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
      end select
    
     ! gas energy
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          word = 'Gas_Energy'
          call OutputGetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
      end select

      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          ! liquid mole fractions
          do i=1,option%nflowspec
            write(word,'(i2)') i
            word = 'Xl(' // trim(adjustl(word)) // ')'
            call OutputGetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          enddo
      end select
  
      select case(option%iflowmode)
        case(MPH_MODE,FLASH2_MODE,G_MODE)
          ! gas mole fractions
          do i=1,option%nflowspec
            write(word,'(i2)') i
            word = 'Xg(' // trim(adjustl(word)) // ')'
            call OutputGetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          enddo
      end select 
      
    case default
  
  end select
  
  if (option%ntrandof > 0) then
    if (associated(reaction)) then
    
      if (reaction%print_pH .and. associated(reaction%species_idx)) then
        if (reaction%species_idx%h_ion_id > 0) then
          call OutputGetVarFromArray(realization,global_vec,PH,reaction%species_idx%h_ion_id)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = 'pH'
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      endif
    
    
      if (reaction%print_total_component) then
        do i=1,reaction%naqcomp
          call OutputGetVarFromArray(realization,global_vec,reaction%print_tot_conc_type,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%primary_species_names(i)) // '_tot_' // trim(tot_mol_char)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word, &
                                      natural_vec,VTK_REAL)
        enddo
      endif
      if (reaction%print_free_ion) then
        do i=1,reaction%naqcomp
          call OutputGetVarFromArray(realization,global_vec,reaction%print_free_conc_type,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%primary_species_names(i)) // '_free_' // trim(free_mol_char)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word, &
                                      natural_vec,VTK_REAL)
        enddo
      endif    
      if (reaction%print_total_bulk) then
        do i=1,reaction%naqcomp
          call OutputGetVarFromArray(realization,global_vec,TOTAL_BULK,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%primary_species_names(i)) // '_total_bulk'
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word, &
                                      natural_vec,VTK_REAL)
        enddo
      endif
      if (reaction%print_act_coefs) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,PRIMARY_ACTIVITY_COEF,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            word = trim(reaction%primary_species_names(i)) // '_gam'
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          endif
        enddo
      endif
      do i=1,reaction%neqcplx
        if (reaction%secondary_species_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,reaction%print_secondary_conc_type,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%secondary_species_names(i)) // &
               '_' // trim(sec_mol_char)
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_VOLUME_FRACTION,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%mineral%kinmnrl_names(i)) // '_vf'
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word, &
                                    natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_RATE,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%mineral%kinmnrl_names(i)) // '_rt'
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word, &
                                    natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%mineral%nmnrl
        if (reaction%mineral%mnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_SATURATION_INDEX,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%mineral%kinmnrl_names(i)) // '_si'
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_density_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_SITE_DENSITY,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%surface_complexation%srfcplxrxn_site_names(i))
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX_FREE,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%surface_complexation%srfcplxrxn_site_names(i))
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%surface_complexation%srfcplx_names(i))
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%surface_complexation%nkinsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
        !TODO(geh): fix
          call OutputGetVarFromArray(realization,global_vec,KIN_SURFACE_CMPLX_FREE,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%surface_complexation%srfcplxrxn_site_names(i))
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      enddo
      do i=1,reaction%surface_complexation%nkinsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,KIN_SURFACE_CMPLX,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = trim(reaction%surface_complexation%srfcplx_names(i))
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      enddo
      if (associated(reaction%kd_print)) then
        do i=1,reaction%naqcomp
          if (reaction%kd_print(i)) then      
            call OutputGetVarFromArray(realization,global_vec,PRIMARY_KD,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            word = trim(reaction%primary_species_names(i)) // '_kd'
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_print)) then
        do i=1,reaction%naqcomp
          if (reaction%total_sorb_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,TOTAL_SORBED,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            word = trim(reaction%primary_species_names(i)) // '_total_sorb'
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_mobile_print)) then
        do i=1,reaction%ncollcomp
          if (reaction%total_sorb_mobile_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,TOTAL_SORBED_MOBILE,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            word = trim(reaction%colloid_species_names(i)) // '_total_sorb_mob'
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          endif
        enddo
      endif
      if (reaction%print_colloid) then
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,COLLOID_MOBILE,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            word = trim(reaction%colloid_names(i)) // '_col_mob_' // &
                 trim(tot_mol_char)
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          endif
        enddo
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,COLLOID_IMMOBILE,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            word = trim(reaction%colloid_names(i)) // '_col_imb_' // &
                 trim(tot_mol_char)
            call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
          endif
        enddo
      endif
      if (reaction%print_age) then
        if (reaction%species_idx%tracer_age_id > 0) then
          call OutputGetVarFromArray(realization,global_vec,AGE, &
            reaction%species_idx%tracer_age_id,reaction%species_idx%tracer_aq_id)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          word = 'Tracer_Age'
          call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
        endif
      endif
    endif
  endif
  
  ! porosity
  if (output_option%print_porosity) then
    word = 'Porosity'
    call OutputGetVarFromArray(realization,global_vec,POROSITY,ZERO_INTEGER)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_INTEGER)
  endif

  ! material id
  word = 'Material_ID'
  call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_INTEGER)

  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)
  
  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)

#if 1
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesVTK(realization)
  endif
#endif
  
#if 0  
  if (output_option%print_tecplot_flux_velocities) then
    if (grid%structured_grid%nx > 1) then
      call OutputFluxVelocitiesVTK(realization,LIQUID_PHASE, &
                                          X_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesVTK(realization,GAS_PHASE, &
                                              X_DIRECTION)
      end select
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesVTK(realization,LIQUID_PHASE, &
                                          Y_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesVTK(realization,GAS_PHASE, &
                                              Y_DIRECTION)
      end select
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesVTK(realization,LIQUID_PHASE, &
                                          Z_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesVTK(realization,GAS_PHASE, &
                                              Z_DIRECTION)
      end select
    endif
  endif
#endif
      
end subroutine OutputVTK

#if 1
! ************************************************************************** !
!
! OutputVelocitiesVTK: Print velocities to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVelocitiesVTK(realization)
 
  use Realization_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  Vec :: global_vec
  Vec :: natural_vec

  PetscReal, pointer :: vec_ptr(:)
  
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  option => realization%option
  output_option => realization%output_option
  discretization => realization%discretization
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-vel.vtk'
  else  
    string = OutputFilenameID(output_option,option)
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-vel-' // trim(string) // '.vtk'
  endif
  
  if (option%myrank == option%io_rank) then
   option%io_buffer = '--> write vtk velocity output file: ' // &
                      trim(filename)
    call printMsg(option)                      
    open(unit=OUTPUT_UNIT,file=filename,action="write")
  
    ! write header
    write(OUTPUT_UNIT,'(''# vtk DataFile Version 2.0'')')
    ! write title
    write(OUTPUT_UNIT,'(''PFLOTRAN output'')')
    write(OUTPUT_UNIT,'(''ASCII'')')
    write(OUTPUT_UNIT,'(''DATASET UNSTRUCTURED_GRID'')')
    
  endif
  
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)    

  ! write out coordinates
  call WriteVTKGrid(OUTPUT_UNIT,realization)

  word = 'Vlx'
  call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)

  word = 'Vly'
  call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)

  word = 'Vlz'
  call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)

  if (option%nphase > 1) then
    word = 'Vgx'
    call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)

    word = 'Vgy'
    call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)

    word = 'Vgz'
    call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_REAL)
  endif

  ! material id
  word = 'Material_ID'
  call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteVTKDataSetFromVec(OUTPUT_UNIT,realization,word,natural_vec,VTK_INTEGER)
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)
  
end subroutine OutputVelocitiesVTK
#endif
! ************************************************************************** !
!
! WriteVTKGrid: Writes a grid in VTK format
! author: Glenn Hammond
! date: 11/05/08
!
! ************************************************************************** !
subroutine WriteVTKGrid(fid,realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscInt :: i, j, k, nx, ny, nz
  PetscReal :: x, y, z
  PetscInt :: nxp1Xnyp1, nxp1, nyp1, nzp1
  PetscInt :: vertex_id

1000 format(es13.6,1x,es13.6,1x,es13.6)
1001 format(i1,8(1x,i8))
  
  call PetscLogEventBegin(logging%event_output_grid_vtk,ierr) 
                              
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  if ((realization%discretization%itype == STRUCTURED_GRID).or. &
        (realization%discretization%itype == STRUCTURED_GRID_MIMETIC))  then

    nx = grid%structured_grid%nx
    ny = grid%structured_grid%ny
    nz = grid%structured_grid%nz
  
    nxp1 = nx+1
    nyp1 = ny+1
    nzp1 = nz+1
  
    if (option%myrank == option%io_rank) then

 1010 format("POINTS",1x,i12,1x,"float")
      write(fid,1010) (nx+1)*(ny+1)*(nz+1)
      do k=0,nz
        if (k > 0) then
          z = z + grid%structured_grid%dz_global(k)
        else
          z = grid%structured_grid%origin(Z_DIRECTION)
        endif
        do j=0,ny
          if (j > 0) then
            y = y + grid%structured_grid%dy_global(j)
          else
            y = grid%structured_grid%origin(Y_DIRECTION)
          endif
          x = grid%structured_grid%origin(X_DIRECTION)
          write(fid,1000) x,y,z
          do i=1,nx
            x = x + grid%structured_grid%dx_global(i)
            write(fid,1000) x,y,z
          enddo
        enddo
      enddo

1020 format('CELLS',1x,i12,1x,i12)
      write(fid,1020) grid%nmax, grid%nmax*9
      nxp1Xnyp1 = nxp1*nyp1
      do k=0,nz-1
        do j=0,ny-1
          do i=0,nx-1
            vertex_id = i+j*nxp1+k*nxp1Xnyp1
            write(fid,1001) 8,vertex_id,vertex_id+1, &
                            vertex_id+nxp1+1,vertex_id+nxp1, &
                            vertex_id+nxp1Xnyp1,vertex_id+nxp1Xnyp1+1, &
                            vertex_id+nxp1Xnyp1+nxp1+1, &
                            vertex_id+nxp1Xnyp1+nxp1
          enddo
        enddo
      enddo

      write(fid,'(a)') ""

1030 format('CELL_TYPES',1x,i12)
      write(fid,1030) grid%nmax
      do i=1,grid%nmax
        write(fid,'(i2)') 12
      enddo

      write(fid,'(a)') ""

    endif
  endif

  call PetscLogEventEnd(logging%event_output_grid_vtk,ierr) 
                            
end subroutine WriteVTKGrid

! ************************************************************************** !
!
! WriteVTKDataSetFromVec: Writes data from a Petsc Vec within a block
!                             of a VTK file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteVTKDataSetFromVec(fid,realization,dataset_name,vec,datatype)

  use Realization_module
  
  implicit none

  PetscInt :: fid
  type(realization_type) :: realization
  Vec :: vec
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: datatype
  
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteVTKDataSet(fid,realization,dataset_name,vec_ptr,datatype, &
                       ZERO_INTEGER) ! 0 implies grid%nlmax
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteVTKDataSetFromVec

! ************************************************************************** !
!
! WriteVTKDataSet: Writes data from an array within a block
!                      of a VTK file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteVTKDataSet(fid,realization,dataset_name,array,datatype, &
                           size_flag)

  use Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscReal :: array(:)
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscInt :: i
  PetscInt :: max_proc, max_proc_prefetch
  PetscMPIInt :: iproc_mpi, recv_size_mpi
  PetscInt :: max_local_size
  PetscMPIInt :: local_size_mpi
  PetscInt :: istart, iend, num_in_array
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)

1001 format(10(es13.6,1x))
1002 format(i3)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call PetscLogEventBegin(logging%event_output_write_vtk,ierr)    

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
      call MPI_Allreduce(grid%nlmax,max_local_size,ONE_INTEGER_MPI, &
                         MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
      max_local_size_saved = max_local_size
      if (OptionPrintToScreen(option)) print *, 'max_local_size_saved: ', &
                                                 max_local_size
    endif
    max_local_size = max_local_size_saved
    local_size_mpi = grid%nlmax
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == VTK_INTEGER) then
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

!    write(fid,'(''CELL_DATA'',i8)') grid%nmax

    if (datatype == VTK_INTEGER) then
      write(fid,'(''SCALARS '',a20,'' int 1'')') dataset_name
    else
      write(fid,'(''SCALARS '',a20,'' float 1'')') dataset_name
    endif
    
    write(fid,'(''LOOKUP_TABLE default'')') 

    if (datatype == VTK_INTEGER) then
      ! This approach makes output files identical, regardless of processor
      ! distribution.  It is necessary when diffing files.
      iend = 0
      do
        istart = iend+1
        if (iend+10 > local_size_mpi) exit
        iend = istart+9
        write(fid,1002) integer_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      integer_data(1:local_size_mpi-iend) = integer_data(iend+1:local_size_mpi)
      num_in_array = local_size_mpi-iend
    else
      iend = 0
      do
        istart = iend+1
        if (iend+10 > local_size_mpi) exit
        iend = istart+9
        write(fid,1001) real_data(istart:iend)
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
      if (datatype == 0) then
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
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(fid,1002) integer_data(istart:iend)
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
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(fid,1001) real_data(istart:iend)
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
    if (datatype == 0) then
      if (num_in_array > 0) &
        write(fid,1002) integer_data(1:num_in_array)
    else
      if (num_in_array > 0) &
        write(fid,1001) real_data(1:num_in_array)
    endif
    write(fid,'(/)')
  else
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        if (option%myrank < max_proc) exit
        call MPI_Bcast(max_proc,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                       option%mycomm,ierr)
      enddo
    endif
#endif    
    if (datatype == VTK_INTEGER) then
      call MPI_Send(integer_data,local_size_mpi,MPIU_INTEGER,option%io_rank, &
                    local_size_mpi,option%mycomm,ierr)
    else
      call MPI_Send(real_data,local_size_mpi,MPI_DOUBLE_PRECISION,option%io_rank, &
                    local_size_mpi,option%mycomm,ierr)
    endif
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        call MPI_Bcast(max_proc,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                       option%mycomm,ierr)
        if (max_proc < 0) exit
      enddo
    endif
#endif
#undef HANDSHAKE
  endif
      
  if (datatype == VTK_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

  call PetscLogEventEnd(logging%event_output_write_vtk,ierr)    

end subroutine WriteVTKDataSet

! ************************************************************************** !
!
! OutputHDF5: Print to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputHDF5(realization)

  use Realization_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  
#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  
  type(realization_type) :: realization

  call printMsg(realization%option,'')
  write(realization%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to &
        &write HDF5 formatted structured grids Darn.")')
  call printErrMsg(realization%option)
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

  type(realization_type) :: realization

#if defined(PARALLELIO_LIB_WRITE)
  integer:: file_id
  integer:: grp_id
  integer:: file_space_id
  integer:: realization_set_id
  integer:: prop_id
  PetscMPIInt :: rank
  integer:: dims(3)
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
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  
  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type

  discretization => realization%discretization
  patch => realization%patch
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option

  if (realization%discretization%itype == UNSTRUCTURED_GRID) then
    call OutputHDF5UGrid(realization)
    return
  endif
  
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
#if defined(PARALLELIO_LIB_WRITE)
  if (.not.first) then
    filename = trim(filename) // CHAR(0)
    call parallelio_open_file(filename, option%iowrite_group_id, &
                              FILE_READWRITE, file_id, ierr)
    if (file_id == -1) first = PETSC_TRUE
  endif
  if (first) then
    filename = trim(filename) // CHAR(0)
    call parallelio_open_file(filename, option%iowrite_group_id, &
                              FILE_CREATE, file_id, ierr)
  endif

#else  ! PARALLELIO_LIB_WRITE is not defined

    ! initialize fortran interface
  call h5open_f(hdf5_err)

#ifdef VAMSI_HDF5_WRITE
  if (mod(option%myrank,option%hdf5_write_group_size) == 0) then 
#endif

    call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
#ifdef VAMSI_HDF5_WRITE
    call h5pset_fapl_mpio_f(prop_id,option%writers,MPI_INFO_NULL,hdf5_err) 
#else
    call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
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
#endif ! PARALLELIO_LIB_WRITE

    if (first) then
      option%io_buffer = '--> creating hdf5 output file: ' // filename
    else
      option%io_buffer = '--> appending to hdf5 output file: ' // filename
    endif
    call printMsg(option)

    if (first) then

      ! create a group for the coordinates data set
#if defined(PARALLELIO_LIB_WRITE)
      string = "Coordinates" // CHAR(0)
      call parallelIO_create_dataset_group(pio_dataset_groupid, string, file_id, &
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

#if defined(PARALLELIO_LIB_WRITE)
      call parallelio_close_dataset_group(pio_dataset_groupid, file_id, &
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
#if defined(PARALLELIO_LIB_WRITE)
    string = trim(string) //CHAR(0)
      ! This opens existing dataset and creates it if needed
    call parallelIO_create_dataset_group(pio_dataset_groupid, string, file_id, &
                                          option%iowrite_group_id, ierr)
    grp_id = file_id
#else
    call h5eset_auto_f(OFF,hdf5_err)
    call h5gopen_f(file_id,string,grp_id,hdf5_err)
    if (hdf5_err /= 0) then
      call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
    endif
    call h5eset_auto_f(ON,hdf5_err)
#endif ! PARALLELIO_LIB_WRITE

#ifdef VAMSI_HDF5_WRITE
  endif
#endif
  
  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)

  ! loop over variables and write to file
  cur_variable => output_option%output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputGetVarFromArray(realization,global_vec,cur_variable%ivar, &
                                cur_variable%isubvar)
    string = cur_variable%name
    if (len_trim(cur_variable%units) > 0) then
      word = cur_variable%units
      call HDF5MakeStringCompabible(word)
      string = trim(string) // ' [' // trim(word) // ']'
    endif
    if (cur_variable%iformat == 0) then
      call HDF5WriteStructDataSetFromVec(string,realization, &
                                         global_vec,grp_id,H5T_NATIVE_DOUBLE)
    else
      call HDF5WriteStructDataSetFromVec(string,realization, &
                                         global_vec,grp_id,H5T_NATIVE_INTEGER)
    endif
    cur_variable => cur_variable%next
  enddo

  if (output_option%print_hdf5_velocities) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)
    call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    if (option%nphase > 1) then
        call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
        string = "Gas X-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
        string = "Gas Y-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
        string = "Gas Z-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)
    endif
  endif

  if (output_option%print_hdf5_flux_velocities) then

    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
        string = "Liquid X-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,X_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas X-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,X_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%ny > 1) then
        string = "Liquid Y-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,Y_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Y-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,Y_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%nz > 1) then
        string = "Liquid Z-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,Z_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Z-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,Z_DIRECTION,grp_id)
        endif
    endif
   
  endif

  call VecDestroy(global_vec,ierr)

#if defined(PARALLELIO_LIB_WRITE)
    call parallelio_close_dataset_group(pio_dataset_groupid, file_id, &
            option%iowrite_group_id, ierr)
    call parallelio_close_file(file_id, option%iowrite_group_id, ierr)
#else
#ifdef VAMSI_HDF5_WRITE
    if (mod(option%myrank,option%hdf5_write_group_size) == 0) then 
#endif
       call h5gclose_f(grp_id,hdf5_err)
       call h5fclose_f(file_id,hdf5_err)
#ifdef VAMSI_HDF5_WRITE
    endif
#endif
     call h5close_f(hdf5_err)
#endif !PARALLELIO_LIB_WRITE
#endif

  hdf5_first = PETSC_FALSE

end subroutine OutputHDF5

! ************************************************************************** !
!
! OutputMAD: Print to HDF5 file for MAD final output
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputMAD(realization)

  use Realization_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
 
#if !defined(PETSC_HAVE_HDF5)
  implicit none
  
  type(realization_type) :: realization

  call printMsg(realization%option,'')
  write(realization%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to ", &
        &"write HDF5 formatted structured grids.")')
  call printErrMsg(realization%option)
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

  type(realization_type) :: realization

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

  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option

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
  call OutputGetVarFromArray(realization,global_vec,LIQUID_PRESSURE,ZERO_INTEGER)
#ifdef ALL
  string = 'Pressure' // trim(option%group_prefix)
#else
  string = 'Pressure'
#endif
  call HDF5WriteStructDataSetFromVec(string,realization,global_vec,file_id,H5T_NATIVE_DOUBLE)

  call VecDestroy(global_vec,ierr)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif
end subroutine OutputMAD

#if defined(PETSC_HAVE_HDF5)
! ************************************************************************** !
!
! WriteHDF5FluxVelocities: Print flux velocities to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5FluxVelocities(name,realization,iphase,direction,file_id)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  use hdf5
  use HDF5_module

  implicit none
  
  character(len=32) :: name
  type(realization_type) :: realization
  PetscInt :: iphase
  PetscInt :: direction
  integer(HID_T) :: file_id

  PetscInt :: i, j, k
  PetscInt :: count, iconn
  PetscInt :: local_id, ghosted_id
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  
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
    
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option  

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
  
#if defined(PARALLELIO_LIB_WRITE)
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
  
  call PetscLogEventBegin(logging%event_output_coordinates_hdf5,ierr) 
#if defined(PARALLELIO_LIB_WRITE)

  name = trim(name) // CHAR(0)
  ! write out grid structure
  rank = 1
  dims = 0
  globaldims = 0
  ! x-direction

  !if (option%myrank == option%io_rank) then
  ! Only process 0 writes coordinates
  if (option%myrank == 0 ) then
     dims(1) = length
     globaldims(1) = length
  else
     dims(1) = 0
     globaldims(1) = length
  endif

 call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
 call parallelio_write_dataset(array, PIO_DOUBLE, rank, globaldims, dims, &
      file_id, name, option%iowrite_group_id, NONUNIFORM_CONTIGUOUS_WRITE, ierr)
 !call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                !hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
 call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

#else
!PARALLELIO_LIB_WRITE is not defined

#ifdef VAMSI_HDF5_WRITE
  if (mod(option%myrank,option%hdf5_write_group_size) == 0) then
#endif
                            
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

#ifdef VAMSI_HDF5_WRITE
  endif
#endif

#endif ! PARALLELIO_LIB_WRITE

  call PetscLogEventEnd(logging%event_output_coordinates_hdf5,ierr) 

end subroutine WriteHDF5Coordinates
#endif

! ************************************************************************** !
!
! GetCoordinates: Extracts coordinates of cells into a PetscVec
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GetCoordinates(grid,vec,direction)

  use Grid_module
  
  implicit none
  
  type(grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  
  if (direction == X_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%x(grid%nL2G(i))
    enddo
  else if (direction == Y_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%y(grid%nL2G(i))
    enddo
  else if (direction == Z_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%z(grid%nL2G(i))
    enddo
  endif
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine GetCoordinates

! ************************************************************************** !
!
! GetVertexCoordinates: Extracts vertex coordinates of cells into a PetscVec
! author: Gautam Bisht
! date: 11/01/2011
!
! ************************************************************************** !
subroutine GetVertexCoordinates(grid,vec,direction,option)

  use Grid_module
  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  type(option_type) :: option
  
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: indices(:)
  PetscReal, allocatable :: values(:)
  
  if (option%mycommsize == 1) then
    call VecGetArrayF90(vec,vec_ptr,ierr)
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%z
        enddo
    end select
    call VecRestoreArrayF90(vec,vec_ptr,ierr)
  else
    ! initialize to -999 to catch bugs
    call VecSet(vec,-999.d0,ierr)
    allocate(values(grid%unstructured_grid%num_vertices_local))
    allocate(indices(grid%unstructured_grid%num_vertices_local))
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%z
        enddo
    end select
    indices(:) = grid%unstructured_grid%vertex_ids_natural(:)-1
    call VecSetValues(vec,grid%unstructured_grid%num_vertices_local, &
                      indices,values,INSERT_VALUES,ierr)
    call VecAssemblyBegin(vec,ierr)
    deallocate(values)
    deallocate(indices)
    call VecAssemblyEnd(vec,ierr)
  endif
  
  
end subroutine GetVertexCoordinates

! ************************************************************************** !
!
! GetCellConnections: This routine returns a vector containing vertex ids
! in natural order of local cells.
! author: Gautam Bisht
! date: 11/01/2011
!
! ************************************************************************** !
subroutine GetCellConnectionsTecplot(grid, vec)

  use Grid_module
  use Unstructured_Grid_Aux_module

  implicit none
  
  type(grid_type) :: grid
  type(unstructured_grid_type),pointer :: ugrid
  Vec :: vec
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  
  ugrid => grid%unstructured_grid
  
  call GridVecGetArrayF90(grid, vec, vec_ptr, ierr)

  ! initialize
  vec_ptr = -999.d0
  do local_id=1, ugrid%nlmax
    ghosted_id = local_id
    select case(ugrid%cell_type(ghosted_id))
      case(HEX_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 8
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case(WEDGE_TYPE)
        offset = (local_id-1)*8
        vec_ptr(offset + 1) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(1,local_id))
        vec_ptr(offset + 2) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(1,local_id))
        vec_ptr(offset + 3) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(4,local_id))
        vec_ptr(offset + 4) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(4,local_id))
        vec_ptr(offset + 5) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(3,local_id))
        vec_ptr(offset + 6) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(2,local_id))
        vec_ptr(offset + 7) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(5,local_id))
        vec_ptr(offset + 8) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(6,local_id))
      case (PYR_TYPE)
        offset = (local_id-1)*8
        ! from Tecplot 360 Data Format Guide
        ! n1=vert1,n2=vert2,n3=vert3,n4=vert4,n5=n6=n7=n8=vert5
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        do ivertex = 5, 8
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(5,local_id))
        enddo
      case (TET_TYPE)
        offset = (local_id-1)*8
        ! from Tecplot 360 Data Format Guide
        ! n1=vert1,n2=vert2,n3=n4=vert3,n5=vert5=n6=n7=n8=vert4
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        vec_ptr(offset + 4) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(3,local_id))
        do ivertex = 5, 8
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(4,local_id))
        enddo
      case (QUAD_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case (TRI_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        ivertex = 4
        vec_ptr(offset + ivertex) = &
          ugrid%vertex_ids_natural(ugrid%cell_vertices(3,local_id))
    end select
  enddo

  call GridVecRestoreArrayF90(grid, vec, vec_ptr, ierr)

end subroutine GetCellConnectionsTecplot

! ************************************************************************** !
!> This routine returns a vector containing vertex ids in natural order of
!! local cells for unstructured grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/31/12
! ************************************************************************** !
subroutine GetCellConnections(grid, vec)

  use Grid_module
  use Unstructured_Grid_Aux_module

  implicit none
  
  type(grid_type) :: grid
  type(unstructured_grid_type),pointer :: ugrid
  Vec :: vec
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  
  ugrid => grid%unstructured_grid
  
  call GridVecGetArrayF90(grid, vec, vec_ptr, ierr)

  ! initialize
  vec_ptr = -999.d0
  do local_id=1, ugrid%nlmax
    ghosted_id = local_id
    select case(ugrid%cell_type(ghosted_id))
      case(HEX_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 8
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case(WEDGE_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 6
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        vec_ptr(offset + 7) = 0
        vec_ptr(offset + 8) = 0
      case (PYR_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 5
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        do ivertex = 6, 8
          vec_ptr(offset + ivertex) = 0
        enddo
      case (TET_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        do ivertex = 5, 8
          vec_ptr(offset + ivertex) = 0
        enddo
      case (QUAD_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case (TRI_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        ivertex = 4
        vec_ptr(offset + 4) = 0
    end select
  enddo

  call GridVecRestoreArrayF90(grid, vec, vec_ptr, ierr)

end subroutine GetCellConnections

! ************************************************************************** !
!
! ConvertArrayToNatural: Converts an array  to natural ordering
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine ConvertArrayToNatural(indices,array,local_size,global_size,option)

  use Option_module
  
  implicit none
  
  PetscInt :: local_size, global_size
  PetscInt :: indices(:)
  PetscReal, pointer :: array(:)
  type(option_type) :: option
  
  Vec :: natural_vec
  PetscInt, allocatable :: indices_zero_based(:)
  PetscReal, pointer :: vec_ptr(:)
  
  call VecCreate(option%mycomm,natural_vec,ierr)
  call VecSetSizes(natural_vec,PETSC_DECIDE,global_size,ierr)
  call VecSetType(natural_vec,VECMPI,ierr)

  allocate(indices_zero_based(local_size))
  indices_zero_based(1:local_size) = indices(1:local_size)-1

  call VecSetValues(natural_vec,local_size,indices_zero_based, &
                    array,INSERT_VALUES,ierr)

  call VecAssemblyBegin(natural_vec,ierr)
  call VecAssemblyEnd(natural_vec,ierr)

  call VecGetLocalSize(natural_vec,local_size,ierr)
  deallocate(array)
  allocate(array(local_size))
  
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)
  array(1:local_size) = vec_ptr(1:local_size)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)

  call VecDestroy(natural_vec,ierr)
  
end subroutine ConvertArrayToNatural

! ************************************************************************** !
!
! OutputGetVarFromArrayAtCoord: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function OutputGetVarFromArrayAtCoord(realization,ivar,isubvar,x,y,z, &
                                      num_cells,ghosted_ids,isubvar1)

  use Realization_module
  use Grid_module
  use Option_module

  implicit none
  
  PetscReal :: OutputGetVarFromArrayAtCoord
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscReal :: x,y,z
  PetscInt :: num_cells
  PetscInt :: ghosted_ids(num_cells)

  type(grid_type), pointer :: grid
  PetscInt :: icell, ghosted_id
  PetscReal :: dx, dy, dz
  PetscReal :: value, sum_value
  PetscReal :: weight, sum_weight, sum_root
  
  sum_value = 0.d0
  sum_weight = 0.d0
  
  grid => realization%patch%grid

  do icell=1, num_cells
    ghosted_id = ghosted_ids(icell)
    dx = x-grid%x(ghosted_id)
    dy = y-grid%y(ghosted_id)
    dz = z-grid%z(ghosted_id)
    sum_root = sqrt(dx*dx+dy*dy+dz*dz)
    value = 0.d0
    value = RealizGetDatasetValueAtCell(realization,ivar,isubvar,ghosted_id, &
      isubvar1)
    if (sum_root < 1.d-40) then ! bail because it is right on this coordinate
      sum_weight = 1.d0
      sum_value = value
      exit
    endif
    weight = 1.d0/sum_root
    sum_weight = sum_weight + weight
    sum_value = sum_value + weight * value
  enddo
  
  OutputGetVarFromArrayAtCoord = sum_value/sum_weight

end function OutputGetVarFromArrayAtCoord

! ************************************************************************** !
!
! OutputGetVarFromArray: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputGetVarFromArray1(realization,vec,ivar,isubvar,isubvar1)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1

  call PetscLogEventBegin(logging%event_output_get_var_from_array,ierr) 
                        
  call RealizationGetDataset(realization,vec,ivar,isubvar,isubvar1)

  call PetscLogEventEnd(logging%event_output_get_var_from_array,ierr) 
  
end subroutine OutputGetVarFromArray1

! ************************************************************************** !
!
! OutputGetCellCenteredVelocities: Computes the cell-centered velocity component 
!                            as an averages of cell face velocities
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputGetCellCenteredVelocities(realization,vec,iphase,direction)

  use Realization_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Patch_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: direction
  PetscInt :: iphase
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscReal :: velocity, area
  PetscReal :: average, sum, max, min, std_dev
  PetscInt :: max_loc, min_loc
  
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, allocatable :: sum_area(:)
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  call PetscLogEventBegin(logging%event_output_get_cell_vel,ierr) 
                            
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
    
  allocate(sum_area(grid%nlmax))
  sum_area(1:grid%nlmax) = 0.d0

  call VecSet(vec,0.d0,ierr)
  call VecGetArrayF90(vec,vec_ptr,ierr)

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
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(direction,iconn)
      velocity = patch%internal_velocities(iphase,sum_connection)* &
                 area
      if (local_id_up > 0) then
        vec_ptr(local_id_up) = vec_ptr(local_id_up) + velocity
        sum_area(local_id_up) = sum_area(local_id_up) + dabs(area)
      endif
      if (local_id_dn > 0) then
        vec_ptr(local_id_dn) = vec_ptr(local_id_dn) + velocity
        sum_area(local_id_dn) = sum_area(local_id_dn) + dabs(area)
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
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(direction,iconn)
      vec_ptr(local_id) = vec_ptr(local_id)+ &
                          patch%boundary_velocities(iphase,sum_connection)* &
                          area
      sum_area(local_id) = sum_area(local_id) + dabs(area)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! divide by total area
  do local_id=1,grid%nlmax
    if (sum_area(local_id) > 0.d0) &
      vec_ptr(local_id) = vec_ptr(local_id)/sum_area(local_id)*output_option%tconv
  enddo
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  deallocate(sum_area)

  call PetscLogEventEnd(logging%event_output_get_cell_vel,ierr) 

end subroutine OutputGetCellCenteredVelocities

! ************************************************************************** !
!
! ComputeFlowMassBalance: 
! author: Glenn Hammond
! date: 03/11/08
!
! ************************************************************************** !
subroutine ComputeFlowMassBalance(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Patch_module
  use Discretization_module

  implicit none
  
  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(discretization_type), pointer :: discretization
  type(output_option_type), pointer :: output_option
  PetscInt :: iconn, i, sum_connection
  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscReal :: flux
  Vec :: global_vec
  Vec :: mass_vec
  Vec :: density_loc
  Vec :: total_mass_vec
  PetscReal :: average, sum, max, min, std_dev
  PetscInt :: max_loc, min_loc
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:), den_loc_p(:)
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  discretization => realization%discretization

  call DiscretizationDuplicateVector(discretization,field%porosity0,total_mass_vec)
  call DiscretizationDuplicateVector(discretization,field%porosity0,mass_vec)
  call DiscretizationDuplicateVector(discretization,field%porosity0,global_vec)
  call DiscretizationDuplicateVector(discretization,field%porosity_loc,density_loc)
  
  call OutputGetVarFromArray(realization,global_vec,LIQUID_DENSITY,ZERO_INTEGER)
  call DiscretizationGlobalToLocal(discretization,global_vec,density_loc,ONEDOF)

  call VecSet(mass_vec,0.d0,ierr)
  call VecGetArrayF90(mass_vec,vec_ptr,ierr)
  call VecGetArrayF90(density_loc,den_loc_p,ierr)

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
      flux = patch%internal_velocities(1,iconn)* &
               cur_connection_set%area(iconn)
      if (patch%internal_velocities(1,iconn) > 0.d0) then
        flux = flux*den_loc_p(ghosted_id_up)
      else
        flux = flux*den_loc_p(ghosted_id_dn)
      endif
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
                          patch%boundary_velocities(1,sum_connection)* &
                          cur_connection_set%area(iconn)* &
                          den_loc_p(grid%nL2G(local_id))
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(mass_vec,vec_ptr,ierr)
  call VecRestoreArrayF90(density_loc,den_loc_p,ierr)

  ! scale by mass of water in cell
  call DiscretizationLocalToGlobal(discretization,density_loc,global_vec,ONEDOF)
  call VecPointWiseMult(total_mass_vec,global_vec,field%volume,ierr) ! global_vec is density
  call DiscretizationLocalToGlobal(discretization,field%porosity_loc,global_vec,ONEDOF)
  call VecPointWiseMult(total_mass_vec,total_mass_vec,global_vec,ierr)
  call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
  call VecPointWiseMult(total_mass_vec,total_mass_vec,global_vec,ierr) ! global_vec is saturation

  ! have to divide through on our own since zero values exist on inactive cells
  call VecGetArrayF90(mass_vec,vec_ptr,ierr)
  call VecGetArrayF90(total_mass_vec,vec2_ptr,ierr)
  do i=1,grid%nlmax
    if (dabs(vec2_ptr(i)) > 1.d-40) then
      vec_ptr(i) = vec_ptr(i)/vec2_ptr(i)
    else
      vec_ptr(i) = 0.d0
    endif
  enddo
  call VecRestoreArrayF90(mass_vec,vec_ptr,ierr)
  call VecRestoreArrayF90(total_mass_vec,vec2_ptr,ierr)

  string = 'mass_balance.tec'
  call OutputVectorTecplot(string,string,realization,mass_vec)

  call VecSum(mass_vec,sum,ierr)
  average = sum/real(grid%nmax)
  call VecSet(global_vec,average,ierr)
  call VecMax(mass_vec,max_loc,max,ierr)
  call VecMin(mass_vec,min_loc,min,ierr)
  call VecAYPX(global_vec,-1.d0,mass_vec,ierr)
  call VecNorm(global_vec,NORM_2,std_dev,ierr)
  string = 'Mass Balance'
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

  call VecDestroy(total_mass_vec,ierr)
  call VecDestroy(mass_vec,ierr)
  call VecDestroy(global_vec,ierr)
  call VecDestroy(density_loc,ierr)

end subroutine ComputeFlowMassBalance

! ************************************************************************** !
!
! OutputMassBalance: Print to Tecplot POINT format
! author: Glenn Hammond
! date: 06/18/08
!
! ************************************************************************** !  
subroutine OutputMassBalance(realization)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Utility_module
  
  use Richards_module
  use Mphase_module
  use Immis_module
  use Miscible_module
  use THC_module
  use THMC_module
  use Reactive_Transport_module
  use General_module
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  
  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option  
  type(coupler_type), pointer :: coupler
  type(global_auxvar_type), pointer :: global_aux_vars_bc_or_ss(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc_or_ss(:)
  
  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word, units
  character(len=MAXSTRINGLENGTH) :: string
  character(len=4) :: strcol
  PetscInt :: fid = 86
  PetscInt :: ios
  PetscInt :: i,icol
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iconn
  PetscInt :: offset
  PetscInt :: iphase, ispec
  PetscInt :: icomp
  PetscReal :: sum_area(4)
  PetscReal :: sum_area_global(4)
  PetscReal :: sum_kg(realization%option%nflowspec,realization%option%nphase)
  PetscReal :: sum_kg_global(realization%option%nflowspec,realization%option%nphase)
  PetscReal :: sum_mol(realization%option%ntrandof,realization%option%nphase)
  PetscReal :: sum_mol_global(realization%option%ntrandof,realization%option%nphase)
  PetscMPIInt :: int_mpi
  PetscBool :: bcs_done
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  reaction => realization%reaction
  output_option => realization%output_option
  
 
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.dat'
  else
    filename = 'mass_balance' // trim(option%group_prefix) // '.dat'
  endif
  
  ! open file
  if (option%myrank == option%io_rank) then

!geh    option%io_buffer = '--> write tecplot mass balance file: ' // trim(filename)
!geh    call printMsg(option)    

    if (output_option%print_column_ids) then
      icol = 1
    else
      icol = -1
    endif
  
    if (mass_balance_first .or. .not.FileExists(filename)) then
      open(unit=fid,file=filename,action="write",status="replace")

      ! write header
      write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'  
      
      header = ''
      if (option%iflowmode > 0) then
        call OutputAppendToHeader(header,'dt_flow',output_option%tunit,'',icol)
      endif
      
      if (option%ntrandof > 0) then
        call OutputAppendToHeader(header,'dt_tran',output_option%tunit,'',icol)
      endif
      write(fid,'(a)',advance="no") trim(header)
      
      header = ''
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          call OutputAppendToHeader(header,'Global Water Mass','[kg]','',icol)
          
        case(THC_MODE)
          call OutputAppendToHeader(header,'Global Water Mass in Liquid Phase', &
                                    '[kg]','',icol)
        case(THMC_MODE)
          call OutputAppendToHeader(header,'Global Water Mass in Liquid Phase', &
                                    '[kg]','',icol)
        case(G_MODE)
          call OutputAppendToHeader(header,'Global Water Mass in Liquid Phase', &
                                    '[mol]','',icol)
          call OutputAppendToHeader(header,'Global Air Mass in Liquid Phase', &
                                    '[mol]','',icol)
          call OutputAppendToHeader(header,'Global Water Mass in Gas Phase', &
                                    '[mol]','',icol)
          call OutputAppendToHeader(header,'Global Air Mass in Gas Phase', &
                                    '[mol]','',icol)
        case(MPH_MODE)
          call OutputAppendToHeader(header,'Global Water Mass in Water Phase', &
                                    '[kmol]','',icol)
          call OutputAppendToHeader(header,'Global CO2 Mass in Water Phase', &
                                    '[kmol]','',icol)
          call OutputAppendToHeader(header,'Global Water Mass in Gas Phase', &
                                    '[kmol]','',icol)
          call OutputAppendToHeader(header,'Global CO2 Mass in Gas Phase', &
                                    '[kmol]','',icol)
        case(IMS_MODE)
          call OutputAppendToHeader(header,'Global Water Mass in Water Phase', &
                                    '[kmol]','',icol)
          call OutputAppendToHeader(header,'Global CO2 Mass in Gas Phase', &
                                    '[kmol]','',icol)
        case(MIS_MODE)
          call OutputAppendToHeader(header,'Global Water Mass in Liquid Phase', &
                                    '[kg]','',icol)
          call OutputAppendToHeader(header,'Global Glycol Mass in Liquid Phase', &
                                    '[kg]','',icol)
      end select
      write(fid,'(a)',advance="no") trim(header)

      if (option%ntrandof > 0) then
        header = ''    
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            string = 'Global ' // trim(reaction%primary_species_names(i))
            call OutputAppendToHeader(header,string,'[mol]','',icol)
          endif
        enddo
        write(fid,'(a)',advance="no") trim(header)
      endif
      
      coupler => patch%boundary_conditions%first
      bcs_done = PETSC_FALSE
      do
        if (.not.associated(coupler)) then
          if (bcs_done) then
            exit
          else
            bcs_done = PETSC_TRUE
            if (associated(patch%source_sinks)) then
              coupler => patch%source_sinks%first
              if (.not.associated(coupler)) exit
            else
              exit
            endif
          endif
        endif

        header = ''
        select case(option%iflowmode)
          case(RICHARDS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,'[kg]','',icol)
            
            units = '[kg/' // trim(output_option%tunit) // ']'
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
          case(THC_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,'[kg]','',icol)
            
            units = '[kg/' // trim(output_option%tunit) // ']'
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
          case(THMC_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,'[kg]','',icol)
            
            units = '[kg/' // trim(output_option%tunit) // ']'
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
          case(MIS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,'[kg]','',icol)
            string = trim(coupler%name) // ' Glycol Mass'
            call OutputAppendToHeader(header,string,'[kg]','',icol)
            
            units = '[kg/' // trim(output_option%tunit) // ']'
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
            string = trim(coupler%name) // ' Glycol Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
          case(G_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,'[mol]','',icol)
            string = trim(coupler%name) // ' Air Mass'
            call OutputAppendToHeader(header,string,'[mol]','',icol)
            
            units = '[mol/' // trim(output_option%tunit) // ']'
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
            string = trim(coupler%name) // ' Air Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
          case(MPH_MODE,IMS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,'[kmol]','',icol)
            string = trim(coupler%name) // ' CO2 Mass'
            call OutputAppendToHeader(header,string,'[kmol]','',icol)
            
            units = '[kmol/' // trim(output_option%tunit) // ']'
            string = trim(coupler%name) // ' Water Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
            string = trim(coupler%name) // ' CO2 Mass'
            call OutputAppendToHeader(header,string,units,'',icol)
        end select
        write(fid,'(a)',advance="no") trim(header)
        
        if (option%ntrandof > 0) then
          header = ''    
          do i=1,reaction%naqcomp
            if (reaction%primary_species_print(i)) then
              string = trim(coupler%name) // ' ' // &
                       trim(reaction%primary_species_names(i))
              call OutputAppendToHeader(header,string,'[mol]','',icol)
            endif
          enddo
          write(fid,'(a)',advance="no") trim(header)
          
          header = ''    
          units = '[mol/' // trim(output_option%tunit) // ']'
          do i=1,reaction%naqcomp
            if (reaction%primary_species_print(i)) then
              string = trim(coupler%name) // ' ' // &
                       trim(reaction%primary_species_names(i))
              call OutputAppendToHeader(header,string,units,'',icol)
            endif
          enddo
          write(fid,'(a)',advance="no") trim(header)
        endif
        coupler => coupler%next
      
      enddo
      
#ifdef COMPUTE_INTERNAL_MASS_FLUX
      do offset = 1, 4
        write(word,'(i6)') offset*100
        select case(option%iflowmode)
          case(RICHARDS_MODE)
            write(fid,'(a)',advance="no") ',"' // &
              trim(adjustl(word)) // 'm Water Mass [kg]"'
          case(THC_MODE)
            write(fid,'(a)',advance="no") ',"' // &
              trim(adjustl(word)) // 'm Water Mass [kg]"'
          case(THMC_MODE)
            write(fid,'(a)',advance="no") ',"' // &
              trim(adjustl(word)) // 'm Water Mass [kg]"'
        end select
        
        if (option%ntrandof > 0) then
          do i=1,reaction%naqcomp
            if (reaction%primary_species_print(i)) then
              write(fid,'(a)',advance="no") ',"' // &
                  trim(adjustl(word)) // 'm ' // &
                  trim(reaction%primary_species_names(i)) // ' [mol]"'
            endif
          enddo
        endif
      enddo
#endif      
      write(fid,'(a)') '' 
    else
      open(unit=fid,file=filename,action="write",status="old",position="append")
    endif 
    
  endif     

100 format(100es16.8)
110 format(100es16.8)

  ! write time
  if (option%myrank == option%io_rank) then
    write(fid,100,advance="no") option%time/output_option%tconv
  endif
  
  if (option%nflowdof > 0) then
    if (option%myrank == option%io_rank) &
      write(fid,100,advance="no") option%flow_dt/output_option%tconv
  endif
  if (option%ntrandof > 0) then
    if (option%myrank == option%io_rank) &
      write(fid,100,advance="no") option%tran_dt/output_option%tconv
  endif
  
! print out global mass balance

  if (option%nflowdof > 0) then
    sum_kg = 0.d0
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        call RichardsComputeMassBalance(realization,sum_kg(1,:))
      case(THC_MODE)
        call THCComputeMassBalance(realization,sum_kg(1,:))
      case(THMC_MODE)
        call THMCComputeMassBalance(realization,sum_kg(1,:))
      case(MIS_MODE)
        call MiscibleComputeMassBalance(realization,sum_kg(:,1))
      case(MPH_MODE)
        call MphaseComputeMassBalance(realization,sum_kg(:,:))
      case(IMS_MODE)
        call ImmisComputeMassBalance(realization,sum_kg(:,1))
      case(G_MODE)
        option%io_buffer = 'Mass balance calculations not yet implemented for General Mode'
        call printErrMsg(option)
        call GeneralComputeMassBalance(realization,sum_kg(1,:))
    end select
    int_mpi = option%nflowspec*option%nphase
    call MPI_Reduce(sum_kg,sum_kg_global, &
                    int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)
                        
    if (option%myrank == option%io_rank) then
      select case(option%iflowmode)
        case(RICHARDS_MODE,MPH_MODE,FLASH2_MODE,IMS_MODE,MIS_MODE,G_MODE, &
             THC_MODE,THMC_MODE)
          do iphase = 1, option%nphase
            do ispec = 1, option%nflowspec
              write(fid,110,advance="no") sum_kg_global(ispec,iphase)
            enddo
          enddo
      end select
    endif
  endif
  
  if (option%ntrandof > 0) then
    sum_mol = 0.d0
    call RTComputeMassBalance(realization,sum_mol)
    int_mpi = option%nphase*option%ntrandof
    call MPI_Reduce(sum_mol,sum_mol_global,int_mpi, &
                    MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)

    if (option%myrank == option%io_rank) then
      do iphase = 1, option%nphase
        do icomp = 1, reaction%naqcomp
          if (reaction%primary_species_print(icomp)) then
            write(fid,110,advance="no") sum_mol_global(icomp,iphase)
          endif
        enddo
      enddo
    endif
  endif
  
  coupler => patch%boundary_conditions%first
  global_aux_vars_bc_or_ss => patch%aux%Global%aux_vars_bc
  if (option%ntrandof > 0) then
    rt_aux_vars_bc_or_ss => patch%aux%RT%aux_vars_bc
  endif    
  bcs_done = PETSC_FALSE
  do 
    if (.not.associated(coupler)) then
      if (bcs_done) then
        exit
      else
        bcs_done = PETSC_TRUE
        if (associated(patch%source_sinks)) then
          coupler => patch%source_sinks%first
          if (.not.associated(coupler)) exit
          global_aux_vars_bc_or_ss => patch%aux%Global%aux_vars_ss
          if (option%ntrandof > 0) then
            rt_aux_vars_bc_or_ss => patch%aux%RT%aux_vars_ss
          endif    
        else
          exit
        endif
      endif
    endif

    offset = coupler%connection_set%offset
    
    if (option%nflowdof > 0) then

#if 0
! compute the total area of the boundary condition
      if (.not.bcs_done) then
        sum_area = 0.d0
        do iconn = 1, coupler%connection_set%num_connections
          sum_area(1) = sum_area(1) + &
            coupler%connection_set%area(iconn)
          if (global_aux_vars_bc_or_ss(offset+iconn)%sat(1) >= 0.5d0) then
            sum_area(2) = sum_area(2) + &
              coupler%connection_set%area(iconn)
          endif
          if (global_aux_vars_bc_or_ss(offset+iconn)%sat(1) > 0.99d0) then
            sum_area(3) = sum_area(3) + &
              coupler%connection_set%area(iconn)
          endif
          sum_area(4) = sum_area(4) + &
            coupler%connection_set%area(iconn)* &
            global_aux_vars_bc_or_ss(offset+iconn)%sat(1)
        enddo

        call MPI_Reduce(sum_area,sum_area_global, &
                        FOUR_INTEGER_MPI,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        option%io_rank,option%mycomm,ierr)
                          
        if (option%myrank == option%io_rank) then
          print *
          write(word,'(es16.6)') sum_area_global(1)
          print *, 'Total area in ' // trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          write(word,'(es16.6)') sum_area_global(2)
          print *, 'Total half-saturated area in '// &
                   trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          write(word,'(es16.6)') sum_area_global(3)
          print *, 'Total saturated area in '// trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          write(word,'(es16.6)') sum_area_global(4)
          print *, 'Total saturation-weighted area [=sum(saturation*area)] in '//&
                     trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          print *
        endif
      endif
#endif

      select case(option%iflowmode)
        case(RICHARDS_MODE)
          ! print out cumulative H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_aux_vars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global
          endif

          ! print out H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_aux_vars_bc_or_ss(offset+iconn)%mass_balance_delta
          enddo
          
          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
          sum_kg = sum_kg*FMWH2O

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global*output_option%tconv
          endif

        case(THC_MODE)
          ! print out cumulative H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_aux_vars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global
          endif

          ! print out H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_aux_vars_bc_or_ss(offset+iconn)%mass_balance_delta
          enddo
          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
          sum_kg = sum_kg*FMWH2O

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global*output_option%tconv
          endif

        case(THMC_MODE)
          ! print out cumulative H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_aux_vars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global
          endif

          ! print out H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_aux_vars_bc_or_ss(offset+iconn)%mass_balance_delta
          enddo
          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
          sum_kg = sum_kg*FMWH2O

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global*output_option%tconv
          endif

        case(MIS_MODE)
          ! print out cumulative mixture flux
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_aux_vars_bc_or_ss(offset+iconn)%mass_balance(icomp,1)
            enddo
            
            if (icomp == 1) then
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O
            else
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWGLYC
            endif
            
            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
          
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)
            endif
          enddo

          ! print out mixture flux
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_aux_vars_bc_or_ss(offset+iconn)%mass_balance_delta(icomp,1)
            enddo
            
        !   mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o/glycol
            if (icomp == 1) then
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O
            else
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWGLYC
            endif

            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)*output_option%tconv
            endif
          enddo

        case(MPH_MODE)
        ! print out cumulative H2O & CO2 fluxes in kmol and kmol/time
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_aux_vars_bc_or_ss(offset+iconn)%mass_balance(icomp,1)
            enddo
            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)
            endif
          enddo
          
        ! print out H2O & CO2 fluxes in kmol and kmol/time
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_aux_vars_bc_or_ss(offset+iconn)%mass_balance_delta(icomp,1)
            enddo

          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
!           sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O ! <<---fix for multiphase!

            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)*output_option%tconv
            endif
          enddo

        case(IMS_MODE)
        ! print out cumulative H2O & CO2 fluxes
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_aux_vars_bc_or_ss(offset+iconn)%mass_balance(icomp,1)
            enddo
            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)
            endif
          enddo
          
        ! print out H2O & CO2 fluxes
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_aux_vars_bc_or_ss(offset+iconn)%mass_balance_delta(icomp,1)
            enddo

          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
!           sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O ! <<---fix for multiphase!

            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)*output_option%tconv
            endif
          enddo
        case(G_MODE)
      end select
    endif
    
    if (option%ntrandof > 0) then

      ! print out cumulative boundary flux
      sum_mol = 0.d0
      do iconn = 1, coupler%connection_set%num_connections
        sum_mol = sum_mol + rt_aux_vars_bc_or_ss(offset+iconn)%mass_balance
      enddo

      int_mpi = option%nphase*option%ntrandof
      call MPI_Reduce(sum_mol,sum_mol_global,int_mpi, &
                      MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)

      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        do iphase = 1, option%nphase
          do icomp = 1, reaction%naqcomp
            if (reaction%primary_species_print(icomp)) then
              write(fid,110,advance="no") -sum_mol_global(icomp,iphase)
            endif
          enddo
        enddo
      endif
    
      ! print out boundary flux
      sum_mol = 0.d0
      do iconn = 1, coupler%connection_set%num_connections
        sum_mol = sum_mol + rt_aux_vars_bc_or_ss(offset+iconn)%mass_balance_delta 
      enddo

      int_mpi = option%nphase*option%ntrandof
      call MPI_Reduce(sum_mol,sum_mol_global,int_mpi, &
                      MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
                      
      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        do iphase = 1, option%nphase
          do icomp = 1, reaction%naqcomp
            if (reaction%primary_species_print(icomp)) then
              write(fid,110,advance="no") -sum_mol_global(icomp,iphase)* &
                                          output_option%tconv
            endif
          enddo
        enddo
      endif
    endif

    coupler => coupler%next
  
  enddo

#ifdef COMPUTE_INTERNAL_MASS_FLUX

  do offset = 1, 4
    iconn = offset*20-1

    if (option%nflowdof > 0) then
      sum_kg = 0.d0
      sum_kg = sum_kg + patch%aux%Global%aux_vars(iconn)%mass_balance

      int_mpi = option%nphase
      call MPI_Reduce(sum_kg,sum_kg_global, &
                      int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
                          
      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        write(fid,110,advance="no") -sum_kg_global
      endif
    endif
    
    if (option%ntrandof > 0) then

      sum_mol = 0.d0
      sum_mol = sum_mol + patch%aux%RT%aux_vars(iconn)%mass_balance

      int_mpi = option%nphase*option%ntrandof
      call MPI_Reduce(sum_mol,sum_mol_global,int_mpi, &
                      MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)

      if (option%myrank == option%io_rank) then
        do iphase = 1, option%nphase
          do icomp = 1, reaction%naqcomp
            if (reaction%primary_species_print(icomp)) then
              ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_mol_global(icomp,iphase)
            endif
          enddo
        enddo
      endif
    endif
  enddo
#endif
  
  if (option%myrank == option%io_rank) then
    write(fid,'(a)') ''
    close(fid)
  endif
  
  mass_balance_first = PETSC_FALSE

end subroutine OutputMassBalance

! ************************************************************************** !
!
! ComputeFlowCellVelocityStats: 
! author: Glenn Hammond
! date: 03/11/08
!
! ************************************************************************** !
subroutine ComputeFlowCellVelocityStats(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Patch_module
  use Discretization_module

  implicit none
  
  type(realization_type) :: realization
  
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
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  discretization => realization%discretization
    
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
subroutine ComputeFlowFluxVelocityStats(realization)
!geh - specifically, the flow velocities at the interfaces between cells
 
  use Realization_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  
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

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
    
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
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
! OutputPermeability: Print vectors for permeability
! author: Glenn Hammond
! date: 08/25/09
!
! ************************************************************************** !
subroutine OutputPermeability(realization)

  use Realization_module
  use Option_module
  use Discretization_module
  use Material_module

  implicit none

  type(realization_type) :: realization
  
  PetscBool :: print_all_three
  PetscInt :: material_property_id
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(material_property_type), pointer :: material_property
  
  option => realization%option

  print_all_three = PETSC_FALSE
  ! check for anisotripic permeabilities  
  do material_property_id = 1, size(realization%material_property_array)
    material_property => &
      realization%material_property_array(material_property_id)%ptr
    if (associated(material_property)) then
      if (.not.material_property%isotropic_permeability) then
        print_all_three = PETSC_TRUE
      endif
    endif
  enddo
  
  if (print_all_three) then
    if (len_trim(option%group_prefix) > 1) then
      string = 'permeabilityX-' // trim(option%group_prefix) // '.tec'
    else
      string = 'permeabilityX.tec'
    endif
    call DiscretizationLocalToGlobal(realization%discretization, &
                                     realization%field%perm_xx_loc, &
                                     realization%field%work,ONEDOF)
    call OutputVectorTecplot(string,string,realization,realization%field%work)
    if (len_trim(option%group_prefix) > 1) then
      string = 'permeabilityY-' // trim(option%group_prefix) // '.tec'
    else
      string = 'permeabilityY.tec'
    endif
    call DiscretizationLocalToGlobal(realization%discretization, &
                                     realization%field%perm_yy_loc, &
                                     realization%field%work,ONEDOF)
    call OutputVectorTecplot(string,string,realization,realization%field%work)
    if (len_trim(option%group_prefix) > 1) then
      string = 'permeabilityZ-' // trim(option%group_prefix) // '.tec'
    else
      string = 'permeabilityZ.tec'
    endif
    call DiscretizationLocalToGlobal(realization%discretization, &
                                     realization%field%perm_zz_loc, &
                                     realization%field%work,ONEDOF)
    call OutputVectorTecplot(string,string,realization,realization%field%work)
  else
    if (len_trim(option%group_prefix) > 1) then
      string = 'permeability-' // trim(option%group_prefix) // '.tec'
    else
      string = 'permeability.tec'
    endif
    call DiscretizationLocalToGlobal(realization%discretization, &
                                     realization%field%perm_xx_loc, &
                                     realization%field%work,ONEDOF)
    call OutputVectorTecplot(string,string,realization,realization%field%work)
  endif
  
end subroutine OutputPermeability

! ************************************************************************** !
!
! OutputPrintCouplers: Prints values of auxilliary variables associated with
!                      couplers (boundary and initial conditions, source
!                      sinks).  Note that since multiple connections for
!                      couplers can exist for a single cell, the latter will
!                      overwrite the former.
! author: Glenn Hammond
! date: 11/02/11
!
! ************************************************************************** !
subroutine OutputPrintCouplers(realization,istep)

  use Realization_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Debug_module
  use Field_module
  use Patch_module
  use Level_module
  use Grid_module
  use Input_module

  type(realization_type) :: realization
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
  
  
  option => realization%option
  flow_debug => realization%debug
  field => realization%field

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
    
    cur_level => realization%level_list%first
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
    call OutputVectorTecplot(string,word,realization,field%work)
      
  enddo

end subroutine OutputPrintCouplers

! ************************************************************************** !
!> This subroutine prints a HDF5 file.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/31/12
! ************************************************************************** !
subroutine OutputHDF5UGrid(realization)

  use Realization_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  
  type(realization_type) :: realization

  call printMsg(realization%option,'')
  write(realization%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to &
        &write HDF5 formatted structured grids Darn.")')
  call printErrMsg(realization%option)
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

  type(realization_type) :: realization

#if defined(PARALLELIO_LIB_WRITE)
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
  
  PetscMPIInt :: hdf5_flag

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
  PetscFortranAddr :: app_ptr
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type

  discretization => realization%discretization
  patch => realization%patch
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option


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

#if defined(PARALLELIO_LIB_WRITE)

  if (.not.first) then
    filename = trim(filename) // CHAR(0)
    call parallelio_open_file(filename, option%iowrite_group_id, &
                              FILE_READWRITE, file_id, ierr)
    if (file_id == -1) first = PETSC_TRUE
  endif
  if (first) then
    filename = trim(filename) // CHAR(0)
    call parallelio_open_file(filename, option%iowrite_group_id, &
                              FILE_CREATE, file_id, ierr)
  endif

  if (first) then
    option%io_buffer = '--> creating hdf5 output file: ' // filename
  else
    option%io_buffer = '--> appending to hdf5 output file: ' // filename
  endif
  call printMsg(option)

  if (first) then
    ! create a group for the coordinates data set
    string = "Domain" // CHAR(0)
    call parallelIO_create_dataset_group(pio_dataset_groupid, string, file_id, &
                                         option%iowrite_group_id, ierr)
    ! set grp_id here
    ! As we already created the group, we will use file_id as group_id
    grp_id = file_id
    call WriteHDF5CoordinatesUGrid(grid,option,grp_id)
    call parallelio_close_dataset_group(pio_dataset_groupid, file_id, &
            option%iowrite_group_id, ierr)
  endif

#else

  !
  !        not(PARALLELIO_LIB_WRITE)
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
    option%io_buffer = '--> creating hdf5 output file: ' // filename
  else
    option%io_buffer = '--> appending to hdf5 output file: ' // filename
  endif
  call printMsg(option)

  if (first) then
    ! create a group for the coordinates data set
    string = "Domain"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
    call WriteHDF5CoordinatesUGrid(grid,option,grp_id)
    call h5gclose_f(grp_id,hdf5_err)
  endif
#endif ! PARALLELIO_LIB_WRITE


    ! create a group for the data set
    write(string,'(''Time:'',es13.5,x,a1)') &
          option%time/output_option%tconv,output_option%tunit
    if (len_trim(output_option%plot_name) > 2) then
      string = trim(string) // ' ' // output_option%plot_name
    endif
#if defined(PARALLELIO_LIB_WRITE)
    string = trim(string) //CHAR(0)
      ! This opens existing dataset and creates it if needed
    call parallelIO_create_dataset_group(pio_dataset_groupid, string, file_id, &
                                          option%iowrite_group_id, ierr)
    grp_id = file_id
#else
    call h5eset_auto_f(OFF,hdf5_err)
    call h5gopen_f(file_id,string,grp_id,hdf5_err)
    if (hdf5_err /= 0) then
      call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
    endif
    call h5eset_auto_f(ON,hdf5_err)
#endif ! PARALLELIO_LIB_WRITE

  ! write out data sets
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)

  select case(option%iflowmode)

    case(FLASH2_MODE,MPH_MODE,THC_MODE,THMC_MODE,MIS_MODE,IMS_MODE, &
         RICHARDS_MODE,G_MODE)

      ! temperature
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          string = "Temperature"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

      ! liquid pressure
      call OutputGetVarFromArray(realization,global_vec,LIQUID_PRESSURE,ZERO_INTEGER)
      string = "Liquid Pressure"
      call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)

      ! gas pressure
      select case(option%iflowmode)
        case (MPH_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_PRESSURE,ZERO_INTEGER)
          string = "Gas Pressure"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,THMC_MODE,RICHARDS_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          string = "Liquid Saturation"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)  
      end select

      ! gas saturation
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE,THC_MODE,THMC_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          string = "Gas Saturation"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

#ifdef ICE
      ! ice saturation
      select case(option%iflowmode)
        case (THC_MODE,THMC_MODE)
          call OutputGetVarFromArray(realization,global_vec,ICE_SATURATION,ZERO_INTEGER)
          string = "Ice Saturation"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

      ! ice density
      select case(option%iflowmode)
        case (THC_MODE,THMC_MODE)
          call OutputGetVarFromArray(realization,global_vec,ICE_DENSITY,ZERO_INTEGER)
          string = "Ice Density"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select
#endif
      
      ! liquid density
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_DENSITY,ZERO_INTEGER)
          string = "Liquid Density"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! gas density
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_DENSITY,ZERO_INTEGER)
          string = "Gas Density"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! liquid viscosity
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_VISCOSITY,ZERO_INTEGER)
          string = "Liquid Viscosity"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! gas viscosity
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_VISCOSITY,ZERO_INTEGER)
          string = "Gas Viscosity"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! liquid mobility
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_MOBILITY,ZERO_INTEGER)
          string = "Liquid Mobility"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! gas mobility
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_MOBILITY,ZERO_INTEGER)
          string = "Gas Mobility"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! liquid energy
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          string = "Liquid Energy"
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! gas energy
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          string = "Gas Energy"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
    
      ! liquid mole fractions
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,THMC_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          do i=1,option%nflowspec
            call OutputGetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i)
            write(string,'(''Liquid Mole Fraction-'',i1)') i
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          enddo
      end select
      
      ! gas mole fractions
      select case(option%iflowmode)
        case (MPH_MODE,FLASH2_MODE,G_MODE)
          do i=1,option%nflowspec
            call OutputGetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i)
            write(string,'(''Gas Mole Fraction-'',i1)') i
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          enddo
      end select
 
      ! phase
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputGetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          string = "Phase"
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,HDF_NATIVE_INTEGER) 
      end select
  
    case default

  end select

  if (option%ntrandof > 0) then
    if (reaction%print_free_conc_type == PRIMARY_MOLALITY) then
      free_mol_char = 'm'
    else
      free_mol_char = 'M'
    endif
    if (reaction%print_tot_conc_type == TOTAL_MOLALITY) then
      tot_mol_char = 'm'
    else
      tot_mol_char = 'M'
    endif
    if (reaction%print_secondary_conc_type == SECONDARY_MOLALITY) then
      tot_mol_char = 'm'
    else
      tot_mol_char = 'M'
    endif
    if (associated(reaction)) then
      if (reaction%print_pH .and. associated(reaction%species_idx)) then
        if (reaction%species_idx%h_ion_id > 0) then
          call OutputGetVarFromArray(realization,global_vec,PH,reaction%species_idx%h_ion_id)
          write(string,'(''pH'')')
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      endif
      if (reaction%print_total_component) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,reaction%print_tot_conc_type,i)
            write(string,'(a,''_tot_'',a)') trim(reaction%primary_species_names(i)), trim(tot_mol_char)
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif
      if (reaction%print_free_ion) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,reaction%print_free_conc_type,i)
            write(string,'(a,''_free_'',a)') trim(reaction%primary_species_names(i)), trim(free_mol_char)
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif
      if (reaction%print_total_bulk) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,TOTAL_BULK,i)
            write(string,'(a,''_total_bulk'')') trim(reaction%primary_species_names(i))
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif
      if (reaction%print_act_coefs) then
        do i=1,reaction%naqcomp
          if (reaction%primary_species_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,PRIMARY_ACTIVITY_COEF,i)
            write(string,'(a)') trim(reaction%primary_species_names(i)) // '_gam'
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif
      do i=1,reaction%neqcplx
        if (reaction%secondary_species_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,reaction%print_secondary_conc_type,i)
          write(string,'(a,a)') trim(reaction%secondary_species_names(i)), trim(sec_mol_char)
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_VOLUME_FRACTION,i)
          write(string,'(a)') trim(reaction%mineral%kinmnrl_names(i)) // '_vf'
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo
      do i=1,reaction%mineral%nkinmnrl
        if (reaction%mineral%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_RATE,i)
          write(string,'(a)') trim(reaction%mineral%kinmnrl_names(i)) // '_rt'
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo
      do i=1,reaction%mineral%nmnrl
        if (reaction%mineral%mnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_SATURATION_INDEX,i)
          write(string,'(a)') trim(reaction%mineral%mineral_names(i)) // '_si'
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_density_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_SITE_DENSITY,i)
          write(string,'(a)') reaction%surface_complexation%srfcplxrxn_site_names(i) // '_den'
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX_FREE,i)
          write(string,'(a)') reaction%surface_complexation%srfcplxrxn_site_names(i)
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo
      do i=1,reaction%surface_complexation%nsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX,i)
          write(string,'(a)') reaction%surface_complexation%srfcplx_names(i)
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo

      ! Kinetic surface complexes

      do i=1,reaction%surface_complexation%nkinsrfcplxrxn
        if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,KIN_SURFACE_CMPLX_FREE,i)
          write(string,'(a)') reaction%surface_complexation%srfcplxrxn_site_names(i)
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo
      do i=1,reaction%surface_complexation%nkinsrfcplx
        if (reaction%surface_complexation%srfcplx_print(i)) then
        !TODO(geh): fix
          call OutputGetVarFromArray(realization,global_vec,KIN_SURFACE_CMPLX,i)
          write(string,'(a)') reaction%surface_complexation%srfcplx_names(i)
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      enddo

     ! Kd

      if (associated(reaction%kd_print)) then
        do i=1,reaction%naqcomp
          if (reaction%kd_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,PRIMARY_KD,i)
            write(string,'(a)') trim(reaction%primary_species_names(i)) // '_kd'
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_print)) then
        do i=1,reaction%naqcomp
          if (reaction%total_sorb_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,TOTAL_SORBED,i)
            write(string,'(a)') trim(reaction%primary_species_names(i)) // '_tot_sorb'
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif
      if (associated(reaction%total_sorb_mobile_print)) then
        do i=1,reaction%ncollcomp
          if (reaction%total_sorb_mobile_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,TOTAL_SORBED_MOBILE,i)
            write(string,'(a)') trim(reaction%colloid_species_names(i)) // '_tot_sorb_mob'
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif
      if (reaction%print_colloid) then
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,COLLOID_MOBILE,i)
            write(string,'(a,''_col_mob_'',a)') trim(reaction%colloid_names(i)), trim(tot_mol_char)
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
        do i=1,reaction%ncoll
          if (reaction%colloid_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,COLLOID_IMMOBILE,i)
            write(string,'(a,''_col_imb_'',a)') trim(reaction%colloid_names(i)), trim(tot_mol_char)
            call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          endif
        enddo
      endif

!     Age
  
      if (reaction%print_age) then
        if (reaction%species_idx%tracer_age_id > 0) then

          call OutputGetVarFromArray(realization,global_vec,AGE, &
            reaction%species_idx%tracer_age_id, &
            reaction%species_idx%tracer_aq_id)
          
          write(string,'(''Tracer_Age'')')
          
          call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        endif
      endif
    endif
  endif
  
  ! porosity
  if (output_option%print_porosity) then
    call OutputGetVarFromArray(realization,global_vec,POROSITY,ZERO_INTEGER)
    string = "Porosity"
    call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec, &
                                       grp_id,H5T_NATIVE_DOUBLE)
  endif
  
  ! material id
  call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
  string = "Material_ID"
  call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id,HDF_NATIVE_INTEGER) 
  
  if (output_option%print_hdf5_velocities) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)
    call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    call OutputGetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id, &
          H5T_NATIVE_DOUBLE)

    if (option%nphase > 1) then
        call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
        string = "Gas X-Velocity"
        call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
        string = "Gas Y-Velocity"
        call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)

        call OutputGetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
        string = "Gas Z-Velocity"
        call HDF5WriteUnstructuredDataSetFromVec(string,realization,global_vec,grp_id, &
            H5T_NATIVE_DOUBLE)
    endif
  endif

  if (output_option%print_hdf5_flux_velocities) then
  
    option%io_buffer = 'HDF5 output of Liquid Flux Velocities for ' // &
      'unstructured grid not supported yet.'

    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
        string = "Liquid X-Flux Velocities"
        !call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,X_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas X-Flux Velocities"
          !call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,X_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%ny > 1) then
        string = "Liquid Y-Flux Velocities"
        !call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,Y_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Y-Flux Velocities"
          !call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,Y_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%nz > 1) then
        string = "Liquid Z-Flux Velocities"
        !call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,Z_DIRECTION,grp_id)
        if (option%nphase > 1) then
          string = "Gas Z-Flux Velocities"
          !call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,Z_DIRECTION,grp_id)
        endif
    endif
   
  endif

  call VecDestroy(global_vec,ierr)

#if defined(PARALLELIO_LIB_WRITE)
!    call parallelio_close_dataset_group(pio_dataset_groupid, file_id, &
!            option%iowrite_group_id, ierr)
!    call parallelio_close_file(file_id, option%iowrite_group_id, ierr)
#else
  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif !PARALLELIO_LIB_WRITE

  hdf5_first = PETSC_FALSE

#endif ! !defined(PETSC_HAVE_HDF5)

end subroutine OutputHDF5UGrid

#if defined(PETSC_HAVE_HDF5)
! ************************************************************************** !
!> This subroutine writes structured coordinates to HDF5 file
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
  
  implicit none
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

#if defined(PARALLELIO_LIB_WRITE)
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

  PetscInt :: local_size
  PetscInt :: i,j
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  PetscInt, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element


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

#if defined(PARALLELIO_LIB_WRITE)
  write(*,*),'PARALLELIO_LIB_WRITE'
  option%io_buffer = 'WriteHDF5CoordinatesUGrid not supported for PARALLELIO_LIB_WRITE'
  call printErrMsg(option)
#else

  !
  !        not(PARALLELIO_LIB_WRITE)
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
#if defined(PARALLELIO_LIB_WRITE)
  write(*,*),'PARALLELIO_LIB_WRITE'
  option%io_buffer = 'WriteHDF5CoordinatesUGrid not supported for PARALLELIO_LIB_WRITE'
  call printErrMsg(option)
#else

  !
  !        not(PARALLELIO_LIB_WRITE)
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
#endif

#ifdef SURFACE_FLOW
! ************************************************************************** !
!> This subroutine is main driver for all output subroutines related to
!! surface flows.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/29/12
! ************************************************************************** !
subroutine Output2(surf_realization,realization,plot_flag,transient_plot_flag)

  use Surface_Realization_module, only : surface_realization_type
  use Realization_module, only : realization_type
  use Option_module, only : OptionCheckTouch, option_type, &
                            output_option_type, printMsg, printErrMsg

  implicit none

  type(surface_realization_type) :: surf_realization
  type(realization_type)         :: realization
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
      option%io_buffer = 'HDF5 output not supported for surface flow'
      call printErrMsg(option)
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

end subroutine Output2

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

  use Surface_Realization_module, only : surface_realization_type
  use Realization_module, only : realization_type
  use Discretization_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Surface_Field_module
  use Patch_module
  
  use Mphase_module
  use Immis_module
  use THC_module
  use THMC_module
  use Richards_module
  use Flash2_module
  use Miscible_module
  use General_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
 
  implicit none

  type(surface_realization_type) :: surf_realization
  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  PetscInt, parameter :: icolumn = -1
  character(len=MAXSTRINGLENGTH) :: filename, string, string2
  character(len=MAXHEADERLENGTH) :: header, header2
  character(len=MAXWORDLENGTH) :: tmp_global_prefix
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(patch_type), pointer :: patch 
  type(reaction_type), pointer :: reaction 
  type(output_option_type), pointer :: output_option
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  Vec :: global_cconn_vec
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  
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
  !call WriteTecplotUGridVertices(OUTPUT_UNIT,surf_realization)
  call WriteTecplotUGridVertices(OUTPUT_UNIT,realization)

  select case(option%iflowmode)
    case(RICHARDS_MODE)
      !pressure
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          call OutputGetVarFromArray(surf_realization,global_vec,LIQUID_PRESSURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(OUTPUT_UNIT,surf_realization,natural_vec,TECPLOT_REAL)
      end select
  end select

  ! material id
  call OutputGetVarFromArray(surf_realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,surf_realization,natural_vec,TECPLOT_INTEGER)

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
subroutine OutputTecplotHeader2(fid,surf_realization,icolumn)

  use Surface_Realization_module
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module

  use Mphase_module
  use Immis_module
  use THC_module
  use THMC_module
  use Richards_module
  use Flash2_module
  use Miscible_module
  use General_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
  use Surface_Flow_Module
  
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
                option%time/output_option%tconv,output_option%tunit

  ! initial portion of header
  header = 'VARIABLES=' // &
            '"X [m]",' // &
            '"Y [m]",' // &
            '"Z [m]"'

  ! write flow variables
  header2 = ''
  select case(option%iflowmode)
    case(RICHARDS_MODE)
      header2 = SurfaceFlowGetTecplotHeader(surf_realization,icolumn)
  end select
  header = trim(header) // trim(header2)

  ! write transport variables
  if (option%ntrandof > 0) then
    string = ''
    !header2 = RTGetTecplotHeader(surf_realization,string,icolumn)
    header = trim(header) // trim(header2)
  endif

  ! add porosity to header
  if (output_option%print_porosity) then
    header = trim(header) // ',"Porosity"'
  endif

  ! write material ids
  header = trim(header) // ',"Material_ID"'

  if (associated(output_option%plot_variables)) then
    do i = 1, size(output_option%plot_variables)
      header = trim(header) // ',"' // &
        trim(PatchGetVarNameFromKeyword(output_option%plot_variables(i), &
                                        option)) // &
        '"'
    enddo
  endif

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

end subroutine OutputTecplotHeader2

! ************************************************************************** !
!> This subroutine prints zone header to Tecplot file.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
function OutputTecplotZoneHeader2(surf_realization,variable_count,tecplot_format)

  use Surface_Realization_module
  use Grid_module
  use Option_module
  
  implicit none

  type(surface_realization_type) :: surf_realization
  PetscInt :: variable_count
  PetscInt :: tecplot_format
  
  character(len=MAXSTRINGLENGTH) :: OutputTecplotZoneHeader2

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
  
  OutputTecplotZoneHeader2 = trim(string) // string2

end function OutputTecplotZoneHeader2

! ************************************************************************** !
!> This subroutine writes unstructure grid vertices for surface grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine WriteTecplotUGridVertices2(fid,surf_realization)

  use Surface_Realization_module
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

end subroutine WriteTecplotUGridVertices2


! ************************************************************************** !
!> This subroutine writes data from an array within a block of a Tecplot file.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine WriteTecplotDataSet2(fid, &
                                surf_realization, &
                                array, &
                                datatype, &
                                size_flag &
                               )

  use Surface_Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(surface_realization_type) :: surf_realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size

  PetscInt, parameter :: num_per_line = 10

  call WriteTecplotDataSetNumPerLine(fid,surf_realization,array,datatype, &
                                     size_flag,num_per_line)
  
end subroutine WriteTecplotDataSet2

! ************************************************************************** !
!> This subroutine writes data from an array within a block of a Tecplot file
!! with a specified number of values per line.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine WriteTecplotDataSetNumPerLine2(fid, &
                                          surf_realization, &
                                          array, &
                                          datatype, &
                                          size_flag, &
                                          num_per_line &
                                         )

  use Surface_Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(surface_realization_type) :: surf_realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  PetscInt :: num_per_line
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscInt :: i
  PetscInt :: max_proc, max_proc_prefetch
  PetscMPIInt :: iproc_mpi, recv_size_mpi
  PetscInt :: max_local_size
  PetscMPIInt :: local_size_mpi
  PetscInt :: istart, iend, num_in_array
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)

  
1000 format(100(i2,1x))
1001 format(100(i4,1x))
1002 format(100(i6,1x))
1003 format(100(i8,1x))
1004 format(100(i10,1x))
1010 format(100(es13.6,1x))

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option

  call PetscLogEventBegin(logging%event_output_write_tecplot,ierr)

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
      call MPI_Allreduce(grid%nlmax,max_local_size,ONE_INTEGER_MPI, &
                         MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
      max_local_size_saved = max_local_size
      write(option%io_buffer,'("max_local_size_saved: ",i9)') max_local_size
      call printMsg(option)
    endif
    max_local_size = max_local_size_saved
    local_size_mpi = grid%nlmax
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

  call PetscLogEventEnd(logging%event_output_write_tecplot,ierr)

end subroutine WriteTecplotDataSetNumPerLine2

! ************************************************************************** !
!> This subroutine writes data from a PETSc Vec within a block of a Tecplot
!! file for surface grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine WriteTecplotDataSetFromVec2( fid, &
                                        surf_realization, &
                                        vec, &
                                        datatype &
                                      )

  use Surface_Realization_module
  
  implicit none

  PetscInt :: fid
  type(surface_realization_type) :: surf_realization
  Vec :: vec
  PetscInt :: datatype
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,surf_realization,vec_ptr,datatype,ZERO_INTEGER) ! 0 implies grid%nlmax
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteTecplotDataSetFromVec2

! ************************************************************************** !
!> This subroutine writes unstructured grid elements for surface grid.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine WriteTecplotUGridElements2(fid, &
                                      surf_realization)

  use Surface_Realization_module
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

end subroutine WriteTecplotUGridElements2

! ************************************************************************** !
!> This subroutine extracts variables indexed by ivar from a multivar array.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/30/12
! ************************************************************************** !
subroutine OutputGetVarFromArray2(surf_realization, &
                                  vec, &
                                  ivar, &
                                  isubvar, &
                                  isubvar1 &
                                  )

  use Surface_Realization_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1

  call PetscLogEventBegin(logging%event_output_get_var_from_array,ierr)

  call SurfaceRealizationGetDataset(surf_realization,vec,ivar,isubvar,isubvar1)

  call PetscLogEventEnd(logging%event_output_get_var_from_array,ierr)

end subroutine OutputGetVarFromArray2

! ************************************************************************** !
!> This routine outputs hydrograph fluxes.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date:
! ************************************************************************** !

subroutine OutputHydrograph(surf_realization)

  use Surface_Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Connection_module
  use Utility_module

  use Surface_Flow_Module
  
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
      write(fid,100,advance="no") option%time/output_option%tconv
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

#endif SURFACE_FLOW

end module Output_module
