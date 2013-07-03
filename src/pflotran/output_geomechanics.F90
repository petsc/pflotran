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

#if 0
  ! write out coordinates
  call WriteTecplotUGridVertices(OUTPUT_UNIT,surf_realization)

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


end module Output_Geomechanics_module

#endif
