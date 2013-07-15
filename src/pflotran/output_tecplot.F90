module Output_Tecplot_module

  use Logging_module 
  use Output_Aux_module
  use Output_Common_module
  
  implicit none

  private

#include "definitions.h"
  PetscInt, parameter, public :: TECPLOT_POINT_FORMAT = 1
  PetscInt, parameter, public :: TECPLOT_BLOCK_FORMAT = 2
  PetscInt, parameter, public :: TECPLOT_FEBRICK_FORMAT = 3
  PetscInt, parameter, public :: TECPLOT_FEQUADRILATERAL_FORMAT = 4  

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petsclog.h"

  public :: OutputTecplotBlock, & 
            OutputTecplotPoint, &
            OutputVelocitiesTecplotBlock, &
            OutputFluxVelocitiesTecplotBlk, &
            OutputVelocitiesTecplotPoint, &
            OutputVectorTecplot, &
            GetCellConnectionsTecplot, &
            WriteTecplotDatasetFromVec, &
            WriteTecplotDatasetNumPerLine, &
            WriteTecplotDataset, &
            OutputPrintExplicitFlowrates, &
            OutputSecondaryContinuumTecplot 

contains

! ************************************************************************** !
!
! OutputTecplotHeader: Print header to Tecplot file
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
subroutine OutputTecplotHeader(fid,realization_base,icolumn)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Patch_module
  
  implicit none

  PetscInt :: fid
  class(realization_base_type) :: realization_base
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
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  output_option => realization_base%output_option

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
  string = OutputTecplotZoneHeader(realization_base,variable_count, &
                                   output_option%tecplot_format)
  write(fid,'(a)') trim(string)

end subroutine OutputTecplotHeader

! ************************************************************************** !
!
! OutputTecplotZoneHeader: Print zone header to Tecplot file
! author: Glenn Hammond
! date: 01/13/12
!
! ************************************************************************** !  
function OutputTecplotZoneHeader(realization_base,variable_count,tecplot_format)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use String_module
  
  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: variable_count
  PetscInt :: tecplot_format
  
  character(len=MAXSTRINGLENGTH) :: OutputTecplotZoneHeader

  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  grid => realization_base%patch%grid
  option => realization_base%option
  output_option => realization_base%output_option

  string = 'ZONE T="' // &
           trim(StringFormatDouble(option%time/output_option%tconv)) // &
           '"'
  string2 = ''
  select case(tecplot_format)
    case (TECPLOT_POINT_FORMAT)
      if ((realization_base%discretization%itype == STRUCTURED_GRID).or. &
          (realization_base%discretization%itype == STRUCTURED_GRID_MIMETIC)) then
        string2 = ', I=' // &
                  trim(StringFormatInt(grid%structured_grid%nx)) // &
                  ', J=' // &
                  trim(StringFormatInt(grid%structured_grid%ny)) // &
                  ', K=' // &
                  trim(StringFormatInt(grid%structured_grid%nz))
      else
        string2 = 'POINT format currently not supported for unstructured'
      endif  
      string2 = trim(string2) // &
              ', DATAPACKING=POINT'
    case default !(TECPLOT_BLOCK_FORMAT,TECPLOT_FEBRICK_FORMAT)
      if ((realization_base%discretization%itype == STRUCTURED_GRID).or. &
          (realization_base%discretization%itype == STRUCTURED_GRID_MIMETIC)) then
        string2 = ', I=' // &
                  trim(StringFormatInt(grid%structured_grid%nx+1)) // &
                  ', J=' // &
                  trim(StringFormatInt(grid%structured_grid%ny+1)) // &
                  ', K=' // &
                  trim(StringFormatInt(grid%structured_grid%nz+1))
      else if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
        string2 = ', N=' // &
                  trim(StringFormatInt(grid%unstructured_grid%num_vertices_global)) // &
                  ', ELEMENTS=' // &
                  trim(StringFormatInt(grid%unstructured_grid%nmax))
        string2 = trim(string2) // ', ZONETYPE=FEBRICK'
      else
        string2 = ', N=' // &
                  trim(StringFormatInt(grid%unstructured_grid%nmax)) // &
                  ', ELEMENTS=' // &
                  trim(StringFormatInt(grid%unstructured_grid%explicit_grid%num_elems))
        string2 = trim(string2) // ', ZONETYPE=FEBRICK'
      endif  
      
      if (grid%itype == EXPLICIT_UNSTRUCTURED_GRID) then
        string3 = ', VARLOCATION=(NODAL)'
      else
        if (variable_count > 4) then
          string3 = ', VARLOCATION=([4-' // &
                    trim(StringFormatInt(variable_count)) // &
                    ']=CELLCENTERED)'
        else
          string3 = ', VARLOCATION=([4]=CELLCENTERED)'
        endif
      endif
      string2 = trim(string2) // trim(string3) // ', DATAPACKING=BLOCK'
    
    end select
  
  OutputTecplotZoneHeader = trim(string) // string2  

end function OutputTecplotZoneHeader

! ************************************************************************** !
!
! OutputTecplotBlock: Print to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !  
subroutine OutputTecplotBlock(realization_base)

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Field_module
  use Patch_module
  
  use Reaction_Aux_module
 
  implicit none

  class(realization_base_type) :: realization_base
  
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
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr
  
  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  
  filename = OutputFilename(output_option,option,'tec','')
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=filename,action="write")
    call OutputTecplotHeader(OUTPUT_UNIT,realization_base,icolumn)
  endif
    
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)  

  ! write out coordinates
  if (realization_base%discretization%itype == STRUCTURED_GRID .or. &
      realization_base%discretization%itype == STRUCTURED_GRID_MIMETIC ) then
    call WriteTecplotStructuredGrid(OUTPUT_UNIT,realization_base)
  else
    call WriteTecplotUGridVertices(OUTPUT_UNIT,realization_base)
  endif

  ! loop over variables and write to file
  cur_variable => output_option%output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputGetVarFromArray(realization_base,global_vec,cur_variable%ivar, &
                                cur_variable%isubvar)
    call DiscretizationGlobalToNatural(discretization,global_vec, &
                                        natural_vec,ONEDOF)
    if (cur_variable%iformat == 0) then
      call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec, &
                                      TECPLOT_REAL)
    else
      call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec, &
                                      TECPLOT_INTEGER)
    endif
    cur_variable => cur_variable%next
  enddo

  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if (realization_base%discretization%itype == UNSTRUCTURED_GRID .and. &
      realization_base%discretization%grid%itype == &
      IMPLICIT_UNSTRUCTURED_GRID)  then
    call WriteTecplotUGridElements(OUTPUT_UNIT,realization_base)
  endif
  
  if (realization_base%discretization%grid%itype ==  &
        EXPLICIT_UNSTRUCTURED_GRID) then
    call WriteTecplotExpGridElements(OUTPUT_UNIT,realization_base)
  endif
    
  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)
  
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesTecplotBlock(realization_base)
  endif
  
  if (output_option%print_tecplot_flux_velocities) then
    if (grid%structured_grid%nx > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization_base,LIQUID_PHASE, &
                                          X_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization_base,GAS_PHASE, &
                                              X_DIRECTION)
      end select
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization_base,LIQUID_PHASE, &
                                          Y_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE, IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization_base,GAS_PHASE, &
                                              Y_DIRECTION)
      end select
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization_base,LIQUID_PHASE, &
                                          Z_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE, IMS_MODE,FLASH2_MODE,G_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization_base,GAS_PHASE, &
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
subroutine OutputVelocitiesTecplotBlock(realization_base)
 
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Field_module
  use Patch_module
  use Variables_module
  
  implicit none

  class(realization_base_type) :: realization_base
  
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
  PetscErrorCode :: ierr

  PetscReal, pointer :: vec_ptr(:)
  
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field
  option => realization_base%option
  output_option => realization_base%output_option
  discretization => realization_base%discretization

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
      string = OutputTecplotZoneHeader(realization_base,TEN_INTEGER, &
                                       TECPLOT_BLOCK_FORMAT)
    else
      string = OutputTecplotZoneHeader(realization_base,SEVEN_INTEGER, &
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
  if (realization_base%discretization%itype == STRUCTURED_GRID .or. &
      realization_base%discretization%itype == STRUCTURED_GRID_MIMETIC)  then
    call WriteTecplotStructuredGrid(OUTPUT_UNIT,realization_base)
  else
    call WriteTecplotUGridVertices(OUTPUT_UNIT,realization_base)
  endif
  
  call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,X_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec,TECPLOT_REAL)

  call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,Y_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec,TECPLOT_REAL)

  call OutputGetCellCenteredVelocities(realization_base,global_vec,LIQUID_PHASE,Z_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec,TECPLOT_REAL)

  if (option%nphase > 1) then
    call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,X_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec,TECPLOT_REAL)

    call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,Y_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec,TECPLOT_REAL)

    call OutputGetCellCenteredVelocities(realization_base,global_vec,GAS_PHASE,Z_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec,TECPLOT_REAL)
  endif

  ! material id
  call OutputGetVarFromArray(realization_base,global_vec,MATERIAL_ID,ZERO_INTEGER)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(OUTPUT_UNIT,realization_base,natural_vec,TECPLOT_INTEGER)
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if (realization_base%discretization%itype == UNSTRUCTURED_GRID .and. &
      realization_base%discretization%grid%itype == &
      IMPLICIT_UNSTRUCTURED_GRID)  then
    call WriteTecplotUGridElements(OUTPUT_UNIT,realization_base)
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
subroutine OutputFluxVelocitiesTecplotBlk(realization_base,iphase, &
                                            direction)
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
  PetscErrorCode :: ierr

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
    
  nullify(array)

  call PetscLogEventBegin(logging%event_output_write_flux_tecplot,ierr) 
                          
  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  
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
  call WriteTecplotDataSet(OUTPUT_UNIT,realization_base,array,TECPLOT_REAL,adjusted_size)
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
  call WriteTecplotDataSet(OUTPUT_UNIT,realization_base,array,TECPLOT_REAL,adjusted_size)
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
  call WriteTecplotDataSet(OUTPUT_UNIT,realization_base,array,TECPLOT_REAL,adjusted_size)
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
  call WriteTecplotDataSet(OUTPUT_UNIT,realization_base,array,TECPLOT_REAL,adjusted_size)
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
subroutine OutputTecplotPoint(realization_base)

  use Realization_Base_class, only : realization_base_type, &
                                     RealizGetDatasetValueAtCell
  use Discretization_module
  use Grid_module
  use Structured_Grid_module
  use Option_module
  use Field_module
  use Patch_module

  use Reaction_Aux_module
 
  implicit none

  class(realization_base_type) :: realization_base
  
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
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: value
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr  
  
  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

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
    call OutputTecplotHeader(OUTPUT_UNIT,realization_base,icolumn)
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
      value = RealizGetDatasetValueAtCell(realization_base,cur_variable%ivar, &
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
    call OutputVelocitiesTecplotPoint(realization_base)
  endif
  
end subroutine OutputTecplotPoint

! ************************************************************************** !
!
! OutputVelocitiesTecplotPoint: Print velocities to Tecplot file in POINT format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVelocitiesTecplotPoint(realization_base)
 
  use Realization_Base_class, only : realization_base_type, &
                                     RealizGetDatasetValueAtCell
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  use Variables_module
  
  implicit none

  class(realization_base_type) :: realization_base
  
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
  PetscErrorCode :: ierr  

  PetscReal, pointer :: vec_ptr_vx(:), vec_ptr_vy(:), vec_ptr_vz(:)
  
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field
  option => realization_base%option
  output_option => realization_base%output_option
  discretization => realization_base%discretization
  
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
  
  call OutputGetCellCenteredVelocities(realization_base,global_vec_vx,LIQUID_PHASE,X_DIRECTION)
  call OutputGetCellCenteredVelocities(realization_base,global_vec_vy,LIQUID_PHASE,Y_DIRECTION)
  call OutputGetCellCenteredVelocities(realization_base,global_vec_vz,LIQUID_PHASE,Z_DIRECTION)

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
    value = RealizGetDatasetValueAtCell(realization_base,MATERIAL_ID, &
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
subroutine OutputVectorTecplot(filename,dataset_name,realization_base,vector)
 
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Field_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Patch_module
  use Variables_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: dataset_name
  class(realization_base_type) :: realization_base
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
  PetscErrorCode :: ierr  

  call PetscLogEventBegin(logging%event_output_vec_tecplot,ierr) 

  option => realization_base%option
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field
  discretization => realization_base%discretization
  
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
    string = OutputTecplotZoneHeader(realization_base,FIVE_INTEGER, &
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

  if (realization_base%discretization%itype == STRUCTURED_GRID .or. &
      realization_base%discretization%itype == STRUCTURED_GRID_MIMETIC)  then
    call WriteTecplotStructuredGrid(fid,realization_base)
  else  
    call WriteTecplotUGridVertices(fid,realization_base)
  endif    

  call DiscretizationGlobalToNatural(discretization,vector,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization_base,natural_vec,TECPLOT_REAL)

  call OutputGetVarFromArray(realization_base,global_vec,MATERIAL_ID,ZERO_INTEGER)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization_base,natural_vec,TECPLOT_INTEGER)
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if (realization_base%discretization%itype == UNSTRUCTURED_GRID .and. &
      realization_base%discretization%grid%itype == &
      IMPLICIT_UNSTRUCTURED_GRID)  then
    call WriteTecplotUGridElements(fid,realization_base)
  endif    

  close(fid)

  call PetscLogEventEnd(logging%event_output_vec_tecplot,ierr) 
                            
end subroutine OutputVectorTecplot


! ************************************************************************** !
!
! WriteTecplotStructuredGrid: Writes structured grid face coordinates 
! author: Glenn Hammond
! date: 02/26/08
!
! ************************************************************************** !
subroutine WriteTecplotStructuredGrid(fid,realization_base)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscInt :: i, j, k, count, nx, ny, nz
  PetscReal :: temp_real
  PetscErrorCode :: ierr  

1000 format(es13.6,1x)
1001 format(10(es13.6,1x))
  
  call PetscLogEventBegin(logging%event_output_str_grid_tecplot,ierr) 
                              
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  
  nx = grid%structured_grid%nx
  ny = grid%structured_grid%ny
  nz = grid%structured_grid%nz
  
  if (option%myrank == option%io_rank) then
    ! x-dir
    count = 0
    do k=1,nz+1
      do j=1,ny+1
        temp_real = realization_base%discretization%origin(X_DIRECTION)
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
      temp_real = realization_base%discretization%origin(Y_DIRECTION)
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
    temp_real = realization_base%discretization%origin(Z_DIRECTION)
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
! WriteTecplotUGridVertices: Writes unstructured grid vertices
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine WriteTecplotUGridVertices(fid,realization_base)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module
  use Variables_module

  implicit none

  PetscInt :: fid
  class(realization_base_type) :: realization_base
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch 
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  PetscInt :: local_size
  PetscErrorCode :: ierr
  PetscInt :: num_cells, icell
  PetscInt :: count
    
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option

1000 format(es13.6,1x)

  if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
    call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                      grid%unstructured_grid%num_vertices_global, &
                      global_vertex_vec,ierr)
    call VecGetLocalSize(global_vertex_vec,local_size,ierr)
    call GetVertexCoordinates(grid, global_vertex_vec,X_COORDINATE,option)
    call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
    call WriteTecplotDataSet(fid,realization_base,vec_ptr,TECPLOT_REAL, &
                             local_size)
    call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

    call GetVertexCoordinates(grid,global_vertex_vec,Y_COORDINATE,option)
    call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
    call WriteTecplotDataSet(fid,realization_base,vec_ptr,TECPLOT_REAL, &
                             local_size)
    call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

    call GetVertexCoordinates(grid,global_vertex_vec, Z_COORDINATE,option)
    call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr)
    call WriteTecplotDataSet(fid,realization_base,vec_ptr,TECPLOT_REAL, &
                             local_size)
    call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr)

    call VecDestroy(global_vertex_vec, ierr)
  else
    if (option%myrank == option%io_rank) then
      if (option%print_explicit_primal_grid) then
      num_cells = grid%unstructured_grid%explicit_grid%num_cells_global
      count = 0
      do icell = 1, num_cells
        write(fid,1000,advance='no') grid%unstructured_grid%explicit_grid% &
                                     vertex_coordinates(icell)%x
        count = count + 1
        if (mod(count,10) == 0) then
          write(fid,'(a)') ""
          count = 0
        endif
      enddo
      if (count /= 0) write(fid,'(a)') ""
      count = 0
      do icell = 1, num_cells
        write(fid,1000,advance='no') grid%unstructured_grid%explicit_grid% &
                                     vertex_coordinates(icell)%y
        count = count + 1
        if (mod(count,10) == 0) then
          write(fid,'(a)') ""
          count = 0
        endif
      enddo
      if (count /= 0) write(fid,'(a)') ""
      count = 0
      do icell = 1, num_cells
        write(fid,1000,advance='no') grid%unstructured_grid%explicit_grid% &
                                     vertex_coordinates(icell)%z
        count = count + 1
        if (mod(count,10) == 0) then
          write(fid,'(a)') ""
          count = 0
        endif
      enddo
      if (count /= 0) write(fid,'(a)') ""      
      elseif (option%print_explicit_dual_grid) then
        write(fid,'(">",/,"Add explicit mesh vertex information here",/,">")')
      else 
        write(fid,'(">",/,"Add explicit mesh vertex information here",/,">")')
      endif
    endif
  endif

end subroutine WriteTecplotUGridVertices


! ************************************************************************** !
!
! WriteTecplotExpGridElements: Writes unstructured explicit grid elements
! author: Satish Karra, LANL
! date: 04/11/13
!
! ************************************************************************** !
subroutine WriteTecplotExpGridElements(fid,realization_base)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module
  
  implicit none

  PetscInt :: fid
  class(realization_base_type) :: realization_base

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch 
  PetscInt, pointer :: temp_int(:)
  PetscInt :: iconn, num_elems, i, num_vertices
  PetscErrorCode :: ierr
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  
  num_elems = grid%unstructured_grid%explicit_grid%num_elems
 
  allocate(temp_int(grid%unstructured_grid%max_nvert_per_cell))
  
  if (option%myrank == option%io_rank) then
    do iconn = 1, num_elems
      num_vertices = grid%unstructured_grid%explicit_grid% &
                       cell_connectivity(0,iconn)
      select case(num_vertices)
        case(EIGHT_INTEGER) ! Hex mesh
          temp_int = grid%unstructured_grid%explicit_grid% &
                       cell_connectivity(1:num_vertices,iconn)
        case(SIX_INTEGER)   ! Wedge 
          temp_int(1) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(1,iconn)         
          temp_int(2) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(1,iconn)  
          temp_int(3) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(4,iconn) 
          temp_int(4) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(4,iconn)
          temp_int(5) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(3,iconn) 
          temp_int(6) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(2,iconn) 
          temp_int(7) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(5,iconn) 
          temp_int(8) = grid%unstructured_grid%explicit_grid% &
                          cell_connectivity(6,iconn) 
        case(FIVE_INTEGER)  ! Pyramid
          do i = 1, 4
            temp_int(i) = grid%unstructured_grid%explicit_grid% &
                            cell_connectivity(i,iconn) 
          enddo
          do i = 5, 8
            temp_int(i) = grid%unstructured_grid%explicit_grid% &
                            cell_connectivity(5,iconn) 
          enddo
        case(FOUR_INTEGER)
          if (grid%unstructured_grid%grid_type == TWO_DIM_GRID) then ! Quad
            do i = 1, 4
              temp_int(i) = grid%unstructured_grid%explicit_grid% &
                              cell_connectivity(i,iconn) 
            enddo
            do i = 5, 8
              temp_int(i) = temp_int(i-4)
            enddo
          else ! Tet
            do i = 1, 3
              temp_int(i) = grid%unstructured_grid%explicit_grid% &
                             cell_connectivity(i,iconn) 
            enddo
            temp_int(4) = temp_int(3)
            do i = 5, 8
              temp_int(i) = grid%unstructured_grid%explicit_grid% &
                              cell_connectivity(4,iconn) 
            enddo
          endif
        case(3) ! Tri
          do i = 1, 3
            temp_int(i) = grid%unstructured_grid%explicit_grid% &
                            cell_connectivity(i,iconn) 
          enddo
          temp_int(4) = temp_int(3)
          do i = 5, 8
            temp_int(i) = temp_int(i-4) 
          enddo
      end select
      write(fid,*) temp_int
    enddo 
  endif
  
  deallocate(temp_int)
   
end subroutine WriteTecplotExpGridElements

! ************************************************************************** !
!
! WriteTecplotUGridElements: Writes unstructured grid elements
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine WriteTecplotUGridElements(fid,realization_base)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Option_module
  use Patch_module
  
  implicit none

  PetscInt :: fid
  class(realization_base_type) :: realization_base

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch 
  Vec :: global_cconn_vec
  type(ugdm_type), pointer :: ugdm_element
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  Vec :: global_vec
  Vec :: natural_vec

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  
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
  call WriteTecplotDataSetNumPerLine(fid,realization_base,vec_ptr, &
                                     TECPLOT_INTEGER, &
                                     grid%unstructured_grid%nlmax*8, &
                                     EIGHT_INTEGER)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)
  call VecDestroy(global_vec,ierr)
  call VecDestroy(natural_vec,ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine WriteTecplotUGridElements

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
  use Unstructured_Cell_module
  
  implicit none
  
  type(grid_type) :: grid
  type(unstructured_grid_type),pointer :: ugrid
  Vec :: vec
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr  
  
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
!
! WriteTecplotDataSetFromVec: Writes data from a Petsc Vec within a block
!                             of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSetFromVec(fid,realization_base,vec,datatype)

  use Realization_Base_class, only : realization_base_type
  
  implicit none

  PetscInt :: fid
  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: datatype
  PetscErrorCode :: ierr  
  
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)                ! 0 implies grid%nlmax
  call WriteTecplotDataSet(fid,realization_base,vec_ptr,datatype,ZERO_INTEGER) 
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteTecplotDataSetFromVec

! ************************************************************************** !
!
! WriteTecplotDataSet: Writes data from an array within a block
!                      of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSet(fid,realization_base,array,datatype,size_flag)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Patch_module

  implicit none

  PetscInt :: fid
  class(realization_base_type) :: realization_base
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size

  PetscInt, parameter :: num_per_line = 10

  call WriteTecplotDataSetNumPerLine(fid,realization_base,array,datatype, &
                                     size_flag,num_per_line) 
  
end subroutine WriteTecplotDataSet

! ************************************************************************** !
!
! WriteTecplotDataSetNumPerLine: Writes data from an array within a block
!                                of a Tecplot file with a specified number
!                                of values per line
! author: Glenn Hammond
! date: 10/25/07, 12/02/11
!
! ************************************************************************** !
subroutine WriteTecplotDataSetNumPerLine(fid,realization_base,array,datatype, &
                                         size_flag,num_per_line)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
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
  PetscErrorCode :: ierr  
  
1000 format(100(i2,1x))
1001 format(100(i4,1x))
1002 format(100(i6,1x))
1003 format(100(i8,1x))
1004 format(100(i10,1x))
1010 format(100(es13.6,1x))

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option

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

end subroutine WriteTecplotDataSetNumPerLine

! ************************************************************************** !
!
! OutputPrintExplicitFlowrates: Prints out the flow rate through a voronoi face
! for explicit grid. This will be used for particle tracking.
! outputs in the format: mag(Res) where the row number is same as connection
! number
! author: Satish Karra, LANL
! date: 04/24/13
!
! ************************************************************************** !
subroutine OutputPrintExplicitFlowrates(realization_base)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  use Output_Common_module  
 
  implicit none

  class(realization_base_type) :: realization_base
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr  
  PetscInt :: iconn
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  filename = OutputFilename(output_option,option,'dat','rates')
  
  call OutputGetExplicitFlowrates(realization_base)
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write rate output file: ' // &
                       trim(filename)
    call printMsg(option)                       
    open(unit=OUTPUT_UNIT,file=filename,action="write")
  endif
  
1000 format(es13.6,1x)
1009 format('')

  call VecGetArrayF90(field%flowrate_inst,vec_ptr,ierr)
  if (option%myrank == option%io_rank) then
    do iconn = 1,size(grid%unstructured_grid%explicit_grid%connections,2)
      write(OUTPUT_UNIT,1000,advance='no') vec_ptr(iconn)
      write(OUTPUT_UNIT,1009) 
    enddo
  endif
  call VecRestoreArrayF90(field%flowrate_inst,vec_ptr,ierr)
    
  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)

end subroutine OutputPrintExplicitFlowrates


! ************************************************************************** !
!
! OutputSecondaryContinuumTecplot: Print secondary continuum variables
! in tecplot format. The output is at a given primary continuum node, 
! and the coordinates in the output are the secondary continuum spatial 
! coordinates 
! author: Satish Karra, LANL
! date: 04/30/2013
!
! ************************************************************************** !  
subroutine OutputSecondaryContinuumTecplot(realization_base)

  use Realization_Base_class, only : realization_base_type, &
                                     RealizGetDatasetValueAtCell
  use Option_module
  use Field_module
  use Patch_module
  use Grid_module
  use Reaction_Aux_module
  use Observation_module
  use Variables_module
  use Secondary_Continuum_Aux_module, only : sec_transport_type, &
                                             sec_heat_type, sec_continuum_type

 
  implicit none

  class(realization_base_type) :: realization_base
  
  PetscInt :: i, comma_count, quote_count
  PetscInt :: icolumn
  character(len=MAXSTRINGLENGTH) :: filename, string, string2
  character(len=MAXSTRINGLENGTH) :: string3
  character(len=MAXHEADERLENGTH) :: header
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch 
  type(output_option_type), pointer :: output_option
  type(observation_type), pointer :: observation
  type(grid_type), pointer :: grid
  type(sec_transport_type), pointer :: rt_sec_tranport_vars(:)
  type(sec_heat_type), pointer :: sec_heat_vars(:)
  type(reaction_type), pointer :: reaction   
  PetscReal :: value
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr  
  PetscInt :: count, icell, fid, sec_id
  PetscInt :: ghosted_id, local_id
  PetscInt :: naqcomp, nkinmnrl
  PetscReal, pointer :: dist(:)
  
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  grid => patch%grid
  output_option => realization_base%output_option

  if (option%use_mc) then
    if (option%ntrandof > 0) then
      rt_sec_tranport_vars => patch%aux%SC_RT%sec_transport_vars
      reaction => realization_base%reaction
    endif
    if (option%iflowmode == TH_MODE .or. option%iflowmode == THC_MODE &
        .or. option%iflowmode == MPH_MODE) then
      sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
    endif
  endif

  ! Here we are assuming that if there are secondary continua for both
  ! heat and reactive transport, then the shape and type of secondary
  ! continua are the same - SK
  if (associated(sec_heat_vars)) then
    dist => sec_heat_vars(1)%sec_continuum%distance
  elseif (associated(rt_sec_tranport_vars)) then
    dist => rt_sec_tranport_vars(1)%sec_continuum%distance
  endif


  ! write points
1000 format(es13.6,1x)
1009 format('')

  count = 0
  observation => patch%observation%first
  do 
    if (.not.associated(observation)) exit
    write(string,'(i6)') option%myrank
    write(string2,'(i6)') count
    string3 = OutputFilenameID(output_option,option)
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-sec-rank' // trim(adjustl(string)) // '-obs' &
               // trim(adjustl(string2)) // '-' // trim(string3) // '.tec'    
    
    if (option%myrank == option%io_rank) then
      option%io_buffer = '--> write tecplot output file: ' // trim(filename)
      call printMsg(option)
    endif
    
    ! open file
    fid = 86
    open(unit=fid,file=filename,action="write")

    ! must initialize icolumn here so that icolumn does not restart with
    ! each observation point
    if (output_option%print_column_ids) then
      icolumn = 1
    else
      icolumn = -1
    endif    
    
    ! write header
    ! write title
    write(fid,'(''TITLE = "'',1es13.5," [",a1,'']"'')') &
              option%time/output_option%tconv,output_option%tunit

    ! initial portion of header
    header = 'VARIABLES=' // &
              '"dist [m]"'
               
    write(fid,'(a)',advance='no') trim(header)
                      
    if (associated(observation%region%coordinates) .and. &
            .not.observation%at_cell_center) then
      option%io_buffer = 'Writing of data at coordinates not ' // &
              'functioning properly for minerals.  Perhaps due to ' // &
              'non-ghosting of vol frac....>? - geh'
      call printErrMsg(option)
      call WriteTecplotHeaderForCoordSec(fid,realization_base, &
                                         observation%region, &
                                         observation% &
                                         print_secondary_data, &
                                         icolumn)
    else
      do icell = 1,observation%region%num_cells
        call WriteTecplotHeaderForCellSec(fid,realization_base, &
                                          observation%region,icell, &
                                          observation% &
                                          print_secondary_data, &
                                          icolumn)
      enddo
    endif

    write(fid,'(a)',advance='yes') ""
    ! write zone header
    write(string,'(''ZONE T= "'',1es13.5,''",'','' I='',i5)') &
                  option%time/output_option%tconv, &
                  option%nsec_cells
    string = trim(string) // ',J=1, K=1, DATAPACKING=POINT'
    write(fid,'(a)',advance='no') trim(string)     
    write(fid,1009)
   
    do sec_id = 1,option%nsec_cells
      write(fid,1000,advance='no') dist(sec_id)  
      do icell = 1,observation%region%num_cells
        local_id = observation%region%cell_ids(icell)
        ghosted_id = grid%nL2G(local_id)
        if (observation%print_secondary_data(1)) then
          write(fid,1000,advance='no') &
          RealizGetDatasetValueAtCell(realization_base,SECONDARY_TEMPERATURE, &
                                      sec_id,ghosted_id)        
        endif
        if (observation%print_secondary_data(2)) then
          do naqcomp = 1, reaction%naqcomp
            write(fid,1000,advance='no') &
            RealizGetDatasetValueAtCell(realization_base, &
                                        SECONDARY_CONCENTRATION, &
                                        sec_id,ghosted_id,naqcomp)
          enddo
        endif
        if (observation%print_secondary_data(3)) then
          do nkinmnrl = 1, reaction%mineral%nkinmnrl
            write(fid,1000,advance='no') &
            RealizGetDatasetValueAtCell(realization_base,SEC_MIN_VOLFRAC, &
                                        sec_id,ghosted_id,nkinmnrl) 
          enddo
        endif        
      enddo
      write(fid,1009)
    enddo         
       
    close(fid)  
    observation => observation%next
    count = count + 1    
  enddo
   
end subroutine OutputSecondaryContinuumTecplot

! ************************************************************************** !
!
! WriteTecplotHeaderForCellSec: Print tecplot header for data at a cell for
! secondary continuum 
! author: Satish Karra, LANL
! date: 04/30/2013
!
! ************************************************************************** !  
subroutine WriteTecplotHeaderForCellSec(fid,realization_base,region,icell, &
                                        print_secondary_data, &
                                        icolumn)

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Output_Aux_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
  PetscInt :: icell
  PetscBool :: print_secondary_data(3)
  PetscInt :: icolumn
  
  PetscInt :: local_id
  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  type(grid_type), pointer :: grid

  grid => realization_base%patch%grid
  
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
  
  call WriteTecplotHeaderSec(fid,realization_base,cell_string, &
                                 print_secondary_data,icolumn)

end subroutine WriteTecplotHeaderForCellSec

! ************************************************************************** !
!
! WriteTecplotHeaderForCoordSec: Print a header for data at a coordinate
! for secondary continuum
! author: Satish Karra, LANL
! date: 04/30/2013
!
! ************************************************************************** !  
subroutine WriteTecplotHeaderForCoordSec(fid,realization_base,region, &
                                         print_secondary_data, &
                                         icolumn)

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
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

  call WriteTecplotHeaderSec(fid,realization_base,cell_string, &
                             print_secondary_data,icolumn)

end subroutine WriteTecplotHeaderForCoordSec

! ************************************************************************** !
!
! WriteTecplotHeaderSec: Print a header for secondary continuum data
! author: Satish Karra, LANL
! date: 04/30/2013
!
! ************************************************************************** !  
subroutine WriteTecplotHeaderSec(fid,realization_base,cell_string, &
                                 print_secondary_data,icolumn)
                                     
  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(reaction_type), pointer :: reaction 
  PetscBool :: print_secondary_data(3)
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscInt :: icolumn
  
  PetscInt :: i,j
  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option  
  
  option => realization_base%option
  output_option => realization_base%output_option
  
  ! add secondary temperature to header
  if (print_secondary_data(1)) then
    select case (option%iflowmode) 
      case (TH_MODE,THC_MODE, MPH_MODE)
        header = ''
        string = 'T'
        call OutputAppendToHeader(header,string,'C',cell_string,icolumn)
      case default
        header = ''
    end select
    write(fid,'(a)',advance="no") trim(header)
  endif
  
  ! add secondary concentrations to header
  if (option%ntrandof > 0) then 
    reaction => realization_base%reaction
    if (print_secondary_data(2)) then
      header = ''
      do j = 1, reaction%naqcomp
        string = 'Free ion ' // trim(reaction%primary_species_names(j))
        call OutputAppendToHeader(header,string,'molal',cell_string, &
                                  icolumn)
      enddo
      write(fid,'(a)',advance="no") trim(header)
    endif
  
  
    ! add secondary mineral volume fractions to header
    if (print_secondary_data(3)) then
      header = ''
      do j = 1, reaction%naqcomp
        string = trim(reaction%mineral%mineral_names(j)) // ' VF'
        call OutputAppendToHeader(header,string,'',cell_string, &
                                  icolumn)
      enddo
      write(fid,'(a)',advance="no") trim(header)
    endif
  endif 
  
end subroutine WriteTecplotHeaderSec

end module Output_Tecplot_module
